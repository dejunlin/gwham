#include <limits>
#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <iterator>
#include <stdio.h>
#include <stdlib.h>
#include <functional>
#include <iomanip>
#include "typedefs.hpp"
#include "hamiltonian.hpp"
#include "ensemble.hpp"
#include "gnarray.hpp"
#include "fileio.hpp"
#include "fileio_utils.hpp"
#include "gmxmdp.hpp"
#include "ensemble_factory.hpp"
#include "mdp_factory.hpp"
#include "gwham.hpp"
#include "metaprog_snippets.hpp"

using namespace std;

typedef gnarray<coordtype, valtype, valtype> narray;
typedef gnarray<coordtype, histcounter, valtype> histogram;

//! read histogram from file and return the total number of
//elements in the histogram
inline
#ifdef __GNUC__
__attribute__((always_inline))
#endif
histcounter readhist(const string fnhist, histogram& hist) {
  fileio fio_hist(fnhist, std::fstream::in, false);
  const auto ndim = hist.getdim();
  fio_hist.skipemptylns();
  histcounter N = 0;
  while(fio_hist.readtsnb()) {
    const auto cols = fio_hist.line2val<string>();
    //the format should be coord, val, counts, ... (the rest are ignored)
    if(cols.size() < ndim*2 + 1) continue;
    const vector<string> coordstr(cols.begin(), cols.begin()+ndim);
    histogram::gridcoord coord;
    strto(coord, coordstr);
    const histogram::gridval n = stoul(cols[ndim*2]);
    hist[coord] = n;
    N += n;
  }
  return N;
}

inline
#ifdef __GNUC__
__attribute__((always_inline))
#endif
void HistogramCV(histogram& hist, histcounter& N, const vector<valtype>& input) {
  if(hist.bin(input) != hist.end()) ++N; 
}

inline
#ifdef __GNUC__
__attribute__((always_inline))
#endif
void HistogramExpECV(vector<histogram>& hists, vector<histcounter>& Ns, vector<uint>::const_iterator& itstate, const linecounter nststates, linecounter& count, vector<valtype>::const_iterator& itpot, const vector<valtype>& input) {
  const auto& state = *itstate;
  ++count;
  if(!(count % nststates)) ++itstate;
  auto& hist = hists[state];
  const auto& pot = *itpot;
  ++itpot;
  vector<valtype> data = input;
  data.push_back(pot);
  HistogramCV(hist, Ns[state], data);
}

inline
#ifdef __GNUC__
__attribute__((always_inline))
#endif
void HistogramExpCV(vector<histogram>& hists, vector<histcounter>& Ns, vector<uint>::const_iterator& itstate, const linecounter nststates, linecounter& count, const vector<valtype>& input) {
  const auto& state = *itstate;
  ++count;
  if(!(count % nststates)) ++itstate;
  auto& hist = hists[state];
  HistogramCV(hist, Ns[state], input);
}

inline
#ifdef __GNUC__
__attribute__((always_inline))
#endif
void HistogramECV(histogram& hist, histcounter& N, vector<valtype>::const_iterator& itpot, const vector<valtype>& input) {
  const auto& pot = *itpot;
  ++itpot;
  vector<valtype> data = input;
  data.push_back(pot);
  HistogramCV(hist, N, data);
}

template < class T, class TS >
inline
#ifdef __GNUC__
__attribute__((always_inline))
#endif
void vecfromTS(vector<T>& ans, TS& ts, uint nrun, string fnprefix) {
  const auto& fnsuffix = ts.fnsuffix;
  for(uint i = 0; i < nrun; ++i) {
    const string fts = fnprefix + "_" + tostr(i) + fnsuffix;
    ts(fts, [&ans] (const vector<valtype>& input) 
            {
              ans.emplace_back(input.front());
            }
      );
  }
}

template < class vTSit, class PMDP, class H >
inline
#ifdef __GNUC__
__attribute__((always_inline))
#endif
void histfromTS (
    vTSit& it,
    vector<histcounter>& Nsamples,
    const PMDP& pmdp, 
    const multimap<PMDP, uint>& pmdp2ihist, 
    vector<H>& hists, 
    const uint& nrun, 
    string fnprefix, 
    const vector<uint>& states, 
    const vector<valtype>& pot) 
{
  using namespace std::placeholders;

  auto itstate = states.cbegin(); 
  auto itpot = pot.cbegin();

  function<void(const vector<valtype>&)> histogrammer;
  linecounter count = 0;
  if(states.size()) {
    //Here we assume the states time-series is the previous one or the one 
    //before if we have to deal with potential energy
    auto prevts = it;
    --prevts;
    if(pot.size()) { --prevts; }
    const linecounter nststates = linecounter(prevts->nst * prevts->fio.ls / (it->nst * it->fio.ls));
    if(pot.size()) {
	histogrammer = std::bind(&HistogramExpECV, ref(hists), ref(Nsamples), ref(itstate), nststates, ref(count), ref(itpot), _1); 
    } else {
	histogrammer = std::bind(&HistogramExpCV, ref(hists), ref(Nsamples), ref(itstate), nststates, ref(count), _1); 
    }
  } else {
    const uint ihist = pmdp2ihist.find(pmdp)->second;
    if(pot.size()) {
	histogrammer = std::bind(&HistogramECV, ref(hists[ihist]), ref(Nsamples[ihist]), ref(itpot), _1); 
    } else {
	histogrammer = std::bind(&HistogramCV, ref(hists[ihist]), ref(Nsamples[ihist]), _1); 
    }
  }

  auto& ts = *it;
  const auto& fnsuffix = ts.fnsuffix;
  for(uint i = 0; i < nrun; ++i) {
    const string fts = fnprefix + "_" + tostr(i) + fnsuffix;
    ts(fts, histogrammer); 
  }
}

int main(int argc, char* argv[]) {
#ifdef MPREALCXX
  mpreal::set_default_prec(mpfr::digits2bits(MPREAL_PRECISION));
#endif
  string cmdline("#Usage: gwham sysname mdpsuffixes nwin nrun rcprint nbins hv lv tol rcvbegin rcvstride rcvend mdp0prefixes fseeds ifmin ifreadhist\n");
  if (argc != 17) {
    cerr << cmdline;
    exit(-1);
  }
  uint k = 1;
  //name of the system
  //(your files should be named as "sysname_WINID_RUNID*.*")
  const string sysname = string(argv[k++]); 
  //file names suffixes for each type of input files ( the number of file types depends on the system mdp setup) 
  const vector<string> mdpsuffixes = split<string>(argv[k++]); 
  //# of replica 
  const uint nwin = atoi(argv[k++]);
  //# of runs performed
  const uint nrun = atoi(argv[k++]);
  //which pull-group/RC you want the PMF on.
  const vector<uint> rcprt = split<uint>(argv[k++]);
  //# of bins in each dimension
  const vector<uint> nbins = split<uint>(argv[k++]);
  //upper bounds in each dimension
  const vector<valtype> hv = split<valtype>(argv[k++]);
  //lower bounds in each dimension
  const vector<valtype> lv = split<valtype>(argv[k++]); 
  //tolerance for WHAM iteration
  const valtype tol = atof(argv[k++]);
  //first time point to read (excluding comment-lines) in the x.xvg files
  //Here I make the assumption that the first datum
  //corresponds to the starting point and should always be skipped
  const linecounter rcvbegin = stoul(argv[k++]);
  //Only read every this number of time points (excluding comment-lines) in the x.xvg files
  const linecounter rcvstride = stoul(argv[k++]); 
  //last time point to read (excluding comment-lines) in the x.xvg files
  const linecounter rcvend = stoul(argv[k++]); 
  //the PMF or any analysis will be done in the ensembles specified by these MDP files
  const vector<string> mdp0prefixes = split<string>(argv[k++]);
  //provide a file which contains seeding dimensionless free energy of 
  //each state in one line seperated by spaces
  //(could be a file that doesn't exist, in which case the program will just seed the free energy to zero)
  const string fseedsstr = string(argv[k++]);
  //if we call minimization routine instead of WHAM iteration
  const bool ifmin = stoul(argv[k++]);
  //if we read the histogram directly instead of binning time-series into histogram
  const bool ifreadhist = stoul(argv[k++]);

  if(is_stdfloat<valtype>::value && ifmin) {
    const string valtypename = typeid(valtype).name();
    throw(General_Exception("Calling minimization routine is likely to result in \
    large fluctuation in free energies, which would be used as exponent of exp() function, \
    and thus numerical unstability and the floating point type '"+valtypename+"' \
    might not have enough precision to withstand this. You should use high precision \
    floating point in order to perform minimization or simply turn off minimization"));
  }
   
  fcout << cmdline;
  fcout << "# ";
  for_each(argv, argv+argc, [] (char* arg) { fcout << "'" << arg << "' "; });
  fcout << endl;

  
  const uint ndim = nbins.size();
  if(ndim != hv.size() || hv.size() != lv.size()) {
    throw(General_Exception("The size of nbins, hv and lv don't match one another"));
  }
  for(uint i = 0; i < hv.size(); ++i) {
    if(hv[i] <= lv[i]) {
      throw(General_Exception("The upper bound of dimension "+tostr(i)+" is no larger than the lower bound"));
    }
  }
  
  if(rcvend  < rcvbegin) {
    throw(General_Exception("rcvend must be larger than rcvbegin"));
  }
   

  //create MDP objects by reading the mdp files
  vector<string> fnprefixes;
  for(uint i = 0; i < nwin; ++i) {
    fnprefixes.emplace_back(sysname+"_"+tostr(i));
  }
  const vector<pMDP> pmdps = CreateMDPs(fnprefixes, mdpsuffixes, suffix2MDP);

  //create ensemble
  multimap<pMDP, uint> pmdp2ipens;
  const auto pens = MDP2Ensemble(pmdps, pmdp2ipens, 1);
  fcout << "#created " << pens.size() << " ensembles\n";
  for(auto& pmdp : pmdps) {
    auto itr = pmdp2ipens.equal_range(pmdp);
    fcout << "#MDP file: " << pmdp->fname << endl;
    pmdp->print();
    fcout << "#contains to the following ensembles: ";
    for(auto& it = itr.first; it != itr.second; ++it) {
      fcout << it->second << " ";
    }
    fcout << endl;
  }

  //create histogram (1 for each ensemble)
  vector<histogram> hists{pens.size(), histogram{ndim, nbins, hv, lv}};
  // total number of valid samples
  vector<histcounter> Nsamples(hists.size(), 0);
  
  if(!ifreadhist) {
    //create the time-series readers
    const bool needpotential = cmpens(pens)&(1<<Ensemble::DTemperature);
    for(const auto& pmdp : pmdps) {
      // if temperature are different among ensembles, we need to read potential energy
      auto vts = pmdp->CreateTimeSeries(needpotential);
  
      fcout << "#require " << vts.size() << " types of time-serie files for " << pmdp->fname << " : ";
      for(const auto& ts : vts) {
        fcout << ts.fnsuffix << " ";
      }
      fcout << endl;
  
      // we overwrite the stride here
      for(auto& ts : vts) {
	const auto& tperline = ts.nst;
	if(rcvbegin != 0 && rcvbegin < tperline) {
          throw(General_Exception("Input arguments rcvbegin ("+tostr(rcvbegin)+") is smaller than the output frequency ("+tostr(tperline)+") in file: "+ts.fnsuffix));
	}
	if(rcvstride < tperline) {
          throw(General_Exception("Input arguments rcvstride ("+tostr(rcvstride)+") is smaller than the output frequency ("+tostr(tperline)+") in file: "+ts.fnsuffix));
	}
	if(rcvend < tperline) {
          throw(General_Exception("Input arguments rcvend ("+tostr(rcvend)+") is smaller than the output frequency ("+tostr(tperline)+") in file: "+ts.fnsuffix));
	}
	//Here I make the assumption that the first datum
	//corresponds to the starting point and should always be skipped
	const linecounter linepert = linecounter(1.0/tperline);
        ts.fio.lb = linecounter(rcvbegin*linepert + 1);
        ts.fio.ls = linecounter(ts.fio.ls * (rcvstride == MAXNLINE ? 1 : rcvstride*linepert));
        ts.fio.le = linecounter(rcvend == MAXNLINE ? rcvend : rcvend*linepert);
      }
  
      const auto tsprefix = getfnfixes(pmdp->fname)[0];
      decltype(vts)::iterator it = vts.begin();
      // For exapnded ensemble, we need to read the state time-series first
      vector<uint> states;
      if(pmdp->isExpandedEnsemble()) {
        // state time-series is always the first one
        vecfromTS(states, *it, nrun, tsprefix); 
        ++it;
      }
  
      vector<valtype> pot;
      if(needpotential) {
        // potential energy time-series 
        vecfromTS(pot, *it, nrun, tsprefix);
        ++it;
      }
  
      //finally we histogram the cv; we first build a histogrammer depends 
      //on whether we are doing expanded ensemble and/or if we need to 
      //deal with potential energy
      histfromTS(it, Nsamples, pmdp, pmdp2ipens, hists, nrun, tsprefix, states, pot); 
    }
  } else {
    cout << "#Reading histograms from files...\n";
    for(uint i = 0; i < pens.size(); ++i) {
      const string fnhist = sysname + "_" + tostr(i) + ".hist";
      Nsamples[i] = readhist(fnhist, hists[i]);
      cout << "#Read " << Nsamples[i] << " sample\n";
    }
  }
  fcout.width(10);
  fcout << "#Number of samples: " << Nsamples << endl;
  for(const auto& h : hists) h.print();
  
  //read the seeding free energy from file
  fileio fio_fseeds(fseedsstr, std::fstream::in, 1, 1, 1, 1);
  vector<valtype> fseeds(0);
  if(fio_fseeds.readaline()) { fseeds = fio_fseeds.line2val(); }

  WHAM<pEnsemble, histogram, narray> wham(hists, pens, Nsamples, tol, fseeds, vector<narray>(0), ifmin);
  
  //calculate the probability density
  //we can have the user specify which ensemble we want to probability to be calculated in
  const vector<pMDP> pmdp0s = CreateMDPs(mdp0prefixes, mdpsuffixes, suffix2MDP);
  multimap<pMDP, uint> pmdp2ipens0;
  auto pE0s = MDP2Ensemble(pmdp0s, pmdp2ipens0, 1);

  for(const auto& pmdp0 : pmdp0s) {
    if(pmdp0->isExpandedEnsemble()) { 
      throw(General_Exception("Can't calculate probability density in a expanded ensemble. Check files in mdp0prefixes"));
    }
    const auto iens = pmdp2ipens0.find(pmdp0)->second;
    const auto pE0 = pE0s[iens];
    auto rho = wham.calrho(rcprt, pE0);
    const auto kB = pmdp0->getkB();
    const auto T = pmdp0->getTs().front();
    //Normalize the probability (just for comparision with other programs)
    printf("#%10s%30s%30s%30s\n","Bin","Vals","PMF","RhoNormalized");
    valtype sum = 0.0;
    narray::iterator itmax = rho.begin();
    for(narray::iterator it = rho.begin(); it != rho.end(); ++it) {
      if(it->second > itmax->second) { itmax = it;}
      sum += it->second;
    }
    for(narray::iterator it = rho.begin(); it != rho.end(); ++it) {
      const coordtype bin = it->first;
      const vector<valtype> val = rho.coord2val(bin);
      const valtype pmf = -kB*T*log(it->second/itmax->second);
      const valtype rhonorm = it->second/sum;
      fcout.width(10);
      fcout << bin;
      fcout.flags(ios::fixed | ios::right);
      fcout.precision(15);
      fcout.width(30);
      fcout << val << pmf;
      fcout.flags(ios::scientific);
      fcout << rhonorm << endl;
    }
  }
}
