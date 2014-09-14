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
#include "typedefs.hpp"
#include "hamiltonian.hpp"
#include "ensemble.hpp"
#include "gnarray.hpp"
#include "fileio.hpp"
#include "fileio_utils.hpp"
#include "gmxmdp.hpp"
#include "ensemble_factory.hpp"
#include "mdp_factory.hpp"

using namespace std;

typedef gnarray<coordtype, valtype, valtype> narray;
typedef gnarray<coordtype, histcounter, valtype> histogram;

void HistogramCV(histogram& hist, linecounter& N, const vector<valtype>& input) {
  if(hist.bin(input) != hist.end()) ++N; 
}

void HistogramExpECV(vector<histogram>& hists, vector<linecounter>& Ns, vector<uint>::const_iterator& itstate, const linecounter nststates, linecounter& count, vector<valtype>::const_iterator& itpot, const vector<valtype>& input) {
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

void HistogramExpCV(vector<histogram>& hists, vector<linecounter>& Ns, vector<uint>::const_iterator& itstate, const linecounter nststates, linecounter& count, const vector<valtype>& input) {
  const auto& state = *itstate;
  ++count;
  if(!(count % nststates)) ++itstate;
  auto& hist = hists[state];
  HistogramCV(hist, Ns[state], input);
}

void HistogramECV(histogram& hist, linecounter& N, vector<valtype>::const_iterator& itpot, const vector<valtype>& input) {
  const auto& pot = *itpot;
  ++itpot;
  vector<valtype> data = input;
  data.push_back(pot);
  HistogramCV(hist, N, data);
}

template < class T, class TS >
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
vector<linecounter> histfromTS(
    vTSit& it, 
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
  vector<linecounter> Nsamples(hists.size(), 0);
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
  return Nsamples;
}

int main(int argc, char* argv[]) {
  string cmdline("#Usage: gwham sysname mdpsuffixes nwin nrun rcprint nbins hv lv tol rcvbegin rcvstride rcvend fseeds\n");
  if (argc != 14) {
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
  const double tol = atof(argv[k++]);
  //first line to read (excluding comment-lines) in the x.xvg files
  const uint rcvbegin = atoi(argv[k++]);
  //Only read every stride lines (excluding comment-lines) in the x.xvg files
  const uint rcvstride = atoi(argv[k++]); 
  //last line to read (excluding comment-lines) in the x.xvg files
  const uint rcvend = atoi(argv[k++]); 
  //provide a file which contains seeding dimensionless free energy of 
  //each state in one line seperated by spaces
  //(could be a file that doesn't exist, in which case the program will just seed the free energy to zero)
  const string fseedsstr = string(argv[k++]);  
  
  cout << cmdline;
  cout << "# ";
  for(int i = 0; i < argc; ++i) {
    cout << "'" << argv[i] << "' ";
  }
  cout << endl;

  
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
  cout << "#created " << pens.size() << " ensembles\n";
  for(auto& pmdp : pmdps) {
    auto itr = pmdp2ipens.equal_range(pmdp);
    cout << "#MDP file: " << pmdp->fname << " contains to the following ensembles: ";
    for(auto& it = itr.first; it != itr.second; ++it) {
      cout << it->second << " ";
    }
    cout << endl;
  }

  //create histogram (1 for each ensemble)
  vector<histogram> hists{pens.size(), histogram{ndim, nbins, hv, lv}};

  //create the time-series readers
  const bool needpotential = cmpens(pens)&(1<<Ensemble::DTemperature);
  for(const auto& pmdp : pmdps) {
    // if temperature are different among ensembles, we need to read potential energy
    auto vts = pmdp->CreateTimeSeries(needpotential);

    cout << "#require " << vts.size() << " types of time-serie files for " << pmdp->fname << " : ";
    for(const auto& ts : vts) {
      cout << ts.fnsuffix << " ";
    }
    cout << endl;

    // we overwrite the stride here
    for(auto& ts : vts) {
      ts.fio.lb = rcvbegin;
      ts.fio.ls *= rcvstride;
      ts.fio.le = rcvend;
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
    const auto Nsamples = histfromTS(it, pmdp, pmdp2ipens, hists, nrun, tsprefix, states, pot); 
  }
  for(const auto& h : hists) h.print();
  
}
