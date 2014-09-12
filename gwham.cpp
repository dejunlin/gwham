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

void HistogramCV(histogram& hist, uint& N, const vector<valtype>& input) {
  if(hist.bin(input) != hist.end()) ++N; 
}

void HistogramCV(vector<histogram>& hists, uint& N, vector<uint>::const_iterator& itstate, vector<valtype>::const_iterator itpot, const vector<valtype>& input) {
  const auto& state = *itstate;
  auto& hist = hists[state];
  const auto& pot = *itpot;
  vector<valtype> data = input;
  data.push_back(pot);
  HistogramCV(hist, N, data);
  ++itstate;
  ++itpot;
}

void HistogramCV(vector<histogram>& hists, uint& N, vector<uint>::const_iterator& itstate, const vector<valtype>& input) {
  const auto& state = *itstate;
  auto& hist = hists[state];
  HistogramCV(hist, N, input);
  ++itstate;
}

void HistogramCV(histogram& hist, uint& N, vector<valtype>::const_iterator itpot, const vector<valtype>& input) {
  const auto& pot = *itpot;
  vector<valtype> data = input;
  data.push_back(pot);
  HistogramCV(hist, N, data);
  ++itpot;
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
  for(const auto& pmdp : pmdps) {
    const bool needpotential = cmpens(pens)&(1<<Ensemble::DTemperature);
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
    auto itstate = states.begin(); 
    if(pmdp->isExpandedEnsemble()) {
      // state time-series is always the first one
      auto& ts = *it;
      const auto tssuffix = ts.fnsuffix;
      for(uint i = 0; i < nrun; ++i) {
	const string fts = tsprefix + "_" + tostr(i) + tssuffix;
	ts(fts, [&states] (const vector<valtype>& input) 
	        {
		  states.emplace_back(input.front());
		}
          );
      }
      ++it;
    }

    vector<valtype> pot;
    auto itpot = pot.begin();
    if(needpotential) {
      // potential energy time-series 
      auto& ts = *it;
      const auto tssuffix = ts.fnsuffix;
      for(uint i = 0; i < nrun; ++i) {
	const string fts = tsprefix + "_" + tostr(i) + tssuffix;
	ts(fts, [&pot] (const vector<valtype>& input) 
	        {
		  pot.emplace_back(input.front());
		}
          );
      }
      ++it;
    }
    
    //finally we histogram the cv; we first build a histogrammer depends 
    //on whether we are doing expanded ensemble and/or if we need to 
    //deal with potential energy
    function<void(const vector<valtype>& input)> cvhistogrammer;
    using namespace std::placeholders;
    uint Nsamples = 0;
    if(states.size()) {
      if(pot.size()) {
	cvhistogrammer = std::bind(static_cast<void(*)(vector<histogram>&, uint&, vector<uint>::const_iterator&, vector<valtype>::const_iterator, const vector<valtype>&)>(&HistogramCV), ref(hists), ref(Nsamples), ref(itstate), ref(itpot), _1); 
      } else {
	cvhistogrammer = std::bind(&HistogramCV, ref(hists), ref(Nsamples), ref(itstate), _1); 
      }
    } else {
      const uint ihist = pmdp2ipens.find(pmdp)->second;
      if(pot.size()) {
	cvhistogrammer = std::bind(&HistogramCV, ref(hists[ihist]), ref(Nsamples), ref(itpot), _1); 
      } else {
	cvhistogrammer = std::bind(&HistogramCV, ref(hists[ihist]), ref(Nsamples), _1); 
      }
    }

    auto& ts = *it;
    const auto tssuffix = ts.fnsuffix;
    for(uint i = 0; i < nrun; ++i) {
      const string fts = tsprefix + "_" + tostr(i) + tssuffix;
      ts(fts, cvhistogrammer); 
    }
  }
}
