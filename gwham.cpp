#include <limits>
#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <iterator>
#include <stdio.h>
#include <stdlib.h>
#include "typedefs.hpp"
#include "hamiltonian.hpp"
#include "ensemble.hpp"
#include "gnarray.hpp"
#include "fileio.hpp"
#include "fileio_utils.hpp"
#include "gmxmdp.hpp"
#include "ensemble_factory.hpp"

using namespace std;

typedef gnarray<coordtype, valtype, valtype> narray;
typedef gnarray<coordtype, histcounter, valtype> histogram;

typedef map<string, function<pMDP(const string&)>> MDP_Factory;
static const MDP_Factory suffix2mdp {
  {"mdp", GMXMDP::CreateMDP},
};


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
   
  for(const auto& mdpsuffix : mdpsuffixes) {
    const MDP_Factory::const_iterator it = suffix2mdp.find(mdpsuffix);
    if(it == suffix2mdp.end()) {
      throw(General_Exception("MDP file with suffix '"+mdpsuffix+"' not supported"));
    }
  }

  //create MDP objects by reading the mdp files
  //first suffix is for MDP file (format is sysname_winid.suffix)
  vector<pMDP> pmdps;
  for(uint i = 0; i < nwin; ++i) {
    for(const auto& mdpsuffix : mdpsuffixes) {
      const string fmdp = sysname+"_"+tostr(i)+"."+mdpsuffix;
      try {
	// check if the file exists
        fileio fio(fmdp, fstream::in, 1, 0, 1, MAXNLINE, ";#@");
      } catch (FILEIO_Exception& fioex) {
	continue;
      }
      pmdps.emplace_back( suffix2mdp.find(mdpsuffix)->second(fmdp) );
    }
  }

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
    // if temperature are different among ensembles, we need to read potential energy
    auto vts = pmdp->CreateTimeSeries(chkens(pens)&(1<<Ensemble::DTemperature));

    cout << "#require " << vts.size() << " types of time-serie files for " << pmdp->fname << " : ";
    for(const auto& ts : vts) {
      cout << ts.fnsuffix << " ";
    }
    cout << endl;

    const auto tsprefix = getfnfixes(pmdp->fname)[0];
    for(auto& ts : vts) {
      const auto tssuffix = ts.fnsuffix;
      for(uint i = 0; i < nrun; ++i) {
	const string fts = tsprefix + "_" + tostr(i) + tssuffix;
        cout << "about to read file: " << fts << endl;
      }
    }

  }

}
