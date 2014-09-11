#include "typedefs.hpp"
#include "hamiltonian.hpp"
#include "ensemble.hpp"
#include "gnarray.hpp"
#include "fileio.hpp"
#include "fileio_utils.hpp"
#include <limits>
#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <iterator>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

typedef gnarray<coordtype, valtype, valtype> narray;
typedef gnarray<coordtype, histcounter, valtype> histogram;

int main(int argc, char* argv[]) {
  string cmdline("#Usage: gwham sysname nwin nrun rcprint nbins hv lv tol rcvbegin rcvstride rcvend fseeds\n");
  if (argc != 13) {
    cerr << cmdline;
    exit(-1);
  }
  uint k = 1;
  //name of the system
  //(your files should be named as "sysname_WINID_RUNID*.*")
  const string sysname = string(argv[k++]); 
  //file names suffixes for each type of input files ( the number of file types depends on the system mdp setup) 
  const vector<string> fnsuffixes = split<string>(argv[k++]); 
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

  const string mdpsuffix = fnsuffixes[0];
  fnsuffixes.erase(fnsuffixes.begin());
  vector<pMDP> pmdps;
  for(uint i = 0; i < nwin; ++i) {
    //first suffix is for MDP file (format is sysname_winid.suffix
    const string fmdp = sysname+"_"+tostr(i)+"."+mdpsuffix;

  }
}
