/*
 * =====================================================================================
 *
 *       Filename:  gwham.cpp
 *
 *    Description:  Generic driver program for generalized WHAM for umbrella sampling
 *
 *        Version:  1.0
 *        Created:  06/02/2014 12:33:18 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dejun Lin 
 *   Organization:  
 *
 * =====================================================================================
 */

#include "typedefs.hpp"
#include "hamiltonian.hpp"
#include "ensemble.hpp"
#include "gwham.hpp"
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
  string cmdline("#Usage: gwham metadata ndim Temperature nbins hv lv tol tsstart tsstride tsend fseeds\n");
  if (argc != 12) {
    cerr << cmdline;
    exit(-1);
  }
  const string datastr = string(argv[1]); //metadata file
  const uint ndim = atoi(argv[2]);  //number of dimension of the histogram
  const valtype T = atof(argv[3]); //temperature
  const string nbinstr = string(argv[4]); //# of bins in each dimension
  const string hvstr = string(argv[5]); //upper bounds in each dimension
  const string lvstr = string(argv[6]); //lower bounds in each dimension
  const double tol = atof(argv[7]); //tolerance for WHAM iteration
  const uint tsstart = atoi(argv[8]); //first line to read (excluding comment-lines) in the x.ts files
  const uint tsstride = atoi(argv[9]); //Only read every stride lines (excluding comment-lines) in the x.ts files
  const uint tsend = atoi(argv[10]); //last line to read (excluding comment-lines) in the x.ts files
  const string fseedsstr = string(argv[11]); //provide a file which contains seeding dimensionless free energy of 
                                             //each state in one line seperated by spaces
  					     //(could be a file that doesn't exist, in which case the program will just seed the free energy to zero) 

  cout << cmdline;
  cout << "# ";
  for(int i = 0; i < argc; ++i) {
    cout << "'" << argv[i] << "' ";
  }
  cout << endl;

  vector<uint> nbins;
  vector<valtype> hv, lv;
  split<uint>(nbins,nbinstr);
  split<valtype>(hv,hvstr);
  split<valtype>(lv,lvstr);
  //Check if the histogram parameters (nbins, hv, lv) make sense
  if(nbins.size() != hv.size() || hv.size() != lv.size()) {
    cerr << "The size of nbins, hv and lv don't match one another" << endl;
    exit(-1);
  }
  for(uint i = 0; i < hv.size(); ++i) {
    if(hv[i] <= lv[i]) {
      cerr << "The upper bound of dimension " << i << " is no larger than the lower bound\n";
      cerr << "hv = " << hv[i] << endl;
      cerr << "lv = " << lv[i] << endl;
      exit(-1);
    }
  }
  
  //Read the metadata file
  //Each line should be tsfile 'rst center' 'rst force', where center and force must be of the same dimension as ndim
  vector<Hamiltonian<RSTXLAMBDAsgl>* > V;
  vector<histogram> hists;
  vector<uint> N;
  fileio fmeta(datastr, std::fstream::in, false, 1, 1, MAXNLINE);
  while(fmeta.readaline()) {
    vector<string> line = fmeta.line2str();
    if(line.size() != 1 + ndim*2) {
      cerr << "Line: ";
      copy(line.begin(), line.end(), ostream_iterator<string>(cout, " "));
      cerr << "\n has wrong format; NOTE, the format should be:\n";
      cerr << "tsfile 'rst center' 'rst force', where center and force must be of the same dimension as ndim\n";
      exit(-1);
    }

    //Set up the hamiltonian parameters
    vector<valtype> rstcenters, rstforces;
    strto<valtype>(rstcenters, vector<string>(line.begin()+1, line.begin()+1+ndim));
    strto<valtype>(rstforces, vector<string>(line.begin()+1+ndim, line.begin()+1+2*ndim));
    vector<valtype> params;
    params.push_back(Boltzmannkcal);
    params.push_back(T);
    params.insert(params.end(), rstforces.begin(), rstforces.end());
    params.insert(params.end(), rstcenters.begin(), rstcenters.end());
    params.push_back(0.0);
    params.push_back(0.0);
    V.push_back(new Hamiltonian<RSTXLAMBDAsgl>(params));
    
    //Read the ts file and generate histogram
    hists.push_back(histogram(ndim,nbins,hv,lv));
    const string tsstr = line[0];
    fileio fts(tsstr, std::fstream::in, false, tsstart, tsstride, tsend);
    const uint lastwin = hists.size()-1;
    histogram& hist = hists[lastwin];
    uint nsamples = 0;
    while(fts.readaline()) {
      //Format of ts file should be: Time 'RC' where RC should be of the same dimension as ndim
      vector<valtype> line = fts.line2val();
      if(hist.bin(vector<valtype>(line.begin()+1, line.end())) != hist.end()) {
	++nsamples;
      }
    }
    N.push_back(nsamples);
  }

  const uint nwin = hists.size();
  cout << "#A total of " << nwin << " ts windows processed\n";
  cout << "#Number of samples in each window: ";
  copy(N.begin(),N.end(),ostream_iterator<uint>(cout," ")); cout << endl;
  //build records of what bins in the histograms are non-zero
  map<coordtype, vector<uint> > record;
  for(uint i = 0; i < nwin; ++i) {
    typename histogram::iterator it;
    for(it = hists[i].begin(); it != hists[i].end(); ++it) {
      const coordtype coord = it->first;
      record[coord].push_back(i);
    }
  }
  cout << "# DBL_MAX supported on this machine is " << DBL_MAX << endl;
  cout << "# DBL_MIN supported on this machine is " << DBL_MIN << endl;
  cout << "# NOTE: WHAM iteration might break if the arguments of a exp() evaluation is ";
  cout << "out of the range [" << MINEXPARG << "," << MAXEXPARG << "]" << endl;
  //read the seeding free energy from file
  fileio fio_fseeds(fseedsstr, std::fstream::in, 1, 1, 1, 1);
  vector<valtype> fseeds(0);
  if(fio_fseeds.readaline()) {
    fseeds = fio_fseeds.line2val();
  }
  //do wham
  WHAM<RSTXLAMBDAsgl, histogram, narray> wham(record,hists, V, N, tol, fseeds);
  //calculate PMF
  vector<valtype> params;
  params.push_back(Boltzmannkcal);
  params.push_back(T);
  params.insert(params.end(), ndim*2, 0.0);
  params.push_back(0.0); //LAMBDAsgl::L
  params.push_back(0.0); //LAMBDAsgl::i
  Hamiltonian<RSTXLAMBDAsgl> V0(params);
  vector<uint> rcprt;
  for(uint i = 0; i < ndim; ++i) { rcprt.push_back(i); }
  narray rho = wham.calrho(rcprt,V0);
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
    const valtype pmf = -Boltzmannkcal*T*log(it->second/itmax->second);
    const valtype rhonorm = it->second/sum;
    for(uint i = 0; i < bin.size(); ++i) { printf("%10d",bin[i]);}
    for(uint i = 0; i < val.size(); ++i) { printf("%30.15lf",val[i]);}
    printf("%30.15lf%30.15lf\n",pmf,rhonorm);
  }
  //calculate the inconsistency between the consensus histogram and the raw historam
  const vector<valtype> eitas = wham.whamvsraw(rho);
  printf("#WHAM inconsistency test\n");
  printf("#%10s%30s%30s\n", "Window", "Center", "Eita");
  for(uint i = 0; i < eitas.size(); ++i) {
    const vector<valtype> rstparams = V[i]->getens()->getRSTparams();
    const vector<valtype> rstcenter(rstparams.end()-ndim, rstparams.end());
    printf("%10d",i);
    for(uint k = 0; k < ndim; ++k) { printf("%30lf", rstcenter[k]); }
    printf("%30lf\n", eitas[i]);
  }
}
