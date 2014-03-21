#include "typedefs.hpp"
#include "hamiltonian.hpp"
#include "ensemble.hpp"
#include "gwham.hpp"
#include "gnarray.hpp"
#include "fileio.hpp"
#include "fileio_utils.hpp"
#include "gmxxvgio.hpp"
#include "gmxmdpio.hpp"
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
  string cmdline("Usage: gwham_gromacs463_umb sysname nrun nwin rcsample rcprint Temperature nbins hv lv tol xvgstart xvgstride xvgend ncolskip\n");
  if (argc != 15) {
    cerr << cmdline;
    exit(-1);
  }
  const string sysname = string(argv[1]); //name of the system 
  				 //(your files must be named as "sysname_RUNIDx_WINID.xvg", "sysname_RUNIDdhdl_WINID.xvg")
  const uint nrun = atoi(argv[2]); //# of runs performed
  const uint nwin = atoi(argv[3]); //# of replica
  const string rcsmpstr = string(argv[4]); //which pull-group/RC you want to be histogramed on, where the pull-group id starts from 0 (1st pull-group not the reference group)
  const string rcprtstr = string(argv[5]); //which pull-group/RC you want the PMF on. 
  					   //Since the mpd2rst and xxvg2hist classes will re-order the RCs so that 
					   //they are the 1st N dimension of the histogram (N is the number of RCs you histogramed on), 
					   //rcrpt must be in the range of [0,N-1]. We could potentially allow this range to extend over N-1 but
					   //that will require we know if all the restrained RCs are the RCs we're interested in
  const valtype T = atof(argv[6]); //what temperature to use to calculate PMF
  const string nbinstr = string(argv[7]); //# of bins in each dimension
  const string hvstr = string(argv[8]); //upper bounds in each dimension
  const string lvstr = string(argv[9]); //lower bounds in each dimension
  const double tol = atof(argv[10]); //tolerance for WHAM iteration
  const uint xvgstart = atoi(argv[11]); //first line to read (excluding comment-lines) in the x.xvg files
  const uint xvgstride = atoi(argv[12]); //Only read every stride lines (excluding comment-lines) in the x.xvg files
  const uint xvgend = atoi(argv[13]); //last line to read (excluding comment-lines) in the x.xvg files
  const uint ncolskip = atoi(argv[14]); //skip the 1st ncolskip elements in any line of a x.xvg file
  
  cout << cmdline;
  cout << "# ";
  copy(argv,argv+argc,ostream_iterator<char*>(cout," "));
  cout << endl;

  vector<uint> rcsmp, rcprt;
  parser<uint>(rcsmp,rcsmpstr);
  parser<uint>(rcprt,rcprtstr);
  vector<uint> nbins;
  vector<valtype> hv, lv;
  parser<uint>(nbins,nbinstr);
  parser<valtype>(hv,hvstr);
  parser<valtype>(lv,lvstr);
  
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

  //Turn rc into a bitmask
  bitset<MAXNRST> rcsmpmask;
  for(uint i = 0; i < rcsmp.size(); ++i) {
    rcsmpmask |= (bitunit << rcsmp[i]);
  }
  //Bitmask for all the pull group that're restrained
  bitset<MAXNRST> rstmask;
  
  const uint ndim = nbins.size();
  //histograms for each window
  vector<histogram> hists(nwin,histogram(ndim,nbins,hv,lv));
  
  //Read the parameters from mdp file, construct the restraint functor and bin the x.xvg file 
  //into histograms
  fileio fmdp(std::fstream::in, false, 1, 1, MAXNLINE);
  mdp2pullpot mdp2rst(rcsmpmask, fmdp);
  vector<Hamiltonian<RSTXLAMBDAsgl>* > V;
  map<uint, vector<umbrella*> > rstfuncts;
  for(uint i = 0; i < nwin; ++i) {
    char winid[MAXNDIGWIN];
    sprintf(winid,"%d",i);
    string mdpfname = sysname + "_0_" + winid + ".mdp";
    mdp2rst(mdpfname, rstfuncts, V);
  }
  
  //Check if all windows have the same temperature, if not, just quit
  for(uint i = 0; i < V.size(); ++i) {
    const RSTXLAMBDAsgl ens = V[i]->getens();
    const valtype ensT = (ens.getparams())[1];
    if( ensT != T ) {
      cerr << "The " << i+1 << "'th ensemble is in temperature " << ensT << " other than the one specified in command line: " << T << endl;
      exit(-1);
    }
  }

  //number of data points for each window
  vector<uint> N(nwin,0);
  fileio fxxvg(std::fstream::in, 1, xvgstart, xvgstride, xvgend);
  xxvg2hist<histogram,umbrella> xvg2hist(ncolskip, rcsmpmask, rstfuncts, fxxvg);
  //Read x.xvg file for each run
  for(uint i = 0; i < nwin; ++i) {
    char winid[MAXNDIGWIN];
    sprintf(winid,"%d",i);
    for(uint j = 0; j < nrun; ++j) {
      char runid[MAXNDIGRUN];
      sprintf(runid,"%d",j);
      string xxvgfname = sysname + "_" + runid + "x_" + winid + ".xvg";
      N[i] += xvg2hist(xxvgfname,hists[i]);
    }
  }
  cout << "#number of samples in each window: ";
  copy(N.begin(),N.end(),ostream_iterator<uint>(cout," ")); cout << endl;
  //print the histograms
  for(uint i = 0; i < nwin; ++i) {
    cout << "#Histogram " << i << endl;
    hists[i].print();
  }

  //build records of what bins in the histograms are non-zero
  map<coordtype, vector<uint> > record;
  for(uint i = 0; i < nwin; ++i) {
    typename histogram::iterator it;
    for(it = hists[i].begin(); it != hists[i].end(); ++it) {
      const coordtype coord = it->first;
      record[coord].push_back(i);
    }
  }
  for(map<coordtype,vector<uint> >::const_iterator it = record.begin(); it != record.end(); ++it) {
    const coordtype coord = it->first;
    const vector<uint> histid = it->second;
    const vector<valtype> vals = hists[histid[0]].coord2val(coord);
    cout << "#Bin: ";
    copy(coord.begin(),coord.end(),ostream_iterator<uint>(cout," "));
    copy(vals.begin(),vals.end(),ostream_iterator<valtype>(cout," "));
    cout << " is contributed from " << histid.size() << " histograms: ";
    copy(histid.begin(),histid.end(),ostream_iterator<uint>(cout," ")); cout << endl;
  }
  
  cout << "# DBL_MAX supported on this machine is " << DBL_MAX << endl;
  cout << "# DBL_MIN supported on this machine is " << DBL_MIN << endl;
  cout << "# NOTE: WHAM iteration might break if the arguments of a exp() evaluation is ";
  cout << "out of the range [" << MINEXPARG << "," << MAXEXPARG << "]" << endl;

  WHAM<RSTXLAMBDAsgl, histogram, narray> wham(record,hists, V, N, tol);
  vector<valtype> params;
  params.push_back(BoltzmannkJ);
  params.push_back(T);
  params.insert(params.end(),rcsmpmask.count()*2,0.0);
  params.push_back(0.0); //LAMBDAsgl::L
  params.push_back(0.0); //LAMBDAsgl::i
  Hamiltonian<RSTXLAMBDAsgl> V0(params);
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
    const valtype pmf = -BoltzmannkJ*T*log(it->second/itmax->second);
    const valtype rhonorm = it->second/sum;
    for(uint i = 0; i < bin.size(); ++i) { printf("%10d",bin[i]);}
    for(uint i = 0; i < val.size(); ++i) { printf("%30.15lf",val[i]);}
    printf("%30.15lf%30.15lf\n",pmf,rhonorm);
  }

}
