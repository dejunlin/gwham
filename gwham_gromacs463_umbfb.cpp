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
  if (argc != 15) {
    cerr << "Usage: gwham_gromacs463_umbfb sysname nrun nwin runmin rcsample rcprint Temperature nbins hv lv tol xvgskip xvgstride ncolskip\n";
    exit(-1);
  }
  const string sysname = string(argv[1]); //name of the system 
  				 //(your files must be named as "sysname_RUNIDx_WINID.xvg", "sysname_RUNIDdhdl_WINID.xvg")
  const uint nrun = atoi(argv[2]); //# of runs performed
  const uint nwin = atoi(argv[3]); //# of replica
  const uint runmin = atoi(argv[4]); //first replica id
  const string rcsmpstr = string(argv[5]); //which pull-group/RC you want to be histogramed on, where the pull-group id starts from 0 (1st pull-group not the reference group)
  const string rcprtstr = string(argv[6]); //which pull-group/RC you want the PMF on. 
  					   //Since the mpd2rst and xxvg2hist classes will re-order the RCs so that 
					   //they are the 1st N dimension of the histogram (N is the number of RCs you histogramed on), 
					   //rcrpt must be in the range of [0,N-1]. We could potentially allow this range to extend over N-1 but
					   //that will require we know if all the restrained RCs are the RCs we're interested in
  const valtype T = atof(argv[7]); //what temperature to use to calculate PMF
  const string nbinstr = string(argv[8]); //# of bins in each dimension
  const string hvstr = string(argv[9]); //upper bounds in each dimension
  const string lvstr = string(argv[10]); //lower bounds in each dimension
  const double tol = atof(argv[11]); //tolerance for WHAM iteration
  const uint xvgskip = atoi(argv[12]); //Skip this many lines from the beginning of the x.xvg files
  const uint xvgstride = atoi(argv[13]); //Only read every stride lines (excluding comment-lines) in the x.xvg files
  const uint ncolskip = atoi(argv[14]); //skip the 1st ncolskip elements of any line of a x.xvg file
  
  cout << "# gwham_gromacs463_umbfb sysname nrun nwin runmin rcsample rcprint Temperature nbins hv lv tol xvgskip xvgstride ncolskip\n";
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

  //Check if rcprt makes sence
  for(uint i = 0; i < rcprt.size(); ++i) {
    if(rcprt[i] >= rcsmp.size()) {
      cerr << "All elements of rcprint must be in the range of [0," << rcsmp.size()-1 << "]\n";
      exit(-1);
    }
  }

  //Turn rcsmp into a bitmask
  bitset<MAXNRST> rcsmpmask;
  for(uint i = 0; i < rcsmp.size(); ++i) {
    rcsmpmask |= (bitunit << rcsmp[i]);
  }
  //Bitmask for all the pull group that're restrained
  bitset<MAXNRST> rstmask;
  
  const uint ndim = nbins.size();
  //histograms for each window
  vector<histogram> hists(nwin,histogram(ndim,nbins,hv,lv));
  cout << "#Histogram is of dimension " << ndim << endl;
  
  //Read the parameters from mdp file, construct the restraint functor and bin the x.xvg file 
  //into histograms
  const mdp2pullpot mdp2rst(rcsmpmask);
  fileio<mdp2pullpot, map<uint,vector<umbrella_fb*> >, vector<Hamiltonian<RST_fbXLAMBDAsgl>* > > fmdp(mdp2rst, std::fstream::in, 0);
  vector<Hamiltonian<RST_fbXLAMBDAsgl>* > V;
  map<uint, vector<umbrella_fb*> > rstfuncts;
  for(uint i = 0; i < nwin; ++i) {
    char winid[MAXNDIGWIN];
    sprintf(winid,"%d",i);
    string mdpfname = sysname + "_0_" + winid + ".mdp";
    fmdp(mdpfname,rstfuncts,V);
  }

  //Check if all windows have the same temperature, if not, just quit
  for(uint i = 0; i < V.size(); ++i) {
    const RST_fbXLAMBDAsgl ens = V[i]->getens();
    const valtype ensT = (ens.getparams())[1];
    if( ensT != T ) {
      cerr << "The " << i+1 << "'th ensemble is in temperature " << ensT << " other than the one specified in command line: " << T << endl;
      exit(-1);
    }
  }

  //number of data points for each window
  vector<uint> N(nwin,0);
  //Read x.xvg file for each run
  for(uint i = 0; i < nwin; ++i) {
    char winid[MAXNDIGWIN];
    sprintf(winid,"%d",i);
    for(uint j = runmin; j < nrun; ++j) {
      char runid[MAXNDIGRUN];
      sprintf(runid,"%d",j);
      string xxvgfname = sysname + "_" + runid + "x_" + winid + ".xvg";
      if(j == 0) {
        const xxvg2hist<histogram,umbrella_fb> xvg2hist(ncolskip, rcsmpmask,rstfuncts,xvgskip,xvgstride);
        fileio<xxvg2hist<histogram,umbrella_fb>,histogram> fxxvg(xvg2hist, std::fstream::in);
        N[i] += fxxvg(xxvgfname,hists[i]);
      } else {
        const xxvg2hist<histogram,umbrella_fb> xvg2hist(ncolskip, rcsmpmask,rstfuncts,0,xvgstride);
        fileio<xxvg2hist<histogram,umbrella_fb>,histogram> fxxvg(xvg2hist, std::fstream::in);
        N[i] += fxxvg(xxvgfname,hists[i]);
      }
    }
  }
  /*cout << "#number of samples in each window: ";
  copy(N.begin(),N.end(),ostream_iterator<uint>(cout," ")); cout << endl;
  //print the histograms
  for(uint i = 0; i < nwin; ++i) {
    cout << "#Histogram " << i << endl;
    hists[i].print();
  }*/

  //build records of what bins in the histograms are non-zero
  map<coordtype, vector<uint> > record;
  for(uint i = 0; i < nwin; ++i) {
    typename histogram::iterator it;
    for(it = hists[i].begin(); it != hists[i].end(); ++it) {
      const coordtype coord = it->first;
      record[coord].push_back(i);
    }
  }
  /*for(map<coordtype,vector<uint> >::const_iterator it = record.begin(); it != record.end(); ++it) {
    const coordtype coord = it->first;
    const vector<uint> histid = it->second;
    const vector<valtype> vals = hists[histid[0]].coord2val(coord);
    cout << "#Bin: ";
    copy(coord.begin(),coord.end(),ostream_iterator<uint>(cout," "));
    copy(vals.begin(),vals.end(),ostream_iterator<valtype>(cout," "));
    cout << " is contributed from histogram: ";
    copy(histid.begin(),histid.end(),ostream_iterator<uint>(cout," ")); cout << endl;
  }*/
  
  cout << "# DBL_MAX supported on this machine is " << DBL_MAX << endl;
  cout << "# DBL_MIN supported on this machine is " << DBL_MIN << endl;
  cout << "# NOTE: WHAM iteration might break if the arguments of a exp() evaluation is ";
  cout << "out of the range [" << MINEXPARG << "," << MAXEXPARG << "]" << endl;

  WHAM<RST_fbXLAMBDAsgl, histogram, narray> wham(record,hists, V, N, tol);
  vector<valtype> params;
  params.push_back(BoltzmannkJ);
  params.push_back(T);
  params.insert(params.end(),rcsmpmask.count()*5,0.0);
  params.push_back(0.0); //LAMBDAsgl::L
  params.push_back(0.0); //LAMBDAsgl::i
  Hamiltonian<RST_fbXLAMBDAsgl> V0(params);
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
