#include "typedefs.hpp"
#include "ensemble.hpp"
#include "gwham.hpp"
#include "gnarray.hpp"
#include "fileio.hpp"
#include "gmxxvgio.hpp"
#include "gmxmdpio.hpp"
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
  if (argc != 12) {
    cerr << "Usage: gwham_gromacs463 sysname nrun nrepl nstate iftemperature ifpressure rc nbins hv lv tol\n";
    exit(-1);
  }
  const string sysname = string(argv[1]); //name of the system 
  				 //(your files must be named as "sysname_RUNIDx_REPLID.xvg", "sysname_RUNIDdhdl_REPLID.xvg")
  const uint nrun = atoi(argv[2]); //# of runs performed
  const uint nrepl = atoi(argv[3]); //# of replica
  const uint nstate = atoi(argv[4]); //# of states
  const bool iftemperature = atoi(argv[5]); //if temperatures are different in different replica
  const bool ifpressure = atoi(argv[6]); //if pressures are different in different replica
  const string rcstr = string(argv[7]); //which pull-group/RC you want the PMF on, where the pull-group id starts from 0 (1st pull-group not the reference group)
  const string nbinstr = string(argv[8]); //# of bins in each dimension
  const string hvstr = string(argv[9]); //upper bounds in each dimension
  const string lvstr = string(argv[10]); //lower bounds in each dimension
  const double tol = atof(argv[11]); //tolerance for WHAM iteration
  
  cout << "#gwham_gromacs463 sysname nrun nrepl nstate iftemperature ifpressure rc nbins hv lv tol\n";
  cout << "#";
  copy(argv,argv+argc,ostream_iterator<char*>(cout," "));
  cout << endl;

  vector<uint> rc;
  parser<uint>(rc,rcstr);
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
  bitset<MAXNRST> rcmask;
  for(uint i = 0; i < rc.size(); ++i) {
    rcmask |= (bitunit << rc[i]);
  }
  
  const uint ndim = nbins.size();
  //histograms for each replica
  vector<histogram> hists(nrepl,histogram(ndim,nbins,hv,lv));
  //number of data points for each replica
  vector<uint> N(nrepl,0);
  
  //Read the parameters from mdp file, construct the restraint functor and bin the x.xvg file 
  //into histograms
  const mdp2pullpot mdp2rst;
  fileio<mdp2pullpot,map<uint,umbrella_fb> > fmdp(mdp2rst, std::fstream::in);
  for(uint i = 0; i < nrepl; ++i) {
    char replid[MAXNDIGWIN];
    sprintf(replid,"%d",i);
    string mdpfname = sysname + "_0_" + replid + ".mdp";
    map<uint,umbrella_fb> rstfunct;
    fmdp(rstfunct,mdpfname);
    //Read x.xvg file for each run
    const xxvg2hist<histogram,umbrella_fb> xvg2hist(rcmask,rstfunct);
    fileio<xxvg2hist<histogram,umbrella_fb>,histogram> fxxvg(xvg2hist, std::fstream::in);
    for(uint j = 0; j < nrun; ++j) {
      char runid[MAXNDIGRUN];
      sprintf(runid,"%d",j);
      string xxvgfname = sysname + "_" + runid + "x_" + replid + ".xvg";
      N[i] += fxxvg(hists[i],xxvgfname);
    }
  }
  cout << "#number of samples in each window: ";
  copy(N.begin(),N.end(),ostream_iterator<uint>(cout," ")); cout << endl;
  /*//print the histograms
  for(uint i = 0; i < nrepl; ++i) {
    cout << "#Histogram " << i << endl;
    hists[i].print();
  }*/

  //build records of what bins in the histograms are non-zero
  map<coordtype, vector<uint> > record;
  for(uint i = 0; i < nrepl; ++i) {
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

  //Construct WHAM
  const valtype kB = 0.0083144621; //kB in kJ/mol/K
  const valtype T = 300; //temperature in kelvin
  vector<Hamiltonian<LAMBDAsgl>* > V;
  for(uint i = 0; i < nrepl; ++i) {
    vector<valtype> params(4,0.0);
    params[0] = kB; 
    params[1] = T; //Temperature
    params[2] = 1;
    params[3] = i;
    V.push_back(new Hamiltonian<LAMBDAsgl>(params));
  }

  cout << "# DBL_MAX supported on this machine is " << DBL_MAX << endl;
  cout << "# DBL_MIN supported on this machine is " << DBL_MIN << endl;
  cout << "# NOTE: WHAM iteration might break if the arguments of a exp() evaluation is ";
  cout << "out of the range [" << MINEXPARG << "," << MAXEXPARG << "]" << endl;

  WHAM<LAMBDAsgl, histogram, narray> wham(record,hists, V, N, tol);
  vector<valtype> params(4,0.0);
  params[0] = kB; 
  params[1] = T; //Temperature
  params[2] = 0;
  params[3] = 0;
  Hamiltonian<LAMBDAsgl> V0(params);
  narray rho = wham.calrho(rc,V0);
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
    for(uint i = 0; i < bin.size(); ++i) { printf("%10d",bin[i]);}
    for(uint i = 0; i < val.size(); ++i) { printf("%30.15lf",val[i]);}
    printf("%30.15lf%30.15lf\n",pmf,rhonorm);
  }
}
