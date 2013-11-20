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
  if (argc != 9) {
    cerr << "Usage: gwham_gromacs463_umb sysname nrun nwin rc nbins hv lv tol\n";
    exit(-1);
  }
  const string sysname = string(argv[1]); //name of the system 
  				 //(your files must be named as "sysname_RUNIDx_WINID.xvg", "sysname_RUNIDdhdl_WINID.xvg")
  const uint nrun = atoi(argv[2]); //# of runs performed
  const uint nwin = atoi(argv[3]); //# of replica
  const string rcstr = string(argv[4]); //which pull-group/RC you want the PMF on, where the pull-group id starts from 0 (1st pull-group not the reference group) 
  const string nbinstr = string(argv[5]); //# of bins in each dimension
  const string hvstr = string(argv[6]); //upper bounds in each dimension
  const string lvstr = string(argv[7]); //lower bounds in each dimension
  const double tol = atof(argv[8]); //tolerance for WHAM iteration
  
  cout << "# gwham_gromacs463_umb sysname nrun nwin rc nbins hv lv tol\n";
  cout << "# ";
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
  //Bitmask for all the pull group that're restrained
  bitset<MAXNRST> rstmask;
  //number of pullgroups in x.xvg file, including those with virtually no restraint on them
  uint npullgrps = 0;
  
  const uint ndim = nbins.size();
  //histograms for each window
  vector<histogram> hists(nwin,histogram(ndim,nbins,hv,lv));
  //number of data points for each replica
  vector<uint> N(nwin,0);
  
  //Read the parameters from mdp file, construct the restraint functor and bin the x.xvg file 
  //into histograms
  const mdp2pullpot mdp2rst;
  fileio<mdp2pullpot,map<uint,umbrella> > fmdp(mdp2rst, std::fstream::in);
  vector<map<uint,umbrella> > rstfuncts;
  for(uint i = 0; i < nwin; ++i) {
    char winid[MAXNDIGWIN];
    sprintf(winid,"%d",i);
    string mdpfname = sysname + "_0_" + winid + ".mdp";
    map<uint,umbrella> rstfunct;
    fmdp(rstfunct,mdpfname);
    rstfuncts.push_back(rstfunct);
    //Read x.xvg file for each run
    const xxvg2hist<histogram,umbrella> xvg2hist(rcmask,rstfunct);
    rstmask = xvg2hist.getcol_rst() >> 2; //2 means 1st column is time; 2nd column is reference group position
    npullgrps = xvg2hist.getnelm();
    fileio<xxvg2hist<histogram,umbrella>,histogram> fxxvg(xvg2hist, std::fstream::in);
    for(uint j = 0; j < nrun; ++j) {
      char runid[MAXNDIGRUN];
      sprintf(runid,"%d",j);
      string xxvgfname = sysname + "_" + runid + "x_" + winid + ".xvg";
      N[i] += fxxvg(hists[i],xxvgfname);
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

  //Construct WHAM
  const valtype kB = 0.0083144621; //kB in kJ/mol/K
  const valtype T = 300; //temperature in kelvin
  vector<Hamiltonian<RST>* > V;
  /*cout << "#rstmask = " << rstmask << endl;
  cout << "#rcmask = " << rcmask << endl;
  cout << "#trajid rstfunct k r" << endl;*/
  for(uint i = 0; i < nwin; ++i) {
    const map<uint,umbrella> funct = rstfuncts[i];
    vector<valtype> k,r;
    for(uint j = 0; j < npullgrps; ++j) { //We need to pack 0 in where there's a RC but not restraint acting on it 
      if( (rstmask & (bitunit << j)).any() ) { //is a restrained RC
        k.push_back(funct.find(j)->second.getk());
        r.push_back(funct.find(j)->second.getinit());
      } else if ( (rcmask & (bitunit << j)).any() ) { //is a RC but not restrained
        k.push_back(0.0); //there's virtually no restraint
	r.push_back(0.0);
      }
    }
    /*cout << "# " << i << " ";
    copy(k.begin(),k.end(),ostream_iterator<valtype>(cout," "));
    copy(r.begin(),r.end(),ostream_iterator<valtype>(cout," "));
    cout << endl;*/

    vector<valtype> params;
    params.push_back(kB);
    params.push_back(T);
    params.insert(params.end(),k.begin(),k.end());
    params.insert(params.end(),r.begin(),r.end());
    V.push_back(new Hamiltonian<RST>(params));
  }
  /*cout << "# Restraint ensemble of each traj: " << endl;
  cout << "# trajid kB T k r" << endl;
  for(uint i = 0; i < V.size(); ++i) {
    cout << "# " << i << " ";
    vector<valtype> params = V[i]->getens().getparams();
    copy(params.begin(),params.end(),ostream_iterator<valtype>(cout," "));
    cout << endl;
  }*/
  
  cout << "# DBL_MAX supported on this machine is " << DBL_MAX << endl;
  cout << "# DBL_MIN supported on this machine is " << DBL_MIN << endl;
  cout << "# NOTE: WHAM iteration might break if the arguments of a exp() evaluation is ";
  cout << "out of the range [" << MINEXPARG << "," << MAXEXPARG << "]" << endl;

  WHAM<RST, histogram, narray> wham(record,hists, V, N, tol);
  vector<valtype> params;
  params.push_back(kB);
  params.push_back(T);
  params.insert(params.end(),ndim*2,0.0);
  Hamiltonian<RST> V0(params);
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