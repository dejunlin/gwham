#include "ensemble.hpp"
#include "typedefs.hpp"
#include "gnarray.hpp"
#include "gwham.hpp"
#include "hamiltonian.hpp"
#include "mc.hpp"
#include "fileio_utils.hpp"
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <string>
#include <map>

int main(int argc, char* argv[]) {
  if (argc != 10) {
    cerr << "Usage: mc dim nbins hv lv nwins nsteps stepsize tol dim_prt\n";
    exit(-1);
  }

  cout << "# mc dim nbins hv lv nwins nsteps stepsize tol dim_prt\n";
  cout << "# ";
  copy(argv,argv+argc,ostream_iterator<char*>(cout," "));
  cout << endl;

  const uint dim = atoi(argv[1]);
  const string nbinstr = string(argv[2]); 
  const string hvstr = string(argv[3]); 
  const string lvstr = string(argv[4]); 
  const string nwinsstr = string(argv[5]); 
  const uint nsteps = atoi(argv[6]);
  const valtype stepsize = atof(argv[7]);
  const valtype tol = atof(argv[8]);
  const string dimprtstr = string(argv[9]);

  vector<uint> nbins, nwins, dimprt;
  vector<valtype> hv, lv;
  parser<uint>(nbins,nbinstr);
  parser<uint>(nwins,nwinsstr);
  parser<uint>(dimprt,dimprtstr);
  parser<valtype>(hv,hvstr);
  parser<valtype>(lv,lvstr);

  if(nbins.size() != hv.size() || hv.size() != lv.size() || lv.size() != nwins.size() || nwins.size() != dim) {
    cerr << "The size of nbins, hv, lv, nwins aren't all " << dim << endl;
    exit(-1);
  }

  for(uint i = 0; i < dimprt.size(); ++i) {
    if(dimprt[i] >= dim) {
      cerr << "The dimension to be printed out must be less than " << dim << endl;
      exit(-1);
    }
  }

  const valtype kB = 0.0019872041; //boltzmann factor in kcal/mol/K
  const valtype T = 300; //temperature in K
  
  vector<valtype> params;
  params.push_back(1.0);
  params.push_back(0.001);
  params.push_back(-0.1);
  const potpoly_dblwell pot(dim,params);

  //const nobias bias(dim);
  const valtype spring = 0.2;
  const vector<valtype> k(dim,spring);
  vector<Hamiltonian<RST>* > V;
  vector<valtype> rstparams;
  rstparams.push_back(kB);
  rstparams.push_back(T);
  rstparams.insert(rstparams.end(),k.begin(),k.end());
  
  narray windows(dim,nwins,hv,lv);
  const ulong nwins_tot = windows.getnelms_tot();
  vector<histogram> hists(nwins_tot,histogram(dim,nbins,hv,lv));
  
  const vector<vector<valtype> > wincentrs = windows.canonical_valseries();
  for(uint i = 0; i < wincentrs.size(); ++i) {
    const vector<valtype> wincentr = wincentrs[i];
    harmonic bias(dim,k,wincentr);
    MC<potpoly_dblwell,harmonic> mc(kB,T,pot,bias,nsteps,stepsize);
    mc(wincentr,hists[i]);
    
    rstparams.resize(2+dim);
    rstparams.insert(rstparams.end(),wincentr.begin(),wincentr.end());
    V.push_back(new Hamiltonian<RST>(rstparams));
  }

  histogram hist_sum(dim,nbins,hv,lv);

  map<coordtype, vector<uint> > record;
  for(uint i = 0; i < nwins_tot; ++i) {
    typename histogram::const_iterator it;
    for(it = hists[i].begin(); it != hists[i].end(); ++it) {
      const coordtype coord = it->first;
      if(hist_sum.find(coord) != hist_sum.end()) {
        hist_sum[coord] += it->second;
      } else { 
        hist_sum[coord] = it->second;
      }
      record[coord].push_back(i);
    }
  }
  //print the histograms
  /*for(uint i = 0; i < nwin; ++i) {
    cout << "#Histogram " << i << endl;
    hists[i].print();
  }
  cout << "#Histogram_tot: " << endl;
  hist_sum.print();*/

  WHAM<RST,histogram,narray> wham(record,hists,V,vector<uint>(nwins_tot,nsteps),tol);
  rstparams.resize(2);
  const vector<valtype> k0(dim,0.0);
  const vector<valtype> r0(dim,0.0);
  rstparams.insert(rstparams.end(),k0.begin(),k0.end());
  rstparams.insert(rstparams.end(),r0.begin(),r0.end());
  Hamiltonian<RST> V0(rstparams);
  /*vector<uint> rhodim(dim,0);
  for(uint i = 0; i < dim; ++i) { rhodim[i] = i; }*/
  narray rho = wham.calrho(dimprt,V0);
  //Normalize the probability (just for comparision with other programs)
  printf("#%10s%30s%30s%30s\n","Bin","Vals","PMF","RhoNormalized");
  valtype sum = 0.0;
  narray::iterator itmax = rho.begin();
  for(narray::iterator it = rho.begin(); it != rho.end(); ++it) {
    if(it->second > itmax->second) { itmax = it;}
    sum += it->second;
  }
  //Output normalized density and PMF
  //Evaluate the difference between WHAM results and analytic result
  vector<valtype> dpmf;
  valtype mindpmf = 9999999999, maxdpmf = -9999999999;
  //potprt is the potential in the dimension we want to print out
  const potpoly_dblwell potprt(dimprt.size(),params);
  for(narray::iterator it = rho.begin(); it != rho.end(); ++it) {
    const coordtype bin = it->first;
    const vector<valtype> val = rho.coord2val(bin);

    const valtype pmf = -kB*T*log(it->second/itmax->second);
    const valtype pmf_analytic = potprt(val);
    const valtype d = pmf-pmf_analytic;
    if( d < mindpmf) { mindpmf = d; }
    else if(d > maxdpmf) { maxdpmf = d;}
    dpmf.push_back(pmf-pmf_analytic);
    const valtype rhonorm = it->second/sum;
    for(uint i = 0; i < bin.size(); ++i) { printf("%10d",bin[i]);}
    for(uint i = 0; i < val.size(); ++i) { printf("%30.15lf",val[i]);}
    printf("%30.15lf%30.15lf\n",pmf,rhonorm);
  }
  //perform WHAM inconsistency test
  vector<valtype> eitas = wham.whamvsraw(rho);
  printf("#%10s%30s%30s\n", "Window", "Vals", "Eita");
  for(uint i = 0; i < eitas.size(); ++i) {
    printf("%10d", i);
    const vector<valtype> rstcenter = wincentrs[i];
    for(uint k = 0; k < rstcenter.size(); ++k) { printf("%30lf", rstcenter[k]); }
    printf("%30lf\n", eitas[i]);
  }

  cout << "# WHAM - Analytic " << endl;
  histogram dpmf_hist(1,vector<uint>(1,100),vector<valtype>(1,maxdpmf),vector<valtype>(1,mindpmf));
  for(uint i = 0; i < dpmf.size(); ++i) {
    dpmf_hist.bin(vector<valtype>(1,dpmf[i]));
  }
  dpmf_hist.print();
  for(uint i = 0; i < wincentrs.size(); ++i) {
    delete V[i];
  }
}
