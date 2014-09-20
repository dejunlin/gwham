#include "ensemble.hpp"
#include "typedefs.hpp"
#include "gnarray.hpp"
#include "gwham.hpp"
#include "hamiltonian.hpp"
#include "mc.hpp"
#include "fileio_utils.hpp"
#include "metaprog_snippets.hpp"
#include <algorithm>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <string>
#include <map>
#include <iomanip>

int main(int argc, char* argv[]) {
  if (argc != 10) {
    cerr << "Usage: mc dim nbins hv lv nwins nsteps stepsize tol dim_prt\n";
    exit(-1);
  }

  cout << "# mc dim nbins hv lv nwins nsteps stepsize tol dim_prt\n";
  cout << "# ";
  for_each(argv, argv+argc, [] (char* arg) { cout << "'" << arg << "' "; });
  cout << endl;

  const uint dim = atoi(argv[1]);
  const string nbinstr = string(argv[2]); 
  const string hvstr = string(argv[3]); 
  const string lvstr = string(argv[4]); 
  const string nwinsstr = string(argv[5]); 
  const linecounter nsteps = stoul(argv[6]);
  const valtype stepsize = atof(argv[7]);
  const valtype tol = atof(argv[8]);
  const string dimprtstr = string(argv[9]);

  vector<uint> nbins, nwins, dimprt;
  vector<valtype> hv, lv;
  split<uint>(nbins,nbinstr);
  split<uint>(nwins,nwinsstr);
  split<uint>(dimprt,dimprtstr);
  split<valtype>(hv,hvstr);
  split<valtype>(lv,lvstr);

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
  vpEnsemble V;
  
  narray windows(dim,nwins,hv,lv);
  const ulong nwins_tot = windows.getnelms_tot();
  vector<histogram> hists(nwins_tot,histogram(dim,nbins,hv,lv));
  
  const vector<vector<valtype> > wincentrs = windows.canonical_valseries();
  for(uint i = 0; i < wincentrs.size(); ++i) {
    const vector<valtype> wincentr = wincentrs[i];
    harmonic bias(dim,k,wincentr);
    MC<potpoly_dblwell,harmonic> mc(kB,T,pot,bias,nsteps,stepsize);
    mc(wincentr,hists[i]);
    
    vFunctVV functors;
    for(uint j = 0; j < dim; ++j) {
      functors.emplace_back(Make_Functor_Wrapper(FunctVV{}, Quadratic<valtype>, k[j]/2, wincentr[j], 0.0));
    }
    V.emplace_back(make_shared<NVT>(kB, Hamiltonian{functors}, T));
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

  WHAM<pEnsemble,histogram,narray> wham(record,hists,V,vector<linecounter>(nwins_tot,nsteps),tol);
  pEnsemble V0 = make_shared<NVT>(kB, Hamiltonian{}, T);
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
    fcout.width(10);
    fcout << bin;
    fcout.flags(ios::fixed | ios::right);
    fcout.precision(15);
    fcout.width(30);
    fcout << val << pmf;
    fcout.flags(ios::scientific);
    fcout << rhonorm << endl;
  }
  //perform WHAM inconsistency test
  vector<valtype> eitas = wham.whamvsraw(rho);
  printf("#%10s%30s%30s\n", "Window", "Vals", "Eita");
  for(uint i = 0; i < eitas.size(); ++i) {
    fcout.width(10);
    fcout << i;
    const vector<valtype> rstcenter = wincentrs[i];
    fcout.width(30);
    fcout.precision(15);
    fcout.flags(ios::fixed);
    fcout << rstcenter << eitas[i] << endl;
  }

  cout << "# WHAM - Analytic " << endl;
  histogram dpmf_hist(1,vector<uint>(1,100),vector<valtype>(1,maxdpmf),vector<valtype>(1,mindpmf));
  for(uint i = 0; i < dpmf.size(); ++i) {
    dpmf_hist.bin(vector<valtype>(1,dpmf[i]));
  }
  dpmf_hist.print();
}
