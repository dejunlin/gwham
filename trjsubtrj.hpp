#if !defined(TRJSUBTRJ_HPP)
#define TRJSUBTRJ_HPP

#include "typedefs.hpp"
#include "fileio_utils.hpp"
#include "ensemble.hpp"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <bitset>
#include <vector>
#include <map>
#include <algorithm>
#include <math.h>

template<class histogram>
class trjsubtrj {
  public:
    typedef typename histogram::hist_coord hist_coord;
    trjsubtrj(const uint& _stride, const bitset<MAXNRST>& _col_rc, const NVT& _nvt);
    linecounter operator() (fstream& fs, histogram& hist, map<hist_coord,vector<valtype> >& subtrj) const; 
  private:
    const uint stride;
    const bitset<MAXNRST> col_rc;
    const NVT nvt;
};

template<class histogram>
trjsubtrj<histogram>::trjsubtrj(const uint& _stride, const bitset<MAXNRST>& _col_rc, const NVT& _nvt):
  stride(_stride),
  col_rc(_col_rc << 1), //assuming the 1st column is time
  nvt(_nvt)
{}

template<class histogram>
linecounter trjsubtrj<histogram>::operator() (fstream& fs, histogram& hist, map<hist_coord,vector<valtype> >& subtrj) const {
  string line;
  linecounter nsamples = 0;
  linecounter lc = 0;
  while(getline(fs,line)) {
    if(line[0] == '@' || line[0] == '#') { continue; }
    ++lc;
    if(lc % stride) { continue; }
    vector<valtype> tmp;
    parser<valtype>(tmp,line);
    const valtype bias = tmp[tmp.size()-1];
    vector<valtype> vals;
    for(uint i = 1; i < tmp.size() - 2; ++i) { //assuming 1st column is time and the last column is bias energy
      if( ((bitunit << i) & col_rc).any() ) { vals.push_back(tmp[i]); }
    }
    
    //Here we postpone the binning until we read all data
    //so that we can sort the bias ascendingly to minimize underflow of exp() calls
    //since we might be dealing with Boltzmann factor of very large bias
    hist_coord coord = hist.val2coord(vals);
    if(coord.size()) { 
      ++nsamples;
      if(subtrj.find(coord) != subtrj.end()) { subtrj[coord].push_back(bias); }
      else { subtrj[coord] = vector<valtype>(1,bias); }
    } 
  }
  typename map<hist_coord,vector<valtype> >::iterator it;
  for(it = subtrj.begin(); it != subtrj.end(); ++it) {
    const hist_coord coord = it->first;
    vector<valtype>& bias = it->second;
    sort(bias.begin(),bias.end());
    //Here we manually bin the bias
    for(vector<valtype>::const_iterator i = bias.begin(); i != bias.end(); ++i) {
      if(hist.find(coord) != hist.end()) { hist[coord] += exp(nvt.ener(*i)); } //TODO: we might need to take out the MAX here to prevent overflow...
      else { hist[coord] = exp(nvt.ener(*i)); }
    }
  }
  return nsamples;
}


#endif
