#if !defined(GMXXVGIO_HPP)
#define GMXXVGIO_HPP

#include "typedefs.hpp"
#include "fileio_utils.hpp"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <bitset>
#include <vector>

using namespace std;

//! Basic gromacs x.xvg file I/O interface
template<class Tdata, class rstfunct> 
class xxvg {
  private:
    //! convert the keys in xxvg::funct into a bitmask
    const bitset<MAXNRST> mapkeys2bitset() const;
  public:
    xxvg(const uint _ncolskip, const bitset<MAXNRST>& _col_rc, const map<uint, vector<rstfunct*> >& _funct, const uint& _stride);
    virtual linecounter operator() (fstream& fs, Tdata& data) const = 0;
    //! Just return col_rst
    const bitset<MAXNRST>& getcol_rst() const;
    //! Shift the col_rst to the left by ncolskip and return it, which is the pullgrp mask
    const bitset<MAXNRST>& getpullgrp_mask() const;
    //! Return xxvg::nelm
    uint getnelm() const;
  protected:
    //!pullgroup id => restraint potential functors for each ensemble (pullgroup id starts with 0)
    const map<uint, vector<rstfunct*> > funct;
    //!how many columns since the 1st one in an x.xvg file we want to skip
    const uint ncolskip;
    //!bitmask what columns we want to process 
    const bitset<MAXNRST> col_rc;
    //!bitmask what columns are restrained RC
    const bitset<MAXNRST> col_rst;
    //!bitmask whatever in col_rst but not in col_rc
    const bitset<MAXNRST> col_pot;
    //!number of element in each data point to bin into histogram, which is number of 1-bit's in col_rc plus that in col_pot
    const uint nelm;
    //!only process every stride lines from the file
    const uint stride;
};

template<class Tdata, class rstfunct>
const bitset<MAXNRST>& xxvg<Tdata, rstfunct>::getcol_rst() const {
  return col_rst;
}

template<class Tdata, class rstfunct>
const bitset<MAXNRST>& xxvg<Tdata, rstfunct>::getpullgrp_mask() const {
  return col_rst >> ncolskip;
}

template<class Tdata, class rstfunct>
uint xxvg<Tdata, rstfunct>::getnelm() const {
  return nelm;
}

template<class Tdata, class rstfunct>
const bitset<MAXNRST> xxvg<Tdata, rstfunct>::mapkeys2bitset() const {
  bitset<MAXNRST> keybit; //we'll need this for the next check
  typename map<uint,vector<rstfunct*> >::const_iterator it;
  for(it = funct.begin(); it != funct.end(); ++it) {
    keybit |= (bitunit << it->first);
  }
  return keybit;
}

template<class Tdata, class rstfunct>
xxvg<Tdata, rstfunct>::xxvg(const uint _ncolskip, const bitset<MAXNRST>& _col_rc, const map<uint, vector<rstfunct*> >& _funct, const uint& _stride) :
  funct(_funct),
  ncolskip(_ncolskip),
  col_rc(_col_rc << ncolskip), //input _col_rc usu. mask pull-group id but here we want column id mask, which is pull-group id mask + ncolskip offset
  col_rst(this->mapkeys2bitset() << ncolskip),
  col_pot(col_rst & (~col_rc)),
  nelm(col_rc.count() + col_pot.count()), //just count how many '1's are there in both bitmasks
  stride(_stride)
{
  /*cout << "col_rc = " << col_rc << endl;
  cout << "col_rst = " << col_rst << endl;
  cout << "col_pot = " << col_pot << endl;*/
}

//! read gromacs x.xvg file and store a trajectory of energy differences between ensemble
/** NOTE that the dE[t][i] produce by this operator is E[i](x[t]), where 
 * x[t] is the coordinates at time t and E[i] is the i'th restraint potential function. 
 * This is different from what's expected by GMBAR class but dE[t][j] - dE[t][i] always 
 * give back E[j](x[t]) - E[i](x[t]) as is used by GMBAR anyway
 */
template<class rstfunct>
class xxvg2dE : public virtual xxvg<vector<vector<valtype> >, rstfunct> {
  using xxvg<vector<vector<valtype> >,rstfunct>::funct;
  using xxvg<vector<vector<valtype> >,rstfunct>::ncolskip;
  using xxvg<vector<vector<valtype> >,rstfunct>::col_rc;
  using xxvg<vector<vector<valtype> >,rstfunct>::col_rst;
  using xxvg<vector<vector<valtype> >,rstfunct>::col_pot;
  using xxvg<vector<vector<valtype> >,rstfunct>::stride;
  public:
    xxvg2dE(const uint _ncolskip, const bitset<MAXNRST>& _col_rc, const map<uint, vector<rstfunct*> >& _funct, const uint& _stride);
    linecounter operator() (fstream& fs, vector<vector<valtype> >& dE) const;
  private:
};

template<class rstfunct>
xxvg2dE<rstfunct>::xxvg2dE(const uint _ncolskip, const bitset<MAXNRST>& _col_rc, const map<uint, vector<rstfunct*> >& _funct, const uint& _stride) :
  xxvg<vector<vector<valtype> >,rstfunct>(_ncolskip, _col_rc, _funct, _stride)
{
}

template<class rstfunct>
linecounter xxvg2dE<rstfunct>::operator() (fstream& fs, vector<vector<valtype> >& dE) const {
  string line;
  linecounter nsamples = 0;
  linecounter lc = 0;
  while(getline(fs,line)) {
    if(line[0] == '@' || line[0] == '#') { continue; }
    ++lc;
    if(lc % stride) { continue; }
    //if(!nsamples) { ++nsamples; continue; } //always skip the 1st line, which is the last line from the last file
    vector<valtype> tmp;
    parser<valtype>(tmp,line);
    vector<valtype> E(0,0);
    for(uint i = ncolskip; i < tmp.size(); ++i) { //we skip the 1st ncolskip column
      if( ((bitunit << i) & col_rst).any() ) {
	vector<rstfunct*> functs = funct.find(i-ncolskip)->second; 
	if(!E.size()) { E.resize(functs.size(),0.0);}
	for(uint j = 0; j < functs.size(); ++j) {
	  E[j] += functs[j]->operator()(tmp[i]);
	}
      }
    }
    //cout <<"#GMXXVG: dE.size() = " << dE.size() << endl;
    //cout <<"#GMXXVG: E.size() = " << E.size() << endl;
    dE.push_back(E);
    ++nsamples;
  }
  return nsamples;
}

//! read gromacs x.xvg file and bin reaction coordinates into histogram
/** rstfunct class should have member 'valtype operator() (const valtype& data)'
 *  all columns in col_rc will be histogramed
 *  all potential energy corresponding to the columns in col_rst will be histogrammed 
 *  (each potential energy usu. corresponds to a restraint from 1 simulation and there are multiple simulations)
 *  unless those which also appear in col_rc, where the RC itself will be histogrammed
 *  A data point is binned into the histogram with the format:
 *  RC1 RC2 ... RCn POT1 POT2 ... POTm
 *  where there're n RC and m pot
 */
template<class histogram, class rstfunct>
class xxvg2hist : public virtual xxvg<histogram,rstfunct> {
  using xxvg<histogram,rstfunct>::funct;
  using xxvg<histogram,rstfunct>::ncolskip;
  using xxvg<histogram,rstfunct>::col_rc;
  using xxvg<histogram,rstfunct>::col_rst;
  using xxvg<histogram,rstfunct>::col_pot;
  using xxvg<histogram,rstfunct>::stride;
  public:
    xxvg2hist(const uint _ncolskip, const bitset<MAXNRST>& _col_rc, const map<uint, vector<rstfunct*> >& _funct, const uint& _stride);
    linecounter operator() (fstream& fs, histogram& hist) const;
  private:
};

template<class histogram, class rstfunct>
xxvg2hist<histogram, rstfunct>::xxvg2hist(const uint _ncolskip, const bitset<MAXNRST>& _col_rc, const map<uint, vector<rstfunct*> >& _funct, const uint& _stride) :
  xxvg<histogram,rstfunct>(_ncolskip, _col_rc, _funct, _stride)
{
}

template<class histogram, class rstfunct>
linecounter xxvg2hist<histogram, rstfunct>::operator() (fstream& fs, histogram& hist) const {
  //check if the histogram is of the same dimension as the data we're about to read
  const uint nelm_hist = hist.getdim();
  const uint nelm_data = col_pot.any() ? (funct.begin()->second).size() + col_rc.count() : col_rc.count();
  if(nelm_hist != nelm_data) { 
    cerr << "Asking to generate data of " << nelm_data << " columns but histogram is expecting " << nelm_hist << " columns." << endl; 
    exit(-1); 
  }
  string line;
  linecounter nsamples = 0;
  linecounter lc = 0;
  while(getline(fs,line)) {
    if(line[0] == '@' || line[0] == '#') { continue; }
    ++lc;
    if(lc % stride) { continue; }
    //if(!nsamples) { ++nsamples; continue; } //always skip the 1st line, which is the last line from the last file
    vector<valtype> tmp;
    parser<valtype>(tmp,line);
    vector<valtype> data(0,0.0);
    vector<valtype> tmpdata(0,0);
    for(uint i = ncolskip; i < tmp.size(); ++i) {//we skip the 1st ncolskip column
      //Save the pot in tmpdata first; after we've processed all the columns, we'll push them in the data
      if( ((bitunit << i) & col_pot).any() ) {
	vector<rstfunct*> functs = funct.find(i-ncolskip)->second; 
	if(!tmpdata.size()) { tmpdata.resize(functs.size(),0.0);}
	for(uint j = 0; j < functs.size(); ++j) {
	  tmpdata[j] += functs[j]->operator()(tmp[i]);
	}
      }
      //We just push the RC in the data first
      else if( ((bitunit << i) & col_rc).any() ) { data.push_back(tmp[i]); } 
    }
    //Here we push the pot in the data
    data.insert(data.end(),tmpdata.begin(),tmpdata.end());
    /*cout << "data = ";
    copy(data.begin(),data.end(),ostream_iterator<valtype>(cout," ")); cout << endl;*/
    if(hist.bin(data) != hist.end()) {++nsamples;}
  }
  if(!nsamples) { cerr << "All data is excluded. Please check the histogram bounds\n"; exit(-1);}
  return nsamples;
}

#endif
