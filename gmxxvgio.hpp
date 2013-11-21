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

//! read gromacs x.xvg file and bin reaction coordinates into histogram
/** rstfunct class should have member 'valtype operator() (const valtype& data)'
 *  all columns in col_rc will be histogramed
 *  all potential energy corresponding to the columns in col_rst will be histogrammed
 *  unless those which also appear in col_rc, where the RC itself will be histogrammed
 */
template<class histogram, class rstfunct>
class xxvg2hist {
  public:
    xxvg2hist(const bitset<MAXNRST>& _col_rc, const map<uint, rstfunct>& _funct);
    linecounter operator() (fstream& fs, histogram& hist) const;
    const bitset<MAXNRST> mapkeys2bitset() const;
    const bitset<MAXNRST>& getcol_rst() const;
    uint getnelm() const;
  private:
    //!pullgroup id => restraint potential functors (pullgroup id starts with 0)
    const map<uint, rstfunct> funct;
    //!bitmask what columns we want to histogram
    const bitset<MAXNRST> col_rc;
    //!bitmask what columns are restrained RC
    const bitset<MAXNRST> col_rst;
    //!bitmask whatever in col_rst but not in col_rc
    const bitset<MAXNRST> col_pot;
    //!number of element to expect, which is number of 1-bit's in col_rc plus that in col_pot
    const uint nelm;
};

template<class histogram, class rstfunct>
const bitset<MAXNRST>& xxvg2hist<histogram, rstfunct>::getcol_rst() const {
  return col_rst;
}

template<class histogram, class rstfunct>
uint xxvg2hist<histogram, rstfunct>::getnelm() const {
  return nelm;
}

template<class histogram, class rstfunct>
const bitset<MAXNRST> xxvg2hist<histogram, rstfunct>::mapkeys2bitset() const {
  bitset<MAXNRST> keybit; //we'll need this for the next check
  typename map<uint,rstfunct>::const_iterator it;
  for(it = funct.begin(); it != funct.end(); ++it) {
    keybit |= (bitunit << it->first);
  }
  return keybit;
}

template<class histogram, class rstfunct>
xxvg2hist<histogram, rstfunct>::xxvg2hist(const bitset<MAXNRST>& _col_rc, const map<uint, rstfunct>& _funct) :
  funct(_funct),
  col_rc(_col_rc << 2), //2 means the 1st column is time and the 2nd column is absolute position of the reference group; TODO: we should extend this class to handle RC in all 3 dimension
  col_rst(this->mapkeys2bitset() << 2),
  col_pot(col_rst & (~col_rc)),
  nelm(col_rc.count() + col_pot.count()) //just count how many '1's are there in both bitmasks
{
  /*cout << "col_rc = " << col_rc << endl;
  cout << "col_rst = " << col_rst << endl;
  cout << "col_pot = " << col_pot << endl;*/
}

template<class histogram, class rstfunct>
linecounter xxvg2hist<histogram, rstfunct>::operator() (fstream& fs, histogram& hist) const {
  //check if the histogram is of the same dimension as the data we're about to read
  const uint nelm_hist = hist.getdim();
  if(nelm_hist != nelm) { 
    cerr << "Asking to read data of " << nelm << " columns but histogram is expecting " << nelm_hist << " columns." << endl; 
    exit(-1); 
  }
  string line;
  linecounter nsamples = 0;
  while(getline(fs,line)) {
    if(line[0] == '@' || line[0] == '#') { continue; }
    //if(!nsamples) { ++nsamples; continue; } //always skip the 1st line, which is the last line from the last file
    vector<valtype> tmp;
    parser<valtype>(tmp,line);
    vector<valtype> data(0,0.0);
    for(uint i = 0; i < tmp.size(); ++i) {
      if( ((bitunit << i) & col_pot).any() ) { data.push_back((funct.find(i-2)->second(tmp[i]))); } //2 means the 1st column is time and the 2nd column is reference group position; TODO
      else if( ((bitunit << i) & col_rc).any() ) { data.push_back(tmp[i]); }
    }
    if(hist.bin(data)) {++nsamples;}
  }
  if(!nsamples) { cerr << "All data is excluded. Please check the histogram bounds\n"; exit(-1);}
  return nsamples;
}
#endif
