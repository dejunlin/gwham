#if !defined(GNARRAY_HPP)
#define GNARRAY_HPP

#include <map>
#include <iterator>
#include <vector>
#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <algorithm>
#include "metaprog_snippets.hpp"

using namespace std;

template<class Tcoord, class Telem, class Tval>
class gnarray {
  public:
    //!just a wrapper around map<Tcoord,Telem>::iterator
    typedef typename map<Tcoord, Telem>::iterator iterator;
    //!just a wrapper around map<Tcoord,Telem>::const_iterator
    typedef typename map<Tcoord, Telem>::const_iterator const_iterator;
    //!just Tcoord
    typedef Tcoord gridcoord;
    //!just Telem
    typedef Telem gridval;
    //!This type
    typedef gnarray<Tcoord, Telem, Tval> ThisType;
    //!empty constructor
    gnarray();
    //!constructor
    gnarray(const uint _dim, const vector<uint>& _nelms, const vector<Tval>& _hv, const vector<Tval>& _lv);
    //!Copy constructor
    gnarray(const gnarray<Tcoord,Telem,Tval>& _gnarr);
    //!constructing by splicing another narray. Note that this will construct an empty narr
    gnarray(const gnarray<Tcoord,Telem,Tval>& _gnarr, const vector<uint>& _dim);
    //!given coordinate in index space, calculate real-space coordinate
    const vector<Tval> coord2val(const Tcoord& coord) const;
    //!given real-space coordinate, calculate coordinate in index space
    const Tcoord val2coord(const vector<Tval>& vals) const;
    //!bin a data point into histogram and increase the value by weight
    iterator bin(const vector<Tval>& data, const Telem& weight = 1);
    //!just a wrapper around map<Tcoord,Telem>::find()
    iterator find(const Tcoord& coord);
    //!just a wrapper around map<Tcoord,Telem>::find()
    const_iterator find(const Tcoord& coord) const;
    //!just a wrapper around map<Tcoord,Telem>::begin()
    iterator begin();
    //!just a wrapper around map<Tcoord,Telem>::begin()
    const_iterator begin() const;
    //!just a wrapper around map<Tcoord,Telem>::end()
    iterator end();
    //!just a wrapper around map<Tcoord,Telem>::end()
    const_iterator end() const;
    //!just a wrapper around map<Tcoord,Telem>::size()
    size_t size() const;

    //! read-only accessor operator
    Telem operator[](const Tcoord& coord) const {
      const_iterator it = narr.find(coord);
      return it->second;
    }
    //! writable accessor operator
    Telem& operator[](const Tcoord& coord) { return narr[coord]; }
    //! assignment operator
    gnarray operator=(const gnarray<Tcoord,Telem,Tval>& rhs);

    //!Return the narr
    const map<Tcoord,Telem>& getnarr() const { return narr;}
    //!Return the number of dimension
    uint getdim() const { return dim; }
    //!Return the nelms
    const vector<uint>& getnelms() const { return nelms;}
    //!Return the nelms_tot
    const ulong& getnelms_tot() const { return nelms_tot;}
    //!Return the hv
    const vector<Tval>& gethv() const { return hv;}
    //!Return the lv
    const vector<Tval>& getlv() const { return lv;}
    //!Return the binsize
    const vector<Tval>& getbinsize() const { return binsize;}
    //!Return the stride_array
    const vector<uint>& getstride_array() const { return stride_array;}

    //!Set the narr
    void setnarr(const map<Tcoord,Telem>& _narr) { narr = _narr;}
    //!Set the number of dimension
    void setdim(const uint& _dim) { dim = _dim; }
    //!Set the nelms
    void setnelms(const vector<uint>& _nelms) { nelms = _nelms;}
    //!Set the hv
    void sethv(const vector<Tval>& _hv) { hv = _hv;}
    //!Set the lv
    void setlv(const vector<Tval>& _lv) { lv = _lv;}
    //!Set the binsize
    void setbinsize(const vector<Tval>& _binsize) { binsize = _binsize;}

    //!Print out all elements
    void print() const;
    //!Print out all elements divided by norm (normalizing factor in histogram)
    void print(const Telem& norm) const;
    //!Output a serialized array of vals in canonical order
    vector<vector<Tval> > canonical_valseries() const; 
    //!Sum of all the elments
    Telem sum() const;

    //!Calculate the overlap between this and a input narray
    /** The overlap is the sum of the lesser elemnts that reside 
     * in the common subset of coordinates of the 2 narrays
     */
    Telem overlap(const ThisType& rhs) const;

  private:
    map<Tcoord, Telem> narr;
    //!number of dimension
    uint dim;
    //!number of elements or bins along each dimension. For histograming only
    vector<uint> nelms;
    //!upper bound in real space. For histograming only
    vector<Tval> hv;
    //!lower bound in real space. For histograming only
    vector<Tval> lv;
    //!binsize in real-space along each dimension. For histograming only
    vector<Tval> binsize;
    //! stride_array[i+1] = nelms[i+1]*nelms[i+2]*...*nelms[dim-1] or 1 (when i=dim-1), the bin-width of dimension 'i+1' ( -1 <= i <= dim-1 ) in unit of regular vector elements
    vector<uint> stride_array;
    //! populate stride_array
    void setstride();
    //! total number of elements
    ulong nelms_tot;
};
template<class Tcoord, class Telem, class Tval>
const Tcoord gnarray<Tcoord,Telem,Tval>::val2coord(const vector<Tval>& vals) const {
  Tcoord coord;
  for(uint i = 0; i < vals.size(); ++i) {
    if(vals[i] >= hv[i] || vals[i] < lv[i] || vals[i] != vals[i]) { //Here we assume the range to be [lv, hv)
      return Tcoord();
    }
    coord.push_back(uint((vals[i]-lv[i])/binsize[i]));
  }
  return coord;
}

template<class Tcoord, class Telem, class Tval>
typename gnarray<Tcoord,Telem,Tval>::iterator gnarray<Tcoord,Telem,Tval>::bin(const vector<Tval>& data, const Telem& weight) {
  const Tcoord coord = val2coord(data);
  if(coord.size()) { //only bin data if they're in bound
    const iterator it = narr.find(coord);
    if(it != narr.end()) { it->second += weight; return it; } //check this to avoid weird uninitialized value
    else { return narr.insert(it,pair<Tcoord,Telem>(coord,weight)); }
  } 
  return narr.end();
}

template<class Tcoord, class Telem, class Tval>
gnarray<Tcoord,Telem,Tval>::gnarray():
  narr(map<Tcoord,Telem>()),
  dim(0),
  nelms(vector<uint>(0,0)),
  hv(vector<Tval>(0,0)),
  lv(vector<Tval>(0,0)),
  binsize(vector<Tval>(0,0)),
  stride_array(vector<uint>(dim+1,0)),
  nelms_tot(0)
{
  setstride();
}

template<class Tcoord, class Telem, class Tval>
gnarray<Tcoord,Telem,Tval>::gnarray(const uint _dim, const vector<uint>& _nelms, const vector<Tval>& _hv, const vector<Tval>& _lv):
  narr(map<Tcoord,Telem>()),
  dim(_dim),
  nelms(_nelms),
  hv(_hv),
  lv(_lv),
  binsize(vector<Tval>(dim,0)),
  stride_array(vector<uint>(dim+1,0)),
  nelms_tot(1)
{
  for(uint i = 0; i < dim; ++i) {
    binsize[i] = (hv[i] - lv[i])/nelms[i];
    nelms_tot *= nelms[i];
  }
  setstride();
}

template<class Tcoord, class Telem, class Tval>
gnarray<Tcoord,Telem,Tval>::gnarray(const gnarray<Tcoord,Telem,Tval>& _gnarr, const vector<uint>& _dim) :
  narr(map<Tcoord,Telem>()),
  dim(_dim.size()),
  nelms(vector<uint>(0,0)),
  hv(vector<Tval>(0,0)),
  lv(vector<Tval>(0,0)),
  binsize(vector<Tval>(0,0)),
  stride_array(vector<uint>(dim+1,0)),
  nelms_tot(1)
{
  for(uint i = 0; i < _dim.size(); ++i) {
    const uint dimindex = _dim[i];
    nelms.push_back(_gnarr.getnelms()[dimindex]);
    hv.push_back(_gnarr.gethv()[dimindex]);
    lv.push_back(_gnarr.getlv()[dimindex]);
    binsize.push_back(_gnarr.getbinsize()[dimindex]);
    nelms_tot *= _gnarr.getnelms()[dimindex];
  }
  setstride();
}

template<class Tcoord, class Telem, class Tval>
gnarray<Tcoord,Telem,Tval>::gnarray(const gnarray<Tcoord,Telem,Tval>& _gnarr) :
  narr(_gnarr.getnarr()),
  dim(_gnarr.getdim()),
  nelms(_gnarr.getnelms()),
  hv(_gnarr.gethv()),
  lv(_gnarr.getlv()),
  binsize(_gnarr.getbinsize()),
  stride_array(_gnarr.getstride_array()),
  nelms_tot(_gnarr.getnelms_tot())
{
  setstride();
}

template<class Tcoord, class Telem, class Tval>
const vector<Tval> gnarray<Tcoord,Telem,Tval>::coord2val(const Tcoord& coord) const {
  vector<Tval> vals;
  for(uint i = 0; i < dim; ++i) {
    vals.push_back(lv[i]+(coord[i]+0.5)*binsize[i]);
  }
  return vals;
}

template<class Tcoord, class Telem, class Tval>
typename gnarray<Tcoord,Telem,Tval>::iterator gnarray<Tcoord,Telem,Tval>::find(const Tcoord& coord) {
  return narr.find(coord);
}

template<class Tcoord, class Telem, class Tval>
typename gnarray<Tcoord,Telem,Tval>::const_iterator gnarray<Tcoord,Telem,Tval>::find(const Tcoord& coord) const {
  return narr.find(coord);
}

template<class Tcoord, class Telem, class Tval>
typename gnarray<Tcoord,Telem,Tval>::iterator gnarray<Tcoord,Telem,Tval>::begin() {
  return narr.begin();
}

template<class Tcoord, class Telem, class Tval>
typename gnarray<Tcoord,Telem,Tval>::const_iterator gnarray<Tcoord,Telem,Tval>::begin() const {
  return narr.begin();
}

template<class Tcoord, class Telem, class Tval>
typename gnarray<Tcoord,Telem,Tval>::iterator gnarray<Tcoord,Telem,Tval>::end() {
  return narr.end();
}

template<class Tcoord, class Telem, class Tval>
typename gnarray<Tcoord,Telem,Tval>::const_iterator gnarray<Tcoord,Telem,Tval>::end() const {
  return narr.end();
}

template<class Tcoord, class Telem, class Tval>
size_t gnarray<Tcoord,Telem,Tval>::size() const {
  return narr.size();
}

template<class Tcoord, class Telem, class Tval>
gnarray<Tcoord,Telem,Tval> gnarray<Tcoord,Telem,Tval>::operator=(const gnarray<Tcoord,Telem,Tval>& rhs) {
  narr = rhs.getnarr();
  dim = rhs.getdim();
  nelms = rhs.getnelms();
  hv = rhs.gethv();
  lv = rhs.getlv();
  binsize = rhs.getbinsize();
  return *this;
}

template<class Tcoord, class Telem, class Tval>
void gnarray<Tcoord,Telem,Tval>::print() const {
  printf("%10s%30s%20s\n","#coord","val","elements");
  for(const_iterator it = narr.begin(); it != narr.end(); ++it) {
    const Tcoord coord = it->first;
    const Telem elem = it->second;
    const vector<Tval> vals = coord2val(coord);
    fcout.width(10);
    fcout << coord;
    fcout.flags(ios::fixed);
    fcout.width(30);
    fcout.precision(15);
    fcout << vals << elem << endl;
  }
}

template<class Tcoord, class Telem, class Tval>
void gnarray<Tcoord,Telem,Tval>::print(const Telem& norm) const {
  printf("%10s%30s%20s\n","#coord","val","elements");
  for(const_iterator it = narr.begin(); it != narr.end(); ++it) {
    const Tcoord coord = it->first;
    const Telem elem = it->second;
    const vector<Tval> vals = coord2val(coord);
    fcout.width(10);
    fcout << coord;
    fcout.flags(ios::fixed);
    fcout.width(30);
    fcout.precision(15);
    fcout << vals << elem/norm << endl;
  }
}

template<class Tcoord, class Telem, class Tval>
vector<vector<Tval> > gnarray<Tcoord,Telem,Tval>::canonical_valseries() const {
  vector<vector<Tval> > vals;
  for(ulong i = 0; i < nelms_tot; ++i) {
    Tcoord coord;
    for(uint j = 0; j < dim; ++j) {
      coord.push_back( (i%stride_array[j]) / stride_array[j+1] );
    }
    vector<Tval> val = coord2val(coord);
    vals.push_back(val);
  }
  return vals;
}

template<class Tcoord, class Telem, class Tval>
Telem gnarray<Tcoord,Telem,Tval>::sum() const {
  Telem ans = 0;
  for(const_iterator it = narr.begin(); it != narr.end(); ++it) {
    ans += it->second;
  }
  return ans;
}

template<class Tcoord, class Telem, class Tval>
Telem gnarray<Tcoord,Telem,Tval>::overlap(const ThisType& rhs) const {
  Telem ans{0};
  for(const_iterator it = narr.begin(); it != narr.end(); ++it) {
    const auto& coord = it->first;
    const auto& val = it->second;
    const auto& itrhs = rhs.find(coord);
    //we only sum the elements whose coordinates are common
    if(itrhs == rhs.end()) continue;
    const auto& valrhs = itrhs->second;
    ans += std::min(val, valrhs);
  }
  return ans;
}

template<class Tcoord, class Telem, class Tval>
void gnarray<Tcoord,Telem,Tval>::setstride() {
  stride_array[dim] = 1;
  if(!dim) { return; }
  for(int i = dim-1; i >= 0; --i) {
    stride_array[i] = nelms[i]*stride_array[i+1];
  }
  return;
}

//! calculate the overlap matrix of a set of narrays
template < class V, class C, class NARRAY >
vector<vector<V> > narrayoverlap(const vector<NARRAY>& hists, const vector<C>& Nsamples) {
  vector<vector<V> > overlap(hists.size(), vector<V>(hists.size(), 0));
  for(uint i = 0; i < hists.size(); ++i) {
    //diagonal is always 1 by definition
    overlap[i][i] = 1;
    const auto& hi = hists[i];
    const auto& Ni = Nsamples[i];
    for(uint j = i+1; j < hists.size(); ++j) {
      const auto ov = hi.overlap(hists[j]);
      //NOTE: there is a chance that V(ov) might 
      //overflow V type in case ov is too large
      const V oij = V(ov)/(Ni+Nsamples[j]);
      overlap[i][j] = oij;
      //just copy the upper half to the lower
      overlap[j][i] = oij;
    }
  }
  return overlap;
}

#endif
