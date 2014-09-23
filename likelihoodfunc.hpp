#ifndef LIKELIHOODFUNC_HPP
#define LIKELIHOODFUNC_HPP
/*
 * =====================================================================================
 *
 *       Filename:  likelyhoodfunc.hpp
 *
 *    Description:  The WHAM likelyhood functor
 *
 *        Version:  1.0
 *        Created:  22/09/14 11:16:36
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dejun Lin (DL), dejun.lin@gmail.com
 *   Organization:  Department of Biochemistry and Biophysics, Medical Center, University of Rochester
 *
 * =====================================================================================
 */
#include <cmath>
using namespace std;

template < class DOS, class NARRAY >
class LikeliHoodFunc {
  typedef typename NARRAY::iterator narrit;
  typedef typename NARRAY::const_iterator narrcit;
  typedef typename DOS::iterator dosit;
  typedef typename DOS::const_iterator doscit;
  typedef decltype( declval<narrit>()->second ) valtype; 
  typedef unsigned int uint;
  typedef unsigned long ulong;
  public:
    LikeliHoodFunc(DOS& _dos, const vector<ulong>& _N, const NARRAY& _C, const vector<NARRAY>& _NgexpmH) 
      : dos(_dos), N(_N), C(_C), NgexpmH(_NgexpmH) 
    {
      if(N.size() != NgexpmH.size()) {
	throw(General_Exception("Size of array N is not the same as array NgexpmH"));
      }
    };

    valtype operator() (const vector<valtype>& deltaf) {
      vector<valtype> f(deltaf.size(), 0.0);
      //here we constraint f[0] to 0.0 so that vector f only 
      //stores the free energy of state 1,2,...,N-1
      vector<valtype> expf(f.size()+1, 0.0);
      expf[0] = 1.0;
      f[0] = deltaf[0];
      expf[1] = exp(f[0]);
      valtype fret = 0.0 - N[1]*f[0];
      for(uint i = 1; i < deltaf.size(); ++i) {
	f[i] = f[i-1] + deltaf[i];
	const auto e = exp(f[i]);
	expf[i+1] = e;
	fret -= N[i+1]*f[i];
      }

      narrcit itC = C.begin();
      vector<narrcit> itsNgexpmH;
      for_each(NgexpmH.begin(), NgexpmH.end(), [&itsNgexpmH](const NARRAY& narr) { itsNgexpmH.emplace_back(narr.begin()); });
      for(dosit itdos = dos.begin(); itdos != dos.end(); ++itdos ) {
	const auto c = itC->second;
	const auto d = dos(itC, itsNgexpmH, expf);
	fret -= c*log(d);
      }
      return fret;
    }
    
    void df (const vector<valtype>& deltaf, vector<valtype>& graddf) {
      vector<valtype> f(deltaf.size(), 0.0);
      f[0] = deltaf[0];
      //here we constraint f[0] to 0.0 so that vector f only 
      //stores the free energy of state 1,2,...,N-1
      vector<valtype> expf(f.size()+1, 0.0);
      expf[0] = 1.0;
      expf[1] = exp(f[0]);
      vector<valtype> gradf(f.size(), 0.0);
      gradf[0] -= N[1];
      for(uint i = 1; i < deltaf.size(); ++i) {
	f[i] = f[i-1] + deltaf[i];
	const auto e = exp(f[i]);
	expf[i+1] = e;
	gradf[i] -= N[i+1];
      }
      narrcit itC = C.begin();
      vector<narrcit> itsNgexpmH;
      for_each(NgexpmH.begin(), NgexpmH.end(), [&itsNgexpmH](const NARRAY& narr) { itsNgexpmH.emplace_back(narr.begin()); });
      for(dosit itdos = dos.begin(); itdos != dos.end(); ++itdos ) {
	vector<narrcit> itsNgexpmH_buff{itsNgexpmH};
	const auto d = dos(itC, itsNgexpmH, expf);
	for(uint i = 0; i < f.size(); ++i) {
	  gradf[i] += d*itsNgexpmH_buff[i+1]->second*expf[i+1];
	  ++itsNgexpmH_buff[i+1];
	}
      }
      graddf.assign(gradf.size(), 0.0);
      for(uint i = 0; i < gradf.size(); ++i) {
	const auto& gf = gradf[i];
	for(uint j = 0; j <= i; ++j) {
	  graddf[j] += gf;
	}
      }
    }
     
  private:
    DOS& dos;
    const vector<ulong>& N;
    const NARRAY& C;
    const vector<NARRAY>& NgexpmH;
};


#endif
