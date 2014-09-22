#include "ensemble.hpp"
#include "typedefs.hpp"
#include "fileio_utils.hpp"
#include <map>
#include <vector>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iterator>

using namespace std;

//histogram class must provide public member function 'coord2val' that maps a coordinate to corresponding value
//histogram class must have member variable 'binsize' which is the size of bin of the histogram
//narray class must have square bracket operator that maps a index (coordinate) to an element
template<class ensemble, class histogram, class narray>
class DensityOfState {
  public:
    DensityOfState();
    DensityOfState(const histogram& hist);

    typedef typename narray::iterator iterator;
    typedef typename narray::const_iterator const_iterator;

    //! calculate density of state given a set of inputs 
    /**
     * @param[in] itC iterator to the numerator array 
     * @param[in] itsNgexpmH iterator to the denumerator N/g*exp(-H) array 
     * @param[in] expf Inverse of the partition function of each state 
     */
    inline
    #ifdef __GNUC__
    __attribute__((always_inline))
    #endif
    auto operator() (
    		      const_iterator& itC, 
                      vector<const_iterator>& itsNgexpmH,  
 		      const vector<valtype>& expf
		    ) -> decltype ( itC->second )
    {
      const auto& num = itC->second;
      ++itC;
      valtype denum = 0.0;
      for(uint k = 0; k < expf.size(); ++k) {
        denum += itsNgexpmH[k]->second * expf[k];
	++itsNgexpmH[k];
      }
      return num/denum;
    }

    //! update density of state given a set of inputs 
    /**
     * @param[in] itC iterator to the numerator array 
     * @param[in] itsNgexpmH iterator to the denumerator N/g*exp(-H) array 
     * @param[in] expf Inverse of the partition function of each state 
     * @param[in] itsexpmH iterator to the denumerator exp(-H) array 
     * @param[in] newexpf updated version of expf 
     */
    inline
    #ifdef __GNUC__
    __attribute__((always_inline))
    #endif
    void operator() (
    		      const_iterator& itC, 
                      vector<const_iterator>& itsNgexpmH,  
 		      const vector<valtype>& expf,
                      vector<const_iterator>& itsexpmH,  
		      vector<valtype>& newexpf
		    )
    {
      newexpf.assign(newexpf.size(), 0.0);
      //Here we loop through all the non-zero histogram values and compute the density of state from the histograms
      for(iterator itDOS = DOS.begin(); itDOS != DOS.end(); ++itDOS ) {
	//This will advance the iterators itC, itsNgexpmH[k]
	const auto dos = this->operator()(itC, itsNgexpmH, expf);
        itDOS->second = dos;
        for(uint l = 0; l < newexpf.size(); ++l) {
          newexpf[l] += dos*itsexpmH[l]->second;
	  ++itsexpmH[l];
        }
      }
      for(auto& x : newexpf) { x = 1/x; }
    }

    //! read only accessor operator 
    valtype operator[](const coordtype& coord) const { return DOS[coord]; };
    //! writable accessor
    valtype& operator[](const coordtype& coord) { return DOS[coord]; };
    //! return DOS
    const narray& getdosarr() const { return DOS;}
    iterator begin() { return DOS.begin(); }
    const_iterator begin() const { return DOS.begin(); }
    iterator end() { return DOS.end(); }
    const_iterator end() const { return DOS.end(); }
    
    void print() const { DOS.print(); }
  private:
    narray DOS;
};

template<class ensemble, class histogram, class narray>
DensityOfState<ensemble,histogram,narray>::DensityOfState():
  DOS(narray())
{}

template<class ensemble, class histogram, class narray>
DensityOfState<ensemble,histogram,narray>::DensityOfState(const histogram& hist)
  : DOS(narray(hist.getdim(),hist.getnelms(),hist.gethv(),hist.getlv())) 
  {}
