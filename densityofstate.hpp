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
    //! update density of state given a set of inputs (assumes g[k][coord] are invariant across all k for any coord)
    /**
     * @param[in] record The 1st vector<uint> is a set of coordinates in at least one of the histograms, whose values is non-zero; and the 2nd vector<uint> keep the index of the histograms
     * @param[in] hists Generic histogram for each trajectory
     * @param[in] V Generic hamiltonian that combine the conserved quantities and the associated parameters. Total number of elements is the total number of states considered
     * @param[in] N N[k] is the number of samples in the k'th trajectory/state. 
     * @param[in] f Free energy of states. Total number of elements is the total number of states considered 
     */

    inline
    #ifdef __GNUC__
    __attribute__((always_inline))
    #endif
    void operator() (
    		      typename narray::const_iterator& itC, 
                      vector<typename narray::const_iterator>& itsNgexpmH,  
                      vector<typename narray::const_iterator>& itsexpmH,  
 		      const vector<valtype>& expf,
		      vector<valtype>& newexpf
		    )
    {
      newexpf.assign(newexpf.size(), 0.0);
      //Here we loop through all the non-zero histogram values and compute the density of state from the histograms
      for(iterator itDOS = DOS.begin(); itDOS != DOS.end(); ++itDOS ) {
        valtype denum = 0.0;
        const auto& num = itC->second;
	++itC;
        for(uint k = 0; k < expf.size(); ++k) {
          denum += itsNgexpmH[k]->second * expf[k];
	  ++itsNgexpmH[k];
        }
        const valtype dos = num/denum;
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
    typedef typename narray::iterator iterator;
    typedef typename narray::const_iterator const_iterator;
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
