#include "ensemble.hpp"
#include "typedefs.hpp"
#include <map>
#include <vector>
#include <math.h>
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
    void operator() (
    		      const narray& C, 
                      const vector<narray>& expmH,  
                      const vector<narray>& Ng,  
 		      const vector<valtype>& expf,
		      vector<valtype>& newexpf
		    );
    //! read only accessor operator 
    valtype operator[](const coordtype& coord) const { return DOS[coord]; };
    //! return DOS
    const narray& getdosarr() const { return DOS;}
    typedef typename narray::iterator iterator;
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

template<class ensemble, class histogram, class narray>
void DensityOfState<ensemble,histogram,narray>::operator () (
    		     	      const narray& C, 
                              const vector<narray>& expmH,  
                              const vector<narray>& Ng,  
         		      const vector<valtype>& expf,
			      vector<valtype>& newexpf
         		                   ) 
{
  newexpf.assign(newexpf.size(), 0.0);
  //Here we loop through all the non-zero histogram values and compute the density of state from the histograms
  for(typename narray::const_iterator it = C.begin(); it != C.end(); ++it) {
    const auto& coord = it->first;
    const auto& num = it->second;
    valtype denum = 0.0;
    for(uint k = 0; k < Ng.size(); ++k) {
      denum += Ng[k][coord] * expf[k] * expmH[k][coord];
    }
    const valtype dos = num/denum;
    DOS[coord] = dos;
    for(uint l = 0; l < newexpf.size(); ++l) {
      newexpf[l] += dos*expmH[l][coord];
    }
  }
  for(auto& x : newexpf) { x = 1/x; }
}
