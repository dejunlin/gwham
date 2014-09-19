#if !defined(GWHAM_HPP)
#define GWHAM_HPP

#include "densityofstate.hpp"
#include "hamiltonian.hpp"
#include "exception.hpp"
#include "fileio_utils.hpp"
#include <stdlib.h>
#include <iterator>
#include <vector>

using namespace std;

//HISTOGRAM class must provide public member function 'coord2val' that maps a coordinate to corresponding value
//HISTOGRAM class must have member variable 'binsize' which is the size of bin along each dimension 
//NARRAY class must have square bracket operator that maps a index (coordinate) to an element
//NARRAY class must have member iterator type
//NARRAY class must have member function find(coordtype) which given a coordinate return the iterator it points at
//NARRAY class must have member function end() which return the past-the-end iterator
template <class PENSEMBLE, class HISTOGRAM, class NARRAY>
class WHAM {
  typedef DensityOfState<PENSEMBLE,HISTOGRAM,NARRAY> DOStype;
  typedef DOStype* DOSptr;
  public:
    //!perform WHAM iteration and update the density of state
    /**
     * @param[in] _record record[i][k] is the index of the k'th histograms that has non-zero value at point whose coordinate is i 
     * @param[in] hists Generic HISTOGRAM for each trajectory
     * @param[in] _g[k][coord] Statistical inefficiency of the k'th trajectory in bin coord
     * @param[in] V[i] is the hamiltonian the i'th state that combine the conserved quantities and the associated parameters
     * @param[in] N N[k][l] is the number of samples the k'th trajectory visiting the l'th state. 
     * @param[in] tol Tolerance for WHAM iteration 
     */
    WHAM(const map<coordtype, vector<uint> >& _record, 
         const vector<HISTOGRAM>& _hists, 
	 const vector<PENSEMBLE>& _V, 
         const vector<linecounter>& _N, 
	 const valtype _tol, 
	 const vector<valtype>& fseeds = vector<valtype>(0),
	 const vector<NARRAY>& _g = vector<NARRAY>{0}
	);

    //!calculate PMF in Hamiltonian _V along the dimension in DOS as specified by dim
    NARRAY calrho(const vector<uint>& dim, const pEnsemble& _V) const;
    //!Calculate the inconsistency between the i'th consensus HISTOGRAM and the corresponding raw HISTOGRAM  
    valtype whamvsrawi(const uint& i, const NARRAY& rho) const;
    //!just loop over WHAM::whamvsrawi for all i
    vector<valtype> whamvsraw(const NARRAY& rho) const;
    //!given coordinate in index space, calculate the corresponding real-space coordinate
    const vector<valtype> coord2val(const coordtype& coord) const;
    //!given coordinate in index space along the specified dimension, calculate the corresponding real-space coordinate
    const vector<valtype> coord2val(const vector<uint>& _dim, const coordtype& coord) const;
    //!shift all the free energies (argument newexpf)by a constant so that newexpf[0] is always zero 
    void shiftf(vector<valtype>& newexpf) const;
    //!print out all free energies
    void printfree() const;
  private:
     //!max number of iterations
     static const ulong MAXIT = 1000000;
     //!check if WHAM iteration should end; also update WHAM::expf if not end
     bool endit(const vector<valtype>& newexpf, const ulong& count); 
     //!record[i][k] is the index of the k'th histograms that has non-zero value at point whose coordinate is i 
     const map<coordtype, vector<uint> >* const record;
     //!Generic HISTOGRAM for each trajectory
     const vector<HISTOGRAM>* const hists;
     //!Statistical inefficiency for each trajectory
     const vector<NARRAY>* const g;
     //!V[i] is the hamiltonian the i'th state that combine the conserved quantities and the associated parameters
     const vector<PENSEMBLE>* const V;
     //!N[k] is the number of samples the k'th trajectory has, assuming each trajectory visits 1 state 
     const vector<linecounter>* const N;
     //!Tolerance for WHAM iteration
     const valtype tol;
     //!bin-size of each dimension of the HISTOGRAM
     const vector<valtype> binsize;
     //!lower bound in real space of the HISTOGRAM
     const vector<valtype> lv;
     //!number of dimension of the HISTOGRAM
     const uint dim;
     //!expf[i] is the exponential of the dimensionless free energy of state i 
     vector<valtype> expf;
     //! Density of state as well as a functor that caclulates it
     DOStype DOS;
};

template <class PENSEMBLE, class HISTOGRAM, class NARRAY>
WHAM<PENSEMBLE,HISTOGRAM,NARRAY>::WHAM(const map<coordtype, vector<uint> >& _record, 
                                      const vector<HISTOGRAM>& _hists,  
	                              const vector<PENSEMBLE>& _V, 
                                      const vector<linecounter>& _N, 
	                              const double _tol, 
                                      const vector<valtype>& fseeds,
                                      const vector<NARRAY>& _g 
	                             ):
				     record(&_record),
				     hists(&_hists),
				     g(_g.size() ? &_g : nullptr),
				     V(&_V),
				     N(&_N),
				     tol(_tol),
				     binsize((*hists)[0].getbinsize()),
				     lv((*hists)[0].getlv()),
				     dim(binsize.size()),
				     expf(vector<valtype>(V->size(),1.0)),
				     DOS({(*hists)[0]})
{
  if(fseeds.size() == expf.size()) {
    for(uint i = 0; i < fseeds.size(); ++i) { expf[i] = exp(fseeds[i]); }
    cout << "# Seeding WHAM iteration with free energies: " << fseeds << endl;
  }
  //first we also precaculate g-weighted number-of-samples here
  vector<HISTOGRAM> h{*hists};
  vector<NARRAY> Ng{hists->size(), DOS.getdosarr()};
  if(g != nullptr) {
    for(uint i = 0; i < h.size(); ++i) {
      auto& hist = h[i];
      for(map<coordtype, vector<uint> >::const_iterator it = record->begin(); it != record->end(); ++it) {
	const auto& coord = it->first;
        const auto& gtmp = (*g)[i][coord];
	const auto ih = hist.find(coord);
	if(ih != hist.end()) { ih->second /= gtmp; }
	Ng[i][coord] = (*N)[i]/gtmp;
      }
    }
  } else {
    for(uint i = 0; i < h.size(); ++i) {
      auto& hist = h[i];
      for(map<coordtype, vector<uint> >::const_iterator it = record->begin(); it != record->end(); ++it) {
	const auto& coord = it->first;
	Ng[i][coord] = (*N)[i];
      }
    }
  }
  //then we cache all the energy 
  vector<NARRAY> expmH{hists->size(), DOS.getdosarr()};
  NARRAY C{DOS.getdosarr()};
  //Again we only loop through non-zero elements of DOS
  for(map<coordtype, vector<uint> >::const_iterator it = record->begin(); it != record->end(); ++it) {
    const coordtype coord = it->first;
    DOS[coord] = 0.0;
    const vector<uint> histids = it->second;
    const vector<valtype> vals = coord2val(coord);
    for(const auto& k : histids) {
      const typename NARRAY::iterator itc = C.find(coord);
      if(itc == C.end()) { C[coord] = h[k][coord]; } 
      else { itc->second += h[k][coord]; }
    }
    for(uint k = 0; k < expmH.size(); ++k) {
      const auto& exparg = -(*V)[k]->ener(vals);
      if(exparg > MAXEXPARG) {
	throw(General_Exception("exp("+tostr(exparg)+") will overflow"));
      } else if(exparg < MINEXPARG) {
        expmH[k][coord] = 0;
	continue;
      }
      expmH[k][coord] = exp(exparg);
    }
  }

  //perform the WHAM iteration
  ulong count = 0;
  vector<double> newexpf(expf);
  do {
    vector<typename NARRAY::const_iterator> itsNg, itsexpmH;
    for_each(Ng.begin(), Ng.end(), [&itsNg](const NARRAY& narr) { itsNg.emplace_back(narr.begin()); });
    for_each(expmH.begin(), expmH.end(), [&itsexpmH](const NARRAY& narr) { itsexpmH.emplace_back(narr.begin()); });
    typename NARRAY::const_iterator itC = C.begin();

    ++count;
    DOS(itC, itsexpmH, itsNg, expf, newexpf);
    shiftf(newexpf);
  } while(!endit(newexpf,count));
}

template <class PENSEMBLE, class HISTOGRAM, class NARRAY>
bool WHAM<PENSEMBLE,HISTOGRAM,NARRAY>::endit(const vector<valtype>& newexpf, const ulong& count) {
  if(count > MAXIT) {
    cerr << "Max number of iteration " << MAXIT << " is reached. Quit\n";
    printfree();
    exit(-1);
  }
  if(count % 100 == 0) { cout << "#At the " << count << "'th iteration\n"; printfree(); } 
  for(uint i = 0; i < newexpf.size(); ++i) {
    if(fabs(log(newexpf[i])/log(expf[i])-1) > tol) {
      expf = newexpf;
      return false;
    }
  }
  cout << "#WHAM converge to " << tol << " after " << count << " iterations" << endl;
  printfree();
  /*cout <<"#Final Density of state" << endl;
  DOS.print();*/
  return true;
}

template <class PENSEMBLE, class HISTOGRAM, class NARRAY>
void WHAM<PENSEMBLE,HISTOGRAM,NARRAY>::shiftf(vector<valtype>& newexpf) const {
  for(int i = newexpf.size()-1; i >= 0; --i) {
    newexpf[i] /= newexpf[0];
  }
}

template <class PENSEMBLE, class HISTOGRAM, class NARRAY>
void WHAM<PENSEMBLE,HISTOGRAM,NARRAY>::printfree() const {
  printf("#%10s%30s\n","State","DimensionlessFreeEnergy");
  for(uint i = 0; i < expf.size(); ++i) {
    printf("#%10d%30.15lf\n",i,log(expf[i]));
  }
}

template <class PENSEMBLE, class HISTOGRAM, class NARRAY>
NARRAY WHAM<PENSEMBLE,HISTOGRAM,NARRAY>::calrho(const vector<uint>& _dim, const pEnsemble& _V) const {

  if(_dim.size() > dim) {
    cerr << "The number of dimension can't be larger than " << dim << endl;
    exit(-1);
  }

  cout << "#Calculating probability density along dimension: ";
  copy(_dim.begin(),_dim.end(),ostream_iterator<uint>(cout," "));
  cout << endl;

  //Construct the total bin-size in all the dimension 
  /*valtype binsize_sub = 1;
  for(uint i = 0; i < binsize.size(); ++i) {
    binsize_sub *= binsize[i];
    cout << "# " << i << " binsize_sub = " << binsize_sub << endl;
  }*/
  //Check if the _dim actually make sense, i.e., if they are in the range of [0, DOS.dim-1]
  //project out all the dimension irrelevent to the PMF BTW
  for(uint i = 0; i < _dim.size(); ++i) {
    if(_dim[i] >= dim) {
      cerr << "The dimension of probability density must be in the range of [0, " << dim-1  << "]" << endl;;
      exit(-1);
    }
    //const uint dimindex = _dim[i];
    //binsize_sub /= binsize[dimindex];
  }

  NARRAY rho(DOS.getdosarr(),_dim);
  //Again we only loop through non-zero elements of DOS
  for(map<coordtype, vector<uint> >::const_iterator it = record->begin(); it != record->end(); ++it) {
    const coordtype coord = it->first;
    const vector<uint> histids = it->second;
    const vector<valtype> vals = coord2val(coord);
    const valtype exparg = -_V->ener(vals);
    //if(exparg > MAXEXPARG ) { cerr << "In WHAM::calrho: exp("<<exparg<<") will overflow!\n"; exit(-1); }
    //else if(exparg < MINEXPARG) { continue; }
    const valtype weight = exp(exparg);

    //construct a vector<uint> as the coordinate of the PMF
    coordtype rhocoord;
    for(uint i = 0; i < _dim.size(); ++i) {
      rhocoord.push_back(coord[_dim[i]]);
    }
    
    //Update the PMF
    //cout << "#DOS[coord] = " << DOS[coord] << " exparg = " << exparg << " weight = " << weight << " binsize_sub = " << binsize_sub << endl;
    const valtype rhov = DOS[coord]*weight; //There should be a binsize factor here but if we only care about relative probability (density), it cancels out in the ratio
                                            //The reason why I didn't multiply with the binsize factor is that it could get very large in multidimensional cases
    const typename NARRAY::iterator rhoit = rho.find(rhocoord);
    if( rhoit != rho.end()) {
      rho[rhocoord] += rhov;
    }
    else {
      rho[rhocoord] = rhov; 
    }
  }
  return rho;
}


template <class PENSEMBLE, class HISTOGRAM, class NARRAY>
valtype WHAM<PENSEMBLE,HISTOGRAM,NARRAY>::whamvsrawi(const uint& i, const NARRAY& rho) const {
  const HISTOGRAM hist = (*hists)[i];
  const valtype sum = hist.sum();
  const valtype Qinv = expf[i]; //this is the inverse of the partition function 
  const pEnsemble Vi = (*V)[i];
  typename HISTOGRAM::const_iterator it;
  valtype eta = 0; //eta in equation 41 from reference DOI: 10.1002/jcc.21989
  for(it = hist.begin(); it != hist.end(); ++it) {
    //this is the raw normalized HISTOGRAM from i'th simulation
    const valtype praw = it->second/sum;
    //next we calculate the consensus HISTOGRAM from the WHAM distribution
    const coordtype coord = it->first;
    const vector<valtype> vals = hist.coord2val(coord);
    const valtype pwham = Qinv * exp(-Vi->ener(vals)) * rho[coord];
    //then we calculate eta, the relative entropy between praw and pwham
    eta += praw * log(praw/pwham);
  }
  return eta;
}

template <class PENSEMBLE, class HISTOGRAM, class NARRAY>
vector<valtype> WHAM<PENSEMBLE,HISTOGRAM,NARRAY>::whamvsraw(const NARRAY& rho) const {
  const uint K = (*hists).size();
  vector<valtype> etas;
  for(uint i = 0; i < K; ++i) {
    etas.push_back(whamvsrawi(i, rho));
  }
  return etas;
}

template <class PENSEMBLE, class HISTOGRAM, class NARRAY>
const vector<valtype> WHAM<PENSEMBLE,HISTOGRAM,NARRAY>::coord2val(const coordtype& coord) const {
  vector<valtype> vals;
  for(uint i = 0; i < dim; ++i) {
    vals.push_back(lv[i]+(coord[i]+0.5)*binsize[i]);
  }
  return vals;
}

template <class PENSEMBLE, class HISTOGRAM, class NARRAY>
const vector<valtype> WHAM<PENSEMBLE,HISTOGRAM,NARRAY>::coord2val(const vector<uint>& _dim, const coordtype& coord) const {
  vector<valtype> vals;
  for(uint i = 0; i < coord.size(); ++i) {
    const uint dimindex = _dim[i];
    vals.push_back(lv[dimindex]+(coord[dimindex]+0.5)*binsize[dimindex]);
  }
  return vals;
}

#endif
