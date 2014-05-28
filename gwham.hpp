#if !defined(GWHAM_HPP)
#define GWHAM_HPP

#include "densityofstate.hpp"
#include "hamiltonian.hpp"
#include <stdlib.h>
#include <iterator>
#include <vector>

using namespace std;

//histogram class must provide public member function 'coord2val' that maps a coordinate to corresponding value
//histogram class must have member variable 'binsize' which is the size of bin along each dimension 
//narray class must have square bracket operator that maps a index (coordinate) to an element
//narray class must have member iterator type
//narray class must have member function find(coordtype) which given a coordinate return the iterator it points at
//narray class must have member function end() which return the past-the-end iterator
template <class ensemble, class histogram, class narray>
class WHAM {
  typedef DensityOfState<ensemble,histogram,narray> DOStype;
  typedef DOStype* DOSptr;
  public:
    //!perform WHAM iteration and update the density of state
    /**
     * @param[in] record record[i][k] is the index of the k'th histograms that has non-zero value at point whose coordinate is i 
     * @param[in] hists Generic histogram for each trajectory
     * @param[in] g Statistical inefficiency for each trajectory
     * @param[in] V[i] is the hamiltonian the i'th state that combine the conserved quantities and the associated parameters
     * @param[in] N N[k][l] is the number of samples the k'th trajectory visiting the l'th state. 
     * @param[in] tol Tolerance for WHAM iteration 
     */
    WHAM(const map<coordtype, vector<uint> >& _record, 
         const vector<histogram>& hists,  
         const vector<narray>& g,  
	 const vector<Hamiltonian<ensemble>* >& V, 
         const vector<vector<uint> >& N, 
	 const valtype _tol 
	);
    //! same as the one above except we assume all elements in g[m][k] are the same for all k given m
    WHAM(const map<coordtype, vector<uint> >& _record, 
         const vector<histogram>& hists,  
	 const vector<Hamiltonian<ensemble>* >& V, 
         const vector<vector<uint> >& N, 
	 const valtype _tol 
	);
    //! same as the one above except we further assume each trajectory only visit 1 state
    /** This means that V.size() == hists.size() 
     * @param[in] N N[k] is the number of samples in the k'th trajectory/state. 
     */
    WHAM(const map<coordtype, vector<uint> >& _record, 
         const vector<histogram>& hists,  
	 const vector<Hamiltonian<ensemble>* >& V, 
         const vector<uint>& N, 
	 const valtype _tol
	);
    //! Based on the density of state, calculate a new set of free energies
    void calnewf(vector<valtype>& f, const vector<Hamiltonian<ensemble>* >& V) const;
    //!calculate PMF in Hamiltonian V along the dimension in DOS as specified by dim
    narray calrho(const vector<uint>& dim, const Hamiltonian<ensemble>& V) const; 
    //!given coordinate in index space, calculate the corresponding real-space coordinate
    const vector<valtype> coord2val(const coordtype& coord) const;
    //!given coordinate in index space along the specified dimension, calculate the corresponding real-space coordinate
    const vector<valtype> coord2val(const vector<uint>& _dim, const coordtype& coord) const;
    //!shift all the free energies (argument newf)by a constant so that newf[0] is always zero 
    void shiftf(vector<valtype>& newf) const;
    //!print out all free energies
    void printfree() const;
  private:
     /*
      * =====================================================================================
      *        Class:  LogLikeFunct
      *  Description:  The log-likelihood function of WHAM
      * =====================================================================================
      */
     class LogLikeFunct
     {
       public:
	 /* ====================  LIFECYCLE     ======================================= */
	 LogLikeFunct (WHAM<ensemble,histogram,narray>* const _wham);                          /* constructor */

	 /* ====================  ACCESSORS     ======================================= */
	 //! call this->(df) and then return F; this is awkward but it's how dlib works:
	 //dlib::find_min() always call getF() and getgradF() in pair such that getF() 
	 //is called first and then is getgradF() -- This is really dumb since in general
	 //you can save a lot of computation by calculating the function and its derivative
	 //together
	 valtype getF(const vector<valtype>& df);
	 //! just return gradF; NOTE that this is expected to be called right after getF() 
	 const vector<valtype> getgradF(const vector<valtype>& df) const;

	 /* ====================  MUTATORS      ======================================= */

	 /* ====================  OPERATORS     ======================================= */

	 /* 
	  * ===  FUNCTION  ======================================================================
	  *         Name:  operator()
	  *  Description:   
	  * ====================================================================================
	  * @param[in] df df[i] is dg_i as in equation 2.20 in the documentation 
          */
	  void 
	   operator() ( 
                      const vector<valtype>& df
	              );

       protected:
	 /* ====================  METHODS       ======================================= */

	 /* ====================  DATA MEMBERS  ======================================= */

       private:
	 /* ====================  METHODS       ======================================= */

	 /* ====================  DATA MEMBERS  ======================================= */
	 //! pointer to the outter object that contains *this 
	 WHAM<ensemble,histogram,narray>* const wham;
	 //! value of the likelihood function, which is updated when LogLikeFunct::operator() is called
	 valtype F;
	 vector<valtype> gradF;

     }; /* -----  end of class LogLikeFunct  ----- */
     friend class LogLikeFunct;

     //!max number of iterations
     static const ulong MAXIT = 1000000;
     //!check if WHAM iteration should end; also update WHAM::f if not end
     bool endit(const vector<valtype>& newf, const ulong& count); 
     //!record[i][k] is the index of the k'th histograms that has non-zero value at point whose coordinate is i
     const map<coordtype, vector<uint> >* const record;
     //!Generic histogram for each trajectory
     const vector<histogram>* const hists;  
     //!Statistical inefficiency for each trajectory in each bin
     const vector<narray>* const g;
     //!V[i] is the hamiltonian functor the i'th state that compute the energy given the values corresponding to a histogram bin m 
     const vector<Hamiltonian<ensemble>* >* const V; 
     //!N[k][l] is the number of samples the k'th trajectory visiting the l'th state
     const vector<vector<uint> >* const N;
     //!N1[k] is the number of samples from the k'th trajectory
     const vector<uint>* const N1;
     //!Tolerance for WHAM iteration
     const valtype tol; 
     //!f[i] is the dimensionless free energy of state i 
     //!bin-size of each dimension of the histogram
     const vector<valtype> binsize;
     //!lower bound in real space of the histogram
     const vector<valtype> lv;
     //!number of dimension of the histogram
     const uint dim;
     //!dimensionless free energy of each state
     vector<valtype> f;
     //!density of state
     DOStype DOS;
     //!the likelihood function to be optimized
     LogLikeFunct funct;
};

template <class ensemble, class histogram, class narray>
WHAM<ensemble,histogram,narray>::LogLikeFunct::LogLikeFunct
(WHAM<ensemble,histogram,narray>* const _wham) :
  wham(_wham),
  F(0),
  gradF(vector<valtype>(wham->f.size()-1, 0.0))
  {};

template <class ensemble, class histogram, class narray>
void WHAM<ensemble,histogram,narray>::LogLikeFunct::operator()
(
  const vector<valtype>& df
) {
  const vector<histogram>& hists = *(wham->hists);  
  const vector<Hamiltonian<ensemble>* >& V = *(wham->V); 
  const vector<uint>& N = *(wham->N1);

  const uint G = df.size();
  const uint K = N.size();
  //TODO: We should check if K == G + 1 here but ignore this for now
  //since we expect WHAM constructor compute df from WHAM::f initially

  // First term in equation 2.22 and 2.23 in doc
  valtype F_part1 = 0.0;
  vector<uint> sN(K,0); //sN[i] = N[i+1] + N[i+2] + ... + N[K-1] is the coefficient of df[i]
  for (  int i = K-2; i >= 0; --i ) {
    sN[i] = sN[i+1] + N[i+1];
  }
  for ( uint i = 0; i < G ; ++i ) {
    F_part1 -= sN[i] * df[i];
    gradF[i] -= sN[i];
  }

  // Second term in equation 2.22 and 2.23 in doc, which is just sum over m: 
  // C_m*ln(C_m*Omega_m) and Omega_m is computed in DOS
  // we need to update WHAM::f from df first
  vector<valtype>& f = wham->f;
  f[0] = 0.0;
  for ( uint i = 0; i < G; ++i) {
    f[i+1] = f[i] + df[i];
  }
  valtype F_part2 = 0.0;
  wham->DOS(wham->record,hists,V,N,f,F_part2,gradF);
  F = F_part1 + F_part2;
}

template <class ensemble, class histogram, class narray>
valtype WHAM<ensemble,histogram,narray>::LogLikeFunct::getF(const vector<valtype>& df) { 
  this->operator()(df);
  return F; 
}

template <class ensemble, class histogram, class narray>
const vector<valtype> WHAM<ensemble,histogram,narray>::LogLikeFunct::getgradF(const vector<valtype>& df) const { return gradF; }

template <class ensemble, class histogram, class narray>
WHAM<ensemble,histogram,narray>::WHAM(const map<coordtype, vector<uint> >& _record, 
                                      const vector<histogram>& _hists,  
                                      const vector<narray>& _g,  
	                              const vector<Hamiltonian<ensemble>* >& _V, 
                                      const vector<vector<uint> >& _N, 
	                              const double _tol 
	                             ):
				     record(&_record),
				     hists(&_hists),
				     g(&_g),
				     V(&_V),
				     N(&_N),
				     N1(NULL),
				     tol(_tol),
				     binsize((*hists)[0].getbinsize()),
				     lv((*hists)[0].getlv()),
				     dim(binsize.size()),
				     f(vector<valtype>(V->size(),0.0)),
				     DOS(DOStype((*hists)[0])),
				     funct(this)
{
  //perform the WHAM iteration
  ulong count = 0;
  vector<double> newf(f);
  do {
    ++count;
    DOS(*record, *hists, *g, *V, *N, f, newf);
    //calnewf(newf,V);
    shiftf(newf);
  } while(!endit(newf,count)); 
}

template <class ensemble, class histogram, class narray>
WHAM<ensemble,histogram,narray>::WHAM(const map<coordtype, vector<uint> >& _record, 
                                      const vector<histogram>& _hists,  
	                              const vector<Hamiltonian<ensemble>* >& _V, 
                                      const vector<vector<uint> >& _N, 
	                              const double _tol 
	                             ):
				     record(&_record),
				     hists(&_hists),
				     g(NULL),
				     V(&_V),
				     N(&_N),
				     N1(NULL),
				     tol(_tol),
				     binsize((*hists)[0].getbinsize()),
				     lv((*hists)[0].getlv()),
				     dim(binsize.size()),
				     f(vector<valtype>(V->size(),0.0)),
				     DOS(DOStype((*hists)[0])),
				     funct(this)
{
  //perform the WHAM iteration
  ulong count = 0;
  vector<double> newf(f);
  do {
    ++count;
    DOS(*record, *hists, *V, *N, f, newf);
    //calnewf(newf,V);
    shiftf(newf);
  } while(!endit(newf,count)); 
}

template <class ensemble, class histogram, class narray>
WHAM<ensemble,histogram,narray>::WHAM(const map<coordtype, vector<uint> >& _record, 
                                      const vector<histogram>& _hists,  
	                              const vector<Hamiltonian<ensemble>* >& _V, 
                                      const vector<uint>& _N, 
	                              const double _tol 
	                             ):
				     record(&_record),
				     hists(&_hists),
				     g(NULL),
				     V(&_V),
				     N(NULL),
				     N1(&_N),
				     tol(_tol),
				     binsize((*hists)[0].getbinsize()),
				     lv((*hists)[0].getlv()),
				     dim(binsize.size()),
				     f(vector<valtype>(V->size(),0.0)),
				     DOS(DOStype((*hists)[0])),
				     funct(this)
{
  //First use a minimizer to get the optimal estimate of WHAM::f
  //TODO: just use dlib L-BFGS here -- we might eventually implement our own

  //perform the WHAM iteration
  ulong count = 0;
  vector<double> newf(f);
  do {
    ++count;
    DOS(*record, *hists, *V, *N1, f, newf);
    //calnewf(newf,V);
    shiftf(newf);
  } while(!endit(newf,count));
}

template <class ensemble, class histogram, class narray>
void WHAM<ensemble,histogram,narray>::calnewf(vector<valtype>& newf, const vector<Hamiltonian<ensemble>* >& _V) const {
  for(uint l = 0; l < newf.size(); ++l) {
    double expmf = 0.0;
    /*cout <<"#calnewf: state ensemble_params" << endl;
    cout << "#" << l << " ";
    const vector<valtype> params = _V[l]->getens().getparams();
    copy(params.begin(),params.end(),ostream_iterator<valtype>(cout," ")); cout << endl;
    cout <<"#calnewf: vals ener" << endl;*/

    /*cout <<"#state RC energy DOS expmf\n";
    cout <<"#be " << l << endl;*/
    //Here we only loop through non-zero elements of DOS
    for(map<coordtype, vector<uint> >::const_iterator it = record->begin(); it != record->end(); ++it) {
      const coordtype coord = it->first;
      const vector<uint> histids = it->second;
      const vector<valtype> vals = coord2val(coord);
      //calculate f[l] for all l
      const valtype exparg = -_V[l]->ener(vals);
      /*cout << "#";
      copy(vals.begin(),vals.end(),ostream_iterator<valtype>(cout," "));
      cout << " " << exparg << endl;*/
      //if(exparg > MAXEXPARG ) { cerr << "In WHAM::calnewf: exp("<<exparg<<") will overflow!\n"; exit(-1); }
      //else if(exparg < MINEXPARG) { continue; }
      expmf += DOS[coord]*exp(exparg); //There should be a bin-size term here but it cancels out with the same term the denumerator of DOS
      //cout << "#" << l << " " << vals[0] << " " << -exparg << " " << DOS[coord] << " " << expmf << endl;
    }
    //cout << "#ed " << newf[l] << endl;
    newf[l] = -log(expmf);
  }
  shiftf(newf);
}

template <class ensemble, class histogram, class narray>
bool WHAM<ensemble,histogram,narray>::endit(const vector<valtype>& newf, const ulong& count) {
  if(count > MAXIT) {
    cerr << "Max number of iteration " << MAXIT << " is reached. Quit\n";
    printfree();
    exit(-1);
  }
  if(count % 100 == 0) { cout << "#At the " << count << "'th iteration\n"; printfree(); } 
  for(uint i = 0; i < newf.size(); ++i) {
    if(fabs((newf[i] - f[i])/f[i]) > tol) {
      f = newf;
      return false;
    }
  }
  cout << "#WHAM converge to " << tol << " after " << count << " iterations" << endl;
  printfree();
  /*cout <<"#Final Density of state" << endl;
  DOS.print();*/
  return true;
}

template <class ensemble, class histogram, class narray>
void WHAM<ensemble,histogram,narray>::shiftf(vector<valtype>& newf) const {
  /*valtype favg = 0.0;
  valtype fmax = -9999999999, fmin = 9999999999;
  for(uint i = 0; i < newf.size(); ++i) {
    favg += newf[i];
    if(newf[i] > fmax) { fmax = newf[i];}
    else if(newf[i] < fmin) { fmin = newf[i]; }
  }
  favg /= newf.size();
  cout << "#Shift all f by " << favg << " with fmax = " << fmax << " and fmin = " << fmin << endl;
  for(uint i = 0; i < newf.size(); ++i) {
    newf[i] -= favg;
  }*/
  for(int i = newf.size()-1; i >= 0; --i) {
    newf[i] -= newf[0];
  }
}

template <class ensemble, class histogram, class narray>
void WHAM<ensemble,histogram,narray>::printfree() const {
  printf("#%10s%30s\n","State","DimensionlessFreeEnergy");
  for(uint i = 0; i < f.size(); ++i) {
    printf("#%10d%30.15lf\n",i,f[i]);
  }
}

template <class ensemble, class histogram, class narray>
narray WHAM<ensemble,histogram,narray>::calrho(const vector<uint>& _dim, const Hamiltonian<ensemble>& _V) const {

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

  narray rho(DOS.getdosarr(),_dim);
  //Again we only loop through non-zero elements of DOS
  for(map<coordtype, vector<uint> >::const_iterator it = record->begin(); it != record->end(); ++it) {
    const coordtype coord = it->first;
    const vector<uint> histids = it->second;
    const vector<valtype> vals = coord2val(coord);
    const valtype exparg = -_V.ener(vals);
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
    const typename narray::iterator rhoit = rho.find(rhocoord);
    if( rhoit != rho.end()) {
      rho[rhocoord] += rhov;
    }
    else {
      rho[rhocoord] = rhov; 
    }
  }
  return rho;
}

template <class ensemble, class histogram, class narray>
const vector<valtype> WHAM<ensemble,histogram,narray>::coord2val(const coordtype& coord) const {
  vector<valtype> vals;
  for(uint i = 0; i < dim; ++i) {
    vals.push_back(lv[i]+(coord[i]+0.5)*binsize[i]);
  }
  return vals;
}

template <class ensemble, class histogram, class narray>
const vector<valtype> WHAM<ensemble,histogram,narray>::coord2val(const vector<uint>& _dim, const coordtype& coord) const {
  vector<valtype> vals;
  for(uint i = 0; i < coord.size(); ++i) {
    const uint dimindex = _dim[i];
    vals.push_back(lv[dimindex]+(coord[dimindex]+0.5)*binsize[dimindex]);
  }
  return vals;
}

#endif
