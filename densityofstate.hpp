#include "hamiltonian.hpp"
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
    //! update density of state given a set of inputs
    /**
     * @param[in] record The 1st vector<uint> is a set of coordinates in at least one of the histograms, whose values is non-zero; and the 2nd vector<uint> keep the index of the histograms
     * @param[in] hists Generic histogram for each trajectory
     * @param[in] g Statistical inefficiency for each trajectory
     * @param[in] V Generic hamiltonian that combine the conserved quantities and the associated parameters. Total number of elements is the total number of states considered
     * @param[in] N N[k][l] is the number of samples the k'th trajectory visiting the l'th state. 
     * @param[in] f Free energy of states. Total number of elements is the total number of states considered 
     */
    void operator() (
    		     const map<coordtype, vector<uint> >& record, 
		     					    
                     const vector<histogram>& hists,  
                     const vector<narray>& g,  
		     const vector<Hamiltonian<ensemble>* >& V, 
		     const vector<vector<uint> >& N, 
		     const vector<valtype>& f,
		     vector<valtype>& newf 
		    );
    //! same as the one above except we assume all elements in g[m][k] are the same for all k given m
    void operator() (
    		     const map<coordtype, vector<uint> >& record, 
                     const vector<histogram>& hists,  
		     const vector<Hamiltonian<ensemble>* >& V, 
		     const vector<vector<uint> >& N, 
		     const vector<valtype>& f,
		     vector<valtype>& newf 
		    ); 
    //! same as the one above except we further assume each trajectory only visit 1 state
    /** This means that V.size() == hists.size() 
     * @param[in] N N[k] is the number of samples in the k'th trajectory/state. 
     */
    void operator() (
    		     const map<coordtype, vector<uint> >& record, 
                     const vector<histogram>& hists,  
		     const vector<Hamiltonian<ensemble>* >& V, 
		     const vector<uint>& N, 
		     const vector<valtype>& f,
		     vector<valtype>& newf
		    );
    //! same as the one above except we don't calculate newf but 
    //instead calculate the 2nd term in equation 2.22 and 2.23 in doc and 
    //save that in F and gradF, respectively
    /** 
     * @param[in] N N[k] is the number of samples in the k'th trajectory/state.
     * @param[out] F is the contribution of the DOS to the WHAM likelihood function
     * @param[out] gradF is the contribution of the DOS to the gradient of the WHAM likelihood function
     */
    void operator() (
    		     const map<coordtype, vector<uint> >& record, 
                     const vector<histogram>& hists,  
		     const vector<Hamiltonian<ensemble>* >& V, 
		     const vector<uint>& N, 
		     const vector<valtype>& q,
		     valtype& F,
		     vector<valtype>& gradF
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
    		     	      const map<coordtype, vector<uint> >& record, 
                              const vector<histogram>& hists,  
                              const vector<narray>& g, 
         		      const vector<Hamiltonian<ensemble>* >& V, 
         		      const vector<vector<uint> >& N, 
         		      const vector<valtype>& f,
			      vector<valtype>& newf
         		                   ) 
{
  //First empty newf array
  newf = vector<valtype>(f.size(),0.0);
  vector<valtype> expnewf(f.size(),0.0);
  //expenergy[k] == exp(f[k]-V[k]->ener(vals))
  vector<valtype> expenergy(f.size(),0.0);
  //Here we loop through all the non-zero histogram values and compute the density of state from the histograms
  for(map<coordtype, vector<uint> >::const_iterator it = record.begin(); it != record.end(); ++it) {
    valtype num = 0.0, denum = 0.0;
    const coordtype coord = it->first;
    const vector<uint> histids = it->second;
    const vector<valtype> vals = hists[0].coord2val(coord);
    for(uint k = 0; k < histids.size(); ++k) {
      const uint histid = histids[k];
      num += hists[histid][coord]/g[histid][coord];  
    }
    for(uint k = 0; k < hists.size(); ++k) {
      valtype denum_part = 0.0;
      for(uint l = 0; l < V.size(); ++l) { //TODO: Can't we move this loop outside the k-loop?
	const valtype exparg = f[l]-V[l]->ener(vals);
	if(exparg > MAXEXPARG) { cerr << "exp(exparg) will overflow!\n"; exit(-1); }
	else if(exparg < MINEXPARG) { continue; }
	expenergy[l] = exp(exparg);
	denum_part += N[k][l]*expenergy[l]; //There should be a bin-size term here but it cancels out with the same term in f[l] 
							 //since we assume all bin-sizes are the same across different histograms
      }
      denum += denum_part/g[k][coord];
    }
    const valtype dos = num/denum;
    DOS[coord] = dos;
    for(uint l = 0; l < newf.size(); ++l) {
      expnewf[l] += dos*expenergy[l]; 
    }
  }
  for(uint l = 0; l < expnewf.size(); ++l) {
    newf[l] = -log(expnewf[l]/exp(f[l])); //exp(f[l]) needs to be taken out because expenergy[l] = exp(f[l]-V[l]->ener(vals))
  }
}

template<class ensemble, class histogram, class narray>
void DensityOfState<ensemble,histogram,narray>::operator () (
    		     	      const map<coordtype, vector<uint> >& record, 
                              const vector<histogram>& hists,  
         		      const vector<Hamiltonian<ensemble>* >& V, 
         		      const vector<vector<uint> >& N, 
         		      const vector<valtype>& f,
			      vector<valtype>& newf
         		                   ) 
{
  //First empty newf array
  newf = vector<valtype>(f.size(),0.0);
  vector<valtype> expnewf(f.size(),0.0);
  //expenergy[k] == exp(f[k]-V[k]->ener(vals))
  vector<valtype> expenergy(f.size(),0.0);
  //Here we loop through all the non-zero histogram values and compute the density of state from the histograms
  for(map<coordtype, vector<uint> >::const_iterator it = record.begin(); it != record.end(); ++it) {
    valtype num = 0.0, denum = 0.0;
    const coordtype coord = it->first;
    const vector<uint> histids = it->second;
    const vector<valtype> vals = hists[0].coord2val(coord);
    for(uint k = 0; k < histids.size(); ++k) {
      const uint histid = histids[k];
      num += hists[histid][coord];  
    }
    for(uint k = 0; k < hists.size(); ++k) {
      valtype denum_part = 0.0;
      for(uint l = 0; l < V.size(); ++l) {
	const valtype exparg = f[l]-V[l]->ener(vals);
	if(exparg > MAXEXPARG) { cerr << "exp("<<exparg<<") will overflow!\n"; exit(-1); }
	else if(exparg < MINEXPARG) { continue; }
	expenergy[l] = exp(exparg);
	denum_part += N[k][l]*expenergy[l]; //There should be a bin-size term here but it cancels out with the same term in f[l] 
							 //since we assume all bin-sizes are the same across different histograms
      }
      denum += denum_part;
    }
    const valtype dos = num/denum;
    DOS[coord] = dos;
    for(uint l = 0; l < newf.size(); ++l) {
      expnewf[l] += dos*expenergy[l]; 
    }
  }
  for(uint l = 0; l < expnewf.size(); ++l) {
    newf[l] = -log(expnewf[l]/exp(f[l])); //exp(f[l]) needs to be taken out because expenergy[l] = exp(f[l]-V[l]->ener(vals))
  }
}

template<class ensemble, class histogram, class narray>
void DensityOfState<ensemble,histogram,narray>::operator () (
    		     	      const map<coordtype, vector<uint> >& record, 
                              const vector<histogram>& hists,  
         		      const vector<Hamiltonian<ensemble>* >& V, 
         		      const vector<uint>& N, 
         		      const vector<valtype>& f,
			      vector<valtype>& newf
         		                   ) 
{
  //First empty newf array
  newf = vector<valtype>(f.size(),0.0);
  vector<valtype> expnewf(f.size(),0.0);
  //expenergy[k] == exp(f[k]-V[k]->ener(vals))
  vector<valtype> expenergy(f.size(),0.0);
  //Here we loop through all the non-zero histogram values and compute the density of state from the histograms
  //printf("#RC k hist_k N_k f_k exparg exp(exparg)\n");
  for(map<coordtype, vector<uint> >::const_iterator it = record.begin(); it != record.end(); ++it) {
    valtype num = 0.0, denum = 0.0;
    const coordtype coord = it->first;
    const vector<uint> histids = it->second;
    const vector<valtype> vals = hists[0].coord2val(coord);
    for(uint k = 0; k < histids.size(); ++k) {
      const uint histid = histids[k];
      num += hists[histid][coord];  
    }
    /*cout << "# ";
    copy(coord.begin(),coord.end(),ostream_iterator<uint>(cout," "));
    cout << num << endl;*/
    for(uint k = 0; k < hists.size(); ++k) {
      /*cout <<"#DOS: state ensemble_params vals ener" << endl;
      cout << k << " ";
      const vector<valtype> params = V[k]->getens().getparams();
      copy(params.begin(),params.end(),ostream_iterator<valtype>(cout," "));
      copy(vals.begin(),vals.end(),ostream_iterator<valtype>(cout," "));
      cout << V[k]->ener(vals) << endl;*/
      valtype denum_part = 0.0;
      const valtype exparg = f[k]-V[k]->ener(vals);
      if(exparg > MAXEXPARG ) { cerr << "exp("<<exparg<<") will overflow!\n"; exit(-1); }
      else if(exparg < MINEXPARG) { continue; }
      expenergy[k] = exp(exparg);
      denum_part += N[k]*expenergy[k]; //There should be a bin-size term here but it cancels out with the same term in f[k] 
                                      //since we assume all bin-sizes are the same across different histograms
      denum += denum_part;

      /*const ulong hist_k = hists[k].find(coord) == hists[k].end() ? 0 : hists[k][coord];
      printf("#%30.15lf%5u%10lu%10u%30.15lf%30.15lf%30.15lf\n",vals[0],k,hist_k,N[k],f[k],exparg,exp(exparg));*/
    }
    const valtype dos = num/denum;
    DOS[coord] = dos;
    //printf("#num = %30.15lf denum = %30.15lf DOS = %30.15lf\n",num,denum,DOS[coord]);
    for(uint l = 0; l < expnewf.size(); ++l) {
      expnewf[l] += dos*expenergy[l]; 
    }
  }
  for(uint l = 0; l < expnewf.size(); ++l) {
    newf[l] = -log(expnewf[l]/exp(f[l])); //exp(f[l]) needs to be taken out because expenergy[l] = exp(f[l]-V[l]->ener(vals))
  }
}

template<class ensemble, class histogram, class narray>
void DensityOfState<ensemble,histogram,narray>::operator () (
    		     	      const map<coordtype, vector<uint> >& record, 
                              const vector<histogram>& hists,  
         		      const vector<Hamiltonian<ensemble>* >& V, 
         		      const vector<uint>& N, 
         		      const vector<valtype>& q,
			      valtype& F,
			      vector<valtype>& gradF
         		                   ) 
{
  const uint G = gradF.size();
  //expenergy[k] == exp(f[k]-V[k]->ener(vals))
  vector<valtype> expenergy(V.size(),0.0);
  //Here we loop through all the non-zero histogram values and compute the density of state from the histograms
  //printf("#RC k hist_k N_k f_k exparg exp(exparg)\n");
  for(map<coordtype, vector<uint> >::const_iterator it = record.begin(); it != record.end(); ++it) {
    valtype num = 0.0, denum = 0.0;
    const coordtype coord = it->first;
    const vector<uint> histids = it->second;
    const vector<valtype> vals = hists[0].coord2val(coord);
    for(uint k = 0; k < histids.size(); ++k) {
      const uint histid = histids[k];
      num += hists[histid][coord];  
    }
    //The contribution of DOS to gradF[i] is a weighted sum of density of state
    //and the weight is just the i'th denum_part 
    vector<valtype> weights(G,0);
    for(uint k = 0; k < hists.size(); ++k) {
      valtype denum_part = 0.0;
      const valtype exparg = -V[k]->ener(vals);
      if(exparg > MAXEXPARG ) { 
	cerr << "exp("<<exparg<<") will overflow!\n"; 
	exit(-1); 
      } else if(exparg < MINEXPARG) { continue; }
      expenergy[k] = exp(exparg);
      denum_part += q[k]*N[k]*expenergy[k]; //There should be a bin-size term here but it cancels out with the same term in f[k] 
                                      //since we assume all bin-sizes are the same across different histograms
      denum += denum_part;
      //update weights for gradF
      weights[k] = denum_part;
    }
    const valtype dos = num/denum;
    DOS[coord] = dos;
    //update F
    F -= num*log(dos);
    //update gradF
    for ( uint i = 0; i < G; ++i) {
      gradF[i] += dos*weights[i];
    }
  }
}
