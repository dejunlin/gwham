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
    //! update density of state given a set of inputs
    /**
     * @param[in] record The 1st vector<uint> is a set of coordinates in at least one of the histograms, whose values is non-zero; and the 2nd vector<uint> keep the index of the histograms
     * @param[in] hists Generic histogram for each trajectory
     * @param[in] g[k][coord] Statistical inefficiency of the k'th trajectory in bin coord
     * @param[in] V Generic hamiltonian that combine the conserved quantities and the associated parameters. Total number of elements is the total number of states considered
     * @param[in] N N[k] is the number of samples in the k'th trajectory/state. 
     * @param[in] f Free energy of states. Total number of elements is the total number of states considered 
     */
    void operator() (
    		     const map<coordtype, vector<uint> >& record, 
                     const vector<histogram>& hists,  
		     const vpEnsemble& V, 
		     const vector<linecounter>& N, 
		     const vector<valtype>& f,
		     vector<valtype>& newf,
                     const vector<narray>& g 
		    );
    //! update density of state given a set of inputs (assumes g[k][coord] are invariant across all k for any coord)
    /**
     * @param[in] record The 1st vector<uint> is a set of coordinates in at least one of the histograms, whose values is non-zero; and the 2nd vector<uint> keep the index of the histograms
     * @param[in] hists Generic histogram for each trajectory
     * @param[in] V Generic hamiltonian that combine the conserved quantities and the associated parameters. Total number of elements is the total number of states considered
     * @param[in] N N[k] is the number of samples in the k'th trajectory/state. 
     * @param[in] f Free energy of states. Total number of elements is the total number of states considered 
     */
    void operator() (
    		     const map<coordtype, vector<uint> >& record, 
                     const vector<histogram>& hists,  
		     const vpEnsemble& V, 
		     const vector<linecounter>& N, 
		     const vector<valtype>& f,
		     vector<valtype>& newf
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
         		      const vpEnsemble& V, 
         		      const vector<linecounter>& N, 
         		      const vector<valtype>& f,
			      vector<valtype>& newf,
                              const vector<narray>& g
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
      num += hists[histid][coord]/g[histid][coord];  
    }
    /*cout << "# ";
    copy(coord.begin(),coord.end(),ostream_iterator<uint>(cout," "));
    cout << num << endl;*/
    for(uint k = 0; k < hists.size(); ++k) {
      /*cout <<"#DOS: state ensemble_params vals ener" << endl;
      cout << k << " ";
      const vector<valtype> params = V[k]->getens()->getparams();
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
      denum += denum_part/g[k][coord];

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
         		      const vpEnsemble& V, 
         		      const vector<linecounter>& N, 
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
      const vector<valtype> params = V[k]->getens()->getparams();
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
