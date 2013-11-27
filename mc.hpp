#include "typedefs.hpp"
#include "gnarray.hpp"
#include "ran2.c"
#include <vector>

typedef gnarray<coordtype, valtype, valtype> narray;
typedef gnarray<coordtype, histcounter, valtype> histogram;

class potpoly_dblwell {
  public:
    //!Constructor
    potpoly_dblwell(const uint& _dim, const vector<valtype>& _params) : 
      dim(_dim),
      params(_params)
    {
    };
    //!Calculate potential energy given coordinates
    valtype operator() (const vector<valtype>& coord) const {
      valtype ans = 0.0;
      for(uint i = 0; i < coord.size(); ++i) {
	const valtype coord2 = coord[i]*coord[i];
	ans += params[0]*coord2*(params[1]*coord2 + params[2]);
      }
      return ans;
    };
    //!access dim
    uint getdim() const { return dim;};
    //!Calculate partition function governed by this potential
    valtype partitionfunct(const valtype& kB, const valtype& T, const vector<uint>& Npts, const vector<valtype>& hv, const vector<valtype>& lv) const {

      vector<uint> stride_array(dim+1,0);
      setstride(Npts,stride_array);

      const valtype beta = 1/(kB*T);

      vector<valtype> dx(dim,0.0);
      ulong Npts_tot = 1;
      valtype dx_tot = 1.0;
      for(uint i = 0; i < dim; ++i) {
	dx[i] = (hv[i] - lv[i])/Npts[i];
	Npts_tot *= Npts[i];
	dx_tot *= dx[i];
      }

      valtype intgl = 0.0;
      for(ulong i = 0; i < Npts_tot; ++i) {
	vector<valtype> coord;
	for(uint j = 0; j < dim; ++j) {
	  const uint index = int( ( i%stride_array[j] ) / stride_array[j+1] );
	  coord.push_back(lv[j]+(index+0.5)*dx[j]);
	}
	const valtype fx = this->operator()(coord);
	intgl += exp(-fx*beta)*dx_tot;
      }
      return intgl;
    };
  private:
    const uint dim;
    const vector<valtype> params;
    //!set-up stride_array
    void setstride(const vector<uint>& Npts, vector<uint>& stride_array) const {
      if(!stride_array.size()) {
	stride_array = vector<uint>(dim+1,0);
      }
      stride_array[dim] = 1; //equivalent to stride(dim-1)
      if(!dim) { return; }
      for(int i = dim-1; i >= 0; --i) {
        stride_array[i] = Npts[i]*stride_array[i+1];
      }
      return;
    };
};

//! This is just an empty bias potential class for debugging
/*class nobias {
  public:
    nobias(const uint& _dim) : dim(_dim) {};
    valtype operator() (const vector<valtype>& coord) const { return 0.0;}
    uint getdim() const { return dim; }
  private:
    const uint dim;
};*/

class harmonic {
  public:
    harmonic(const uint& _dim, const vector<valtype>& _k, const vector<valtype>& _r) :
      dim(_dim),
      k(_k),
      r(_r)
    {
      if(dim != k.size() || dim != r.size()) {
	cerr << "Number of dimension of input parameters for harmonic doesn't match what's specified\n";
	exit(-1);
      }
    }
    valtype operator() (const vector<valtype>& coord) const {
      valtype ans = 0.0;
      for(uint i = 0; i < coord.size(); ++i) {
	const valtype dev = coord[i] - r[i];
	ans += 0.5*k[i]*dev*dev;
      }
      return ans;
    }
    valtype ener(const vector<valtype>& coord) const {
      return this->operator()(coord);
    }
    uint getdim() const { return dim;};
  private:
    const uint dim;
    const vector<valtype> k;
    const vector<valtype> r;
};

template<class potener, class bias>
class MC {
  public:
    MC(const valtype& _kB, const valtype& _T, const potener& _U, const bias& _V, const ulong& _nsteps, const valtype& _stepsize) :
      kB(_kB),
      T(_T),
      beta(1/(kB*T)),
      U(_U),
      V(_V),
      nsteps(_nsteps),
      stepsize(_stepsize)
    {
      if(U.getdim() != V.getdim()) { cerr << "Mismatched dimension in potential and bias\n"; exit(-1); }
    }
    //! propagate the state and generate histogram
    /** Return the acceptance ratio
     */
    valtype operator() (const vector <valtype>& init_coord, histogram& hist) const {
      //Initialize random seed
      //srand48(time(NULL));
      srandom(time(NULL));
      long idum = -random();
      ran2(&idum);

      const valtype stephalf = stepsize/2;
      vector<valtype> coord  = init_coord;
      const uint dim = coord.size();
      ulong naccept = 0;

      for(ulong t = 0; t <= nsteps; ++t) {
	vector<valtype> newcoord = coord;
	for(uint i = 0; i < dim; ++i) {
	  newcoord[i] += ran2(&idum)*stepsize - stephalf;
	}
	const valtype E = U(coord) + V(coord);
	/*cout << "#coord = ";
	copy(coord.begin(),coord.end(),ostream_iterator<valtype>(cout," "));
	cout << endl;
	cout << "U(coord) = " << U(coord) << " Bias(coord) = " << V(coord) << endl;*/
	const valtype newE = U(newcoord) + V(newcoord);
	const valtype dE = newE - E;
	if(dE < 0) { 
	  coord = newcoord;
	  ++naccept;
	} else if( exp(-dE*beta) > ran2(&idum)) {
	    coord = newcoord;
	    ++naccept;
	}
	hist.bin(coord);
      }
      return naccept/nsteps;
    };
  private:
    const valtype kB;
    const valtype T; 
    const valtype beta;
    const potener U;
    const bias V;
    const ulong nsteps;
    const valtype stepsize;
};
