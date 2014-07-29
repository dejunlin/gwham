/*
 * =====================================================================================
 *
 *       Filename:  gmxmdp.hpp
 *
 *    Description:  GROMACS mdp file reader and Hamiltonian factory
 *
 *        Version:  1.0
 *        Created:  03/07/14 16:33:51
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dejun Lin (DL), dejun.lin@gmail.com
 *   Organization:  Department of Biochemistry and Biophysics, Medical Center, University of Rochester
 *
 * =====================================================================================
 */

#include "typedefs.hpp"
#include "fileio.hpp"
#include "fileio_utils.hpp"
#include "exception.hpp"
#include <string>
#include <string.h>

using namespace std;

/*
 * =====================================================================================
 *        Class:  GMXMDP
 *  Description:  
 * =====================================================================================
 */
template < class THamiltonian >
class GMXMDP
{
  public:
    /* ====================  LIFECYCLE     ======================================= */
    GMXMDP (const string& fname);    /* constructor */

    /* ====================  ACCESSORS     ======================================= */

    /* ====================  MUTATORS      ======================================= */

    /* ====================  OPERATORS     ======================================= */
    
    THamiltonian genH() const;

  protected:
    
    /*
     * =====================================================================================
     *        Class:  FEP
     *  Description:  Free Energy Perturbation parameters
     * =====================================================================================
     */
    class FEP
    {
      public:
	/* ====================  LIFECYCLE     ======================================= */
	FEP () ;                             /* constructor */

	/* ====================  ACCESSORS     ======================================= */

	/* ====================  MUTATORS      ======================================= */

	/* ====================  OPERATORS     ======================================= */
	//parse the mdp options for FEP
        bool operator()(const string& opt, const string& optval) throw(GMXMDP_Exception, FILEIO_Exception);
	/* ====================  DATA MEMBERS  ======================================= */
	vector<double> Lbond, Lmass, Lvdw, Lcoul, Lrst, Lcnt1, Lcnt2, Lcnt3;
	int nstdhdl;
	int Linit;
	enum FEPType {Yes, Expanded, NFEPTypes} fepT;

      protected:
	/* ====================  METHODS       ======================================= */

	/* ====================  DATA MEMBERS  ======================================= */

      private:
	/* ====================  METHODS       ======================================= */

	/* ====================  DATA MEMBERS  ======================================= */

    } fep; /* -----  end of class FEP  ----- */

    
    /*
     * =====================================================================================
     *        Class:  PULL
     *  Description:  Pulling parameters
     * =====================================================================================
     */
    class PULL
    {
      public:
	/* ====================  LIFECYCLE     ======================================= */
	PULL ();                             /* constructor */

	/* ====================  ACCESSORS     ======================================= */

	/* ====================  MUTATORS      ======================================= */

	/* ====================  OPERATORS     ======================================= */
	//parse the mdp options for pull
        bool operator()(const string& opt, const string& optval) throw(GMXMDP_Exception, FILEIO_Exception);
	/* ====================  DATA MEMBERS  ======================================= */
	enum PullType {Umbrella, Constraint, ConstantForce, NPullTypes} pullT;
	enum PullGeom {Distance, Direction, DirectionPeriodic, Cylinder, Position, NGeomTypes} geomT;
	short dim;
	int nstx, nstf;
	int npgrps;
	/*
	 * =====================================================================================
	 *        Class:  PULLGRP
	 *  Description:  pull-group parameters
	 * =====================================================================================
	 */
	class PULLGRP
	{
	  public:
	    /* ====================  LIFECYCLE     ======================================= */
	    PULLGRP ();                             /* constructor */

	    /* ====================  ACCESSORS     ======================================= */

	    /* ====================  MUTATORS      ======================================= */

	    /* ====================  OPERATORS     ======================================= */
	    //parse the mdp options for pull-group
            bool operator()(const string& opt, const string& optval, const string& i) throw(GMXMDP_Exception, FILEIO_Exception);

	  protected:
	    /* ====================  METHODS       ======================================= */

	    /* ====================  DATA MEMBERS  ======================================= */
	    vector<double> vec, init, initB;
	    double rate, k, kB;

	  private:
	    /* ====================  METHODS       ======================================= */

	    /* ====================  DATA MEMBERS  ======================================= */

	}; /* -----  end of class PULLGRP  ----- */
	vector<PULLGRP> pullgrps;


      protected:
	/* ====================  METHODS       ======================================= */

	/* ====================  DATA MEMBERS  ======================================= */

      private:
	/* ====================  METHODS       ======================================= */

	/* ====================  DATA MEMBERS  ======================================= */

    } pull; /* -----  end of class PULL  ----- */


    /* ====================  METHODS       ======================================= */

    /* ====================  DATA MEMBERS  ======================================= */

  private:
    /* ====================  METHODS       ======================================= */

    /* ====================  DATA MEMBERS  ======================================= */

}; /* ----------  end of template class GMXMDP  ---------- */

template < class THamiltonian >
GMXMDP<THamiltonian>::GMXMDP(const string& fname):
  fep(FEP()),
  pull(PULL())
{
  fileio fio(fname, fstream::in, 1, 0, 1, MAXNLINE);
  while(fio.readaline()) {
    cout << "Reading line " << fio.line << endl;
    const vector<string> tmpstr = fio.line2str("=");
    if(!tmpstr.size()) { continue; } //skip empty lines
    const string& opt = tmpstr[0];
    const string& optval = tmpstr[1];
    cout << opt << " ===== " << optval << endl;
    try {
      const bool p = fep(opt, optval) || pull(opt, optval);
    } catch(GMXMDP_Exception& gmxmdpex) {
      cerr << "Error understanding mdp file: " << fname << ": " << gmxmdpex.what() << endl;
      terminate();
    } catch(FILEIO_Exception& fioex) {
      cerr << "Error reading mdp file: " << fname << ": " << fioex.what() << endl;
      terminate();
    }
  }
}

template < class THamiltonian >
GMXMDP<THamiltonian>::FEP::FEP() :
   Lbond(vector<double>(0, 0.0)),
   Lmass(vector<double>(0, 0.0)),
   Lvdw(vector<double>(0, 0.0)),
   Lcoul(vector<double>(0, 0.0)),
   Lrst(vector<double>(0, 0.0)),
   Lcnt1(vector<double>(0, 0.0)),
   Lcnt2(vector<double>(0, 0.0)),
   Lcnt3(vector<double>(0, 0.0)),
   nstdhdl(-1),
   Linit(-1),
   fepT(NFEPTypes)
{
}

template < class THamiltonian >
bool GMXMDP<THamiltonian>::FEP::operator()(const string& opt, const string& optval) throw(GMXMDP_Exception, FILEIO_Exception)
{
  if(matchkey(opt, "free-energy")) {
    if(matchkey(optval,"yes")) {
	fepT = Yes; 
    } else if(matchkey(optval, "expanded")) {
	fepT = Expanded;
    } else {
      throw(GMXMDP_Exception("Unknown FEP type: " + optval));
    }
  } else if(matchkey(opt, "bonded-lambdas")) {
    setvec(optval, Lbond);
  } else if(matchkey(opt, "mass-lambdas")) {
    setvec(optval, Lmass);
  } else if(matchkey(opt, "vdw-lambdas")) {
    setvec(optval, Lvdw);
  } else if(matchkey(opt, "coul-lambdas")) {
    setvec(optval, Lcoul);
  } else if(matchkey(opt, "restraint-lambdas")) {
    setvec(optval, Lrst);
  } else if(matchkey(opt, "umbrella-contact1-lambdas")) {
    setvec(optval, Lcnt1);
  } else if(matchkey(opt, "umbrella-contact2-lambdas")) {
    setvec(optval, Lcnt2);
  } else if(matchkey(opt, "umbrella-contact3-lambdas")) {
    setvec(optval, Lcnt3);
  } else if(matchkey(opt, "init-lambda-state")) {
    setsc(optval, Linit);
  } else if(matchkey(opt, "nstdhdl")) {
    setsc(optval, nstdhdl);
  } else { 
    return false;
  }
  return true;
}

template < class THamiltonian >
GMXMDP<THamiltonian>::PULL::PULL() :
  pullT(NPullTypes),
  geomT(NGeomTypes),
  dim(0),
  nstx(-1),
  nstf(-1),
  npgrps(0),
  pullgrps(vector<PULLGRP>(npgrps, PULLGRP()))
{
}

template < class THamiltonian >
bool GMXMDP<THamiltonian>::PULL::operator()(const string& opt, const string& optval) throw(GMXMDP_Exception, FILEIO_Exception)
{
  if(matchkey(opt, "pull")) {
    if(matchkey(optval,"umbrella")) {
      pullT = Umbrella; 
    } else if(matchkey(optval, "constraint")) {
      pullT = Constraint;
    } else if(matchkey(optval, "constant-force")) {
      pullT = ConstantForce;
    } else {
      throw(GMXMDP_Exception("Unknown Pull type: " + optval));
    }
  } else if(matchkey(opt, "pull_geometry") || matchkey(opt, "pull-geometry")) {
    if(matchkey(optval,"distance")) {
      geomT = Distance; 
    } else if(matchkey(optval, "direction")) {
      geomT = Direction; 
    } else if(matchkey(optval, "direction-periodic")) {
      geomT = DirectionPeriodic; 
    } else if(matchkey(optval, "cylinder")) {
      geomT = Cylinder; 
    } else if(matchkey(optval, "position")) {
      geomT = Position; 
    } else {
      throw(GMXMDP_Exception("Unknown Pull-geometry type: " + optval));
    }
  } else if(matchkey(opt, "pull_dim") || matchkey(opt, "pull-dim")) {
    vector<string> dims;
    setvec(optval, dims);
    if(dims.size() != 3) { throw(GMXMDP_Exception("pull dimension must be <= 3 but we got: " + tostr(dims.size()))); }
    for(uint i = 0; i < 3; ++i) {
      if(cmatchkey(dims[i], "Y")) {
	dim |= 1 << i;
      } else if(cmatchkey(dims[i], "N")) {
	const short mask = numeric_limits<short>::max() ^ (1 << i) ;
	dim |= mask;
      } else {
        throw(GMXMDP_Exception("Unknown Pull-dim: " + optval));
      }
    }
  } else if(matchkey(opt, "pull_nstxout") || matchkey(opt, "pull-nstxout")) {
    setsc(optval, nstx);
  } else if(matchkey(opt, "pull_start") || matchkey(opt, "pull-start")) {
    return true; //not really need this but just to avoid additional if-else
  } else if(matchkey(opt, "pull_nstfout") || matchkey(opt, "pull-nstfout")) {
    setsc(optval, nstf);
  } else if(matchkey(opt, "pull_ngroups") || matchkey(opt, "pull-ngroups")) {
    setsc(optval, npgrps);
    pullgrps.resize(npgrps);
  } else if(matchkey(opt, "pull_group0") || matchkey(opt, "pull-group0")) {
    return true;
  } else {
    for(uint i = 0; i < npgrps; ++i) {
      if(pullgrps[i](opt, optval, tostr(i+1))) { return true; }
    }
    return false;
  }
  return true;
}

template < class THamiltonian >
GMXMDP<THamiltonian>::PULL::PULLGRP::PULLGRP() :
  vec(vector<double>(0, 0.0)),
  init(vector<double>(0, 0.0)),
  initB(vector<double>(0, 0.0)),
  rate(0.0),
  k(0.0),
  kB(k)
{
}

template < class THamiltonian >
bool GMXMDP<THamiltonian>::PULL::PULLGRP::operator()(const string& opt, const string& optval, const string& i) throw(GMXMDP_Exception, FILEIO_Exception)
{
  if(matchkey(opt, "pull_group" + i) || matchkey(opt, "pull-group" + i)) {
    return true; //not really need this
  } else if(matchkey(opt, "pull_vec" + i) || matchkey(opt, "pull-vec" + i)) {
    setvec(optval, vec);
  } else if(matchkey(opt, "pull_init" + i) || matchkey(opt, "pull-init" + i)) {
    setvec(optval, init);
  } else if(matchkey(opt, "pull_initB" + i) || matchkey(opt, "pull-initB" + i)) {
    setvec(optval, initB);
  } else if(matchkey(opt, "pull_rate" + i) || matchkey(opt, "pull-rate" + i)) {
    setsc(optval, rate);
  } else if(matchkey(opt, "pull_k" + i) || matchkey(opt, "pull-k" + i)) {
    setsc(optval, k);
  } else if(matchkey(opt, "pull_kB" + i) || matchkey(opt, "pull-kB" + i)) {
    setsc(optval, kB);
  } else {
    return false;
  }
  return true;
}
