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
    void print() const; /*  print out all the parameters read */

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
        void print() const; /*  print out all the parameters read */

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
        void print() const; /*  print out all the parameters read */

	/* ====================  MUTATORS      ======================================= */

	/* ====================  OPERATORS     ======================================= */
	//parse the mdp options for pull
        bool operator()(const string& opt, const string& optval) throw(GMXMDP_Exception, FILEIO_Exception);
	/* ====================  DATA MEMBERS  ======================================= */
	enum PullType {Umbrella, Constraint, ConstantForce, Contact, NPullTypes} pullT;
	enum PullGeom {Distance, Direction, DirectionPeriodic, Cylinder, Position, NGeomTypes} geomT;
	//bit mask of the pulled dimension in left-to-right order (e.g., 001 <=> Y N N, meaning only x dimention is pulled)
	short dim;
	int nstx, nstf;
	int npgrps;
	int ncntgrps;
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
            void print() const; /*  print out all the parameters read */

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

	/*
	 * =====================================================================================
	 *        Class:  PULLCNTGRP
	 *  Description:  contact-group parameters
	 * =====================================================================================
	 */
	class PULLCNTGRP
	{
	  public:
	    /* ====================  LIFECYCLE     ======================================= */
	    PULLCNTGRP ();                             /* constructor */

	    /* ====================  ACCESSORS     ======================================= */
            void print() const; /*  print out all the parameters read */

	    /* ====================  MUTATORS      ======================================= */

	    /* ====================  OPERATORS     ======================================= */
	    //parse the mdp options for contact-group
            bool operator()(const string& opt, const string& optval, const string& i) throw(GMXMDP_Exception, FILEIO_Exception);

	  protected:
	    /* ====================  METHODS       ======================================= */

	    /* ====================  DATA MEMBERS  ======================================= */
	    double rate, nc, ncB, k, kB;

	  private:
	    /* ====================  METHODS       ======================================= */

	    /* ====================  DATA MEMBERS  ======================================= */

	}; /* -----  end of class PULLCNTGRP  ----- */
	vector<PULLCNTGRP> pullcntgrps;


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
  fileio fio(fname, fstream::in, 1, 0, 1, MAXNLINE, ";");
  while(fio.readaline()) {
    const vector<string> tmpstr = fio.line2str("=");
    const string& opt = trimltcm(tmpstr[0], ";", " \t");
    const string& optval = trimltcm(tmpstr[1], ";", " \t");
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
void GMXMDP<THamiltonian>::print() const {
  fep.print();
  pull.print();
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
void GMXMDP<THamiltonian>::FEP::print() const {
  printf("#%20s = %-5d\n", "free-energy", fepT);
  printf("#%20s = ", "bonded-lambdas");
  copy(Lbond.begin(), Lbond.end(), ostream_iterator<double>(cout, " "));
  cout << endl;
  printf("#%20s = ", "mass-lambdas");
  copy(Lmass.begin(), Lmass.end(), ostream_iterator<double>(cout, " "));
  cout << endl;
  printf("#%20s = ", "vdw-lambdas");
  copy(Lvdw.begin(), Lvdw.end(), ostream_iterator<double>(cout, " "));
  cout << endl;
  printf("#%20s = ", "coul-lambdas");
  copy(Lcoul.begin(), Lcoul.end(), ostream_iterator<double>(cout, " "));
  cout << endl;
  printf("#%20s = ", "rst-lambdas");
  copy(Lrst.begin(), Lrst.end(), ostream_iterator<double>(cout, " "));
  cout << endl;
  printf("#%20s = ", "cnt1-lambdas");
  copy(Lcnt1.begin(), Lcnt1.end(), ostream_iterator<double>(cout, " "));
  cout << endl;
  printf("#%20s = ", "cnt2-lambdas");
  copy(Lcnt2.begin(), Lcnt2.end(), ostream_iterator<double>(cout, " "));
  cout << endl;
  printf("#%20s = ", "cnt3-lambdas");
  copy(Lcnt3.begin(), Lcnt3.end(), ostream_iterator<double>(cout, " "));
  cout << endl;
}

template < class THamiltonian >
bool GMXMDP<THamiltonian>::FEP::operator()(const string& opt, const string& optval) throw(GMXMDP_Exception, FILEIO_Exception)
{
  if(nocmmatchkey(opt, "free-energy")) {
    if(nocmmatchkey(optval,"yes")) {
	fepT = Yes; 
    } else if(nocmmatchkey(optval, "expanded")) {
	fepT = Expanded;
    } else {
      throw(GMXMDP_Exception("Unknown FEP type: " + optval));
    }
  } else if(nocmmatchkey(opt, "bonded-lambdas")) {
    setvec(optval, Lbond);
  } else if(nocmmatchkey(opt, "mass-lambdas")) {
    setvec(optval, Lmass);
  } else if(nocmmatchkey(opt, "vdw-lambdas")) {
    setvec(optval, Lvdw);
  } else if(nocmmatchkey(opt, "coul-lambdas")) {
    setvec(optval, Lcoul);
  } else if(nocmmatchkey(opt, "restraint-lambdas")) {
    setvec(optval, Lrst);
  } else if(nocmmatchkey(opt, "umbrella-contact1-lambdas")) {
    setvec(optval, Lcnt1);
  } else if(nocmmatchkey(opt, "umbrella-contact2-lambdas")) {
    setvec(optval, Lcnt2);
  } else if(nocmmatchkey(opt, "umbrella-contact3-lambdas")) {
    setvec(optval, Lcnt3);
  } else if(nocmmatchkey(opt, "init-lambda-state")) {
    setsc(optval, Linit);
  } else if(nocmmatchkey(opt, "nstdhdl")) {
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
  ncntgrps(0),
  pullgrps(vector<PULLGRP>(npgrps, PULLGRP())),
  pullcntgrps(vector<PULLCNTGRP>(ncntgrps, PULLCNTGRP()))
{
}

template < class THamiltonian >
void GMXMDP<THamiltonian>::PULL::print() const {
  printf("#%20s = %-5d\n", "Pull", pullT);
  printf("#%20s = %-5d\n", "Pull-geometry", geomT);
  printf("#%20s = %-5d\n", "Pull-dim", dim);
  printf("#%20s = %-5d\n", "Pull-nstxout", nstx);
  printf("#%20s = %-5d\n", "Pull-nstfout", nstf);
  printf("#%20s = %-5d\n", "Pull-ngroups", npgrps);
  printf("#%20s = %-5d\n", "Pull-ncontactgroups", ncntgrps);
  if(pullT == Contact) {
    for(uint i = 0; i < ncntgrps; ++i) {
      printf("#%20s%-5d\n", "Pull-contactgroup", i);
      pullcntgrps[i].print();
    }
  } else {
    for(uint i = 0; i < npgrps; ++i) {
      printf("#%20s%-5d\n", "Pullgroup", i+1);
      pullgrps[i].print();
    }
  }
}

template < class THamiltonian >
bool GMXMDP<THamiltonian>::PULL::operator()(const string& opt, const string& optval) throw(GMXMDP_Exception, FILEIO_Exception)
{
  if(nocmmatchkey(opt, "pull")) {
    if(nocmmatchkey(optval,"umbrella")) {
      pullT = Umbrella; 
    } else if(nocmmatchkey(optval, "constraint")) {
      pullT = Constraint;
    } else if(nocmmatchkey(optval, "constant-force")) {
      pullT = ConstantForce;
    } else if(nocmmatchkey(optval, "contact")) {
      pullT = Contact;
    } else {
      throw(GMXMDP_Exception("Unknown Pull type: '" + optval + "'"));
    }
  } else if(nocmmatchkey(opt, "pull_geometry") || nocmmatchkey(opt, "pull-geometry")) {
    if(nocmmatchkey(optval,"distance")) {
      geomT = Distance; 
    } else if(nocmmatchkey(optval, "direction")) {
      geomT = Direction; 
    } else if(nocmmatchkey(optval, "direction-periodic")) {
      geomT = DirectionPeriodic; 
    } else if(nocmmatchkey(optval, "cylinder")) {
      geomT = Cylinder; 
    } else if(nocmmatchkey(optval, "position")) {
      geomT = Position; 
    } else {
      throw(GMXMDP_Exception("Unknown Pull-geometry type: '" + optval + "'"));
    }
  } else if(nocmmatchkey(opt, "pull_dim") || nocmmatchkey(opt, "pull-dim")) {
    vector<string> dims;
    setvec(optval, dims);
    if(dims.size() != 3) { throw(GMXMDP_Exception("pull dimension must be <= 3 but we got: " + tostr(dims.size()))); }
    for(uint i = 0; i < 3; ++i) {
      if(nocmcmatchkey(dims[i], "Y")) {
	dim |= 1 << i;
      } else if(nocmcmatchkey(dims[i], "N")) {
	const short mask = numeric_limits<short>::max() ^ (1 << i) ;
	dim &= mask;
      } else {
        throw(GMXMDP_Exception("Unknown Pull-dim: " + optval));
      }
    }
  } else if(nocmmatchkey(opt, "pull_nstxout") || nocmmatchkey(opt, "pull-nstxout")) {
    setsc(optval, nstx);
  } else if(nocmmatchkey(opt, "pull_start") || nocmmatchkey(opt, "pull-start")) {
    return true; //not really need this but just to avoid additional if-else
  } else if(nocmmatchkey(opt, "pull_nstfout") || nocmmatchkey(opt, "pull-nstfout")) {
    setsc(optval, nstf);
  } else if(nocmmatchkey(opt, "pull_ngroups") || nocmmatchkey(opt, "pull-ngroups")) {
    setsc(optval, npgrps);
    pullgrps.resize(npgrps);
  } else if(nocmmatchkey(opt, "pull_ncontactgroups") || nocmmatchkey(opt, "pull-ncontactgroups")) {
    setsc(optval, ncntgrps);
    pullcntgrps.resize(ncntgrps);
  } else if(nocmmatchkey(opt, "pull_group0") || nocmmatchkey(opt, "pull-group0")) {
    return true;
  } else {
    if(pullT == Contact) {
      for(uint i = 0; i < ncntgrps; ++i) {
	if(pullcntgrps[i](opt, optval, tostr(i))) { return true; }
      }
    } else {
      for(uint i = 0; i < npgrps; ++i) {
        if(pullgrps[i](opt, optval, tostr(i+1))) { return true; }
      } 
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
void GMXMDP<THamiltonian>::PULL::PULLGRP::print() const {
  printf("#  %20s = ", "pull_vec");
  copy(vec.begin(), vec.end(), ostream_iterator<double>(cout, " "));
  cout << endl;
  printf("#  %20s = ", "pull_init");
  copy(init.begin(), init.end(), ostream_iterator<double>(cout, " "));
  cout << endl;
  printf("#  %20s = ", "pull_initB");
  copy(initB.begin(), initB.end(), ostream_iterator<double>(cout, " "));
  cout << endl;
  printf("#  %20s = %-10.5lf\n", "pull_rate", rate);
  printf("#  %20s = %-10.5lf\n", "pull_k", k);
  printf("#  %20s = %-10.5lf\n", "pull_kB", kB);
}

template < class THamiltonian >
bool GMXMDP<THamiltonian>::PULL::PULLGRP::operator()(const string& opt, const string& optval, const string& i) throw(GMXMDP_Exception, FILEIO_Exception)
{
  if(nocmmatchkey(opt, "pull_group" + i) || nocmmatchkey(opt, "pull-group" + i)) {
    return true; //not really need this
  } else if(nocmmatchkey(opt, "pull_vec" + i) || nocmmatchkey(opt, "pull-vec" + i)) {
    setvec(optval, vec);
  } else if(nocmmatchkey(opt, "pull_init" + i) || nocmmatchkey(opt, "pull-init" + i)) {
    setvec(optval, init);
  } else if(nocmmatchkey(opt, "pull_initB" + i) || nocmmatchkey(opt, "pull-initB" + i)) {
    setvec(optval, initB);
  } else if(nocmmatchkey(opt, "pull_rate" + i) || nocmmatchkey(opt, "pull-rate" + i)) {
    setsc(optval, rate);
  } else if(nocmmatchkey(opt, "pull_k" + i) || nocmmatchkey(opt, "pull-k" + i)) {
    setsc(optval, k);
  } else if(nocmmatchkey(opt, "pull_kB" + i) || nocmmatchkey(opt, "pull-kB" + i)) {
    setsc(optval, kB);
  } else {
    return false;
  }
  return true;
}

template < class THamiltonian >
GMXMDP<THamiltonian>::PULL::PULLCNTGRP::PULLCNTGRP() :
  rate(0.0),
  nc(0.0),
  ncB(0.0),
  k(0.0),
  kB(k)
{
}

template < class THamiltonian >
void GMXMDP<THamiltonian>::PULL::PULLCNTGRP::print() const {
  printf("#  %20s = %-10.5lf\n", "rate", rate);
  printf("#  %20s = %-10.5lf\n", "nc", nc);
  printf("#  %20s = %-10.5lf\n", "ncB", ncB);
  printf("#  %20s = %-10.5lf\n", "k", k);
  printf("#  %20s = %-10.5lf\n", "kB", kB);
}

template < class THamiltonian >
bool GMXMDP<THamiltonian>::PULL::PULLCNTGRP::operator()(const string& opt, const string& optval, const string& i) throw(GMXMDP_Exception, FILEIO_Exception)
{
  const string cntgrpnm = "pull_contactgroup" + i;
  if(nocmmatchkey(opt, cntgrpnm + "_pgrp1")) {
    return true; //not really need this
  } else if(nocmmatchkey(opt, cntgrpnm + "_pgrp2")) {
    return true; //not really need this
  } else if(nocmmatchkey(opt, cntgrpnm + "_ncrate")) {
    setsc(optval, rate);
  } else if(nocmmatchkey(opt, cntgrpnm + "_k")) {
    setsc(optval, k);
  } else if(nocmmatchkey(opt, cntgrpnm + "_kB")) {
    setsc(optval, kB);
  } else if(nocmmatchkey(opt, cntgrpnm + "_nc")) {
    setsc(optval, nc);
  } else if(nocmmatchkey(opt, cntgrpnm + "_ncB")) {
    setsc(optval, ncB);
  } else {
    return false;
  }
  return true;
}
