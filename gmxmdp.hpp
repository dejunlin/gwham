#if !defined(GMXMDP_HPP)
#define GMXMDP_HPP
/*
 * =====================================================================================
 *
 *       Filename:  gmxmdp.hpp
 *
 *    Description:  GROMACS mdp file reader 
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
#include "mdp.hpp"
#include "functor.hpp"
#include <string>
#include <string.h>

using namespace std;

/*
 * =====================================================================================
 *        Class:  GMXMDP
 *  Description:  
 * =====================================================================================
 */
class GMXMDP : public MDP
{
  protected:
    /*
     * =====================================================================================
     *        Class:  GMXGENERIC
     *  Description:  
     * =====================================================================================
     */
    class GMXGENERIC : public GENERIC
    {
      public:
	/* ====================  LIFECYCLE     ======================================= */
	GMXGENERIC () {};                             /* constructor */

	/* ====================  ACCESSORS     ======================================= */
        virtual void print() const; /*  print out all the parameters read */

	/* ====================  MUTATORS      ======================================= */
	virtual void doublechk() throw(MDP_Exception); /*  check and combine parameters */

	/* ====================  OPERATORS     ======================================= */
        //parse the generic parameters
        virtual bool operator()(const string& opt, const string& optval) throw(MDP_Exception, FILEIO_Exception);

      protected:
	/* ====================  METHODS       ======================================= */

	/* ====================  DATA MEMBERS  ======================================= */

      private:
	/* ====================  METHODS       ======================================= */

	/* ====================  DATA MEMBERS  ======================================= */

    }; /* -----  end of class GMXGENERIC  ----- */

    
    /*
     * =====================================================================================
     *        Class:  FEP
     *  Description:  Free Energy Perturbation parameters
     * =====================================================================================
     */
    class GMXFEP : public FEP
    {
      public:
	/* ====================  LIFECYCLE     ======================================= */
	GMXFEP () ;                             /* constructor */

	/* ====================  ACCESSORS     ======================================= */
        virtual void print() const; /*  print out all the parameters read */
	virtual vector<valtype> getAllLambdas() const;
	virtual vector<valtype> getInitLambdas() const;

	/* ====================  MUTATORS      ======================================= */
	virtual void doublechk() throw(MDP_Exception); /*  check and combine parameters */

	/* ====================  OPERATORS     ======================================= */
	//parse the mdp options for FEP
        virtual bool operator()(const string& opt, const string& optval) throw(MDP_Exception, FILEIO_Exception);
	/* ====================  DATA MEMBERS  ======================================= */
	vector<valtype> Lcnt1, Lcnt2, Lcnt3;
	int nstdhdl, nstexpanded;
	valtype Tmc;

      protected:
	/* ====================  METHODS       ======================================= */

	/* ====================  DATA MEMBERS  ======================================= */

      private:
	/* ====================  METHODS       ======================================= */

	/* ====================  DATA MEMBERS  ======================================= */

    }; /* -----  end of class FEP  ----- */

    
    /*
     * =====================================================================================
     *        Class:  PULL
     *  Description:  Pulling parameters
     * =====================================================================================
     */
    class GMXPULL : public PULL
    {
      public:
	enum PullType {Umbrella, UmbrellaFlatBottom, Constraint, ConstantForce, Contact, NPullTypes} pullT;
	enum PullGeom {Distance, Direction, DirectionPeriodic, Cylinder, Position, NGeomTypes} geomT;
      public:
	/* ====================  LIFECYCLE     ======================================= */
	GMXPULL ();                             /* constructor */
	virtual ~GMXPULL() { };

	/* ====================  ACCESSORS     ======================================= */
        virtual void print() const; /*  print out all the parameters read */

	/* ====================  MUTATORS      ======================================= */
	virtual void doublechk() throw(MDP_Exception); /*  check and combine parameters */

	/* ====================  OPERATORS     ======================================= */
	//parse the mdp options for pull
        virtual bool operator()(const string& opt, const string& optval) throw(MDP_Exception, FILEIO_Exception);
	/* ====================  DATA MEMBERS  ======================================= */
	int nstx, nstf; /* output frequency in steps for RC and forces*/
	/*
	 * =====================================================================================
	 *        Class:  PULLGRP
	 *  Description:  regular gromacs pull-group parameters
	 * =====================================================================================
	 */
	class GMXPULLGRP : public PULLGRP
	{
	  public:
	    /* ====================  LIFECYCLE     ======================================= */
	    GMXPULLGRP ();                             /* constructor */
	    GMXPULLGRP (const vector<valtype>& _vec, const vector<valtype>& _init, const vector<valtype>& _initB, const valtype& _k, const valtype& _kB, const valtype& _r0,const valtype& _r0B,const valtype& _r1,const valtype& _r1B,const valtype& _k0,const valtype& _k0B,const valtype& _k1,const valtype& _k1B,const vector<valtype>& _Lrst);                             /* constructor */
	    virtual ~GMXPULLGRP() {  };

	    /* ====================  ACCESSORS     ======================================= */
            virtual void print() const; /*  print out all the parameters read */

	    /* ====================  MUTATORS      ======================================= */
	    virtual void doublechk() throw(MDP_Exception); /*  check and combine parameters */

	    /* ====================  OPERATORS     ======================================= */
            virtual bool operator()(const string& opt, const string& optval, const string& i) throw(MDP_Exception, FILEIO_Exception); /* parse the mdp options for pull-group */

	    vector<valtype> vec, init, initB;
	    valtype k, kB, r0, r0B, r1, r1B, k0, k0B, k1, k1B;
	    vector<valtype> Lrst; /* restraint-lambdas */
	  protected:
	    /* ====================  METHODS       ======================================= */

	    /* ====================  DATA MEMBERS  ======================================= */

	  private:
	    /* ====================  METHODS       ======================================= */

	    /* ====================  DATA MEMBERS  ======================================= */

	}; /* -----  end of class PULLGRP  ----- */

	/*
	 * =====================================================================================
	 *        Class:  GMXPULLCNTGRP
	 *  Description:  gromacs contact-group parameters (my own hack-up)
	 * =====================================================================================
	 */
	class GMXPULLCNTGRP : public PULLGRP
	{
	  public:
	    /* ====================  LIFECYCLE     ======================================= */
	    GMXPULLCNTGRP ();                             /* constructor */
	    virtual ~GMXPULLCNTGRP() { };

	    /* ====================  ACCESSORS     ======================================= */
            virtual void print() const; /*  print out all the parameters read */

	    /* ====================  MUTATORS      ======================================= */
	    virtual void doublechk() throw(MDP_Exception); /*  check and combine parameters */


	    /* ====================  OPERATORS     ======================================= */
	    //parse the mdp options for contact-group
            virtual bool operator()(const string& opt, const string& optval, const string& i) throw(MDP_Exception, FILEIO_Exception);

	    valtype nc, ncB, k, kB;
	    vector<valtype> Lcnt; /* contact-lambdas */
	  protected:
	    /* ====================  METHODS       ======================================= */

	    /* ====================  DATA MEMBERS  ======================================= */

	  private:
	    /* ====================  METHODS       ======================================= */

	    /* ====================  DATA MEMBERS  ======================================= */

	}; /* -----  end of class GMXPULLCNTGRP  ----- */

      protected:
	/* ====================  METHODS       ======================================= */

	/* ====================  DATA MEMBERS  ======================================= */

      private:
	/* ====================  METHODS       ======================================= */

	/* ====================  DATA MEMBERS  ======================================= */

    }; /* -----  end of class PULL  ----- */
    typedef GMXPULL::GMXPULLGRP GMXPULLGRP;
    typedef GMXPULL::GMXPULLCNTGRP GMXPULLCNTGRP;
  public:
    /* ====================  LIFECYCLE     ======================================= */
    GMXMDP (const string& fname);    /* constructor */
    virtual ~GMXMDP() { 
    };

    /* ====================  ACCESSORS     ======================================= */
    virtual void print() const; /*  print out all the parameters read */
    virtual uint cmp(const MDP& mdp) const throw (MDP_Exception);

    /* ====================  MUTATORS      ======================================= */
    virtual void doublechk() throw(MDP_Exception);
    void setLPull(); 

    /* ====================  OPERATORS     ======================================= */
    

    /* ====================  METHODS       ======================================= */

    /* ====================  DATA MEMBERS  ======================================= */

  private:
    /* ====================  METHODS       ======================================= */

    /* ====================  DATA MEMBERS  ======================================= */

}; /* ----------  end of template class GMXMDP  ---------- */

#endif
