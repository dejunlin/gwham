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
  public:
    /* ====================  LIFECYCLE     ======================================= */
    GMXMDP (const string& fname);    /* constructor */

    /* ====================  ACCESSORS     ======================================= */
    virtual void print() const; /*  print out all the parameters read */
    virtual Qt cmp(const MDP& mdp) const; /*  compare with another mdp object */
    virtual bool hasTemperature() const; /*  if has temperature */
    virtual bool hasPressure() const; /*  if has pressure */
    virtual bool hasLambda() const; /*  if has Lambdas */
    virtual bool hasRestraint() const; /*  if has restraint */
    virtual valtype getTemperature() const; /*  get temperature */
    virtual valtype getPressure() const; /*  get pressure */
    virtual vector<valtype> getFEPLambda() const; /*  get FEPLambdas */
    virtual vector<Functor> getRestraint() const; /*  get restraint */

    /* ====================  MUTATORS      ======================================= */

    /* ====================  OPERATORS     ======================================= */
    
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
	GMXGENERIC ();                             /* constructor */

	/* ====================  ACCESSORS     ======================================= */
        virtual void print() const; /*  print out all the parameters read */
        virtual Qt cmp(const GENERIC& mdp) const; /*  compare with another mdp object */
        virtual valtype getT() const { return T; }; 
        virtual valtype getP() const { return P; }; 

	/* ====================  MUTATORS      ======================================= */

	/* ====================  OPERATORS     ======================================= */
        //parse the generic parameters
        virtual bool operator()(const string& opt, const string& optval) throw(GMXMDP_Exception, FILEIO_Exception);

      protected:
	/* ====================  METHODS       ======================================= */

	/* ====================  DATA MEMBERS  ======================================= */

      private:
	/* ====================  METHODS       ======================================= */

	/* ====================  DATA MEMBERS  ======================================= */

    } generic; /* -----  end of class GMXGENERIC  ----- */

    
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

	/* ====================  MUTATORS      ======================================= */

	/* ====================  OPERATORS     ======================================= */
	//parse the mdp options for FEP
        virtual bool operator()(const string& opt, const string& optval) throw(GMXMDP_Exception, FILEIO_Exception);
	/* ====================  DATA MEMBERS  ======================================= */
	vector<valtype> Lcnt1, Lcnt2, Lcnt3;
	int nstdhdl, nstexpanded;
	int Linit;
	enum FEPType {Yes, Expanded, NFEPTypes} fepT;
	valtype Tmc;

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
    class GMXPULL : public PULL
    {
      public:
	/* ====================  LIFECYCLE     ======================================= */
	GMXPULL ();                             /* constructor */

	/* ====================  ACCESSORS     ======================================= */
        virtual void print() const; /*  print out all the parameters read */

	/* ====================  MUTATORS      ======================================= */

	/* ====================  OPERATORS     ======================================= */
	//parse the mdp options for pull
        virtual bool operator()(const string& opt, const string& optval) throw(GMXMDP_Exception, FILEIO_Exception);
	/* ====================  DATA MEMBERS  ======================================= */
	enum PullType {Umbrella, Constraint, ConstantForce, Contact, NPullTypes} pullT;
	enum PullGeom {Distance, Direction, DirectionPeriodic, Cylinder, Position, NGeomTypes} geomT;
	int nstx, nstf;
	uint ncntgrps;
	/*
	 * =====================================================================================
	 *        Class:  PULLGRP
	 *  Description:  pull-group parameters
	 * =====================================================================================
	 */
	class GMXPULLGRP : public PULLGRP
	{
	  public:
	    /* ====================  LIFECYCLE     ======================================= */
	    GMXPULLGRP ();                             /* constructor */

	    /* ====================  ACCESSORS     ======================================= */
            virtual void print() const; /*  print out all the parameters read */

	    /* ====================  MUTATORS      ======================================= */

	    /* ====================  OPERATORS     ======================================= */
	    //parse the mdp options for pull-group
            virtual bool operator()(const string& opt, const string& optval, const string& i) throw(GMXMDP_Exception, FILEIO_Exception);

	  protected:
	    /* ====================  METHODS       ======================================= */

	    /* ====================  DATA MEMBERS  ======================================= */
	    vector<valtype> vec;

	  private:
	    /* ====================  METHODS       ======================================= */

	    /* ====================  DATA MEMBERS  ======================================= */

	}; /* -----  end of class PULLGRP  ----- */
	vector<GMXPULLGRP> pullgrps;

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
	    valtype rate, nc, ncB, k, kB;

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

#endif
