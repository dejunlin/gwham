#if !defined(MDP_HPP)
#define MDP_HPP
/*
 * =====================================================================================
 *
 *       Filename:  mdp.hpp
 *
 *    Description:  Interfacial classes for MD-parameters reader
 *
 *        Version:  1.0
 *        Created:  31/07/14 09:36:54
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dejun Lin (DL), dejun.lin@gmail.com
 *   Organization:  Department of Biochemistry and Biophysics, Medical Center, University of Rochester
 *
 * =====================================================================================
 */
#include <vector>
#include "typedefs.hpp"
#include "functor.hpp"

using namespace std;

/*
 * =====================================================================================
 *        Class:  MDP
 *  Description:  Interfacial class for MD-parameters reader
 * =====================================================================================
 */

//The thermodynamic quantities we care about
enum Qt {Temperature, Pressure, Lambdas, Restraints};

class MDP
{
  public:
    /* ====================  LIFECYCLE     ======================================= */
    MDP () {};                             /* constructor */

    /* ====================  ACCESSORS     ======================================= */
    virtual void print() const = 0; /*  print out all the parameters read */
    virtual uint cmp(const MDP& mdp) const = 0; /*  compare with another mdp object */
    virtual bool hasTemperature() const = 0; /*  if has temperature */
    virtual bool hasPressure() const = 0; /*  if has pressure */
    virtual bool hasLambda() const = 0; /*  if has Lambdas */
    virtual bool hasRestraint() const = 0; /*  if has restraint */
    virtual valtype getTemperature() const = 0; /*  get temperature */
    virtual valtype getPressure() const = 0; /*  get pressure */
    virtual vector<valtype> getFEPLambda() const = 0; /*  get FEPLambdas */
    virtual string getRestraintType() const = 0; /*  get the type of restraint */
    virtual vector<valtype> getRestraint() const = 0; /*  get restraint */

    /* ====================  MUTATORS      ======================================= */

    /* ====================  OPERATORS     ======================================= */

  protected:
    /* ====================  METHODS       ======================================= */

    /* ====================  DATA MEMBERS  ======================================= */
    
    /*
     * =====================================================================================
     *        Class:  GENERIC
     *  Description:  Generic parameters
     * =====================================================================================
     */
    class GENERIC
    {
      public:
        /* ====================  LIFECYCLE     ======================================= */
        GENERIC () {};                             /* constructor */
    
        /* ====================  ACCESSORS     ======================================= */
        virtual void print() const = 0; /*  print out all the parameters read */
        virtual valtype getT() const { return T; }; 
        virtual valtype getP() const { return P; }; 
    
        /* ====================  MUTATORS      ======================================= */
    
        /* ====================  OPERATORS     ======================================= */
        //parse the generic parameters
        virtual bool operator()(const string& opt, const string& optval) throw(GMXMDP_Exception, FILEIO_Exception) = 0;
    
      protected:
        /* ====================  METHODS       ======================================= */
    
        /* ====================  DATA MEMBERS  ======================================= */
        //reference temperature and pressure NOTE these parameters might not be relevant
        //to the WHAM calculation since some FEP parameters might overwrite these
        valtype T, P;
    
      private:
        /* ====================  METHODS       ======================================= */
    
        /* ====================  DATA MEMBERS  ======================================= */
    }; /* -----  end of class GENERIC  ----- */

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
	FEP () {};                             /* constructor */

	/* ====================  ACCESSORS     ======================================= */
        void print() const = 0; /*  print out all the parameters read */
	virtual vector<valtype> getL() const {
	  vector<valtype> allL(Lbond.begin(), Lbond.end());
	  allL.insert(allL.end(), Lmass.begin(), Lmass.end());
	  allL.insert(allL.end(), Lvdw.begin(), Lvdw.end());
	  allL.insert(allL.end(), Lcoul.begin(), Lcoul.end());
	  allL.insert(allL.end(), Lrst.begin(), Lrst.end());
	  allL.insert(allL.end(), Ltemp.begin(), Ltemp.end());
	  return allL;
	}

	/* ====================  MUTATORS      ======================================= */

	/* ====================  OPERATORS     ======================================= */
	//parse the mdp options for FEP
        virtual bool operator()(const string& opt, const string& optval) throw(GMXMDP_Exception, FILEIO_Exception) = 0;
	/* ====================  DATA MEMBERS  ======================================= */
	vector<valtype> Lbond, Lmass, Lvdw, Lcoul, Lrst, Ltemp;
	int Linit;

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
    class PULL
    {
      public:
	/* ====================  LIFECYCLE     ======================================= */
	PULL () {};                             /* constructor */

	/* ====================  ACCESSORS     ======================================= */
        void print() const = 0; /*  print out all the parameters read */

	/* ====================  MUTATORS      ======================================= */

	/* ====================  OPERATORS     ======================================= */
	//parse the mdp options for pull
        bool operator()(const string& opt, const string& optval) throw(GMXMDP_Exception, FILEIO_Exception) = 0;
	/* ====================  DATA MEMBERS  ======================================= */
	//bit mask of the pulled dimension in left-to-right order (e.g., 001 <=> Y N N, meaning only x dimention is pulled)
	short dim;
	uint npgrps;
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
	    PULLGRP () {};                             /* constructor */

	    /* ====================  ACCESSORS     ======================================= */
            virtual void print() const = 0; /*  print out all the parameters read */

	    /* ====================  MUTATORS      ======================================= */

	    /* ====================  OPERATORS     ======================================= */
	    //parse the mdp options for pull-group
            virtual bool operator()(const string& opt, const string& optval, const string& i) throw(GMXMDP_Exception, FILEIO_Exception) = 0;

	  protected:
	    /* ====================  METHODS       ======================================= */

	    /* ====================  DATA MEMBERS  ======================================= */
	    vector<valtype> init, initB;
	    valtype rate, k, kB;

	  private:
	    /* ====================  METHODS       ======================================= */

	    /* ====================  DATA MEMBERS  ======================================= */

	}; /* -----  end of class PULLGRP  ----- */

      protected:
	/* ====================  METHODS       ======================================= */

	/* ====================  DATA MEMBERS  ======================================= */

      private:
	/* ====================  METHODS       ======================================= */

	/* ====================  DATA MEMBERS  ======================================= */

    }; /* -----  end of class PULL  ----- */

  private:
    /* ====================  METHODS       ======================================= */

    /* ====================  DATA MEMBERS  ======================================= */

}; /* -----  end of class MDP  ----- */

#endif
