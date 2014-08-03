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
  protected:
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
        GENERIC () :
          T(-1),
          P(-1)
	  {};                             /* constructor */
	virtual ~GENERIC() {};
    
        /* ====================  ACCESSORS     ======================================= */
        virtual void print() const = 0; /*  print out all the parameters read */
        valtype getT() const { return T; }; 
        valtype getP() const { return P; }; 
    
        /* ====================  MUTATORS      ======================================= */
    
        /* ====================  OPERATORS     ======================================= */
        //parse the generic parameters
        virtual bool operator()(const string& opt, const string& optval) throw(MDP_Exception, FILEIO_Exception) = 0;
    
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
	enum FEPType {Yes, Expanded, NFEPTypes};
	FEPType fepT;

      public:
	/* ====================  LIFECYCLE     ======================================= */
	FEP () : 
	  fepT(NFEPTypes),  
          Lbond(vector<valtype>(0, 0.0)),
          Lmass(vector<valtype>(0, 0.0)),
          Lvdw(vector<valtype>(0, 0.0)),
          Lcoul(vector<valtype>(0, 0.0)),
          Lrst(vector<valtype>(0, 0.0)),
          Ltemp(vector<valtype>(0, 0.0)),
          Linit(-1)
          {};                             /* constructor */
	virtual ~FEP() {};

	/* ====================  ACCESSORS     ======================================= */
        virtual void print() const = 0; /*  print out all the parameters read */
	virtual FEPType getFEPT() const { return fepT; }
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
        virtual bool operator()(const string& opt, const string& optval) throw(MDP_Exception, FILEIO_Exception) = 0;
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
	/*
	 * =====================================================================================
	 *        Class:  PULLGRP
	 *  Description:  pull-group parameters
	 * =====================================================================================
	 */
	class PULLGRP
	{
	  public:
	    enum RestraintType { Quad, QuadFlat, NRestraintTypes };
	    RestraintType rstT;
	    valtype rate;
	    const Functor<valtype, valtype>* rstfunct;

	  public:
	    /* ====================  LIFECYCLE     ======================================= */
	    PULLGRP () :
	      rstT(NRestraintTypes),
              rate(0.0),
	      rstfunct(NULL)
	      {};                             /* constructor */

	    virtual ~PULLGRP() {
	      if(rstfunct) { delete rstfunct; }
	    }

	    /* ====================  ACCESSORS     ======================================= */
            virtual void print() const = 0; /*  print out all the parameters read */
	    RestraintType getRSTT() const { return rstT; }
	    const Functor<valtype, valtype>* getRSTF() const { return rstfunct; }

	    /* ====================  MUTATORS      ======================================= */
	    void setRSTT(const RestraintType& _rstT) { rstT = _rstT; }

	    /* ====================  OPERATORS     ======================================= */
	    //parse the mdp options for pull-group
            virtual bool operator()(const string& opt, const string& optval, const string& i) throw(MDP_Exception, FILEIO_Exception) = 0;

	  protected:
	    /* ====================  METHODS       ======================================= */

	    /* ====================  DATA MEMBERS  ======================================= */

	  private:
	    /* ====================  METHODS       ======================================= */

	    /* ====================  DATA MEMBERS  ======================================= */


	}; /* -----  end of class PULLGRP  ----- */
      public:
	/* ====================  LIFECYCLE     ======================================= */
	PULL () : 
	  dim(0),  
	  npgrps(0),
	  ppullgrps(vector<PULLGRP*>(npgrps, NULL))
	  {};                             /* constructor */

	virtual ~PULL() {
	  for(uint i = 0; i < ppullgrps.size(); ++i) {
	    if(ppullgrps[i]) {
	      delete ppullgrps[i];
	    }
	  }
	}

	/* ====================  ACCESSORS     ======================================= */
        virtual void print() const = 0; /*  print out all the parameters read */
	uint getNPG() const { return npgrps; }
	const vector<PULLGRP*>& getPPGRPS() const { return ppullgrps; }

	/* ====================  MUTATORS      ======================================= */

	/* ====================  OPERATORS     ======================================= */
	//parse the mdp options for pull
        virtual bool operator()(const string& opt, const string& optval) throw(MDP_Exception, FILEIO_Exception) = 0;
	/* ====================  DATA MEMBERS  ======================================= */
	//bit mask of the pulled dimension in left-to-right order (e.g., 001 <=> Y N N, meaning only x dimention is pulled)
	short dim;
	uint npgrps;

      protected:
	/* ====================  METHODS       ======================================= */

	/* ====================  DATA MEMBERS  ======================================= */
	vector<PULLGRP*> ppullgrps;

      private:
	/* ====================  METHODS       ======================================= */

	/* ====================  DATA MEMBERS  ======================================= */

    }; /* -----  end of class PULL  ----- */

  public:
    /* ====================  LIFECYCLE     ======================================= */
    MDP (const string& _fname) : 
      fname(_fname),
      pgeneric(NULL), 
      pfep(NULL), 
      ppull(NULL) 
      {};                             /* constructor */

    virtual ~MDP() {
      if(pgeneric) delete pgeneric;
      if(pfep) delete pfep;
      if(ppull) delete ppull;
    }

    /* ====================  ACCESSORS     ======================================= */
    virtual void print() const = 0; /*  print out all the parameters read */
    virtual uint cmp(const MDP& mdp) const throw (MDP_Exception) {
      uint QtMask = 0;
      if( mdp.hasTemperature() != this->hasTemperature() ) {
	throw(MDP_Exception("Comparing a constant-T MDP with a nonconstant-T MDP"));
      } else if( mdp.hasPressure() !=  this->hasPressure() ) {
	throw(MDP_Exception("Comparing a constant-P MDP with a nonconstant-P MDP"));
      } else if( mdp.hasFEPLambda() !=  this->hasFEPLambda() ) {
	throw(MDP_Exception("Comparing a FEP MDP with a non-FEP MDP"));
      } else if( mdp.getFEPType() !=  this->getFEPType() ) {
	throw(MDP_Exception("Comparing a MDP of FEP Type '" + tostr(mdp.getFEPType()) + "' to another of FEP type '" + tostr(this->getFEPType()) + "'"));
      } else if( mdp.getFEPLambda().size() !=  this->getFEPLambda().size() ) {
	throw(MDP_Exception("Comparing two MDP's of incompatible FEP Lambdas (of different dimensions)"));
      } else if( mdp.hasRestraint() !=  this->hasRestraint() ) {
	throw(MDP_Exception("Comparing a Restrained MDP with a non-Restrained MDP"));
      } else if( mdp.getRestraintType().size() !=  this->getRestraintType().size() ) {
	throw(MDP_Exception("Comparing two Restrained MDP's of incompatible restraint types (of different dimensions)"));
      } else {
	QtMask |= uint(mdp.getTemperature() != this->getTemperature()) << Temperature;
	QtMask |= uint(mdp.getPressure() != this->getPressure()) << Pressure;
	QtMask |= uint(mdp.getFEPLambda() != this->getFEPLambda()) << Lambdas;
	vector<const Functor<valtype, valtype>* > myfuncts = this->getRestraintFunctor();
	vector<const Functor<valtype, valtype>* > theirfuncts = mdp.getRestraintFunctor();
	for(uint i = 0; i < myfuncts.size(); ++i) {
	  if(*myfuncts[i] != *theirfuncts[i]) { 
	    QtMask |= 1 << Restraints; 
	    break; 
	  }
	}
      }
      return QtMask;
    }; /*  compare with another mdp object */

    bool hasTemperature() const { return pgeneric->getT() > 0; }; 
    bool hasPressure() const { return pgeneric->getP() > 0; };; 
    bool hasFEPLambda() const { return pfep->getFEPT() != FEP::NFEPTypes; }; 
    bool hasRestraint() const { return ppull->getNPG() != 0; }; 

    valtype getTemperature() const { return pgeneric->getT(); }; 
    valtype getPressure() const { return pgeneric->getP(); };
    FEP::FEPType getFEPType() const { return pfep->getFEPT(); };
    vector<valtype> getFEPLambda() const  { return pfep->getL(); };

    typedef PULL::PULLGRP PULLGRP;
    typedef PULLGRP::RestraintType RSTType;
    vector<RSTType> getRestraintType() const {
      const vector<PULLGRP*>& ppgrps = ppull->getPPGRPS();
      vector<RSTType> rsts;
      for(uint i = 0; i < ppgrps.size(); ++i) {
	rsts.push_back(ppgrps[i]->getRSTT());
      }
      return rsts;
    };

    vector<const Functor<valtype, valtype>* > getRestraintFunctor() const {
      const vector<PULLGRP*>& ppgrps = ppull->getPPGRPS();
      vector<const Functor<valtype, valtype>* > rstfuncts;
      for(uint i = 0; i < ppgrps.size(); ++i) {
	rstfuncts.push_back(ppgrps[i]->getRSTF());
      }
      return rstfuncts;
    };

    /* ====================  MUTATORS      ======================================= */

    /* ====================  OPERATORS     ======================================= */

  protected:
    /* ====================  METHODS       ======================================= */

    /* ====================  DATA MEMBERS  ======================================= */
    const string fname;
    GENERIC* pgeneric;
    FEP* pfep;
    PULL* ppull;

  private:
    /* ====================  METHODS       ======================================= */

    /* ====================  DATA MEMBERS  ======================================= */

}; /* -----  end of class MDP  ----- */

#endif
