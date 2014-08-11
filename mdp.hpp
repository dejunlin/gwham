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
enum Qt {Temperature, Pressure, Restraints, BondLambdas, MassLambdas, VdwLambdas, CoulLambdas, RestraintLambdas, TemperatureLambdas};

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
	  kB(-1),
          T(-1),
          P(-1)
	  {};                             /* constructor */
	virtual ~GENERIC() {};
    
        /* ====================  ACCESSORS     ======================================= */
        virtual void print() const = 0; /*  print out all the parameters read */
	valtype getkB() const { return kB; };
        valtype getT() const { return T; }; 
        valtype getP() const { return P; }; 
    
        /* ====================  MUTATORS      ======================================= */
	virtual void doublechk() throw(MDP_Exception) = 0; /*  check and combine parameters */
    
        /* ====================  OPERATORS     ======================================= */
        //parse the generic parameters
        virtual bool operator()(const string& opt, const string& optval) throw(MDP_Exception, FILEIO_Exception) = 0;
    
      protected:
        /* ====================  METHODS       ======================================= */
    
        /* ====================  DATA MEMBERS  ======================================= */
        //! Boltzman constant, reference temperature and pressure 
	/* NOTE these parameters might not be relevant 
	 * to the WHAM calculation since some FEP parameters 
	 * might overwrite these
	 */
	 valtype kB, T, P;
    
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
	virtual vector<valtype> getAllLambdas() const {
	  vector<valtype> allL(Lbond.begin(), Lbond.end());
	  allL.insert(allL.end(), Lmass.begin(), Lmass.end());
	  allL.insert(allL.end(), Lvdw.begin(), Lvdw.end());
	  allL.insert(allL.end(), Lcoul.begin(), Lcoul.end());
	  allL.insert(allL.end(), Lrst.begin(), Lrst.end());
	  allL.insert(allL.end(), Ltemp.begin(), Ltemp.end());
	  return allL;
	}
	virtual vector<valtype> getInitLambdas() const {
	  vector<valtype> initL;
	  initL.push_back(Lbond[Linit]);
	  initL.push_back(Lmass[Linit]);
	  initL.push_back(Lvdw[Linit]);
	  initL.push_back(Lcoul[Linit]);
	  initL.push_back(Lrst[Linit]);
	  initL.push_back(Ltemp[Linit]);
	  return initL;
	}

	/* ====================  MUTATORS      ======================================= */
	virtual void doublechk() throw(MDP_Exception) = 0; /*  check and combine parameters */

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
	    vPFunct rstfunct;

	  public:
	    /* ====================  LIFECYCLE     ======================================= */
	    PULLGRP () :
	      rstT(NRestraintTypes),
              rate(0.0),
	      rstfunct(vPFunct(0, NULL))
	      {};                             /* constructor */

	    virtual ~PULLGRP() {
	      for(uint i = 0; i < rstfunct.size(); ++i) {
	        if(rstfunct[i]) { delete rstfunct[i]; }
	      }
	    }

	    /* ====================  ACCESSORS     ======================================= */
            virtual void print() const = 0; /*  print out all the parameters read */
	    RestraintType getRSTT() const { return rstT; }
	    const vPFunct& getRSTF() const { return rstfunct; }

	    /* ====================  MUTATORS      ======================================= */
	    void setRSTT(const RestraintType& _rstT) { rstT = _rstT; }
	    virtual void doublechk() throw(MDP_Exception) = 0; /*  check and combine parameters */

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

	/* ====================  MUTATORS      ======================================= */
	virtual void doublechk() throw(MDP_Exception) = 0; /*  check and combine parameters */

	/* ====================  OPERATORS     ======================================= */
	//parse the mdp options for pull
        virtual bool operator()(const string& opt, const string& optval) throw(MDP_Exception, FILEIO_Exception) = 0;
	/* ====================  DATA MEMBERS  ======================================= */
	//bit mask of the pulled dimension in left-to-right order (e.g., 001 <=> Y N N, meaning only x dimention is pulled)
	short dim;
	uint npgrps;
	vector<PULLGRP*> ppullgrps;

      protected:
	/* ====================  METHODS       ======================================= */

	/* ====================  DATA MEMBERS  ======================================= */

      private:
	/* ====================  METHODS       ======================================= */

	/* ====================  DATA MEMBERS  ======================================= */

    }; /* -----  end of class PULL  ----- */

  public:
    typedef PULL::PULLGRP::vPFunct vPFunct;
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
      } else if( mdp.hasRestraint() !=  this->hasRestraint() ) {
	throw(MDP_Exception("Comparing a Restrained MDP with a non-Restrained MDP. You might want to add a dummy potential to the non-restrained one."));
      } else if( mdp.getNPullGroups() !=  this->getNPullGroups() ) {
	throw(MDP_Exception("Comparing two Restrained MDP's of different number of pull-groups. You might want to add a dummy potential to the non-restrained one"));
      } else {
	QtMask |= uint(mdp.getTemperature() != this->getTemperature()) << Temperature;
	QtMask |= uint(mdp.getPressure() != this->getPressure()) << Pressure;
	//whenenver FEP is activated for both MDP, we assume that we have to histogram on the perturbation energy anyway...
	if( mdp.hasFEPLambda() && this->hasFEPLambda() ) {
	  const uint myLsize = this->getLbond().size(); /*  assume we already make sure all lambdas have the same dimension */
          const uint theirLsize = mdp.getLbond().size(); 
          const vector<valtype> myzeros(myLsize, 0); 
          const vector<valtype> theirzeros(theirLsize, 0);
	  QtMask |= uint( mdp.getLbond() != this->getLbond() && !(mdp.getLbond() != theirzeros && this->getLbond() != myzeros) ) << BondLambdas;
	  QtMask |= uint( mdp.getLmass() != this->getLmass() && !(mdp.getLmass() != theirzeros && this->getLmass() != myzeros) ) << MassLambdas;
	  QtMask |= uint( mdp.getLvdw() != this->getLvdw() && !(mdp.getLvdw() != theirzeros && this->getLvdw() != myzeros) ) << VdwLambdas;
	  QtMask |= uint( mdp.getLcoul() != this->getLcoul() && !(mdp.getLcoul() != theirzeros && this->getLcoul() != myzeros) ) << CoulLambdas;
	  QtMask |= uint( mdp.getLrst() != this->getLrst() && !(mdp.getLrst() != theirzeros && this->getLrst() != myzeros) ) << RestraintLambdas;
	  QtMask |= uint( mdp.getLtemp() != this->getLtemp() && !(mdp.getLtemp() != theirzeros && this->getLtemp() != myzeros) ) << TemperatureLambdas;
	}
	vPFunct myfuncts = this->getRestraintFunctor();
	vPFunct theirfuncts = mdp.getRestraintFunctor();
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
    bool hasBondLambda() const { 
      const vector<valtype>& Lbond = this->getLbond();
      const vector<valtype> zeros(Lbond.size(), 0.0);
      return Lbond != zeros;
    }
    bool hasMassLambda() const { 
      const vector<valtype>& Lmass = this->getLmass();
      const vector<valtype> zeros(Lmass.size(), 0.0);
      return Lmass != zeros;
    }
    bool hasVdwLambda() const { 
      const vector<valtype>& Lvdw = this->getLvdw();
      const vector<valtype> zeros(Lvdw.size(), 0.0);
      return Lvdw != zeros;
    }
    bool hasCoulLambda() const { 
      const vector<valtype>& Lcoul = this->getLcoul();
      const vector<valtype> zeros(Lcoul.size(), 0.0);
      return Lcoul != zeros;
    }
    bool hasTempLambda() const { 
      const vector<valtype>& Ltemp = this->getLtemp();
      const vector<valtype> zeros(Ltemp.size(), 0.0);
      return Ltemp != zeros;
    }
    
    valtype getkB() const { return pgeneric->getkB(); }
    valtype getTemperature() const { return pgeneric->getT(); }; 
    valtype getPressure() const { return pgeneric->getP(); };
    FEP::FEPType getFEPType() const { return pfep->getFEPT(); };
    const vector<valtype>& getLbond() const { return pfep->Lbond; }
    const vector<valtype>& getLmass() const { return pfep->Lmass; }
    const vector<valtype>& getLvdw() const { return pfep->Lvdw; }
    const vector<valtype>& getLcoul() const { return pfep->Lcoul; }
    const vector<valtype>& getLrst() const { return pfep->Lrst; }
    const vector<valtype>& getLtemp() const { return pfep->Ltemp; }
    valtype getLinit() const { return pfep->Linit; }
    vector<valtype> getFEPAllLambdas() const  { return pfep->getAllLambdas(); };
    vector<valtype> getFEPInitLambdas() const  { return pfep->getInitLambdas(); };

    typedef PULL::PULLGRP PULLGRP;
    typedef PULLGRP::RestraintType RSTType;
    
    uint getNPullGroups() const {
      return ppull->npgrps;
    }
    
    //! Return a vector of restraint functors from the pull-groups
    vPFunct getRestraintFunctor() const {
      const vector<PULLGRP*>& ppgrps = ppull->ppullgrps;
      vPFunct rstfuncts(0);
      for(uint i = 0; i < ppgrps.size(); ++i) {
	const vPFunct& functs = ppgrps[i]->getRSTF();
	rstfuncts.insert(rstfuncts.end(), functs.begin(), functs.end());
      }
      return rstfuncts;
    };

    /* ====================  MUTATORS      ======================================= */
    virtual void doublechk() throw(MDP_Exception) = 0; /*  check and combine parameters; should be called at the end of the constructor */

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
