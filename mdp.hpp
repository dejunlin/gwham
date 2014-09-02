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
#include "exception.hpp"
#include "functor.hpp"
#include "metaprog_snippets.hpp"

using namespace std;
//! True if the vector is a vector of zeros
template < class T >
bool iszero(const vector<T>& op) {
  return op == vector<T>(op.size(), T(0));
}

//! True if 2 vector of type T are not the same and they're both non-zero
template < class T >
bool nonzero_diff(const vector<T>& op1, const vector<T>& op2) {
  return op1 != op2 && !( iszero(op1) && iszero(op2) );
}

/*
 * =====================================================================================
 *        Class:  MDP
 *  Description:  Interfacial class for MD-parameters reader
 * =====================================================================================
 */

class MDP
{
  public:
    //The thermodynamic quantities we care about
    enum Qt {Temperature, Pressure, Restraints, BondLambdas, MassLambdas, VdwLambdas, CoulLambdas, RestraintLambdas, TemperatureLambdas};
    //! File name 
    const string fname;
  protected:
    //! Boltzman constant
    valtype kB = 0;
    //! Reference temperature
    valtype T = -1;
    //! Reference pressure 
    valtype P = -1;
    //! Free-energy-perturbation type
    enum FEPType {No, Yes, Expanded, NFEPTypes} fepT = NFEPTypes;
    //! Initial lambda state id
    int Linit = -1;
    //! All the lambda vectors
    vector<valtype> Lbond, Lmass, Lvdw, Lcoul, Lrst, Ltemp;
    //! Restraint functors
    vector<vFunctVV> rstfuncts;

  public:
    /* ====================  LIFECYCLE     ======================================= */
    //! constructor
    MDP (const string& _fname) : 
      fname(_fname)
      {};
    //! copy constructor
    MDP (const MDP& src) : 
      MDP(src.fname)
      {
	*this = src;
      };
    //! copy assignment
    MDP& operator=(const MDP& src) {
      kB = src.getkB();
      T = src.getTemperature();
      P = src.getPressure();
      fepT = src.getfepT();
      Linit = src.getLinit();
      Lbond = src.getLbond();
      Lmass = src.getLmass();
      Lvdw = src.getLvdw();
      Lcoul = src.getLcoul();
      Lrst = src.getLrst();
      Ltemp = src.getLtemp();

      for(const auto& rstfunct : src.getrstfuncts()) {
	rstfuncts.emplace_back(rstfunct);
      }

      return *this;
    };
    /* ====================  ACCESSORS     ======================================= */
    bool hasTemperature() const { return T > 0; }; 
    bool hasPressure() const { return P > 0; }; 
    bool hasFEPLambda() const { return fepT != NFEPTypes && fepT != No; }; 
    bool hasLbond() const { return iszero(Lbond); }
    bool hasLmass() const { return iszero(Lmass); }
    bool hasLvdw() const { return iszero(Lvdw); }
    bool hasLcoul() const { return iszero(Lcoul); }
    bool hasLrst() const { return iszero(Lrst); }
    bool hasLtemp() const { return iszero(Ltemp); }
    bool isExpandedEnsemble() const { return fepT == Expanded; };
    bool hasRestraint() const { return rstfuncts.size() != 0; };
    uint NRestraints() const { return rstfuncts.size(); }
    const valtype& getkB() const { return kB; };
    const valtype& getTemperature() const { return T; }; 
    const valtype& getPressure() const { return P; };
    const FEPType& getfepT() const { return fepT; };
    const int& getLinit() const { return Linit; }
    const vector<valtype>& getLbond() const { return Lbond; }
    const vector<valtype>& getLmass() const { return Lmass; }
    const vector<valtype>& getLvdw() const { return Lvdw; }
    const vector<valtype>& getLcoul() const { return Lcoul; }
    const vector<valtype>& getLrst() const { return Lrst; }
    const vector<valtype>& getLtemp() const { return Ltemp; }
    const vector<vFunctVV>& getrstfuncts() const { return rstfuncts; }
    //! print out all the parameters read 
    virtual void print() const = 0; 
    //! compare 2 MDP objects and return a mask indicating what thermodynamic quantities are different
    virtual uint cmp(const MDP& mdp) const {
      uint QtMask = 0;
      if( mdp.hasTemperature() != this->hasTemperature() ) {
	throw(MDP_Exception("Comparing a constant-T MDP with a nonconstant-T MDP"));
      } else if( mdp.hasPressure() !=  this->hasPressure() ) {
	throw(MDP_Exception("Comparing a constant-P MDP with a nonconstant-P MDP"));
      } else if( mdp.hasFEPLambda() !=  this->hasFEPLambda() ) {
	throw(MDP_Exception("Comparing a FEP MDP with a non-FEP MDP"));
      } else if( mdp.hasRestraint() !=  this->hasRestraint() ) {
	throw(MDP_Exception("Comparing a Restrained MDP with a non-Restrained MDP."));
      } else if( mdp.NRestraints() !=  this->NRestraints() ) {
	throw(MDP_Exception("Comparing two Restrained MDP's of different number of restraint potential."));
      } else {
	QtMask |= uint(mdp.getTemperature() != this->getTemperature()) << Temperature;
	QtMask |= uint(mdp.getPressure() != this->getPressure()) << Pressure;
	//whenenver FEP is activated for both MDP, we assume that we have to histogram on the perturbation energy anyway...
	if( mdp.hasFEPLambda() && this->hasFEPLambda() ) {
	  // assume we already make sure all lambdas have the same dimension 
	  QtMask |= uint( nonzero_diff(mdp.getLbond(), this->getLbond()) ) << BondLambdas;
	  QtMask |= uint( nonzero_diff(mdp.getLmass(), this->getLmass()) ) << MassLambdas;
	  QtMask |= uint( nonzero_diff(mdp.getLvdw(), this->getLvdw()) ) << VdwLambdas;
	  QtMask |= uint( nonzero_diff(mdp.getLcoul(), this->getLcoul()) ) << CoulLambdas;
	  QtMask |= uint( nonzero_diff(mdp.getLrst(), this->getLrst()) ) << RestraintLambdas;
	  QtMask |= uint( nonzero_diff(mdp.getLtemp(), this->getLtemp()) ) << TemperatureLambdas;
	}
	QtMask |= uint( mdp.getrstfuncts() != this->getrstfuncts() ) << Restraints;
      }
      return QtMask;
    }; 

    //! Return a list of lambdas
    virtual vector<valtype> getlambdas() const {
      vector<valtype> ans;
      ans.insert(ans.end(), Lbond.begin(), Lbond.end());
      ans.insert(ans.end(), Lmass.begin(), Lmass.end());
      ans.insert(ans.end(), Lvdw.begin(), Lvdw.end());
      ans.insert(ans.end(), Lcoul.begin(), Lcoul.end());
      ans.insert(ans.end(), Lrst.begin(), Lrst.end());
      ans.insert(ans.end(), Ltemp.begin(), Ltemp.end());
      return ans;
    };

    /* ====================  MUTATORS      ======================================= */

    /* ====================  OPERATORS     ======================================= */

  protected:
    /* ====================  METHODS       ======================================= */

    /* ====================  DATA MEMBERS  ======================================= */

  private:
    /* ====================  METHODS       ======================================= */
    //! check and combine parameters; should be called at the end of the constructor
    virtual void doublechk() = 0; 

    /* ====================  DATA MEMBERS  ======================================= */

}; /* -----  end of class MDP  ----- */


#endif
