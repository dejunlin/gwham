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
#include "hamiltonian.hpp"
#include "timeseries.hpp"
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
    valtype kB = Boltzmannkcal;
    /** The following are expanded ensemble parameters.
     * The size of these vectors should be the same and
     * this size represents the number accessible states 
     * to the system.
     * NOTE that we only support restraints and temperate 
     * now. TODO: For other FEP parameters, we need to have 
     * a vector of potential energy functors for each type 
     * of FEP lambdas and add them to MDP::Hs via MDP::setexpand()
     */
    uint Nstates = 0;
    //! lambda parameters in Ls array
    enum FEPLambdas { Lbond, Lmass, Lvdw, Lcoul, Lrst, Ltemp, Lpress, NFEPLambdas };
    //! All the FEP lambdas
    array<vector<valtype>, NFEPLambdas> Ls;
    //! Restraint functors
    /** rstfuncts[i][j] is the j'th functor of the i'th state
     * NOTE that each state could have multiple functors, usually 
     * with each corresponding to 1 pull-group
     */
    vector<vFunctVV> rstfuncts;
    //! Temperatures
    vector<valtype> Ts;
    //! Presssures
    vector<valtype> Ps;
    //! Extra Hamiltonian
    vector<Hamiltonian> Hs;

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
      Ts = src.getTs();
      Ps = src.getPs();
      Hs = src.getHs();
      Nstates = src.getNstates();
      Ls = src.getLs(); 
      rstfuncts.clear();
      for(const auto& rstfunct : src.getrstfuncts()) {
	rstfuncts.emplace_back(rstfunct);
      }
      return *this;
    };
    /* ====================  ACCESSORS     ======================================= */
    bool hasTemperature() const { return Ts.size() != 0; }; 
    bool hasPressure() const { return !iszero(Ps); }; 
    bool hasFEPLambda() const { return getlambdas().size() != 0; }; 
    bool hasRestraint() const { return rstfuncts.size() != 0; };
    bool hasLbond() const { return !iszero(Ls[Lbond]); }
    bool hasLmass() const { return !iszero(Ls[Lmass]); }
    bool hasLvdw() const { return !iszero(Ls[Lvdw]); }
    bool hasLcoul() const { return !iszero(Ls[Lcoul]); }
    bool hasLrst() const { return !iszero(Ls[Lrst]); }
    bool hasLtemp() const { return !iszero(Ls[Ltemp]); }
    bool hasLpress() const { return !iszero(Ls[Lpress]); }
    bool isExpandedEnsemble() const { return Nstates > 1; };
    const valtype& getkB() const { return kB; };
    const uint& getNstates() const { return Nstates; };
    const vector<valtype>& getTs() const { return Ts; }
    const vector<valtype>& getPs() const { return Ps; }
    const vector<Hamiltonian>& getHs() const { return Hs; }
    const array<vector<valtype>, NFEPLambdas>& getLs() const { return Ls; }
    const vector<valtype>& getLbond() const { return Ls[Lbond]; }
    const vector<valtype>& getLmass() const { return Ls[Lmass]; }
    const vector<valtype>& getLvdw() const { return Ls[Lvdw]; }
    const vector<valtype>& getLcoul() const { return Ls[Lcoul]; }
    const vector<valtype>& getLrst() const { return Ls[Lrst]; }
    const vector<valtype>& getLtemp() const { return Ls[Ltemp]; }
    const vector<vFunctVV>& getrstfuncts() const { return rstfuncts; }
    //! print out all the parameters read 
    virtual void print() const = 0; 

    //! Return a list of lambdas
    virtual vector<valtype> getlambdas() const {
      vector<valtype> ans;
      for(const auto& L : Ls) ans.insert(ans.end(), L.begin(), L.end());
      return ans;
    };

    /* ====================  MUTATORS      ======================================= */

    /* ====================  OPERATORS     ======================================= */
    //! Time-series factory
    virtual vector<TimeSeries> gents() const = 0;

  protected:
    /* ====================  METHODS       ======================================= */

    /* ====================  DATA MEMBERS  ======================================= */

  private:
    /* ====================  METHODS       ======================================= */
    //! check and combine parameters; should be called at the end of the constructor
    virtual void doublechk() = 0;
    //! set the lambda vectors
    virtual void setLs() = 0;
    //! set the expanded ensemble states
    virtual void setexpand() = 0;

    /* ====================  DATA MEMBERS  ======================================= */

}; /* -----  end of class MDP  ----- */


#endif
