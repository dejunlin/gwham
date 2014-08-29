#if !defined(ENSEMBLE_FACTORY_HPP)
#define ENSEMBLE_FACTORY_HPP
/*
 * =====================================================================================
 *
 *       Filename:  ensemble_factory.hpp
 *
 *    Description:  Factory utility that generates instances of ensemble class
 *
 *        Version:  1.0
 *        Created:  31/07/14 13:23:24
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dejun Lin (DL), dejun.lin@gmail.com
 *   Organization:  Department of Biochemistry and Biophysics, Medical Center, University of Rochester
 *
 * =====================================================================================
 */
#include "mdp.hpp"
#include "exception.hpp"
#include <vector>

/*
 * =====================================================================================
 *        Class:  Ensemble_Factory
 *  Description:  
 * =====================================================================================
 */
class Ensemble_Factory
{
  public:
    /* ====================  LIFECYCLE     ======================================= */
    Ensemble_Factory (const vector<MDP*>& mdps); 

    /* ====================  ACCESSORS     ======================================= */

    /* ====================  MUTATORS      ======================================= */

    /* ====================  OPERATORS     ======================================= */
    const vPEns& getEnsemble() const;

  protected:
    /* ====================  METHODS       ======================================= */

    /* ====================  DATA MEMBERS  ======================================= */
    //! Thermodynamic quantities masks 
    /* These bits (left to right) tell if the following 
     * thermodynamic quantites (in the same order defined 
     * in MDP::Qt) are different among different mdps: 
     * temperatures, pressures, lambdas, restraints
     */
    uint QtMask;
    //! vector of ensembles functor
    vPEns ensembles; 
  private:
    /* ====================  METHODS       ======================================= */

    /* ====================  DATA MEMBERS  ======================================= */

}; /* -----  end of class Ensemble_Factory  ----- */

#endif
