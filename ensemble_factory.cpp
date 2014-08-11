/*
 * =====================================================================================
 *
 *       Filename:  ensemble_factory.cpp
 *
 *    Description:  Factory utility that generates instances of ensemble class
 *
 *        Version:  1.0
 *        Created:  31/07/14 15:00:21
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dejun Lin (DL), dejun.lin@gmail.com
 *   Organization:  Department of Biochemistry and Biophysics, Medical Center, University of Rochester
 *
 * =====================================================================================
 */

#include "ensemble_factory.hpp"
#include "exception.hpp"
#include "mdp.hpp"

Ensemble_Factory::Ensemble_Factory(const vector<MDP*>& mdps) throw(Ensemble_Factory_Exception):
  QtMask(0),
  ensembles(vPEns(0, NULL))
{
  //we only check each one against the 1st mdp since we only need to know if the are different
  const MDP* const mdp0 = mdps[0];
  for(uint i = 1; i < mdps.size(); ++i) {
    QtMask |= mdp0->cmp(*mpds[i]);
  }
  //we only handle restraint for now
  if(QtMask != (1 << Restraints)) {
    throw(Ensemble_Factory_Exception("we can only handle MDPs differing in restraints"));
  }

  //construct the ensemble functors. We pretend to handle the general
  //case, e.g., the ensemble can differ in any thermodynamic quantites, 
  //so that the developers know how they can add things they need here 
  if(QtMask & (1 << Restraints)) {
    for(uint i = 0; i < mdps.size(); ++i) {
      //we only need the NVT ensemble for this
      //NOTE that the restraint functors rstfuncts is a vector of 
      //pointer to functors, which are created in the MDP class 
      //and will be destroyed once the MDP object goes out of scope.
      //We should use smart-pointer for the functor pointer
      const MDP* const mdp = mdps[i];
      const vPFunct& rstfuncts = mdp->getRestraintFunctor();
      vPEns.push_back(new NVT(mdp->getkB(), mdp->getTemperature(), mdp->getPressure(), Hamiltonian(rstfuncts))); 
    }
  }
}
