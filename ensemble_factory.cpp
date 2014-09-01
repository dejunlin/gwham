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

Ensemble_Factory::Ensemble_Factory(const vpMDP& mdps) :
  QtMask(0),
  ensembles(vPEns(0, NULL))
{
  vector<MDP> _mdps;
  //check for expanded ensemble
  for(const auto& pmdp : mdps) {
    _mdps.push_back(*pmdp);
    if(pmdp->isExpandedEnsemble()) {  continue; }
     
  }
  
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
      const MDP* const mdp = mdps[i];
      vPEns.push_back(new NVT(mdp->getkB(), mdp->getTemperature(), mdp->getPressure(), Hamiltonian(mdp->getrstfuncts()))); 
    }
  }
}

const vPEns& Ensemble_Factory::getEnsemble() const { 
  return ensembles;
}
