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
#include "mdp.hpp"

Ensemble_Factory::Ensemble_Factory(const vector<MDP*>& mdps) :
  QtMask(0)
{
  //we only check each one against the 1st mdp since we only need to know if the are different
  const MDP* const mdp0 = mdps[0];
  for(uint i = 1; i < mdps.size(); ++i) {
    QtMask |= mdp0->cmp(*mpds[i]);
  }
}
