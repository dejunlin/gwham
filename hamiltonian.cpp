/*
 * =====================================================================================
 *
 *       Filename:  hamiltonian.cpp
 *
 *    Description:  Hamiltonian of various ensembles 
 *
 *        Version:  1.0
 *        Created:  11/08/14 13:55:02
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dejun Lin (DL), dejun.lin@gmail.com
 *   Organization:  Department of Biochemistry and Biophysics, Medical Center, University of Rochester
 *
 * =====================================================================================
 */

#include "hamiltonian.hpp"

Hamiltonian::Hamiltonian(const vFunctVV& _potentials) : 
  potentials(_potentials),
  Npot(potentials.size())
  {};

Hamiltonian::Hamiltonian(vFunctVV&& _potentials) :
  potentials(std::move(_potentials)),
  Npot(potentials.size())
  {};

Hamiltonian::Hamiltonian(const Hamiltonian& _H) :
  Hamiltonian(_H.getPotentialFuncts())
  {};

Hamiltonian::Hamiltonian(Hamiltonian&& _H) :
  Hamiltonian(std::move(_H.getPotentialFuncts()))
  {};

const vFunctVV& Hamiltonian::getPotentialFuncts() const {
  return potentials;
}

valtype Hamiltonian::ener(const vector<valtype>& vals) const {
  valtype ans = 0.0;
  for(uint i = 0; i < Npot; ++i) {
    ans += potentials[i](vals[i]);
  }
  return ans;
}
