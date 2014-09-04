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

Hamiltonian::Hamiltonian() :
  Hamiltonian(vFunctVV{})
  {};

Hamiltonian::Hamiltonian(const Hamiltonian& _H) :
  Hamiltonian(_H.getPotentialFuncts())
  {};

Hamiltonian::Hamiltonian(Hamiltonian&& _H) :
  Hamiltonian(std::move(_H.getPotentialFuncts()))
  {};

Hamiltonian& Hamiltonian::operator=(const Hamiltonian& src) {
  potentials.clear();
  const auto& srcpots = src.getPotentialFuncts();
  for(const auto& srcpot : srcpots) potentials.emplace_back(srcpot);
  Npot = potentials.size();
  return *this;
}

vFunctVV& Hamiltonian::getPotentialFuncts() {return potentials; }
const vFunctVV& Hamiltonian::getPotentialFuncts() const {return potentials; }

valtype Hamiltonian::ener(const vector<valtype>& vals) const {
  valtype ans = 0.0;
  for(uint i = 0; i < Npot; ++i) {
    ans += potentials[i](vals[i]);
  }
  return ans;
}

bool Hamiltonian::operator==(const Hamiltonian& src) const {
  return (this == &src) || (this->potentials == src.getPotentialFuncts()); 
};

bool Hamiltonian::operator!=(const Hamiltonian& src) const {
  return !(*this == src);
}
