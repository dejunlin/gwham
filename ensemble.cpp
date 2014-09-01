#include "typedefs.hpp"
#include "ensemble.hpp"
#include <iostream>
#include <iterator>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

Ensemble::Ensemble():
  kB(0.0)
  {}

Ensemble::Ensemble(const valtype _kB):
  kB(_kB)
  {}

NVE::NVE(const valtype _kB, const Hamiltonian& _H) :
  Ensemble(_kB),
  H(_H)
  {}

valtype NVE::ener(const vector<valtype>& vals) const {
  // Note that we don't do bound-check here
  // and we don't have temperature here
  return H.ener(vals);
}

NVT::NVT(const valtype _kB, const valtype _T, const Hamiltonian& _H):
  NVE(_kB, _H),
  T(_T)
  {}

valtype NVT::ener(const vector<valtype>& vals) const {
  // Note that we don't do bound-check here
  return H.ener(vals)/(kB*T);
}

NPT::NPT(const valtype _kB, const valtype _T, const Hamiltonian& _H, const valtype _P):
  Ensemble(_kB),
  NVT(_kB, _T, _H),
  P(_P)
  {}

valtype NPT::ener(const vector<valtype>& vals) const {
  // Note that we don't do bound-check here
  // Also, we don't slice the input vector
  // by just assuming the hamiltonian functor
  // will use and only use as many as it needs
  // so that the last element is the volumn term
  return (H.ener(vals) + P*vals.back())/(kB*T);
}

