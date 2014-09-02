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

uint Ensemble::cmp(const Ensemble& src) const {
  return uint(kB != src.getkB()) << DkB;
};

bool Ensemble::operator==(const Ensemble& src) const {
  return !(this->cmp(src));
};

bool Ensemble::operator!=(const Ensemble& src) const {
  return !(*this == src);
}

NVE::NVE(const valtype _kB, const Hamiltonian& _H) :
  Ensemble(_kB),
  H(_H)
  {}

uint NVE::cmp(const Ensemble& src) const {
  uint Qt = Ensemble::cmp(src);
  if(this == &src) { return Qt; }
  try {
    const NVE& _src = dynamic_cast<const NVE&>(src);
    Qt |= uint(this->H != _src.getH()) << DHamiltonian ;
    return Qt;
  } catch(const bad_cast& bcex) {
    return Qt;
  };
}

valtype NVE::ener(const vector<valtype>& vals) const {
  // Note that we don't do bound-check here
  // and we don't have temperature here
  return H.ener(vals);
}

NVT::NVT(const valtype _kB, const Hamiltonian& _H, const valtype _T):
  Ensemble(_kB),
  NVE(_kB, _H),
  T(_T)
  {}

uint NVT::cmp(const Ensemble& src) const {
  uint Qt = NVE::cmp(src);
  if(this == &src) { return Qt; }
  try {
    const NVT& _src = dynamic_cast<const NVT&>(src);
    Qt |= uint(this->T != _src.getT()) << DTemperature;
    return Qt;
  } catch (const bad_cast& bcex) {
    return Qt;
  }
}

valtype NVT::ener(const vector<valtype>& vals) const {
  // Note that we don't do bound-check here
  return H.ener(vals)/(kB*T);
}

NPT::NPT(const valtype _kB, const Hamiltonian& _H, const valtype _T, const valtype _P):
  Ensemble(_kB),
  NVE(_kB, _H),
  NVT(_kB, _H, _T),
  P(_P)
  {}

uint NPT::cmp(const Ensemble& src) const {
  uint Qt = NVT::cmp(src);
  if(this == &src) { return Qt; }
  try {
    const NPT& _src = dynamic_cast<const NPT&>(src);
    Qt |= (this->P != _src.getP());
    return Qt;
  } catch (const bad_cast& bcex) {
    return Qt;
  }
}

valtype NPT::ener(const vector<valtype>& vals) const {
  // Note that we don't do bound-check here
  // Also, we don't slice the input vector
  // by just assuming the hamiltonian functor
  // will use and only use as many as it needs
  // so that the last element is the volumn term
  return (H.ener(vals) + P*vals.back())/(kB*T);
}

