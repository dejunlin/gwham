#if !defined(HAMILTONIAN_HPP)
#define HAMILTONIAN_HPP
#include "typedefs.hpp"

template <class ensemble>
class Hamiltonian {
  public:
    Hamiltonian();
    Hamiltonian(const vector<double>& _params); 
    Hamiltonian(const Hamiltonian<ensemble>& H); 
    valtype ener(const vector<valtype>& vals) const; //take a set of conserved quantity and return the energy
    ensemble getens() const;
  private:
    const ensemble ens;
};

template <class ensemble>
Hamiltonian<ensemble>::Hamiltonian():
  ens(ensemble())
  {}

template <class ensemble>
Hamiltonian<ensemble>::Hamiltonian(const vector<valtype>& _params):
  ens(ensemble(_params))
  {}

template <class ensemble>
Hamiltonian<ensemble>::Hamiltonian(const Hamiltonian<ensemble>& H):
  ens(H.getens())
  {}

template <class ensemble>
double Hamiltonian<ensemble>::ener(const vector<valtype>& vals) const {
  return ens.ener(vals); 
}

template <class ensemble>
ensemble Hamiltonian<ensemble>::getens() const {
  return ens; 
}
#endif
