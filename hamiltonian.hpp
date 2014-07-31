#if !defined(HAMILTONIAN_HPP)
#define HAMILTONIAN_HPP
#include "typedefs.hpp"

template <class ensemble>
class Hamiltonian {
  public:
    Hamiltonian();
    Hamiltonian(const vector<double>& _params); 
    Hamiltonian(const Hamiltonian<ensemble>& H);
    ~Hamiltonian();
    valtype ener(const vector<valtype>& vals) const; //take a set of conserved quantity and return the energy
    const ensemble* getens() const;
  private:
    const ensemble* const ens;
};

template <class ensemble>
Hamiltonian<ensemble>::Hamiltonian():
  ens(new ensemble())
  {}

template <class ensemble>
Hamiltonian<ensemble>::Hamiltonian(const vector<valtype>& _params):
  ens(new ensemble(_params))
  {}

template <class ensemble>
Hamiltonian<ensemble>::Hamiltonian(const Hamiltonian<ensemble>& H):
  ens(H.getens()) //TODO: we want a deep copy here by have ens return a list of parameters
  {}

  template <class ensemble>
Hamiltonian<ensemble>::~Hamiltonian() {
  delete ens;
}

template <class ensemble>
double Hamiltonian<ensemble>::ener(const vector<valtype>& vals) const {
  return ens->ener(vals); 
}

template <class ensemble>
const ensemble* Hamiltonian<ensemble>::getens() const {
  return ens; 
}
#endif
