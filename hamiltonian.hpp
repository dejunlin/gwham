#if !defined(HAMILTONIAN_HPP)
#define HAMILTONIAN_HPP
/*
 * =====================================================================================
 *
 *       Filename:  hamiltonian.hpp
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
#include "typedefs.hpp"
#include "functor.hpp"

class Hamiltonian {
  public:
    Hamiltonian(const vPFunct& _potentials);
    Hamiltonian(const Hamiltonian& _H);
    //! return the functors
    const vPFunct& getPotentialFuncts() const;
    //! return the total energy (only potential energy for now)
    valtype ener(const vector<valtype>& vals) const; 
  private:
    //! potential energy functors
    const vPFunct potentials; 
};

#endif
