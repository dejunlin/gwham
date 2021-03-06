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
    //! Target constructor
    Hamiltonian(const vFunctVV& _potentials);
    //! Target move constructor
    Hamiltonian(vFunctVV&& _potentials);
    //! Empty constructor
    Hamiltonian();
    //! Copy constructor
    Hamiltonian(const Hamiltonian& _H);
    //! Move constructor
    Hamiltonian(Hamiltonian&& _H);
    //! Assignment
    Hamiltonian& operator=(const Hamiltonian& src);
    //! access the functors
    vFunctVV& getPotentialFuncts();
    const vFunctVV& getPotentialFuncts() const;
    //! return the total energy (only potential energy for now)
    valtype ener(const vector<valtype>& vals) const;
    //! check if the energy functors are the same
    bool operator==(const Hamiltonian& src) const;
    bool operator!=(const Hamiltonian& src) const;
  private:
    //! potential energy functors
    vFunctVV potentials;
    //! Number of potential energy functors
    uint Npot = 0;
};

#endif
