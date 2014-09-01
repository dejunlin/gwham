#if !defined(ENSEMBLE_HPP)
#define ENSEMBLE_HPP

#include <vector>
#include "typedefs.hpp"
#include "hamiltonian.hpp"

using namespace std;

class Ensemble {
  public:
    Ensemble();
    Ensemble(const valtype _kB);
    //! return the dimensionless total energy (in the unit of kB*T)
    /*  NOTE that there's no bound-check nor any guarantee that 
     *  the user necessarily supplies the right number of parameters
     *  here. Under any circumstance we assume the corresponding subclass
     *  can use as many as elements from the input vector vals, 
     *  i.e., there's always enough elements in it for this function to work
     *  */
    virtual valtype ener(const vector<valtype>& vals) const = 0;
    virtual ~Ensemble() {};
  protected:
    //Boltzmann factor in whatever energy unit the user supplies
    const valtype kB; 
};

class NVE : public virtual Ensemble {
  public:
    NVE(const valtype _kB, const Hamiltonian& _H);
    //! NVE::ener() returns the energy not weighted by kB*T
    virtual valtype ener(const vector<valtype>& vals) const; 
    virtual ~NVE() {};
  protected:
    //! Hamiltonian of the system
    const Hamiltonian H; 
};

class NVT : public virtual NVE {
  public:
    NVT(const valtype _kB, const valtype _T, const Hamiltonian& _H);
    virtual valtype ener(const vector<valtype>& vals) const; 
    virtual ~NVT() {};
  protected:
    //! Temperature
    const valtype T;
};

class NPT : public virtual NVT {
  public:
    NPT(const valtype _kB, const valtype _T, const Hamiltonian& _H, const valtype _P);
    valtype ener(const vector<valtype>& vals) const; 
    virtual ~NPT() {};
  protected:
    //! Pressure
    const valtype P; 
};

#endif
