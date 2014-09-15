#if !defined(ENSEMBLE_HPP)
#define ENSEMBLE_HPP

#include <vector>
#include "typedefs.hpp"
#include "hamiltonian.hpp"

using namespace std;

class Ensemble {
  public:
    enum Qt {DkB, DHamiltonian, DTemperature, DPressure};
    Ensemble();
    Ensemble(const valtype _kB);
    Ensemble(const Ensemble& src);
    //! return the dimensionless total energy (in the unit of kB*T)
    /*  NOTE that there's no bound-check nor any guarantee that 
     *  the user necessarily supplies the right number of parameters
     *  here. Under any circumstance we assume the corresponding subclass
     *  can use as many as elements from the input vector vals, 
     *  i.e., there's always enough elements in it for this function to work
     *  */
    virtual valtype ener(const vector<valtype>& vals) const = 0;
    //! Compare the parameters of two ensembles 
    /* NOTE: we are essentially comparing the units of kB here
     */ 
    virtual uint cmp(const Ensemble& src) const;
    //! true iff. this->cmp(src) == 0
    virtual bool operator==(const Ensemble& src) const;
    //! true iff. this->cmp(src) != 0
    virtual bool operator!=(const Ensemble& src) const;
    //! return kB
    const valtype& getkB() const { return kB; };
    virtual ~Ensemble() {};
  protected:
    //Boltzmann factor in whatever energy unit the user supplies
    const valtype kB; 
};

class NVE : public virtual Ensemble {
  public:
    NVE(const valtype _kB, const Hamiltonian& _H);
    NVE(const NVE& src);
    //! NVE::ener() returns the energy not weighted by kB*T
    virtual valtype ener(const vector<valtype>& vals) const;
    //! Compare the parameters of two ensembles
    virtual uint cmp(const Ensemble& src) const;
    //! access the Hamiltonian
    Hamiltonian& getH() { return H; };
    const Hamiltonian& getH() const { return H; };
    virtual ~NVE() {};
  protected:
    //! Hamiltonian of the system
    Hamiltonian H; 
};

class NVT : public virtual NVE {
  public:
    NVT(const valtype _kB, const Hamiltonian& _H, const valtype _T);
    NVT(const NVT& src);
    virtual valtype ener(const vector<valtype>& vals) const; 
    //! Compare the parameters of two ensembles
    virtual uint cmp(const Ensemble& src) const;
    //! return the temperature
    const valtype& getT() const { return T; }
    virtual ~NVT() {};
  protected:
    //! Temperature
    const valtype T;
};

class NPT : public virtual NVT {
  public:
    NPT(const valtype _kB, const Hamiltonian& _H, const valtype _T,  const valtype _P);
    NPT(const NPT& src);
    //! Compare the parameters of two ensembles
    virtual uint cmp(const Ensemble& src) const;
    //! return the pressure
    const valtype& getP() const { return P; }
    valtype ener(const vector<valtype>& vals) const; 
    virtual ~NPT() {};
  protected:
    //! Pressure
    const valtype P; 
};

#endif
