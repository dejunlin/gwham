#if !defined(ENSEMBLE_HPP)
#define ENSEMBLE_HPP

#include <vector>
#include "typedefs.hpp"

using namespace std;

class Ensemble {
  public:
    Ensemble();
    Ensemble(const double _kB);
    virtual double ener(const vector<valtype>& vals) const = 0; //there's no need to constrain the dimension of the vector here 
                                                           //since each derived class just access from the beginning however 
							   //many terms they need and there's no guarantee that the user must 
							   //specify the right parameters anyway. 
							   //You need to make sure the API provides the bound-check
  protected:
    const double kB; //Boltzmann factor in whatever energy unit the user supplies
};

class NVT : public virtual Ensemble {
  public:
    NVT(const double _kB, const double _T);
    NVT(const vector<double>& params);
    double ener(const vector<valtype>& vals) const; //vals[0] -- the fundamental hamiltonian
  protected:
    const double T; //temperature
};

class NPT : public virtual NVT {
  public:
    NPT(const double _kB, const double _T, const double _P);
    NPT(const vector<double>& params);
    double ener(const vector<valtype>& vals) const; //vals[0] -- the fundamental hamiltonian
                                               //vals[1] -- the volume
  protected:
    const double P; //pressure
};

//Just lambda parameters
//Diference between LAMBDA and NVTL should be that LAMBDA::T is static const
class LAMBDA : public virtual Ensemble {
  public:
    LAMBDA();
    LAMBDA(const double _kB, const double _T, const vector<double>& _L);
    LAMBDA(const vector<double>& params);
    double ener(const vector<valtype>& vals) const; 
  protected:
    const double T; //temperature
    const vector<double> L; //Lambda parameters
};

//Same as LAMBDA but a single point lambda
class LAMBDAsgl : public virtual Ensemble {
  public:
    LAMBDAsgl();
    LAMBDAsgl(const double _kB, const double _T, const double& _L, const uint& _i);
    LAMBDAsgl(const vector<double>& params);
    double ener(const vector<valtype>& vals) const; 
    vector<double> getparams() const;
  protected:
    const double T; //temperature
    const double L; //Lambda parameters
    const uint i; //which val to use in LAMBDAsgl::ener()
};

//Harmonic restraint potential which takes a vector of reaction coordinates and calculate the restraint energy
class RST : public virtual Ensemble {
  public:
    RST();
    RST(const double& _kB, const double& _T, const vector<double>& _k, const vector<double>& _r);
    RST(const vector<double>& params);
    double ener(const vector<valtype>& vals) const;
    vector<double> getparams() const;
  protected:
    const double T;
    const vector<double> k;
    const vector<double> r;
};

class NVTL : public virtual NVT {
  public:
    NVTL(const double _kB, const double _T, const vector<double>& _L);
    NVTL(const vector<double>& params);
    double ener(const vector<valtype>& vals) const; 
  protected:
    const vector<double> L; //Lambda parameters
};

class NPTL : public virtual NPT, public virtual NVTL {
  public:
    NPTL(const double _kB, const double _T, const double _P, const vector<double>& _L);
    NPTL(const vector<double>& params);
    double ener(const vector<valtype>& vals) const;
};

#endif
