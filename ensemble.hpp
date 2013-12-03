#if !defined(ENSEMBLE_HPP)
#define ENSEMBLE_HPP

#include <vector>
#include "typedefs.hpp"

using namespace std;

class Ensemble {
  public:
    Ensemble();
    Ensemble(const valtype _kB);
    virtual valtype ener(const vector<valtype>& vals) const = 0; //there's no need to constrain the dimension of the vector here 
                                                           //since each derived class just access from the beginning however 
							   //many terms they need and there's no guarantee that the user must 
							   //specify the right parameters anyway. 
							   //You need to make sure the API provides the bound-check
  protected:
    const valtype kB; //Boltzmann factor in whatever energy unit the user supplies
};

class NVT : public virtual Ensemble {
  public:
    NVT(const valtype _kB, const valtype _T);
    NVT(const vector<valtype>& params);
    valtype ener(const vector<valtype>& vals) const; //vals[0] -- the fundamental hamiltonian
    valtype ener(const valtype& pot) const; //vals[0] -- the fundamental hamiltonian
  protected:
    const valtype T; //temperature
};

class NPT : public virtual NVT {
  public:
    NPT(const valtype _kB, const valtype _T, const valtype _P);
    NPT(const vector<valtype>& params);
    valtype ener(const vector<valtype>& vals) const; //vals[0] -- the fundamental hamiltonian
                                               //vals[1] -- the volume
  protected:
    const valtype P; //pressure
};

//Just lambda parameters
//Diference between LAMBDA and NVTL should be that LAMBDA::T is static const
class LAMBDA : public virtual Ensemble {
  public:
    LAMBDA();
    LAMBDA(const valtype _kB, const valtype _T, const vector<valtype>& _L);
    LAMBDA(const vector<valtype>& params);
    valtype ener(const vector<valtype>& vals) const; 
  protected:
    const valtype T; //temperature
    const vector<valtype> L; //Lambda parameters
};

//Same as LAMBDA but a single point lambda
class LAMBDAsgl : public virtual Ensemble {
  public:
    LAMBDAsgl();
    LAMBDAsgl(const valtype _kB, const valtype _T, const valtype& _L, const uint& _i);
    LAMBDAsgl(const vector<valtype>& params);
    valtype ener(const vector<valtype>& vals) const; 
    vector<valtype> getparams() const;
  protected:
    const valtype T; //temperature
    const valtype L; //Lambda parameters
    const uint i; //which val to use in LAMBDAsgl::ener()
};

//Harmonic restraint potential which takes a vector of reaction coordinates and calculate the restraint energy
class RST : public virtual Ensemble {
  public:
    RST();
    RST(const valtype& _kB, const valtype& _T, const vector<valtype>& _k, const vector<valtype>& _r);
    RST(const vector<valtype>& params);
    valtype ener(const vector<valtype>& vals) const;
    vector<valtype> getparams() const;
  protected:
    const valtype T;
    const vector<valtype> k;
    const vector<valtype> r;
};

//Flat-bottom harmonic restraint potential
class RST_fb : public virtual Ensemble {
  public:
    RST_fb();
    RST_fb(const valtype& _kB, const valtype& _T, const vector<valtype>& _k0, const vector<valtype>& _k1, const vector<valtype>& _init, const vector<valtype>& _r0, const vector<valtype>& _r1);
    RST_fb(const vector<valtype>& params);
    valtype ener(const vector<valtype>& vals) const;
    vector<valtype> getparams() const;
  protected:
    const valtype T;
    const vector<valtype> k0;
    const vector<valtype> k1;
    const vector<valtype> init;
    const vector<valtype> r0;
    const vector<valtype> r1;
};

class RSTXLAMBDAsgl : public virtual RST, public virtual LAMBDAsgl {
  public:
    RSTXLAMBDAsgl();
    RSTXLAMBDAsgl(const valtype& _kB, const valtype& _T, const vector<valtype>& _k, const vector<valtype>& _r, const valtype& _L, const uint& _i);
    RSTXLAMBDAsgl(const vector<valtype>& params);
    //! This assume that elements of vals are ordered so that the first N are for RST while LAMBDAsgl::i is the index of the element in vals
    valtype ener(const vector<valtype>& vals) const;
    vector<valtype> getparams() const;
  private:
};

class RST_fbXLAMBDAsgl : public virtual RST_fb, public virtual LAMBDAsgl {
  public:
    RST_fbXLAMBDAsgl();
    RST_fbXLAMBDAsgl(const valtype& _kB, const valtype& _T, const vector<valtype>& _k0, const vector<valtype>& _k1, const vector<valtype>& _init, const vector<valtype>& _r0, const vector<valtype>& _r1, const valtype& _L, const uint& _i);
    RST_fbXLAMBDAsgl(const vector<valtype>& params);
    //! This assume that elements of vals are ordered so that the first N are for RST_fb while LAMBDAsgl::i is the index of the element in vals
    valtype ener(const vector<valtype>& vals) const;
    vector<valtype> getparams() const;
  private:
};

class NVTL : public virtual NVT {
  public:
    NVTL(const valtype _kB, const valtype _T, const vector<valtype>& _L);
    NVTL(const vector<valtype>& params);
    valtype ener(const vector<valtype>& vals) const; 
  protected:
    const vector<valtype> L; //Lambda parameters
};

class NPTL : public virtual NPT, public virtual NVTL {
  public:
    NPTL(const valtype _kB, const valtype _T, const valtype _P, const vector<valtype>& _L);
    NPTL(const vector<valtype>& params);
    valtype ener(const vector<valtype>& vals) const;
};

#endif
