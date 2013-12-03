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

NVT::NVT(const valtype _kB, const valtype _T):
  Ensemble(_kB),
  T(_T)
  {}

NVT::NVT(const vector<valtype>& params):
  Ensemble(params[0]),
  T(params[1])
  {}

valtype NVT::ener(const valtype& pot) const {
  return pot/(kB*T);
}

valtype NVT::ener(const vector<valtype>& vals) const {
  return vals[0]/(kB*T);
}

NPT::NPT(const valtype _kB, const valtype _T, const valtype _P):
  Ensemble(_kB),
  NVT(_kB, _T),
  P(_P)
  {}

NPT::NPT(const vector<valtype>& params):
  Ensemble(params[0]),
  NVT(params[1],params[2]),
  P(params[3])
  {}

valtype NPT::ener(const vector<valtype>& vals) const {
  return (vals[0] + P*vals[1])/(kB*T);
}

LAMBDA::LAMBDA() :
  Ensemble(),
  T(0.0),
  L(vector<valtype>(0,0.0))
  {}

LAMBDA::LAMBDA(const valtype _kB, const valtype _T, const vector<valtype>& _L) :
  Ensemble(_kB),
  T(_T),
  L(_L)
  {}

LAMBDA::LAMBDA(const vector<valtype>& params) :
  Ensemble(params[0]),
  T(params[1]),
  L(vector<valtype>(params.begin()+2,params.end()))
  {}

valtype LAMBDA::ener(const vector<valtype>& vals) const {
  valtype V = 0.0;
  for(uint i = 0; i < L.size(); ++i) {
    V += L[i]*vals[i]; 
  }
  return V/(kB*T);
}

LAMBDAsgl::LAMBDAsgl() :
  Ensemble(),
  T(0.0),
  L(0.0),
  i(0)
  {}

LAMBDAsgl::LAMBDAsgl(const valtype _kB, const valtype _T, const valtype& _L, const uint& _i) :
  Ensemble(_kB),
  T(_T),
  L(_L),
  i(_i)
  {
    cout << "# LAMBDAsgl parameters: " << kB << " " << T << " " << L << " " << i << endl;
  }

LAMBDAsgl::LAMBDAsgl(const vector<valtype>& params) :
  Ensemble(params[0]),
  T(params[1]),
  L(params[2]),
  i(params[3])
  {
    cout << "# LAMBDAsgl parameters: " << kB << " " << T << " " << L << " " << i << endl;
  }

valtype LAMBDAsgl::ener(const vector<valtype>& vals) const {
  return vals[i]*L/(kB*T);
}

vector<valtype> LAMBDAsgl::getparams() const {
  vector<valtype> params;
  params.push_back(kB);
  params.push_back(T);
  params.push_back(L);
  params.push_back(i);
  return params;
}

RST::RST() :
  Ensemble(),
  T(0.0),
  k(vector<valtype>(0,0.0)),
  r(vector<valtype>(0,0.0))
  {}

RST::RST(const valtype& _kB, const valtype& _T, const vector<valtype>& _k, const vector<valtype>& _r) :
  Ensemble(_kB),
  T(_T),
  k(_k),
  r(_r)
  {}

RST::RST(const vector<valtype>& params) :
  Ensemble(params[0]),
  T(params[1]),
  k(vector<valtype>(params.begin()+2,params.begin()+2+uint((params.size()-2)/2))),
  r(vector<valtype>(params.begin()+2+uint((params.size()-2)/2),params.end()))
  {
    cout << "# RST parameters: " << kB << " " << T << " ";
    copy(k.begin(),k.end(),ostream_iterator<valtype>(cout," "));
    copy(r.begin(),r.end(),ostream_iterator<valtype>(cout," "));
    cout << endl;
    if(k.size() != r.size()) {
      cerr << "The size of k and r vectors don't match in RST initialization" << endl;
      cerr << "k.size() = " << k.size() << endl;
      cerr << "r.size() = " << r.size() << endl;
      exit(-1);
    }
  }

valtype RST::ener(const vector<valtype>& vals) const {
  valtype V = 0.0;
  for(uint i = 0; i < vals.size(); ++i) {
    V += 0.5*k[i]*(vals[i]-r[i])*(vals[i]-r[i]);
  }
  //printf("#val = %30.15lf v = %30.15lf Vbeta = %30.15lf\n",vals[0],V,V/(kB*T));
  return V/(kB*T);
}

vector<valtype> RST::getparams() const {
  vector<valtype> params;
  params.push_back(kB);
  params.push_back(T);
  params.insert(params.end(),k.begin(),k.end());
  params.insert(params.end(),r.begin(),r.end());
  return params;
}

RST_fb::RST_fb() :
  Ensemble(),
  T(1.0),
  k0(vector<valtype>(0,0.0)),
  k1(vector<valtype>(0,0.0)),
  init(vector<valtype>(0,0.0)),
  r0(vector<valtype>(0,0.0)),
  r1(vector<valtype>(0,0.0))
  {}

RST_fb::RST_fb(const valtype& _kB, const valtype& _T, const vector<valtype>& _k0, const vector<valtype>& _k1, const vector<valtype>& _init, const vector<valtype>& _r0, const vector<valtype>& _r1) :
  Ensemble(_kB),
  T(_T),
  k0(_k0),
  k1(_k1),
  init(_init),
  r0(_r0),
  r1(_r1)
  {}

RST_fb::RST_fb(const vector<valtype>& params) :
  Ensemble(params[0]),
  T(params[1]),
  k0(vector<valtype>(params.begin()+2, params.begin()+2+uint((params.size()-2)/5))),
  k1(vector<valtype>(params.begin()+2+k0.size(), params.begin()+2+2*k0.size())),
  init(vector<valtype>(params.begin()+2+2*k0.size(), params.begin()+2+3*k0.size())),
  r0(vector<valtype>(params.begin()+2+3*k0.size(), params.begin()+2+4*k0.size())),
  r1(vector<valtype>(params.begin()+2+4*k0.size(), params.begin()+2+5*k0.size()))
  {
    cout << "# RST_fb parameters: " << kB << " " << T << " " << endl;
    cout << "#k0 = "; copy(k0.begin(),k0.end(),ostream_iterator<valtype>(cout," ")); cout << endl;
    cout << "#k1 = "; copy(k1.begin(),k1.end(),ostream_iterator<valtype>(cout," ")); cout << endl;
    cout << "#init = "; copy(init.begin(),init.end(),ostream_iterator<valtype>(cout," ")); cout << endl;
    cout << "#r0 = "; copy(r0.begin(),r0.end(),ostream_iterator<valtype>(cout," ")); cout << endl;
    cout << "#r1 = "; copy(r1.begin(),r1.end(),ostream_iterator<valtype>(cout," ")); cout << endl;
    cout << endl;
    if(k0.size() != k1.size() || k1.size() != init.size() || init.size() != r0.size() || r0.size() != r1.size()) {
      cerr << "The size of k0, k1, init, r0 or r1 vectors don't match in RST_fb initialization" << endl;
      cerr << "k0.size() = " << k0.size() << endl;
      cerr << "k1.size() = " << k1.size() << endl;
      cerr << "init.size() = " << init.size() << endl;
      cerr << "r0.size() = " << r0.size() << endl;
      cerr << "r1.size() = " << r1.size() << endl;
      exit(-1);
    }
  }

valtype RST_fb::ener(const vector<valtype>& vals) const {
  valtype V = 0.0;
  for(uint i = 0; i < vals.size(); ++i) {
    const valtype val = vals[i];
    const valtype dev = val - init[i];
    if( dev >= r0[i] && dev <= r1[i] ) { continue; }
    else if ( dev < r0[i] ) { V += 0.5*k0[i]*(dev-r0[i])*(dev-r0[i]); }
    else { V += 0.5*k1[i]*(dev-r1[i])*(dev-r1[i]); }
  }
  return V/(kB*T);
}

vector<valtype> RST_fb::getparams() const {
  vector<valtype> params;
  params.push_back(kB);
  params.push_back(T);
  params.insert(params.end(),k0.begin(),k0.end());
  params.insert(params.end(),k1.begin(),k1.end());
  params.insert(params.end(),init.begin(),init.end());
  params.insert(params.end(),r0.begin(),r0.end());
  params.insert(params.end(),r1.begin(),r1.end());
  return params;
}

RST_fbXLAMBDAsgl::RST_fbXLAMBDAsgl() : 
  Ensemble(),
  RST_fb(),
  LAMBDAsgl()
  {}

RST_fbXLAMBDAsgl::RST_fbXLAMBDAsgl(const valtype& _kB, const valtype& _T, const vector<valtype>& _k0, const vector<valtype>& _k1, const vector<valtype>& _init, const vector<valtype>& _r0, const vector<valtype>& _r1, const valtype& _L, const uint& _i) : 
  Ensemble(_kB),
  RST_fb(_kB,_T,_k0,_k1,_init,_r0,_r1),
  LAMBDAsgl(_kB,_T,_L,_i)
  {}

RST_fbXLAMBDAsgl::RST_fbXLAMBDAsgl(const vector<valtype>& params) :
  Ensemble(params[0]),
  RST_fb(vector<valtype>(params.begin(),params.end()-2)),
  LAMBDAsgl(params[0],params[1],params[params.size()-2],params[params.size()-1])
  {}

valtype RST_fbXLAMBDAsgl::ener(const vector<valtype>& vals) const {
  //First are vals for RST_fb
  const uint RST_size = k0.size();
  valtype V = RST_fb::ener(vector<valtype>(vals.begin(),vals.begin()+RST_size));
  //Then are for LAMBDAsgl
  //NOTE that there's no bound check here since we assume that LAMBDAsgl::i and LAMBDAsgl::L
  //will be 0 if there's virtually no elements in vals that are supposed to be processed by LAMBDAsgl
  V += LAMBDAsgl::ener(vals);
  //V is already weighted by kB*T in both RST_fb and LAMBDAsgl
  return V;
}

vector<valtype> RST_fbXLAMBDAsgl::getparams() const {
  vector<valtype> RST_fbparams = RST_fb::getparams();
  vector<valtype> LAMBDAsglparams = LAMBDAsgl::getparams();
  vector<valtype> params(RST_fbparams);
  params.insert(params.end(),LAMBDAsglparams.begin(),LAMBDAsglparams.end());
  return params;
}

RSTXLAMBDAsgl::RSTXLAMBDAsgl() : 
  Ensemble(),
  RST(),
  LAMBDAsgl()
  {}

RSTXLAMBDAsgl::RSTXLAMBDAsgl(const valtype& _kB, const valtype& _T, const vector<valtype>& _k, const vector<valtype>& _r, const valtype& _L, const uint& _i) : 
  Ensemble(_kB),
  RST(_kB,_T,_k,_r),
  LAMBDAsgl(_kB,_T,_L,_i)
  {}

RSTXLAMBDAsgl::RSTXLAMBDAsgl(const vector<valtype>& params) :
  Ensemble(params[0]),
  RST(vector<valtype>(params.begin(),params.end()-2)),
  LAMBDAsgl(params[0],params[1],params[params.size()-2],params[params.size()-1])
  {}

valtype RSTXLAMBDAsgl::ener(const vector<valtype>& vals) const {
  //First are vals for RST
  const uint RST_size = k.size();
  valtype V = RST::ener(vector<valtype>(vals.begin(),vals.begin()+RST_size));
  //Then are for LAMBDAsgl
  //NOTE that there's no bound check here since we assume that LAMBDAsgl::i and LAMBDAsgl::L
  //will be 0 if there's virtually no elements in vals that are supposed to be processed by LAMBDAsgl
  V += LAMBDAsgl::ener(vals);
  //V is already weighted by kB*T in both RST and LAMBDAsgl
  return V;
}

vector<valtype> RSTXLAMBDAsgl::getparams() const {
  vector<valtype> RSTparams = RST::getparams();
  vector<valtype> LAMBDAsglparams = LAMBDAsgl::getparams();
  vector<valtype> params(RSTparams);
  params.insert(params.end(),LAMBDAsglparams.begin(),LAMBDAsglparams.end());
  return params;
}

NVTL::NVTL(const valtype _kB, const valtype _T, const vector<valtype>& _L) :
  Ensemble(_kB),
  NVT(_kB,_T),
  L(_L)
  {}

NVTL::NVTL(const vector<valtype>& params) :
  Ensemble(params[0]),
  NVT(params[1],params[2]),
  L(vector<valtype>(params.begin()+3,params.end()))
  {}

valtype NVTL::ener(const vector<valtype>& vals) const {
  valtype V = vals[0];
  for(uint i = 0; i < L.size(); ++i) {
    V += L[i]*vals[i+1]; 
  }
  return V/(kB*T);
}

NPTL::NPTL(const valtype _kB, const valtype _T, const valtype _P, const vector<valtype>& _L) :
  Ensemble(_kB),
  NVT(_kB, _T),
  NPT(_kB, _T, _P),
  NVTL(_kB, _T, _L)
  {}

NPTL::NPTL(const vector<valtype>& params) :
  Ensemble(params[0]),
  NVT(params[0],params[1]),
  NPT(params[0],params[1],params[2]),
  NVTL(params[0],params[1],vector<valtype>(params.begin()+3,params.end()))
  {}

valtype NPTL::ener(const vector<valtype>& vals) const {
  valtype V = vals[0] + P*vals[1];
  for(uint i = 0; i < L.size(); ++i) {
    V += L[i]*vals[i+2]; 
  }
  return V/(kB*T);
}
