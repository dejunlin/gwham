#include "ensemble.hpp"
#include <iostream>
#include <iterator>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

Ensemble::Ensemble():
  kB(0.0)
  {}

Ensemble::Ensemble(const double _kB):
  kB(_kB)
  {}

NVT::NVT(const double _kB, const double _T):
  Ensemble(_kB),
  T(_T)
  {}

NVT::NVT(const vector<double>& params):
  Ensemble(params[0]),
  T(params[1])
  {}

double NVT::ener(const vector<valtype>& vals) const {
  return vals[0]/(kB*T);
}

NPT::NPT(const double _kB, const double _T, const double _P):
  Ensemble(_kB),
  NVT(_kB, _T),
  P(_P)
  {}

NPT::NPT(const vector<double>& params):
  Ensemble(params[0]),
  NVT(params[1],params[2]),
  P(params[3])
  {}

double NPT::ener(const vector<valtype>& vals) const {
  return (vals[0] + P*vals[1])/(kB*T);
}

LAMBDA::LAMBDA() :
  Ensemble(),
  T(0.0),
  L(vector<double>(0,0.0))
  {}

LAMBDA::LAMBDA(const double _kB, const double _T, const vector<double>& _L) :
  Ensemble(_kB),
  T(_T),
  L(_L)
  {}

LAMBDA::LAMBDA(const vector<double>& params) :
  Ensemble(params[0]),
  T(params[1]),
  L(vector<double>(params.begin()+2,params.end()))
  {}

double LAMBDA::ener(const vector<valtype>& vals) const {
  double V = 0.0;
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

LAMBDAsgl::LAMBDAsgl(const double _kB, const double _T, const double& _L, const uint& _i) :
  Ensemble(_kB),
  T(_T),
  L(_L),
  i(_i)
  {}

LAMBDAsgl::LAMBDAsgl(const vector<double>& params) :
  Ensemble(params[0]),
  T(params[1]),
  L(params[2]),
  i(params[3])
  {}

double LAMBDAsgl::ener(const vector<valtype>& vals) const {
  return vals[i]*L/(kB*T);
}

vector<double> LAMBDAsgl::getparams() const {
  vector<double> params;
  params.push_back(kB);
  params.push_back(T);
  params.push_back(L);
  params.push_back(i);
  return params;
}

RST::RST() :
  Ensemble(),
  T(0.0),
  k(vector<double>(0,0.0)),
  r(vector<double>(0,0.0))
  {}

RST::RST(const double& _kB, const double& _T, const vector<double>& _k, const vector<double>& _r) :
  Ensemble(_kB),
  T(_T),
  k(_k),
  r(_r)
  {}

RST::RST(const vector<double>& params) :
  Ensemble(params[0]),
  T(params[1]),
  k(vector<double>(params.begin()+2,params.begin()+2+uint((params.size()-2)/2))),
  r(vector<double>(params.begin()+2+uint((params.size()-2)/2),params.end()))
  {
    cout << "# RST parameters: " << kB << " " << T << " ";
    copy(k.begin(),k.end(),ostream_iterator<double>(cout," "));
    copy(r.begin(),r.end(),ostream_iterator<double>(cout," "));
    cout << endl;
    if(k.size() != r.size()) {
      cerr << "The size of k and r vectors don't match in RST initialization" << endl;
      cerr << "k.size() = " << k.size() << endl;
      cerr << "r.size() = " << r.size() << endl;
      exit(-1);
    }
  }

double RST::ener(const vector<valtype>& vals) const {
  double V = 0.0;
  for(uint i = 0; i < vals.size(); ++i) {
    V += 0.5*k[i]*(vals[i]-r[i])*(vals[i]-r[i]);
  }
  //printf("#val = %30.15lf v = %30.15lf Vbeta = %30.15lf\n",vals[0],V,V/(kB*T));
  return V/(kB*T);
}

vector<double> RST::getparams() const {
  vector<double> params;
  params.push_back(kB);
  params.push_back(T);
  params.insert(params.end(),k.begin(),k.end());
  params.insert(params.end(),r.begin(),r.end());
  return params;
}

NVTL::NVTL(const double _kB, const double _T, const vector<double>& _L) :
  Ensemble(_kB),
  NVT(_kB,_T),
  L(_L)
  {}

NVTL::NVTL(const vector<double>& params) :
  Ensemble(params[0]),
  NVT(params[1],params[2]),
  L(vector<double>(params.begin()+3,params.end()))
  {}

double NVTL::ener(const vector<valtype>& vals) const {
  double V = vals[0];
  for(uint i = 0; i < L.size(); ++i) {
    V += L[i]*vals[i+1]; 
  }
  return V/(kB*T);
}

NPTL::NPTL(const double _kB, const double _T, const double _P, const vector<double>& _L) :
  Ensemble(_kB),
  NVT(_kB, _T),
  NPT(_kB, _T, _P),
  NVTL(_kB, _T, _L)
  {}

NPTL::NPTL(const vector<double>& params) :
  Ensemble(params[0]),
  NVT(params[0],params[1]),
  NPT(params[0],params[1],params[2]),
  NVTL(params[0],params[1],vector<double>(params.begin()+3,params.end()))
  {}

double NPTL::ener(const vector<valtype>& vals) const {
  double V = vals[0] + P*vals[1];
  for(uint i = 0; i < L.size(); ++i) {
    V += L[i]*vals[i+2]; 
  }
  return V/(kB*T);
}
