#include "typedefs.hpp"
#include "gmbar.hpp"
#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

const valtype GMBAR::kB = 0.0083144621;
const valtype GMBAR::T = 300;


GMBAR::GMBAR(const vector<vector<vector<valtype> > >& dEs, const vector<uint>& N, const vector<valtype>& _f, const valtype& _tol):
  f(_f),
  tol(_tol)
{
  ulong count = 0;
  vector<valtype> newf(f);
  do {
    ++count;
    calnewf(dEs,N,newf);
  } while(!endit(newf,count));
}

void GMBAR::calnewf(const vector<vector<vector<valtype> > >& dEs, const vector<uint>& N, vector<valtype>& newf) const {
  const valtype beta = 1/(kB*T);
  //First reset all newf to zero
  newf = vector<valtype>(f.size(),0.0);
  vector<valtype> expnewf = vector<valtype>(f.size(),0.0); //expnewf[i] == exp(-newf[i])
  //Perform an update on newf
  for(uint k = 0; k < dEs.size(); ++k) { //Loop through all trajs
    const vector<vector<valtype> >& traj = dEs[k];
    for(uint t = 0; t < traj.size(); ++t) { //Loop through each frame of traj
      const vector<valtype>& dE = traj[t]; //dE here should have the size of f
                                           //dE[m] == dEs[k][t][m] == E[m](x[k][t]) - E[k](x[k][t])
					   //where x[k][t] is the configuration at time t from simulation k
					   //and E[m] is the Hamiltonian of simulation m
      //cout << "#Size of dE = " << dE.size() << endl;
      for(uint i = 0; i < expnewf.size(); ++i) {
        valtype sum_denom = 0.0;
	for(uint j = 0; j < expnewf.size(); ++j) {
	  const valtype exparg = f[j] - dE[j]*beta; 
          if(exparg > MAXEXPARG ) { cerr << k << " " << t << " " << i << " " << j << " " << f[j] << " " << dE[j] << " " << dE[i] << " exp("<<exparg<<") will overflow!\n"; exit(-1); }
          else if(exparg < MINEXPARG) { continue; }
	  sum_denom += N[j]*exp(exparg);
	}
	  const valtype exparg = -dE[i]*beta; 
          if(exparg > MAXEXPARG ) { cerr << "exp("<<exparg<<") will overflow!\n"; exit(-1); }
          else if(exparg < MINEXPARG) { continue; }
	expnewf[i] += exp(exparg)/sum_denom;
      }
    }
  }
  for(uint i = 0; i < newf.size(); ++i) {
    newf[i] = -log(expnewf[i]);
  }
  shiftf(newf);
}

const vector<valtype>& GMBAR::getf() const {
  return f;
}

bool GMBAR::endit(const vector<valtype>& newf, const ulong& count) {
  if(count > MAXIT) {
    cerr << "Max number of iteration " << MAXIT << " is reached. Quit\n";
    printfree();
    exit(-1);
  }
  if(count % 10 == 0) { cout << "#At the " << count << "'th iteration\n"; printfree(); } 
  for(uint i = 0; i < newf.size(); ++i) {
    if(fabs((newf[i] - f[i])/f[i]) > tol) {
      f = newf;
      return false;
    }
  }
  cout << "#MBAR converge to " << tol << " after " << count << " iterations" << endl;
  printfree();
  return true;
}

void GMBAR::printfree() const {
  printf("#%10s%30s\n","State","DimensionlessFreeEnergy");
  for(uint i = 0; i < f.size(); ++i) {
    printf("#%10d%30.15lf\n",i,f[i]);
  }
}

void GMBAR::shiftf(vector<valtype>& newf) const {
  for(int i = newf.size()-1; i >= 0; --i) {
    newf[i] -= newf[0];
  }
}
