/*
 * =====================================================================================
 *
 *       Filename:  gmxmdp.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  08/24/2014 10:41:25 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dejun Lin (), dejun.lin@gmail.com
 *   Organization:  Department of Biochemistry and Biophysics, Medical Center, University of Rochester
 *
 * =====================================================================================
 */

#include "gmxmdp.hpp"
#include "fileio.hpp"
#include "fileio_utils.hpp"
#include "typedefs.hpp"
#include <regex>

using namespace std;

const map<string, GMXMDP::FEPType> GMXMDP::str2fepT = 
    {
      {"no", GMXMDP::NoFEP},
      {"yes", GMXMDP::Yes},
      {"expanded", GMXMDP::Expanded},
    };

const map<string, GMXMDP::PullType> GMXMDP::str2pullT = 
{
  {"no", GMXMDP::NoPull},
  {"umbrella", GMXMDP::Umbrella},
  {"umbrella-flat-bottom", GMXMDP::UmbrellaFlatBottom},
  {"constraint", GMXMDP::Constraint},
  {"constant-force", GMXMDP::ConstantForce},
  {"contact", GMXMDP::Contact},
};

const map<string, GMXMDP::PullGeom> GMXMDP::str2geomT = 
{
  {"distance", GMXMDP::Distance},
  {"direction", GMXMDP::Direction},
  {"direction-periodic", GMXMDP::DirectionPeriodic},
  {"cylinder", GMXMDP::Cylinder},
  {"position", GMXMDP::Position},
};

const map<string, bool> GMXMDP::str2dim = 
{
  {"Y", 1},
  {"N", 0},
};

GMXMDP::GMXMDP(const string& fname) : MDP(fname), ifinit(this->initopts()) {
  fileio fio(fname, fstream::in, 1, 0, 1, MAXNLINE, ";");

  /** The followings the expected regex pattern for scalar and vector options in mdp
   * regex opt_sc(R"(^\s*(\S*)\s*=\s*([-+]?\d+[.]?\d*[eE]?\d*\s*)(;.*)?$)");
   * regex opt_vec(R"(^\s*(\S*)\s*=\s*([-+]?\d+[.]?\d*[eE]?\d*)(\s+[-+]?\d+[.]?\d*[eE]?\d*){1,}(\s*;.*)?$)");
   */

  const regex spliter(R"(\s*=\s*|(\s*;.*$)|\s+)");
  try 
  {
    while(fio.readaline()) {
      const string& line = fio.line;
      // Here we split all the mdp options by '=', ';' or empty spaces
      auto i = sregex_token_iterator(line.begin(), line.end(), spliter, -1);
      auto ie = sregex_token_iterator();
      setMDPopt(i, ie, key2val) || 
      setMDPopt(i, ie, key2str) ||
      setMDPopt(i, ie, key2int) ||
      setMDPopt(i, ie, key2vvec) ||
      setMDPopt(i, ie, key2vecstr);
      setenum();
      setmask();
      setpgrps();
      for(uint k = 0; k < pgrps.size(); ++k) {
        setMDPopt(i, ie, pgrps[k].key2val) || 
        setMDPopt(i, ie, pgrps[k].key2vvec);
      }
    }
    doublechk();
  } catch (MDP_Exception& gmxmdpex) {
    cerr << "Error understanding mdp file: " << fname << ": " << gmxmdpex.what() << endl;
    terminate();
  } catch (FILEIO_Exception& fioex) {
    cerr << "Error reading mdp file: " << fname << ": " << fioex.what() << endl;
    terminate();
  }

};

vector<valtype> GMXMDP::getlambdas() const {
  auto ans = MDP::getlambdas();
  ans.insert(ans.end(), Lcnt1.begin(), Lcnt1.end());
  ans.insert(ans.end(), Lcnt2.begin(), Lcnt2.end());
  ans.insert(ans.end(), Lcnt3.begin(), Lcnt3.end());
  return ans;
}

void GMXMDP::checkgeneric() {
  if(T != T) {
    throw(MDP_Exception("Temperature is not set"));
  }
  if(P != P) {
    throw(MDP_Exception("Pressure is not set"));
  }
}

void GMXMDP::checkfep() {
  if(fepT != NFEPTypes) {
    if(getlambdas().size() == 0) {
      throw(MDP_Exception("FEP is activated but no lambda is set"));
    }
    if(Linit < 0) {
      throw(MDP_Exception("FEP initial lambda must be a non-negative integer"));
    }
    vector<reference_wrapper<vector<valtype>>> _Ls = {
      Ls[Lbond], Ls[Lmass], Ls[Lvdw], Ls[Lcoul], Ls[Ltemp], Ls[Lpress],
      Lcnt1, Lcnt2, Lcnt3
    };

    uint NL = 0;
    for(auto& _L : _Ls) {
      auto& L = _L.get();
      if(NL && L.size() && NL != L.size()) {
	throw(MDP_Exception("Sizes of all FEP lambdas should be the same"));
      } else if(!NL && L.size()) { NL = L.size(); }
    }
    if(NL) {
      for(auto& _L : _Ls) { 
        auto& L = _L.get();
	if(!L.size()) {	L.resize(NL, 0); }
      } 
    }
  }
}

void GMXMDP::checkpull() {
  if(pullT != NPullTypes) {
    if(npgrps < 0) {
      throw(MDP_Exception("Pull is activated but no pull-group is defined"));
    }
    if(dim == 0) {
      throw(MDP_Exception("Pull is activated but pull dimension is 0"));
    }
    if(pullT == Contact && pgrps.size() > 3) {
      throw(MDP_Exception("Gromacs MDP can only have <= 3 contact groups for now"));
    }
  }
}

void GMXMDP::setAB() {
  //First set the temperature to ref-t
  setMDPABpair(key2val, "ref-t", "sim-temp-high");
  //Do the same thing to pressure
  setMDPABpair(key2val, "ref-p", "sim-press-high");
  //Now sync. the high and low pressure/temperature
  setMDPABpair(key2val, "high", "low");
  //Then are the pull-group parameters
  for(auto& pgrp : pgrps) {
    setMDPABpair(pgrp.key2val);
    setMDPABpair(pgrp.key2vvec);
  }
}

void GMXMDP::setpgrpLrst() {
  vector<vector<valtype> > lambdas;
  switch(pullT) {
    case Umbrella:
    case UmbrellaFlatBottom: {
      lambdas.emplace_back(Ls[Lrst]);
      break;
    }
    case Contact: {
      lambdas.emplace_back(Lcnt1);
      lambdas.emplace_back(Lcnt2);
      lambdas.emplace_back(Lcnt3);
      break;
    }
    default:
      break;
  }
  vector<valtype> lambda; 
  for(GMXPGRP& pgrp : pgrps) {
    if(lambdas.size()) { 
      lambda = lambdas.front(); 
      lambdas.erase(lambdas.begin()); 
    }
    pgrp.Lrst = lambda;
  }
}

void GMXMDP::expandpgrps() {
  //here we treat the position restraint as if there's a functor for each dimension
  if((pullT != Umbrella && pullT != UmbrellaFlatBottom) || geomT != Position) {
    return; 
  }
  vector<uint> bits(0);
  for(uint i = 0; i < 3; ++i) { if(dim & (i << 1)) { bits.push_back(i); } }
  for(vector<GMXPGRP>::iterator i = pgrps.begin(); i != pgrps.end(); ++i) {
    const GMXPGRP& pgrp = *i;
    const vector<valtype> init = pgrp.init;
    const vector<valtype> initB = pgrp.initB;
    //handle j = 0 case
    i->init = vector<valtype>(1, init[bits[0]]);
    i->initB = vector<valtype>(1, initB[bits[0]]);
    for(uint j = 1; j < bits.size(); ++j) {
      i = pgrps.emplace(i+1, pgrp);
      const uint bit = bits[j];
      //! now i points to the new pgrp
      i->init = vector<valtype>(1, init[bit]);
      i->initB = vector<valtype>(1, initB[bit]);
    }
  }
}

void GMXMDP::setpgrprstfuncts() {
  setpgrpLrst();
  expandpgrps();
  for(GMXPGRP& pgrp : pgrps) {
    switch(pgrp.rstT) {
      case GMXPGRP::Quad: {
	for(const valtype& L : pgrp.Lrst) {
	  const valtype ref = (1-L)*pgrp.init[0] + L*pgrp.initB[0];
	  const valtype fconst = ((1-L)*pgrp.k + L*pgrp.kB)/2;
	  pgrp.rstfuncts.push_back( Make_Functor_Wrapper(FunctVV{}, Quadratic<valtype>, fconst, ref, 0.0) );
	}
	break;
      }
      case GMXPGRP::QuadFlat: {
	for(const valtype& L : pgrp.Lrst) {
	  const valtype ref = (1-L)*pgrp.init[0] + L*pgrp.initB[0];
	  const valtype Lr = (1-L)*pgrp.r0 + L*pgrp.r0B;
	  const valtype Rr = (1-L)*pgrp.r1 + L*pgrp.r1B;
	  const valtype Lk = ((1-L)*pgrp.k0 + L*pgrp.k0B)/2;
	  const valtype Rk = ((1-L)*pgrp.k1 + L*pgrp.k1B)/2;
	  pgrp.rstfuncts.push_back( Make_Functor_Wrapper(FunctVV{}, QuadraticFlat<valtype>, ref, Lr, Rr, Lk, Rk, 0.0) );
	}
	break;
      }
      default:
        break;
    }
    if(Nstates && Nstates != pgrp.Lrst.size()) { 
	throw(MDP_Exception("The number of rstfunctors should be the same of the number of FEP states in expanded ensemble"));
    } else if(!Nstates && pgrp.Lrst.size()) { Nstates = pgrp.Lrst.size(); }
  }
}

void GMXMDP::doublechk() {
  checkgeneric();
  checkfep();
  checkpull();
  setexpand();
};

void GMXMDP::setLs() {
  switch(fepT) {
    case Expanded:
      Nstates = Ls[0].size();
      break;
    case Yes:
      for(auto& L : Ls) L = { L[Linit] };
      Lcnt1 = { Lcnt1[Linit] };
      Lcnt2 = { Lcnt2[Linit] };
      Lcnt3 = { Lcnt3[Linit] };
      Nstates = 1;
      break;
    default:
      for(auto& L : Ls) L = { 0 };
      Lcnt1 = { 0 };
      Lcnt2 = { 0 };
      Lcnt3 = { 0 };
      Nstates = 1;
      break;  
  }
}

void GMXMDP::setexpand() {
  setAB();
  setLs();
  setpgrprstfuncts();
  //First temperature; only support linear scaling simulated tempering now
  Ts.clear();
  for(const auto& L : Ls[Ltemp]) {Ts.emplace_back(Tmin + L*(Tmax-Tmin));}
  //Then pressure; only linear scaling
  Ps.clear();
  for(const auto& L : Ls[Lpress]) {Ps.emplace_back(Pmin + L*(Pmax-Pmin));}
  //Then restraints
  rstfuncts.resize(Nstates);
  for(const auto& pgrp : pgrps) {
    const auto& pgrprst = pgrp.rstfuncts;
    for(uint i = 0; i < Nstates; ++i) {
      rstfuncts[i].emplace_back(pgrprst[i]);
    }
  }
}

void GMXMDP::print() const {
  printMDPopt(key2val);
  printMDPopt(key2int);
  printMDPopt(key2str);
  printMDPopt(key2vvec);
  printMDPopt(key2vecstr);
  for(const GMXPGRP& pgrp : pgrps) {
    printf("#%27s%-3d\n", "Pull-group",pgrp.gid);
    printMDPopt(pgrp.key2val);
    printMDPopt(pgrp.key2vvec);
    uint i = 0;
    for(const FunctVV& rstfunct : pgrp.rstfuncts) {
      printf("#%27s%-3d = ", "Restraint Functor",i);
      cout << rstfunct.getParams() << endl;
      ++i;
    }
  }
};

void GMXMDP::setpgrps() {
  if((npgrps < 0 && ncntgrps < 0) || pgrps.size()) { return; }
  switch(pullT) {
    case NPullTypes:
      throw(MDP_Exception("Pull-type must be set before pull-group"));
      break;
    case Contact: {
      for(int k = 0; k < ncntgrps; ++k) { pgrps.emplace_back(k); }
      break; 
    }
    default: {
      for(int k = 0; k < npgrps; ++k) { pgrps.emplace_back(k+1); }
      break;
    }
  }
};

void GMXMDP::setenum() {
  if(fepT == NFEPTypes) { setMDPenum(str2fepT, fepTstr, fepT); }
  if(pullT == NPullTypes) { setMDPenum(str2pullT, pullTstr, pullT); }
  if(geomT == NGeomTypes) { setMDPenum(str2geomT, geomTstr, geomT); }
  if(pgrps.size() && pgrps[0].rstT == GMXPGRP::NRestraintTypes) {
    GMXPGRP::RestraintType rstT = GMXPGRP::NRestraintTypes;
    switch(pullT) {
      case Umbrella:
      case Contact: {
        rstT = GMXPGRP::Quad;
        break;
      }
      case UmbrellaFlatBottom : {
        rstT = GMXPGRP::QuadFlat;
        break;
      }
      default:
	break;
    }
    for(GMXPGRP& pgrp : pgrps) { pgrp.rstT = rstT; }
  }
}

void GMXMDP::setmask() {
  if(dim != 0) { return; }
  for(uint i = 0; i < dimstr.size(); ++i) {
    const string& bitstr = dimstr[i];
    const map<string, bool>::const_iterator it = str2dim.find(bitstr);
    if(it != str2dim.end()) {
      if(it->second) { dim |= 1 << i; } 
      else { dim &= (numeric_limits<int>::max() ^ (1 << i)); }
    }
  }
}
