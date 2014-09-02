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

const map<string, MDP::FEPType> GMXMDP::str2fepT = 
    {
      {"no", MDP::No},
      {"yes", MDP::Yes},
      {"expanded", MDP::Expanded},
    };

const map<string, GMXMDP::PullType> GMXMDP::str2pullT = 
{
  {"no", GMXMDP::No},
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

uint GMXMDP::cmp(const MDP& mdp) const {
  uint QtMask = MDP::cmp(mdp);
  //TODO: we need to handle if only Temperature-lambda 
  //are non-zero and expanded ensemble is turned on...
  //Turn off Restraint-Lambdas when Restraint is on. 
  //The catch here is that position restraint 
  //can also be affected by Restraint-Lambdas. 
  //TODO: we need a better way to find out if 
  //position restraint interactions are affected.
  if((QtMask & Restraints) && (QtMask & RestraintLambdas)) {
    cout << "# In file " << fname << ": Restraint and FEP restraint-lambda are both turned on but we can only handle pulling interaction not position restraint, the latter of which is affected by restraint-lambda also; we manually turn off restraint-lambda here\n";
    QtMask ^= (1 << RestraintLambdas);
  }
  return QtMask;

}

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
    uint NL = 0;
    if(Lbond.size()) { NL = Lbond.size(); } 
    if(Lmass.size()) { NL = Lmass.size(); } 
    if(Lvdw.size()) { NL = Lvdw.size(); } 
    if(Lcoul.size()) { NL = Lcoul.size(); } 
    if(Lrst.size()) { NL = Lrst.size(); } 
    if(Ltemp.size()) { NL = Ltemp.size(); } 
    if(Lcnt1.size()) { NL = Lcnt1.size(); } 
    if(Lcnt2.size()) { NL = Lcnt2.size(); } 
    if(Lcnt3.size()) { NL = Lcnt3.size(); }
    if(!Lbond.size()) { Lbond = vector<valtype>(NL, 0); }
    else if(Lbond.size() != NL) { throw(MDP_Exception("FEP bond-lambda must be of the same size of other lambdas")); }
    if(!Lmass.size()) { Lmass = vector<valtype>(NL, 0); }
    else if(Lmass.size() != NL) { throw(MDP_Exception("FEP mass-lambda must be of the same size of other lambdas")); }
    if(!Lvdw.size()) { Lvdw = vector<valtype>(NL, 0); }
    else if(Lvdw.size() != NL) { throw(MDP_Exception("FEP vdw-lambda must be of the same size of other lambdas")); }
    if(!Lcoul.size()) { Lcoul = vector<valtype>(NL, 0); }
    else if(Lcoul.size() != NL) { throw(MDP_Exception("FEP coul-lambda must be of the same size of other lambdas")); }
    if(!Lrst.size()) { Lrst = vector<valtype>(NL, 0); }
    else if(Lrst.size() != NL) { throw(MDP_Exception("FEP rst-lambda must be of the same size of other lambdas")); }
    if(!Ltemp.size()) { Ltemp = vector<valtype>(NL, 0); }
    else if(Ltemp.size() != NL) { throw(MDP_Exception("FEP temp-lambda must be of the same size of other lambdas")); }
    if(!Lcnt1.size()) { Lcnt1 = vector<valtype>(NL, 0); }
    else if(Lcnt1.size() != NL) { throw(MDP_Exception("FEP cnt1-lambda must be of the same size of other lambdas")); }
    if(!Lcnt2.size()) { Lcnt2 = vector<valtype>(NL, 0); }
    else if(Lcnt2.size() != NL) { throw(MDP_Exception("FEP cnt2-lambda must be of the same size of other lambdas")); }
    if(!Lcnt3.size()) { Lcnt3 = vector<valtype>(NL, 0); }
    else if(Lcnt3.size() != NL) { throw(MDP_Exception("FEP cnt3-lambda must be of the same size of other lambdas")); }
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
  //First is the temperature
  setMDPABpair(key2val, "ref-t", "sim-temp-high");
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
      lambdas.emplace_back(Lrst);
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
  uint Nstates = 0;
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
	throw(MDP_Exception("The number of rstfunctors for each pull-group are not consistent"));
    }
    Nstates = pgrp.Lrst.size();
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
      break;
    case Yes:
      Lbond = vector<valtype>(1, Lbond[Linit]);
      Lmass = vector<valtype>(1, Lmass[Linit]);
      Lvdw = vector<valtype>(1, Lvdw[Linit]);
      Lcoul = vector<valtype>(1, Lcoul[Linit]);
      Lrst = vector<valtype>(1, Lrst[Linit]);
      Ltemp = vector<valtype>(1, Ltemp[Linit]);
      Lcnt1 = vector<valtype>(1, Lcnt1[Linit]);
      Lcnt2 = vector<valtype>(1, Lcnt2[Linit]);
      Lcnt3 = vector<valtype>(1, Lcnt3[Linit]);
      break;
    default:
      Lbond = vector<valtype>(1, 0);
      Lmass = vector<valtype>(1, 0);
      Lvdw = vector<valtype>(1, 0);
      Lcoul = vector<valtype>(1, 0);
      Lrst = vector<valtype>(1, 0);
      Ltemp = vector<valtype>(1, 0);
      Lcnt1 = vector<valtype>(1, 0);
      Lcnt2 = vector<valtype>(1, 0);
      Lcnt3 = vector<valtype>(1, 0);
      break;  
  }
}

void GMXMDP::setexpand() {
  setAB();
  setLs();
  setpgrprstfuncts();
  //First temperature; only support linear scaling simulated tempering now
  Ts.clear();
  for(const auto& L : Ltemp) {Ts.emplace_back(Tmin + L*(Tmax-Tmin));}
  //Then restraints
  const auto Nstates = Lrst.size();
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
