/*
 * =====================================================================================
 *
 *       Filename:  gmxmdp.cpp
 *
 *    Description:  GROMACS mdp file reader 
 *
 *        Version:  1.0
 *        Created:  31/07/14 08:46:23
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dejun Lin (DL), dejun.lin@gmail.com
 *   Organization:  Department of Biochemistry and Biophysics, Medical Center, University of Rochester
 *
 * =====================================================================================
 */
#include "gmxmdp.hpp"
#include "typedefs.hpp"

GMXMDP::GMXMDP(const string& _fname):
  MDP(_fname)
{
  pgeneric = new GMXGENERIC;
  pfep = new GMXFEP;
  ppull = new GMXPULL;

  fileio fio(fname, fstream::in, 1, 0, 1, MAXNLINE, ";");
  while(fio.readaline()) {
    const vector<string> tmpstr = fio.line2str("=");
    const string& opt = trimltcm(tmpstr[0], ";", " \t");
    const string& optval = trimltcm(tmpstr[1], ";", " \t");
    try {
      (*pgeneric)(opt, optval) || (*pfep)(opt, optval) || (*ppull)(opt, optval);
    } catch(MDP_Exception& gmxmdpex) {
      cerr << "Error understanding mdp file: " << fname << ": " << gmxmdpex.what() << endl;
      terminate();
    } catch(FILEIO_Exception& fioex) {
      cerr << "Error reading mdp file: " << fname << ": " << fioex.what() << endl;
      terminate();
    }
  }
  this->doublechk();
}


void GMXMDP::print() const {
  pgeneric->print();
  pfep->print();
  ppull->print();
}

uint GMXMDP::cmp(const MDP& mdp) const throw (MDP_Exception) {
  uint QtMask = MDP::cmp(mdp);
  //TODO: we need to handle if only Temperature-lambda are non-zero and expanded ensemble is turned on...
  //Turn off Restraint-Lambdas when Restraint is on. NOTE that position restraint 
  //can also be affected by Restraint-Lambdas. TODO: we need a better way to find out if position restraint interactions are affected.
  if((QtMask & Restraints) && (QtMask & RestraintLambdas)) {
    cout << "# In file " << fname << ": Restraint and FEP restraint-lambda are both turned on but we can only handle pulling interaction not position restraint, the latter of which is affected by restraint-lambda also; we manually turn off restraint-lambda here\n";
    QtMask ^= (1 << RestraintLambdas);
  }
  return QtMask;
}

void GMXMDP::doublechk() throw(MDP_Exception) {
  pgeneric->doublechk();
  pfep->doublechk();

  //we need to inform the pullgrp object about the FEP parameters
  this->setLPull();
  ppull->doublechk(); /*  this will also initialize the restraint functors */
}

void GMXMDP::setLPull() {
  vector<PULLGRP*>& ppgrps = ppull->ppullgrps; 
  switch (dynamic_cast<GMXPULL*>(ppull)->pullT) {
    case GMXPULL::Umbrella:
    case GMXPULL::UmbrellaFlatBottom: {
      vector<valtype>& Lrst = pfep->Lrst;
      switch (pfep->fepT) {
	case FEP::Yes: {
	  for(uint i = 0; i < ppgrps.size(); ++i) {
	    GMXPULLGRP* ppgrp = dynamic_cast<GMXPULLGRP*>(ppgrps[i]);
	    ppgrp->Lrst = vector<valtype>(1, Lrst[pfep->Linit]);
	  }
	  break;
	}
	case FEP::Expanded: {
	  for(uint i = 0; i < ppgrps.size(); ++i) {
	    GMXPULLGRP* ppgrp = dynamic_cast<GMXPULLGRP*>(ppgrps[i]);
	    ppgrp->Lrst = Lrst;
	  }
	  break;
	}
	default:
	  break;
      }
      break;
    }
    case GMXPULL::Contact: {
      if(ppgrps.size() > 3) {
        throw(MDP_Exception("Gromacs MDP can only have <= 3 contact groups for now"));
      }
      GMXFEP* pgmxfep = dynamic_cast<GMXFEP*>(pfep);
      vector<vector<valtype> > Lcnts;
      Lcnts.push_back(pgmxfep->Lcnt1);
      Lcnts.push_back(pgmxfep->Lcnt2);
      Lcnts.push_back(pgmxfep->Lcnt3);
      switch (pfep->fepT) {
	case FEP::Yes:
	  for(uint i = 0; i < ppgrps.size(); ++i) {
	    GMXPULLCNTGRP* ppgrp = dynamic_cast<GMXPULLCNTGRP*>(ppgrps[i]);
	    ppgrp->Lcnt = vector<valtype>(1, Lcnts[i][pfep->Linit]);
	  }
	  break;
	case FEP::Expanded:
	  for(uint i = 0; i < ppgrps.size(); ++i) {
	    GMXPULLCNTGRP* ppgrp = dynamic_cast<GMXPULLCNTGRP*>(ppgrps[i]);
	    ppgrp->Lcnt = Lcnts[i];
	  }
	  break;
	default:
	  break;
      }
      break;
    }
    default:
      break;
  }
}

void GMXMDP::GMXGENERIC::print() const {
  printf("#%20s = %-10.5f\n", "Temperature", T);
  printf("#%20s = %-10.5f\n", "Pressure", P);
}

void GMXMDP::GMXGENERIC::doublechk() throw(MDP_Exception) {
  if(T < 0) {
    throw(MDP_Exception("Temperature is negative"));
  }
  if(P < 0) {
    throw(MDP_Exception("Pressure is negative"));
  }
}

bool GMXMDP::GMXGENERIC::operator()(const string& opt, const string& optval) throw(MDP_Exception, FILEIO_Exception) {
  //first set the boltzman constant in GROMACS unit (kJ/(mol*K))
  kB = 0.0083144621;
  if(nocmmatchkey(opt, "ref-t") || nocmmatchkey(opt, "ref_t")) {
    setsc(optval, T);
  } else if(nocmmatchkey(opt, "ref-p") || nocmmatchkey(opt, "ref_p")) {
    setsc(optval, P);
  } else {
    return false;
  }
  return true;
}

GMXMDP::GMXFEP::GMXFEP() :
   Lcnt1(vector<valtype>(0, 0.0)),
   Lcnt2(vector<valtype>(0, 0.0)),
   Lcnt3(vector<valtype>(0, 0.0)),
   nstdhdl(-1),
   nstexpanded(-1),
   Tmc(-1)
{
}

void GMXMDP::GMXFEP::print() const {
  printf("#%20s = %-5d\n", "free-energy", fepT);
  printf("#%20s = ", "bonded-lambdas");
  copy(Lbond.begin(), Lbond.end(), ostream_iterator<valtype>(cout, " "));
  cout << endl;
  printf("#%20s = ", "mass-lambdas");
  copy(Lmass.begin(), Lmass.end(), ostream_iterator<valtype>(cout, " "));
  cout << endl;
  printf("#%20s = ", "vdw-lambdas");
  copy(Lvdw.begin(), Lvdw.end(), ostream_iterator<valtype>(cout, " "));
  cout << endl;
  printf("#%20s = ", "coul-lambdas");
  copy(Lcoul.begin(), Lcoul.end(), ostream_iterator<valtype>(cout, " "));
  cout << endl;
  printf("#%20s = ", "rst-lambdas");
  copy(Lrst.begin(), Lrst.end(), ostream_iterator<valtype>(cout, " "));
  cout << endl;
  printf("#%20s = ", "temperature-lambdas");
  copy(Ltemp.begin(), Ltemp.end(), ostream_iterator<valtype>(cout, " "));
  cout << endl;
  printf("#%20s = ", "cnt1-lambdas");
  copy(Lcnt1.begin(), Lcnt1.end(), ostream_iterator<valtype>(cout, " "));
  cout << endl;
  printf("#%20s = ", "cnt2-lambdas");
  copy(Lcnt2.begin(), Lcnt2.end(), ostream_iterator<valtype>(cout, " "));
  cout << endl;
  printf("#%20s = ", "cnt3-lambdas");
  copy(Lcnt3.begin(), Lcnt3.end(), ostream_iterator<valtype>(cout, " "));
  cout << endl;
  printf("#  %20s = %-5d\n", "Initial L-State", Linit);
  printf("#  %20s = %-5d\n", "nstdhdl", nstdhdl);
  printf("#  %20s = %-5d\n", "nstexpanded", nstexpanded);
  printf("#  %20s = %-10.5lf\n", "mc-temperature", Tmc);
}

void GMXMDP::GMXFEP::doublechk() throw(MDP_Exception) {
  if(fepT != FEP::NFEPTypes) {
    if(getAllLambdas().size() == 0) {
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

vector<valtype> GMXMDP::GMXFEP::getAllLambdas() const {
  vector<valtype> Lcnts(Lcnt1.begin(), Lcnt1.end());
  Lcnts.insert(Lcnts.end(), Lcnt2.begin(), Lcnt2.end());
  Lcnts.insert(Lcnts.end(), Lcnt3.begin(), Lcnt3.end());
  vector<valtype> allL = FEP::getAllLambdas();
  allL.insert(allL.end(), Lcnts.begin(), Lcnts.end());
  return allL;
}

vector<valtype> GMXMDP::GMXFEP::getInitLambdas() const {
  vector<valtype> initL = FEP::getInitLambdas();
  initL.push_back(Lcnt1[Linit]);
  initL.push_back(Lcnt2[Linit]);
  initL.push_back(Lcnt3[Linit]);
  return initL;
}

bool GMXMDP::GMXFEP::operator()(const string& opt, const string& optval) throw(MDP_Exception, FILEIO_Exception)
{
  if(nocmmatchkey(opt, "free-energy")) {
    if(nocmmatchkey(optval,"yes")) {
	fepT = Yes; 
    } else if(nocmmatchkey(optval, "expanded")) {
	fepT = Expanded;
    } else {
      throw(MDP_Exception("Unknown GMXFEP type: " + optval));
    }
  } else if(nocmmatchkey(opt, "bonded-lambdas")) {
    setvec(optval, Lbond);
  } else if(nocmmatchkey(opt, "mass-lambdas")) {
    setvec(optval, Lmass);
  } else if(nocmmatchkey(opt, "vdw-lambdas")) {
    setvec(optval, Lvdw);
  } else if(nocmmatchkey(opt, "coul-lambdas")) {
    setvec(optval, Lcoul);
  } else if(nocmmatchkey(opt, "restraint-lambdas")) {
    setvec(optval, Lrst);
  } else if(nocmmatchkey(opt, "temperature-lambdas")) {
    setvec(optval, Ltemp);
  } else if(nocmmatchkey(opt, "umbrella-contact1-lambdas")) {
    setvec(optval, Lcnt1);
  } else if(nocmmatchkey(opt, "umbrella-contact2-lambdas")) {
    setvec(optval, Lcnt2);
  } else if(nocmmatchkey(opt, "umbrella-contact3-lambdas")) {
    setvec(optval, Lcnt3);
  } else if(nocmmatchkey(opt, "init-lambda-state")) {
    setsc(optval, Linit);
  } else if(nocmmatchkey(opt, "nstdhdl")) {
    setsc(optval, nstdhdl);
  } else if(nocmmatchkey(opt, "nstexpanded")) {
    setsc(optval, nstexpanded);
  } else if(nocmmatchkey(opt, "mc-temperature")) {
    setsc(optval, Tmc);
  } else {
    return false;
  }
  return true;
}

GMXMDP::GMXPULL::GMXPULL() :
  pullT(NPullTypes),
  geomT(NGeomTypes),
  nstx(-1),
  nstf(-1)
{
}

void GMXMDP::GMXPULL::print() const {
  printf("#%20s = %-5d\n", "Pull", pullT);
  printf("#%20s = %-5d\n", "Pull-geometry", geomT);
  printf("#%20s = %-5d\n", "Pull-dim", dim);
  printf("#%20s = %-5d\n", "Pull-nstxout", nstx);
  printf("#%20s = %-5d\n", "Pull-nstfout", nstf);
  if(pullT == Contact) {
    printf("#%20s = %-5d\n", "Pull-ncontactgroups", npgrps);
    for(uint i = 0; i < npgrps; ++i) {
      printf("#%20s%-5d\n", "Pull-contactgroup", i);
      ppullgrps[i]->print();
    }
  } else {
    printf("#%20s = %-5d\n", "Pull-ngroups", npgrps);
    for(uint i = 0; i < npgrps; ++i) {
      printf("#%20s%-5d\n", "Pullgroup", i+1);
      ppullgrps[i]->print();
    }
  }
}

void GMXMDP::GMXPULL::doublechk() throw(MDP_Exception) {
  if(pullT != NPullTypes) {
    if(npgrps == 0) {
      throw(MDP_Exception("Pull is activated but no pull-group is defined"));
    }
    if(dim == 0) {
      throw(MDP_Exception("Pull is activated but pull dimension is 0"));
    } else {
      //here we treat the position restraint as if there's a functor for each dimension
      if((pullT == Umbrella || pullT == UmbrellaFlatBottom) && geomT == Position) {
	vector<uint> bits(0);
	for(uint i = 0; i < 2; ++i) { 
	  if(dim & (i << 1)) {
	    bits.push_back(i);
	  }
	}
	for(vector<PULLGRP*>::iterator i = ppullgrps.begin(); i != ppullgrps.end(); ++i) {
	  GMXPULLGRP* pi = dynamic_cast<GMXPULLGRP*>(*i);
	  const vector<valtype> init = pi->init;
	  const vector<valtype> initB = pi->initB;
	  //handle j = 0 case
	  pi->init = vector<valtype>(1, init[bits[0]]);
	  pi->initB = vector<valtype>(1, initB[bits[0]]);
	  for(uint j = 1; j < bits.size(); ++j) {
	    const uint bit = bits[j];
	    const vector<valtype> newinit(1, pi->init[bit]);
            const vector<valtype> newinitB(1, pi->initB[bit]);
	    i = ppullgrps.insert(i+1, new GMXPULLGRP(pi->vec, newinit, newinitB, pi->k, pi->kB, pi->r0, pi->r0B, pi->r1, pi->r1B, pi->k0, pi->k0B, pi->k1, pi->k1B, pi->Lrst)); 
	  }
	}
      }
    }

    for(vector<PULLGRP*>::iterator i = ppullgrps.begin(); i != ppullgrps.end(); ++i) {
      (*i)->doublechk();
    }
  }
}

bool GMXMDP::GMXPULL::operator()(const string& opt, const string& optval) throw(MDP_Exception, FILEIO_Exception)
{
  if(nocmmatchkey(opt, "pull")) {
    if(nocmmatchkey(optval,"umbrella")) {
      pullT = Umbrella; 
    } else if(nocmmatchkey(optval, "umbrella-flat-bottom")) {
      pullT = UmbrellaFlatBottom;
    } else if(nocmmatchkey(optval, "constraint")) {
      pullT = Constraint;
    } else if(nocmmatchkey(optval, "constant-force")) {
      pullT = ConstantForce;
    } else if(nocmmatchkey(optval, "contact")) {
      pullT = Contact;
    } else {
      throw(MDP_Exception("Unknown Pull type: '" + optval + "'"));
    }
  } else if(nocmmatchkey(opt, "pull_geometry") || nocmmatchkey(opt, "pull-geometry")) {
    if(nocmmatchkey(optval,"distance")) {
      geomT = Distance; 
    } else if(nocmmatchkey(optval, "direction")) {
      geomT = Direction; 
    } else if(nocmmatchkey(optval, "direction-periodic")) {
      geomT = DirectionPeriodic; 
    } else if(nocmmatchkey(optval, "cylinder")) {
      geomT = Cylinder; 
    } else if(nocmmatchkey(optval, "position")) {
      geomT = Position; 
    } else {
      throw(MDP_Exception("Unknown Pull-geometry type: '" + optval + "'"));
    }
  } else if(nocmmatchkey(opt, "pull_dim") || nocmmatchkey(opt, "pull-dim")) {
    vector<string> dims;
    setvec(optval, dims);
    if(dims.size() != 3) { throw(MDP_Exception("pull dimension must be <= 3 but we got: " + tostr(dims.size()))); }
    for(uint i = 0; i < 3; ++i) {
      if(nocmcmatchkey(dims[i], "Y")) {
	dim |= 1 << i;
      } else if(nocmcmatchkey(dims[i], "N")) {
	const short mask = numeric_limits<short>::max() ^ (1 << i) ;
	dim &= mask;
      } else {
        throw(MDP_Exception("Unknown Pull-dim: " + optval));
      }
    }
  } else if(nocmmatchkey(opt, "pull_nstxout") || nocmmatchkey(opt, "pull-nstxout")) {
    setsc(optval, nstx);
  } else if(nocmmatchkey(opt, "pull_start") || nocmmatchkey(opt, "pull-start")) {
    return true; //not really need this but just to avoid additional if-else
  } else if(nocmmatchkey(opt, "pull_nstfout") || nocmmatchkey(opt, "pull-nstfout")) {
    setsc(optval, nstf);
  } else if(nocmmatchkey(opt, "pull_ngroups") || nocmmatchkey(opt, "pull-ngroups")) {
    if(pullT == Contact) { return true; }
    setsc(optval, npgrps);
    for(uint i = 0; i < ppullgrps.size(); ++i) { if(ppullgrps[i]) { delete ppullgrps[i]; } }
    ppullgrps.clear();
    for(uint i = 0; i < npgrps; ++i) { ppullgrps.push_back(new GMXPULLGRP); }
  } else if(nocmmatchkey(opt, "pull_ncontactgroups") || nocmmatchkey(opt, "pull-ncontactgroups")) {
    setsc(optval, npgrps);
    for(uint i = 0; i < ppullgrps.size(); ++i) { if(ppullgrps[i]) { delete ppullgrps[i]; } }
    ppullgrps.clear();
    for(uint i = 0; i < npgrps; ++i) { ppullgrps.push_back(new GMXPULLCNTGRP);}
  } else if(nocmmatchkey(opt, "pull_group0") || nocmmatchkey(opt, "pull-group0")) {
    return true;
  } else {
    for(uint i = 0; i < npgrps; ++i) {
      const string gid = pullT == Contact ? tostr(i) : tostr(i+1);
      if(ppullgrps[i] && (*ppullgrps[i])(opt, optval, gid)) { return true; }
    }
    return false;
  } 
  //make sure we initialize pullgrps properly
  if(npgrps) {
    switch (pullT) {
      case Umbrella: {
        for(uint i = 0; i < npgrps; ++i) { 
	  ppullgrps[i]->setRSTT(PULLGRP::Quad);
	}
	break;
      }
      case UmbrellaFlatBottom: {
        for(uint i = 0; i < npgrps; ++i) { 
	  ppullgrps[i]->setRSTT(PULLGRP::QuadFlat);
	}
        break;
      }
      case Contact: {
        for(uint i = 0; i < npgrps; ++i) { 
	  ppullgrps[i]->setRSTT(PULLGRP::Quad);
	}
        break;
      }
      default:
	break;
    }
  }
  return true;
}

GMXMDP::GMXPULL::GMXPULLGRP::GMXPULLGRP() :
  vec(vector<valtype>(0, 0.0)),
  init(vector<valtype>(0, 0.0)),
  initB(vector<valtype>(0, 0.0)),
  k(-1),
  kB(-1),
  r0(-1),
  r0B(-1),
  r1(-1),
  r1B(-1),
  k0(-1),
  k0B(-1),
  k1(-1),
  k1B(-1),
  Lrst(vector<valtype>(0, 0.0))
{
}

GMXMDP::GMXPULL::GMXPULLGRP::GMXPULLGRP(
  const vector<valtype>& _vec, 
  const vector<valtype>& _init,
  const vector<valtype>& _initB,
  const valtype& _k,
  const valtype& _kB,
  const valtype& _r0,
  const valtype& _r0B,
  const valtype& _r1,
  const valtype& _r1B,
  const valtype& _k0,
  const valtype& _k0B,
  const valtype& _k1,
  const valtype& _k1B,
  const vector<valtype>& _Lrst) :
  vec(_vec),
  init(_init),
  initB(_initB),
  k(_k),
  kB(_kB),
  r0(_r0),
  r0B(_r0B),
  r1(_r1),
  r1B(_r1B),
  k0(_k0),
  k0B(_k0B),
  k1(_k1),
  k1B(_k1B),
  Lrst(_Lrst)
{
}

void GMXMDP::GMXPULL::GMXPULLGRP::print() const {
  printf("#  %20s = ", "pull_vec");
  copy(vec.begin(), vec.end(), ostream_iterator<valtype>(cout, " "));
  cout << endl;
  printf("#  %20s = ", "pull_init");
  copy(init.begin(), init.end(), ostream_iterator<valtype>(cout, " "));
  cout << endl;
  printf("#  %20s = ", "pull_initB");
  copy(initB.begin(), initB.end(), ostream_iterator<valtype>(cout, " "));
  cout << endl;
  printf("#  %20s = %-10.5lf\n", "pull_rate", rate);
  printf("#  %20s = %-10.5lf\n", "pull_k", k);
  printf("#  %20s = %-10.5lf\n", "pull_kB", kB);
  printf("#  %20s = %-10.5lf\n", "pull_r0", r0);
  printf("#  %20s = %-10.5lf\n", "pull_r0B", r0B);
  printf("#  %20s = %-10.5lf\n", "pull_r1", r1);
  printf("#  %20s = %-10.5lf\n", "pull_r1B", r1B);
  printf("#  %20s = %-10.5lf\n", "pull_k0", k0);
  printf("#  %20s = %-10.5lf\n", "pull_k0B", k0B);
  printf("#  %20s = %-10.5lf\n", "pull_k1", k1);
  printf("#  %20s = %-10.5lf\n", "pull_k1B", k1B);
  printf("#  %20s = %-5lu\n", "Nfunctors", rstfunct.size());
  for(uint i = 0; i < rstfunct.size(); ++i) {
    printf("#    %20s%-5u = ", "Functor", i);
    const vector<valtype> params = rstfunct[i].getParams();
    copy(params.begin(), params.end(), ostream_iterator<valtype>(cout, " "));
    cout << endl;
  }
}

void GMXMDP::GMXPULL::GMXPULLGRP::doublechk() throw(MDP_Exception) {
  switch (rstT) {
    case Quad: {
      for(uint i = 0; i < Lrst.size(); ++i) {
	const valtype& L = Lrst[i];
	const valtype ref = init[0] + (initB[0] - init[0])*L;
	const valtype fconst = (k + (kB - k)*L)/2;
	rstfunct.push_back(Quadratic(fconst, ref, 0));
      }
      break;
    }
    case QuadFlat: {
      for(uint i = 0; i < Lrst.size(); ++i) {
	const valtype& L = Lrst[i];
	const valtype ref = init[0] + (initB[0] - init[0])*L;
	const valtype Lr = r0 + (r0B - r0)*L;
	const valtype Rr = r1 + (r1B - r1)*L;
	const valtype Lk = (k0 + (k0B - k0)*L)/2;
	const valtype Rk = (k1 + (k1B - k1)*L)/2;
	rstfunct.push_back(QuadraticFlat(ref, Lr, Rr, Lk, Rk, 0));
      }
      break;
    }
    default:
      break;
  }
}

bool GMXMDP::GMXPULL::GMXPULLGRP::operator()(const string& opt, const string& optval, const string& i) throw(MDP_Exception, FILEIO_Exception)
{
  if(nocmmatchkey(opt, "pull_group" + i) || nocmmatchkey(opt, "pull-group" + i)) {
    return true; //not really need this
  } else if(nocmmatchkey(opt, "pull_vec" + i) || nocmmatchkey(opt, "pull-vec" + i)) {
    setvec(optval, vec);
  } else if(nocmmatchkey(opt, "pull_init" + i) || nocmmatchkey(opt, "pull-init" + i)) {
    setvec(optval, init);
  } else if(nocmmatchkey(opt, "pull_initB" + i) || nocmmatchkey(opt, "pull-initB" + i)) {
    setvec(optval, initB);
  } else if(nocmmatchkey(opt, "pull_rate" + i) || nocmmatchkey(opt, "pull-rate" + i)) {
    setsc(optval, rate);
  } else if(nocmmatchkey(opt, "pull_k" + i) || nocmmatchkey(opt, "pull-k" + i)) {
    setsc(optval, k);
  } else if(nocmmatchkey(opt, "pull_kB" + i) || nocmmatchkey(opt, "pull-kB" + i)) {
    setsc(optval, kB);
  } else if(nocmmatchkey(opt, "pull_umb_r0" + i)) {
    setsc(optval, r0);
  } else if(nocmmatchkey(opt, "pull_umb_r0B" + i)) {
    setsc(optval, r0B);
  } else if(nocmmatchkey(opt, "pull_umb_r1" + i)) {
    setsc(optval, r1);
  } else if(nocmmatchkey(opt, "pull_umb_r1B" + i)) {
    setsc(optval, r1B);
  } else if(nocmmatchkey(opt, "pull_umb_k0" + i)) {
    setsc(optval, k0);
  } else if(nocmmatchkey(opt, "pull_umb_k0B" + i)) {
    setsc(optval, k0B);
  } else if(nocmmatchkey(opt, "pull_umb_k1" + i)) {
    setsc(optval, k1);
  } else if(nocmmatchkey(opt, "pull_umb_k1B" + i)) {
    setsc(optval, k1B);
  } else {
    return false;
  }

  return true;
}

GMXMDP::GMXPULL::GMXPULLCNTGRP::GMXPULLCNTGRP() :
  nc(-1),
  ncB(-1),
  k(-1),
  kB(-1),
  Lcnt(vector<valtype>(0, 0.0))
  {}

void GMXMDP::GMXPULL::GMXPULLCNTGRP::print() const {
  printf("#  %20s = %-10.5lf\n", "rate", rate);
  printf("#  %20s = %-10.5lf\n", "nc", nc);
  printf("#  %20s = %-10.5lf\n", "ncB", ncB);
  printf("#  %20s = %-10.5lf\n", "k", k);
  printf("#  %20s = %-10.5lf\n", "kB", kB);
  printf("#  %20s = %-5lu\n", "Nfunctors", rstfunct.size());
  for(uint i = 0; i < rstfunct.size(); ++i) {
    printf("#    %20s%-5u = ", "Functor", i);
    const vector<valtype> params = rstfunct[i].getParams();
    copy(params.begin(), params.end(), ostream_iterator<valtype>(cout, " "));
    cout << endl;
  }
}

void GMXMDP::GMXPULL::GMXPULLCNTGRP::doublechk() throw(MDP_Exception) {
  switch (rstT) {
    case Quad: {
      for(uint i = 0; i < Lcnt.size(); ++i) {
	const valtype& L = Lcnt[i];
	const valtype ref = nc + (ncB - nc)*L;
	const valtype fconst = (k + (kB - k)*L)/2;
	rstfunct.push_back(Quadratic(fconst, ref, 0));
      }
      break;
    }
    default:
      break;
  }
}

bool GMXMDP::GMXPULL::GMXPULLCNTGRP::operator()(const string& opt, const string& optval, const string& i) throw(MDP_Exception, FILEIO_Exception)
{
  const string cntgrpnm = "pull_contactgroup" + i;
  if(nocmmatchkey(opt, cntgrpnm + "_pgrp1")) {
    return true; //not really need this
  } else if(nocmmatchkey(opt, cntgrpnm + "_pgrp2")) {
    return true; //not really need this
  } else if(nocmmatchkey(opt, cntgrpnm + "_ncrate")) {
    setsc(optval, rate);
  } else if(nocmmatchkey(opt, cntgrpnm + "_k")) {
    setsc(optval, k);
  } else if(nocmmatchkey(opt, cntgrpnm + "_kB")) {
    setsc(optval, kB);
  } else if(nocmmatchkey(opt, cntgrpnm + "_nc0")) {
    setsc(optval, nc);
  } else if(nocmmatchkey(opt, cntgrpnm + "_nc0B")) {
    setsc(optval, ncB);
  } else {
    return false;
  }
  return true;
}
