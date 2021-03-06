#if !defined(ENSEMBLE_FACTORY_HPP)
#define ENSEMBLE_FACTORY_HPP
/*
 * =====================================================================================
 *
 *       Filename:  ensemble_factory.hpp
 *
 *    Description:  Factory utility that generates instances of ensemble class
 *
 *        Version:  1.0
 *        Created:  31/07/14 13:23:24
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dejun Lin (DL), dejun.lin@gmail.com
 *   Organization:  Department of Biochemistry and Biophysics, Medical Center, University of Rochester
 *
 * =====================================================================================
 */
#include <vector>
#include <bitset>
#include <limits>
#include "typedefs.hpp"
#include "mdp.hpp"
#include "exception.hpp"
#include "ensemble.hpp"
#include "hamiltonian.hpp"

template <class PMDP>
typename std::enable_if<std::is_pointer<PMDP>::value, void>::type
chkmdp ( const PMDP& pmdp )
{
  //For simplicity, we only support restraints and temperature
  //for now
  if(pmdp->hasFEPLambda() && ( 
     pmdp->hasLbond() || pmdp->hasLmass() || pmdp->hasLvdw() || 
     pmdp->hasLcoul() )
    ) {
    throw(MDP_Exception("Can't use GWHAM to handle FEP for now"));
  }
}		/* -----  end of template function chkmdp  ----- */

template < class PMDP >
typename std::enable_if<std::is_pointer<PMDP>::value, void>::type
genens(const PMDP& pmdp, vpEnsemble& ens, multimap<PMDP, uint>& pmdp2ipens, bool combinestates) {
  //TODO: when we support expanded ensemble on other FEP parameter
  //set, we need to retrieve the corresponding vector of 
  //functors/parameters from the MDP object
  const uint Nstates = pmdp->getNstates();
  const auto& Ts = pmdp->getTs();
  const auto& Ps = pmdp->getPs();
  const vector<Hamiltonian>& Hs = pmdp->getHs();

  for(uint i = 0; i < Nstates; ++i) {
    pEnsemble newpens;
    if(pmdp->hasTemperature() && pmdp->hasPressure()) {
      //NPT
      newpens = make_shared<NPT>(pmdp->getkB(), Hs[i], Ts[i], Ps[i]);
    } else if(pmdp->hasTemperature()) {
      //NVT
      newpens = make_shared<NVT>(pmdp->getkB(), Hs[i], Ts[i]);
    } else if(pmdp->hasPressure()) {
      //NPE
      throw(MDP_Exception("Ensemble NPE not supported"));
    }  else {
      //NVE
      newpens = make_shared<NVE>(pmdp->getkB(), Hs[i]);
    }
    //if true, we'll check if this new ensemble is equivalent to 
    //the ensemble we've already created -- if yes, we skip pushing
    //it into ens
    if(combinestates) {
	bool skip = false;
	for(uint j = 0; j < ens.size(); ++j) {
	  const auto& oldpens = ens[j];
	  if(*newpens == *oldpens) {
	    pmdp2ipens.insert(make_pair(pmdp, j));
	    skip = true; 
	    break; 
	  }
	}
	if(skip) { continue; }
    }
    ens.emplace_back(newpens);
    pmdp2ipens.insert(make_pair(pmdp, ens.size()-1));
  }
}

template < class pE >
typename std::enable_if<is_pointer<pE>::value, uint>::type
cmpens(const vector<pE>& pens) {
  uint Qt = 0;
  for(auto& pen : pens) {
    if(typeid(*pen) != typeid(*pens[0])) {
      throw(Ensemble_Factory_Exception("Compared ensembles are not the same ensemble"));
    }
    Qt |= pen->cmp(*pens[0]);
  }
  return Qt;
}

void adjustens(vpEnsemble& ens) {
  const uint Qt = cmpens(ens);
  //! If pressures are consistent, try down-casting the pointer to NVT 
  if(!(Qt & (1<<Ensemble::DPressure))) {
    for(auto& pens : ens) { 
      auto _pens = dynamic_pointer_cast<NVT>(pens);
      if(_pens == nullptr) {
	//reaching here means the pens is pointer to upper class of NVT
	break;
      }
      pens.reset();
      pens = make_shared<NVT>(_pens->getkB(), _pens->getH(), _pens->getT()); 
    }
  }
  //! If temperautre are diff., add potential energy
  if(Qt & (1<<Ensemble::DTemperature)) {
    for(auto& pens : ens) { 
      auto _pens = dynamic_pointer_cast<NVE>(pens);
      _pens->getH().getPotentialFuncts().emplace_back(Make_Functor_Wrapper(FunctVV{}, Linear<valtype>, 1.0, 0.0)); 
    }
  }
}

template <class PMDP>
typename std::enable_if<std::is_pointer<PMDP>::value, vpEnsemble>::type
MDP2Ensemble (const vector<PMDP>& pmdps, multimap<PMDP, uint>& pmdp2ipens, const bool combinestates)
{
  vpEnsemble ans;
  for(const auto& pmdp : pmdps) {
    chkmdp(pmdp);   
    genens(pmdp, ans, pmdp2ipens, combinestates);
  }
  //! adjust ensembles if we need to take care of temperature and pressure
  adjustens(ans);

  return ans;
}		/* -----  end of template function generate_ensemble  ----- */

#endif
