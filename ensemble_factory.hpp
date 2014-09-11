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
genens(const PMDP& pmdp, vpEnsemble& ens, bool combinestates) {
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
    //TODO: in case we have redundant ensembles, we need to create 
    //a map to tell which mdp coresponds to which ensembles
    if(combinestates) {
	bool skip = false;
	for(const auto& oldpens : ens) {
	  if(*newpens == *oldpens) { skip = true; break; }
	}
	if(skip) { continue; }
    }
    ens.emplace_back(newpens);
  }
}

void chkens(vpEnsemble& ens) {
  uint Qt = 0;
  for(auto& pens : ens) {
    if(typeid(*pens) != typeid(*ens[0])) {
      throw(MDP_Exception("Input MDP files must refer to systems of the same ensemble"));
    }
    Qt |= pens->cmp(*ens[0]);
  }

  //! If pressures are consistent, down cast the pointer to NVT 
  if(!(Qt & (1<<Ensemble::DPressure))) {
    for(auto& pens : ens) { 
      auto _pens = dynamic_pointer_cast<NVT>(pens);
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
MDP2Ensemble (const vector<PMDP>& pmdps, const bool combinestates)
{
  vpEnsemble ans;
  for(const auto& pmdp : pmdps) {
    chkmdp(pmdp);   
    genens(pmdp, ans, combinestates);
  }
  //! Check if we need to take care of temperature and pressure
  chkens(ans);

  return ans;
}		/* -----  end of template function generate_ensemble  ----- */

#endif
