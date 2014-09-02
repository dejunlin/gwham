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

template <class MDP>
vpEnsemble MDP2Ensemble (const vector<MDP>& mdps, const bool combinestates)
{
  vpEnsemble ans;
  for(const auto& mdp : mdps) {
    //For simplicity, we only support restraints for now
    //TODO: to support FEP, all we need to do is add more 
    //potential energy functor into the Hamiltonian in the
    //lines followed
    if(mdp.hasFEPLambda() && ( 
       mdp.hasLbond() || mdp.hasLmass() || mdp.hasLvdw() || 
       mdp.hasLcoul() || mdp.hasLtemp() )
      ) {
      throw(MDP_Exception("Can't use GWHAM to handle FEP for now"));
    }

    //In case of expanded ensemble where a simulation can jump 
    //between different restraint potential, there should be 
    //Nstates functors defined for each pull-group and
    //mdp.rstfuncts[i][j] is the j'th functor of the i'th pull-group
    const vector<vFunctVV>& pgrprsts = mdp.getrstfuncts();
    uint Nstates = 0;
    for(const auto& pgrprst : pgrprsts) {
      if(Nstates && Nstates != pgrprst.size()) {  
	throw(MDP_Exception("The number of rstfunctors for each pull-group are not consistent"));
      };
      Nstates = pgrprst.size();
    }
    //Here rstfuncts[i][j] is the j'th functor of the i'th ensemble
    vector<vFunctVV> rstfuncts(Nstates);
    
    for(const auto& pgrprst : pgrprsts) {
      for(uint i = 0; i < pgrprst.size(); ++i) {
	rstfuncts[i].emplace_back(pgrprst[i]);
      }
    }

    //TODO: when we support expanded ensemble on temperature and/or 
    //pressure space, we need to extract the vector of temperatures
    //and states out of the MDP object and construct vectors of 
    //temperature/pressure corresponding to each state
    vector<valtype> Ts(Nstates, mdp.getTemperature());
    vector<valtype> Ps(Nstates, mdp.getPressure());
    for(uint i = 0; i < Nstates; ++i) {
      auto newpens = 
          make_shared<NPT>(mdp.getkB(), 
                           Hamiltonian(rstfuncts[i]),
			   Ts[i],
			   Ps[i]
                          );
      //if true, we'll check if this new ensemble is equivalent to 
      //the ensemble we've already created -- if yes, we skip pushing
      //it into ans
      //TODO: in case we have redundant ensembles, we need to create 
      //a map to tell which mdp coresponds to which ensembles
      if(combinestates) {
	bool skip = false;
	for(const auto& oldpens : ans) {
	  if(*newpens == *oldpens) { skip = true; break; }
	}
	if(skip) { continue; }
      }
      ans.emplace_back(newpens);
    }
  }
  //! Check if we need to take care of temperature and pressure
  uint Qt = 0;
  for(auto& pens : ans) {
    Qt |= pens->cmp(*ans[0]);
  }
/*   bitset<numeric_limits<uint>::digits> bits(Qt);
 *   cout << bits << endl;
 */

  //! If pressures are consistent, down cast the pointer to NVT 
  if(!(Qt & (1<<Ensemble::DPressure))) {
    for(auto& pens : ans) { 
      auto _pens = dynamic_pointer_cast<NVT>(pens);
      pens.reset();
      pens = make_shared<NVT>(_pens->getkB(), _pens->getH(), _pens->getT()); 
    }
  }
  //! If temperautre are diff., add potential energy
  if(Qt & (1<<Ensemble::DTemperature)) {
    for(auto& pens : ans) { 
      auto _pens = dynamic_pointer_cast<NVE>(pens);
      _pens->getH().getPotentialFuncts().emplace_back(Make_Functor_Wrapper(FunctVV{}, Linear<valtype>, 1.0, 0.0)); 
    }
  } 
  
  cout << "Constructed " << ans.size() << " ensembles\n";

  return ans;
}		/* -----  end of template function generate_ensemble  ----- */


#endif
