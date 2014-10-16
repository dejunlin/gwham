#ifndef MDP_FACTORY_HPP
#define MDP_FACTORY_HPP
/*
 * =====================================================================================
 *
 *       Filename:  mdp_factory.hpp
 *
 *    Description:  MDP factory function
 *
 *        Version:  1.0
 *        Created:  12/09/14 10:05:43
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dejun Lin (DL), dejun.lin@gmail.com
 *   Organization:  Department of Biochemistry and Biophysics, Medical Center, University of Rochester
 *
 * =====================================================================================
 */
#include <map>
#include <functional>
#include "typedefs.hpp"
#include "mdp.hpp"
#include "gmxmdp.hpp"
#include "fileio.hpp"

using namespace std;

typedef map<string, function<pMDP(const string&)>> Tsuffix2MDP;
static const Tsuffix2MDP suffix2MDP {
  {"mdp", GMXMDP::CreateMDP},
};

//! This is the MDP factory function that creates pointers to new 
// MDP instances based on the combination of file name prefixes and 
// suffixes. The supported suffixes must exist in the map type 
// suffix2mdp
template < class PMDP >
typename std::enable_if<is_pointer<PMDP>::value, vector<PMDP>>::type
CreateMDPs(const vector<string>& MDPprefixes,
           const vector<string>& MDPsuffixes,
           const map<string, function<PMDP(const string&)>>& suffix2mdp) 
{
  using Tmap = const map<string, function<PMDP(const string&)>>;
  using Tmapciterator = typename Tmap::const_iterator;
  //Do sanity check to see if all the input types of MDP supported
  for(const auto& MDPsuffix : MDPsuffixes) {
    const Tmapciterator it = suffix2mdp.find(MDPsuffix);
    if(it == suffix2mdp.end()) {
      throw(MDP_Factory_Exception("MDP file with suffix '"+MDPsuffix+"' not supported"));
    }
  }

  vector<PMDP> ans;
  //Here we loop over all the combination of prefixes and suffixes and 
  //create MDPs
  for(const auto& MDPprefix : MDPprefixes) {
    for(const auto& MDPsuffix : MDPsuffixes) {
      const string fmdp = MDPprefix+"."+MDPsuffix;
      try {
	// check if the file exists
	// if not, a FILEIO_Exception will be fired
        fileio fio(fmdp, fstream::in, false, 0, 1, MAXNLINE, ";#@");
      } catch (FILEIO_Exception& fioex) {
	continue;
      }
      ans.emplace_back( suffix2mdp.find(MDPsuffix)->second(fmdp) );
    }
  }
  return ans;
}

#endif
