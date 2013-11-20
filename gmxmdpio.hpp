#if !defined(GMXMDPIO_HPP)
#define GMXMDPIO_HPP

#include "typedefs.hpp"
#include "gmxpullpot.hpp"
#include <fstream>
#include <map>

using namespace std;

//! read gromacs mdp file and construct the restraint potential functor from the pull group potential
class mdp2pullpot {
  public:
    mdp2pullpot() {};
    linecounter operator() (fstream& fs, map<uint,umbrella>& funct) const;
    linecounter operator() (fstream& fs, map<uint,umbrella_fb>& funct) const;
  private:
    //! Given a directive keyword in the COM-pull section in mdp file, extract the pullgrp id
    /** 'str' is of the syntax 'pull_.*[0-9]'
     */
    uint getpullgrpid(const string& directive, const string& str) const;
};

#endif
