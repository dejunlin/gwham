#if !defined(GMXMDPIO_HPP)
#define GMXMDPIO_HPP

#include "typedefs.hpp"
#include "gmxpullpot.hpp"
#include "hamiltonian.hpp"
#include "ensemble.hpp"
#include <string>
#include <map>
#include "fileio.hpp"

using namespace std;

//! read gromacs mdp file and construct the restraint potential functor from the pull group potential
/** This class also instructs what other gromacs file classes will do
 */
class mdp2pullpot {
  public:
    mdp2pullpot(const bitset<MAXNRST>& _rcmask, const fileio& _fio);
    linecounter operator() (const string& fname, map<uint,vector<umbrella*> >& funct, vector<Hamiltonian<RSTXLAMBDAsgl>* >& V);
    linecounter operator() (const string& fname, map<uint,vector<umbrella_fb*> >& funct, vector<Hamiltonian<RST_fbXLAMBDAsgl>* >& V);
  private:
    //! Given a directive keyword in the COM-pull section in mdp file, extract the pullgrp id
    /** 'str' is of the syntax 'pull_.*[0-9]'
     */
    uint getpullgrpid(const string& directive, const string& str) const;
    uint getcontactgrpid(const string& str) const;
    //! Reaction Coordinate mask that indicates what pull groups are interesting
    const bitset<MAXNRST> rcmask;
    fileio fio;
};

#endif
