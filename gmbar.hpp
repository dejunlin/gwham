#if !defined(GMBAR_HPP)
#define GMBAR_HPP

#include "typedefs.hpp"
#include <vector>

using namespace std;

//! MBAR equation as in J. Chem. Phys. 129, 124105, 2008 equation 11

class GMBAR {
  public:
    GMBAR(const vector<vector<vector<valtype> > >& dEs, const vector<uint>& N, const vector<valtype>& _f, const valtype& _tol);
    const vector<valtype>& getf() const;
  private:
    static const valtype kB = 0.0083144621;
    static const valtype T = 300;
    //!max number of iterations
    static const ulong MAXIT = 1000000;
    //! free energy of each state
    vector<valtype> f;
    //! tolerance in the iteration
    const valtype tol;
    //! Check if the iteration should end; also update GMBAR::f if not end
    bool endit(const vector<valtype>& newf, const ulong& count);
    //!shift all the free energies (argument newf)by a constant so that newf[0] is always zero 
    void shiftf(vector<valtype>& newf) const;
    //! print out all free energies
    void printfree() const;
    //! calculate newf based on dEs array
    void calnewf(const vector<vector<vector<valtype> > >& dEs, const vector<uint>& N, vector<valtype>& newf) const;
};

#endif
