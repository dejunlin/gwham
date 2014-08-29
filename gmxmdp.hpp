#if !defined(GMXMDP_HPP)
#define GMXMDP_HPP
/*
 * =====================================================================================
 *
 *       Filename:  gmxmdp.hpp
 *
 *    Description:  GROMACS mdp file reader 
 *
 *        Version:  1.0
 *        Created:  03/07/14 16:33:51
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dejun Lin (DL), dejun.lin@gmail.com
 *   Organization:  Department of Biochemistry and Biophysics, Medical Center, University of Rochester
 *
 * =====================================================================================
 */

#include "typedefs.hpp"
#include "fileio.hpp"
#include "fileio_utils.hpp"
#include "exception.hpp"
#include "mdp.hpp"
#include "functor.hpp"
#include <string>
#include <string.h>
#include <map>
#include <vector>
#include <iterator>

using namespace std;
/*
 * =====================================================================================
 *        Class:  GMXMDP
 *  Description:  
 * =====================================================================================
 */
class GMXMDP : public MDP
{
  public:
    class GMXPGRP;
    /* ====================  LIFECYCLE     ======================================= */
    GMXMDP (const string& fname);    /* constructor */

    /* ====================  ACCESSORS     ======================================= */
    //! print out all the parameters read 
    virtual void print() const; 
    //! compare 2 MDP objects and return a mask indicating what thermodynamic quantities are different
    virtual uint cmp(const MDP& mdp) const;

    /* ====================  MUTATORS      ======================================= */

    /* ====================  OPERATORS     ======================================= */
    

    /* ====================  METHODS       ======================================= */

   /* ====================  DATA MEMBERS  ======================================= */
    //! FEP type
    static const map<string, FEPType> str2fepT; 
    string fepTstr = "";
    //! Pull type
    enum PullType {
      No, Umbrella, UmbrellaFlatBottom, Constraint, ConstantForce, Contact, NPullTypes
    } pullT = NPullTypes;
    static const map<string, PullType> str2pullT; 
    string pullTstr = "";
    //! Pull geometry
    enum PullGeom {
      Distance, Direction, DirectionPeriodic, Cylinder, Position, NGeomTypes
    } geomT = NGeomTypes;
    static const map<string, PullGeom> str2geomT;
    string geomTstr = "";
    //! Output frequency in steps for dhdl file
    int nstdhdl = -1;
    //! Output frequency in steps for doing expanded ensemble
    int nstexpanded = -1;
    //! Temperature for doing Monte Carlo jump between FEP states 
    valtype Tmc = T;
    //! Contact lambdas
    vector<valtype> Lcnt1, Lcnt2, Lcnt3;
    //! Pull dimension
    vector<string> dimstr;
    static const map<string, bool> str2dim;
    //! Pull dimension
    int dim = 0;
    //! Output frequency in steps for pulled quantity 
    int nstx = -1;
    //! Output frequency in steps for pulling forces
    int nstf = -1;
    //! Number of pull groups
    int npgrps = -1;
    //! Number of contact groups
    int ncntgrps = -1;
    //! Pull groups
    vector<GMXPGRP> pgrps;

    template < class T > using rw = reference_wrapper<T>;
    typedef rw<valtype> valref;
    typedef rw<int> intref;
    typedef rw<string> strref;
    typedef rw<vector<valtype> > vvref;
    typedef rw<vector<string> > vsref;
    typedef map<string, valref> dblopts;
    typedef map<string, intref> intopts;
    typedef map<string, strref> stropts;
    typedef map<string, vvref> vecopts;
    typedef map<string, vsref> vecstropts;

    //! double-type options
    dblopts key2val;  
    //! string-type options
    stropts key2str; 
    //! vector<string>-type options -- usually convertible to masks
    vecstropts key2vecstr;
    //! int-type options
    intopts key2int;
    //! vector<double>-type options
    vecopts key2vvec;
    //! if the MDP key maps have been initialized
    const bool ifinit;
  private:
    /* ====================  METHODS       ======================================= */
    //! populate the pgrps vector
    void setpgrps();
    //! convert string into enum
    void setenum();
    //! convert string into mask
    void setmask();
    //! integrate the input info. and make final processing
    virtual void doublechk();
    //! check temperature, pressure ... etc.
    void checkgeneric();
    //! check FEP parameters
    void checkfep();
    //! check pull parameters
    void checkpull();
    //! convert contact parameters
    void setcntgrp();
    //! set GMXPGRP::Lrst
    void setpgrpLrst();
    //! we treate each dimension of pulling as an individual pull group
    void expandpgrps();
    //! populate GMXPGRP::rstfuncts
    void setrstfuncts();
    //! Initialize the MDP key maps
    bool initopts() {
      key2val =
      {
        {"ref-t", T},
        {"ref_t", T},
        {"ref-p", P},
        {"ref_p", P},
        {"mc-temperature", Tmc},
      };  
      
      key2str = 
      {
        {"free-energy", fepTstr},
        {"pull", pullTstr},
        {"pull_geometry", geomTstr},
        {"pull-geometry", geomTstr},
      };
  
      key2vecstr = 
      {
        {"pull_dim", dimstr},
        {"pull-dim", dimstr},
      };
  
      key2int =
      {
        {"init-lambda-state", Linit},
        {"nstdhdl", nstdhdl},
        {"nstexpanded", nstexpanded},
        {"pull_nstxout", nstx},
        {"pull-nstxout", nstx},
        {"pull_nstfout", nstf},
        {"pull-nstfout", nstf},
        {"pull-ngroups", npgrps},
        {"pull_ngroups", npgrps},
        {"pull_ncontactgroups", ncntgrps},
      };
  
      key2vvec =
      {
        {"bonded-lambdas", Lbond},
        {"mass-lambdas", Lmass},
        {"vdw-lambdas", Lvdw},
        {"coul-lambdas", Lcoul},
        {"restraint-lambdas", Lrst},
        {"temperature-lambdas", Ltemp},
        {"umbrella-contact1-lambdas", Lcnt1},
        {"umbrella-contact2-lambdas", Lcnt2},
        {"umbrella-contact3-lambdas", Lcnt3},
      };
      return true;
    };

    /* ====================  DATA MEMBERS  ======================================= */
}; /* ----------  end of template class GMXMDP  ---------- */

class GMXMDP::GMXPGRP {
 public:
  //! target constructor
  GMXPGRP(const uint& _gid) : 
    gid(_gid), gidstr(tostr(gid)),
    ifinit(this->initopts())
  {};
  //! copy constructor
  GMXPGRP(const GMXPGRP& src) :
    GMXPGRP(src.gid)
  {
    this->operator=(src);
  };
  //! copy assignment
  GMXPGRP& operator=(const GMXPGRP& src) {
    gid = src.gid;
    gidstr = src.gidstr;
    rstT = src.rstT;
    init = src.init;
    initB = src.initB;
    vec = src.vec;
    k = src.k;
    kB = src.kB;
    r0 = src.r0;
    r0B = src.r0B;
    r1 = src.r1;
    r1B = src.r1B;
    k0 = src.k0;
    k0B = src.k0B;
    k1 = src.k1;
    k1B = src.k1B;
    nc0 = src.nc0;
    nc0B = src.nc0B;
    Lrst = src.Lrst;

    rstfuncts.clear();
    for(const FunctVV& rstfunct: src.rstfuncts) {
      rstfuncts.emplace_back(rstfunct);
    }

    return *this;
  }
  //! Pull group id (for pull-group, this is zero initial; for contact, 1 initial)
  uint gid;
  //! Pull group id (for pull-group, this is zero initial; for contact, 1 initial)
  string gidstr;
  //! Reference position in state A and B
  vector<valtype> init{}, initB{};
  //! Velocity of the reference position
  vector<valtype> vec;
  //! Coefficient of the harmonic potential in state A and B
  valtype k = 0, kB = 0;
  //! Parameters for the flat-bottom harmonic potential
  valtype r0 = 0, r0B = r0, r1 = 0, r1B = r1, k0 = 0, k0B = k0, k1 = k0, k1B = k1;
  //! Reference contact number
  valtype nc0 = 0, nc0B = 0;
  //! Restraint lambdas
  vector<valtype> Lrst;
  //! Restraint type
  enum RestraintType { Quad, QuadFlat, NRestraintTypes } rstT = NRestraintTypes;
  //! double-type options
  dblopts key2val; 
  //! vector<double>-type options
  vecopts key2vvec;
  //! Restraint functors
  vFunctVV rstfuncts; 
  //! if all the option maps are initialized
  const bool ifinit;
 private:
  //! Initialize the MDP key maps
  bool initopts() {
    key2val =  
    {
      {"pull_k" + gidstr, k},
      {"pull_kB" + gidstr, kB},
      {"pull-k" + gidstr, k},
      {"pull-kB" + gidstr, kB},
      {"pull_umb_r0" + gidstr, r0},
      {"pull_umb_r0B" + gidstr, r0B},
      {"pull_umb_r1" + gidstr, r1},
      {"pull_umb_r1B" + gidstr, r1B},
      {"pull_umb_k0" + gidstr, k0},
      {"pull_umb_k0B" + gidstr, k0B},
      {"pull_umb_k1" + gidstr, k1},
      {"pull_umb_k1B" + gidstr, k1B},
      {"pull_contactgroup" + gidstr + "_k", k},
      {"pull_contactgroup" + gidstr + "_kB", kB},
    };
    key2vvec =  
    {
      {"pull_vec" + gidstr, vec},
      {"pull_init" + gidstr, init},
      {"pull_initB" + gidstr, initB},
      {"pull-vec" + gidstr, vec},
      {"pull-init" + gidstr, init},
      {"pull-initB" + gidstr, initB},
      {"pull_contactgroup" + gidstr + "_nc0", init},
      {"pull_contactgroup" + gidstr + "_nc0B", initB},
    };
    return true;
  }
};

//! This function handles MDP option with option value that 
//can be interpreted as a type T which has 
//member function void emplace_back(T::value_type&)
//which usually means the type T is a dynamic-size container
template < class T >
typename enable_if<check_if<T, has_memfn_emplace_back>::value, void>::type 
setMDPoptval(const string& val, T& container) {
  typedef typename T::value_type value_type;
  container.emplace_back(stosc<value_type>(val));
}

//! This function set MDP option for scalar value 
template < class T >
typename enable_if<std::is_scalar<T>::value || is_string<T>::value, void>::type 
setMDPoptval(const string& val, T& scalar) {
  setsc(val, scalar);
}

template < class InIterator, class T >
bool setMDPopt(const InIterator& itbegin, const InIterator& itend, map<string,reference_wrapper<T>>& key2opt) {
  typedef typename map<string,reference_wrapper<T>>::iterator OIterator;
  const string key(*itbegin);
  OIterator oit = key2opt.find(key);
  if(oit == key2opt.end()) { return false; }
  T& optdata = oit->second.get();
  InIterator iit = itbegin;
  for(++iit; iit != itend; ++iit) {
    setMDPoptval(string(*iit), optdata);
  }
  return true;
}

template < class T >
typename enable_if<check_if<T, can_be_streamed<ostream>>::value, void>::type
printMDPopt(const map<string, reference_wrapper<T>>& key2opt) {
  typename map<string, reference_wrapper<T>>::const_iterator it = key2opt.begin();
  for(; it != key2opt.end(); ++it) {
    const string& key = it->first;
    const T& data = it->second.get();
    printf ("#%27s = ", key.c_str());
    cout << data << endl;
  }
}

template < class T >
typename enable_if<std::is_enum<T>::value, void>::type
setMDPenum(const map<string, T>& str2enum, const string& str, T& optval) {
  const typename map<string,T>::const_iterator it = str2enum.find(str);
  if(it != str2enum.end()) { optval = it->second; }
};

#endif
