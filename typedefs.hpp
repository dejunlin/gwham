#if !defined(TYPEDEFS_HPP)
#define TYPEDEFS_HPP

#include <vector>
#include <math.h>
#include <float.h>
#include <limits>
#include <climits>
#include <bitset>
#include <memory>
#include "functor.hpp"

using namespace std;

typedef unsigned int uint;
typedef unsigned long ulong;

//! Number of bits in unsigned int
static const uint Nbits_UINT = numeric_limits<uint>::digits; 
//! Number of bits in unsigned long
static const uint Nbits_ULONG = numeric_limits<ulong>::digits;
//! Max number of restraint potential
static const uint MAXNRST = 1000;
//! Max number of digits for number of windows
static const uint MAXNDIGWIN = 50;
//! Max number of digits for number of run-cycles
static const uint MAXNDIGRUN = 50;
//!just a unit bit
const bitset<MAXNRST> bitunit = bitset<MAXNRST>(1);

typedef vector<uint> coordtype; //type of the coordinate for histogram and narray
typedef double valtype; //the data type we histogram
typedef ulong histcounter; //the histogram counter
typedef ulong linecounter; //the file line counter
static const ulong MAXNLINE = ULONG_MAX; //maximum number of lines that can be counted 

static constexpr enable_if<numeric_limits<valtype>::has_quiet_NaN, valtype>::type
NaN = numeric_limits<valtype>::quiet_NaN();

const valtype MAXEXPARG = log(DBL_MAX); //max of the arguments to std::exp() call
const valtype MINEXPARG = log(DBL_MIN); //min of the arguments to std::exp() call

//!kB in kJ/mol/K
const valtype BoltzmannkJ = 0.0083144621; 
const valtype Boltzmannkcal = 0.0019872041;

//! short-hand definition for scalar double -> scalar double function
typedef Functor<const function<valtype(const valtype&)>, const vector<valtype>, const void> FunctVV;

//! short-hand for a vector of scalar functors
typedef vector< FunctVV > vFunctVV;

//! pointer to MDP
class MDP;
typedef shared_ptr<MDP> pMDP;
typedef vector<pMDP> vpMDP;

//! pointer to Ensemble
class Ensemble;
class NVE;
class NVT;
class NPT;
typedef shared_ptr<Ensemble> pEnsemble;
typedef vector<pEnsemble> vpEnsemble;
typedef shared_ptr<NVE> pNVE;
typedef shared_ptr<NVT> pNVT;
typedef shared_ptr<NPT> pNPT;

#endif
