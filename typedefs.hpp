#if !defined(TYPEDEFS_HPP)
#define TYPEDEFS_HPP

#include <vector>
#include <math.h>
#include <float.h>
#include <limits>
#include <bitset>

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

const valtype MAXEXPARG = log(DBL_MAX); //max of the arguments to std::exp() call
const valtype MINEXPARG = log(DBL_MIN); //min of the arguments to std::exp() call

#endif
