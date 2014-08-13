#if !defined(FUNCTOR_HPP)
#define FUNCTOR_HPP
/*
 * =====================================================================================
 *
 *       Filename:  functor.hpp
 *
 *    Description:  All kinds of functors
 *
 *        Version:  1.0
 *        Created:  31/07/14 16:04:54
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dejun Lin (DL), dejun.lin@gmail.com
 *   Organization:  Department of Biochemistry and Biophysics, Medical Center, University of Rochester
 *
 * =====================================================================================
 */
#include <vector>
#include <typeinfo>
#include <functional>
#include "typedefs.hpp"
#include "exception.hpp"

/*
 * =====================================================================================
 *        Class:  Functor
 *  Description:  Wrapper to all kinds of functors with basic comparison operator
 *                Note that the comparison only works if the rhs provides implementation 
 *                of getParams() member function
 * =====================================================================================
 */
template < class Params, class Output, class ... Input >
class Functor
{
  public:
    typedef function<Output(Input...)> Tf; 
    /* ====================  LIFECYCLE     ======================================= */
    //! Default constructor
    Functor() = default;
    //! Delegating constructor
    /** Only works if Tfunctor::getParams() is implemented
     */
    template < class Tfunctor > Functor (const Tfunctor& _f) : Functor(Tf(_f), _f.getParams())  {};

    //! Target constructor
    Functor (const Tf& _funct, const Params& _params) : funct(_funct), params(_params) {};
    //! Target constructor enabling move semantics in Params
    Functor (const Tf& _funct, const Params&& _params) : funct(_funct), params(move(_params)) {};

    virtual ~Functor() {};
    /* ====================  ACCESSORS     ======================================= */
    virtual const Params& getParams() const { return this->params; };

    /* ====================  MUTATORS      ======================================= */

    /* ====================  OPERATORS     ======================================= */
    virtual Output operator() (const Input&... in) const {
      return funct(in...);
    };

    template < class TheirType >
    bool operator== (const TheirType& rhs) const {
      try {
	const TheirType& lhs = dynamic_cast<const TheirType&>(*this);
	return &lhs == &rhs || lhs.getParams() == rhs.getParams();
      } catch (bad_cast& bcex) {
	return false;
      }
    };

    template < class TheirType >
    bool operator!= (const TheirType& rhs) const {
      return !this->operator==(rhs);
    };

  protected:
    /* ====================  METHODS       ======================================= */

    /* ====================  DATA MEMBERS  ======================================= */
    const Tf funct;
    const Params params;

  private:
    /* ====================  METHODS       ======================================= */

    /* ====================  DATA MEMBERS  ======================================= */

}; /* ----------  end of template class Functor  ---------- */

//! short-hand definition for scalar double -> scalar double function
typedef Functor<vector<valtype>, valtype, valtype> FunctVV;
typedef vector< FunctVV > vFunctVV; 

/*
 * =====================================================================================
 *        Class:  Quadratic
 *  Description:  1-dimensional quadratic functor
 * =====================================================================================
 */
class Quadratic 
{
  public:
    enum ParamEnum { K, R0, C };
    /* ====================  LIFECYCLE     ======================================= */
    Quadratic() : Quadratic(vector<valtype>({0,0,0})) {};
    Quadratic(const Quadratic& src) : Quadratic(src.getParams()) {};
    Quadratic (const valtype& _k, const valtype& _r0, const valtype& _c=0) :
      Quadratic(vector<valtype>({_k, _r0, _c}))
    {};

    Quadratic(const vector<valtype>&& _params) :
      params(move(_params)),
      k(params[K]), r0(params[R0]), c(params[C]) 
    {};

    Quadratic(const vector<valtype>& _params) :
      params(_params),
      k(params[K]), r0(params[R0]), c(params[C]) 
    {};

    virtual ~Quadratic() {};

    /* ====================  ACCESSORS     ======================================= */
    virtual vector<valtype> getParams() const {
      return params;
    };

    /* ====================  MUTATORS      ======================================= */

    /* ====================  OPERATORS     ======================================= */
    virtual valtype operator() (const valtype& r) const {
      const valtype dr = r - r0;
      return k*dr*dr + c;
    };

  protected:
    /* ====================  METHODS       ======================================= */

    /* ====================  DATA MEMBERS  ======================================= */
    const vector<valtype> params;
    const valtype& k, r0, c; 

  private:
    /* ====================  METHODS       ======================================= */

    /* ====================  DATA MEMBERS  ======================================= */

}; /* -----  end of class Quadratic  ----- */

/*
 * =====================================================================================
 *        Class:  QuadraticFlat
 *  Description:  1-dimensional Flat-bottom quadratic functor
 * =====================================================================================
 */
class QuadraticFlat 
{
  public:
    /* ====================  LIFECYCLE     ======================================= */
    QuadraticFlat (const valtype& _ref, 
                   const valtype& _r0, 
		   const valtype& _r1, 
		   const valtype& _k0, 
		   const valtype& _k1, 
		   const valtype& _c) :
                   ref(_ref),
                   r0(_r0),
		   r1(_r1),
		   k0(_k0),
		   k1(_k1),
		   c(_c) 
		   {};                             /* constructor */
    virtual ~QuadraticFlat() {};

    /* ====================  ACCESSORS     ======================================= */
    virtual vector<valtype> getParams() const {
      vector<valtype> ans;
      ans.push_back(ref);
      ans.push_back(r0);
      ans.push_back(r1);
      ans.push_back(k0);
      ans.push_back(k1);
      ans.push_back(c);
      return ans;
    };

    /* ====================  MUTATORS      ======================================= */

    /* ====================  OPERATORS     ======================================= */
    virtual valtype operator() (const valtype& r) const {
      const valtype dev = r - ref;
      if(dev >= r0 && dev <= r1) { return c; }
      else if(dev < r0) { 
	const valtype ldev = dev - r0;
	return k0*ldev*ldev + c;
      } else {
	const valtype rdev = dev - r1;
	return k1*rdev*rdev + c;
      }
    }

  protected:
    /* ====================  METHODS       ======================================= */

    /* ====================  DATA MEMBERS  ======================================= */
    const valtype ref, r0, r1, k0, k1, c;

  private:
    /* ====================  METHODS       ======================================= */

    /* ====================  DATA MEMBERS  ======================================= */

}; /* -----  end of class QuadraticFlat  ----- */

#endif
