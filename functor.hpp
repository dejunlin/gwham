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
#include "typedefs.hpp"
#include "exception.hpp"

/*
 * =====================================================================================
 *        Class:  Functor
 *  Description:  Interfacial class for all kinds of functors
 * =====================================================================================
 */
template < class Input, class Output >
class Functor
{
  public:
    /* ====================  LIFECYCLE     ======================================= */
    Functor () {};                             /* constructor */

    /* ====================  ACCESSORS     ======================================= */
    virtual vector<Input> getParams() const = 0;

    /* ====================  MUTATORS      ======================================= */

    /* ====================  OPERATORS     ======================================= */
    virtual Output operator() (const Input& in) const = 0;

    template < class TheirType >
    bool operator== (const TheirType& rhs) const {
      try {
	const TheirType& lhs = dynamic_cast<const TheirType&>(*this);
	return &lhs == &rhs || lhs.getParams() == rhs.getParams();
      } catch (bad_cast& bcex) {
	return false;
      }
    };

    virtual bool operator!= (const Functor<Input, Output>& rhs) const {
      return !this->operator==(rhs);
    };
    virtual ~Functor() {};

  protected:
    /* ====================  METHODS       ======================================= */

    /* ====================  DATA MEMBERS  ======================================= */

  private:
    /* ====================  METHODS       ======================================= */

    /* ====================  DATA MEMBERS  ======================================= */

}; /* ----------  end of template class Functor  ---------- */

//! short-hand definition for a vector of functors
/** We can use smart-pointer here if the performance 
 * is good
 */
typedef const Functor<valtype, valtype>* cpFunctorVV;
typedef vector< cpFunctorVV > vPFunct; 

/*
 * =====================================================================================
 *        Class:  Quadratic
 *  Description:  1-dimensional quadratic functor
 * =====================================================================================
 */
class Quadratic : public Functor<valtype, valtype> 
{
  public:
    /* ====================  LIFECYCLE     ======================================= */
    Quadratic() : k(0), r0(0), c(0) {};
    Quadratic (const valtype _k, const valtype _r0, const valtype _c=0) : k(_k), r0(_r0), c(_c) {};                             /* constructor */
    Quadratic(const Quadratic& src) :
      k(src.getParams()[0]),
      r0(src.getParams()[1]),
      c(src.getParams()[2]) 
    {};
    virtual ~Quadratic() {};

    /* ====================  ACCESSORS     ======================================= */
    virtual vector<valtype> getParams() const {
      vector<valtype> ans;
      ans.push_back(k);
      ans.push_back(r0);
      ans.push_back(c);
      return ans;
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
    const valtype k, r0, c; 

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
class QuadraticFlat : public Functor<valtype, valtype>
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
