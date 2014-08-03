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

    /* ====================  MUTATORS      ======================================= */

    /* ====================  OPERATORS     ======================================= */
    virtual Output operator() (const Input& in) const = 0;
    virtual bool operator== (const Functor<Input, Output>& rhs) const = 0;
    virtual bool operator!= (const Functor<Input, Output>& rhs) const {
      return !this->operator==(rhs);
    };
    virtual ~Functor() = 0;

  protected:
    /* ====================  METHODS       ======================================= */

    /* ====================  DATA MEMBERS  ======================================= */

  private:
    /* ====================  METHODS       ======================================= */

    /* ====================  DATA MEMBERS  ======================================= */

}; /* ----------  end of template class Functor  ---------- */

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

    virtual bool operator== (const Functor<valtype, valtype>& rhs) const {
      if(this == &rhs) { return true; }
      try {
	const Quadratic& rrhs = dynamic_cast<const Quadratic&>(rhs);
        return this->getParams() == rrhs.getParams();
      } catch (bad_cast& bcex) {
	return false;
      }
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
 *        Class:  QuadraticNwrap
 *  Description:  Wrapper functor of N quadratic
 * =====================================================================================
 */
//class QuadraticNwrap : Functor<vector<valtype>, valtype>
//{
//  public:
//    /* ====================  LIFECYCLE     ======================================= */
//    QuadraticNwrap (const vector<valtype>& _ks, const vector<valtype>& _r0s, const vector<valtype>& _cs) throw(Functor_Exception)
//    : Qfuncts(vector<Quadratic>(0))  
//    {
//      if(_ks.size() != _r0s.size() || _r0s.size() != _cs.size()) {
//	throw(Functor_Exception("Size of QuadraticNwrap parameters don't match each other"));
//      }
//      for(uint i = 0; i < _ks.size(); ++i) {
//	Qfuncts.push_back(Quadratic(_ks[i], _r0s[i], _cs[i]));
//      }
//    };                             /* constructor */
//    QuadraticNwrap (const vector<valtype>& _ks, const vector<valtype>& _r0s) throw(Functor_Exception)
//    : Qfuncts(vector<Quadratic>(0))  
//    {
//      if(_ks.size() != _r0s.size()) {
//	throw(Functor_Exception("Sizes of QuadraticNwrap parameters don't match each other"));
//      }
//      for(uint i = 0; i < _ks.size(); ++i) {
//	Qfuncts.push_back(Quadratic(_ks[i], _r0s[i], 0));
//      }
//    };                             /* constructor */
//    virtual ~QuadraticNwrap() {};
//
//    /* ====================  ACCESSORS     ======================================= */
//    const vector<Quadratic>& getQfuncts() const {
//      return this->Qfuncts;
//    }
//    /* ====================  MUTATORS      ======================================= */
//
//    /* ====================  OPERATORS     ======================================= */
//    virtual valtype operator() (const vector<valtype>& r) const throw(Functor_Exception) {
//      if(r.size() != Qfuncts.size()) {
//	throw(Functor_Exception("Size of QuadraticNwrap arguments doesn't match the number of functors"));
//      }
//      valtype ans = 0.0;
//      for(uint i = 0; i < r.size(); ++i) {
//	ans += Qfuncts[i](r[i]);
//      }
//      return ans;
//    };
//
//    virtual bool operator== (const Functor<vector<valtype>, valtype>& rhs) const {
//      if(this == &rhs) { return true; }
//      try {
//	const QuadraticNwrap& rrhs = dynamic_cast<const QuadraticNwrap&>(rhs);
//        return this->getQfuncts() == rrhs.getQfuncts();
//      } catch (bad_cast& bcex) {
//	return false;
//      }
//    };
//  protected:
//    /* ====================  METHODS       ======================================= */
//
//    /* ====================  DATA MEMBERS  ======================================= */
//    vector<Quadratic> Qfuncts;
//
//  private:
//    /* ====================  METHODS       ======================================= */
//
//    /* ====================  DATA MEMBERS  ======================================= */
//
//}; /* -----  end of class QuadraticNwrap  ----- */

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

    virtual bool operator== (const Functor<valtype, valtype>& rhs) const {
      if(this == &rhs) { return true; }
      try {
	const QuadraticFlat& rrhs = dynamic_cast<const QuadraticFlat&>(rhs);
        return this->getParams() == rrhs.getParams();
      } catch (bad_cast& bcex) {
	return false;
      }
    };
  protected:
    /* ====================  METHODS       ======================================= */

    /* ====================  DATA MEMBERS  ======================================= */
    const valtype ref, r0, r1, k0, k1, c;

  private:
    /* ====================  METHODS       ======================================= */

    /* ====================  DATA MEMBERS  ======================================= */

}; /* -----  end of class QuadraticFlat  ----- */

#endif
