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
#include <stdexcept>
#include <array>
#include <utility>
#include "metaprog_snippets.hpp"

using namespace std;

/*
 * =====================================================================================
 *        Class:  FunctorWrapper
 *  Description:  Wrapper to all kinds of functors with basic comparison operator
 *                Note that the comparison only works if the rhs provides implementation 
 *                of getParams() member function
 * =====================================================================================
 */
template < class Params, class Output, class ... Input >
class FunctorWrapper
{
  public:
    typedef function<Output(Input...)> Tf; 
    /* ====================  LIFECYCLE     ======================================= */
    //! Default constructor
    FunctorWrapper() = default;

    //! Delegating constructor
    /** Only works if Tfunctor::getParams() is implemented
     */
    template < class Tfunctor > FunctorWrapper (const Tfunctor& _f) : FunctorWrapper(Tf(_f), _f.getParams())  {};

    //! Target constructor
    FunctorWrapper (const Tf& _funct, const Params& _params) : funct(_funct), params(_params) {};

    //! Target constructor enabling move semantics in Params
    FunctorWrapper (const Tf& _funct, const Params&& _params) : funct(_funct), params(move(_params)) {};

    //! Destructor
    virtual ~FunctorWrapper() {};

    /* ====================  ACCESSORS     ======================================= */
    virtual const Params& getParams() const { return this->params; };

    /* ====================  MUTATORS      ======================================= */

    /* ====================  OPERATORS     ======================================= */
    //! calculate the function value
    virtual Output operator() (const Input... in) const {
      return funct(forward<Input...>(in...));
    };

    //! == with another object
    template < class TheirType >
    bool operator== (const TheirType& rhs) const {
      try {
	const TheirType& lhs = dynamic_cast<const TheirType&>(*this);
	return &lhs == &rhs || lhs.getParams() == rhs.getParams();
      } catch (bad_cast& bcex) {
	return false;
      }
    };

    //! =! with another object
    template < class TheirType >
    bool operator!= (const TheirType& rhs) const {
      return !this->operator==(rhs);
    };

  protected:
    /* ====================  METHODS       ======================================= */

    /* ====================  DATA MEMBERS  ======================================= */
    //! The virtual functor that does the caclulation
    const Tf funct;
    //! The parameters associated with the virtual functor
    const Params params;

  private:
    /* ====================  METHODS       ======================================= */

    /* ====================  DATA MEMBERS  ======================================= */

}; /* ----------  end of template class FunctorWrapper  ---------- */

//! short-hand definition for scalar double -> scalar double function
typedef FunctorWrapper<vector<valtype>, valtype, valtype&> FunctVV;

//! short-hand for a vector of scalar functors
typedef vector< FunctVV > vFunctVV; 

/*
 * =====================================================================================
 *        Class:  Functor
 *  Description:  Generic functor class where the user can just plug in a custom 
 *                function (or function pointer) and can treat the arguments as 
 *                either parameters or arguments for this->operator()
 * =====================================================================================
 */
//! Primary template
template < class ... T > class Functor; 
//! Partially specialized for a functor using compiler-time fixed-size container as parameter holder
template < 
	   template <class, class...> class Funct,
	   class Output,
	   class ... Arg,
           template <class,size_t> class ParamsContainer,
           class ParamsContained,
	   class ParamsOutput,
	   class ... Params
         >
class Functor < 
                Funct<Output(Arg...)>,
                ParamsContainer<ParamsContained, sizeof...(Params)>, 
	        ParamsOutput,
	        Params... 
	      >
{
  public:
    //! functor type
    typedef Funct<Output(Arg...)> Tf;
    //! Number of parameters
    constexpr static uint Nparams = sizeof...(Params);
    //! parameters type
    typedef ParamsContainer<ParamsContained, Nparams> Tp;
    //! Here we make a index sequence for the parameters
    typedef ArithmeticSeq<IndexSeq<>, IndexSeq<> > IndexSeqGen;
    typedef StopAtMaxN<IndexSeq<>, IndexSeq<Nparams> > IndexSeqStop;
    typedef typename make_index_seq<IndexSeqStop::decision, IndexSeq<>,  IndexSeqGen, IndexSeqStop>::type ParamIndex;

    /* ====================  LIFECYCLE     ======================================= */
    //! Construct from a functor and a list of parameters
    Functor (const Tf& _funct, const Params... params) : Functor(_funct, Tp{params...}, ParamIndex()) {};

    /* ====================  ACCESSORS     ======================================= */
    //! Get the parameters
    const Tp& _getParams() const { return params; };
    //! Output the parameters (to other wrapper/interface) 
    const ParamsOutput& getParams() const { return paramsout; };

    /* ====================  OPERATORS     ======================================= */
    //! Evaluate the functor
    /** Note that we really need perfect forwarding here because 'arg' is 
     * lvalue in this function body but Tf::operator() could take rvalue
     */
    Output operator()(Arg... arg) const { return funct(forward<Arg...>(arg...)); };
  protected:
    /* ====================  METHODS       ======================================= */
    //! Construct from a functor and a container of parameters -- Target constructor
    /** Note the last argument has to be provided at the called site since the 
     * template parameters can't be deduced from default function argument
     */
    template < size_t ... indices >
    Functor (const Tf& _funct, const Tp& _params, const IndexSeq<indices...> paramid) : 
      funct(_funct), params(_params), paramsout{params[indices]...} {};

    /* ====================  DATA MEMBERS  ======================================= */
    const Tf funct;
    Tp params;
    const ParamsOutput paramsout;
}; /* ----------  end of template class Functor  ---------- */

/*
 * =====================================================================================
 *        Class:  Quadratic
 *  Description:  1-dimensional quadratic functor
 * =====================================================================================
 */
class Quadratic 
{
  public:
    //! The elements in the Quadratic::params array
    enum ParamEnum { K, R0, C, NPARAM };
    typedef array<valtype, NPARAM> Params;
    /* ====================  LIFECYCLE     ======================================= */
    //! construct an empty object
    Quadratic() : Quadratic(0,0,0) {};

    //! copy constructor
    Quadratic(const Quadratic& src) : Quadratic(src._getParams()) {};

    //! construct from a list of parameters
    Quadratic (const valtype& _k, const valtype& _r0, const valtype& _c=0) :
      Quadratic(Params{_k, _r0, _c})
    {};
    
    //! construct from a vector (rvalue)
    Quadratic(const Params&& _params) : 
      params(move(_params)),
      k(params.at(K)), r0(params.at(R0)), c(params.at(C)) 
    {};

    //! consttruct from a vector (lvalue)
    Quadratic(const Params& _params) :
      params(_params),
      k(params.at(K)), r0(params.at(R0)), c(params.at(C)) 
    {};

    //! destructor
    virtual ~Quadratic() {};

    /* ====================  ACCESSORS     ======================================= */
    //! return the parameters (vector interface)
    virtual vector<valtype> getParams() const { return vector<valtype>(params.begin(), params.end()); };
    //! return the parameters (array interface)
    virtual const Params& _getParams() const { return params; }

    /* ====================  MUTATORS      ======================================= */

    /* ====================  OPERATORS     ======================================= */
    //! calculate the value of the quadratic function
    virtual valtype operator() (const valtype& r) const {
      const valtype dr = r - r0;
      return k*dr*dr + c;
    };

  protected:
    /* ====================  METHODS       ======================================= */

    /* ====================  DATA MEMBERS  ======================================= */
    //! a vector of parameters
    const Params params;
    //! reference to the elements in params
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
