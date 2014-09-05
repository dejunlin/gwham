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
#include <array>
#include <utility>
#include "metaprog_snippets.hpp"

using namespace std;

/*
 * =====================================================================================
 *        Class:  Functor
 *  Description:  Generic functor class that provide access to the functor parameters
 *                Users can just plug in a custom function (or function pointer) and 
 *                can treat the arguments as either parameters or functor arguments 
 * =====================================================================================
 */
//! Primary template
template < class ... T > class Functor; 
//! Partially specialized for a functor using compiler-time fixed-size container as parameter holder
template < 
	   template <class, class...> class Funct,
	   class Output,
	   class ... Arg,
           template <class...> class ParamsContainer,
           class ParamsContained,
	   class ... ParamsContainerOther,
	   class ... Params
         >
class Functor < 
                const Funct<Output(Arg...)>,
                const ParamsContainer<ParamsContained, ParamsContainerOther...>, 
	        Params... 
	      >
{
  public:
    //! functor type
    typedef const Funct<Output(Arg...)> Tf;
    //! parameters type
    typedef const ParamsContainer<ParamsContained, ParamsContainerOther...> Tp;
    //! This functor type
    typedef Functor<const Funct<Output(Arg...)>,const ParamsContainer<ParamsContained, ParamsContainerOther...>,Params...> ThisType;

    /* ====================  LIFECYCLE     ======================================= */
    //! Empty constructor
    Functor() = default;

    //! Construct from a functor and a parameter initializer list
    /** Note that the parameter type Tp must be brace-initializible 
     */
    Functor (Tf _funct, initializer_list<ParamsContained> lparam) : Functor(forward<Tf>(_funct), Tp{lparam}) {};

    //! Construct from a functor and a list of parameters
    Functor (Tf _funct, Params... params) : Functor(forward<Tf>(_funct), Tp{forward<Params>(params)...}) {};

    //! Copy construction
    Functor (const ThisType& src) : Functor(src.funct, src.params) {};

    //! Construct from a functor and a container of parameters -- Target constructor
    Functor (Tf _funct, Tp _params) : 
      funct(_funct), params(_params) {};
    /* ====================  ACCESSORS     ======================================= */
    //! Get the parameters
    Tp& getParams() const { return params; };

    /* ====================  OPERATORS     ======================================= */
    //! Evaluate the functor
    /** Note that we really need perfect forwarding here because 'arg' is 
     * lvalue in this function body but Tf::operator() could take rvalue
     */
    Output operator()(Arg... arg) const { return funct(forward<Arg...>(arg...)); };

    //! Compare two functors
    /** Note that this only compare the parameter sets against each other, which is 
     * in fact not enought to tell whether the 2 functors are different or not
     */
    bool operator==(const ThisType& rhs) const { return this == &rhs || this->getParams() == rhs.getParams(); }
    bool operator!=(const ThisType& rhs) const { return !((*this)==rhs); }
  protected:
    /* ====================  METHODS       ======================================= */

    /* ====================  DATA MEMBERS  ======================================= */
    Tf funct;
    Tp params;
}; /* ----------  end of template class Functor  ---------- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Make_Functor 
 *  Description:  bind any arbituray function with any arbituary parameters into a 
 *                std::function object. 
 * =====================================================================================
 */

template < size_t i > struct myplaceholder{};
namespace std {
  template < size_t i >
  struct is_placeholder< ::myplaceholder<i> > : public integral_constant<size_t, i>{};
}

//! This is the actual binder function that takes an list of index and parameters
//and bind the parameters to the function by std::bind. The number of parameters 
//is arbituray but must match the size of the index -- the trick is call std::bind
//using parameter pack expansion of the index list
template < size_t ... indices, template <size_t...> class I, class F, class ... Params >
auto Make_Functor_Impl(I<indices...>, const F& f, Params&&... params) 
  -> decltype ( bind(f, forward<Params>(params)..., myplaceholder<indices+1>{}...) )
{
  return bind(f, forward<Params>(params)..., myplaceholder<indices+1>{}...);
}

//! This function determines the number of parameters and make an 
//index list mapping each parameter to a integral constant so that 
//Make_Functor_Impl() can unroll the index list when calling std::bind
template < class Output, class ... Input, class ... Params >
auto Make_Functor ( Output (*const f)(Input...), Params&&... params ) 
  -> decltype ( Make_Functor_Impl(container_index<sizeof...(Input) - sizeof...(Params)>{}, f, forward<Params>(params)...) )
{
  constexpr size_t Nph = sizeof...(Input) - sizeof...(Params);
  using Index = container_index<Nph>;
  return Make_Functor_Impl(Index{}, f, forward<Params>(params)...);
}		/* -----  end of template function Make_Functor  ----- */

//! First argument is the return type; second argument is the pointer to a function
//The rest are the parameters you want to bind into the resultant functor --
//the number of parameters should be no larger than the number of arguments actually
//taken by the function pointer
template < class Functor, class Output, class ... Input, class ... Params >
auto Make_Functor_Wrapper (Functor dummy, Output (*const f)(Input...), Params&&... params) 
  -> decltype (dummy)
{
  return { Make_Functor(f, forward<Params>(params)...), {forward<Params>(params)...} };
}

//! This is for binding member function of type T. It's the same 
//as Make_Functor<Output, Input..., Params...> except that we have an 
//extra argument for the class type that contains the function
template < class T, class Output, class ... Input, class ... Params >
auto Make_Functor ( Output (T::*const f)(Input...), T* pobj, Params&&... params ) 
  -> decltype ( Make_Functor_Impl(container_index<sizeof...(Input) - sizeof...(Params)>{}, f, pobj, forward<Params>(params)...) )
{
  constexpr size_t Nph = sizeof...(Input) - sizeof...(Params);
  using Index = container_index<Nph>;
  return Make_Functor_Impl(Index{}, f, pobj, forward<Params>(params)...);
}		/* -----  end of template function Make_Functor  ----- */

//! This is for binding member function of type T. It's the same 
//as Make_Functor<Output, Input..., Params...> except that we have an 
//extra argument for the class type that contains the function
template < class Functor, class T, class Output, class ... Input, class ... Params >
auto Make_Functor_Wrapper (Functor dummy, T* pobj, Output (T::*const f)(Input...), Params&&... params) 
  -> decltype (dummy)
{
  return { Make_Functor(f, pobj, forward<Params>(params)...), {forward<Params>(params)...} };
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Quadratic
 *  Description:  1-dimensional quadratic function
 * =====================================================================================
 */
template < class T >
T Quadratic ( const T k, const T r0, const T c, const T& r )
{
  const T dr = r - r0;
  return k*dr*dr + c;
}		/* -----  end of function Quadratic  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  QuadraticFlat
 *  Description:  1-dimensional Flat-bottom quadratic function
 * =====================================================================================
 */
template < class T >
T QuadraticFlat (const T ref, 
               const T r0, 
               const T r1, 
               const T k0, 
               const T k1, 
               const T c,
	       const T& r) 
{
  const T dev = r - ref;
  if(dev >= r0 && dev <= r1) { return c; }
  else if(dev < r0) { 
	const T ldev = dev - r0;
	return k0*ldev*ldev + c;
  } else {
	const T rdev = dev - r1;
	return k1*rdev*rdev + c;
  }
}		/* -----  end of function QuadraticFlat  ----- */


template <class T>
T Linear ( const T a, const T b, const T& r )
{
  return a*r + b;
}		/* -----  end of template function Linear  ----- */

#endif
