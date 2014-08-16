/*
 * =====================================================================================
 *
 *       Filename:  metaprog_snippet.hpp
 *
 *    Description:  template metaprogramming snippets
 *
 *        Version:  1.0
 *        Created:  08/16/2014 12:14:35 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dejun Lin (), dejun.lin@gmail.com
 *   Organization:  Department of Biochemistry and Biophysics, Medical Center, University of Rochester
 *
 * =====================================================================================
 */

/*
 * =====================================================================================
 *        Class:  condition dependent type 
 *  Description:  The nested type depends on the condition and can only be either 
 *                TrueType or FalseType
 * =====================================================================================
 */
template < bool Condition, class TrueType, class FalseType > struct if_;

template < class TrueType, class FalseType > struct if_<true, TrueType, FalseType> {
  typedef TrueType result;
};

template < class TrueType, class FalseType > struct if_<false, TrueType, FalseType> {
  typedef FalseType result;
};

/*
 * =====================================================================================
 *        Class:  get the last element from parameter pack
 *  Description:  
 * =====================================================================================
 */
template < size_t ... i > struct Last_Index{};

template < size_t a, size_t ... rest > struct Last_Index<a, rest...> {
  constexpr static size_t value = Last_Index<rest...>::value;
};

template < size_t z > struct Last_Index<z> {
  constexpr static size_t value = z;
};

/*
 * =====================================================================================
 *        Class:  Extract the specified subset of elements from a parameter pack Superset
 *  Description:  The subset is defined by yet another sequence in Subset
 * =====================================================================================
 */
//! Extract the N'th element
template < class Superset, size_t N > struct Extract_One_Element{};
template < template <size_t...> class I, size_t n, size_t ... indices, size_t N > 
struct Extract_One_Element<I<n, indices...>, N> {
  constexpr static size_t value = Extract_One_Element<I<indices...>, N-1>::value;
};
template < template <size_t...> class I, size_t n, size_t ... indices> 
struct Extract_One_Element<I<n, indices...>, 1> {
  constexpr static size_t value = n;
};

//! Extract a sequence of element
template < class Superset, class SubsetIDs, class Output  > struct Extract_Elements{};
template < size_t ... indices, size_t id1, size_t ... subsetids, size_t ... output, template <size_t...> class I >
struct Extract_Elements< I<indices...>, I<id1, subsetids...>, I<output...> >  {
  typedef I<indices...> Superset;
  constexpr static size_t newelement = Extract_One_Element<Superset, id1>::value;
  typedef typename Extract_Elements<Superset, I<subsetids...>, I<output..., newelement> >::type type;  
};
template < size_t ... indices, size_t ... output, template <size_t...> class I >
struct Extract_Elements< I<indices...>, I<>, I<output...> >  {
  typedef I<output...> type;  
};


/*
 * =====================================================================================
 *        Class:  Indices generator for constexpr integer 
 *  Description:  These helper classes define a compile-time sequence 
 *                Indices<N1, N2, ..., > for a given N1, a propagator template P<N> 
 *                where Ni+1 == P<Ni>::value and a termination condition T<Indices<N1...> >
 *                that determine the end of the sequence.
 *                The idea is to recursively define a nested type in the helper template
 *                class make_indices where a new Ni is pushed back in Indices<N1..., Ni-1> 
 *                at each recursion step until the termination condistion is satisfied
 * =====================================================================================
 */
//! These are the propagators
//! Propagate the sequence indices...,Ni by appending a new element Ni+dN
template < size_t dN > struct IncrdN {
  /** This nested template is needed in order to hide the parameter dN
   * from the user. E.g., if another template class takes a templated class
   * with one type parameter as parameter, we should plug in the alias:
   * template < class T > using Incr = typename IncrdN<N>::template Incr<T>;
   */
  template < class T > struct Incr{};
  template < size_t ... indices, template <size_t...> class I>
  struct Incr<I<indices...> > {
    constexpr static size_t Ni = Last_Index<indices...>::value;
    typedef I<indices..., Ni+dN> type;
  };
};
//! Propagate the Fibonacci sequence
template < class T > struct Fibonacci{};
template < size_t ... indices, template <size_t...> class I > 
struct Fibonacci<I<indices...> > {
  typedef I<indices...> Indices;
  constexpr static size_t j = sizeof...(indices);
  constexpr static size_t i = j-1;
  constexpr static size_t k = Extract_One_Element<Indices, i>::value + Extract_One_Element<Indices, j>::value;
  typedef I<indices..., k> type;
};

template <template <size_t...> class I >
struct Fibonacci<I<> > {
  typedef typename Fibonacci<I<1, 1> >::type type;
};

template <template <size_t...> class I >
struct Fibonacci<I<1> > {
  typedef typename Fibonacci<I<1, 1> >::type type;
};

//! These are the termination condition
//! Terminate if a maximum number of elements is reached in the sequence
template < size_t Max > struct StopIfMax {
  template < class T > struct StopIfMaxN {};
  template < size_t ... indices, template <size_t...> class I> 
  struct StopIfMaxN<I<indices...> > {
    constexpr static bool decision = sizeof...(indices) >= Max ? true : false;
  };
};

//! This is the target type to be generated
template < size_t ... N > struct Index {};

//! Primary template
template < bool Terminate, class Indices, class Propagator, class Terminator > 
struct make_indices{};
//! Propagate if Terminate == false
template < class Indices, 
	   template <class> class Propagator, 
	   template <class> class Terminator 
	 >
struct make_indices< false, Indices, Propagator<Indices>, Terminator<Indices> > {
  typedef typename Propagator<Indices>::type newIndices;
  typedef Propagator<newIndices> newPropagator;
  typedef Terminator<newIndices> newTerminator;
  constexpr static bool newdecision = newTerminator::decision;
  typedef typename make_indices<newdecision, newIndices, newPropagator, newTerminator>::type type; 
};

//! Stop if Terminate == true 
template < class Indices, class Propagator, class Terminator > 
struct make_indices<true, Indices, Propagator, Terminator> {
  typedef Indices type;
};
