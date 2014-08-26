#if !defined(METAPROG_SNIPPETS_HPP)
#define METAPROG_SNIPPETS_HPP
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
#include <iostream>
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
 *                Indices<N1, N2, ..., > for a given N1, a propagator template P<Indices<N1,...> > 
 *                where Ni+1 == P<Indices<N1,N2,...,Ni> >::value and a termination 
 *                condition T<Indices<N1...> > that determine the end of the sequence.
 *                The idea is to recursively define a nested type in the helper template
 *                class make_index_seq where a new Ni is pushed back in Indices<N1..., Ni-1> 
 *                at each recursion step until the termination condistion is satisfied
 * =====================================================================================
 */
//! These are the propagators
//! Propagate the arithmetic sequence with a commone difference of dN and initial term Ni
template < class T, class Arg > struct ArithmeticSeq{};
template < size_t ... indices, template <size_t...> class I, size_t Ni, size_t dN>
struct ArithmeticSeq<I<indices...>, I<Ni, dN> > {
  typedef I<indices..., Ni+dN> type;
  typedef I<Ni+dN, dN> seed;
};
//! Specialize for index sequence
template < template <size_t...> class I>
struct ArithmeticSeq<I<>, I<> > {
  typedef I<0> type;
  typedef I<1, 1> seed;
};

//! Propagate the Fibonacci sequence
template < class T, class Arg > struct Fibonacci{};
template < size_t ... indices, template <size_t...> class I, size_t Ni, size_t Nj > 
struct Fibonacci<I<indices...>, I<Ni,Nj> > {
  typedef I<indices...> Indices;
  typedef I<indices..., Ni+Nj> type;
  typedef I<Nj, Ni+Nj> seed; 
};

template <template <size_t...> class I >
struct Fibonacci<I<>, I<> > {
  typedef typename Fibonacci<I<1, 1>, I<1,1> >::type type;
  typedef I<1, 2> seed;
};

template <template <size_t...> class I >
struct Fibonacci<I<1>, I<> > {
  typedef typename Fibonacci<I<1, 1>, I<1,1> >::type type;
  typedef I<1, 2> seed;
};

//! These are the termination condition
//! Terminate if a maximum number of elements is reached in the sequence
template < class T, class Arg > struct StopAtMaxN {};
template < size_t ... indices, template <size_t...> class I, size_t Max> 
struct StopAtMaxN<I<indices...>, I<Max> > {
  constexpr static bool decision = sizeof...(indices) >= Max ? true : false;
};

//! This is the target type to be generated
template < size_t ... N > struct IndexSeq {};

//! Primary template
template < bool Terminate, class Indices, class Propagator, class Terminator > 
struct make_index_seq{};
//! Propagate if Terminate == false
template < class Indices, 
	   template <class, class> class Propagator,
	   class PropagatorArg, 
	   template <class, class> class Terminator, 
	   class TerminatorArg 
	 >
struct make_index_seq< false, Indices, Propagator<Indices, PropagatorArg>, Terminator<Indices, TerminatorArg> > {
  typedef typename Propagator<Indices, PropagatorArg>::type newIndices;
  typedef typename Propagator<Indices, PropagatorArg>::seed newPropagatorArg;
  typedef Propagator<newIndices,newPropagatorArg> newPropagator;
  typedef Terminator<newIndices,TerminatorArg> newTerminator;
  constexpr static bool newdecision = newTerminator::decision;
  typedef typename make_index_seq<newdecision, newIndices, newPropagator, newTerminator>::type type; 
};

//! Stop if Terminate == true 
template < class Indices, class Propagator, class Terminator > 
struct make_index_seq<true, Indices, Propagator, Terminator> {
  typedef Indices type;
};

template < size_t N >
using make_container_index = 
  make_index_seq< StopAtMaxN<IndexSeq<>, IndexSeq<N>>::decision,
                  IndexSeq<>,
		  ArithmeticSeq<IndexSeq<>, IndexSeq<>>,
		  StopAtMaxN<IndexSeq<>, IndexSeq<N>> >;

template < size_t N >
using container_index = typename make_container_index<N>::type;

/*
 * =====================================================================================
 *        Class:  Container Copier 
 *  Description:  These helper classes copy the subset of a container of type Src into 
 *                another container of Type Des where the size of Src is known at 
 *                compilation time. The subset to be copied is identified by a list of 
 *                indices. 
 * =====================================================================================
 */
//! Primary template
template < class Src, class Des, class ... Sub > struct ContainerCopier{};
//! Copy a subset
template < class Src, class Des, template <size_t...> class I, size_t ... indices > 
struct ContainerCopier<Src, Des, I<indices...> > {
  ContainerCopier(Src src) : value({src[indices]...,}) {
  };
  Des value;
};


/*
 * =====================================================================================
 *        Class:  check_if
 *  Description:  check if a class T satisify a condition 
 *                The idea is to make the compiler decide which member function template
 *                'f' to overload based on whether the parameter 'Condition' has a 
 *                valid definition of nested template class 'type' -- if yes, then 
 *                TMached f(Condition::template type<C>*) is seen because interpreting
 *                '0' in f<T>(0) as Condition::template type<C>* is more specific than
 *                the vararg list '...'.
 *                The implementation of class Condition for checking if the class T has 
 *                a member function 'Mfn' could look like this:
 *                template <class Ret, class ... Arg>
 *                struct Condition {
 *                  template < class T, Ret (T::*)(Arg&&...) = &T::Mfn > struct type {};
 *                };
 * =====================================================================================
 */
template < class T, class Condition >
struct check_if
{
  typedef char TMatched;
  typedef long TUnmatched;
  
  template < class C > static TMatched f(const typename Condition::template type<C>*) {};
  template < class C > static TUnmatched f(...) {};

  constexpr static bool value = (sizeof(f<T>(0)) == sizeof(TMatched));
}; /* ----------  end of template class check_if  ---------- */

#endif
