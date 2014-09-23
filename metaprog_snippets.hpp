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
#include <memory>
#include <type_traits>
#include <utility>
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
//! Extract the N'th index 
template < class Superset, size_t N > struct Extract_One_Index{};
template < template <size_t...> class I, size_t n, size_t ... indices, size_t N > 
struct Extract_One_Index<I<n, indices...>, N> {
  constexpr static size_t value = Extract_One_Index<I<indices...>, N-1>::value;
};
template < template <size_t...> class I, size_t n, size_t ... indices> 
struct Extract_One_Index<I<n, indices...>, 1> {
  constexpr static size_t value = n;
};

//! Extract a sequence of indices 
template < class Superset, class SubsetIDs, class Output  > struct Extract_Indices{};
template < size_t ... indices, size_t id1, size_t ... subsetids, size_t ... output, template <size_t...> class I >
struct Extract_Indices< I<indices...>, I<id1, subsetids...>, I<output...> >  {
  typedef I<indices...> Superset;
  constexpr static size_t newelement = Extract_One_Index<Superset, id1>::value;
  typedef typename Extract_Indices<Superset, I<subsetids...>, I<output..., newelement> >::type type;  
};
template < size_t ... indices, size_t ... output, template <size_t...> class I >
struct Extract_Indices< I<indices...>, I<>, I<output...> >  {
  typedef I<output...> type;  
};

//! Extract the N'th type
template < size_t N, class ... T > struct Extract_One_Type {
  static_assert(N <= sizeof...(T), "Extract_One_Type<T..., N> index out of bound");
  typedef void type;
};
template < size_t N, class T1, class ... T > 
struct Extract_One_Type<N, T1, T...> {
  typedef Extract_One_Type<N-1, T...> NextIteration;
  typedef typename NextIteration::type type;
  static inline const type& get(const T1&, const T&...args) {
    return NextIteration::get(args...); 
  }
};
template < class T1, class ... T>
struct Extract_One_Type<1, T1, T...> {
  typedef T1 type;
  static inline const type& get(const T1& a, const T&...) {
    return a;
  }
};

template <>
struct Extract_One_Type<1> {
  typedef void type;
  static inline type get() {
    return;
  }
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
 *                std::true_type f(int*) is seen because interpreting '0' in f<T>(0) 
 *                as int* is more specific than the vararg list '...'.
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
  template < class C > 
  static auto f(int*)
    ->
    decltype(
      //NOTE: the first expression (left operand in comma operator) is discarded
      std::declval<typename Condition::template type<C>>(), 
      //NOTE: std::true_type is used if the left operand is valid
      std::true_type() 
    ); 

  template < class C > static auto f(...) -> std::false_type;

  constexpr static bool value = decltype(f<T>(0))::value;
}; /* ----------  end of template class check_if  ---------- */

//Here are some examples of the Condition class
//! This class check if a class T has member function 
// T::size_type size() const 
struct has_memfn_size {
  template < class T, typename T::size_type (T::*)() const = &T::size >
  struct type {};
};

//! This class check if a class T has member function 
// void emplace_back(T::value_type&) 
struct has_memfn_emplace_back {
  template < class T, void (T::*)(typename T::value_type&) = &T::emplace_back >
  struct type {};
};

//! This class check if a class T has member function 
// T::type get() const 
struct has_memfn_get {
  template < class T, typename T::type& (T::*)() const = &T::get >
  struct type {};
};

//!This class check if a class T can be iterated using nested const_iterator 
struct can_be_const_iterate {
  template < class T, 
  class C = typename T::const_iterator,
  C (T::*)() const = &T::begin,
  C (T::*)() const = &T::end
           >
  struct type {};
};

//! This class check if a class T has unerlying elements
//that can be retrieved by std::get<0>
struct can_be_get {
  template < class T, typename T::value_type& (*)(T&) = &std::get<0> >
  struct type {};
};

template <class S>
struct can_be_streamed {
  template < class T >
  using type = std::remove_reference<decltype(std::declval<S&>() << std::declval<T>())>;
};

/*
 * =====================================================================================
 *        Class:  Type trait of std::string 
 *  Description:  
 * =====================================================================================
 */
template < class T > struct is_string : std::false_type {};
template <typename charT, typename traits, typename Alloc>
struct is_string<std::basic_string<charT, traits, Alloc> > : std::true_type {};

/*
 * =====================================================================================
 *        Class:  Type trait of the default floating point number, i.e., float, double or long double 
 *  Description:  
 * =====================================================================================
 */
template < class T > struct is_stdfloat : std::false_type {};
template <>
struct is_stdfloat<float> : std::true_type {};
template <>
struct is_stdfloat<double> : std::true_type {};
template <>
struct is_stdfloat<long double> : std::true_type {};

/*
 * =====================================================================================
 *        Class:  Type trait of shared_ptr 
 *  Description:  
 * =====================================================================================
 */
namespace std {
template < class T >
struct is_pointer<std::shared_ptr<T>> : std::true_type {};
}

/*
 * =====================================================================================
 *        Class:  FormatStream 
 *  Description:  stream any object with a certain format
 * =====================================================================================
 */
template < class S >
class FormatStream {
  using ThisType = FormatStream<S>;
  public:
    //! Constructor with default setting
    FormatStream(S& _stream) : FormatStream(_stream, 0, 6) {};

    //! Target constructor
    FormatStream(S& _stream, const uint& _width, const uint& _precision) :
      stream(_stream), w(_width), prec(_precision) {};
    
    //! specialize for string
    ThisType& operator<< (const std::string& input) {
      stream.width(w);
      stream.precision(prec);
      stream << input;
      return *this;
    }
    

    //! stream the basic scalar type to the contained object stream
    template < class T >
    typename std::enable_if<
      check_if<T, can_be_streamed<S>>::value
    ,ThisType&>::type
    operator<< (const T& input) {
      stream.width(w);
      stream.precision(prec);
      stream << input;
      return *this;
    }
    
    //! stream a type which can be const-iterated through and each element can be streamed
    template < class T >
    typename std::enable_if<
      check_if<T, can_be_const_iterate>::value &&
      check_if<decltype( *(std::declval<typename T::const_iterator>()) ), can_be_streamed<ThisType>>::value
    ,ThisType&>::type
    operator<< (const T& input) {
      for(auto& i : input) {
        *this << i; 
      }
      return *this;
    }
    
    //! stream a type that has member function get
    template < class T >
    typename std::enable_if<
      check_if<T, has_memfn_get>::value && 
      check_if<decltype( std::declval<T>().get() ), can_be_streamed<ThisType>>::value
    ,ThisType&>::type
    operator<< (const T& input) {
      *this << input.get();
      return *this;
    }
    
    //! just handle something like std::endl
    ThisType& operator<< (S& (*f)(S&)) {
      f(stream);
      return *this;
    }
    
    //! set the width
    void width(const uint _width) { w = _width; }

    //! set the precision
    void precision(const uint _precision) { prec = _precision; }

    //! set whatever flags the S object has
    template < class F >
    void flags(const F& flags) { stream.flags(flags); }

  private:
    S& stream;
    uint w;
    uint prec;
};

static FormatStream<decltype(std::cout)> fcout{std::cout}; 


/*
 * =====================================================================================
 *     Function: Greatest common divisor 
 *  Description: binary method
 * =====================================================================================
 */

//This function get the number of trailing zeros of an unsigned intgral 
template <class T>
typename std::enable_if<std::is_integral<T>::value && std::is_unsigned<T>::value, T>::type
ctz(T a) 
{
  T shift = 0;
  for(; a & 1 == 0; ++shift) a >> 1;
}

#ifdef __GNUC__
#define CTZ(x) __builtin_ctz(x) 
#else
#define CTZ(x) ctz(x)
#endif 

template <class T>
typename std::enable_if<std::is_integral<T>::value && std::is_unsigned<T>::value, T>::type
GCD ( T a, T b )
{
  const T zero = T(0);
  if(a == zero) return a;
  if(b == zero) return b;
  const T shift = CTZ(a|b);
  a >>= CTZ(a);
  do {
    b >>= CTZ(b);
    if(a>b) { std::swap(a,b); }
    b -= a;
  } while ( b != 0 );				/* -----  end do-while  ----- */
  return a << shift;
}		/* -----  end of template function GCD  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  LCM
 *  Description:  Largets common multiple
 * =====================================================================================
 */
template <class T>
typename std::enable_if<std::is_integral<T>::value && std::is_unsigned<T>::value, T>::type
LCM ( T a, T b )
{
  return a*b/GCD(a,b);
}		/* -----  end of function LCM  ----- */

#endif
