/*
  Nelder-Meade Simplex Optimizer
  (c) 2009 Tod D. Romo, Grossfield Lab, URMC
  
  Based loosely on the NRC implementation (Numerical Recipes in C (1996):411)
*/

#if !defined(LOOS_SIMPLEX_HPP)
#define LOOS_SIMPLEX_HPP

#include <iostream>
#include <vector>
#include "../typedefs.hpp"

//! Nelder-Meade Simplex Optimizer (based loosely on the NRC (1996) implementation)
/**
 * The Simplex class is templated so it can work with doubles, floats,
 * etc.  What is actually optimized is a functor that takes a vector
 * containing the parameters.
 *
 * For example, suppose you're trying to optimize the 2D Rosenbrock
 * function, you could define a functor like,
\code
struct Rosenbrock {
  valtype operator()(const vector<valtype>& v) {
    valtype x1 = 1.0 - v[0];
    valtype x2 = v[1] - v[0] * v[0];
    valtype x = x1 * x1 + 105.0 * x2 * x2;
    return(x);
  }
};
\endcode
 * 
 * Now to optimize this, create a Simplex instance with the same type,
 * and appropriate number of dimensions,
\code
Simplex<valtype> opt(2);
\endcode
 * 
 * This create a new 2-D Simplex optimizer.  The optimizer needs to
 * have a starting guess (seed) as well as an initial step size
 * (length).  For the Rosenbrock, let's make the initial guess (0,0)
 * and length of 2 in both dimensions,
\code
vector<valtype> seeds(2, 0.0);
vector<valtype> lengths(2, 2.0);

opt.seedLengths(lengths);
vector<valtype> final = opt.optimize(seeds, Rosenbrock());
valtype final_value = opt.finalValue();
\endcode
 *
 * Now \a final is a vector containing the best fit parameters and \a
 * final_value is the value of the Rosenbrock function at that
 * location.
 */
template<typename T = valtype>
class Simplex {
  typedef std::vector< std::vector<T> >     VVectors;


  template<class C>
  T modify(T factor, C& ftor) {
    int j;
    T f1, f2, val;
    
    f1 = (1.0 - factor) / ndim;
    f2 = f1 - factor;
    
    // Compute the new simplex vertex
    
    for (j=0; j<ndim; ++j)
      trial[j] = simpsum[j] * f1 - simplex[worst][j] * f2;
    
    val = ftor(trial);
    
    // Use this vertex in the simplex if it's better...
    
    if (val < values[worst]) {
      values[worst] = val;
      for (j=0; j<ndim; j++) {
        simpsum[j] += trial[j] - simplex[worst][j];
        simplex[worst][j] = trial[j];
      }
    }
    
    return(val);
  }



  // The core of the simplex optimizer...
  template<class C>
  void core(C& ftor) {
    int i, next_worst, j, mpts;
    T sum, saved, val, num, den;

    mpts = ndim + 1;

    for (j=0; j<ndim; j++) {
      for (sum = 0.0, i = 0; i<mpts; i++)
        sum += simplex[i][j];
      simpsum[j] = sum;
    }

    n_evals = 0;

    while (n_evals <= maxiters) {
      best = 1;
      if (values[0] > values[1]) {
        worst = 0;
        next_worst = 1;
      } else {
        worst = 1;
        next_worst = 0;
      }


      // Find candidate vertices in the simplex...

      for (i=0; i<mpts; i++) {
        if (values[i] <= values[best])
          best = i;
        if (values[i] > values[worst]) {
          next_worst = worst;
          worst = i;
        } else if (values[i] > values[next_worst] && i != worst)
          next_worst = i;
      }

      // Check for convergence...
      // Some conditions may have the numerator and denominator equal (or both 0)
      // which can cause problems below, hence the special test.

      num = fabs(values[worst] - values[best]);
      den = fabs(values[worst]) + fabs(values[best]);
      rtol = 2.0 * num / den;
      if (rtol < tol || num == den)
        return;

      // Now try reflecting, contracting, expanding the simplex...

      n_evals += 2;
      val = modify(-1.0, ftor);
      if (val <= values[best])
        val = modify(2.0, ftor);
      else if (val >= values[next_worst]) {

        saved = values[worst];
        val = modify(0.5, ftor);

        if (val >= saved) {
          for (i=0; i<mpts; i++) {
            if (i != best) {
              for (j=0; j<ndim; j++)
                simplex[i][j] = simpsum[j] = 0.5 * (simplex[i][j] + simplex[best][j]);
              values[i] = ftor(simpsum);
            }
          }
          n_evals += ndim;

          for (j=0; j<ndim; j++) {
            for (sum=0.0, i=0; i<mpts; i++)
              sum += simplex[i][j];
            simpsum[j] = sum;
          }
        }

      } else
        --n_evals;

    }

  }


  void allocateSpace(const int n) {
    q = std::vector<T>(n+1, 0.0);
    qq = std::vector<T>(n+1, 0.0);
    simpsum = std::vector<T>(n, 0.0);
    values = std::vector<T>(n+1, 0.0);
    trial = std::vector<T>(n, 0.0);
    simplex.clear();
    for (int i=0; i<n+1; ++i)
      simplex.push_back(std::vector<T>(n, 0.0));
  }



public:

  Simplex(const int n) : tol(1e-3), ndim(n), maxiters(2000), best(-1), worst(-1), n_evals(0) {
    allocateSpace(n);
  }

  //! Set the number of dimensions
  void dim(const int n) { ndim = n; allocateSpace(n); }

  //! Initial guess
  void seedLengths(const std::vector<T>& seeds) { characteristics = seeds; }

  //! Convergence criterion
  void tolerance(const valtype d) { tol = d; }

  //! Limit on the number of function evaluations to perform
  void maximumIterations(const int n) { maxiters = n; }

  int maximumIterations() const { return maxiters; }
 
  //! Return the number of iterations that it actually went through 
  int numberOfIterations() const { return n_evals; }

  //! Retrieve the final (best fit) parameters
  std::vector<T> finalParameters(void) const {
    if (best < 0)
      throw(std::logic_error("Simplex has not been optimized"));
    return(simplex[best]);
  }

  //! Final (best) value
  T finalValue(void) const { return(values[best]); }
 
  //! Optimize via a functor
  template<class C>
  std::vector<T> optimize(std::vector<T>& f, C& ftor) {
    int i, j, n;

    if (characteristics.size() != f.size())
      throw(std::logic_error("Invalid seed"));

    n = ndim + 1;

    // Initial simplex is based on Nelder-Meade's construction...

    for (i=0; i<ndim; i++) {
      q[i] = characteristics[i] * ( (sqrtf(n) + ndim) / (n * sqrt(2.0)) );
      qq[i] = characteristics[i] * ( (sqrtf(n) - 1.0) / (n * sqrt(2.0)) );
    }

    for (j=0; j<n; j++)
      for (i=0; i<ndim; i++)
        if (j == i+1)
          simplex[j][i] = f[i] + qq[i];
        else
          simplex[j][i] = f[i] + q[i];

    for (j=0; j<n; ++j) {
      values[j] = ftor(simplex[j]);
    }

    core(ftor);
    return(simplex[best]);
  }

private:
  valtype tol;
  int ndim, maxiters, best, worst;
  valtype rtol;
  int n_evals;

  std::vector<T> characteristics;
  std::vector<T> simpsum;
  std::vector<T> values;
  std::vector<T> q;
  std::vector<T> qq;
  std::vector<T> trial;
  VVectors simplex;

};




#endif
