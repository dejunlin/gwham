/*
 * =====================================================================================
 *
 *       Filename:  Optimizer.hpp
 *
 *    Description:  Interface to various optimization routines   
 *
 *        Version:  1.0
 *        Created:  05/28/2014 09:27:50 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dejun Lin, 
 *   Organization:  
 *
 * =====================================================================================
 */
#include "dlib/optimization.h"
#include <algorithm>
#include "typedefs.hpp"
#include <climits>

namespace dlibOPT {
/*
 * =====================================================================================
 *        Class:  Funct
 *  Description:  Wrapper to kernel which calculates the value of the function
 *  NOTE that this should be used in conjunction with gradFunct since it requires
 *  Funct and gradFunct are called in pair such that Funct is called first and then 
 *  is gradFunct. The Kernel::operator()() is the function that calculates both 
 *  function value and gradient, which are stored as member variables of Kernel
 * =====================================================================================
 */
template < class Kernel,                        /* The kernel that actually do the calculation */
	   class iKT,                           /*  Input type for Kernel::operator() */
	   class iT,                            /*  Input type of argument to whatever optimziation routine*/
	   class oT                             /*  Output type of argument to whatever optimziation routine*/
	 >
class Funct
{
  public:
    /* ====================  LIFECYCLE     ======================================= */
    Funct (Kernel* const _pKernel) : pKernel(_pKernel) {};                             /* constructor */

    /* ====================  ACCESSORS     ======================================= */

    /* ====================  MUTATORS      ======================================= */

    /* ====================  OPERATORS     ======================================= */
    oT operator()(const iT& x) const {
      const iKT ikx(x.begin(), x.end());
      pKernel->operator()(ikx);
      return pKernel->getF();
    };

  protected:
    /* ====================  METHODS       ======================================= */

    /* ====================  DATA MEMBERS  ======================================= */

  private:
    /* ====================  METHODS       ======================================= */

    /* ====================  DATA MEMBERS  ======================================= */
    Kernel* const pKernel;

}; /* ----------  end of template class Funct  ---------- */


/*
 * =====================================================================================
 *        Class:  gradFunct
 *  Description:  Wrapper to kernel which calculates the gradient of the function
 *  NOTE that this should be used in conjunction with Funct since it requires
 *  Funct and gradFunct are called in pair such that Funct is called first and then 
 *  is gradFunct. The Kernel::operator()() is the function that calculates both 
 *  function value and gradient, which are stored as member variables of Kernel
 * =====================================================================================
 */
template < class Kernel,                        /* The kernel that actually do the calculation */
	   class iKT,                           /*  Input type for Kernel::operator() */
	   class iT,                            /*  Input type of argument to whatever optimziation routine*/
	   class oT                             /*  Output type of argument to whatever optimziation routine*/
	 >
class gradFunct
{
  public:
    /* ====================  LIFECYCLE     ======================================= */
    gradFunct (Kernel* const _pKernel) : pKernel(_pKernel) {};                             /* constructor */

    /* ====================  ACCESSORS     ======================================= */

    /* ====================  MUTATORS      ======================================= */

    /* ====================  OPERATORS     ======================================= */
    oT operator()(const iT& x) const {
      // assume the return type of Kernel::getgradF() is iKT
      const iKT& gradF = pKernel->getgradF();
      const int N = gradF.size();
      oT rgradF(N);
      for(int i = 0; i < N; ++i) { rgradF(i) = gradF[i]; } 
      return rgradF; 
    };

  protected:
    /* ====================  METHODS       ======================================= */

    /* ====================  DATA MEMBERS  ======================================= */

  private:
    /* ====================  METHODS       ======================================= */

    /* ====================  DATA MEMBERS  ======================================= */
    Kernel* const pKernel;

}; /* ----------  end of template class gradFunct  ---------- */


template<class Kernel, 
	 class iT
        >
void optimize_dlib_lbfgs(Kernel& funct, iT& startingpoint, const double& tol) {
  typedef dlib::matrix<double,0,1> cvect;
  Funct<Kernel, iT, cvect, double> ft(&funct);
  gradFunct<Kernel, iT, cvect, cvect> gradft(&funct);
  cvect sp(startingpoint.size());
  cvect upbound(startingpoint.size());
  cvect lowbound(startingpoint.size());
  for(uint i = 0; i < startingpoint.size(); ++i) {
    sp(i) = startingpoint[i];
    upbound(i) = 99999999;
    lowbound(i) = DBL_MIN;
  }
  /*  dlib::find_min(dlib::lbfgs_search_strategy(10),  // The 10 here is basically a measure of how much memory L-BFGS will use.
           dlib::objective_delta_stop_strategy(tol).be_verbose(),  // Adding be_verbose() causes a message to be 
                                                                    // printed for each iteration of optimization.
            ft, gradft, sp, -1);
  */
  dlib::find_min_box_constrained(dlib::lbfgs_search_strategy(10),  // The 10 here is basically a measure of how much memory L-BFGS will use.
           dlib::objective_delta_stop_strategy(tol).be_verbose(),  // Adding be_verbose() causes a message to be 
                                                                    // printed for each iteration of optimization.
            ft, gradft, sp, lowbound, upbound);
}

}

