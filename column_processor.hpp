#if !defined(COLUMNPROCESSOR_HPP)
#define COLUMNPROCESSOR_HPP

/*
 * =====================================================================================
 *
 *       Filename:  column_processor.hpp
 *
 *    Description:  Generic columnated data processor
 *
 *        Version:  1.0
 *        Created:  12/08/14 10:12:26
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dejun Lin (DL), dejun.lin@gmail.com
 *   Organization:  Department of Biochemistry and Biophysics, Medical Center, University of Rochester
 *
 * =====================================================================================
 */
#include "typedefs.hpp"
#include <vector>
#include <iterator>

using namespace std;
/*
 * =====================================================================================
 *        Class:  ColumnProcessor
 *  Description:  Generic columnated data processor
 * =====================================================================================
 */
class ColumnProcessor
{
  public:
    /* ====================  LIFECYCLE     ======================================= */
    ColumnProcessor (const uint& _Ncols, const uint& _Mask) : Ncols(_Ncols), Mask(_Mask) {};                             /* constructor */

    /* ====================  ACCESSORS     ======================================= */

    /* ====================  MUTATORS      ======================================= */

    /* ====================  OPERATORS     ======================================= */
    //! read from the iterator, point it to the right location for the next processor
    void operator()(vector<valtype>::const_iterator& it, vector<valtype>& output) const {
      for(uint i = 0; i < Ncols; ++i) {
	if(Masks & (1<<i)) { ouptput.push_back(*it); }
	++it;
      }
    };

  protected:
    /* ====================  METHODS       ======================================= */

    /* ====================  DATA MEMBERS  ======================================= */
    //! Total number of columns expected
    const uint Ncols;
    //! Which columns to be used
    const uint Masks;

  private:
    /* ====================  METHODS       ======================================= */

    /* ====================  DATA MEMBERS  ======================================= */

}; /* -----  end of class ColumnProcessor  ----- */

#endif
