#if !defined(TIMESERIES_HPP)
#define TIMESERIES_HPP

/*
 * =====================================================================================
 *
 *       Filename:  timeseries.hpp
 *
 *    Description: Interfacial time-series-file reader class
 *
 *        Version:  1.0
 *        Created:  08/11/2014 08:14:27 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dejun Lin (), dejun.lin@gmail.com
 *   Organization:  Department of Biochemistry and Biophysics, Medical Center, University of Rochester
 *
 * =====================================================================================
 */
#include "typedefs.hpp"
#include "fileio.hpp"
#include "exception.hpp"
#include <vector>
#include "column_processor.hpp"

/*
 * =====================================================================================
 *        Class:  TimeSeries
 *  Description:  Interfacial time-series-file reader class
 * =====================================================================================
 */
template < class Output >
class TimeSeries
{
  public:
    /* ====================  LIFECYCLE     ======================================= */
    TimeSeries (fileio& _fio, const vector<column_processor>& _cps) : fio(_fio), cps(_cps) {};                             /* constructor */

    /* ====================  ACCESSORS     ======================================= */

    /* ====================  MUTATORS      ======================================= */

    /* ====================  OPERATORS     ======================================= */
    virtual linecounter operator()(const string& fname, Output& data) = 0;

  protected:
    /* ====================  METHODS       ======================================= */

    /* ====================  DATA MEMBERS  ======================================= */
    fileio fio;
    const vector<column_processor> cps;

  private:
    /* ====================  METHODS       ======================================= */

    /* ====================  DATA MEMBERS  ======================================= */

}; /* ----------  end of template class TimeSeries  ---------- */

#endif
