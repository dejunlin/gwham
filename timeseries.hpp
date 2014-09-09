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
#include "fileio_utils.hpp"
#include "exception.hpp"
#include <vector>
#include <functional>
#include <limits>
/*
 * =====================================================================================
 *        Class:  TimeSeries
 *  Description:  Interfacial time-series-file reader class
 * =====================================================================================
 */
template < class INPUT >
class TimeSeries
{
  public:
    typedef TimeSeries<INPUT> ThisType;
    typedef std::reference_wrapper<INPUT> rINPUT;
    //! the file reader
    fileio fio;
    //! Mask indicates which columns are used; 
    const ulong Mask = 0;
    //! Total number of columns in input (should be <= the number of digits in Mask)
    const ulong iNcol = 0;
    //! Total number of columns where are interested in
    const ulong oNcol = 0;
  protected:
    /* ====================  METHODS       ======================================= */

    /* ====================  DATA MEMBERS  ======================================= */
    vector<INPUT> vInput;
    vector<rINPUT> vrInput;

  public:
    typedef decltype(TimeSeries<INPUT>::Mask) MaskT;
    constexpr static int Nbit = numeric_limits<MaskT>::digits;
    /* ====================  LIFECYCLE     ======================================= */
    TimeSeries (fileio&& _fio, const ulong _Mask, const ulong _iNcol, const ulong _oNcol) :
      fio(std::move(_fio)), 
      Mask(_Mask),
      iNcol(_iNcol),
      oNcol(_oNcol),
      vInput(vector<INPUT>(iNcol))
      { init(); };

    TimeSeries (const fileio& _fio, const ulong _Mask, const ulong _iNcol, const ulong _oNcol) : 
      fio(_fio), 
      Mask(_Mask),
      iNcol(_iNcol),
      oNcol(_oNcol),
      vInput(vector<INPUT>(iNcol))
      { init(); };
    
    TimeSeries (const ThisType& src) :
      TimeSeries(src.fio, src.Mask, src.iNcol, src.oNcol) {}; 
    /* ====================  ACCESSORS     ======================================= */

    /* ====================  MUTATORS      ======================================= */

    /* ====================  OPERATORS     ======================================= */
    //! Read a file and process them using the DATAPROCESSOR object
    template < class DATAPROCESSOR >
    linecounter operator()(const string& fname, const DATAPROCESSOR& dproc) {
      linecounter Nl = 0;
      if(!fio.fopen(fname)) { return Nl; }
      vector<INPUT> out(vrInput.size());
      while(fio.readaline()) {
	const vector<double> cols(fio.line2val());
	if(cols.size() != iNcol) { 
	  throw(TimeSeries_Exception("Wrong number of columns in time-series file: "+fname+" at line: '"+fio.line+"'. Expected "+tostr(iNcol)+" but got "+tostr(cols.size())));
	}
	for(ulong i = 0; i < vInput.size(); ++i) vInput[i] = cols[i];
	for(ulong i = 0; i < out.size(); ++i) out[i] = vrInput[i].get();
	dproc(out);
	++Nl;
      }
    }

  private:
    /* ====================  METHODS       ======================================= */
    void init() {
      if(iNcol > Nbit) {
	  throw(TimeSeries_Exception("Number of column in time-series file (" + tostr(iNcol) +") exceeds the limit (" + tostr(Nbit) + ")"));
      }

      //This initialize the references to the columns in the input data
      for(ulong i = 0; i < iNcol; ++i) {
	if(Mask & (1ul<<i)) {
	  vrInput.emplace_back(vInput[i]);
	}
      }
    }

    /* ====================  DATA MEMBERS  ======================================= */

}; /* ----------  end of template class TimeSeries  ---------- */

template < class INPUT>
constexpr int TimeSeries<INPUT>::Nbit;
#endif
