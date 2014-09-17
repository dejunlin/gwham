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
    //! file name suffix
    const string fnsuffix;
    //! Mask indicates which columns are used; 
    const ulong Mask = 0;
    //! Total number of columns in input (should be <= the number of digits in Mask)
    const ulong iNcol = 0;
    //! Total number of columns where are interested in
    const ulong oNcol = 0;
    //! How many units of time per line of this time-series correspond to
    const ulong nst = 0;
  protected:
    /* ====================  METHODS       ======================================= */

    /* ====================  DATA MEMBERS  ======================================= */
    vector<uint> colids;

  public:
    typedef decltype(TimeSeries<INPUT>::Mask) MaskT;
    constexpr static int Nbit = numeric_limits<MaskT>::digits;
    /* ====================  LIFECYCLE     ======================================= */
    TimeSeries (fileio&& _fio, const string& _fnsuffix, const ulong _Mask, const ulong _iNcol, const ulong _oNcol, const ulong _nst) :
      fio(std::move(_fio)),
      fnsuffix(_fnsuffix),
      Mask(_Mask),
      iNcol(_iNcol),
      oNcol(_oNcol),
      nst(_nst)
      { init(); };

    TimeSeries (const fileio& _fio, const string& _fnsuffix, const ulong _Mask, const ulong _iNcol, const ulong _oNcol, const ulong _nst) : 
      fio(_fio), 
      fnsuffix(_fnsuffix),
      Mask(_Mask),
      iNcol(_iNcol),
      oNcol(_oNcol),
      nst(_nst)
      { init(); };
    
    TimeSeries (const ThisType& src) :
      TimeSeries(src.fio, src.fnsuffix, src.Mask, src.iNcol, src.oNcol, src.nst) {}; 
    /* ====================  ACCESSORS     ======================================= */

    /* ====================  MUTATORS      ======================================= */

    /* ====================  OPERATORS     ======================================= */
    //! Read a file and process them using the DATAPROCESSOR object
    /** This is the same as the TimeSeries::operator() except that 
     * the number of input columns will be checked for every line,
     * which is slow
     */
    template < class DATAPROCESSOR >
    linecounter read(const string& fname, const DATAPROCESSOR& dproc) {
      linecounter Nl = 0;
      if(!fio.fopen(fname)) { return Nl; }
      vector<INPUT> out(colids.size(),INPUT{0});
      while(fio.readaline()) {
	const vector<double> cols(fio.line2val());
	if(cols.size() != iNcol) { 
	  throw(TimeSeries_Exception("Wrong number of columns in time-series file: "+fname+" at line: '"+fio.line+"'. Expected "+tostr(iNcol)+" but got "+tostr(cols.size())));
	}
	for(uint i = 0; i < colids.size(); ++i) { out[i] = cols[colids[i]]; }
	dproc(out);
	++Nl;
      }
      return Nl;
    }

    //! Read a file and process them using the DATAPROCESSOR object
    template < class DATAPROCESSOR >
    linecounter operator()(const string& fname, const DATAPROCESSOR& dproc) {
      linecounter Nl = 0;
      if(!fio.fopen(fname)) { return Nl; }
      vector<INPUT> out(colids.size(),INPUT{0});
      if(fio.lb == 0 && fio.ls == 1 && fio.le >= MAXNLINE) {
        fio.skipemptylns();
        do {
          const vector<double> cols(fio.line2val());
          for(uint i = 0; i < colids.size(); ++i) { out[i] = cols[colids[i]]; }
          dproc(out);
          ++Nl;
        } while(fio.readeveryline());
      } else {
        while(fio.readaline()) {
          const vector<double> cols(fio.line2val());
          for(uint i = 0; i < colids.size(); ++i) { out[i] = cols[colids[i]]; }
          dproc(out);
          ++Nl;
        }
      }
      return Nl;
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
	  colids.emplace_back(i);
	}
      }
    }

    /* ====================  DATA MEMBERS  ======================================= */

}; /* ----------  end of template class TimeSeries  ---------- */

template < class INPUT>
constexpr int TimeSeries<INPUT>::Nbit;
#endif
