#if !defined(FILEIO_HPP)
#define FILEIO_HPP

#include "typedefs.hpp"
#include <fstream>
#include <vector>
#include "fileio_utils.hpp"

using namespace std;

class fileio {
  public:
    fileio(const string& _fname, const ios_base::openmode& mode, const bool& perm=true, const linecounter _lb=0, const linecounter _ls=1, const linecounter _le=MAXNLINE, const string& cm="#@;");
    fileio(const ios_base::openmode& mode, const bool& perm=true, const linecounter _lb=0, const linecounter _ls=1, const linecounter _le=MAXNLINE, const string& cm="#@;");
    fileio(const fileio& _fio);

    //! if fileio::line contains only whitespaces
    inline
    #ifdef __GNUC__
    __attribute__((always_inline))
    #endif
    bool emptyline() const { return emptystr(line, comments); }

    //! skip empty lines 
    /** this will read the first non-empty line into fileio::line
     */
    inline
    #ifdef __GNUC__
    __attribute__((always_inline))
    #endif
    void skipemptylns() {
      while(getline(fs, line)) {
        if(emptyline()) { continue; }
	else { break; }
      }
    }

    //! Read a time-series style file without line number check
    /** Assume in the input file all the empty lines are 
     * in the head of the files and are already skipped
     * and no such empty lines should occur for the rest 
     * and that you need every line from the file, 
     * i.e., lb:ls:le = 0:1:MAXNLINE
     */
    inline
    #ifdef __GNUC__
    __attribute__((always_inline))
    #endif
    bool readtsnb() {
      getline(fs, line);
      ++lc;
      return !line.empty();
    }

    //! Read a time-series style file
    /** Assume in the input file all the empty lines are 
     * in the head of the files and are already skipped
     * and no such empty lines should occur for the rest
     */
    inline
    #ifdef __GNUC__
    __attribute__((always_inline))
    #endif
    bool readts() {
      while(getline(fs, line)) {
        ++lc;
        if(lc < lb) { line.clear(); continue; }
        else if(lc > le) { line.clear(); break; }
        else if(lc % ls) { line.clear(); continue; }
        break;
      }
      return !line.empty();
    }

    //! read one line at a time and cache it into fileio::line
    /** This assumes that you need every line from the file, 
     * i.e., lb:ls:le = 0:1:MAXNLINE
     */
    inline
    #ifdef __GNUC__
    __attribute__((always_inline))
    #endif
    bool readalinenb() {
      while(getline(fs, line)) {
        if(emptyline()) { continue; }
	else { break; }
      };
      ++lc;
      return !line.empty();
    }

    //! read one line at a time and cache it into fileio::line
    inline
    #ifdef __GNUC__
    __attribute__((always_inline))
    #endif
    bool readaline() {
      line.clear();
      while(getline(fs, line)) {
        if(emptyline()) { continue; }
        ++lc;
        if(lc < lb) { line.clear(); continue; }
        else if(lc > le) { line.clear(); break; }
        else if(lc % ls) { line.clear(); continue; }
        break;
      }
      return !line.empty();
    }

    //! read all lines and cached them into fileio::lines
    void readall();
    //! just open a file handle with the mode specified at instantiation
    bool fopen(const string& _fname);
    //! just open a file handle
    bool fopen(const string& _fname, const ios_base::openmode& mode);
    //! break line into vector of valtype 
    inline
    #ifdef __GNUC__
    __attribute__((always_inline))
    #endif
    vector<valtype> line2val() const {
      vector<valtype> ans;
      split<valtype>(ans, line);
      return ans;
    };
    //! break line into vector of string using specified deliminator 
    inline
    #ifdef __GNUC__
    __attribute__((always_inline))
    #endif
    vector<string> line2str(const string& delims=" \t") const {
      vector<string> ans;
      split<string>(ans, line, delims);
      return ans;
    }
    //! return fileio::line for read-only
    const string& getln() const;
    //! return fileio::lines for read-only
    const vector<string>& getlines() const;
    //! reset lc to 0
    void resetlc();
    string fname; //filename
    ios_base::openmode iomode; //read or write
    const bool ifperm; //if true, won't exit even though a file can't be opened
    linecounter lb; //beginning (from 1) line of input file to be actually parsed 
    linecounter ls; //parse the input file every this many of lines, starting from fileio::lb
    linecounter le; //ending line of input file to be actually parsed
    linecounter lc; //number of lines actually read, excluding comments but including lines before fileio::lb
    const string comments;
    string line; //used by fileio::readbyline()
    vector<string> lines; //used by fileio::readall()
  private:
    //! open a new file handle
    /**
     * this function will test if the old filehandle is open, if so close it;
     * it also clear fileio::line, fileio::lines
     */
    bool fopen();
    fstream fs; //fstream handle
};

#endif
