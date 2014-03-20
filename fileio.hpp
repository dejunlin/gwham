#if !defined(FILEIO_HPP)
#define FILEIO_HPP

#include "typedefs.hpp"
#include <fstream>
#include <vector>
#include "fileio_utils.hpp"

using namespace std;

class fileio {
  public:
    fileio(const string& _fname, const ios_base::openmode& mode, const bool& perm, const linecounter _lb, const linecounter _ls, const linecounter _le);
    fileio(const ios_base::openmode& mode, const bool& perm, const linecounter _lb, const linecounter _ls, const linecounter _le);
    fileio(const fileio& _fio);
    //! read one line at a time and cache it into fileio::line 
    bool readaline(); 
    //! read all lines and cached them into fileio::lines
    void readall();
    //! just open a file handle with the mode specified at instantiation
    bool fopen(const string& _fname);
    //! just open a file handle
    bool fopen(const string& _fname, const ios_base::openmode& mode);
    //! break line into vector of valtype 
    vector<valtype> line2val() const;
    //! break line into vector of string 
    vector<string> line2str() const;
    //! return fileio::lines for read-only
    const vector<string>& getlines() const;
    string fname; //filename
    ios_base::openmode iomode; //read or write
    const bool ifperm; //if true, won't exit even though a file can't be opened
    const linecounter lb; //beginning (from 1) line of input file to be actually parsed 
    const linecounter ls; //parse the input file every this many of lines, starting from fileio::lb
    const linecounter le; //ending line of input file to be actually parsed
    linecounter lc; //number of lines actually read, excluding comments but including lines before fileio::lb
    string line; //used by fileio::readbyline()
    vector<string> lines; //used by fileio::readall()
  private:
    //! just open a file handle
    bool fopen();
    fstream fs; //fstream handle
};

#endif
