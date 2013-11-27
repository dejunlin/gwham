#if !defined(FILEIO_HPP)
#define FILEIO_HPP

#include "typedefs.hpp"
#include "iostream"
#include <fstream>
#include <stdlib.h>

using namespace std;

//class Taction should provide operator() (fstream& fs, Tdata& data) 
template<class Taction, class Tdata1, class Tdata2 = int>
class fileio {
  public:
    fileio(const Taction& act, const ios_base::openmode& mode);
    //! just open a file handle
    void fopen(const string& fname, fstream& fs) const;
    //! write data to fstream or read data from fstream and return the number line processed
    linecounter operator() (const string& fname, Tdata1& data1) const; 
    //! write data to fstream or read data from fstream and return the number line processed
    linecounter operator() (const string& fname, Tdata1& data1, Tdata2& data2) const; 
  private:
    const Taction funct; //what action to take (how to read/write file in what format)
    const ios_base::openmode iomode; //read or write 
};


template<class Taction, class Tdata1, class Tdata2>
fileio<Taction, Tdata1, Tdata2>::fileio(const Taction& act, const ios_base::openmode& mode) :
  funct(act),
  iomode(mode)
{}

template<class Taction, class Tdata1, class Tdata2>
void fileio<Taction, Tdata1, Tdata2>::fopen(const string& fname, fstream& fs) const {
  fs.open(fname.c_str(),iomode);
  if(!fs.good()) {
    cerr << "Error opening file: " << fname << endl;
    exit(-1);
  }
  else { cout << "#Opening file: " << fname << endl; }
}

template<class Taction, class Tdata1, class Tdata2>
linecounter fileio<Taction, Tdata1, Tdata2>::operator()(const string& fname, Tdata1& data1) const {
  fstream fs;
  fopen(fname,fs);
  const ulong nlines = funct(fs,data1);
  fs.close();
  return nlines;
}

template<class Taction, class Tdata1, class Tdata2>
linecounter fileio<Taction, Tdata1, Tdata2>::operator()(const string& fname, Tdata1& data1, Tdata2& data2) const {
  fstream fs;
  fopen(fname,fs);
  const ulong nlines = funct(fs,data1,data2);
  fs.close();
  return nlines;
}

#endif
