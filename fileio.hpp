#if !defined(FILEIO_HPP)
#define FILEIO_HPP

#include "typedefs.hpp"
#include "iostream"
#include <fstream>
#include <stdlib.h>

using namespace std;

//class Taction should provide operator() (fstream& fs, Tdata& data) 
template<class Taction, class Tdata>
class fileio {
  public:
    fileio(const Taction& act, const ios_base::openmode& mode);
    //! write data to fstream or read data from fstream and return the number line processed
    linecounter operator() (Tdata& data, const string& fname) const; 
  private:
    const Taction funct; //what action to take (how to read/write file in what format)
    const ios_base::openmode iomode; //read or write 
};


template<class Taction, class Tdata>
fileio<Taction, Tdata>::fileio(const Taction& act, const ios_base::openmode& mode) :
  funct(act),
  iomode(mode)
{}

template<class Taction, class Tdata>
linecounter fileio<Taction, Tdata>::operator()(Tdata& data, const string& fname) const {
  fstream fs(fname.c_str(),iomode);
  if(!fs.good()) {
    cerr << "Error opening file: " << fname << endl;
    exit(-1);
  }
  else { cout << "#Opening file: " << fname << endl; }
  const ulong nlines = funct(fs,data);
  fs.close();
  return nlines;
}

#endif
