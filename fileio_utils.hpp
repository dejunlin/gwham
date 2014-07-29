#if !defined(FILEIO_UTILS_HPP)
#define FILEIO_UTILS_HPP

#include <sstream>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <iterator>
#include <typeinfo>
#include "exception.hpp"

using namespace std;

template <class T>
string tostr(const T& input) {
  stringstream ss;
  ss << input;
  return ss.str();
}

template <class Tp> 
void strconverter(vector<Tp>& output, const vector<string>& input) throw(FILEIO_Exception)  {
  for(uint i = 0; i < input.size(); ++i) {
    istringstream iss(input[i]);
    vector<Tp> tmp((istream_iterator<Tp>(iss)), istream_iterator<Tp>());
    if(tmp.size() == 0) {
      throw(FILEIO_Exception("Failed converting string at: '" + input[i] + "' to type '" + typeid(Tp).name() + "'"));
    } else if(tmp.size() > 1) {
      throw(FILEIO_Exception("Convertion of string at: '" + input[i] + "' to type '" + typeid(Tp).name() + "' will result in " + tostr(tmp.size()) + " elements (probably due to the presence of whitespaces)"));
    }
    output.push_back(tmp[0]);
  }
}

//split a string into an array by delims
template <class Tp> 
void parser(vector<Tp>& output, const string& str, string delims=" \t") throw(FILEIO_Exception) {
  //first split str into a array of std::string
  vector<string> tmp;
  parser<string>(tmp, str, delims);

  //then parse each element into output
  try {
    strconverter(output, tmp);
  } catch(FILEIO_Exception& fioex) {
    fioex.prepend("Error parsing string: '" + str + "' into a vector of type '" + typeid(Tp).name() + "' ");
    throw fioex;
  }
}

//specialized parser for array of string -- NOTE that delims can still be defaulted to the one in the general template
template <>
void parser<string>(vector<string>& output, const string& str, string delims) throw(FILEIO_Exception);

//case-sensitive comparison; NOTE std::string::compare() return 0 if and only if the two string match
template <class T>
bool matchkey(const string& entry, const T& key) {
  const string& keystr = tostr(key);
  return entry.compare(keystr) ? false : true; 
};

//case-insensitive comparison
template <class T>
bool cmatchkey(const string& entry, const T& key) {
  const string& keystr = tostr(key);
  return strcasecmp(entry.c_str(), keystr.c_str()) ? false : true; 
};

template <class T>
void setsc(const string& entry, T& output) {
  vector<T> outvec;
  parser(outvec, entry);
  output = outvec[0]; 
};

template <class T>
void setvec(const string& entry, vector<T>& output) throw(FILEIO_Exception) { 
  try { parser(output, entry); } 
  catch(FILEIO_Exception& fioex) { throw fioex; } 
};

#endif
