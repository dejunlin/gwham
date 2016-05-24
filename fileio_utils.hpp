#if !defined(FILEIO_UTILS_HPP)
#define FILEIO_UTILS_HPP

#include <sstream>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <array>
#include <iterator>
#include <typeinfo>
#include "exception.hpp"

using namespace std;
using uint = unsigned int;

template <class T>
string tostr(const T& input) {
  stringstream ss;
  ss << input;
  return ss.str();
}

template <class Tp> 
void strto(vector<Tp>& output, const vector<string>& input) {
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

//split a string into an array by empty spaces
template <class Tp> 
void split(vector<Tp>& output, const string& str) {
  istringstream iss(str);
  Tp val;
  while (iss >> val) {
    output.emplace_back(val);
  }
}

//split a string into an array by delims
template <class Tp> 
void split(vector<Tp>& output, const string& str, string delims) {
  //first split str into a array of std::string
  vector<string> tmp;
  split<string>(tmp, str, delims);

  //then parse each element into output
  try {
    strto(output, tmp);
  } catch(FILEIO_Exception& fioex) {
    fioex.prepend("Error parsing string: '" + str + "' into a vector of type '" + typeid(Tp).name() + "' ");
    throw fioex;
  }
}

//specialized split for array of string -- NOTE that delims can still be defaulted to the one in the general template
template <>
void split<string>(vector<string>& output, const string& str, string delims);

template <class Tp>
vector<Tp> split(const string& str, string delims=" \t") {
  vector<Tp> ans;
  split<Tp>(ans, str, delims);
  return ans;
}

template <class T>
void setsc(const string& entry, T& output) {
  stringstream ss(entry);
  if(!(ss >> output)) {
    throw(FILEIO_Exception("Error parsing string '" + entry + "' into type '" + typeid(T).name() + "'"));
  }
};

template <class T>
void setvec(const string& entry, vector<T>& output) {
  split(output, entry); 
};

template <class T>
T stosc(const string& input) {
  T ans;
  stringstream ss(input);
  if(!(ss >> ans)) {
    throw(FILEIO_Exception("Error parsing string '" + input + "' into type '" + typeid(T).name() + "'"));
  }
  return ans;
};

//! Get the prefix and suffix of a file name
array<string, 2> getfnfixes(const string& entry);

//! check if the input string is a comment or empty string
bool emptystr(const string& input, const string& comments);

#endif
