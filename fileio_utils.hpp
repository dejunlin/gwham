#if !defined(FILEIO_UTILS_HPP)
#define FILEIO_UTILS_HPP

#include <sstream>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <vector>
#include "exception.hpp"

using namespace std;

template <class Tp> 
void parser(vector<Tp>& output, const string& str, string delims=" \t") throw(FILEIO_Exception) {
	char * str_c = new char[str.size()+1];
	str_c[str.size()] = '\0';
	memcpy(str_c,str.c_str(),str.size());

	char * pch = strtok(str_c,delims.c_str());
	while(pch !=NULL) {
		istringstream iss(pch);
		Tp tok;
		if(!(iss>>tok)) {
		  delete[] str_c;
		  throw(FILEIO_Exception("Error parsing string: '" + str + "'"));
		}
		output.push_back(tok);
		pch = strtok(NULL,delims.c_str());
	}
	delete[] str_c;
}

template <class Tp> 
void parser(vector<Tp>& output, const vector<string>& input) {
  for(uint i = 0; i < input.size(); ++i) {
    output.push_back(atof(input[i].c_str()));
  }
}

template <class T>
string tostr(const T& input) {
  stringstream ss;
  ss << input;
  return ss.str();
}

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
