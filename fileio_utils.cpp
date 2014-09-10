/*
 * =====================================================================================
 *
 *       Filename:  fileio_utils.cpp
 *
 *    Description:  fileio utility functions
 *
 *        Version:  1.0
 *        Created:  29/07/14 10:47:57
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dejun Lin (DL), dejun.lin@gmail.com
 *   Organization:  Department of Biochemistry and Biophysics, Medical Center, University of Rochester
 *
 * =====================================================================================
 */
#include "fileio_utils.hpp"
#include <regex>
#include <array>

//specialized split for array of string -- NOTE that delims can still be defaulted to the one in the general template
template <>
void split<string>(vector<string>& output, const string& str, string delims) throw(FILEIO_Exception) {
  char * str_c = new char[str.size()+1];
  str_c[str.size()] = '\0';
  memcpy(str_c,str.c_str(),str.size());

  char * pch = strtok(str_c, delims.c_str());
  while(pch !=NULL) {
    istringstream iss(pch);
    output.push_back(iss.str());
    pch = strtok(NULL, delims.c_str());
  }
  delete[] str_c;
}

array<string, 2> getfnfixes(const string& entry) {
  const regex fnamebase(R"(^(.*)[.](.*?)$)");
  auto i = sregex_token_iterator(entry.begin(), entry.end(), fnamebase, {1,2});
  array<string, 2> ans;
  ans[0] = *i++;
  ans[1] = *i;
  return ans;
}

bool emptystr(const string& input, const string& comments) {
  const regex emptyregex(R"(^\s*[)"+comments+R"(]+.*\n*$|^\s*\n*$)");
  return regex_match(input, emptyregex);
}
