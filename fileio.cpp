#include "typedefs.hpp"
#include "iostream"
#include <fstream>
#include <stdlib.h>
#include <limits.h>
#include <vector>
#include "fileio.hpp"
#include "fileio_utils.hpp"

fileio::fileio(
  const string& _fname,
  const ios_base::openmode& mode, 
  const bool& perm = true, 
  const linecounter _lb=0, 
  const linecounter _ls=1, 
  const linecounter _le=MAXNLINE) :
  fname(_fname),
  iomode(mode),
  ifperm(perm),
  lb(_lb),
  ls(_ls),
  le(_le),
  lc(0),
  line(string()),
  lines(vector<string>(0))
{
  fopen();
}

fileio::fileio(
  const ios_base::openmode& mode, 
  const bool& perm = true, 
  const linecounter _lb=0, 
  const linecounter _ls=1, 
  const linecounter _le=MAXNLINE) :
  fname(string()),
  iomode(mode),
  ifperm(perm),
  lb(_lb),
  ls(_ls),
  le(_le),
  lc(0),
  line(string()),
  lines(vector<string>(0))
{
}

fileio::fileio(const fileio& _fio) :
  fname(_fio.fname),
  iomode(_fio.iomode),
  ifperm(_fio.ifperm),
  lb(_fio.lb),
  ls(_fio.ls),
  le(_fio.le),
  lc(_fio.lc),
  line(_fio.line),
  lines(_fio.lines)
{
}

bool fileio::fopen(const string& _fname) {
  if(fs.is_open()) { fs.close(); }
  if(!line.empty()) { line.clear(); }
  if(!lines.empty()) { lines.clear(); }
  fname = _fname;
  return fopen();
}

bool fileio::fopen(const string& _fname, const ios_base::openmode& mode) {
  if(fs.is_open()) { fs.close(); }
  if(!line.empty()) { line.clear(); }
  if(!lines.empty()) { lines.clear(); }
  fname = _fname;
  iomode = mode;
  return fopen();
}

bool fileio::fopen() {
  fs.open(fname.c_str(),iomode);
  if(!fs.good()) {
    if(!ifperm) {
      cerr << "Error opening file: " << fname << endl;
      exit(-1);
    } else {
      return false;
    }
  }
  else {
    cout << "#Opening file: " << fname << endl;
    return true;
  }
}

bool fileio::readaline() {
  line.clear();
  while(getline(fs, line)) {
    if(line[0] == '@' || line[0] == '#' || line[0] == ';' || line.empty()) { continue; }
    ++lc;
    if(lc < lb) { line.clear(); continue; }
    if(lc % ls) { line.clear(); continue; }
    if(lc > le) { line.clear(); break; }
    break;
  }
  return !line.empty();
}

void fileio::readall() {
  line.clear();
  while(getline(fs, line)) {
    if(line[0] == '@' || line[0] == '#' || line[0] == ';' || line.empty()) { continue; }
    ++lc;
    if(lc < lb) { continue; }
    if(lc % ls) { continue; }
    if(lc > le) { break; }
    lines.push_back(line);
    line.clear();
  }
}

vector<valtype> fileio::line2val() const {
  vector<valtype> ans;
  parser<valtype>(ans, line);
  return ans;
}

vector<string> fileio::line2str() const {
  vector<string> ans;
  parser<string>(ans, line);
  return ans;
}

const vector<string>& fileio::getlines() const { return lines; }
