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
  const bool& perm, 
  const linecounter _lb, 
  const linecounter _ls, 
  const linecounter _le,
  const string& cm
  ) :
  fname(_fname),
  iomode(mode),
  ifperm(perm),
  lb(_lb),
  ls(_ls),
  le(_le),
  lc(0),
  comments(cm),
  line(string()),
  lines(vector<string>(0))
{
  fopen();
}

fileio::fileio(
  const ios_base::openmode& mode, 
  const bool& perm, 
  const linecounter _lb, 
  const linecounter _ls, 
  const linecounter _le,
  const string& cm
  ) :
  fname(string()),
  iomode(mode),
  ifperm(perm),
  lb(_lb),
  ls(_ls),
  le(_le),
  lc(0),
  comments(cm),
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
  comments(_fio.comments),
  line(_fio.line),
  lines(_fio.lines)
{
}

bool fileio::fopen(const string& _fname) {
  fname = _fname;
  return fopen();
}

bool fileio::fopen(const string& _fname, const ios_base::openmode& mode) {
  fname = _fname;
  iomode = mode;
  return fopen();
}

bool fileio::fopen() {
  if(fs.is_open()) { fs.close(); }
  if(!line.empty()) { line.clear(); }
  if(!lines.empty()) { lines.clear(); }
  fs.open(fname.c_str(),iomode);
  if(!fs.good()) {
    if(!ifperm) {
      throw(FILEIO_Exception("Error opening file: "+fname));
    } else {
      return false;
    }
  }
  else {
    cout << "#Opening file: " << fname << endl;
    return true;
  }
}

bool fileio::emptyline() const {
  return emptystr(line, comments); 
}

bool fileio::readaline() {
  line.clear();
  while(getline(fs, line)) {
    if(emptyline()) { continue; }
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
    if(emptyline()) { continue; }
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

vector<string> fileio::line2str(const string& delims) const {
  vector<string> ans;
  parser<string>(ans, line, delims);
  return ans;
}

const vector<string>& fileio::getlines() const { return lines; }
const string& fileio::getln() const { return line; }

void fileio::resetlc() { lc = 0; }
