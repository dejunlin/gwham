#ifndef EXCEPTION_HPP
#define EXCEPTION_HPP

#include <exception>
#include <errno.h>
#include <string>

using namespace std;

class General_Exception : public std::exception {
  public:
    explicit General_Exception(const string& _msg) : msg(_msg) {};
    explicit General_Exception(const string& _msg, const int& _errno) : msg(_msg + " " + string(strerror(_errno))) {};
    explicit General_Exception(const int& _errno) : msg("general exception " + string(strerror(_errno))) {};
    explicit General_Exception() : msg("genearl execption") {};
    virtual ~General_Exception() throw() {};
    virtual const char* what() const throw() { return msg.c_str(); }
  protected:
    const string msg;
};

class FILEIO_Exception : public General_Exception {
  public:
    FILEIO_Exception(const string& _msg) : General_Exception("FILEIO exception: " + _msg) {}; 
};

class GMXMDP_Exception : public General_Exception {
  public:
    GMXMDP_Exception(const string& _msg) : General_Exception("GMXMDP exception: " + _msg) {}; 
};

#endif
