#!/bin/bash
export CC="gcc"
export CXX="g++"
export CPPFLAGS="-std=c++11 -g"
export LDFLAGS="-std=c++11 -g"
cmake .. -DCMAKE_BUILD_TYPE=Release 
