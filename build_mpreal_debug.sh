#!/bin/bash
export CC="gcc"
export CXX="g++"
cmake .. -DCMAKE_BUILD_TYPE=Debug -DMPREALCXX=30
