#!/bin/bash
export CC="gcc"
export CXX="g++"
cmake .. -DCMAKE_BUILD_TYPE=Release -DMPREALCXX=20
