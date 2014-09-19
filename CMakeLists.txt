cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

option(MPREALCXX "Use mpfrc++" 50)

project("GWHAM")
set(CMAKE_CXX_FLAGS                "-Wall -std=c++11")
set(CMAKE_CXX_FLAGS_DEBUG          "-O0 -g")
set(CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE        "-O4 -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")

if(MPREALCXX)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lgmp -lmpfr")
  add_definitions( -DMPREALCXX=${MPREALCXX} )
endif(MPREALCXX)

add_executable(gwham gwham.cpp fileio.cpp fileio_utils.cpp gmxmdp.cpp hamiltonian.cpp ensemble.cpp) 
add_executable(mc mc.cpp fileio_utils.cpp hamiltonian.cpp ensemble.cpp)
