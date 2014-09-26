cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

option(MPREALCXX "Use mpfrc++" 50)

project("GWHAM")
set(CMAKE_CXX_FLAGS                "-Wall -std=c++11")
set(CMAKE_CXX_FLAGS_DEBUG          "-O0 -g")
set(CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE        "-O4 -DNDEBUG -ffast-math")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")

add_executable(gwham gwham.cpp fileio.cpp fileio_utils.cpp gmxmdp.cpp hamiltonian.cpp ensemble.cpp) 
add_executable(mc mc.cpp fileio_utils.cpp hamiltonian.cpp ensemble.cpp)

if(MPREALCXX)
  add_definitions( -DMPREALCXX=${MPREALCXX} )
  target_link_libraries(gwham gmp mpfr)
  target_link_libraries(mc gmp mpfr)
endif(MPREALCXX)
