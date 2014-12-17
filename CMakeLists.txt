cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)

project("GWHAM")

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.9)
        message(FATAL_ERROR "GCC version must be at least 4.9!")
  endif()
endif ()

#if we use mpfrc++ for multiple precision math
option(MPREALCXX "Use mpfrc++" 50)

set(CMAKE_CXX_FLAGS                "-Wall")
set(CMAKE_CXX_FLAGS_DEBUG          "-O0 -g")
set(CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE        "-O4 -DNDEBUG -ffast-math")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")

add_executable(gwham gwham.cpp fileio.cpp fileio_utils.cpp gmxmdp.cpp hamiltonian.cpp ensemble.cpp) 
add_executable(mc mc.cpp fileio_utils.cpp hamiltonian.cpp ensemble.cpp)

#make sure we use c++11 here
set_property(TARGET gwham mc PROPERTY CXX_STANDARD 11)
set_property(TARGET gwham mc PROPERTY CXX_STANDARD_REQUIRED ON)

if(MPREALCXX)
  add_definitions( -DMPREALCXX=${MPREALCXX} )
  target_link_libraries(gwham gmp mpfr)
  target_link_libraries(mc gmp mpfr)
endif(MPREALCXX)
