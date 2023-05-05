cmake_minimum_required (VERSION 2.8)
project (icgca) 

option(STATIC "STATIC" OFF)
set(CMAKE_VERBOSE_MAKEFILE on)

if(STATIC)
set(CMAKE_EXE_LINKER_FLAGS -static)
set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
set(CMAKE_EXE_LINK_DYNAMIC_C_FLAGS)
set(CMAKE_EXE_LINK_DYNAMIC_CXX_FLAGS)
set(CMAKE_SHARED_LIBRARY_C_FLAGS)
set(CMAKE_SHARED_LIBRARY_CXX_FLAGS)
set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS)
set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)
set(BUILD_SHARED_LIBS OFF)
endif(STATIC)

set(CMAKE_CXX_FLAGS "-std=c++17 -O3 -Wall -Wextra -Wshadow -Wfatal-errors -Wno-sign-compare -Wno-literal-suffix -g")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.2 -mpclmul -march=native -funroll-loops -Wstrict-overflow -Wstrict-aliasing -pedantic")


SET(CMAKE_C_FLAGS "-std=c99 -Wall")
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse4.2 -mpclmul -march=native -funroll-loops -Wstrict-overflow -Wstrict-aliasing -Wextra -pedantic -Wshadow")

SET(SOURCES 
    src/graph.cpp
    src/search.cpp
	src/clhash/clhash.c
)
add_executable (read src/main.cpp ${SOURCES})
if(STATIC)
SET_TARGET_PROPERTIES (asppp PROPERTIES LINK_SEARCH_START_STATIC 1)
endif(STATIC)

