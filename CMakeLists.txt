cmake_minimum_required (VERSION 3.1)

project (TLDC) 

set(CMAKE_BUILD_TYPE Release)

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

add_subdirectory("src/clingo")

set(CMAKE_CXX_FLAGS "-std=c++17 -O3 -fopenmp -Wall -Wextra -Wshadow -Wfatal-errors -Wno-sign-compare -Wno-literal-suffix")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.2 -mpclmul -march=native -funroll-loops -Wstrict-overflow -Wstrict-aliasing -pedantic")

SET(CMAKE_C_FLAGS "-std=c11 -O3 -fopenmp -Wall")
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse4.2 -mpclmul -march=native -funroll-loops -Wstrict-overflow -Wstrict-aliasing -Wextra -pedantic -Wshadow")

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_library(GMP_LIB gmp)
find_library(GMPXX_LIB gmpxx)
add_library(nauty STATIC IMPORTED)
set_target_properties(nauty PROPERTIES
    IMPORTED_LOCATION ${PROJECT_SOURCE_DIR}/src/nauty2_8_6/nauty.a
)

SET(SOURCES 
    src/graph.cpp
    src/popen2.cpp
    src/Decomposer.cpp
    src/search.cpp
    src/parallel_search.cpp
    src/treewidth_search.cpp
    src/nauty_pathwidth_search.cpp
    src/clhash/clhash.c
)
add_executable (tldc src/main.cpp ${SOURCES})
if(STATIC)
SET_TARGET_PROPERTIES (tldc PROPERTIES LINK_SEARCH_START_STATIC 1)
endif(STATIC)
target_link_libraries (tldc libclingo nauty ${GMP_LIB} ${GMPXX_LIB})
