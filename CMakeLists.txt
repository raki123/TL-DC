cmake_minimum_required (VERSION 3.12)

project (TLDC) 

set(CMAKE_BUILD_TYPE ReleaseWithDebugInfo)
add_compile_definitions(VERBOSE)

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

add_subdirectory("src/ankerl")

set(CMAKE_CXX_FLAGS "-O3 -fopenmp -Wall -Wextra -Wshadow -Wfatal-errors -Wno-sign-compare -Wno-literal-suffix")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.2 -mpclmul -march=native -funroll-loops -Wstrict-overflow -Wstrict-aliasing -pedantic")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")


SET(CMAKE_C_FLAGS "-O3 -fopenmp -Wall")
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse4.2 -mpclmul -march=native -funroll-loops -Wstrict-overflow -Wstrict-aliasing -Wextra -pedantic -Wshadow")
# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address")

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory("src/pathwidth")
add_subdirectory("src/dominators/Static/src")

find_library(GMP_LIB gmp)
find_library(GMPXX_LIB gmpxx)
add_library(nauty STATIC IMPORTED)
set_target_properties(nauty PROPERTIES
    IMPORTED_LOCATION ${PROJECT_SOURCE_DIR}/src/nauty2_8_6/nauty.a
)
add_library(jemalloc STATIC IMPORTED)
set_target_properties(jemalloc PROPERTIES
    IMPORTED_LOCATION ${PROJECT_SOURCE_DIR}/src/jemalloc/lib/libjemalloc.a
)


SET(SOURCES 
    src/graph.cpp
    src/digraph.cpp
    src/branch_decomposition.cpp
    src/tree_decomposition.cpp
    src/annotated_decomposition.cpp
    src/types.cpp
    src/clhasher.cpp
    src/meet_in_the_middle_search.cpp
    src/meet_in_the_middle_directed_search.cpp
    src/parallel_nauty_search.cpp
    src/parallel_directed_search.cpp
    src/parallel_search.cpp
    src/treewidth_search.cpp
    src/treewidth_directed_search.cpp
    src/nauty_pathwidth_search.cpp
    src/clhash/clhash.c
    src/canonizer.cpp
)
add_executable (tldc src/main.cpp ${SOURCES})
if(STATIC)
    SET_TARGET_PROPERTIES (tldc PROPERTIES LINK_SEARCH_START_STATIC 1)
endif(STATIC)
target_include_directories(tldc PUBLIC
                          "${PROJECT_SOURCE_DIR}/src"
                          "${PROJECT_SOURCE_DIR}/src/pathwidth"
                          "${PROJECT_SOURCE_DIR}/src/ankerl/include"
                          )
target_link_libraries (tldc nauty jemalloc ${GMP_LIB} ${GMPXX_LIB} libpathwidth libstaticdominators)
