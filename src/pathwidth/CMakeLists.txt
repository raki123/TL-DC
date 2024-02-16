add_library(libpathwidth heuristics.cpp)
target_include_directories(libpathwidth PUBLIC .)


SET(PD_SOURCES 
    heuristics.cpp
    path_decomposition_main.cpp
)
add_executable (path_decompose ${PD_SOURCES})
if(STATIC)
    SET_TARGET_PROPERTIES (path_decompose PROPERTIES LINK_SEARCH_START_STATIC 1)
endif(STATIC)