SET(SOURCES 
    dgraph_iter.cpp
    dgraph_lt.cpp
    dgraph_sdom.cpp
    dgraph_slt.cpp
    dgraph_snca.cpp
    dgraph.cpp
    rfw_timer.cpp
)
add_library(libstaticdominators ${SOURCES})
target_include_directories(libstaticdominators PUBLIC .)