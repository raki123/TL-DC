#include "graph.h"
#include "parallel_search.h"

#include <iostream>

int main() {
    fpc::Graph initial_graph(std::cin);
    initial_graph.print_stats();
    initial_graph.preprocess();
    initial_graph.normalize();
    initial_graph.print_stats();
    fpc::ParallelSearch search(initial_graph.to_canon_nauty(), initial_graph.max_length(), OMP_NUM_THREADS);
    auto res = search.search();
    fpc::Edge_weight final_result = 0;
    fpc::Edge_length min_idx = initial_graph.is_all_pair() ? 3 : 0;
    for(fpc::Edge_length l = min_idx; l < res.size(); l++) {
        if(res[l] != 0) {
            std::cerr << res[l] << " paths of length " << static_cast<size_t>(l) << std::endl;
        }
        final_result += res[l];
    }
    res = initial_graph.extra_paths();
    for(fpc::Edge_length l = min_idx; l < res.size(); l++) {
        if(res[l] != 0) {
            std::cerr << res[l] << " extra paths of length " << static_cast<size_t>(l) << std::endl;
        }
        final_result += res[l];
    }
    if(initial_graph.is_all_pair()) {
        final_result /= 2;  
    }
    std::cout << final_result << std::endl;
    search.print_stats();
    return 0;
}