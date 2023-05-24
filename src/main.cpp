#include "graph.h"
#include "parallel_search.h"

#include <iostream>

int main() {
    Graph graph(std::cin);
    graph.print_stats();
    graph.preprocess();
    graph.normalize();
    graph.print_stats();
    ParallelSearch search(graph);
    auto res = search.search();
    Edge_weight final_result = 0;
    Edge_length min_idx = graph.is_all_pair() ? 3 : 0;
    for(Edge_length l = min_idx; l < res.size(); l++) {
        if(res[l] != 0) {
            std::cerr << res[l] << " paths of length " << static_cast<size_t>(l) << std::endl;
        }
        final_result += res[l];
    }
    res = graph.extra_paths();
    for(Edge_length l = min_idx; l < res.size(); l++) {
        if(res[l] != 0) {
            std::cerr << res[l] << " extra paths of length " << static_cast<size_t>(l) << std::endl;
        }
        final_result += res[l];
    }
    if(graph.is_all_pair()) {
        final_result /= 2;  
    }
    std::cout << final_result << std::endl;
    search.print_stats();
    return 0;
}