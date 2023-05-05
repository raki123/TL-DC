#include "graph.h"
#include "search.h"

#include <iostream>

int main() {
    Graph graph(std::cin);
    // graph.preprocess();
    std::cerr << "Extra paths: " << graph.extra_paths() << std::endl;
    graph.normalize();
    Search search(graph);
    auto res = search.search();
    Edge_weight final_result = graph.extra_paths();
    for(Edge_length l = 0; l < res.size(); l++) {
        std::cerr << res[l] << " paths of length " << l << std::endl;
        final_result += res[l];
    }
    std::cout << final_result << std::endl;
    search.print_stats();
    return 0;
}