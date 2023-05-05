#include "graph.h"
#include "search.h"

#include <iostream>

int main(int argc, const char **argv) {
    Graph graph(std::cin);
    graph.preprocess();
    std::cerr << "Extra paths: " << graph.extra_paths() << std::endl;
    graph.normalize();
    Search search(graph);
    auto res = search.search();
    for(Edge_length l = 0; l < res.size(); l++) {
        std::cout << res[l] << " paths of length " << l << std::endl;
    }
    return 0;
}