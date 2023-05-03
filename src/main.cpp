#include "graph.h"

#include <iostream>

int main(int argc, const char **argv) {
    Graph graph(std::cin);
    graph.preprocess();
    std::cerr << "Extra paths: " << graph.extra_paths() << std::endl;
    // graph.encode_unary(std::cout);
    return 0;
}