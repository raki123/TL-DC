#include "graph.h"

#include <iostream>

int main(int argc, const char **argv) {
    Graph graph(std::cin);
    graph.preprocess();
    graph.encode_unary(std::cout);
    return 0;
}