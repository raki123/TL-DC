#include "graph.h"

#include <iostream>

int main(int argc, const char **argv) {
    Graph graph(std::cin);
    graph.preprocess();
    return 0;
}