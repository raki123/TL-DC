#include "graph.h"
#include "treewidth_search.h"

#include <iostream>
#include "Decomposer.hpp"

int main() {
    fpc::Graph initial_graph(std::cin);
    // initial_graph.preprocess();
    // initial_graph.normalize();
    initial_graph.print_stats();

    Decomposer d;
    std::vector<std::pair<Edge, std::vector<vertex_t>>> r = std::move(d.path_decompose(initial_graph));
    auto r2 = std::move(d.tree_decompose(initial_graph));


    //std::vector<std::pair<Edge, std::vector<vertex_t>>> r = std::move(d.tree_decompose(initial_graph));


    // fpc::AnnotatedDecomposition example = {
    //     {Edge(0,1), {0,1}},
    //     {Edge(1,2), {0,1,2}},
    //     {Edge(0,3), {0,1,2,3}},
    //     {Edge(1,4), {1,2,3,4}},
    //     {Edge(2,5), {2,3,4,5}},
    //     {Edge(3,4), {3,4,5}},
    //     {Edge(4,5), {3,4,5}},
    //     {Edge(3,6), {3,4,5,6}},
    //     {Edge(4,7), {4,5,6,7}},
    //     {Edge(5,8), {5,6,7,8}},
    //     {Edge(6,7), {6,7,8}},
    //     {Edge(7,8), {7,8}}
    // };
    // size_t invalid = size_t(-1);
    // fpc::AnnotatedDecomposition example = {
    //     {2, fpc::NodeType::LEAF, Edge(0,1), std::make_pair(invalid, invalid), {0,1}},
    //     {2, fpc::NodeType::LEAF, Edge(0,2), std::make_pair(invalid, invalid), {0,2}},
    //     {3, fpc::NodeType::JOIN, Edge(invalid,invalid), std::make_pair(1, 0), {0,1,2}},
    //     {4, fpc::NodeType::PATH_LIKE, Edge(1,3), std::make_pair(3, invalid), {1,2,3}},
    //     {invalid, fpc::NodeType::PATH_LIKE, Edge(2,3), std::make_pair(2, invalid), {2,3}},
    // };

    fpc::TreewidthSearch search(initial_graph, r2, 4);
    auto res = search.search();
    search.print_stats();
    fpc::Edge_weight final_result = 0;
    for(fpc::Edge_length l = 0; l < res.size(); l++) {
        if(res[l] != 0) {
            std::cerr << res[l] << " paths of length " << static_cast<size_t>(l) << std::endl;
        }
        final_result += res[l];
    }
    res = initial_graph.extra_paths();
    for(fpc::Edge_length l = 0; l < res.size(); l++) {
        if(res[l] != 0) {
            std::cerr << res[l] << " extra paths of length " << static_cast<size_t>(l) << std::endl;
        }
        final_result += res[l];
    }
    std::cout << final_result << std::endl;
    return 0;
}
