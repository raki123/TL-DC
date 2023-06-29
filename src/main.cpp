#include "graph.h"
#include "treewidth_search.h"

#include <iostream>
#include "Decomposer.hpp"

int main() {
    fpc::Graph initial_graph(std::cin);
    // initial_graph.normalize();
    initial_graph.print_stats();

    std::vector<std::pair<Edge, std::vector<vertex_t>>> r;
    r = {
        {Edge(0,1), {0,1}},
        {Edge(1,2), {0,1,2}},
        {Edge(0,3), {0,1,2,3}},
        {Edge(1,4), {1,2,3,4}},
        {Edge(2,5), {2,3,4,5}},
        {Edge(3,4), {3,4,5}},
        {Edge(4,5), {3,4,5}},
        {Edge(3,6), {3,4,5,6}},
        {Edge(4,7), {4,5,6,7}},
        {Edge(5,8), {5,6,7,8}},
        {Edge(6,7), {6,7,8}},
        {Edge(7,8), {7,8}}
    };
	// Decomposer d;
	// auto r = std::move(d.decompose(initial_graph));
	for (auto it = r.begin(); it != r.end(); ++it)
	{
		std::cout << "edge " << it->first.first << "," << it->first.second << "; bag ";
		for (auto jt = it->second.begin(); jt != it->second.end(); jt++)
			std::cout << *jt << " ";
		std::cout << std::endl;
	}

    fpc::TreewidthSearch search(initial_graph, r);
    auto res = search.search();
    fpc::Edge_weight final_result = 0;
    fpc::Edge_length min_idx = initial_graph.is_all_pair() ? 0 : 0;
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
    return 0;
}
