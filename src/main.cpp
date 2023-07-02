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
	auto td = std::move(d.decompose(initial_graph));
	std::cerr << "decomposed, root " << td.first << std::endl;

	for (auto it = td.second.first.begin(); it != td.second.first.end(); ++it)
		std::cerr << "td edge " << it->first << " -> " << it->second << std::endl;

	for (auto it = td.second.second.begin(); it != td.second.second.end(); ++it)
	{	
		std::cerr << "bag " << it->first << std::endl;
		std::cerr << "edges { ";
		for (auto jt = it->second.first.begin(); jt != it->second.first.end(); ++jt)
			std::cerr << jt->first << "," << jt->second << "; ";
		std::cerr << "}" << std::endl;
	
		std::cerr << "bag { ";
		for (auto jt = it->second.second.begin(); jt != it->second.second.end(); ++jt)
			std::cerr << *jt << "; ";
		std::cerr << "}" << std::endl;
	}		

    std::vector<std::pair<Edge, std::vector<vertex_t>>> r;
    auto actual_td = td.second;
    int cur = td.first;
    while(actual_td.first.count(cur) != 0) {
        auto edges = actual_td.second[cur].first;
        auto bag = actual_td.second[cur].second;
        for(auto edge : edges) {
            r.push_back(std::make_pair(edge, bag));
        }
        cur = actual_td.first[cur];
    }
    auto edges = actual_td.second[cur].first;
    auto bag = actual_td.second[cur].second;
    for(auto edge : edges) {
        r.push_back(std::make_pair(edge, bag));
    }
    // r = {
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
    // r = {
    //     {Edge(0,1), {0,1}},
    //     {Edge(0,2), {0,1,2}},
    //     {Edge(1,3), {1,2,3}},
    //     {Edge(2,3), {2,3}}
    // };

    fpc::TreewidthSearch search(initial_graph, r, 4);
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
