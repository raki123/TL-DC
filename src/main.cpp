#include "graph.h"
#include "treewidth_search.h"
#include "parallel_search.h"

#include <iostream>
#include <chrono>
#include "Decomposer.hpp"

void details(size_t bags, size_t nr_bags, size_t joinch, size_t max_join, size_t nr_joins)
{
	std::cerr << "max bag size " << bags << ", nr bags " << nr_bags << ", max join child bags " << joinch << ", max join bag size " << max_join << ", nr of joins " << nr_joins << std::endl;
}

int main() {
    auto start = std::chrono::system_clock::now();
    fpc::Graph initial_graph(std::cin);
    initial_graph.preprocess();
    initial_graph.normalize();
    initial_graph.print_stats();


    size_t max_bagsize, nr_bags;
    size_t t_max_bagsize, max_join_child, max_join, nr_joins, t_nr_bags;

    Decomposer d;
    // std::vector<std::pair<Edge, std::vector<vertex_t>>> r = std::move(d.path_decompose(initial_graph));
    //path
    auto rp = std::move(d.tree_decompose(initial_graph, true, &max_bagsize, &nr_bags));

    std::cerr << "PATH: ";
    details(max_bagsize, nr_bags, 0, 0, 0);

    //tree
    auto rt = std::move(d.tree_decompose(initial_graph, false, &t_max_bagsize, &t_nr_bags, &max_join_child, &max_join, &nr_joins));

    std::cerr << "TREE: ";
    details(t_max_bagsize, t_nr_bags, max_join_child, max_join, nr_joins);


    AnnotatedDecomposition* r = &rt;

    if (true || max_bagsize <= max_join_child + 1)
    {
        r = &rp;
	std::cerr << "fallback to PD " << std::endl;
    }


    auto &r2 = *r;

    std::vector<Edge_weight> res;
    bool use_treewidth = initial_graph.is_all_pair();
    if(true || use_treewidth) {
        fpc::TreewidthSearch search(initial_graph, r2, 1);
        res = search.search();
        search.print_stats();
    } else {
        fpc::ParallelSearch search(initial_graph.to_canon_nauty(), initial_graph.max_length(), 4);
        res = search.search();
        search.print_stats();
    }
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
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now() - start);
    std::cout << final_result << " " << duration.count() << std::endl;
    return 0;
}
