#include "graph.h"
#include "treewidth_search.h"
#include "nauty_pathwidth_search.h"
#include "parallel_search.h"

#include <iostream>
#include "Decomposer.hpp"

void details(size_t bags, size_t nr_bags, size_t joinch, size_t max_join, size_t nr_joins)
{
	std::cerr << "max bag size " << bags << ", nr bags " << nr_bags << ", max join child bags " << joinch << ", max join bag size " << max_join << ", nr of joins " << nr_joins << std::endl;
}

int main() {
    fpc::Graph initial_graph(std::cin);
    initial_graph.preprocess();
    initial_graph.normalize();
    initial_graph.print_stats();



    // portfolio parameters
    double nr_automporphisms = initial_graph.nr_automorphisms();
    Edge_length max_length = initial_graph.max_length();
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


    bool is_all_pair = initial_graph.is_all_pair(), use_pw = false;
    
    AnnotatedDecomposition* r = &rt;

    if ((!is_all_pair && max_bagsize < 19) || (is_all_pair && max_bagsize <= 22))
    {
    	use_pw = true;
        r = &rp;
	std::cerr << "fallback to PD " << std::endl;
    }


    auto &r2 = *r;


    std::vector<Edge_weight> res;
    if((!is_all_pair && max_bagsize >= 16 && max_bagsize <= 25 && nr_automporphisms >= 1e+11 && nr_automporphisms <= 1e+45) || 
    	(is_all_pair && use_pw && max_bagsize >= 21)) {
        //fpc::NautyPathwidthSearch search(initial_graph, r2, 4);
        fpc::NautyPathwidthSearch search(initial_graph, rp, 4);
        res = search.search();
        search.print_stats();
    }
    else if (is_all_pair || use_pw || (max_length > 59 && std::max(max_join_child,t_max_bagsize) < 23) || std::max(max_join_child,t_max_bagsize) < 21) {
        fpc::TreewidthSearch search(initial_graph, r2, 4);
        res = search.search();
        search.print_stats();
    } else {
        fpc::ParallelSearch search(initial_graph.to_canon_nauty(true), initial_graph.max_length(), 4);
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
    std::cout << final_result << std::endl;
    return 0;
}
