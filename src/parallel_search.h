#pragma once
#include "graph.h"
#include "search.h"
#include <unordered_map>
#include <utility>
#include <omp.h>

namespace fpc {

typedef sparsegraph PCacheKey;

struct sg_hash {
public:
    size_t operator()(PCacheKey const& key) const {
        return hasher__((const char *)key.v, sizeof(size_t)*key.nv + sizeof(int)*(key.nv + key.nde));
    }  
};
struct sg_equal {
public:
    bool operator()(PCacheKey const& lhs, PCacheKey const& rhs) const {
        return lhs.nv == rhs.nv && lhs.nde == rhs.nde && std::memcmp(lhs.v, rhs.v, sizeof(size_t)*lhs.nv + sizeof(int)*(lhs.nv + lhs.nde)) == 0;
    }
};

class ParallelSearch {
    public:
    ParallelSearch(sparsegraph input, Edge_length max_length);

    std::vector<Edge_weight> search();

    void print_stats();
    private:
    size_t nthreads_;
    Edge_length max_length_;
    Edge_length invalid_;

    sparsegraph initial_;

    std::vector<std::unordered_map<PCacheKey, std::vector<Edge_weight>, sg_hash, sg_equal>> cache_; 

    std::vector<Edge_weight> result_;
    std::vector<std::vector<Edge_weight>> thread_local_result_;

    // // recursive search
    // std::vector<Edge_weight> search(Vertex start, Edge_length budget);
    // more efficient search if budget is equal to length of shortest path
    Edge_weight dag_search(sparsegraph const& sg, Vertex start, std::vector<Edge_length> const& distance_to_goal);



    void prune_articulation(sparsegraph const& sg, Vertex start, std::vector<Edge_length>& distance);
    bool ap_util(sparsegraph const& sg, Vertex u, std::vector<char>& visited, std::vector<Vertex>& disc, std::vector<Vertex>& low, int& time, int parent, Vertex start, std::vector<Edge_length>& distance);
    void prune_util(sparsegraph const& sg, Vertex u, std::vector<Edge_length>& distance);
    // void component_util(Vertex u, std::vector<Vertex>& disc);

    // // ap datastructures
    // std::vector<Vertex> ap_disc_;
    // std::vector<Vertex> ap_low_;
    // std::vector<char> ap_visited_;

    // // for splitting based on articulation points between start and goal
    // Vertex last_ap_;
    // std::vector<std::pair<Vertex, Vertex>> ap_start_goal_; 
    // std::vector<std::vector<Vertex>> ap_components_;

    // helper functions
    void pruning_dijkstra(sparsegraph const& sg, Vertex prune, std::vector<Edge_length>& distance, Edge_length budget);
    void reverse_pruning_dijkstra(sparsegraph const& sg, Vertex prune, std::vector<Edge_length>& distance, std::vector<Edge_length> const& forward_distance, Edge_length budget);

    // stats
    std::vector<size_t> pos_hits_;
    std::vector<size_t> neg_hits_;

    // size_t splits = 0;

    std::vector<size_t> edges_;
    std::vector<size_t> propagations_;
    std::vector<size_t> dags_;

};

} // namespace fpc