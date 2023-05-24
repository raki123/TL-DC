#pragma once
#include "graph.h"
#include "search.h"
#include <unordered_map>
#include <utility>
#include <omp.h>

typedef std::vector<Edge_length> PCacheKey;

struct pvector_hash {
public:
    size_t operator()(const PCacheKey &key) const {
        return hasher__(key);
    }  
};

class ParallelSearch {
    public:
    ParallelSearch(Graph& input);

    std::vector<Edge_weight> search();

    void print_stats();
    private:
    size_t nthreads_;
    bool enable_dag_;
    Edge_length max_length_;
    std::vector<Vertex> terminals_;
    std::vector<std::vector<Vertex>> neighbors_;
    std::vector<std::vector<std::vector<std::pair<Edge_length, Edge_weight>>>> adjacency_; 

    Edge_length invalid_;
    std::vector<std::vector<Edge_length>> distance_;

    // // if vertex v is included, then all vertices in exclude_[v] must be excluded
    // std::vector<std::vector<Vertex>> exclude_;

    std::vector<char> visited_;

    std::vector<std::vector<std::unordered_map<PCacheKey, std::vector<Edge_weight>, pvector_hash>>> cache_; 

    std::vector<Edge_weight> result_;
    std::vector<std::vector<Edge_weight>> thread_local_result_;

    // // recursive search
    // std::vector<Edge_weight> search(Vertex start, Edge_length budget);
    // // more efficient search if budget is equal to length of shortest path
    // std::vector<Edge_weight> dag_search(Vertex start, Edge_length budget);



    void prune_articulation(Vertex start, std::vector<Edge_length>& distance);
    bool ap_util(Vertex u, std::vector<char>& visited, std::vector<Vertex>& disc, std::vector<Vertex>& low, int& time, int parent, Vertex start, std::vector<Edge_length>& distance);
    void prune_util(Vertex u, std::vector<Edge_length>& distance);
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
    std::vector<Vertex> neighbors(Vertex v) { assert(v >= 0 && v < neighbors_.size()); return neighbors_[v]; };
    void dijkstra(Vertex start, std::vector<Edge_length>& distance);
    void pruning_dijkstra(Vertex start, Vertex prune, std::vector<Edge_length>& distance, std::vector<Edge_length> const& old_distance, Edge_length budget);

    // stats
    std::vector<size_t> pos_hits_;
    std::vector<size_t> neg_hits_;

    // size_t splits = 0;

    std::vector<size_t> edges_;
    std::vector<size_t> propagations_;
    std::vector<size_t> dags_;

};