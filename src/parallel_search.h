#pragma once
#include "graph.h"
#include "clhash/clhash.h"
#include <unordered_map>
#include <utility>

typedef std::vector<char> CacheKey;

extern clhasher hasher__;


struct vector_hash {
public:
    size_t operator()(const CacheKey &key) const {
        return hasher__(key);
    }  
};


class ParallelSearch {
    public:
    ParallelSearch(Graph& input);

    std::vector<Edge_weight> search();

    // void print_stats();
    private:
    bool enable_dag_;
    Edge_length max_length_;
    std::vector<Vertex> terminals_;
    std::vector<std::vector<Vertex>> neighbors_;
    std::vector<std::vector<std::vector<std::pair<Edge_length, Edge_weight>>>> adjacency_; 

    Edge_length invalid_;
    std::vector<Edge_length> distance_to_goal_;
    std::vector<std::vector<Edge_length>> distance_;

    // // if vertex v is included, then all vertices in exclude_[v] must be excluded
    // std::vector<std::vector<Vertex>> exclude_;

    std::vector<char> visited_;

    std::vector<std::vector<std::unordered_map<CacheKey, std::vector<Edge_weight>, vector_hash>>> cache_; 

    // // recursive search
    // std::vector<Edge_weight> search(Vertex start, Edge_length budget);
    // // more efficient search if budget is equal to length of shortest path
    // std::vector<Edge_weight> dag_search(Vertex start, Edge_length budget);



    void prune_articulation(Vertex start);
    bool ap_util(Vertex u, std::vector<char>& visited, std::vector<Vertex>& disc, std::vector<Vertex>& low, int& time, int parent, Vertex start);
    void prune_util(Vertex u, std::vector<char>& visited);
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
    void pruning_dijkstra(Vertex start, Vertex prune, std::vector<Edge_length>& distance, std::vector<char>& visited, Edge_length budget);

    // // stats
    // size_t pos_hits = 0;
    // size_t neg_hits = 0;

    // size_t splits = 0;

    // size_t edges = 0;
    // size_t propagations = 0;
    // size_t dags = 0;


};