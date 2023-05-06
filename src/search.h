#pragma once
#include "graph.h"
#include "clhash/clhash.h"
#include <unordered_map>

typedef std::vector<Edge_length> CacheKey;

extern clhasher hasher__;

namespace std {
    template<>
    class hash<CacheKey> {
    public:
        size_t operator()(const CacheKey &key) const {
            return hasher__(key);
        }  
    };
}

class Search {
    public:
    Search(Graph& input);

    std::vector<Edge_weight> search() { return search(terminals_[0], max_length_); };

    void print_stats();
    private:
    Edge_length max_length_;
    std::vector<Vertex> terminals_;
    std::vector<std::set<Vertex>> neighbors_;
    std::vector<std::vector<std::map<Edge_length, Edge_weight>>> adjacency_; 

    Edge_length invalid_;
    std::vector<Edge_length> distance_to_goal_;
    std::vector<std::vector<Edge_length>> distance_;

    std::vector<bool> visited_;

    std::vector<std::unordered_map<CacheKey, std::pair<Edge_length, std::vector<Edge_weight>>>> cache_; 

    std::vector<Edge_weight> search(Vertex start, Edge_length budget);
    std::vector<Edge_weight> dag_search(Vertex start, Edge_length budget);

    // helper functions
    std::set<Vertex> neighbors(Vertex v) { assert(v >= 0 && v < neighbors_.size()); return neighbors_[v]; };
    void dijkstra(Vertex start, std::vector<Edge_length>& distance, Edge_length budget);
    void pruning_dijkstra(Vertex start, Vertex prune, std::vector<Edge_length>& distance, Edge_length budget);

    // stats
    size_t pos_hits = 0;
    size_t neg_hits = 0;

    size_t edges = 0;
    size_t dags = 0;


};