#pragma once
#include "graph.h"

class Search {
    public:
    Search(Graph& input);

    std::vector<Edge_weight> search() { return search(terminals_[0], max_length_); };
    private:
    Edge_length max_length_;
    std::vector<Vertex> terminals_;
    std::vector<std::set<Vertex>> neighbors_;
    std::vector<std::vector<std::map<Edge_length, Edge_weight>>> adjacency_; 

    Edge_length invalid_;
    std::vector<Edge_length> distance_to_goal_;


    std::vector<Edge_weight> search(Vertex start, Edge_length budget);

    // helper functions
    std::set<Vertex> neighbors(Vertex v) { assert(v >= 0 && v < neighbors_.size()); return neighbors_[v]; };
    void dijkstra(Vertex start, std::vector<Edge_length>& distance);


};