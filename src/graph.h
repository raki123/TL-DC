#pragma once

#include <assert.h>

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <unordered_set>
#include <queue>

typedef uint16_t Vertex;
typedef std::pair<Vertex,Vertex> Edge;
typedef uint16_t Edge_length;
typedef uint64_t Edge_weight;
typedef std::pair<Edge_length, Edge_weight> Weight;
typedef std::priority_queue<std::pair<Edge_length, Vertex>, std::vector<std::pair<Edge_length, Vertex>>, std::greater<std::pair<Edge_length, Vertex>>> DijkstraQueue;

class Graph {

    public:
    friend class Search;
    Graph(std::istream &input);

    void preprocess();

    void print_stats();

    std::vector<Edge_weight> extra_paths() { return extra_paths_; }

    void normalize();

    bool is_all_pair() { return all_pair_; }
    private:
    Edge_length max_length_;
    std::vector<std::set<Vertex>> neighbors_;
    std::vector<std::vector<std::map<Edge_length, Edge_weight>>> adjacency_; 
    std::vector<Vertex> terminals_;

    // if vertex v is included, then all vertices in exclude_[v] must be excluded
    std::vector<std::vector<Vertex>> exclude_;

    std::vector<Edge_weight> extra_paths_;

    bool all_pair_ = false;

    void add_edge(Edge edge, Weight weight);
    void add_exclude(Vertex v, Vertex w);
    void remove_edge(Edge edge);
    void remove_vertex(Vertex v);
    std::set<Vertex> neighbors(Vertex v);
    // check returns true if a vertex is in an exclusion constraint
    bool fixed(Vertex v) { return !exclude_[v].empty(); }

    Graph(Vertex n);
    Graph subgraph(std::vector<Vertex> restrict_to);

    // utility functions
    void dijkstra(Vertex start, std::vector<Edge_length>& distance, const std::set<Vertex>& forbidden);
    std::vector<std::vector<Vertex>> components(const std::set<Vertex>& forbidden);
    std::vector<Vertex> find_separator(size_t size, size_t min_component_size, bool terminals_in_same = true);
    
    // preprocessing subroutines
    Vertex preprocess_start_goal_edges();
    Vertex preprocess_isolated();
    Vertex preprocess_forwarder();
    Vertex preprocess_unreachable();
    Vertex preprocess_twins();
    Vertex preprocess_unusable_edge();
    Vertex preprocess_position_determined();
    Vertex preprocess_two_separator();


};