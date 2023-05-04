#pragma once

#include <assert.h>

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <unordered_set>

typedef uint16_t Vertex;
typedef std::pair<Vertex,Vertex> Edge;
typedef uint16_t Edge_length;
typedef uint64_t Edge_weight;
typedef std::pair<Edge_length, Edge_weight> Weight;

class Graph {

    public:
    Graph(std::istream &input);

    void preprocess();

    // encodings
    void encode_unary(std::ostream& out); 
    void encode_lenghtless(std::ostream& output);

    Edge_weight extra_paths() { return extra_paths_; }
    private:
    std::vector<std::vector<std::map<Edge_length, Edge_weight>>> adjacency_; 
    std::vector<std::set<Vertex>> neighbors_;
    Edge_length max_length_;
    std::vector<Vertex> terminals_;

    Edge_weight extra_paths_ = 0;

    void add_edge(Edge edge, Weight weight);
    void remove_edge(Edge edge);
    void remove_vertex(Vertex v);

    // helper functions
    std::set<Vertex> neighbors(Vertex v);
    void dijkstra(Vertex start, std::vector<Edge_length>& distance, const std::set<Vertex>& forbidden);

    // preprocessing subroutines
    Vertex preprocess_start_goal_edges();
    Vertex preprocess_isolated();
    Vertex preprocess_forwarder();
    Vertex preprocess_unreachable();
    Vertex preprocess_twins();
    Vertex preprocess_unusable_edge();
    Vertex preprocess_position_determined();

};