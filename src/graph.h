#pragma once

#include <iostream>
#include <vector>
#include <unordered_map>

typedef uint16_t Vertex;
typedef std::pair<Vertex,Vertex> Edge;
typedef uint16_t Edge_length;
typedef uint64_t Edge_weight;

class Graph {
    using Weight = std::pair<Edge_length, Edge_weight>;

    public:
    Graph(std::istream &input);

    void preprocess();

    private:
    std::vector<std::unordered_map<Vertex, std::vector<Weight>>> adjacency_; 
    Edge_length max_length_;
    std::vector<Vertex> terminals_;

    void add_edge(Edge edge, Weight weight);
    void remove_vertex(Vertex v);

    std::vector<Vertex> neighbors(Vertex v);

    Vertex preprocess_isolated();

};