// TLDC "Too long; Didn't Count" A length limited path counter.
// Copyright (C) 2023 Rafael Kiesel, Markus Hecher

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

#pragma once
extern "C" {
#include "nauty2_8_6/nausparse.h"
}
#include <assert.h>

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <unordered_set>
#include <queue>


class Decomposer;

namespace fpc {

typedef uint16_t Vertex;
typedef std::pair<Vertex,Vertex> Edge;
typedef uint8_t Edge_length;
typedef uint64_t Edge_weight;
typedef std::pair<Edge_length, Edge_weight> Weight;
typedef std::priority_queue<std::pair<Edge_length, Vertex>, std::vector<std::pair<Edge_length, Vertex>>, std::greater<std::pair<Edge_length, Vertex>>> DijkstraQueue;

class Graph {

    public:
    friend class Search;
    friend class ParallelSearch;
    friend class TreewidthSearch;
    friend class NautyPathwidthSearch;
    friend class ::Decomposer;
    Graph(std::istream &input);

    void preprocess();

    void print_stats();

    std::vector<Edge_weight> extra_paths() { return extra_paths_; }

    void normalize(bool reorder = false);

    bool is_all_pair() { return all_pair_; }

    Edge_length max_length() { return max_length_; }
    Edge_length min_length();

    std::vector<vertex_t> terminals() { return terminals_; }

    size_t nr_vertices();
    size_t nr_edges();

    sparsegraph to_canon_nauty(bool reorder);
    double nr_automorphisms();
    private:
    Edge_length max_length_;
    std::vector<std::set<Vertex>> neighbors_;
    std::vector<std::vector<std::map<Edge_length, Edge_weight>>> adjacency_; 
    std::vector<Vertex> terminals_;

    // if vertex v is included, then all vertices in exclude_[v] must be excluded
    std::vector<std::set<Vertex>> exclusion_classes_;
    std::vector<size_t> exclude_;

    std::vector<Edge_weight> extra_paths_;

    bool all_pair_ = false;

    void add_edge(Edge edge, Weight weight);
    void add_exclude(Vertex v, Vertex w);
    void remove_edge(Edge edge);
    void remove_vertex(Vertex v);
    std::set<Vertex> neighbors(Vertex v);
    // check returns true if a vertex is in an exclusion constraint
    bool fixed(Vertex v) { return exclusion_classes_[exclude_[v]].size() > 1; }

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
    Vertex preprocess_three_separator();
    Vertex limit_max_length();

    // preprocessing stats
    Vertex isolated_removed = 0;
    Vertex forwarder_removed = 0;
    Vertex position_determined_removed = 0;
    Vertex twin_edges_removed = 0;
    Vertex unreachable_removed = 0;
    Vertex unusable_edge_removed = 0;
    Vertex two_sep_removed = 0;
    Vertex three_sep_removed = 0;
    Vertex max_length_decrease = 0;


};

} // namespace fpc
