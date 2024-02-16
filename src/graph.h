// TLDC "Too long; Didn't Count" A length limited path counter.
// Copyright (C) 2023-2024 Rafael Kiesel, Markus Hecher

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
#include <ranges>
#include <set>
#include <unordered_set>
#include <vector>

#include "count_structures.h"
#include "types.h"

namespace srg = std::ranges;

namespace fpc {

class Graph {

public:
  friend class Search;
  friend class ParallelNautySearch;
  Graph(std::istream &input);
  Graph(Vertex n);

  void preprocess();
  void weighted_preprocess();
  void reduce_modulo_equivalence();

  void print_stats() const;

  void print_ordered(std::ostream &os,
                     std::vector<std::size_t> const &order) const;

  Edge_weight const &extra_paths() const { return extra_paths_; }

  void normalize(bool reorder = false);

  bool is_all_pair() const { return all_pair_; }

  bool is_length_limited() const;

  Edge_length max_length() const { return max_length_; }
  Edge_length min_length() const;

  // how many digits the number of paths needs at most in base 2.
  // I.e., #paths < 2^(bit_upper_bound)
  std::size_t bit_upper_bound() const;

  std::vector<Vertex> const &terminals() const { return terminals_; }

  size_t nr_vertices() const {
    assert(nr_vertices_ == count_vertices());
    return nr_vertices_;
  };
  size_t nr_edges() const {
    assert(nr_edges_ == count_edges());
    return nr_edges_;
  };

  std::vector<Vertex> const &neighbors(Vertex v) const {
    assert(v < neighbors_.size());
    return neighbors_[v];
  };

  void add_edge(Vertex v, Vertex w) {
    assert(v != w);
    nr_edges_++;
    assert(!has_edge(v, w));
    if (neighbors_[v].empty()) {
      nr_vertices_++;
    }
    if (neighbors_[w].empty()) {
      nr_vertices_++;
    }
    adjacency_[v][w] = true;
    adjacency_[w][v] = true;
    neighbors_[v].push_back(w);
    neighbors_[w].push_back(v);
  }
  bool has_edge(Vertex v, Vertex w) const { return adjacency_[v][w]; }
  void remove_edge(Vertex v, Vertex w) {
    assert(has_edge(v, w));
    nr_edges_--;
    adjacency_[v][w] = false;
    adjacency_[w][v] = false;
    auto it = srg::find(neighbors_[v], w);
    assert(it != neighbors_[v].end());
    neighbors_[v].erase(it);
    it = srg::find(neighbors_[w], v);
    assert(it != neighbors_[w].end());
    neighbors_[w].erase(it);
    if (neighbors_[v].empty()) {
      nr_vertices_--;
    }
    if (neighbors_[w].empty()) {
      nr_vertices_--;
    }
    if (v > w) {
      std::swap(v, w);
    }
    edge_weights_.erase(Edge(v, w));
  };

  Edge_length edge_length(Vertex v, Vertex w) const {
    if (auto it = find_weight(v, w); it != edge_weights_.end()) {
      return it->second.offset();
    }
    return 1;
  }

  std::optional<Limited_count<mpz_class>> edge_weight(Vertex v, Vertex w) {
    if (auto it = find_weight(v, w); it != edge_weights_.end()) {
      return it->second;
    }
    return {};
  }

  void remove_vertex(Vertex v) {
    auto neighs = neighbors(v);
    for (auto neigh : neighs) {
      assert(neigh != v);
      remove_edge(v, neigh);
    }
    assert(neighbors_[v].empty());
  };

  Graph subgraph(std::vector<Vertex> restrict_to) const;
  Graph copy() const;

  sparsegraph to_canon_nauty(bool reorder);
  std::vector<sparsegraph> all_pair_nauty(bool reorder);
  double nr_automorphisms();

  // get vector of vectors containing the non-trivial equivalence classes w.r.t.
  // neighborhood diversity
  std::vector<std::vector<Vertex>> equivalence_classes() const;

  // utility functions
  void dijkstra(Vertex start, std::vector<Edge_length> &distance,
                const std::set<Vertex> &forbidden) const;

private:
  size_t count_vertices() const {
    size_t rv = 0;
    for (Vertex v = 0; v < adjacency_.size(); v++) {
      if (!neighbors(v).empty()) {
        rv++;
      }
    }
    return rv;
  };
  size_t count_edges() const {
    size_t rv = 0;
    for (Vertex v = 0; v < adjacency_.size(); v++) {
      rv += neighbors(v).size();
    }
    assert(rv % 2 == 0);
    return rv / 2;
  };

  std::unordered_map<Edge, Limited_count<mpz_class>>::iterator
  find_weight(Vertex v, Vertex w) {
    if (v > w) {
      std::swap(v, w);
    }
    return edge_weights_.find(Edge(v, w));
  }

  std::unordered_map<Edge, Limited_count<mpz_class>>::const_iterator
  find_weight(Vertex v, Vertex w) const {
    if (v > w) {
      std::swap(v, w);
    }
    return edge_weights_.find(Edge(v, w));
  }

  std::unordered_set<Vertex>
  must_use_subgraph(Vertex must_use,
                    std::vector<Edge_length> const &distance_from_start,
                    std::vector<Edge_length> const &distance_to_goal,
                    Edge_length back_ward_allowance) const;

  size_t nr_vertices_;
  size_t nr_edges_;
  Edge_length max_length_;
  std::vector<std::vector<Vertex>> neighbors_;
  std::vector<std::vector<char>> adjacency_;
  std::vector<Vertex> terminals_;
  std::unordered_map<Edge, Limited_count<mpz_class>> edge_weights_;

  Edge_weight extra_paths_;
  std::optional<bool> is_length_limited_ = std::nullopt;
  std::optional<std::size_t> bit_upper_bound_ = std::nullopt;

  bool all_pair_ = false;

  void contract_edges(Vertex v);
  void contract_terminal_edges(Vertex v, Vertex terminal);
  // preprocessing subroutines
  Vertex preprocess_start_goal_edges();
  Vertex preprocess_start_goal_forwarders();
  Vertex preprocess_isolated();
  Vertex preprocess_forwarder();
  Vertex preprocess_unreachable();
  Vertex preprocess_unusable_edge();
  Vertex preprocess_articulation_points();
  std::pair<bool, Vertex> ap_util(Vertex u, std::vector<char> &visited,
                                  std::vector<Vertex> &disc,
                                  std::vector<Vertex> &low, int &time,
                                  int parent);
  Vertex prune_util(Vertex ap, Vertex prune);
  Vertex preprocess_position_determined(Edge_length budget);

  // weight inducing preprocessing
  Vertex weighted_preprocess_forwarder();

  // preprocessing stats
  Vertex isolated_removed = 0;
  Vertex forwarder_removed = 0;
  Vertex position_determined_removed = 0;
  Vertex articulation_point_removed = 0;
  Vertex unreachable_removed = 0;
  Vertex unusable_edge_removed = 0;
  Vertex max_length_decrease = 0;
};

} // namespace fpc
