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
#include <assert.h>

#include <iostream>
#include <optional>
#include <ranges>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "count_structures.h"
#include "types.h"

namespace srg = std::ranges;

namespace fpc {

class Digraph {

public:
  Digraph(std::istream &input);
  Digraph(Vertex n);

  void preprocess();
  void weighted_preprocess();
  void reduce_modulo_equivalence();

  void print_stats() const;

  void print_ordered(std::ostream &os,
                     std::vector<std::size_t> const &order) const;

  Edge_weight const &extra_paths() const { return extra_paths_; }

  void normalize(bool reorder = false);

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
  size_t nr_arcs() const {
    assert(nr_arcs_ == count_arcs());
    return nr_arcs_;
  };

  std::vector<Vertex> const &in_neighbors(Vertex v) const {
    assert(v < in_neighbors_.size());
    return in_neighbors_[v];
  };

  std::vector<Vertex> const &out_neighbors(Vertex v) const {
    assert(v < out_neighbors_.size());
    return out_neighbors_[v];
  };

  bool has_arc(Vertex from, Vertex to) const {
    auto it = srg::find(in_neighbors_[to], from);
    return it != in_neighbors_[to].end();
  }
  void remove_arc(Vertex from, Vertex to) {
    assert(has_arc(from, to));
    nr_arcs_--;
    auto it = srg::find(in_neighbors_[to], from);
    assert(it != in_neighbors_[to].end());
    std::swap(*it, in_neighbors_[to].back());
    in_neighbors_[to].pop_back();
    it = srg::find(out_neighbors_[from], to);
    assert(it != out_neighbors_[from].end());
    std::swap(*it, out_neighbors_[from].back());
    out_neighbors_[from].pop_back();
    if (in_neighbors_[to].empty() && out_neighbors_[to].empty()) {
      nr_vertices_--;
    }
    if (in_neighbors_[from].empty() && out_neighbors_[from].empty()) {
      nr_vertices_--;
    }
    arc_weights_.erase(Edge(from, to));
  }

  Edge_length arc_length(Vertex from, Vertex to) const {
    if (auto it = arc_weights_.find(Edge(from, to)); it != arc_weights_.end()) {
      return it->second.offset();
    }
    return 1;
  }

  std::optional<Limited_count<mpz_class>> arc_weight(Vertex from, Vertex to) {
    if (auto it = arc_weights_.find(Edge(from, to)); it != arc_weights_.end()) {
      return it->second;
    }
    return {};
  }

  Digraph subgraph(std::vector<Vertex> const &restrict_to) const;
  Digraph copy() const;

  // utility functions
  void dijkstra(Vertex start, std::vector<Edge_length> &distance, bool outgoing,
                std::optional<Vertex> forbidden = std::nullopt) const;
  std::vector<std::vector<Vertex>> strong_connected_components() const;

  // get vector of vectors containing the non-trivial equivalence classes w.r.t.
  // neighborhood diversity
  std::vector<std::vector<Vertex>> equivalence_classes() const;

  mpz_class naive_dfs() const;
  void print_dfvs_instance() const;

private:
  void scc_dfs(Vertex start, std::vector<Vertex> &disc,
               std::vector<Vertex> &low, std::vector<char> &on_stack,
               std::vector<Vertex> &stack, Vertex &time,
               std::vector<std::vector<Vertex>> &sccs) const;

  std::unordered_set<Vertex>
  must_use_subgraph(Vertex must_use,
                    std::vector<Edge_length> const &distance_from_start,
                    std::vector<Edge_length> const &distance_to_goal,
                    Edge_length back_ward_allowance) const;

  void add_arc(Vertex from, Vertex to) {
    assert(from != to);
    nr_arcs_++;
    if (in_neighbors_[to].empty() && out_neighbors_[to].empty()) {
      nr_vertices_++;
    }
    if (in_neighbors_[from].empty() && out_neighbors_[from].empty()) {
      nr_vertices_++;
    }
    in_neighbors_[to].push_back(from);
    out_neighbors_[from].push_back(to);
  };
  void contract_arcs(Vertex v, bool unsafe = false);
  void remove_vertex(Vertex v) {
    auto in_neighs = in_neighbors(v);
    for (auto neigh : in_neighs) {
      assert(neigh != v);
      remove_arc(neigh, v);
    }
    assert(in_neighbors_[v].empty());

    auto out_neighs = out_neighbors(v);
    for (auto neigh : out_neighs) {
      assert(neigh != v);
      remove_arc(v, neigh);
    }
    assert(out_neighbors_[v].empty());
  }

  size_t count_vertices() const {
    size_t rv = 0;
    for (Vertex v = 0; v < in_neighbors_.size(); v++) {
      if (!in_neighbors(v).empty() || !out_neighbors(v).empty()) {
        rv++;
      }
    }
    return rv;
  };
  size_t count_arcs() const {
    size_t rv = 0;
    for (Vertex v = 0; v < in_neighbors_.size(); v++) {
      rv += out_neighbors(v).size();
    }
    return rv;
  };

  size_t nr_vertices_;
  size_t nr_arcs_;
  Edge_length max_length_;
  std::vector<std::vector<Vertex>> in_neighbors_;
  std::vector<std::vector<Vertex>> out_neighbors_;
  std::vector<Vertex> terminals_;
  std::unordered_map<Edge, Limited_count<mpz_class>> arc_weights_;

  Edge_weight extra_paths_;
  std::optional<bool> is_length_limited_ = std::nullopt;
  std::optional<std::size_t> bit_upper_bound_ = std::nullopt;

  // preprocessing subroutines
  Vertex preprocess_start_goal_edges();
  Vertex preprocess_isolated();
  Vertex preprocess_unreachable();
  Vertex preprocess_forwarder();
  Vertex preprocess_dominator_arcs();
  Vertex preprocess_unusable_arcs();
  Vertex preprocess_position_determined(Edge_length budget);

  // weight inducing preprocessing
  Vertex weighted_preprocess_forwarder();

  // preprocessing stats
  Vertex isolated_removed = 0;
  Vertex forwarder_removed = 0;
  Vertex position_determined_removed = 0;
  Vertex unreachable_removed = 0;
  Vertex dominator_arcs_removed = 0;
  Vertex unusable_arcs_removed = 0;
  Vertex max_length_decrease = 0;
};

} // namespace fpc
