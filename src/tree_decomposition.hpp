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
#include "assert.h"
#include "types.h"
#include <ranges>
#include <vector>

namespace fpc {

class Graph;
class Digraph;

class TreeDecomposition {
public:
  enum class Strategy { HTD, FLOW_CUTTER, TAMAKI };
  std::size_t nr_bags() const { return nr_bags_; }
  std::size_t width() const { return width_; }

  static TreeDecomposition from_graph(Graph const &graph, Strategy strategy);
  static TreeDecomposition from_graph(Digraph const &graph, Strategy strategy);

  void set_root(std::size_t root);
  std::size_t get_root() const { return root_; }

  template <typename fun> void foreach_post_order(fun const &handle) const {
    // {neighbor index, bag index}
    std::vector<std::pair<std::size_t, std::size_t>> stack;
    stack.emplace_back(0, root_);
    while (!stack.empty()) {
      auto [neigh_idx, bag_idx] = stack.back();
      if (neigh_idx < children_[bag_idx].size()) {
        // go down
        stack.back().first++;
        stack.emplace_back(0, children_[bag_idx][neigh_idx]);
      } else {
        // done with everything below, handle this
        stack.pop_back();
        handle(bag_idx);
      }
    }
  }
  template <typename fun> void foreach_pre_order(fun const &handle) const {
    // {neighbor index, bag index}
    std::vector<std::pair<std::size_t, std::size_t>> stack;
    stack.emplace_back(0, root_);
    handle(root_);
    while (!stack.empty()) {
      auto [neigh_idx, bag_idx] = stack.back();
      if (neigh_idx < children_[bag_idx].size()) {
        // go down
        stack.back().first++;
        stack.emplace_back(0, children_[bag_idx][neigh_idx]);
        // first handle current
        handle(children_[bag_idx][neigh_idx]);
      } else {
        // go up until we find something new
        stack.pop_back();
      }
    }
  }

  std::vector<Vertex> const &operator[](std::size_t idx) const {
    return node_label_[idx];
  }
  std::vector<std::size_t> const &children(std::size_t idx) const {
    assert(root_ != -1);
    return children_[idx];
  }
  std::size_t parent(std::size_t idx) const {
    assert(root_ != -1);
    return parent_[idx];
  }
  void enforce_child_limit(std::size_t max_children);

private:
  struct Pipes {
    int in;
    int out;
  };
  static Pipes open_htd(int iterations);
  static Pipes open_flow_cutter();
  static Pipes open_tamaki();
  static void write_instance(Graph const &graph, int in);
  static void write_instance(Digraph const &graph, int in);
  static TreeDecomposition from_file(int out);

  TreeDecomposition(std::size_t nr_bags, std::size_t width)
      : nr_bags_(nr_bags), width_(width), root_(-1), neighbors_(nr_bags_),
        children_(nr_bags_), parent_(nr_bags_, -1), node_label_(nr_bags_) {}

  void set_label(std::size_t bag, std::vector<Vertex> label) {
    assert(node_label_[bag].empty());
    node_label_[bag] = label;
  }
  void add_edge(std::size_t v, std::size_t w) {
    neighbors_[v].push_back(w);
    neighbors_[w].push_back(v);
  }
  bool is_leaf(std::size_t bag) const { return neighbors_[bag].size() == 1; }

  friend std::ostream &operator<<(std::ostream &out,
                                  const TreeDecomposition &td);

  std::size_t nr_bags_;
  std::size_t width_;
  std::size_t root_;
  std::vector<std::vector<std::size_t>> neighbors_;
  std::vector<std::vector<std::size_t>> children_;
  std::vector<std::size_t> parent_;
  std::vector<std::vector<Vertex>> node_label_;
};

std::ostream &operator<<(std::ostream &out, const TreeDecomposition &td);
} // namespace fpc