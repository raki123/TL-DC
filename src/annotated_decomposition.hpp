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
#include "logging.h"
#include "types.h"
#include <optional>

namespace fpc {

enum NodeType { LEAF, PATH_LIKE, JOIN };

struct AnnotatedNode {
  size_t parent; // parent of the node or -1 if root
  NodeType type;
  Edge edge;                          // if LEAF or PATH_LIKE
  std::pair<size_t, size_t> children; // if JOIN the indices of the child nodes
                                      // in the annotated decomposition
  std::vector<Vertex> bag;

  void stats() const;
};

class TreeDecomposition;
class BranchDecomposition;
class Graph;
class Digraph;

// contains annotated nodes
// path like nodes continuously from leaf to last node before join
// join node occurs after both children nodes
class AnnotatedDecomposition {
public:
  std::size_t size() const { return nr_bags_; }
  std::size_t width() const { return width_; }
  std::size_t join_width() const { return max_join_; }

  Vertex get_root() const { return root_; }

  static AnnotatedDecomposition invalid() {
    return AnnotatedDecomposition(0, -1);
  }

  static AnnotatedDecomposition
  from_vertex_order(std::vector<std::size_t> const &order, Graph const &graph);
  static AnnotatedDecomposition
  from_vertex_order(std::vector<std::size_t> const &order,
                    Digraph const &graph);

  static AnnotatedDecomposition
  from_tree_decomposition(TreeDecomposition const &tree_decomposition,
                          Graph const &graph);
  static AnnotatedDecomposition
  from_tree_decomposition(TreeDecomposition const &tree_decomposition,
                          Digraph const &graph);

  static AnnotatedDecomposition
  from_branch_decomposition(BranchDecomposition const &branch_decomposition,
                            Graph const &graph);
  static AnnotatedDecomposition
  from_branch_decomposition(BranchDecomposition const &branch_decomposition,
                            Digraph const &graph);

  template <typename fun>
  void
  foreach_post_order(fun const &handle,
                     std::optional<Vertex> at_or_below = std::nullopt) const {
    // {neighbor index, bag index}
    std::vector<std::pair<Vertex, Vertex>> stack;
    stack.emplace_back(0, at_or_below.value_or(root_));
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
  template <typename fun>
  void
  foreach_pre_order(fun const &handle,
                    std::optional<Vertex> at_or_below = std::nullopt) const {
    // {neighbor index, bag index}
    std::vector<std::pair<Vertex, Vertex>> stack;
    stack.emplace_back(0, at_or_below.value_or(root_));
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

  AnnotatedNode const &operator[](std::size_t idx) const {
    return node_label_[idx];
  }
  AnnotatedNode &operator[](std::size_t idx) { return node_label_[idx]; }

  std::pair<Vertex, Vertex> children(Vertex idx) const {
    assert(root_ != -1);
    assert(children_[idx].size() <= 2);
    if (children_[idx].empty()) {
      return {-1, -1};
    } else if (children_[idx].size() == 1) {
      return {children_[idx].front(), -1};
    }
    return {children_[idx].front(), children_[idx].back()};
  }
  Vertex parent(Vertex idx) const {
    assert(root_ != -1);
    return parent_[idx];
  }

  void print_stats() const {
    LOG << "max bag size " << width_ << ", max join child bags " << nr_max_join_
        << ", max join bag size " << max_join_ << ", nr of joins " << nr_joins_
        << std::endl;
  }

private:
  AnnotatedDecomposition(std::size_t nr_bags, std::size_t width)
      : nr_bags_(nr_bags), width_(width), nr_max_join_(0), max_join_(0),
        nr_joins_(0), root_(-1), neighbors_(nr_bags_), children_(nr_bags_),
        parent_(nr_bags_, -1), node_label_(nr_bags_) {}

  AnnotatedDecomposition(std::vector<AnnotatedNode> nodes);

  void set_root(Vertex root);

  void set_label(std::size_t bag, AnnotatedNode label) {
    assert(node_label_[bag].bag.empty());
    node_label_[bag] = label;
  }
  void add_edge(Vertex v, Vertex w) {
    neighbors_[v].push_back(w);
    neighbors_[w].push_back(v);
  }

  friend std::ostream &operator<<(std::ostream &out,
                                  const AnnotatedDecomposition &ad);

  std::size_t nr_bags_;
  std::size_t width_;
  std::size_t nr_max_join_;
  std::size_t max_join_;
  std::size_t nr_joins_;
  Vertex root_;
  std::vector<std::vector<Vertex>> neighbors_;
  std::vector<std::vector<Vertex>> children_;
  std::vector<Vertex> parent_;
  std::vector<AnnotatedNode> node_label_;
};

} // namespace fpc
