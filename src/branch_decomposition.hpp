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
#include "annotated_decomposition.hpp"
#include "graph.h"
#include <ranges>
#include <vector>

namespace fpc {
class BranchDecomposition {
public:
  enum class Strategy { HICKS };
  std::size_t nr_edges() const { return nr_edges_; }
  std::size_t nr_bags() const { return nr_bags_; }
  std::size_t width() const { return width_; }

  static BranchDecomposition from_graph(Graph const &graph, Strategy strategy);

  template <typename fun> void foreach_post_order(fun const &handle) const {
    std::vector<Vertex> stack;
    Vertex cur = root_;
    bool down = true;
    while (cur != root_ || down) {
      if (down) {
        if (is_leaf(cur)) {
          down = false;
          handle(cur);
        } else {
          stack.push_back(cur);
          cur = left_child(cur);
        }
      } else {
        auto last = cur;
        cur = stack.back();
        stack.pop_back();
        if (left_child(cur) == last) {
          stack.push_back(cur);
          cur = right_child(cur);
          down = true;
        } else {
          handle(cur);
        }
      }
    }
  }

  bool is_leaf(Vertex bag) const {
    assert(children_[bag].empty() || leaf_label_[bag] == Edge(-1, -1));
    return children_[bag].empty();
  }
  Vertex left_child(Vertex bag) const {
    assert(children_[bag].size() == 2);
    return children_[bag].front();
  }
  Vertex right_child(Vertex bag) const {
    assert(children_[bag].size() == 2);
    return children_[bag].back();
  }
  Vertex parent(Vertex bag) const { return parent_[bag]; }

  Edge const &label(Vertex bag) const { return leaf_label_[bag]; }

private:
  static BranchDecomposition hicks(Graph const &graph);

  BranchDecomposition(std::size_t nr_bags, std::size_t nr_edges,
                      std::size_t width)
      : nr_bags_(nr_bags), nr_edges_(nr_edges), width_(width), root_(-1),
        children_(nr_bags_), parent_(nr_bags_, std::size_t(-1)),
        leaf_label_(nr_bags_, Edge(-1, -1)), edge_count_(nr_edges_) {}

  void set_label(std::size_t bag, Edge label) {
    assert(leaf_label_[bag] == Edge(-1, -1));
    assert(label.first < nr_edges_);
    assert(label.second < nr_edges_);
    leaf_label_[bag] = label;
    edge_count_[label.first]++;
    edge_count_[label.second]++;
  }
  void set_root(Vertex root) {
    assert(root_ == Vertex(-1));
    root_ = root;
  }
  void add_arc(Vertex parent, Vertex child) {
    assert(root_ != child);
    assert(parent_[child] == std::size_t(-1));
    children_[parent].push_back(child);
    parent_[child] = parent;
  }

  std::size_t nr_bags_;
  std::size_t nr_edges_;
  std::size_t width_;
  Vertex root_;
  std::vector<std::vector<Vertex>> children_;
  std::vector<std::size_t> parent_;
  std::vector<Edge> leaf_label_;
  std::vector<std::size_t> edge_count_;
};
} // namespace fpc