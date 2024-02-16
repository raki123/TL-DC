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

#include "annotated_decomposition.hpp"
#include "branch_decomposition.hpp"
#include "digraph.h"
#include "graph.h"
#include "tree_decomposition.hpp"
#include <map>

namespace fpc {

void AnnotatedNode::stats() const {
  LOG << "TYPE " << type << ", parent: " << parent << ", edge: (" << edge.first
      << "," << edge.second << ") children: (" << children.first << ","
      << children.second << ") bag {";
  for (auto it = bag.begin(); it != bag.end(); ++it) {
    LOG << *it << ",";
  }
  LOG << "}" << std::endl;
}

AnnotatedDecomposition::AnnotatedDecomposition(std::vector<AnnotatedNode> nodes)
    : nr_bags_(nodes.size()), width_(0), nr_max_join_(0), max_join_(0),
      nr_joins_(0), root_(-1), neighbors_(nr_bags_), children_(nr_bags_),
      parent_(nr_bags_, -1), node_label_(nr_bags_) {

  std::map<size_t, size_t> idx_remap;
  size_t invalid_td_idx = -1;
  idx_remap[invalid_td_idx] = invalid_td_idx;
  size_t cur = nodes.size() - 1;

  size_t root = 0;
  while (root < nodes.size() && nodes[root].parent != invalid_td_idx) {
    root++;
  }

  std::vector<size_t> remap_stack = {root};
  std::vector<AnnotatedNode> ordered;
  while (!remap_stack.empty()) {
    size_t top = remap_stack.back();
    idx_remap[top] = cur--;
    remap_stack.pop_back();
    auto &top_node = nodes[top];
    switch (top_node.type) {
    case NodeType::LEAF:
      ordered.push_back(top_node);
      break;
    case NodeType::PATH_LIKE:
      assert(top > 0);
      ordered.push_back(top_node);
      remap_stack.push_back(top_node.children.first);
      break;
    case NodeType::JOIN:
      ordered.push_back(top_node);
      remap_stack.push_back(top_node.children.first);
      remap_stack.push_back(top_node.children.second);
      break;
    }
  }
  nodes = std::vector<AnnotatedNode>(ordered.rbegin(), ordered.rend());
  for (auto &node : nodes) {
    node.parent = idx_remap[node.parent];
    node.children.first = idx_remap[node.children.first];
    node.children.second = idx_remap[node.children.second];
  }

  for (std::size_t i = 0; i < nodes.size(); i++) {
    auto const &node = nodes[i];
    if (node.children.first != -1) {
      add_edge(i, node.children.first);
    }
    if (node.children.second != -1) {
      assert(node.type == JOIN);
      add_edge(i, node.children.second);
    }
    width_ = std::max(width_, node.bag.size());
    if (node.type == JOIN) {
      nr_joins_++;
      auto join_size = nodes[node.children.first].bag.size() +
                       nodes[node.children.second].bag.size();
      if (join_size > max_join_) {
        max_join_ = join_size;
        nr_max_join_ = 1;
      } else if (join_size == max_join_) {
        nr_max_join_++;
      }
    }
    set_label(i, node);
  }

  set_root(root);
}

void AnnotatedDecomposition::set_root(Vertex root) {
  std::vector<std::tuple<Vertex, Vertex>> stack;
  root_ = root;
  stack.emplace_back(-1, root);
  while (!stack.empty()) {
    auto [parent, cur] = stack.back();
    stack.pop_back();
    children_[cur].clear();
    parent_[cur] = parent;
    for (auto neighbor : neighbors_[cur]) {
      if (neighbor == parent) {
        continue;
      }
      children_[cur].push_back(neighbor);
      stack.emplace_back(cur, neighbor);
    }

    node_label_[cur].children = children(cur);
    node_label_[cur].parent = parent_[cur];
  }
}

AnnotatedDecomposition
AnnotatedDecomposition::from_vertex_order(std::vector<std::size_t> const &order,
                                          Graph const &graph) {
  std::vector<AnnotatedNode> nodes;
  std::size_t width = 0;
  std::set<Vertex> active;
  std::vector<char> taken(graph.nr_vertices(), false);
  for (auto v : order) {
    active.insert(v);
    for (auto neigh : graph.neighbors(v)) {
      if (taken[neigh]) {
        continue;
      }
      active.insert(neigh);
      width = std::max(width, active.size());
      AnnotatedNode node;
      node.parent = nodes.size() + 1;
      node.edge = Edge(v, neigh);
      node.bag = {active.begin(), active.end()};
      node.children.first = nodes.size() - 1;
      node.children.second = -1;
      node.type = (nodes.empty() ? LEAF : PATH_LIKE);
      nodes.push_back(node);
    }
    taken[v] = true;
    active.erase(v);
  }
  nodes.back().parent = std::size_t(-1);

  AnnotatedDecomposition rv(std::move(nodes));

  return rv;
}

AnnotatedDecomposition
AnnotatedDecomposition::from_vertex_order(std::vector<std::size_t> const &order,
                                          Digraph const &graph) {
  std::vector<AnnotatedNode> nodes;
  std::size_t width = 0;
  std::set<Vertex> active;
  std::vector<char> taken(graph.nr_vertices(), false);
  for (auto v : order) {
    active.insert(v);
    for (auto neigh : graph.in_neighbors(v)) {
      if (taken[neigh]) {
        continue;
      }
      active.insert(neigh);
      width = std::max(width, active.size());
      AnnotatedNode node;
      node.parent = nodes.size() + 1;
      node.edge = Edge(neigh, v);
      node.bag = {active.begin(), active.end()};
      node.children.first = nodes.size() - 1;
      node.children.second = -1;
      node.type = (nodes.empty() ? LEAF : PATH_LIKE);
      nodes.push_back(node);
    }
    for (auto neigh : graph.out_neighbors(v)) {
      if (taken[neigh]) {
        continue;
      }
      active.insert(neigh);
      width = std::max(width, active.size());
      AnnotatedNode node;
      node.parent = nodes.size() + 1;
      node.edge = Edge(v, neigh);
      node.bag = {active.begin(), active.end()};
      node.children.first = nodes.size() - 1;
      node.children.second = -1;
      node.type = (nodes.empty() ? LEAF : PATH_LIKE);
      nodes.push_back(node);
    }
    taken[v] = true;
    active.erase(v);
  }
  nodes.back().parent = std::size_t(-1);

  AnnotatedDecomposition rv(std::move(nodes));

  return rv;
}

AnnotatedDecomposition AnnotatedDecomposition::from_tree_decomposition(
    TreeDecomposition const &tree_decomposition, Graph const &graph) {
  std::vector<AnnotatedNode> nodes;
  // vertices are eliminated in the bag in which they last occur
  std::vector<Vertex> last_occurrence(graph.nr_vertices(),
                                      tree_decomposition.get_root());
  tree_decomposition.foreach_post_order([&](auto bag_idx) {
    for (auto v : tree_decomposition[bag_idx]) {
      last_occurrence[v] = bag_idx;
    }
  });

  std::vector<char> taken(graph.nr_vertices(), false);
  std::vector<std::size_t> bag_to_last_annotated(tree_decomposition.nr_bags(),
                                                 -1);
  std::unordered_set<std::size_t> used_as_child;
  tree_decomposition.foreach_post_order([&](auto bag_idx) {
    // determine the type
    NodeType type = [&]() {
      auto nr_children = srg::count_if(
          tree_decomposition.children(bag_idx), [&](auto child_idx) {
            return bag_to_last_annotated[child_idx] != -1;
          });
      if (nr_children == 0) {
        return LEAF;
      }
      if (nr_children == 1) {
        return PATH_LIKE;
      }
      assert(nr_children == 2);
      return JOIN;
    }();

    // for joins we first add the join node before we can eliminate
    // vertices
    if (type == JOIN) {
      auto first_idx = tree_decomposition.children(bag_idx).front();
      auto first_annotated = bag_to_last_annotated[first_idx];
      auto second_idx = tree_decomposition.children(bag_idx).back();
      auto second_annotated = bag_to_last_annotated[second_idx];
      assert(first_annotated != -1);
      assert(second_annotated != -1);
      assert(first_annotated != second_annotated);
      // otherwise we did not eliminate anything in one of the branches
      AnnotatedNode node;
      node.parent = nodes.size() + 1;
      node.edge = Edge(-1, -1);
      node.bag = tree_decomposition[bag_idx];
      node.children.first = first_annotated;
      auto [_, inserted] = used_as_child.insert(first_annotated);
      assert(inserted);
      nodes[first_annotated].parent = nodes.size();
      node.children.second = second_annotated;
      auto [_2, inserted2] = used_as_child.insert(second_annotated);
      assert(inserted2);
      nodes[second_annotated].parent = nodes.size();
      node.type = type;
      nodes.push_back(node);
      bag_to_last_annotated[bag_idx] = nodes.size() - 1;

      // afterwards we may elminate edges in a path like manner
      type = PATH_LIKE;
    }
    auto const &bag = tree_decomposition[bag_idx];
    for (auto v : bag) {
      if (last_occurrence[v] != bag_idx) {
        // v is eliminated later
        continue;
      }

      for (auto neigh : graph.neighbors(v)) {
        if (taken[neigh]) {
          continue;
        }
        AnnotatedNode node;
        node.parent = nodes.size() + 1;
        node.edge = Edge(v, neigh);
        node.bag = bag;
        node.children.first = type == LEAF ? -1 : nodes.size() - 1;
        if (type != LEAF) {
          auto [_, inserted] = used_as_child.insert(nodes.size() - 1);
          assert(inserted);
        }
        node.children.second = -1;
        node.type = type;
        nodes.push_back(node);
        // only the first edge that is eliminated is at a leaf node
        type = PATH_LIKE;
        bag_to_last_annotated[bag_idx] = nodes.size() - 1;
      }
      taken[v] = true;
    }
    if (bag_to_last_annotated[bag_idx] == -1) {
      std::size_t first_annotated = -1;
      std::size_t second_annotated = -1;
      if (tree_decomposition.children(bag_idx).size() > 0) {
        auto first_idx = tree_decomposition.children(bag_idx).front();
        first_annotated = bag_to_last_annotated[first_idx];
      }
      if (tree_decomposition.children(bag_idx).size() > 1) {
        auto second_idx = tree_decomposition.children(bag_idx).back();
        second_annotated = bag_to_last_annotated[second_idx];
      }
      assert(first_annotated == -1 || second_annotated == -1);
      if (first_annotated != -1) {
        bag_to_last_annotated[bag_idx] = first_annotated;
      }
      if (second_annotated != -1) {
        bag_to_last_annotated[bag_idx] = second_annotated;
      }
    }
  });

  for (auto t : taken) {
    assert(t);
  }

  nodes.back().parent = std::size_t(-1);

  AnnotatedDecomposition rv(std::move(nodes));
  return rv;
}

AnnotatedDecomposition AnnotatedDecomposition::from_tree_decomposition(
    TreeDecomposition const &tree_decomposition, Digraph const &graph) {
  std::vector<AnnotatedNode> nodes;
  // vertices are eliminated in the bag in which they last occur
  std::vector<Vertex> last_occurrence(graph.nr_vertices(),
                                      tree_decomposition.get_root());
  tree_decomposition.foreach_post_order([&](auto bag_idx) {
    for (auto v : tree_decomposition[bag_idx]) {
      last_occurrence[v] = bag_idx;
    }
  });

  std::vector<char> taken(graph.nr_vertices(), false);
  std::vector<std::size_t> bag_to_last_annotated(tree_decomposition.nr_bags(),
                                                 -1);
  std::unordered_set<std::size_t> used_as_child;
  tree_decomposition.foreach_post_order([&](auto bag_idx) {
    // determine the type
    NodeType type = [&]() {
      auto nr_children = srg::count_if(
          tree_decomposition.children(bag_idx), [&](auto child_idx) {
            return bag_to_last_annotated[child_idx] != -1;
          });
      if (nr_children == 0) {
        return LEAF;
      }
      if (nr_children == 1) {
        return PATH_LIKE;
      }
      assert(nr_children == 2);
      return JOIN;
    }();

    // for joins we first add the join node before we can eliminate
    // vertices
    if (type == JOIN) {
      auto first_idx = tree_decomposition.children(bag_idx).front();
      auto first_annotated = bag_to_last_annotated[first_idx];
      auto second_idx = tree_decomposition.children(bag_idx).back();
      auto second_annotated = bag_to_last_annotated[second_idx];
      assert(first_annotated != -1);
      assert(second_annotated != -1);
      assert(first_annotated != second_annotated);
      // otherwise we did not eliminate anything in one of the branches
      AnnotatedNode node;
      node.parent = nodes.size() + 1;
      node.edge = Edge(-1, -1);
      node.bag = tree_decomposition[bag_idx];
      node.children.first = first_annotated;
      auto [_, inserted] = used_as_child.insert(first_annotated);
      assert(inserted);
      nodes[first_annotated].parent = nodes.size();
      node.children.second = second_annotated;
      auto [_2, inserted2] = used_as_child.insert(second_annotated);
      assert(inserted2);
      nodes[second_annotated].parent = nodes.size();
      node.type = type;
      nodes.push_back(node);
      bag_to_last_annotated[bag_idx] = nodes.size() - 1;

      // afterwards we may elminate edges in a path like manner
      type = PATH_LIKE;
    }
    auto const &bag = tree_decomposition[bag_idx];
    for (auto v : bag) {
      if (last_occurrence[v] != bag_idx) {
        // v is eliminated later
        continue;
      }

      for (auto neigh : graph.in_neighbors(v)) {
        if (taken[neigh]) {
          continue;
        }
        AnnotatedNode node;
        node.parent = nodes.size() + 1;
        node.edge = Edge(neigh, v);
        node.bag = bag;
        node.children.first = type == LEAF ? -1 : nodes.size() - 1;
        if (type != LEAF) {
          auto [_, inserted] = used_as_child.insert(nodes.size() - 1);
          assert(inserted);
        }
        node.children.second = -1;
        node.type = type;
        nodes.push_back(node);
        // only the first edge that is eliminated is at a leaf node
        type = PATH_LIKE;
        bag_to_last_annotated[bag_idx] = nodes.size() - 1;
      }
      for (auto neigh : graph.out_neighbors(v)) {
        if (taken[neigh]) {
          continue;
        }
        AnnotatedNode node;
        node.parent = nodes.size() + 1;
        node.edge = Edge(v, neigh);
        node.bag = bag;
        node.children.first = type == LEAF ? -1 : nodes.size() - 1;
        if (type != LEAF) {
          auto [_, inserted] = used_as_child.insert(nodes.size() - 1);
          assert(inserted);
        }
        node.children.second = -1;
        node.type = type;
        nodes.push_back(node);
        // only the first edge that is eliminated is at a leaf node
        type = PATH_LIKE;
        bag_to_last_annotated[bag_idx] = nodes.size() - 1;
      }
      taken[v] = true;
    }
    if (bag_to_last_annotated[bag_idx] == -1) {
      std::size_t first_annotated = -1;
      std::size_t second_annotated = -1;
      if (tree_decomposition.children(bag_idx).size() > 0) {
        auto first_idx = tree_decomposition.children(bag_idx).front();
        first_annotated = bag_to_last_annotated[first_idx];
      }
      if (tree_decomposition.children(bag_idx).size() > 1) {
        auto second_idx = tree_decomposition.children(bag_idx).back();
        second_annotated = bag_to_last_annotated[second_idx];
      }
      assert(first_annotated == -1 || second_annotated == -1);
      if (first_annotated != -1) {
        bag_to_last_annotated[bag_idx] = first_annotated;
      }
      if (second_annotated != -1) {
        bag_to_last_annotated[bag_idx] = second_annotated;
      }
    }
  });

  for (auto t : taken) {
    assert(t);
  }

  nodes.back().parent = std::size_t(-1);

  AnnotatedDecomposition rv(std::move(nodes));
  return rv;
}

AnnotatedDecomposition
AnnotatedDecomposition::from_branch_decomposition(BranchDecomposition const &bd,
                                                  Graph const &graph) {

  std::vector<AnnotatedNode> nodes(bd.nr_bags());
  std::vector<std::unordered_map<Vertex, std::size_t>> used_edges(bd.nr_bags());
  auto handle = [&](Vertex bag) {
    if (bd.is_leaf(bag)) {
      auto const &edge = bd.label(bag);
      nodes[bag].edge = edge;
      nodes[bag].children = {-1, -1};
      nodes[bag].type = LEAF;
      used_edges[bag][edge.first]++;
      used_edges[bag][edge.second]++;
      std::erase_if(used_edges[bag], [&](auto const &item) {
        auto const &[v, used] = item;
        return used == graph.neighbors(v).size();
      });
    } else {
      nodes[bag].edge = Edge(-1, -1);
      nodes[bag].children = {bd.left_child(bag), bd.right_child(bag)};
      nodes[bag].type = JOIN;
      used_edges[bag] = used_edges[bd.left_child(bag)];
      for (auto [v, used] : used_edges[bd.right_child(bag)]) {
        used_edges[bag][v] += used;
      }
      std::erase_if(used_edges[bag], [&](auto const &item) {
        auto const &[v, _] = item;
        assert(used_edges[bd.left_child(bag)][v] +
                   used_edges[bd.right_child(bag)][v] <=
               graph.neighbors(v).size());
        return used_edges[bd.left_child(bag)][v] == graph.neighbors(v).size() ||
               used_edges[bd.right_child(bag)][v] == graph.neighbors(v).size();
      });
    }
    nodes[bag].parent = bd.parent(bag);
    std::ranges::transform(used_edges[bag],
                           std::back_insert_iterator(nodes[bag].bag),
                           [](auto const &item) { return item.first; });
  };
  bd.foreach_post_order(handle);
  nodes.back().parent = -1;

  AnnotatedDecomposition rv(nodes);
  return rv;
}

AnnotatedDecomposition
AnnotatedDecomposition::from_branch_decomposition(BranchDecomposition const &bd,
                                                  Digraph const &graph) {
  throw std::runtime_error("Not implemented");

  std::vector<AnnotatedNode> nodes(bd.nr_bags());
  std::vector<std::unordered_map<Vertex, std::size_t>> used_edges(bd.nr_bags());

  auto arc_count = [&](Vertex v) {
    return graph.in_neighbors(v).size() + graph.out_neighbors(v).size();
  };

  auto handle = [&](Vertex bag) {
    if (bd.is_leaf(bag)) {
      auto const &edge = bd.label(bag);
      NodeType type = LEAF;
      if (graph.has_arc(edge.first, edge.second)) {
        nodes[bag].edge = edge;
        nodes[bag].children = {-1, -1};
        nodes[bag].type = type;
        used_edges[bag][edge.first]++;
        used_edges[bag][edge.second]++;
        type = PATH_LIKE;
      }
      auto swapped_edge = Edge(edge.second, edge.first);
      if (graph.has_arc(swapped_edge.first, swapped_edge.second)) {
        nodes[bag].edge = swapped_edge;
        if (type == PATH_LIKE) {
          nodes[bag].children = {-1, -1};
        }
        nodes[bag].type = type;
        used_edges[bag][swapped_edge.first]++;
        used_edges[bag][swapped_edge.second]++;
      }
      std::erase_if(used_edges[bag], [&](auto const &item) {
        auto const &[v, used] = item;
        return used == arc_count(v);
      });
    } else {
      nodes[bag].edge = Edge(-1, -1);
      nodes[bag].children = {bd.left_child(bag), bd.right_child(bag)};
      nodes[bag].type = JOIN;
      used_edges[bag] = used_edges[bd.left_child(bag)];
      for (auto [v, used] : used_edges[bd.right_child(bag)]) {
        used_edges[bag][v] += used;
      }
      std::erase_if(used_edges[bag], [&](auto const &item) {
        auto const &[v, _] = item;
        assert(used_edges[bd.left_child(bag)][v] +
                   used_edges[bd.right_child(bag)][v] <=
               arc_count(v));
        return used_edges[bd.left_child(bag)][v] == arc_count(v) ||
               used_edges[bd.right_child(bag)][v] == arc_count(v);
      });
    }
    nodes[bag].parent = bd.parent(bag);
    std::ranges::transform(used_edges[bag],
                           std::back_insert_iterator(nodes[bag].bag),
                           [](auto const &item) { return item.first; });
  };
  bd.foreach_post_order(handle);
  nodes.back().parent = -1;

  AnnotatedDecomposition rv(nodes);
  return rv;
}

std::ostream &operator<<(std::ostream &out, const AnnotatedDecomposition &ad) {
  ad.foreach_post_order([&](auto bag_idx) {
    for (auto v : ad[bag_idx].bag) {
      out << "bag(" << bag_idx + 1 << "," << v << ")." << std::endl;
    }
  });
  return out;

  out << "ad " << ad.size() << " " << ad.width() << std::endl;
  ad.foreach_post_order([&](auto bag_idx) {
    out << "b " << bag_idx << " ";
    for (auto v : ad[bag_idx].bag) {
      out << v << " ";
    }
    out << std::endl;
    for (auto child : ad.children_[bag_idx]) {
      out << bag_idx << " " << child << std::endl;
    }
  });
  return out;
}

} // namespace fpc
