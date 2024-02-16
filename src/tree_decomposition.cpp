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

#include "tree_decomposition.hpp"
#include "digraph.h"
#include "graph.h"
#include "popen2.h"
#include <fstream>
#include <ios>
#include <iostream>
#include <ranges>
#include <sstream>
#include <unordered_map>

namespace fpc {
TreeDecomposition TreeDecomposition::from_graph(Graph const &graph,
                                                Strategy strategy) {
  switch (strategy) {
  case Strategy::HTD: {
    int iterations = 100'000'000 / graph.nr_edges();
    auto [in, out] = open_htd(iterations);
    write_instance(graph, in);
    return from_file(out);
  }
  case Strategy::FLOW_CUTTER: {
    auto [in, out] = open_flow_cutter();
    write_instance(graph, in);
    return from_file(out);
  }
  case Strategy::TAMAKI: {
    auto [in, out] = open_tamaki();
    write_instance(graph, in);
    return from_file(out);
  }
  default:
    assert(false);
  }
}

TreeDecomposition TreeDecomposition::from_graph(Digraph const &graph,
                                                Strategy strategy) {
  switch (strategy) {
  case Strategy::HTD: {
    int iterations = 100'000'000 / graph.nr_arcs();
    auto [in, out] = open_htd(iterations);
    write_instance(graph, in);
    return from_file(out);
  }
  case Strategy::FLOW_CUTTER: {
    auto [in, out] = open_flow_cutter();
    write_instance(graph, in);
    return from_file(out);
  }
  case Strategy::TAMAKI: {
    auto [in, out] = open_tamaki();
    write_instance(graph, in);
    return from_file(out);
  }
  default:
    assert(false);
  }
}

TreeDecomposition::Pipes TreeDecomposition::open_htd(int iterations) {
  int in, out;
  std::string iters = std::to_string(iterations);
  if (popen2("/usr/bin/timeout", &in, &out, "15", "./utils/htd",
             "--child-limit", "2", "--seed", "1", "--opt", "width",
             "--iterations", "0", "--patience", "-1", "--strategy", "min-fill",
             nullptr) <= 0) {
    throw std::runtime_error("popen failed");
  }
  return {in, out};
}

TreeDecomposition::Pipes TreeDecomposition::open_flow_cutter() {
  int in, out;
  if (popen2("/usr/bin/timeout", &in, &out, "15", "./utils/flow_cutter",
             nullptr) <= 0) {
    throw std::runtime_error("popen failed");
  }
  return {in, out};
}

TreeDecomposition::Pipes TreeDecomposition::open_tamaki() {
  int in, out;
  if (popen2("/usr/bin/timeout", &in, &out, "15", "./utils/tamaki-heuristic",
             nullptr) <= 0) {
    throw std::runtime_error("popen failed");
  }
  return {in, out};
}

void TreeDecomposition::write_instance(Graph const &graph, int in) {
  std::size_t nr_edges = graph.nr_edges();
  std::stringstream s;
  s << "p tw " << graph.nr_vertices() << " " << nr_edges << std::endl;
  if (write(in, s.str().c_str(), strlen(s.str().c_str())) == -1) {
    throw std::runtime_error("write failed");
  }
  for (Vertex v = 0; v < graph.nr_vertices(); v++) {
    for (auto neigh : graph.neighbors(v)) {
      s.str("");
      s.clear();
      if (v < neigh) {
        s << (v + 1) << " " << neigh + 1 << std::endl;
        if (write(in, s.str().c_str(), strlen(s.str().c_str())) == -1) {
          throw std::runtime_error("write failed");
        }
      }
    }
  }
  close(in);
}

void TreeDecomposition::write_instance(Digraph const &graph, int in) {
  std::size_t nr_arcs = graph.nr_arcs();
  std::stringstream s;
  s << "p tw " << graph.nr_vertices() << " " << nr_arcs << std::endl;
  if (write(in, s.str().c_str(), strlen(s.str().c_str())) == -1) {
    throw std::runtime_error("write failed");
  }
  for (Vertex v = 0; v < graph.nr_vertices(); v++) {
    for (auto neigh : graph.out_neighbors(v)) {
      s.str("");
      s.clear();
      s << (v + 1) << " " << neigh + 1 << std::endl;
      if (write(in, s.str().c_str(), strlen(s.str().c_str())) == -1) {
        throw std::runtime_error("write failed");
      }
    }
  }
  close(in);
}

TreeDecomposition TreeDecomposition::from_file(int out) {
  TreeDecomposition rv(0, 0);
  FILE *fout = fdopen(out, "r");

  char *line = NULL;
  int read = 0;
  std::size_t len;

  std::string tmp;
  std::size_t bags_left = -1;
  std::size_t edges_left = -1;
  bool found_header = false;
  while (bags_left > 0 || edges_left > 0) {
    read = getline(&line, &len, fout);
    if (read == -1) {
      if (errno != 0) {
        throw std::runtime_error("getline() failed.");
      } else {
        break;
      }
    }
    std::stringstream ss(line);
    ss >> tmp;
    if (tmp == "c") {
    } else if (tmp == "") {
    } else if (tmp == "s") {
      ss >> tmp;
      assert(tmp == "td");
      int width, nr_vertices;
      ss >> bags_left >> width >> nr_vertices;
      edges_left = bags_left > 0 ? bags_left - 1 : 0;
      rv = TreeDecomposition(bags_left, width);
      found_header = true;
    } else if (tmp == "b") {
      int bid;
      ss >> bid;
      std::vector<Vertex> label;
      Vertex v;
      while (ss >> v) {
        label.push_back(v - 1);
      }
      rv.set_label(bid - 1, label);
      bags_left--;
    } else {
      int a = stoi(tmp);
      int b;
      ss >> b;
      rv.add_edge(a - 1, b - 1);
      edges_left--;
    }
  }
  fclose(fout);
  close(out);
  free(line);
  if (!found_header) {
    throw std::runtime_error(
        "Tree Decomposition Solver did not finish successfully.");
  }
  auto root_child = [&rv]() {
    for (std::size_t v = 0; v < rv.nr_bags(); v++) {
      if (rv.neighbors_[v].size() <= 1) {
        return v;
      }
    }
    assert(false);
  }();
  rv.set_root(root_child);
  return rv;
}

void TreeDecomposition::set_root(std::size_t root) {
  std::vector<std::tuple<std::size_t, std::size_t>> stack;
  root_ = root;
  stack.emplace_back(-1, root);
  while (!stack.empty()) {
    auto [parent, cur] = stack.back();
    width_ = std::max(width_, node_label_[cur].size());
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
  }
}

void TreeDecomposition::enforce_child_limit(std::size_t max_children) {
  assert(max_children >= 2);
  std::vector<std::vector<std::size_t>> new_neighbors;
  std::vector<std::vector<Vertex>> new_label;
  std::size_t cur_idx = 1;

  struct Data {
    std::size_t new_bag;
    std::size_t neighbor;
    std::size_t bag;
  };
  std::vector<Data> stack;
  stack.emplace_back(0, 0, root_);
  new_neighbors.push_back({});
  new_label.push_back(node_label_[root_]);
  while (!stack.empty()) {
    auto [new_bag_idx, neigh_idx, bag_idx] = stack.back();
    if (neigh_idx + 2 < children_[bag_idx].size()) {
      // more than two children left

      // add the copy of the current bag
      new_neighbors.push_back({});
      new_neighbors[new_bag_idx].push_back(cur_idx);
      new_neighbors[cur_idx].push_back(new_bag_idx);
      new_label.push_back(node_label_[bag_idx]);
      stack.back().new_bag = cur_idx;
      stack.back().neighbor++;
      cur_idx++;

      // add the child
      new_neighbors.push_back({});
      new_neighbors[new_bag_idx].push_back(cur_idx);
      new_neighbors[cur_idx].push_back(new_bag_idx);
      auto next = children_[bag_idx][neigh_idx];
      new_label.push_back(node_label_[next]);
      stack.emplace_back(cur_idx, 0, next);
      cur_idx++;
    } else if (neigh_idx < children_[bag_idx].size()) {
      // two children or less left

      stack.back().neighbor++;
      // add the child
      new_neighbors.push_back({});
      new_neighbors[new_bag_idx].push_back(cur_idx);
      new_neighbors[cur_idx].push_back(new_bag_idx);
      auto next = children_[bag_idx][neigh_idx];
      new_label.push_back(node_label_[next]);
      stack.emplace_back(cur_idx, 0, next);
      cur_idx++;
    } else {
      // nothing left. self has been handled before being added to stack
      stack.pop_back();
    }
  }
  assert(cur_idx == new_neighbors.size());
  assert(cur_idx == new_label.size());
  root_ = 0;
  nr_bags_ = cur_idx;
  neighbors_ = new_neighbors;
  children_ = std::vector<std::vector<std::size_t>>(nr_bags_);
  parent_ = std::vector<std::size_t>(nr_bags_, -1);
  node_label_ = new_label;
  set_root(root_);
}

std::ostream &operator<<(std::ostream &out, const TreeDecomposition &td) {
  out << "s td " << td.nr_bags() << " " << td.width() << " ??\n";
  td.foreach_post_order([&](auto bag_idx) {
    out << "b " << bag_idx << " ";
    for (auto v : td[bag_idx]) {
      out << v << " ";
    }
    out << "\n";
  });
  td.foreach_post_order([&](auto bag_idx) {
    for (auto child : td.children(bag_idx)) {
      out << bag_idx << " " << child << "\n";
    }
  });

  return out;
}

} // namespace fpc
