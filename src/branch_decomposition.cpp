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

#include "branch_decomposition.hpp"
#include "popen2.h"
#include <ios>
#include <ranges>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

namespace fpc {
BranchDecomposition BranchDecomposition::from_graph(Graph const &graph,
                                                    Strategy strategy) {
  switch (strategy) {
  case Strategy::HICKS:
    return hicks(graph);
  default:
    assert(false);
  }
}

BranchDecomposition BranchDecomposition::hicks(Graph const &graph) {
  Vertex nr_edges = graph.nr_edges();
  std::vector<Edge> id_to_edge(nr_edges);
  std::size_t edge_id = 0;
  int in, out;
  BranchDecomposition rv(0, 0, 0);
  Graph decomp_tree(2 * nr_edges - 1);
  if (popen2("./utils/bw", &in, &out, nullptr) <= 0) {
    throw std::runtime_error("popen failed");
  }
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
        id_to_edge[edge_id++] = Edge(v, neigh);
        s << (v + 1) << " " << neigh + 1 << std::endl;
        if (write(in, s.str().c_str(), strlen(s.str().c_str())) == -1) {
          throw std::runtime_error("write failed");
        }
      }
    }
  }
  close(in);
  FILE *fout = fdopen(out, "r");

  char *line = NULL;
  int read = 0;
  std::size_t len;

  std::string tmp;
  std::size_t bags_left = nr_edges;
  std::size_t edges_left = 2 * nr_edges - 3;
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
      assert(tmp == "bd");
      int es, bs, w, des;
      ss >> es >> bs >> w >> des;
      assert(es == nr_edges);
      assert(des == edges_left);
      rv = BranchDecomposition(bs + 1, es, w);
    } else if (tmp == "b") {
      int bid;
      ss >> bid;
      int e;
      ss >> e;
      rv.set_label(bid - 1, id_to_edge[e - 1]);
      bags_left--;
    } else {
      int a = stoi(tmp);
      int b;
      ss >> b;
      decomp_tree.add_edge(a - 1, b - 1);
      edges_left--;
    }
  }
  fclose(fout);
  close(out);
  free(line);

  auto root_child = [&decomp_tree]() {
    for (Vertex v = 0; v < decomp_tree.nr_vertices(); v++) {
      if (decomp_tree.neighbors(v).size() == 1) {
        return v;
      }
    }
    assert(false);
  }();
  Vertex root = 2 * nr_edges - 2;
  assert(decomp_tree.neighbors(root).empty());
  assert(rv.nr_bags() == root + 1);
  auto root_neighbor = decomp_tree.neighbors(root_child).front();
  decomp_tree.remove_edge(root_child, root_neighbor);
  decomp_tree.add_edge(root, root_child);
  decomp_tree.add_edge(root, root_neighbor);
  rv.set_root(root);
  std::vector<Vertex> stack = {root};
  std::unordered_set<Vertex> seen = {root};
  while (!stack.empty()) {
    auto cur = stack.back();
    stack.pop_back();
    for (auto child : decomp_tree.neighbors(cur)) {
      auto [_, inserted] = seen.emplace(child);
      if (!inserted) {
        continue;
      }
      rv.add_arc(cur, child);
      stack.push_back(child);
      seen.insert(child);
    }
  }
  return rv;
}

} // namespace fpc
