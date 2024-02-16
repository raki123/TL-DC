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
#include "digraph.h"
#include "graph.h"
#include "pathwidth/heuristics.hpp"
#include "pathwidth/path_decomposer.hpp"
#include "tree_decomposition.hpp"
#include <chrono>
#include <vector>

namespace fpc {
class AnnotatedDecomposition;

template <typename GraphType>
AnnotatedDecomposition make_via_path_decomposition(std::size_t nr_threads,
                                                   GraphType const &graph) {

  std::vector<std::vector<std::size_t>> neighbors;
  for (Vertex v = 0; v < graph.nr_vertices(); v++) {
    std::vector<std::size_t> neighs;
    if constexpr (std::same_as<GraphType, Graph>) {
      for (auto neigh : graph.neighbors(v)) {
        neighs.push_back(neigh);
      }
    } else {
      static_assert(std::same_as<GraphType, Digraph>);
      for (auto neigh : graph.out_neighbors(v)) {
        neighs.push_back(neigh);
      }
      for (auto neigh : graph.in_neighbors(v)) {
        neighs.push_back(neigh);
      }
      srg::sort(neighs);
      auto [end, last] = srg::unique(neighs);
      neighs.erase(end, last);
    }
    neighbors.push_back(neighs);
  }

  std::vector<std::size_t> thread_local_min_width(nr_threads, std::size_t(-1));
  std::vector<std::vector<std::size_t>> thread_local_best_order(nr_threads);

#pragma omp parallel for schedule(static, 1) num_threads(nr_threads)
  for (auto id = 0; id < nr_threads; id++) {
    fpt::PathDecomposer<fpt::NaturalOrderHeuristic> natural(neighbors);
    thread_local_min_width[id] = natural.make_order();
    thread_local_best_order[id] = natural.get_order();
    fpt::PathDecomposer<fpt::H1> path_decomposer(neighbors);
    auto start_time = std::chrono::high_resolution_clock::now();
    auto now = std::chrono::high_resolution_clock::now();
    auto duration =
        std::chrono::duration_cast<std::chrono::seconds>(now - start_time)
            .count();
    while (duration < 15) {
      auto width = path_decomposer.make_probabilistic_order();
      if (width < thread_local_min_width[id]) {
        thread_local_min_width[id] = width;
        thread_local_best_order[id] = path_decomposer.get_order();
      }
      now = std::chrono::high_resolution_clock::now();
      duration =
          std::chrono::duration_cast<std::chrono::seconds>(now - start_time)
              .count();
    }
  }
  std::size_t min_width = std::size_t(-1);
  std::vector<std::size_t> best_order;
  for (auto id = 0; id < nr_threads; id++) {
    if (thread_local_min_width[id] < min_width) {
      min_width = thread_local_min_width[id];
      best_order = std::move(thread_local_best_order[id]);
    }
  }

  auto decomposition =
      AnnotatedDecomposition::from_vertex_order(best_order, graph);
  decomposition.foreach_post_order([&](auto bag_idx) {
    assert(decomposition[bag_idx].bag.size() <= min_width);
  });

  return decomposition;
}

template <typename GraphType>
AnnotatedDecomposition make_via_tree_decomposition(std::size_t,
                                                   GraphType const &graph) {

  auto td =
      TreeDecomposition::from_graph(graph, TreeDecomposition::Strategy::TAMAKI);
  td.enforce_child_limit(2);

  auto decomposition =
      AnnotatedDecomposition::from_tree_decomposition(td, graph);
  decomposition.foreach_post_order([&](auto bag_idx) {
    assert(decomposition[bag_idx].bag.size() <= td.width());
  });

  return decomposition;
}
} // namespace fpc