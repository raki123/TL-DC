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
#include "clhasher.h"
#include "graph.h"
#include <omp.h>
#include <unordered_map>
#include <utility>

namespace fpc {

class ParallelNautySearch {
  using Cache = std::unordered_map<NPCacheKey, std::vector<Edge_weight>,
                                   sg_hash, sg_equal>;

public:
  ParallelNautySearch(sparsegraph input, Edge_length max_length,
                      size_t nthreads);
  ~ParallelNautySearch() {
    for (size_t i = 0; i < nthreads_; i++) {
      free(thread_local_lab_[i]);
      free(thread_local_sg_[i].v);
    }
    free(initial_.v);
  }

  void add_to_initial(std::vector<sparsegraph> &to_add);

  std::vector<Edge_weight> search();

  void print_stats();

private:
  size_t nthreads_;
  Edge_length max_length_;
  Edge_length invalid_;

  sparsegraph initial_;

  std::vector<Cache> cache_;

  std::vector<Edge_weight> result_;
  std::vector<std::vector<Edge_weight>> thread_local_result_;

  // reusable data structures for nauty calls
  std::vector<sparsegraph> thread_local_sg_;
  std::vector<int *> thread_local_lab_;
  std::vector<int *> thread_local_ptn_;
  std::vector<int *> thread_local_orbits_;

  // more efficient search if budget is equal to length of shortest path
  Edge_weight dag_search(sparsegraph const &sg, Vertex start,
                         std::vector<Edge_length> const &distance_to_goal);

  void prune_articulation(sparsegraph const &sg, Vertex start,
                          std::vector<Edge_length> &distance);
  bool ap_util(sparsegraph const &sg, Vertex u, std::vector<char> &visited,
               std::vector<Vertex> &disc, std::vector<Vertex> &low, int &time,
               int parent, Vertex start, std::vector<Edge_length> &distance);
  void prune_util(sparsegraph const &sg, Vertex u,
                  std::vector<Edge_length> &distance);

  // // ap datastructures
  // std::vector<Vertex> ap_disc_;
  // std::vector<Vertex> ap_low_;
  // std::vector<char> ap_visited_;

  // for merging bridges between start and goal
  std::vector<std::vector<std::pair<Vertex, Vertex>>> thread_local_bridges_;

  // helper functions
  void pruning_dijkstra(sparsegraph const &sg,
                        std::vector<Edge_length> &distance, Edge_length budget);
  void reverse_pruning_dijkstra(
      sparsegraph const &sg, Vertex prune, std::vector<Edge_length> &distance,
      std::vector<Edge_length> const &forward_distance, Edge_length budget);

  // stats
  std::vector<size_t> pos_hits_;
  std::vector<size_t> neg_hits_;

  // size_t splits = 0;

  std::vector<size_t> edges_;
  std::vector<size_t> propagations_;
  std::vector<size_t> dags_;
  std::vector<size_t> bridges_;
};

} // namespace fpc