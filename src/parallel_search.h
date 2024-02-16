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
#include "base_solver.h"
#include "canonizer.h"
#include "clhasher.h"
#include "graph.h"
#include <omp.h>
#include <optional>
#include <unordered_map>
#include <utility>

namespace fpc {

template <template <typename> typename Count_structure, typename count_t>
class ParallelSearch : public Solver {
  struct CacheKey {
    Vertex start;
    std::vector<Edge_length> distance_to_goal;
  };
  struct key_hash {
    size_t operator()(CacheKey const &key) const {
      return key.start + hasher__(key.distance_to_goal);
    }
  };
  struct key_eq {
    size_t operator()(CacheKey const &lhs, CacheKey const &rhs) const {
      return lhs.start == rhs.start &&
             lhs.distance_to_goal == rhs.distance_to_goal;
    }
  };
  using Partial_result = Count_structure<count_t>;
  using Cache = std::unordered_map<CacheKey, Partial_result, key_hash, key_eq>;

public:
  ParallelSearch(Graph const &input, std::optional<Vertex> must_use,
                 size_t nthreads);

  virtual mpz_class search() override;

  virtual void print_stats() const override;

private:
  size_t nthreads_;
  bool enable_dag_;
  Edge_length max_length_;
  std::vector<Vertex> terminals_;
  std::vector<std::vector<Vertex>> neighbors_;
  std::optional<Vertex> must_use_;

  Edge_length invalid_;
  std::vector<std::vector<Edge_length>> distance_;

  Canonizer canonizer_;
  std::vector<Cache> cache_;

  std::vector<count_t> thread_local_result_;

  // more efficient search if budget is equal to length of shortest path
  Edge_weight dag_search(Vertex start,
                         std::vector<Edge_length> const &distance_to_goal,
                         bool used_must_use);

  void prune_articulation(Vertex start, std::vector<Edge_length> &distance);
  bool ap_util(Vertex u, std::vector<char> &visited, std::vector<Vertex> &disc,
               std::vector<Vertex> &low, int &time, int parent, Vertex start,
               std::vector<Edge_length> &distance);
  void prune_util(Vertex u, std::vector<Edge_length> &distance);
  // void component_util(Vertex u, std::vector<Vertex>& disc);

  // // ap datastructures
  // std::vector<Vertex> ap_disc_;
  // std::vector<Vertex> ap_low_;
  // std::vector<char> ap_visited_;

  // // for splitting based on articulation points between start and goal
  // Vertex last_ap_;
  // std::vector<std::pair<Vertex, Vertex>> ap_start_goal_;
  // std::vector<std::vector<Vertex>> ap_components_;

  // helper functions
  std::vector<Vertex> const &neighbors(Vertex v) {
    assert(v < neighbors_.size());
    return neighbors_[v];
  };
  void dijkstra(Vertex start, std::vector<Edge_length> &distance);
  void pruning_dijkstra(Vertex start, Vertex prune,
                        std::vector<Edge_length> &distance,
                        std::vector<Edge_length> const &old_distance,
                        Edge_length budget);

  // stats
  std::vector<size_t> pos_hits_;
  std::vector<size_t> neg_hits_;

  // size_t splits = 0;

  std::vector<size_t> edges_;
  std::vector<size_t> propagations_;
  std::vector<size_t> dags_;
};

template class ParallelSearch<Limited_count, uint64_t>;
template class ParallelSearch<Limited_count, uint128_t>;
template class ParallelSearch<Limited_count, uint256_t>;
template class ParallelSearch<Limited_count, uint512_t>;
template class ParallelSearch<Limited_count, uint1024_t>;
template class ParallelSearch<Limited_count, mpz_class>;

template class ParallelSearch<Unlimited_count, uint64_t>;
template class ParallelSearch<Unlimited_count, uint128_t>;
template class ParallelSearch<Unlimited_count, uint256_t>;
template class ParallelSearch<Unlimited_count, uint512_t>;
template class ParallelSearch<Unlimited_count, uint1024_t>;
template class ParallelSearch<Unlimited_count, mpz_class>;
} // namespace fpc