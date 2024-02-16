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
#include "base_solver.h"
#include "clhasher.h"
#include "count_structures.h"
#include "graph.h"
#include <limits>
#include <omp.h>
#include <unordered_map>

namespace fpc {

template <template <typename> typename Count_structure, typename count_t>
class NautyPathwidthSearch : public Solver {
  using Partial_result = Count_structure<count_t>;
  using Cache =
      std::unordered_map<TWCacheKey, Partial_result, twc_hash, twc_equal>;

public:
  NautyPathwidthSearch(Graph const &input, AnnotatedDecomposition decomposition,
                       size_t nthreads);

  virtual mpz_class search() override;
  virtual void print_stats() const override;

  const frontier_index_t invalid_index_ =
      std::numeric_limits<frontier_index_t>::max();
  const frontier_index_t no_edge_index_ =
      std::numeric_limits<frontier_index_t>::max() - 1;
  const frontier_index_t two_edge_index_ =
      std::numeric_limits<frontier_index_t>::max() - 2;

private:
  size_t nthreads_;
  Graph graph_;
  Edge_length max_length_;
  bool is_all_pair_;
  std::vector<Vertex> terminals_;
  AnnotatedDecomposition decomposition_;
  std::vector<std::vector<Vertex>> remaining_edges_after_this_;

  std::vector<std::vector<frontier_index_t>> bag_local_idx_map_;
  std::vector<std::vector<Vertex>> bag_local_vertex_map_;
  Edge_length invalid_distance_ = std::numeric_limits<Edge_length>::max();
  std::vector<std::vector<std::vector<Edge_length>>> bag_local_distance_;
  std::vector<sparsegraph> sparsegraph_after_this_;

  std::vector<count_t> thread_local_result_;

  std::vector<std::pair<Cache, Cache>> cache_;

  void propagateLoop(Frontier &frontier, size_t bag_idx, size_t last_idx,
                     Partial_result &partial_results, bool takeable,
                     bool skippable, size_t thread_id);

  void includeSolutions(Frontier const &frontier, size_t bag_idx,
                        Partial_result const &partial_result);

  bool canTake(Frontier &frontier, size_t bag_idx,
               Partial_result const &partial_result);
  bool canSkip(Frontier &frontier, size_t bag_idx,
               Partial_result const &partial_result);

  bool distancePrune(Frontier &frontier,
                     std::vector<frontier_index_t> const &paths,
                     std::vector<frontier_index_t> const &cut_paths,
                     size_t bag_idx, size_t offset);

  void take(Frontier &frontier, size_t bag_idx);
  void skip(Frontier &frontier, size_t bag_idx);

  void advance(Frontier &frontier, size_t bag_idx);

  sparsegraph construct_sparsegraph(Frontier const &frontier, size_t last_idx);

  // stats
  std::vector<size_t> pos_hits_;
  std::vector<size_t> neg_hits_;

  std::vector<size_t> edges_;
  std::vector<size_t> propagations_;
};

template class NautyPathwidthSearch<Limited_count, uint64_t>;
template class NautyPathwidthSearch<Limited_count, uint128_t>;
template class NautyPathwidthSearch<Limited_count, uint256_t>;
template class NautyPathwidthSearch<Limited_count, uint512_t>;
template class NautyPathwidthSearch<Limited_count, uint1024_t>;
template class NautyPathwidthSearch<Limited_count, mpz_class>;

template class NautyPathwidthSearch<Unlimited_count, uint64_t>;
template class NautyPathwidthSearch<Unlimited_count, uint128_t>;
template class NautyPathwidthSearch<Unlimited_count, uint256_t>;
template class NautyPathwidthSearch<Unlimited_count, uint512_t>;
template class NautyPathwidthSearch<Unlimited_count, uint1024_t>;
template class NautyPathwidthSearch<Unlimited_count, mpz_class>;
} // namespace fpc
