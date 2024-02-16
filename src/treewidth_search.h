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

#include "ankerl/unordered_dense.h"
#include "annotated_decomposition.hpp"
#include "base_solver.h"
#include "clhasher.h"
#include "count_structures.h"
#include "graph.h"
#include <array>
#include <limits>
#include <omp.h>
#include <unordered_map>

namespace fpc {

template <template <typename> typename Count_structure, typename count_t>
class TreewidthSearch : public Solver {
  using Partial_result = Count_structure<count_t>;
  using Cache =
      ankerl::unordered_dense::map<nFrontier, Partial_result, nFrontierHash>;

public:
  TreewidthSearch(Graph const &input, AnnotatedDecomposition decomposition,
                  size_t nthreads);

  virtual mpz_class search() override;
  virtual void print_stats() const override;

private:
  size_t nthreads_;
  Graph graph_;
  std::size_t max_width_;
  Edge_length max_length_;
  bool is_all_pair_;
  std::vector<Vertex> terminals_;
  std::array<frontier_index_t, 2> terminals_idx_;
  AnnotatedDecomposition decomposition_;
  std::vector<std::optional<Partial_result>> edge_weights_;
  std::vector<std::vector<Vertex>> remaining_edges_after_this_;
  std::vector<Vertex> remaining_vertices_after_this_;
  std::vector<std::vector<frontier_index_t>> eliminated_after_this_;

  frontier_index_t invalid_index_ =
      std::numeric_limits<frontier_index_t>::max();
  frontier_index_t no_edge_index_ =
      std::numeric_limits<frontier_index_t>::max() - 1;
  frontier_index_t two_edge_index_ =
      std::numeric_limits<frontier_index_t>::max() - 2;
  Edge_length invalid_distance_ = std::numeric_limits<Edge_length>::max();
  std::vector<std::vector<std::vector<Edge_length>>> bag_local_distance_;

  std::vector<count_t> thread_local_result_;

  std::size_t mod_val_;
  std::vector<omp_lock_t> mod_lock_;
  std::vector<std::pair<Cache, Cache>> cur_cache_;
  std::vector<Cache> next_cache_;
  std::vector<std::vector<std::pair<Cache, Cache>>> join_cache_;

  void cache(nFrontier &frontier, size_t new_idx,
             Partial_result &partial_results, size_t thread_id);

  void includeSolutions(nFrontier const &frontier,
                        Partial_result const &partial_result, size_t bag_idx);

  bool canTake(nFrontier &frontier, size_t bag_idx,
               Partial_result &partial_result);
  bool canSkip(nFrontier &frontier, size_t bag_idx,
               Partial_result const &partial_result);

  bool distancePrune(nFrontier const &frontier,
                     std::vector<frontier_index_t> const &paths,
                     std::vector<frontier_index_t> const &cut_paths,
                     size_t bag_idx, size_t offset);

  bool merge(nFrontier &left, nFrontier const &right, size_t bag_idx,
             Partial_result &left_result, Partial_result const &right_result);

  void advance(nFrontier &frontier, size_t bag_idx);

  // stats
  std::vector<size_t> pos_hits_;
  std::vector<size_t> neg_hits_;

  // size_t splits = 0;

  std::vector<size_t> edges_;
  std::vector<size_t> propagations_;
  std::vector<size_t> merges_;
  std::vector<size_t> unsuccessful_merges_;
};

template class TreewidthSearch<Limited_count, uint64_t>;
template class TreewidthSearch<Limited_count, uint128_t>;
template class TreewidthSearch<Limited_count, uint256_t>;
template class TreewidthSearch<Limited_count, uint512_t>;
template class TreewidthSearch<Limited_count, uint1024_t>;
template class TreewidthSearch<Limited_count, mpz_class>;

template class TreewidthSearch<Unlimited_count, uint64_t>;
template class TreewidthSearch<Unlimited_count, uint128_t>;
template class TreewidthSearch<Unlimited_count, uint256_t>;
template class TreewidthSearch<Unlimited_count, uint512_t>;
template class TreewidthSearch<Unlimited_count, uint1024_t>;
template class TreewidthSearch<Unlimited_count, mpz_class>;

} // namespace fpc
