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
#include "base_solver.h"
#include "clhasher.h"
#include "count_structures.h"
#include "digraph.h"
#include <omp.h>
#include <utility>

namespace fpc {

// partial result is not actually relevant but we have it so that make_solver
// works with it
template <template <typename> typename Count_structure, typename count_t>
class MITMDirectedSearch : public Solver {
  using CacheKey = std::vector<Vertex>;
  using Cache = ankerl::unordered_dense::map<CacheKey, Limited_count<count_t>,
                                             vector_hash<Vertex>>;

public:
  MITMDirectedSearch(Digraph const &input, size_t nthreads);

  virtual mpz_class search() override;

  virtual void print_stats() const override;

  std::pair<double, double> percent_done_after(std::size_t paths) const;

private:
  Vertex nr_vertices_;
  size_t nthreads_;
  Edge_length max_length_;
  std::vector<Vertex> terminals_;
  Digraph const &graph_;
  std::vector<std::vector<Vertex>> in_neighbors_;
  std::vector<std::vector<Vertex>> out_neighbors_;

  Edge_length invalid_;
  std::vector<Edge_length> distance_from_start_;
  std::vector<Edge_length> distance_to_goal_;

  std::vector<Cache> cache_;

  // result per middle vertex
  std::vector<count_t> local_result_;

  // helper functions
  std::vector<Vertex> const &in_neighbors(Vertex v) const {
    assert(v < in_neighbors_.size());
    return in_neighbors_[v];
  }

  std::vector<Vertex> const &out_neighbors(Vertex v) const {
    assert(v < out_neighbors_.size());
    return out_neighbors_[v];
  }

  void insert_all_subsets(std::vector<Vertex> const &path, Vertex mid,
                          Edge_length len);
  void enumerate_ms_paths(Vertex mid, Edge_length ms_max_len,
                          Edge_length mt_max_len);

  void inclusion_exclusion(std::vector<Vertex> const &path, Vertex mid,
                           Edge_length len);
  void enumerate_mt_paths(Vertex mid, Edge_length ms_max_len,
                          Edge_length mt_max_len);
};

template class MITMDirectedSearch<Limited_count, uint64_t>;
template class MITMDirectedSearch<Limited_count, uint128_t>;
template class MITMDirectedSearch<Limited_count, uint256_t>;
template class MITMDirectedSearch<Limited_count, uint512_t>;
template class MITMDirectedSearch<Limited_count, uint1024_t>;
template class MITMDirectedSearch<Limited_count, mpz_class>;

template class MITMDirectedSearch<Unlimited_count, uint64_t>;
template class MITMDirectedSearch<Unlimited_count, uint128_t>;
template class MITMDirectedSearch<Unlimited_count, uint256_t>;
template class MITMDirectedSearch<Unlimited_count, uint512_t>;
template class MITMDirectedSearch<Unlimited_count, uint1024_t>;
template class MITMDirectedSearch<Unlimited_count, mpz_class>;
} // namespace fpc