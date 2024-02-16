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

#include "count_structures.h"
#include "types.h"
#include <memory>

namespace fpc {
class Solver {
public:
  virtual ~Solver() = default;
  virtual mpz_class search() = 0;
  virtual void print_stats() const = 0;
};

typedef std::unique_ptr<Solver> SolverPtr;

namespace {

template <template <template <typename> typename, typename> typename SolverType,
          template <typename> typename CountStructure, typename GraphType,
          typename... Args>
inline SolverPtr make_bit_bounded_solver(GraphType const &graph, Args... args) {
  auto bit_bound = graph.bit_upper_bound();
  if (bit_bound <= 64) {
    return SolverPtr(new SolverType<CountStructure, uint64_t>(graph, args...));
  } else if (bit_bound <= 128) {
    return SolverPtr(new SolverType<CountStructure, uint128_t>(graph, args...));
  } else if (bit_bound <= 256) {
    return SolverPtr(new SolverType<CountStructure, uint256_t>(graph, args...));
  } else if (bit_bound <= 512) {
    return SolverPtr(new SolverType<CountStructure, uint512_t>(graph, args...));
  } else if (bit_bound <= 1024) {
    return SolverPtr(
        new SolverType<CountStructure, uint1024_t>(graph, args...));
  } else {
    return SolverPtr(new SolverType<CountStructure, mpz_class>(graph, args...));
  }
}
} // namespace

template <template <template <typename> typename, typename> typename SolverType,
          typename GraphType, typename... Args>
inline SolverPtr make_solver(GraphType const &graph, Args... args) {
  if (graph.is_length_limited()) {
    return make_bit_bounded_solver<SolverType, Limited_count>(graph, args...);
  } else {
    return make_bit_bounded_solver<SolverType, Unlimited_count>(graph, args...);
  }
}

} // namespace fpc