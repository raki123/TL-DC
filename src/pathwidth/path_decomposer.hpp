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
#include <algorithm>
#include <assert.h>
#include <iostream>
#include <map>
#include <optional>
#include <random>
#include <set>
#include <vector>

namespace fpt {

template <typename T>
concept is_heuristic = requires(T t) {
  { t.operator()(std::size_t()) } -> std::same_as<double>;
  { t.take(std::size_t()) } -> std::same_as<void>;
  { t.reset() } -> std::same_as<void>;
};

template <typename Heuristic>
requires is_heuristic<Heuristic>
class PathDecomposer {
public:
  PathDecomposer(std::vector<std::vector<std::size_t>> const &neighbors,
                 std::optional<std::size_t> first = std::nullopt,
                 std::optional<std::size_t> last = std::nullopt)
      : heuristic_(neighbors), first_(first), last_(last), n_(neighbors.size()),
        neighbors_(neighbors), taken_(n_, false), remaining_neighbors_(n_, 0) {
    assert(std::all_of(neighbors_.begin(), neighbors_.end(),
                       [](auto const &neighs) { return !neighs.empty(); }));
  }

  std::size_t make_order() {
    heuristic_.reset();
    if (n_ == 0) {
      return {};
    }
    std::fill(taken_.begin(), taken_.end(), false);
    for (std::size_t i = 0; i < n_; i++) {
      remaining_neighbors_[i] = neighbors_[i].size();
    }
    active_.clear();
    order_.clear();
    width_ = 0;
    std::size_t remaining = n_;
    take(first_ ? *first_ : choose_first());
    while (--remaining) {
      take(choose_next());
    }
    for (std::size_t i = 0; i < n_; i++) {
      assert(taken_[i]);
    }
    assert(order_.size() == n_);
    return width_;
  }

  std::size_t make_probabilistic_order() {
    heuristic_.reset();
    if (n_ == 0) {
      return {};
    }
    std::fill(taken_.begin(), taken_.end(), false);
    for (std::size_t i = 0; i < n_; i++) {
      remaining_neighbors_[i] = neighbors_[i].size();
    }
    active_.clear();
    order_.clear();
    width_ = 0;
    std::size_t remaining = n_;
    take(first_ ? *first_ : probabilistic_choose_first());
    while (--remaining) {
      take(probabilistic_choose_next());
    }
    for (std::size_t i = 0; i < n_; i++) {
      assert(taken_[i]);
    }
    assert(order_.size() == n_);
    return width_;
  }

  std::vector<std::size_t> get_order() const {
    assert(order_.size() == n_);
    return order_;
  }

private:
  void take(std::size_t v) {
    assert(v < n_);
    assert(!taken_[v]);
    heuristic_.take(v);
    taken_[v] = true;
    if (remaining_neighbors_[v]) {
      active_.insert(v);
      width_ = std::max(width_, active_.size());
    } else {
      width_ = std::max(width_, active_.size() + 1);
    }
    for (auto neigh : neighbors_[v]) {
      if (--remaining_neighbors_[neigh] == 0 && taken_[neigh]) {
        active_.erase(neigh);
        if (!first_ && !last_) {
          order_.push_back(neigh);
        }
      }
    }
    if (!remaining_neighbors_[v] || first_ || last_) {
      order_.push_back(v);
    }
  }

  std::size_t choose_first() {
    std::size_t best = -1;
    double best_val = []() {
      if constexpr (Heuristic::MaxSemantics) {
        return -std::numeric_limits<double>::infinity();
      } else {
        return std::numeric_limits<double>::infinity();
      }
    }();
    for (std::size_t i = 0; i < n_; i++) {
      auto val = heuristic_(i);
      if constexpr (Heuristic::MaxSemantics) {
        if (val > best_val) {
          best = i;
          best_val = val;
        }
      } else {
        if (val < best_val) {
          best = i;
          best_val = val;
        }
      }
    }
    assert(best != std::size_t(-1));
    return best;
  }

  std::size_t probabilistic_choose_first() {
    double sum = 0;
    std::vector<double> values_(neighbors_.size(), 0);
    for (std::size_t i = 0; i < n_; i++) {
      values_[i] = heuristic_(i);
      if constexpr (!Heuristic::MaxSemantics) {
        values_[i] = 1 / values_[i];
      }
      sum += values_[i];
    }

    assert(sum > 0);
    std::uniform_real_distribution<double> distribution(0.0, sum);
    auto chosen = distribution(rng_);
    sum = 0;
    for (std::size_t i = 0; i < n_; i++) {
      sum += values_[i];
      if (sum > chosen) {
        return i;
      }
    }
    assert(false);
  }

  std::size_t choose_next() {
    std::size_t best = -1;
    double best_val = []() {
      if constexpr (Heuristic::MaxSemantics) {
        return -std::numeric_limits<double>::infinity();
      } else {
        return std::numeric_limits<double>::infinity();
      }
    }();
    for (auto v : active_) {
      for (auto neigh : neighbors_[v]) {
        if (taken_[neigh]) {
          continue;
        }
        if (last_ && neigh == *last_) {
          continue;
        }
        auto val = heuristic_(neigh);
        if constexpr (Heuristic::MaxSemantics) {
          if (val > best_val) {
            best = neigh;
            best_val = val;
          }
        } else {
          if (val < best_val) {
            best = neigh;
            best_val = val;
          }
        }
      }
    }
    if (best == std::size_t(-1)) {
      assert(last_);
      return *last_;
    }
    return best;
  }

  std::size_t probabilistic_choose_next() {
    double sum = 0;
    std::map<std::size_t, double> values_;
    for (auto v : active_) {
      for (auto neigh : neighbors_[v]) {
        if (taken_[neigh]) {
          continue;
        }
        if (last_ && neigh == *last_) {
          continue;
        }
        auto val = heuristic_(neigh);
        if constexpr (!Heuristic::MaxSemantics) {
          val = 1 / val;
        }
        auto [it, inserted] = values_.try_emplace(neigh, val);
        if (inserted) {
          sum += it->second;
        }
      }
    }

    std::uniform_real_distribution<double> distribution(0.0, sum);
    auto chosen = distribution(rng_);
    sum = 0;
    for (auto [neigh, val] : values_) {
      sum += val;
      if (sum > chosen) {
        return neigh;
      }
    }
    assert(last_);
    return *last_;
  }
  Heuristic heuristic_;
  std::optional<std::size_t> first_;
  std::optional<std::size_t> last_;
  // static data
  std::size_t const n_;
  std::vector<std::vector<std::size_t>> const neighbors_;
  // dynamic data
  std::size_t width_;
  std::vector<char> taken_;
  std::vector<std::size_t> remaining_neighbors_;
  std::set<std::size_t> active_;
  std::vector<std::size_t> order_;

  // random
  std::mt19937 rng_;
};
} // namespace fpt
