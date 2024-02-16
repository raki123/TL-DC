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
#include <set>
#include <vector>

namespace fpt {

// https://arxiv.org/abs/1702.05710
struct H1 {
  static constexpr bool MaxSemantics = true;
  H1(std::vector<std::vector<std::size_t>> const &neighbors)
      : neighbors_(neighbors), taken_(neighbors.size(), false),
        S_val_(neighbors.size(), 0) {
    auto min_val = neighbors_.size();
    for (std::size_t i = 0; i < neighbors_.size(); i++) {
      S_val_[i] = neighbors_[i].size();
      if (S_val_[i] < min_val) {
        min_val = S_val_[i];
        Q_.clear();
      }
      if (S_val_[i] == min_val) {
        Q_.insert(i);
      }
    }
  };
  double operator()(std::size_t vertex) const;
  void take(std::size_t vertex);
  void reset();

private:
  std::vector<std::vector<std::size_t>> const &neighbors_;
  std::vector<char> taken_;
  std::set<std::size_t> active_;
  std::vector<std::size_t> S_val_;
  std::set<std::size_t> S_;
  std::set<std::size_t> P_;
  std::set<std::size_t> Q_;
};

struct DegreeHeuristic {
  static constexpr bool MaxSemantics = false;
  DegreeHeuristic(std::vector<std::vector<std::size_t>> const &neighbors)
      : neighbors_(neighbors), taken_(neighbors.size(), false) {
    max_degree_ = 0;
    for (auto &neighs : neighbors_) {
      max_degree_ = std::max(max_degree_, neighs.size());
    }
  };
  double operator()(std::size_t vertex) const;
  void take(std::size_t vertex);
  void reset();

private:
  std::vector<std::vector<std::size_t>> const &neighbors_;
  std::vector<char> taken_;
  std::size_t max_degree_;
};

struct MinFillHeuristic {
  static constexpr bool MaxSemantics = false;
  MinFillHeuristic(std::vector<std::vector<std::size_t>> const &neighbors)
      : neighbors_(neighbors), taken_(neighbors.size(), false){};
  double operator()(std::size_t vertex) const;
  void take(std::size_t vertex);
  void reset();

private:
  std::vector<std::vector<std::size_t>> const &neighbors_;
  std::vector<char> taken_;
};

struct NaturalOrderHeuristic {
  static constexpr bool MaxSemantics = false;
  NaturalOrderHeuristic(std::vector<std::vector<std::size_t>> const &){};
  double operator()(std::size_t vertex) const { return 0.01 + vertex; }
  void take(std::size_t) {}
  void reset() {}
};
} // namespace fpt