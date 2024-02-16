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

#include "heuristics.hpp"
#include <algorithm>
#include <assert.h>
#include <iostream>
#include <ranges>

namespace srg = std::ranges;

namespace fpt {

double H1::operator()(std::size_t vertex) const { return Q_.contains(vertex); }

void H1::take(std::size_t vertex) {
  taken_[vertex] = true;

  if (S_val_[vertex]) {
    active_.insert(vertex);
  }
  for (auto neigh : neighbors_[vertex]) {
    if (--S_val_[neigh] == 0 && taken_[neigh]) {
      active_.erase(neigh);
    }
  }

  auto min_val = neighbors_.size();
  for (auto v : active_) {
    if (S_val_[v] < min_val) {
      min_val = S_val_[v];
      S_.clear();
    }
    if (S_val_[v] == min_val) {
      S_.insert(v);
    }
  }
  min_val = 0;
  for (auto v : active_) {
    for (auto neigh : neighbors_[v]) {
      if (taken_[neigh])
        continue;

      auto val = srg::count_if(neighbors_[neigh],
                               [&](auto other) { return S_.contains(other); });
      if (val > min_val) {
        min_val = val;
        P_.clear();
      }
      if (val == min_val) {
        P_.insert(neigh);
      }
    }
  }

  min_val = neighbors_.size();
  for (auto v : P_) {
    if (S_val_[v] < min_val) {
      min_val = S_val_[v];
      Q_.clear();
    }
    if (S_val_[v] == min_val) {
      Q_.insert(v);
    }
  }
}

void H1::reset() {
  srg::fill(taken_, false);
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
  active_.clear();
}

double DegreeHeuristic::operator()(std::size_t vertex) const {
  return 0.01 + srg::count_if(neighbors_[vertex], [this](std::size_t neigh) {
           return !taken_[neigh];
         });
}

void DegreeHeuristic::take(std::size_t vertex) { taken_[vertex] = true; }

void DegreeHeuristic::reset() { srg::fill(taken_, false); }

double MinFillHeuristic::operator()(std::size_t vertex) const {
  int fill = 0;
  for (auto neigh : neighbors_[vertex]) {
    for (auto other_neigh : neighbors_[vertex]) {
      if (taken_[neigh] && taken_[other_neigh]) {
        continue;
      }
      if (std::find(neighbors_[neigh].begin(), neighbors_[neigh].end(),
                    other_neigh) == neighbors_[neigh].end()) {
        ++fill;
      }
    }
  }
  return 0.01 + fill;
}

void MinFillHeuristic::take(std::size_t vertex) { taken_[vertex] = true; }

void MinFillHeuristic::reset() { srg::fill(taken_, false); }
} // namespace fpt