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
#include "graph.h"

namespace fpc {

class Canonizer {
public:
  Canonizer(std::vector<std::vector<Vertex>> equivalence_classes,
            Edge_length invalid, Vertex end)
      : equivalence_classes_(std::move(equivalence_classes)),
        invalid_(invalid) {
    for (auto &eq : equivalence_classes_) {
      if (auto it = std::find(eq.begin(), eq.end(), end); it != eq.end()) {
        std::swap(*it, eq.back());
        eq.pop_back();
      }
    }
    for (auto i = 0; i < equivalence_classes_.size(); ++i) {
      if (equivalence_classes_[i].size() == 1) {
        std::swap(equivalence_classes_.back(), equivalence_classes_[i]);
        equivalence_classes_.pop_back();
        --i;
      }
    }
  }

  static Canonizer noop_canonizer() {
    return Canonizer({}, std::numeric_limits<Edge_length>::max(), 0);
  }

  void operator()(Vertex &start,
                  std::vector<Edge_length> &distance_to_goal) const;

private:
  std::vector<std::vector<Vertex>> equivalence_classes_;
  Edge_length invalid_;
};
} // namespace fpc