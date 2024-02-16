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

#include "canonizer.h"

namespace fpc {

void Canonizer::operator()(Vertex &start,
                           std::vector<Edge_length> &distance) const {
  for (auto const &eq : equivalence_classes_) {
    Vertex valid = 0;
    bool found_start = false;
    Edge_length value = invalid_;
    for (Vertex v : eq) {
      assert(v != start || distance[v] == invalid_);
      if (distance[v] != invalid_) {
        valid++;
        assert(value == invalid_ || value == distance[v]);
        value = distance[v];
        distance[v] = invalid_;
      }
      found_start = found_start || start == v;
    }
    if (found_start) {
      // since start has distance invalid_, this does not go out of bounds
      assert(valid < eq.size());
      start = eq[valid];
    }
    while (valid > 0) {
      --valid;
      distance[eq[valid]] = value;
    }
  }
}
} // namespace fpc