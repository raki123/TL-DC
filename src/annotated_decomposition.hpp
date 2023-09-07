// TLDC "Too long; Didn't Count" A length limited path counter.
// Copyright (C) 2023 Rafael Kiesel, Markus Hecher

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

enum NodeType { LEAF, PATH_LIKE, JOIN };

struct AnnotatedNode {
public:
  size_t parent; // parent of the node or -1 if root
  NodeType type;
  Edge edge;                          // if LEAF or PATH_LIKE
  std::pair<size_t, size_t> children; // if JOIN the indices of the child nodes
                                      // in the annotated decomposition
  std::vector<vertex_t> bag;

  void stats() const;
};

// contains annotated nodes
// path like nodes continuously from leaf to last node before join
// join node occurs after both children nodes
typedef std::vector<AnnotatedNode> AnnotatedDecomposition;

} // namespace fpc
