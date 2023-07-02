#pragma once

#include "graph.h"

namespace fpc {

enum NodeType {
    LEAF,
    PATH_LIKE,
    JOIN
};

struct NodeAnnotation {
  public:
    NodeType type;
    Edge edge; // if LEAF or PATH_LIKE
    std::pair<size_t, size_t> children; // if JOIN the indices of the child nodes in the annotated decomposition
};


// contains annotated nodes
// path like nodes continuously from leaf to last node before join
// join node occurs after both children nodes
typedef std::vector<std::pair<NodeAnnotation, std::vector<vertex_t>>> AnnotatedDecomposition;

} // namespace