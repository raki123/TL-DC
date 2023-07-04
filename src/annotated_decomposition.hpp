#pragma once

#include "graph.h"

namespace fpc {

enum NodeType {
    LEAF,
    PATH_LIKE,
    JOIN
};

struct AnnotatedNode {
  public:
    size_t parent; // parent of the node or -1 if root
    NodeType type;
    Edge edge; // if LEAF or PATH_LIKE
    std::pair<size_t, size_t> children; // if JOIN the indices of the child nodes in the annotated decomposition
    std::vector<vertex_t> bag;

    void stats() const;
};


// contains annotated nodes
// path like nodes continuously from leaf to last node before join
// join node occurs after both children nodes
typedef std::vector<AnnotatedNode> AnnotatedDecomposition;

} // namespace
