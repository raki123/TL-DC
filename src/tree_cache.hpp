#pragma once

#include "treewidth_search.h"

namespace fpc {

struct TreeNode {
  public:
    std::vector<size_t> cache_idx;
    size_t children[4] = {size_t(-1), size_t(-1), size_t(-1), size_t(-1)};

    static TreeNode construct_tree(
        std::vector<TreeNode>& nodes,
        std::vector<std::pair<Frontier, std::vector<Edge_weight>>>::const_iterator entries_begin,
        std::vector<std::pair<Frontier, std::vector<Edge_weight>>>::const_iterator entries_end
    );
};

} // namespace fpc