#include "tree_cache.hpp"

namespace fpc {

TreeNode TreeNode::construct_tree(
        std::vector<TreeNode>& nodes,
        std::vector<std::pair<Frontier, std::vector<Edge_weight>>>::const_iterator entries_begin,
        std::vector<std::pair<Frontier, std::vector<Edge_weight>>>::const_iterator entries_end) {
    nodes.push_back(TreeNode());
    TreeNode root = nodes.back();
    size_t root_idx = nodes.size() - 1;
    size_t cache_idx = 0;
    while(entries_begin != entries_end) {
        auto const& frontier = (*entries_begin).first;
        size_t cur = root_idx;
        for(size_t idx = 0; idx < frontier.size(); idx++) {
            size_t type = std::max(frontier[idx], frontier_index_t(252)) - 252;
            assert(type >= 0 && type < 4);
            if(nodes[cur].children[type] == size_t(-1)) {
                nodes[cur].children[type] = nodes.size();
                nodes.push_back(TreeNode());
            }
            cur = nodes[cur].children[type];
        }
        nodes[cur].cache_idx.push_back(cache_idx++);
        ++entries_begin;
    }
    return root;
}

} //namespace fpc