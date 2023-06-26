#include "treewidth_search.h"

namespace fpc {

TreewidthSearch::TreewidthSearch(Graph& input, std::vector<std::pair<Edge, std::vector<vertex_t>>> path_decomposition) 
    :   nthreads_(0),
        graph_(input),
        max_length_(graph_.max_length()),
        is_all_pair_(graph_.is_all_pair()),
        terminals_(graph_.terminals()),
        path_decomposition_(path_decomposition),
        remaining_edges_after_this_(path_decomposition_.size()),
        bag_local_idx_map_(
            path_decomposition_.size(),
            std::vector<frontier_index_t>(graph_.adjacency_.size(), invalid_index_)),
        result_(max_length_ + 1, 0),
        thread_local_result_(
            nthreads_,
            std::vector<Edge_weight>(max_length_ + 1, 0)
        ),
        cache_(
            path_decomposition_.size(),
            std::unordered_map<Frontier, std::vector<Edge_weight>, vec_hash>()
        ),
        pos_hits_(nthreads_, 0),
        neg_hits_(nthreads_, 0),
        edges_(nthreads_, 0),
        propagations_(nthreads_, 0) {
    
    std::vector<char> seen(graph_.adjacency_.size(), false);
    std::vector<vertex_t> last_remaining;
    std::vector<frontier_index_t> last_idx;
    if(!is_all_pair_) {
        last_remaining = {graph_.neighbors(terminals_[0]).size(), graph_.neighbors(terminals_[1]).size()};
        last_idx = std::vector<frontier_index_t>(graph_.adjacency_.size(), invalid_index_);
        last_idx[terminals_[0]] = start_index_;
        last_idx[terminals_[1]] = goal_index_;
        seen[terminals_[0]] = true;
        seen[terminals_[1]] = true;
    }
    for(size_t bag_idx = 0; bag_idx < path_decomposition_.size(); bag_idx++) {
        auto &bag = path_decomposition_[bag_idx].second;
        auto &idx = bag_local_idx_map_[bag_idx];
        auto &remaining = remaining_edges_after_this_[bag_idx];
        if(!is_all_pair_) {
            idx[terminals_[0]] = start_index_;
            idx[terminals_[1]] = goal_index_;
            remaining.resize(bag.size() + 2);
        } else {
            remaining.resize(bag.size());
        }
        frontier_index_t cur_idx = is_all_pair_?0:2;
        for(auto v : bag) {
            if(is_all_pair_ || (terminals_[0] != v && terminals_[1] != v)) {
                assert(cur_idx != invalid_index_);
                idx[v] = cur_idx++;
            }
            if(seen[v]) {
                frontier_index_t new_idx = idx[v];
                frontier_index_t old_idx = last_idx[v];
                remaining[new_idx] = last_remaining[old_idx];
            } else {
                seen[v] = true;
                remaining[idx[v]] = graph_.neighbors(v).size();
            }
        }
        Edge edge = path_decomposition_[bag_idx].first;
        remaining[idx[edge.first]]--;
        assert(remaining[idx[edge.first]] >= 0);
        remaining[idx[edge.second]]--;
        assert(remaining[idx[edge.second]] >= 0);
        last_remaining = remaining;
        last_idx = idx;
    }
    Frontier initial_frontier;
    if(!is_all_pair_) {
        initial_frontier = {no_edge_index_, no_edge_index_};
    }
    for(auto _ : path_decomposition_[0].second) {
        initial_frontier.push_back(no_edge_index_);
    }
    // cached vectors are {offset, results}
    std::vector<Edge_weight> initial_result = {0, 1};
    cache_[0][initial_frontier] = initial_result;
    includeSolutions(initial_frontier, 0, initial_result);
}

std::vector<Edge_weight> TreewidthSearch::search() {
    size_t thread_id = 0;
    for(size_t bag_idx = 0; bag_idx < path_decomposition_.size(); bag_idx++) {
        for(size_t bucket = 0; bucket < cache_[bag_idx].bucket_count(); bucket++) {
            for(auto task_it = cache_[bag_idx].begin(bucket); task_it != cache_[bag_idx].end(bucket); ++task_it) {
                auto const& old_frontier = task_it->first;
                auto const& result = task_it->second;
                Edge_length budget = max_length_ - result[0];
                for(auto [takeable, skippable] : {std::make_pair(true, false), std::make_pair(false, true)}) {
                    Frontier new_frontier = old_frontier;
                    frontier_index_t new_idx = bag_idx;
                    std::vector<Edge_weight> new_result = result;
                    // propagate while only one of the two is possible
                    while(takeable ^ skippable) {
                        if(takeable) {
                            take(new_frontier, new_idx);
                            new_result[0]++;
                        } else {
                            skip(new_frontier, new_idx);
                        }
                        new_idx++;
                        takeable = canTake(new_frontier, new_idx);
                        skippable = canSkip(new_frontier, new_idx);
                        includeSolutions(new_frontier, new_idx, new_result);
                    }
                    // both are possible, so we have a new decision edge
                    // put it into the cache
                    if(takeable && skippable) {
                        auto ins = cache_[new_idx].insert(
                            std::make_pair(new_frontier, new_result)
                        );
                        if(!ins.second) {
                            pos_hits_[thread_id]++;
                            // there is already an element with that key
                            // instead increase the partial result for that key
                            auto old_offset = ins.first->second[0];
                            auto new_offset = new_result[0];
                            if(old_offset > new_offset) {
                                // we reached this state with a smaller offset
                                // add the paths we currently cached to the new result
                                if(new_result.size() + new_offset < ins.first->second.size() + old_offset) {
                                    new_result.resize(ins.first->second.size() + old_offset - new_offset);
                                }
                                for(Edge_length res_length = 1; res_length < ins.first->second.size(); res_length++) {
                                    new_result[old_offset - new_offset + res_length] += ins.first->second[res_length];
                                }
                                ins.first->second = new_result;
                            } else {
                                if(ins.first->second.size() + old_offset < new_result.size() + new_offset) {
                                    ins.first->second.resize(new_result.size() + new_offset - old_offset);
                                }
                                for(Edge_length res_length = 1; res_length < new_result.size(); res_length++) {
                                    ins.first->second[new_offset - old_offset + res_length] += new_result[res_length];
                                }
                            }
                        } else {
                            neg_hits_[thread_id]++;
                        }
                    }
                }
            }
        }
        cache_[bag_idx] = {};
        // std::cerr << bag_idx << std::endl;
        // print_stats();
    }
    for(Edge_length length = 0; length <= max_length_; length++) {
        for(size_t id = 0; id < nthreads_; id++) {
            result_[length] += thread_local_result_[id][length];
        }
    }
    return result_;
}

void TreewidthSearch::includeSolutions(Frontier const& frontier, size_t bag_idx, std::vector<Edge_weight> const& partial_result) {
    assert(partial_result[0] + 1 <= max_length_);
    Edge edge = path_decomposition_[bag_idx].first;
    auto v_idx = bag_local_idx_map_[bag_idx][edge.first];
    auto w_idx = bag_local_idx_map_[bag_idx][edge.second];
    // must not be closing a path to a loop
    if(frontier[v_idx] == w_idx) {
        assert(frontier[w_idx] == v_idx);
        return;
    }
    assert(frontier[w_idx] != v_idx);

    size_t thread_id = 0;
    size_t number_paths = 0;
    for(frontier_index_t idx = 0; idx < frontier.size(); idx++) {
        if(frontier[idx] != no_edge_index_ && frontier[idx] != two_edge_index_) {
            number_paths++;
        }
    }
    assert(number_paths % 2 == 0);
    number_paths /= 2;
    if(is_all_pair_) {
        // no way this leads to a complete path
        if(number_paths > 2) {
            return;
        }
        // this can lead to a path if the first path and the second path connect
        if(number_paths == 2) {
            // cannot connect if they arent ends of paths
            if( frontier[v_idx] == no_edge_index_ 
                || frontier[v_idx] == two_edge_index_ 
                || frontier[w_idx] == no_edge_index_ 
                || frontier[w_idx] == two_edge_index_) {
                return;
            }
            // we connected them
            // add the partial result
            for(Edge_length res_length = 1; res_length < partial_result.size() && partial_result[0] + res_length <= max_length_; res_length++) {
                // current offset + 1 for the additional edge
                thread_local_result_[thread_id][partial_result[0] + res_length] += partial_result[res_length];
            }
            return;
        }
        // this can lead to a path if we continue the existing path
        if(number_paths == 1) {
            // cannot connect if neither is the end of the path
            if(     (frontier[v_idx] == no_edge_index_ 
                    || frontier[v_idx] == two_edge_index_)
                &&  (frontier[w_idx] == no_edge_index_ 
                    || frontier[w_idx] == two_edge_index_)) {
                return;
            }
            // we connected them
            // add the partial result
            for(Edge_length res_length = 1; res_length < partial_result.size() && partial_result[0] + res_length <= max_length_; res_length++) {
                // current offset + 1 for the additional edge
                thread_local_result_[thread_id][partial_result[0] + res_length] += partial_result[res_length];
            }
            return;
        }
        assert(number_paths == 0);
        assert(partial_result[0] == 0);
        thread_local_result_[thread_id][1] += 1;
        return;
    } else {
        // no way this leads to a complete path
        if(number_paths > 2) {
            return;
        }
        // this can lead to a path if the first path and the second path connect
        if(number_paths == 2) {
            // leads to a path if start and goal are both part of the path
            if(     (frontier[start_index_] != w_idx 
                    && frontier[start_index_] != v_idx)
                ||  (frontier[goal_index_] != w_idx 
                    && frontier[goal_index_] != v_idx)) {
                return;
            }
            assert(frontier[goal_index_] != frontier[start_index_]);
            // add the partial result
            for(Edge_length res_length = 1; res_length < partial_result.size() && partial_result[0] + res_length <= max_length_; res_length++) {
                // current offset + 1 for the additional edge
                thread_local_result_[thread_id][partial_result[0] + res_length] += partial_result[res_length];
            }
            return;
        }
        // this can lead to a path if we continue the existing path
        if(number_paths == 1) {
            // can leads to a path only if start or goal are part of the path
            if(frontier[start_index_] != w_idx 
                && frontier[start_index_] != v_idx
                && frontier[goal_index_] != w_idx 
                && frontier[goal_index_] != v_idx) {
                return;
            }
            // now the other index must be the missing of goal and end
            if(v_idx != start_index_
                && v_idx != goal_index_
                && w_idx != start_index_
                && w_idx != goal_index_) {
                return;
            }
            // we connected them
            // add the partial result
            for(Edge_length res_length = 1; res_length < partial_result.size() && partial_result[0] + res_length <= max_length_; res_length++) {
                // current offset + 1 for the additional edge
                thread_local_result_[thread_id][partial_result[0] + res_length] += partial_result[res_length];
            }
            return;
        }
        // we should already have covered the edges between start and goal during preprocessing
        assert(false);
    }
}

bool TreewidthSearch::canTake(Frontier& frontier, size_t bag_idx, std::vector<Edge_weight> const& partial_result) {
    assert(partial_result[0] + 1 <= max_length_);
    Edge edge = path_decomposition_[bag_idx].first;
    auto v_idx = bag_local_idx_map_[bag_idx][edge.first];
    auto w_idx = bag_local_idx_map_[bag_idx][edge.second];
    // must not be closing a path to a loop
    if(frontier[v_idx] == w_idx) {
        assert(frontier[w_idx] == v_idx);
        return false;
    }
    assert(frontier[w_idx] != v_idx);
    if(!is_all_pair_) {
        // start and goal are already connected
        if(frontier[start_idx_] == goal_idx_) {
            assert(frontier[goal_idx_] == start_idx_);
            return false;
        }
        // make sure start goal have only one outgoing edge
        if((v_idx == start_idx_ || v_idx == goal_idx_) && frontier[v_idx] != no_edge_index_) {
            return false;
        }
        if((w_idx == start_idx_ || w_idx == goal_idx_) && frontier[w_idx] != no_edge_index_) {
            return false;
        }

    }
    // make sure we took less than two edges at both ends of the new edge
    if(frontier[v_idx] == two_edge_index_ || frontier[w_idx] == two_edge_index_) {
        return false;
    }
    // length based pruning
    auto old_v = frontier[v_idx];
    auto old_w = frontier[w_idx];
    std::vector<std::pair<frontier_index_t>> restore({std::make_pair(v_idx, old_v), std::make_pair(w_idx, old_w)});
    if(old_v == no_edge_index_ && old_w == no_edge_idx_) {
        frontier[v_idx] = w_idx;
        frontier[w_idx] = v_idx;
    } else if(old_v != no_edge_index_ && old_w == no_edge_idx_) {
        restore.push_back(std::make_pair(old_v, frontier[old_v]));
        frontier[old_v] = w_idx;
        frontier[v_idx] = two_edge_index_;
        frontier[w_idx] = old_v;
    } else if(old_v == no_edge_index_ && old_w != no_edge_idx_) {
        restore.push_back(std::make_pair(old_w, frontier[old_w]));
        frontier[old_w] = v_idx;
        frontier[w_idx] = two_edge_index_;
        frontier[v_idx] = w_idx;
    } else {
        restore.push_back(std::make_pair(old_v, frontier[old_v]));
        restore.push_back(std::make_pair(old_w, frontier[old_w]));
        frontier[old_w] = v_idx;
        frontier[old_v] = w_idx;
        frontier[w_idx] = two_edge_index_;
        frontier[v_idx] = two_edge_index_;
    }
    size_t number_paths = 0;
    for(frontier_index_t idx = 0; idx < frontier.size(); idx++) {
        if(frontier[idx] != no_edge_index_ && frontier[idx] != two_edge_index_) {
            number_paths++;
        }
    }
    assert(number_paths % 2 == 0);
    number_paths /= 2;
    for(auto [idx, rest] : restore) {
        frontier[idx] = rest;
    }
    return number_paths + partial_result[0] <= max_length_;
}

bool TreewidthSearch::canSkip(Frontier& frontier, size_t bag_idx, std::vector<Edge_weight> const& partial_result) {
    assert(partial_result[0] + 1 <= max_length_);
    Edge edge = path_decomposition_[bag_idx].first;
    auto v_idx = bag_local_idx_map_[bag_idx][edge.first];
    auto w_idx = bag_local_idx_map_[bag_idx][edge.second];
    auto v_remaining = remaining_edges_after_this_[bag_idx][v_idx];
    auto w_remaining = remaining_edges_after_this_[bag_idx][w_idx];
    if(!is_all_pair_) {
        // start and goal are already connected
        if(frontier[start_idx_] == goal_idx_) {
            assert(frontier[goal_idx_] == start_idx_);
            return false;
        }
        // make sure start goal have at least one outgoing edge
        if((v_idx == start_idx_ || v_idx == goal_idx_) && frontier[v_idx] == no_edge_index_ && v_remaining == 0) {
            return false;
        }
        if((w_idx == start_idx_ || w_idx == goal_idx_) && frontier[w_idx] == no_edge_index_ && w_remaining == 0) {
            return false;
        }
    }
    // make sure we took exactly two or zero edges at both ends of the new edge
    if(v_remaining == 0 && frontier[v_idx] != two_edge_index_ && frontier[v_idx] != no_edge_index_) {
        return false;
    }
    if(w_remaining == 0 && frontier[w_idx] != two_edge_index_ && frontier[w_idx] != no_edge_index_) {
        return false;
    }
    // length based pruning
    size_t number_paths = 0;
    for(frontier_index_t idx = 0; idx < frontier.size(); idx++) {
        if(frontier[idx] != no_edge_index_ && frontier[idx] != two_edge_index_) {
            number_paths++;
        }
    }
    assert(number_paths % 2 == 0);
    number_paths /= 2;
    // + 1 because to connect number_paths, we need number_paths - 1 edges at least
    // but number_paths may be zero and we do not want to underflow (if thats what its called)
    return number_paths + partial_result[0] <= max_length_ + 1;
}

void TreewidthSearch::take(Frontier& frontier, size_t bag_idx) {
    
}

void TreewidthSearch::skip(Frontier& frontier, size_t bag_idx) {
    
}

}