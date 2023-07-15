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

#include "treewidth_search.h"
#include <algorithm>
namespace fpc {

TreewidthSearch::TreewidthSearch(Graph& input, AnnotatedDecomposition decomposition, size_t nthreads) 
    :   nthreads_(nthreads),
        graph_(input),
        max_length_(graph_.max_length()),
        is_all_pair_(graph_.is_all_pair()),
        terminals_(graph_.terminals()),
        decomposition_(decomposition),
        remaining_edges_after_this_(decomposition_.size()),
        bag_local_idx_map_(
            decomposition_.size(),
            std::vector<frontier_index_t>(graph_.adjacency_.size(), invalid_index_)),
        bag_local_vertex_map_(decomposition_.size()),
        bag_local_distance_(decomposition_.size()),
        result_(max_length_ + 1, 0),
        thread_local_result_(
            nthreads_,
            std::vector<Edge_weight>(max_length_ + 1, 0)
        ),
        cache_(decomposition_.size()),
        pos_hits_(nthreads_, 0),
        neg_hits_(nthreads_, 0),
        edges_(nthreads_, 0),
        propagations_(nthreads_, 0),
        merges_(nthreads_, 0),
        unsuccessful_merges_(nthreads_, 0) {
    
    omp_set_num_threads(nthreads_);

    if(!is_all_pair_) {
        for(auto &node : decomposition_) {
            auto &bag = node.bag;
            if(std::find(bag.begin(), bag.end(), terminals_[0]) == bag.end()) {
                bag.push_back(terminals_[0]);
            }
            if(std::find(bag.begin(), bag.end(), terminals_[1]) == bag.end()) {
                bag.push_back(terminals_[1]);
            }
        }
    }

    std::map<size_t, size_t> idx_remap;
    size_t invalid_td_idx = -1;
    idx_remap[invalid_td_idx] = invalid_td_idx;
    size_t cur = decomposition_.size() - 1;
    size_t root = 0;
    while(root < decomposition_.size() && decomposition_[root].parent != invalid_td_idx) {
        root++;
    }
    assert(root < decomposition_.size());
    std::vector<size_t> remap_stack = {root};
    AnnotatedDecomposition ordered;
    while(!remap_stack.empty()) {
        size_t top = remap_stack.back();
        idx_remap[top] = cur--;
        remap_stack.pop_back();
        auto &top_node = decomposition_[top];
        switch (top_node.type) {
        case NodeType::LEAF:
            ordered.push_back(top_node);
            break;
        case NodeType::PATH_LIKE:
            assert(top > 0);
            ordered.push_back(top_node);
            remap_stack.push_back(top_node.children.first);
            break;
        case NodeType::JOIN:
            ordered.push_back(top_node);
            remap_stack.push_back(top_node.children.first);
            remap_stack.push_back(top_node.children.second);
            break;
        }
    }
    decomposition_ = AnnotatedDecomposition(ordered.rbegin(), ordered.rend());
    for(auto &node : decomposition_) {
        node.parent = idx_remap[node.parent];
        node.children.first = idx_remap[node.children.first];
        node.children.second = idx_remap[node.children.second];
    }

    std::set<vertex_t> empty;
    std::vector<vertex_t> all;
    for(vertex_t i = 0; i < graph_.adjacency_.size(); i++) {
        all.push_back(i);
    }
    Graph cur_graph = graph_.subgraph(all);;
    for(size_t bag_idx = 0; bag_idx < decomposition_.size(); bag_idx++) {
        auto &node = decomposition_[bag_idx];
        auto &idx = bag_local_idx_map_[bag_idx];
        auto &vertex = bag_local_vertex_map_[bag_idx];
        auto &remaining = remaining_edges_after_this_[bag_idx];
        switch(node.type) {
        case NodeType::LEAF:
            cur_graph = graph_.subgraph(all);
            break;
        case NodeType::PATH_LIKE:
            break;
        case NodeType::JOIN:
            cur_graph = graph_.subgraph(all);
            std::vector<size_t> stack = {bag_idx};
            while(!stack.empty()) {
                size_t top = stack.back();
                stack.pop_back();
                auto &top_node = decomposition_[top];
                switch (top_node.type) {
                case NodeType::LEAF:
                    cur_graph.remove_edge(top_node.edge);
                    break;
                case NodeType::PATH_LIKE:
                    cur_graph.remove_edge(top_node.edge);
                    assert(top > 0);
                    stack.push_back(top - 1);
                    break;
                case NodeType::JOIN:
                    stack.push_back(top_node.children.first);
                    stack.push_back(top_node.children.second);
                    break;
                }
            }
            break;
        }
        for(size_t i = 0; i < node.bag.size(); i++) {
            if(!is_all_pair_ && (node.bag[i] == terminals_[0] || node.bag[i] == terminals_[1])) {
                continue;
            }
            if(node.type != NodeType::JOIN && cur_graph.neighbors(node.bag[i]).size() == 0) {
                    std::swap(node.bag[i], node.bag.back());
                    node.bag.pop_back();
                    i--;
            } else if(cur_graph.neighbors(node.bag[i]).size() == graph_.neighbors(node.bag[i]).size() 
                && node.bag[i] != node.edge.first && node.bag[i] != node.edge.second) {
                std::swap(node.bag[i], node.bag.back());
                node.bag.pop_back();
                i--;
            } else if(node.type == NodeType::JOIN) {
                auto v = node.bag[i];
                auto left_child = node.children.first;
                auto right_child = node.children.second;
                auto left_idx = bag_local_idx_map_[left_child][v];
                auto right_idx = bag_local_idx_map_[right_child][v];
                if(right_idx != invalid_index_ && remaining_edges_after_this_[right_child][right_idx] == 0) {
                    assert(left_idx == invalid_index_);
                    std::swap(node.bag[i], node.bag.back());
                    node.bag.pop_back();
                    i--;
                } else if(left_idx != invalid_index_ && remaining_edges_after_this_[left_child][left_idx] == 0) {
                    assert(right_idx == invalid_index_);
                    std::swap(node.bag[i], node.bag.back());
                    node.bag.pop_back();
                    i--;
                }
            }
        }
        std::sort(node.bag.begin(), node.bag.end());
        bag_local_distance_[bag_idx].resize(node.bag.size());
        remaining.resize(node.bag.size());
        if(node.type != NodeType::JOIN) {
            cur_graph.remove_edge(node.edge);
        }
        frontier_index_t cur_idx = 0;
        for(auto v : node.bag) {
            assert(cur_idx != invalid_index_);
            idx[v] = cur_idx++;
            vertex.push_back(v);
            remaining[idx[v]] = cur_graph.neighbors(v).size();
        }
        for(auto v : node.bag) {
            std::vector<Edge_length> distance(graph_.adjacency_.size(), invalid_distance_);
            cur_graph.dijkstra(v, distance, empty);
            std::vector<Edge_length> local_distance(node.bag.size(), invalid_distance_);
            for(auto other : node.bag) {
                if(other != v) {
                    local_distance[idx[other]] = distance[other];
                }
            }
            bag_local_distance_[bag_idx][idx[v]] = local_distance;
        }
    }
    for(vertex_t v = 0; v < graph_.adjacency_.size(); v++) {
        if(cur_graph.neighbors(v).size() != 0) {
            std::cerr << v << std::endl;
        }
        assert(cur_graph.neighbors(v).size() == 0);
    }
    for(size_t bag_idx = 0; bag_idx < decomposition_.size(); bag_idx++) {
        auto &node = decomposition_[bag_idx];
        if(node.type == NodeType::LEAF) {
            Frontier initial_frontier(decomposition_[bag_idx].bag.size(), no_edge_index_);
            if(!is_all_pair_) {
                assert(bag_local_idx_map_[bag_idx][terminals_[0]] != invalid_index_);
                assert(bag_local_idx_map_[bag_idx][terminals_[1]] != invalid_index_);
                initial_frontier[bag_local_idx_map_[bag_idx][terminals_[0]]] = invalid_index_;
                initial_frontier[bag_local_idx_map_[bag_idx][terminals_[1]]] = invalid_index_;
            }
            // cached vectors are {offset, results}
            std::vector<Edge_weight> initial_result = {0, 1};
            cache_[bag_idx].first[initial_frontier] = initial_result;
            includeSolutions(initial_frontier, bag_idx, initial_result);
        }
    }
}

std::vector<Edge_weight> TreewidthSearch::search() {
    for(size_t bag_idx = 0; bag_idx < decomposition_.size(); bag_idx++) {
        std::cerr << "\r" << bag_idx << " / " << decomposition_.size() - 1 << " of type ";
        if(decomposition_[bag_idx].type == NodeType::LEAF) {
            std::cerr << "leaf, with " << cache_[bag_idx].first.size() << " entries.";
        } else if(decomposition_[bag_idx].type == NodeType::PATH_LIKE) {
            std::cerr << "path, with " << cache_[bag_idx].first.size() << " entries.";
        } else if(decomposition_[bag_idx].type == NodeType::JOIN) {
            std::cerr << "join, with (" << cache_[bag_idx].first.size() << "," << cache_[bag_idx].second.size() << ") entries.";
        }
        if(decomposition_[bag_idx].type == LEAF || decomposition_[bag_idx].type == PATH_LIKE) {
            // PATH_LIKE/LEAF
            #pragma omp parallel for default(shared) shared(edges_) shared(propagations_) shared(bag_idx) shared(max_length_) shared(cache_) shared(decomposition_) shared(thread_local_result_) shared(pos_hits_) shared(neg_hits_) shared(bag_local_idx_map_) shared(bag_local_vertex_map_)
            for(size_t bucket = 0; bucket < cache_[bag_idx].first.bucket_count(); bucket++) {
                size_t thread_id = omp_get_thread_num();
                for(auto task_it = cache_[bag_idx].first.begin(bucket); task_it != cache_[bag_idx].first.end(bucket); ++task_it) {
                    auto const& old_frontier = task_it->first;
                    auto const& result = task_it->second;
                    for(auto pp : {std::make_pair(true, false), std::make_pair(false, true)}) {
                        bool takeable = pp.first;
                        bool skippable = pp.second;
                        Frontier new_frontier = old_frontier;
                        size_t last_idx = -1;
                        size_t new_idx = bag_idx;
                        std::vector<Edge_weight> new_result = result;
                        edges_[thread_id]++;
                        propagations_[thread_id]--;
                        // propagate while only one of the two is possible
                        while((takeable ^ skippable) && new_idx + 1 < decomposition_.size() && decomposition_[new_idx].type != JOIN) {
                            propagations_[thread_id]++;
                            // Edge edge = decomposition_[new_idx].edge;
                            // auto v_idx = bag_local_idx_map_[new_idx][edge.first];
                            // auto w_idx = bag_local_idx_map_[new_idx][edge.second];
                            // std::cerr << "Edge: (" << size_t(v_idx) << "," << size_t(w_idx) << ") Idx:" << size_t(new_idx) << std::endl;
                            // for(auto idx : new_frontier) {
                            //     std::cerr << size_t(idx) << " ";
                            // }
                            // std::cerr << std::endl;
                            // std::cerr << "Remaining: ";
                            // for(auto i = 0; i < new_frontier.size(); i++) {
                            //     std::cerr << size_t(remaining_edges_after_this_[new_idx][i]) << " ";
                            // }
                            // std::cerr << std::endl;
                            // std::cerr << takeable << " " << skippable << std::endl;
                            if(takeable) {
                                take(new_frontier, new_idx);
                                new_result[0]++;
                            } else {
                                skip(new_frontier, new_idx);
                            }
                            last_idx = new_idx;
                            new_idx = decomposition_[new_idx].parent;
                            if(new_idx < decomposition_.size() && decomposition_[new_idx].type != JOIN) {
                                // Edge edge = decomposition_[new_idx].first;
                                // auto v_idx = bag_local_idx_map_[new_idx][edge.first];
                                // auto w_idx = bag_local_idx_map_[new_idx][edge.second];
                                // std::cerr << "Edge: (" << size_t(v_idx) << "," << size_t(w_idx) << ") Idx:" << size_t(new_idx) << std::endl;
                                // for(auto idx : new_frontier) {
                                //     std::cerr << size_t(idx) << " ";
                                // }
                                // std::cerr << std::endl;
                                // std::cerr << "Remaining: ";
                                // for(auto i = 0; i < new_frontier.size(); i++) {
                                //     std::cerr << size_t(remaining_edges_after_this_[new_idx][i]) << " ";
                                // }
                                // std::cerr << std::endl;
                                // std::cerr << "Result: ";
                                // for(auto i = 0; i < new_result.size(); i++) {
                                //     std::cerr << new_result[i] << " ";
                                // }
                                // std::cerr << std::endl;
                                takeable = canTake(new_frontier, new_idx, new_result);
                                skippable = canSkip(new_frontier, new_idx, new_result);
                                includeSolutions(new_frontier, new_idx, new_result);
                                // std::cerr << takeable << " " << skippable << std::endl;
                            }
                        }
                        // both are possible, so we have a new decision edge
                        // put it into the cache
                        if(new_idx < decomposition_.size() 
                            && (
                                    (takeable && skippable)
                                ||  (decomposition_[new_idx].type == JOIN && last_idx == decomposition_[new_idx].children.first))) {
                            #pragma omp critical 
                            {
                                auto ins = cache_[new_idx].first.insert(
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
                        } else if(new_idx < decomposition_.size() 
                            && decomposition_[new_idx].type == JOIN) {
                            assert(last_idx == decomposition_[new_idx].children.second);
                            #pragma omp critical 
                            {
                                auto ins = cache_[new_idx].second.insert(
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
            }
        } else {
            // JOIN
            // FIXME: check if its better to take either bag as the outer one
            #pragma omp parallel for default(shared) shared(merges_) shared(unsuccessful_merges_) shared(edges_) shared(propagations_) shared(bag_idx) shared(max_length_) shared(cache_) shared(decomposition_) shared(thread_local_result_) shared(pos_hits_) shared(neg_hits_) shared(bag_local_idx_map_) shared(bag_local_vertex_map_)
            for(size_t bucket = 0; bucket < cache_[bag_idx].first.bucket_count(); bucket++) {
                size_t thread_id = omp_get_thread_num();
                for(auto task_it = cache_[bag_idx].first.begin(bucket); task_it != cache_[bag_idx].first.end(bucket); ++task_it) {
                    // FIXME: can we just modify the object somehow?
                    for(auto const&[right_frontier, right_result] : cache_[bag_idx].second) {
                        auto left_frontier = task_it->first;
                        auto left_result = task_it->second;
                        // for(auto idx : left_frontier) {
                        //     std::cerr << size_t(idx) << " ";
                        // }
                        // std::cerr << std::endl;
                        // for(auto idx : right_frontier) {
                        //     std::cerr << size_t(idx) << " ";
                        // }
                        // std::cerr << std::endl;
                        // std::cerr << "Remaining: ";
                        // for(auto i = 0; i < left_frontier.size(); i++) {
                        //     std::cerr << size_t(remaining_edges_after_this_[bag_idx][i]) << " ";
                        // }
                        // std::cerr << std::endl;
                        merges_[thread_id]++;
                        if(!merge(left_frontier, right_frontier, bag_idx, left_result, right_result)) {
                            // std::cerr << "Not merged" << std::endl;
                            unsuccessful_merges_[thread_id]++;
                            continue;
                        }
                        // std::cerr << "Merged" << std::endl;
                        size_t last_idx = bag_idx;
                        size_t new_idx = decomposition_[bag_idx].parent;
                        assert(new_idx != size_t(-1));
                        bool takeable = false;
                        bool skippable = false;
                        if(decomposition_[new_idx].type != JOIN) {
                            // Edge edge = decomposition_[new_idx].edge;
                            // auto v_idx = bag_local_idx_map_[new_idx][edge.first];
                            // auto w_idx = bag_local_idx_map_[new_idx][edge.second];
                            // std::cerr << "Edge: (" << size_t(v_idx) << "," << size_t(w_idx) << ") Idx:" << size_t(new_idx) << std::endl;
                            // for(auto idx : left_frontier) {
                            //     std::cerr << size_t(idx) << " ";
                            // }
                            // std::cerr << std::endl;
                            // std::cerr << "Remaining: ";
                            // for(auto i = 0; i < left_frontier.size(); i++) {
                            //     std::cerr << size_t(remaining_edges_after_this_[new_idx][i]) << " ";
                            // }
                            // std::cerr << std::endl;
                            // std::cerr << "Result: ";
                            // for(auto i = 0; i < left_result.size(); i++) {
                            //     std::cerr << left_result[i] << " ";
                            // }
                            // std::cerr << std::endl;
                            takeable = canTake(left_frontier, new_idx, left_result);
                            skippable = canSkip(left_frontier, new_idx, left_result);
                            includeSolutions(left_frontier, new_idx, left_result);
                            // std::cerr << takeable << " " << skippable << std::endl;
                            if(takeable ^ skippable) {
                                edges_[thread_id]++;
                                propagations_[thread_id]--;
                            }
                        }
                        // propagate while only one of the two is possible
                        while((takeable ^ skippable) && new_idx + 1 < decomposition_.size() && decomposition_[new_idx].type != JOIN) {
                            propagations_[thread_id]++;
                            // Edge edge = decomposition_[new_idx].edge;
                            // auto v_idx = bag_local_idx_map_[new_idx][edge.first];
                            // auto w_idx = bag_local_idx_map_[new_idx][edge.second];
                            // std::cerr << "Edge: (" << size_t(v_idx) << "," << size_t(w_idx) << ") Idx:" << size_t(new_idx) << std::endl;
                            // for(auto idx : left_frontier) {
                            //     std::cerr << size_t(idx) << " ";
                            // }
                            // std::cerr << std::endl;
                            // std::cerr << "Remaining: ";
                            // for(auto i = 0; i < left_frontier.size(); i++) {
                            //     std::cerr << size_t(remaining_edges_after_this_[new_idx][i]) << " ";
                            // }
                            // std::cerr << std::endl;
                            // std::cerr << takeable << " " << skippable << std::endl;
                            if(takeable) {
                                take(left_frontier, new_idx);
                                left_result[0]++;
                            } else {
                                skip(left_frontier, new_idx);
                            }
                            last_idx = new_idx;
                            new_idx = decomposition_[new_idx].parent;
                            if(new_idx < decomposition_.size() && decomposition_[new_idx].type != JOIN) {
                                // Edge edge = decomposition_[new_idx].edge;
                                // auto v_idx = bag_local_idx_map_[new_idx][edge.first];
                                // auto w_idx = bag_local_idx_map_[new_idx][edge.second];
                                // std::cerr << "Edge: (" << size_t(v_idx) << "," << size_t(w_idx) << ") Idx:" << size_t(new_idx) << std::endl;
                                // for(auto idx : left_frontier) {
                                //     std::cerr << size_t(idx) << " ";
                                // }
                                // std::cerr << std::endl;
                                // std::cerr << "Remaining: ";
                                // for(auto i = 0; i < left_frontier.size(); i++) {
                                //     std::cerr << size_t(remaining_edges_after_this_[new_idx][i]) << " ";
                                // }
                                // std::cerr << std::endl;
                                // std::cerr << "Result: ";
                                // for(auto i = 0; i < left_result.size(); i++) {
                                //     std::cerr << left_result[i] << " ";
                                // }
                                // std::cerr << std::endl;
                                takeable = canTake(left_frontier, new_idx, left_result);
                                skippable = canSkip(left_frontier, new_idx, left_result);
                                includeSolutions(left_frontier, new_idx, left_result);
                                // std::cerr << takeable << " " << skippable << std::endl;
                            }
                        }
                        // both are possible, so we have a new decision edge
                        // put it into the cache
                        if(new_idx < decomposition_.size() 
                            && (
                                    (takeable && skippable)
                                ||  (decomposition_[new_idx].type == JOIN && last_idx == decomposition_[new_idx].children.first))) {
                            #pragma omp critical 
                            {
                                auto ins = cache_[new_idx].first.insert(
                                    std::make_pair(left_frontier, left_result)
                                );
                                if(!ins.second) {
                                    pos_hits_[thread_id]++;
                                    // there is already an element with that key
                                    // instead increase the partial result for that key
                                    auto old_offset = ins.first->second[0];
                                    auto new_offset = left_result[0];
                                    if(old_offset > new_offset) {
                                        // we reached this state with a smaller offset
                                        // add the paths we currently cached to the new result
                                        if(left_result.size() + new_offset < ins.first->second.size() + old_offset) {
                                            left_result.resize(ins.first->second.size() + old_offset - new_offset);
                                        }
                                        for(Edge_length res_length = 1; res_length < ins.first->second.size(); res_length++) {
                                            left_result[old_offset - new_offset + res_length] += ins.first->second[res_length];
                                        }
                                        ins.first->second = left_result;
                                    } else {
                                        if(ins.first->second.size() + old_offset < left_result.size() + new_offset) {
                                            ins.first->second.resize(left_result.size() + new_offset - old_offset);
                                        }
                                        for(Edge_length res_length = 1; res_length < left_result.size(); res_length++) {
                                            ins.first->second[new_offset - old_offset + res_length] += left_result[res_length];
                                        }
                                    }
                                } else {
                                    neg_hits_[thread_id]++;
                                }
                            }
                        } else if(new_idx < decomposition_.size() 
                            && decomposition_[new_idx].type == JOIN) {
                            assert(last_idx == decomposition_[new_idx].children.second);
                            #pragma omp critical 
                            {
                                auto ins = cache_[new_idx].second.insert(
                                    std::make_pair(left_frontier, left_result)
                                );
                                if(!ins.second) {
                                    pos_hits_[thread_id]++;
                                    // there is already an element with that key
                                    // instead increase the partial result for that key
                                    auto old_offset = ins.first->second[0];
                                    auto new_offset = left_result[0];
                                    if(old_offset > new_offset) {
                                        // we reached this state with a smaller offset
                                        // add the paths we currently cached to the new result
                                        if(left_result.size() + new_offset < ins.first->second.size() + old_offset) {
                                            left_result.resize(ins.first->second.size() + old_offset - new_offset);
                                        }
                                        for(Edge_length res_length = 1; res_length < ins.first->second.size(); res_length++) {
                                            left_result[old_offset - new_offset + res_length] += ins.first->second[res_length];
                                        }
                                        ins.first->second = left_result;
                                    } else {
                                        if(ins.first->second.size() + old_offset < left_result.size() + new_offset) {
                                            ins.first->second.resize(left_result.size() + new_offset - old_offset);
                                        }
                                        for(Edge_length res_length = 1; res_length < left_result.size(); res_length++) {
                                            ins.first->second[new_offset - old_offset + res_length] += left_result[res_length];
                                        }
                                    }
                                } else {
                                    neg_hits_[thread_id]++;
                                }
                            }
                        }
                    }
                }
            }
        }

        cache_[bag_idx] = {};
    }
    for(Edge_length length = 0; length <= max_length_; length++) {
        for(size_t id = 0; id < nthreads_; id++) {
            result_[length] += thread_local_result_[id][length];
        }
    }
    std::cerr << std::endl;
    return result_;
}

void TreewidthSearch::includeSolutions(Frontier const& frontier, size_t bag_idx, std::vector<Edge_weight> const& partial_result) {
    assert(partial_result[0] + 1 <= max_length_);
    Edge edge = decomposition_[bag_idx].edge;
    auto v_idx = bag_local_idx_map_[bag_idx][edge.first];
    auto w_idx = bag_local_idx_map_[bag_idx][edge.second];
    // must not be closing a path to a loop
    if(frontier[v_idx] == w_idx) {
        assert(frontier[w_idx] == v_idx);
        return;
    }
    assert(frontier[w_idx] != v_idx);

    if(frontier[w_idx] == two_edge_index_ || frontier[v_idx] == two_edge_index_) {
        return;
    }

    size_t thread_id = omp_get_thread_num();
    size_t number_paths = 0;
    size_t number_cut_paths = 0;
    for(frontier_index_t idx = 0; idx < frontier.size(); idx++) {
        if(frontier[idx] != no_edge_index_ && frontier[idx] != two_edge_index_) {
            number_paths++;
        }
        if(frontier[idx] == invalid_index_) {
            number_cut_paths++;
            number_paths--;
        }
    }
    assert(number_paths % 2 == 0);
    number_paths /= 2;
    number_paths += number_cut_paths;
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
}

bool TreewidthSearch::canTake(Frontier& frontier, size_t bag_idx, std::vector<Edge_weight> const& partial_result) {
    assert(partial_result[0] + 1 <= max_length_);
    Edge edge = decomposition_[bag_idx].edge;
    auto v_idx = bag_local_idx_map_[bag_idx][edge.first];
    auto w_idx = bag_local_idx_map_[bag_idx][edge.second];
    auto v_remaining = remaining_edges_after_this_[bag_idx][v_idx];
    auto w_remaining = remaining_edges_after_this_[bag_idx][w_idx];

    // must not be closing a path to a loop
    if(frontier[v_idx] == w_idx) {
        assert(frontier[w_idx] == v_idx);
        return false;
    }

    int takeable_idx_v = std::max(frontier[v_idx], frontier_index_t(252)) - 252;
    if(v_remaining) {
        takeable_idx_v += 4;
    }
    int takeable_idx_w = std::max(frontier[w_idx], frontier_index_t(252)) - 252;
    if(w_remaining) {
        takeable_idx_w += 4;
    }

    // make sure we took less than two edges at both ends of the new edge
    if(!takeable_[takeable_idx_v][takeable_idx_w]) {
        return false;
    }

    assert(w_idx < frontier.size());
    assert(v_idx < frontier.size());
    assert(frontier[w_idx] != v_idx);
    // length based pruning
    auto old_v = frontier[v_idx];
    auto old_w = frontier[w_idx];
    std::vector<std::pair<frontier_index_t, frontier_index_t>> restore({std::make_pair(v_idx, old_v), std::make_pair(w_idx, old_w)});
    if(old_v == no_edge_index_) {
        if(old_w == no_edge_index_) {
            frontier[v_idx] = w_idx;
            frontier[w_idx] = v_idx;
        } else if(old_w == invalid_index_) {
            frontier[w_idx] = two_edge_index_;
            frontier[v_idx] = invalid_index_;
        } else {
            restore.push_back(std::make_pair(old_w, frontier[old_w]));
            frontier[old_w] = v_idx;
            frontier[w_idx] = two_edge_index_;
            frontier[v_idx] = old_w;
        }
    } else if(old_v == invalid_index_) {
        if(old_w == no_edge_index_) {
            frontier[v_idx] = two_edge_index_;
            frontier[w_idx] = invalid_index_;
        } else if(old_w == invalid_index_) {
            assert(false);
        } else {
            restore.push_back(std::make_pair(old_w, frontier[old_w]));
            frontier[old_w] = invalid_index_;
            frontier[w_idx] = two_edge_index_;
            frontier[v_idx] = two_edge_index_;
        }
    } else {
        restore.push_back(std::make_pair(old_v, frontier[old_v]));
        if(old_w == no_edge_index_) {
            frontier[old_v] = w_idx;
            frontier[v_idx] = two_edge_index_;
            frontier[w_idx] = old_v;
        } else if(old_w == invalid_index_) {
            frontier[old_v] = invalid_index_;
            frontier[v_idx] = two_edge_index_;
            frontier[w_idx] = two_edge_index_;
        } else {
            restore.push_back(std::make_pair(old_w, frontier[old_w]));
            frontier[old_w] = old_v;
            frontier[old_v] = old_w;
            frontier[w_idx] = two_edge_index_;
            frontier[v_idx] = two_edge_index_;
        }
    }
    if(v_remaining == 0) {
        assert(frontier[v_idx] != invalid_index_);
        assert(frontier[v_idx] != no_edge_index_);
        if(frontier[v_idx] != two_edge_index_) {
            frontier[frontier[v_idx]] = invalid_index_;
            frontier[v_idx] = two_edge_index_;
        }
    }
    if(w_remaining == 0) {
        assert(frontier[w_idx] != invalid_index_);
        assert(frontier[w_idx] != no_edge_index_);
        if(frontier[w_idx] != two_edge_index_) {
            frontier[frontier[w_idx]] = invalid_index_;
            frontier[w_idx] = two_edge_index_;
        }
    }
    std::vector<frontier_index_t> paths;
    std::vector<frontier_index_t> cut_paths;
    for(frontier_index_t idx = 0; idx < frontier.size(); idx++) {
        // std::cerr << size_t(frontier[idx]) << " ";
        if(frontier[idx] != no_edge_index_ && frontier[idx] != two_edge_index_) {
            if(frontier[idx] == invalid_index_) {
                cut_paths.push_back(idx);
            } else {
                paths.push_back(idx);
            }
        }
    }
    // std::cerr << std::endl;
    if(cut_paths.size() > 2) {
        for(auto [idx, rest] : restore) {
            frontier[idx] = rest;
        }
        return false;
    }
    assert(paths.size() % 2 == 0);

    for(auto [idx, rest] : restore) {
        frontier[idx] = rest;
    }
    // + 1 - 1 since were taking the edge
    if(paths.size()/2 + cut_paths.size() > 1) {
        if(paths.size()/2 + cut_paths.size() + partial_result[0] > max_length_) {
            return false;
        }
    }
     else if(paths.size()/2 + cut_paths.size() + partial_result[0]  + 1 > max_length_) {
        return false;
    }

    return distancePrune(frontier, paths, cut_paths, bag_idx, partial_result[0] + 1);
}

bool TreewidthSearch::canSkip(Frontier& frontier, size_t bag_idx, std::vector<Edge_weight> const& partial_result) {
    assert(partial_result[0] + 1 <= max_length_);
    Edge edge = decomposition_[bag_idx].edge;
    auto v_idx = bag_local_idx_map_[bag_idx][edge.first];
    auto w_idx = bag_local_idx_map_[bag_idx][edge.second];
    auto v_remaining = remaining_edges_after_this_[bag_idx][v_idx];
    auto w_remaining = remaining_edges_after_this_[bag_idx][w_idx];

    int skippable_idx_v = std::max(frontier[v_idx], frontier_index_t(252)) - 252;
    if(v_remaining) {
        skippable_idx_v += 4;
    }
    int skippable_idx_w = std::max(frontier[w_idx], frontier_index_t(252)) - 252;
    if(w_remaining) {
        skippable_idx_w += 4;
    }
    if(!skippable_[skippable_idx_v][skippable_idx_w]) {
        return false;
    }
    if(v_remaining == 0 && w_remaining == 0
        && frontier[v_idx] == w_idx) {
        return false;
    }
    // length based pruning
    std::vector<frontier_index_t> paths;
    std::vector<frontier_index_t> cut_paths;
    for(frontier_index_t idx = 0; idx < frontier.size(); idx++) {
        // std::cerr << size_t(frontier[idx]) << " ";
        if(frontier[idx] != no_edge_index_ && frontier[idx] != two_edge_index_) {
            if(frontier[idx] == invalid_index_) {
                cut_paths.push_back(idx);
            } else {
                paths.push_back(idx);
            }
        }
    }
    assert(cut_paths.size() <= 2);
    assert(paths.size() % 2 == 0);
    bool last_and_incomplete_v = v_remaining == 0 && frontier[v_idx] != two_edge_index_ && frontier[v_idx] != no_edge_index_;
    bool last_and_incomplete_w = w_remaining == 0 && frontier[w_idx] != two_edge_index_ && frontier[w_idx] != no_edge_index_;
    if(last_and_incomplete_v && last_and_incomplete_w) {
        if(cut_paths.size() > 0) {
            return false;
        }
    } else if(last_and_incomplete_v || last_and_incomplete_w) {
        if(cut_paths.size() > 1) {
            return false;
        }
    }

    return distancePrune(frontier, paths, cut_paths, bag_idx, partial_result[0]);
}

bool TreewidthSearch::distancePrune(
        Frontier& frontier, 
        std::vector<frontier_index_t> const& paths, 
        std::vector<frontier_index_t> const& cut_paths, 
        size_t bag_idx, 
        size_t offset) {
    // advanced length based pruning
    auto &distance = bag_local_distance_[bag_idx];
    if(cut_paths.size() == 2 && paths.size() == 0) {
        // we have to connect the cut paths
        return distance[cut_paths[0]][cut_paths[1]] + offset <= max_length_;
    }
    if(cut_paths.size() + paths.size() <= 1) {
        return offset + 1 <= max_length_;
    }

    // there is at least one other path that we have to incorporate
    // connecting the cut paths is not an option
    size_t min_dist = 0;
    for(auto v : cut_paths) {
        Edge_length cur_min_dist = invalid_distance_;
        for(auto other : paths) {
            cur_min_dist = std::min(cur_min_dist, distance[v][other]);
        }
        if(cur_min_dist == invalid_distance_) {
            // we cannot connect this cut path
            return false;
        }
        min_dist += cur_min_dist;
    }
    if(min_dist + offset > max_length_) {
        return false;
    }
    size_t path_min_dist = 0;
    Edge_length path_max_dist = 0;
    Edge_length path_snd_max_dist = 0;
    for(auto v : paths) {
        Edge_length cur_min_dist = invalid_distance_;
        for(auto cut_end : cut_paths) {
            cur_min_dist = std::min(cur_min_dist, distance[v][cut_end]);
        }
        for(auto other : paths) {
            if(frontier[v] != other) {
                cur_min_dist = std::min(cur_min_dist, distance[v][other]);
            }
        }
        if(cur_min_dist >= path_max_dist) {
            path_snd_max_dist = path_max_dist;
            path_max_dist = cur_min_dist;
        } else if(cur_min_dist > path_snd_max_dist) {
            path_snd_max_dist = cur_min_dist;
        }
        path_min_dist += cur_min_dist;
    }
    assert(path_min_dist >= path_max_dist + path_snd_max_dist);
    path_min_dist -= path_max_dist;
    path_min_dist -= path_snd_max_dist;
    path_min_dist /= 2;
    if(path_min_dist + min_dist + offset >  max_length_) {
        return false;
    }
    return true;
}

void TreewidthSearch::take(Frontier& frontier, size_t bag_idx) {
    Edge edge = decomposition_[bag_idx].edge;
    auto v_idx = bag_local_idx_map_[bag_idx][edge.first];
    auto w_idx = bag_local_idx_map_[bag_idx][edge.second];
    auto old_v = frontier[v_idx];
    auto old_w = frontier[w_idx];

    if(old_v == no_edge_index_) {
        if(old_w == no_edge_index_) {
            frontier[v_idx] = w_idx;
            frontier[w_idx] = v_idx;
        } else if(old_w == invalid_index_) {
            frontier[w_idx] = two_edge_index_;
            frontier[v_idx] = invalid_index_;
        } else {
            frontier[old_w] = v_idx;
            frontier[w_idx] = two_edge_index_;
            frontier[v_idx] = old_w;
        }
    } else if(old_v == invalid_index_) {
        if(old_w == no_edge_index_) {
            frontier[v_idx] = two_edge_index_;
            frontier[w_idx] = invalid_index_;
        } else if(old_w == invalid_index_) {
            assert(false);
        } else {
            frontier[old_w] = invalid_index_;
            frontier[w_idx] = two_edge_index_;
            frontier[v_idx] = two_edge_index_;
        }
    } else {
        if(old_w == no_edge_index_) {
            frontier[old_v] = w_idx;
            frontier[v_idx] = two_edge_index_;
            frontier[w_idx] = old_v;
        } else if(old_w == invalid_index_) {
            frontier[old_v] = invalid_index_;
            frontier[v_idx] = two_edge_index_;
            frontier[w_idx] = two_edge_index_;
        } else {
            frontier[old_w] = old_v;
            frontier[old_v] = old_w;
            frontier[w_idx] = two_edge_index_;
            frontier[v_idx] = two_edge_index_;
        }
    }

    advance(frontier, bag_idx);
}

void TreewidthSearch::skip(Frontier& frontier, size_t bag_idx) {
    advance(frontier, bag_idx);
}

bool TreewidthSearch::merge(Frontier& left, Frontier const& right, size_t bag_idx, std::vector<Edge_weight>& left_result, std::vector<Edge_weight> const& right_result) {
    if(left_result[0] + right_result[0] > max_length_) {
        return false;
    }    
    bool left_empty = left_result[0] == 0;
    bool right_empty = right_result[0] == 0;
    // merge the frontiers
    bool found_solution = false;
    for(size_t idx = 0; idx < right.size(); idx++) {
        // start goal case
        if(!is_all_pair_ && (bag_local_idx_map_[bag_idx][terminals_[0]] == idx || bag_local_idx_map_[bag_idx][terminals_[1]] == idx)) {
            assert(left[idx] == invalid_index_ || left[idx] == two_edge_index_);
            assert(right[idx] == invalid_index_ || right[idx] == two_edge_index_);
            if(left[idx] == invalid_index_) {
                left[idx] = right[idx];
            } else if(right[idx] == two_edge_index_) {
                return false;
            }
            continue;
        }
        // rest
        if(right[idx] == no_edge_index_) {
            continue;
        } else if(right[idx] == two_edge_index_) {
            if(left[idx] != no_edge_index_) {
                return false;
            }
            left[idx] = two_edge_index_;
        } else if(right[idx] == invalid_index_) {
            if(left[idx] == no_edge_index_) {
                left[idx] = invalid_index_;
            } else if(left[idx] == two_edge_index_) {
                return false;
            } else if(left[idx] == invalid_index_) {
                if(found_solution) {
                    return false;
                }
                left[idx] = two_edge_index_;
                found_solution = true;
            } else {
                assert(left[left[idx]] == idx);
                left[left[idx]] = invalid_index_;
                left[idx] = two_edge_index_;
            }
        } else if(right[idx] > idx) {
            assert(right[right[idx]] == idx);
            if(left[idx] == no_edge_index_) {
                if(left[right[idx]] == no_edge_index_) {
                    left[idx] = right[idx];
                    left[right[idx]] = idx;
                } else if(left[right[idx]] == two_edge_index_) {
                    return false;
                } else if(left[right[idx]] == invalid_index_) {
                    left[idx] = invalid_index_;
                    left[right[idx]] = two_edge_index_;
                } else {
                    left[idx] = left[right[idx]];
                    left[left[right[idx]]] = idx;
                    left[right[idx]] = two_edge_index_;
                }
            } else if(left[idx] == two_edge_index_) {
                return false;
            } else if(left[idx] == invalid_index_) {
                if(left[right[idx]] == no_edge_index_) {
                    left[idx] = two_edge_index_;
                    left[right[idx]] = invalid_index_;
                } else if(left[right[idx]] == two_edge_index_) {
                    return false;
                } else if(left[right[idx]] == invalid_index_) {
                    if(found_solution) {
                        return false;
                    }
                    left[idx] = two_edge_index_;
                    left[right[idx]] = two_edge_index_;
                    found_solution = true;
                } else {
                    left[idx] = two_edge_index_;
                    left[left[right[idx]]] = invalid_index_;
                    left[right[idx]] = two_edge_index_;
                }
            } else {
                if(left[right[idx]] == no_edge_index_) {
                    left[left[idx]] = right[idx];
                    left[right[idx]] = left[idx];
                    left[idx] = two_edge_index_;
                } else if(left[right[idx]] == two_edge_index_) {
                    return false;
                } else if(left[right[idx]] == invalid_index_) {
                    left[left[idx]] = invalid_index_;
                    left[idx] = two_edge_index_;
                    left[right[idx]] = two_edge_index_;
                } else {
                    if(left[right[idx]] == idx) {
                        assert(left[idx] == right[idx]);
                        assert(right[left[idx]] == idx);
                        // found a loop
                        return false;
                    }
                    left[left[idx]] = left[right[idx]];
                    left[left[right[idx]]] = left[idx];
                    left[idx] = two_edge_index_;
                    left[right[idx]] = two_edge_index_;
                }
            }
        }
    }
    assert(!left_empty || !right_empty || !found_solution);
    // mark out of scope edge ends and populate paths, cut_paths
    std::vector<frontier_index_t> paths;
    std::vector<frontier_index_t> cut_paths;
    for(frontier_index_t idx = 0; idx < right.size(); idx++) {
        // start goal case
        if(!is_all_pair_ && (bag_local_idx_map_[bag_idx][terminals_[0]] == idx || bag_local_idx_map_[bag_idx][terminals_[1]] == idx)) {
            if(left[idx] == invalid_index_) {
                if(!remaining_edges_after_this_[bag_idx][idx]) {
                    return false;
                }
                cut_paths.push_back(idx);
            }
            continue;
        }
        // rest
        if(left[idx] == invalid_index_) {
            if(!remaining_edges_after_this_[bag_idx][idx]) {
                if(found_solution) {
                    return false;
                }
                left[idx] = two_edge_index_;
                found_solution = true;
            } else {
                cut_paths.push_back(idx);
            }
        } else if(left[idx] != no_edge_index_ && left[idx] != two_edge_index_ && idx > left[idx]) {
            assert(left[left[idx]] == idx);
            if(!remaining_edges_after_this_[bag_idx][idx] && !remaining_edges_after_this_[bag_idx][left[idx]]) {
                if(found_solution) {
                    return false;
                }
                left[left[idx]] = two_edge_index_;
                left[idx] = two_edge_index_;
                found_solution = true;
            } else if(!remaining_edges_after_this_[bag_idx][idx] && remaining_edges_after_this_[bag_idx][left[idx]]) {
                left[left[idx]] = invalid_index_;
                cut_paths.push_back(left[idx]);
                left[idx] = two_edge_index_;
            } else if(remaining_edges_after_this_[bag_idx][idx] && !remaining_edges_after_this_[bag_idx][left[idx]]) {
                left[left[idx]] = two_edge_index_;
                left[idx] = invalid_index_;
                cut_paths.push_back(idx);
            } else {
                paths.push_back(idx);
                paths.push_back(left[idx]);
            }
        }
    }

    assert(paths.size() % 2 == 0);
    // adapt result
    // FIXME: optimize this
    std::vector<Edge_weight> new_result(left_result.size() + right_result.size() - 2, 0);
    new_result[0] = left_result[0] + right_result[0];
    for(size_t l_length = 1; l_length < left_result.size(); l_length++) {
        for(size_t r_length = 1; r_length < right_result.size(); r_length++) {
            new_result[l_length + r_length - 1] += left_result[l_length]*right_result[r_length];
        }
    }
    left_result = new_result;
    left_result.resize(std::min(left_result.size(), max_length_ - left_result[0] + 2));

    // possibly include solutions
    size_t thread_id = omp_get_thread_num();
    if(found_solution) {
        if(left_empty || right_empty) {
            return false;
        }
        // we cannot continue this either way but if there are not other partial paths we have to include the solution
        if(paths.size() + cut_paths.size() == 0) {
            for(Edge_length res_length = 1; res_length < left_result.size() && left_result[0] + res_length - 1 <= max_length_; res_length++) {
                thread_local_result_[thread_id][left_result[0] + res_length - 1] += left_result[res_length];
            }
        }
        return false;
    }
    if(!left_empty && !right_empty && paths.size()/2 + cut_paths.size() == 1) {
        // if either is empty then we already counted the solutions
        // this way there is currently exactly one path and we have not seen it yet
        for(Edge_length res_length = 1; res_length < left_result.size() && left_result[0] + res_length - 1 <= max_length_; res_length++) {
            thread_local_result_[thread_id][left_result[0] + res_length - 1] += left_result[res_length];
        }
    }

    if(cut_paths.size() > 2) {
        return false;
    }

    // - 1 since have to take one edge less
    if(paths.size()/2 + cut_paths.size() > 1) {
        if(paths.size()/2 + cut_paths.size() + left_result[0] > max_length_ + 1) {
            return false;
        }
    } else if(paths.size()/2 + cut_paths.size() + left_result[0] > max_length_) {
        return false;
    }

    if(!distancePrune(left, paths, cut_paths, bag_idx, left_result[0])) {
        return false;
    }
    advance(left, bag_idx);
    return true;
}

void TreewidthSearch::advance(Frontier& frontier, size_t bag_idx) {
    Frontier old = frontier;
    size_t next_idx = decomposition_[bag_idx].parent;
    auto &bag = decomposition_[next_idx].bag;
    frontier.resize(bag.size());
    std::fill(frontier.begin(), frontier.end(), no_edge_index_);
    auto &old_idx = bag_local_idx_map_[bag_idx];
    auto &new_idx = bag_local_idx_map_[next_idx];
    auto &old_vertex = bag_local_vertex_map_[bag_idx];
    int found_two = 0;
    int found_invalid = 0;
    int found_path = 0;
    for(auto v : bag) {
        if(old_idx[v] != invalid_index_) {
            if(old[old_idx[v]] == no_edge_index_) {
                frontier[new_idx[v]] = no_edge_index_;
            } else if(old[old_idx[v]] == two_edge_index_) {
                frontier[new_idx[v]] = two_edge_index_;
                found_two += 1;
            } else if(old[old_idx[v]] == invalid_index_) {
                frontier[new_idx[v]] = invalid_index_;
                found_invalid += 1;
            } else {
                frontier[new_idx[v]] = new_idx[old_vertex[old[old_idx[v]]]];
                if(frontier[new_idx[v]] != invalid_index_) {
                    found_path += 1;
                } else {
                    found_path += 2;
                }
            }
            assert(frontier[new_idx[v]] == no_edge_index_ || frontier[new_idx[v]] == two_edge_index_ || remaining_edges_after_this_[bag_idx][old_idx[v]] > 0);
        }
    }
    assert(found_path % 2 == 0);
    assert(!found_two || found_invalid || found_path);
    assert(found_invalid <= 2);
    if(found_invalid == 2) {
        for(frontier_index_t idx = 0; idx < old.size(); idx++) {
            if(old[idx] == no_edge_index_ && remaining_edges_after_this_[bag_idx][idx] == 1) {
                auto advanced_idx = new_idx[old_vertex[idx]];
                frontier[advanced_idx] = two_edge_index_;
            }
        }
    }
}


void TreewidthSearch::print_stats() {
    size_t pos_hits = 0, neg_hits = 0;
    for(size_t i = 0; i < nthreads_; i++) {
        pos_hits += pos_hits_[i];
        neg_hits += neg_hits_[i];
    }
    size_t edges = 0, propagations = 0;
    for(size_t i = 0; i < nthreads_; i++) {
        edges += edges_[i];
        propagations += propagations_[i];
    }
    size_t merges = 0, unsuccessful_merges = 0;
    for(size_t i = 0; i < nthreads_; i++) {
        merges += merges_[i];
        unsuccessful_merges += unsuccessful_merges_[i];
    }
    std::cerr << "Cache hit rate: " << 100*pos_hits/(double)(pos_hits + neg_hits) << "% (" << pos_hits << "/" << pos_hits + neg_hits << ")" << std::endl;
    std::cerr << "#Edges: " << edges << " #Propagations: " << propagations << std::endl; 
    std::cerr << "#Merges: " << merges << " #Unsuccessful merges: " << unsuccessful_merges << std::endl;
}


}
