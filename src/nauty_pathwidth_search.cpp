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

#include "nauty_pathwidth_search.h"
#include <algorithm>
namespace fpc {

void prettyPrint(Frontier frontier) {
    for(auto idx : frontier) {
        if(idx <= 252) {
            std::cerr << size_t(idx) << " ";
        } else if(idx == 253) {
            std::cerr << "-" << " ";
        } else if(idx == 254) {
            std::cerr << "#" << " ";
        } else  {
            std::cerr << "*" << " ";
        }
    }
    std::cerr << std::endl;
}

NautyPathwidthSearch::NautyPathwidthSearch(Graph& input, AnnotatedDecomposition decomposition, size_t nthreads) 
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
        sparsegraph_after_this_(decomposition_.size()),
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
	//std::cerr << bag_idx << " ";
	//node.stats();
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
        auto &sg = sparsegraph_after_this_[bag_idx];
        SG_INIT(sg);
        size_t nr_vertices = cur_graph.nr_vertices();
        for(auto v : node.bag) {
            if(cur_graph.neighbors(v).empty()) {
                nr_vertices++;
            }
        }
        nr_vertices += node.bag.size();
        size_t nr_edges = 2*cur_graph.nr_edges();
        nr_edges += 3*node.bag.size();
        std::vector<size_t> new_name(graph_.adjacency_.size(), size_t(-1));
        std::vector<size_t> reverse(nr_vertices, size_t(-1));
        for(auto v : node.bag) {
            new_name[v] = bag_local_idx_map_[bag_idx][v] + node.bag.size();
            reverse[new_name[v]] = v;
        }
        size_t cur_name = 2*node.bag.size();
        for(Vertex v = 0; v < graph_.adjacency_.size(); v++) {
            if(cur_graph.neighbors(v).empty()) {// && v != node.edge.first && v != node.edge.second) {
                continue;
            }         
            if(new_name[v] == size_t(-1)) {
                new_name[v] = cur_name++;
                reverse[new_name[v]] = v;
            }
        }
        sg.v = (edge_t *)calloc(sizeof(edge_t)*(nr_vertices + nr_edges) + sizeof(degree_t)*nr_vertices,1);
        sg.d = (degree_t *)(sg.v + nr_vertices);
        sg.e = (edge_t *)(sg.d + nr_vertices);
        sg.nv = nr_vertices;
        sg.nde = nr_edges - 3*node.bag.size();
        sg.vlen = nr_vertices;
        sg.dlen = nr_vertices;
        sg.elen = nr_edges;
        size_t cur_e_idx = 0;
        for(size_t i = 0; i < node.bag.size(); i++) {
            sg.v[i] = cur_e_idx;
            sg.d[i] = 0;
            cur_e_idx += 2;
        }
        for(size_t i = node.bag.size(); i < nr_vertices; i++) {
            size_t v = reverse[i];
            assert(v != size_t(-1));
            sg.v[i] = cur_e_idx;
            for(auto w : cur_graph.neighbors(v)) {
                sg.e[cur_e_idx++] = new_name[w];
            }
            sg.d[i] = cur_e_idx - sg.v[i];
            if(i < 2*node.bag.size()) {
                cur_e_idx++;
            }
        }
        assert(cur_e_idx == nr_edges);
    }
    for(vertex_t v = 0; v < graph_.adjacency_.size(); v++) {
        if(cur_graph.neighbors(v).size() != 0) {
            std::cerr << v << std::endl;
        }
        assert(cur_graph.neighbors(v).size() == 0);
    }
    for(size_t bag_idx = 0; bag_idx < decomposition_.size(); bag_idx++) {
        auto &node = decomposition_[bag_idx];
        // std::cerr << "<" << bag_idx << ">"; node.stats();
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
            auto sg = construct_sparsegraph(initial_frontier, size_t(-1));
            cache_[bag_idx].first[std::make_pair(sg, initial_frontier)] = initial_result;
        }
    }
}

std::vector<Edge_weight> NautyPathwidthSearch::search() {
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
            #pragma omp parallel for default(shared) //shared(edges_) shared(propagations_) shared(bag_idx) shared(max_length_) shared(cache_) shared(decomposition_) shared(thread_local_result_) shared(pos_hits_) shared(neg_hits_) shared(bag_local_idx_map_) shared(bag_local_vertex_map_)
            for(size_t bucket = 0; bucket < cache_[bag_idx].first.bucket_count(); bucket++) {
                size_t thread_id = omp_get_thread_num();
                for(auto task_it = cache_[bag_idx].first.begin(bucket); task_it != cache_[bag_idx].first.end(bucket); ++task_it) {
                    auto copy_frontier = task_it->first.second;
                    auto copy_result = task_it->second;
                    propagateLoop(copy_frontier, bag_idx, -1, copy_result, true, false, thread_id);
                    copy_frontier = task_it->first.second;
                    copy_result = task_it->second;
                    propagateLoop(copy_frontier, bag_idx, -1, copy_result, false, true, thread_id);
                    free(task_it->first.first.v);
                }
            }
        } else {
            if(cache_[bag_idx].first.size() > cache_[bag_idx].second.size()) {
                std::swap(cache_[bag_idx].first, cache_[bag_idx].second);
            }
            std::vector<std::pair<Frontier, std::vector<Edge_weight>>> right_vector;
            std::transform(
                cache_[bag_idx].second.begin(),
                cache_[bag_idx].second.end(), 
                std::back_inserter(right_vector),
                [](std::pair<TWCacheKey, std::vector<Edge_weight>> const& key) {
                    free(key.first.first.v);
                    return std::make_pair(key.first.second, key.second);
            });
            std::sort(right_vector.begin(), right_vector.end());
            // std::cerr << "Remaining: ";
            // for(auto i = 0; i < remaining_edges_after_this_[bag_idx].size(); i++) {
            //     std::cerr << size_t(remaining_edges_after_this_[bag_idx][i]) << " ";
            // }
            // std::cerr << std::endl;
            // std::cerr << std::endl;
            // for(size_t bucket = 0; bucket < cache_[bag_idx].first.bucket_count(); bucket++) {
            //     for(auto task_it = cache_[bag_idx].first.begin(bucket); task_it != cache_[bag_idx].first.end(bucket); ++task_it) {
            //         for(auto [frontier, result] : right_vector) {
            //             prettyPrint(task_it->first);
            //             prettyPrint(frontier);
            //             auto left_copy = task_it->first;
            //             auto right_copy = frontier;
            //             auto left_result_copy = task_it->second;
            //             auto right_result_copy = result;
            //             if(merge(left_copy, right_copy, bag_idx, left_result_copy, right_result_copy)) {
            //                 std::cerr << "Merged" << std::endl;
            //             } else {
            //                 std::cerr << "Not Merged" << std::endl;
            //             }
            //             std::vector<std::pair<Frontier, std::vector<Edge_weight>>> single = {std::make_pair(right_copy, right_result_copy)};
            //             auto begin = single.begin();
            //             auto end = single.end();
            //             std::vector<frontier_index_t> cut_paths;
            //             std::vector<frontier_index_t> paths;
            //             bool found_solution = false;
            //             auto left = task_it->first;
            //             auto left_result = task_it->second;
            //             mergeStep(left, bag_idx, 0, found_solution, cut_paths, paths, left_result, begin, end);
            //         }
            //     }
            // }
            // std::cerr << std::endl;
            // continue;
            // JOIN
            // FIXME: check if its better to take either bag as the outer one
            #pragma omp parallel for default(shared) shared(merges_) shared(unsuccessful_merges_) shared(edges_) shared(propagations_) shared(bag_idx) shared(max_length_) shared(cache_) shared(decomposition_) shared(thread_local_result_) shared(pos_hits_) shared(neg_hits_) shared(bag_local_idx_map_) shared(bag_local_vertex_map_)
            for(size_t bucket = 0; bucket < cache_[bag_idx].first.bucket_count(); bucket++) {
                size_t thread_id = omp_get_thread_num();
                for(auto task_it = cache_[bag_idx].first.begin(bucket); task_it != cache_[bag_idx].first.end(bucket); ++task_it) {
                    free(task_it->first.first.v);
                    if(right_vector.size() == 0) {
                        continue;
                    }
                    auto begin = right_vector.begin();
                    auto end = right_vector.end();
                    std::vector<frontier_index_t> cut_paths;
                    std::vector<frontier_index_t> paths;
                    auto right = right_vector[0].first;
                    bool found_solution = false;
                    auto left = task_it->first.second;
                    auto left_result = task_it->second;
                    if(right[0] <= 252) {
                        // open interval constructor
                        Frontier next(right.begin(), right.begin());
                        std::vector<Edge_weight> empty;
                        next.push_back(two_edge_index_);
                        auto middle = std::lower_bound(begin, end, std::make_pair(next, empty));
                        auto value_it = begin;
                        while(value_it != middle) {
                            begin = value_it;
                            frontier_index_t cur = begin->first[0];
                            assert(cur < 252);
                            next.back() = cur + 1;
                            value_it = std::lower_bound(begin, middle, std::make_pair(next, empty));
                            mergeStep(left, bag_idx, 0, found_solution, cut_paths, paths, left_result, begin, value_it);
                        }
                        next.back() = no_edge_index_;
                        begin = middle;
                        middle = std::lower_bound(middle, end, std::make_pair(next, empty));
                        mergeStep(left, bag_idx, 0, found_solution, cut_paths, paths, left_result, begin, middle);
                        next.back() = invalid_index_;
                        begin = middle;
                        middle = std::lower_bound(middle, end, std::make_pair(next, empty));
                        mergeStep(left, bag_idx, 0, found_solution, cut_paths, paths, left_result, begin, middle);
                        mergeStep(left, bag_idx, 0, found_solution, cut_paths, paths, left_result, middle, end);
                    } else if(right[0] == two_edge_index_) {
                        // open interval constructor
                        Frontier next(right.begin(), right.begin());
                        std::vector<Edge_weight> empty;
                        next.push_back(no_edge_index_);
                        auto middle = std::lower_bound(begin, end, std::make_pair(next, empty));
                        mergeStep(left, bag_idx, 0, found_solution, cut_paths, paths, left_result, begin, middle);
                        next.back() = invalid_index_;
                        begin = middle;
                        middle = std::lower_bound(middle, end, std::make_pair(next, empty));
                        mergeStep(left, bag_idx, 0, found_solution, cut_paths, paths, left_result, begin, middle);
                        mergeStep(left, bag_idx, 0, found_solution, cut_paths, paths, left_result, middle, end);
                    } else if(right[0] == no_edge_index_) {
                        // open interval constructor
                        Frontier next(right.begin(), right.begin() + 0);
                        std::vector<Edge_weight> empty;
                        next.push_back(invalid_index_);
                        auto middle = std::lower_bound(begin, end, std::make_pair(next, empty));
                        mergeStep(left, bag_idx, 0, found_solution, cut_paths, paths, left_result, begin, middle);
                        mergeStep(left, bag_idx, 0, found_solution, cut_paths, paths, left_result, middle, end);
                    } else if(right[0] == invalid_index_) {
                        mergeStep(left, bag_idx, 0, found_solution, cut_paths, paths, left_result, begin, end);
                    }
                    continue;
                    // std::vector<std::pair<size_t, frontier_index_t>> stack;
                    // stack.push_back(std::make_pair(0, 0));
                    // while(!stack.empty()) {
                    //     auto [cur_node, cur_idx] = stack.back();
                    //     stack.pop_back();
                    //     if(cur_node == size_t(-1)) {
                    //         continue;
                    //     }
                    //     if(cur_idx < task_it->first.size()) {
                    //         int mergeable_idx = std::max(task_it->first[cur_idx], frontier_index_t(252)) - 252;
                    //         if(!is_all_pair_ && (bag_local_idx_map_[bag_idx][terminals_[0]] == cur_idx || bag_local_idx_map_[bag_idx][terminals_[1]] == cur_idx)) {
                    //             if(mergeable_idx == 1) {
                    //                 stack.push_back(std::make_pair(nodes[cur_node].children[3], cur_idx + 1));
                    //             } else {
                    //                 assert(mergeable_idx == 3);
                    //                 stack.push_back(std::make_pair(nodes[cur_node].children[1], cur_idx + 1));
                    //                 stack.push_back(std::make_pair(nodes[cur_node].children[3], cur_idx + 1));
                    //             }
                    //         } else {
                    //             if(mergeable_idx == 0) { // path
                    //                 stack.push_back(std::make_pair(nodes[cur_node].children[0], cur_idx + 1));
                    //                 stack.push_back(std::make_pair(nodes[cur_node].children[2], cur_idx + 1));
                    //                 stack.push_back(std::make_pair(nodes[cur_node].children[3], cur_idx + 1));
                    //             } else if(mergeable_idx == 1) { // two edges
                    //                 stack.push_back(std::make_pair(nodes[cur_node].children[2], cur_idx + 1));
                    //             } else if(mergeable_idx == 2) { // no edge
                    //                 stack.push_back(std::make_pair(nodes[cur_node].children[0], cur_idx + 1));
                    //                 stack.push_back(std::make_pair(nodes[cur_node].children[1], cur_idx + 1));
                    //                 stack.push_back(std::make_pair(nodes[cur_node].children[2], cur_idx + 1));
                    //                 stack.push_back(std::make_pair(nodes[cur_node].children[3], cur_idx + 1));
                    //             } else if(mergeable_idx == 3) { // cut path
                    //                 stack.push_back(std::make_pair(nodes[cur_node].children[0], cur_idx + 1));
                    //                 stack.push_back(std::make_pair(nodes[cur_node].children[2], cur_idx + 1));
                    //                 stack.push_back(std::make_pair(nodes[cur_node].children[3], cur_idx + 1));
                    //             }
                    //         }
                    //         continue;
                    //     }
                    //     for(auto cache_idx : nodes[cur_node].cache_idx) {
                    //         auto const&[right_frontier, right_result] = right_vector[cache_idx];
                    //         if(right_result[0] + task_it->second[0] > max_length_) {
                    //             break;
                    //         }
                    //         auto left_frontier = task_it->first;
                    //         auto left_result = task_it->second;
                    //         // for(auto idx : left_frontier) {
                    //         //     std::cerr << size_t(idx) << " ";
                    //         // }
                    //         // std::cerr << std::endl;
                    //         // for(auto idx : right_frontier) {
                    //         //     std::cerr << size_t(idx) << " ";
                    //         // }
                    //         // std::cerr << std::endl;
                    //         // std::cerr << "Remaining: ";
                    //         // for(auto i = 0; i < left_frontier.size(); i++) {
                    //         //     std::cerr << size_t(remaining_edges_after_this_[bag_idx][i]) << " ";
                    //         // }
                    //         // std::cerr << std::endl;
                    //         merges_[thread_id]++;
                    //         if(!merge(left_frontier, right_frontier, bag_idx, left_result, right_result)) {
                    //             // std::cerr << "Not merged" << std::endl;
                    //             unsuccessful_merges_[thread_id]++;
                    //             continue;
                    //         }
                    //         // std::cerr << "Merged" << std::endl;
                    //     }
                    // }
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
    for(size_t bag_idx = 0; bag_idx < decomposition_.size(); bag_idx++) {
        free(sparsegraph_after_this_[bag_idx].v);
    }
    std::cerr << std::endl;
    return result_;
}

void NautyPathwidthSearch::propagateLoop(Frontier &frontier, size_t bag_idx, size_t last_idx, std::vector<Edge_weight>& partial_results, bool takeable, bool skippable, size_t thread_id) {
    size_t new_idx = bag_idx;
    if(decomposition_[new_idx].type != JOIN && (takeable ^ skippable)) {
        edges_[thread_id]++;
        propagations_[thread_id]--;
    }
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
            take(frontier, new_idx);
            partial_results[0]++;
            size_t paths = 0;
            for(frontier_index_t idx : frontier) {
                if(idx <= 252) {
                    paths++;
                } else if(idx == invalid_index_) {
                    paths += 2;
                }
            }
            if(paths/2 == 1) {
                for(Edge_length res_length = 1; res_length < partial_results.size() && partial_results[0] + res_length - 1 <= max_length_; res_length++) {
                    thread_local_result_[thread_id][partial_results[0] + res_length - 1] += partial_results[res_length];
                }
            } 
        } else {
            skip(frontier, new_idx);
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
            takeable = canTake(frontier, new_idx, partial_results);
            skippable = canSkip(frontier, new_idx, partial_results);
            if(!takeable) {
                includeSolutions(frontier, new_idx, partial_results);
            }
            // std::cerr << takeable << " " << skippable << std::endl;
        }
    }
    // both are possible, so we have a new decision edge
    // put it into the cache
    if(new_idx < decomposition_.size() 
        && (
                (takeable && skippable)
            ||  (decomposition_[new_idx].type == JOIN && last_idx == decomposition_[new_idx].children.first))) {
        auto sg = construct_sparsegraph(frontier, last_idx);
        #pragma omp critical 
        {
            auto ins = cache_[new_idx].first.insert(
                std::make_pair(std::make_pair(sg,frontier), partial_results)
            );
            if(!ins.second) {
                free(sg.v);
                pos_hits_[thread_id]++;
                // there is already an element with that key
                // instead increase the partial result for that key
                auto old_offset = ins.first->second[0];
                auto new_offset = partial_results[0];
                if(old_offset > new_offset) {
                    // we reached this state with a smaller offset
                    // add the paths we currently cached to the new result
                    if(partial_results.size() + new_offset < ins.first->second.size() + old_offset) {
                        partial_results.resize(ins.first->second.size() + old_offset - new_offset);
                    }
                    for(Edge_length res_length = 1; res_length < ins.first->second.size(); res_length++) {
                        partial_results[old_offset - new_offset + res_length] += ins.first->second[res_length];
                    }
                    ins.first->second = partial_results;
                } else {
                    if(ins.first->second.size() + old_offset < partial_results.size() + new_offset) {
                        ins.first->second.resize(partial_results.size() + new_offset - old_offset);
                    }
                    for(Edge_length res_length = 1; res_length < partial_results.size(); res_length++) {
                        ins.first->second[new_offset - old_offset + res_length] += partial_results[res_length];
                    }
                }
            } else {
                neg_hits_[thread_id]++;
            }
        }
    } else if(new_idx < decomposition_.size() 
        && decomposition_[new_idx].type == JOIN) {
        assert(last_idx == decomposition_[new_idx].children.second);
        auto sg = construct_sparsegraph(frontier, last_idx);
        #pragma omp critical 
        {
            auto ins = cache_[new_idx].second.insert(
                std::make_pair(std::make_pair(sg,frontier), partial_results)
            );
            if(!ins.second) {
                free(sg.v);
                pos_hits_[thread_id]++;
                // there is already an element with that key
                // instead increase the partial result for that key
                auto old_offset = ins.first->second[0];
                auto new_offset = partial_results[0];
                if(old_offset > new_offset) {
                    // we reached this state with a smaller offset
                    // add the paths we currently cached to the new result
                    if(partial_results.size() + new_offset < ins.first->second.size() + old_offset) {
                        partial_results.resize(ins.first->second.size() + old_offset - new_offset);
                    }
                    for(Edge_length res_length = 1; res_length < ins.first->second.size(); res_length++) {
                        partial_results[old_offset - new_offset + res_length] += ins.first->second[res_length];
                    }
                    ins.first->second = partial_results;
                } else {
                    if(ins.first->second.size() + old_offset < partial_results.size() + new_offset) {
                        ins.first->second.resize(partial_results.size() + new_offset - old_offset);
                    }
                    for(Edge_length res_length = 1; res_length < partial_results.size(); res_length++) {
                        ins.first->second[new_offset - old_offset + res_length] += partial_results[res_length];
                    }
                }
            } else {
                neg_hits_[thread_id]++;
            }
        }
    }
}

void NautyPathwidthSearch::includeSolutions(Frontier const& frontier, size_t bag_idx, std::vector<Edge_weight> const& partial_result) {
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

bool NautyPathwidthSearch::canTake(Frontier& frontier, size_t bag_idx, std::vector<Edge_weight> const& partial_result) {
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

bool NautyPathwidthSearch::canSkip(Frontier& frontier, size_t bag_idx, std::vector<Edge_weight> const& partial_result) {
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

bool NautyPathwidthSearch::distancePrune(
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

void NautyPathwidthSearch::take(Frontier& frontier, size_t bag_idx) {
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

void NautyPathwidthSearch::skip(Frontier& frontier, size_t bag_idx) {
    advance(frontier, bag_idx);
}

void NautyPathwidthSearch::restoreStep(
        Frontier &left, 
        std::vector<frontier_index_t>& cut_paths,
        std::vector<frontier_index_t>& paths,
        std::vector<std::pair<frontier_index_t, frontier_index_t>> restore,
        size_t cut_paths_size,
        size_t paths_size) {
    for(size_t i = restore.size(); i-- > 0;) {
        auto [idx, value] = restore[i];
        left[idx] = value;
    }
    cut_paths.resize(cut_paths_size);
    paths.resize(paths_size);
}

void NautyPathwidthSearch::mergeStep(
        Frontier &left,
        size_t bag_idx,
        frontier_index_t idx,
        bool found_solution,
        std::vector<frontier_index_t>& cut_paths,
        std::vector<frontier_index_t>& paths,
        std::vector<Edge_weight>& left_result,
        block_iter begin,
        block_iter end) {
    if(begin == end) {
        return;
    }
    // FIXME: dont restore if unnecessary
    // remember the previous values before change and restore before returning
    std::vector<std::pair<frontier_index_t, frontier_index_t>> restore;
    size_t cut_paths_size = cut_paths.size();
    size_t paths_size = paths.size();
    auto const& right = begin->first;
    // std::cerr << size_t(idx) << std::endl;
    // prettyPrint(left);
    // prettyPrint(right);
    // start goal case
    if(!is_all_pair_ && (bag_local_idx_map_[bag_idx][terminals_[0]] == idx || bag_local_idx_map_[bag_idx][terminals_[1]] == idx)) {
        assert(left[idx] == invalid_index_ || left[idx] == two_edge_index_);
        assert(right[idx] == invalid_index_ || right[idx] == two_edge_index_);
        // merge
        if(left[idx] == invalid_index_) {
            restore.push_back({idx, left[idx]});
            left[idx] = right[idx];
        } else if(right[idx] == two_edge_index_) {
            restoreStep(left, cut_paths, paths, restore, cut_paths_size, paths_size);
            return;
        }
        // paths
        if(left[idx] == invalid_index_) {
            if(!remaining_edges_after_this_[bag_idx][idx]) {
                restoreStep(left, cut_paths, paths, restore, cut_paths_size, paths_size);
                return;
            }
            cut_paths.push_back(idx);
        }
    } else {
        // rest
        if(right[idx] == no_edge_index_) {
            // standard path tracking
        } else if(right[idx] == two_edge_index_) {
            if(left[idx] != no_edge_index_) {
                restoreStep(left, cut_paths, paths, restore, cut_paths_size, paths_size);
                return;
            }
            restore.push_back({idx, left[idx]});
            left[idx] = two_edge_index_;
            // no path tracking
        } else if(right[idx] == invalid_index_) {
            if(left[idx] == no_edge_index_) {
                restore.push_back({idx, left[idx]});
                left[idx] = invalid_index_;
                // standard path tracking
            } else if(left[idx] == two_edge_index_) {
                restoreStep(left, cut_paths, paths, restore, cut_paths_size, paths_size);
                return;
            } else if(left[idx] == invalid_index_) {
                if(found_solution) {
                    restoreStep(left, cut_paths, paths, restore, cut_paths_size, paths_size);
                    return;
                }
                restore.push_back({idx, left[idx]});
                left[idx] = two_edge_index_;
                found_solution = true;
                // no path tracking
            } else {
                assert(left[left[idx]] == idx);
                if(left[idx] < idx) {
                    // time travel path tracking!
                    frontier_index_t aux_idx = left[idx]; 
                    restore.push_back({left[idx], left[left[idx]]});
                    left[left[idx]] = invalid_index_;
                    restore.push_back({idx, left[idx]});
                    left[idx] = two_edge_index_;
                    if(!remaining_edges_after_this_[bag_idx][aux_idx]) {
                        if(found_solution) {
                            restoreStep(left, cut_paths, paths, restore, cut_paths_size, paths_size);
                            return;
                        }
                        restore.push_back({aux_idx, left[aux_idx]});
                        left[aux_idx] = two_edge_index_;
                        found_solution = true;
                    } else {
                        cut_paths.push_back(aux_idx);
                    }
                } else {
                    // otherwise standard path tracking
                    restore.push_back({left[idx], left[left[idx]]});
                    left[left[idx]] = invalid_index_;
                    restore.push_back({idx, left[idx]});
                    left[idx] = two_edge_index_;
                }
            }
        } else if(right[idx] > idx) {
            assert(right[right[idx]] == idx);
            if(left[idx] == no_edge_index_) {
                if(left[right[idx]] == no_edge_index_) {
                    restore.push_back({idx, left[idx]});
                    left[idx] = right[idx];
                    restore.push_back({right[idx], left[right[idx]]});
                    left[right[idx]] = idx;
                    // standard path tracking
                } else if(left[right[idx]] == two_edge_index_) {
                    restoreStep(left, cut_paths, paths, restore, cut_paths_size, paths_size);
                    return;
                } else if(left[right[idx]] == invalid_index_) {
                    restore.push_back({idx, left[idx]});
                    left[idx] = invalid_index_;
                    restore.push_back({right[idx], left[right[idx]]});
                    left[right[idx]] = two_edge_index_;
                    // standard path tracking
                } else {
                    restore.push_back({idx, left[idx]});
                    left[idx] = left[right[idx]];
                    restore.push_back({left[right[idx]], left[left[right[idx]]]});
                    left[left[right[idx]]] = idx;
                    restore.push_back({right[idx], left[right[idx]]});
                    left[right[idx]] = two_edge_index_;
                    // standard path tracking
                }
            } else if(left[idx] == two_edge_index_) {
                restoreStep(left, cut_paths, paths, restore, cut_paths_size, paths_size);
                return;
            } else if(left[idx] == invalid_index_) {
                if(left[right[idx]] == no_edge_index_) {
                    restore.push_back({idx, left[idx]});
                    left[idx] = two_edge_index_;
                    restore.push_back({right[idx], left[right[idx]]});
                    left[right[idx]] = invalid_index_;
                    // standard path tracking
                } else if(left[right[idx]] == two_edge_index_) {
                    restoreStep(left, cut_paths, paths, restore, cut_paths_size, paths_size);
                    return;
                } else if(left[right[idx]] == invalid_index_) {
                    if(found_solution) {
                        restoreStep(left, cut_paths, paths, restore, cut_paths_size, paths_size);
                        return;
                    }
                    restore.push_back({idx, left[idx]});
                    left[idx] = two_edge_index_;
                    restore.push_back({right[idx], left[right[idx]]});
                    left[right[idx]] = two_edge_index_;
                    found_solution = true;
                    // no path tracking
                } else {
                    if(left[right[idx]] < idx) {
                        // time travel path tracking!
                        frontier_index_t aux_idx = left[right[idx]]; 
                        restore.push_back({idx, left[idx]});
                        left[idx] = two_edge_index_;
                        restore.push_back({left[right[idx]], left[left[right[idx]]]});
                        left[left[right[idx]]] = invalid_index_;
                        restore.push_back({right[idx], left[right[idx]]});
                        left[right[idx]] = two_edge_index_;
                        if(!remaining_edges_after_this_[bag_idx][aux_idx]) {
                            if(found_solution) {
                                restoreStep(left, cut_paths, paths, restore, cut_paths_size, paths_size);
                                return;
                            }
                            restore.push_back({aux_idx, left[aux_idx]});
                            left[aux_idx] = two_edge_index_;
                            found_solution = true;
                        } else {
                            cut_paths.push_back(aux_idx);
                        }
                    } else {
                        // otherwise standard path tracking
                        restore.push_back({idx, left[idx]});
                        left[idx] = two_edge_index_;
                        restore.push_back({left[right[idx]], left[left[right[idx]]]});
                        left[left[right[idx]]] = invalid_index_;
                        restore.push_back({right[idx], left[right[idx]]});
                        left[right[idx]] = two_edge_index_;
                    }
                }
            } else {
                if(left[right[idx]] == no_edge_index_) {
                    restore.push_back({left[idx], left[left[idx]]});
                    left[left[idx]] = right[idx];
                    restore.push_back({right[idx], left[right[idx]]});
                    left[right[idx]] = left[idx];
                    restore.push_back({idx, left[idx]});
                    left[idx] = two_edge_index_;
                    // standard path tracking
                } else if(left[right[idx]] == two_edge_index_) {
                    restoreStep(left, cut_paths, paths, restore, cut_paths_size, paths_size);
                    return;
                } else if(left[right[idx]] == invalid_index_) {
                    if(left[idx] < idx) {
                        // time travel path tracking!
                        frontier_index_t aux_idx = left[idx]; 
                        restore.push_back({left[idx], left[left[idx]]});
                        left[left[idx]] = invalid_index_;
                        restore.push_back({idx, left[idx]});
                        left[idx] = two_edge_index_;
                        restore.push_back({right[idx], left[right[idx]]});
                        left[right[idx]] = two_edge_index_;
                        if(!remaining_edges_after_this_[bag_idx][aux_idx]) {
                            if(found_solution) {
                                restoreStep(left, cut_paths, paths, restore, cut_paths_size, paths_size);
                                return;
                            }
                            restore.push_back({aux_idx, left[aux_idx]});
                            left[aux_idx] = two_edge_index_;
                            found_solution = true;
                        } else {
                            cut_paths.push_back(aux_idx);
                        }
                    } else {
                        // otherwise standard path tracking
                        restore.push_back({left[idx], left[left[idx]]});
                        left[left[idx]] = invalid_index_;
                        restore.push_back({idx, left[idx]});
                        left[idx] = two_edge_index_;
                        restore.push_back({right[idx], left[right[idx]]});
                        left[right[idx]] = two_edge_index_;
                    }
                } else {
                    if(left[right[idx]] == idx) {
                        assert(left[idx] == right[idx]);
                        assert(right[left[idx]] == idx);
                        // found a loop
                        restoreStep(left, cut_paths, paths, restore, cut_paths_size, paths_size);
                        return;
                    }
                    if(left[right[idx]] < idx && left[idx] < idx) {
                        frontier_index_t aux_idx = left[idx];
                        restore.push_back({left[idx], left[left[idx]]});
                        left[left[idx]] = left[right[idx]];
                        restore.push_back({left[right[idx]], left[left[right[idx]]]});
                        left[left[right[idx]]] = left[idx];
                        restore.push_back({idx, left[idx]});
                        left[idx] = two_edge_index_;
                        restore.push_back({right[idx], left[right[idx]]});
                        left[right[idx]] = two_edge_index_;
                        // time travel path tracking!
                        if(!remaining_edges_after_this_[bag_idx][aux_idx] && !remaining_edges_after_this_[bag_idx][left[aux_idx]]) {
                            if(found_solution) {
                                restoreStep(left, cut_paths, paths, restore, cut_paths_size, paths_size);
                                return;
                            }
                            restore.push_back({left[aux_idx], left[left[aux_idx]]});
                            left[left[aux_idx]] = two_edge_index_;
                            restore.push_back({aux_idx, left[aux_idx]});
                            left[aux_idx] = two_edge_index_;
                            found_solution = true;
                        } else if(!remaining_edges_after_this_[bag_idx][aux_idx] && remaining_edges_after_this_[bag_idx][left[aux_idx]]) {
                            restore.push_back({left[aux_idx], left[left[aux_idx]]});
                            left[left[aux_idx]] = invalid_index_;
                            cut_paths.push_back(left[aux_idx]);
                            restore.push_back({aux_idx, left[aux_idx]});
                            left[aux_idx] = two_edge_index_;
                        } else if(remaining_edges_after_this_[bag_idx][aux_idx] && !remaining_edges_after_this_[bag_idx][left[aux_idx]]) {
                            restore.push_back({left[aux_idx], left[left[aux_idx]]});
                            left[left[aux_idx]] = two_edge_index_;
                            restore.push_back({aux_idx, left[aux_idx]});
                            left[aux_idx] = invalid_index_;
                            cut_paths.push_back(aux_idx);
                        } else {
                            paths.push_back(aux_idx);
                            paths.push_back(left[aux_idx]);
                        }
                    } else {
                        restore.push_back({left[idx], left[left[idx]]});
                        left[left[idx]] = left[right[idx]];
                        restore.push_back({left[right[idx]], left[left[right[idx]]]});
                        left[left[right[idx]]] = left[idx];
                        restore.push_back({idx, left[idx]});
                        left[idx] = two_edge_index_;
                        restore.push_back({right[idx], left[right[idx]]});
                        left[right[idx]] = two_edge_index_;
                        // standard path tracking
                    }
                }
            }
        }
        // standard path tracking
        if(left[idx] == invalid_index_) {
            if(!remaining_edges_after_this_[bag_idx][idx]) {
                if(found_solution) {
                    restoreStep(left, cut_paths, paths, restore, cut_paths_size, paths_size);
                    return;
                }
                restore.push_back({idx, left[idx]});
                left[idx] = two_edge_index_;
                found_solution = true;
            } else {
                cut_paths.push_back(idx);
            }
        } else if(left[idx] != no_edge_index_ && left[idx] != two_edge_index_ && idx > left[idx]) {
            assert(left[left[idx]] == idx);
            if(!remaining_edges_after_this_[bag_idx][idx] && !remaining_edges_after_this_[bag_idx][left[idx]]) {
                if(found_solution) {
                    restoreStep(left, cut_paths, paths, restore, cut_paths_size, paths_size);
                    return;
                }
                restore.push_back({left[idx], left[left[idx]]});
                left[left[idx]] = two_edge_index_;
                restore.push_back({idx, left[idx]});
                left[idx] = two_edge_index_;
                found_solution = true;
            } else if(!remaining_edges_after_this_[bag_idx][idx] && remaining_edges_after_this_[bag_idx][left[idx]]) {
                restore.push_back({left[idx], left[left[idx]]});
                left[left[idx]] = invalid_index_;
                cut_paths.push_back(left[idx]);
                restore.push_back({idx, left[idx]});
                left[idx] = two_edge_index_;
            } else if(remaining_edges_after_this_[bag_idx][idx] && !remaining_edges_after_this_[bag_idx][left[idx]]) {
                restore.push_back({left[idx], left[left[idx]]});
                left[left[idx]] = two_edge_index_;
                restore.push_back({idx, left[idx]});
                left[idx] = invalid_index_;
                cut_paths.push_back(idx);
            } else {
                paths.push_back(idx);
                paths.push_back(left[idx]);
            }
        }
    }
    if(found_solution && (left_result[0] == 0 || paths.size() + cut_paths.size() != 0)) {
        restoreStep(left, cut_paths, paths, restore, cut_paths_size, paths_size);
        return;
    }

    if(cut_paths.size() > 2) {
        restoreStep(left, cut_paths, paths, restore, cut_paths_size, paths_size);
        return;
    }

    // - 1 since have to take one edge less
    if(paths.size()/2 + cut_paths.size() > 1) {
        if(paths.size()/2 + cut_paths.size() + left_result[0] > max_length_ + 1) {
            restoreStep(left, cut_paths, paths, restore, cut_paths_size, paths_size);
            return;
        }
    } else if(paths.size()/2 + cut_paths.size() + left_result[0] > max_length_) {
        restoreStep(left, cut_paths, paths, restore, cut_paths_size, paths_size);
        return;
    }

    // prettyPrint(left);
    if(idx + 1 < right.size()) {
        // proceed recursively
        if(right[idx + 1] <= 252) {
            // open interval constructor
            Frontier next(right.begin(), right.begin() + idx + 1);
            std::vector<Edge_weight> empty;
            next.push_back(two_edge_index_);
            auto middle = std::lower_bound(begin, end, std::make_pair(next, empty));
            auto value_it = begin;
            while(value_it != middle) {
                begin = value_it;
                frontier_index_t cur = begin->first[idx + 1];
                assert(cur < 252);
                next.back() = cur + 1;
                value_it = std::lower_bound(begin, middle, std::make_pair(next, empty));
                mergeStep(left, bag_idx, idx + 1, found_solution, cut_paths, paths, left_result, begin, value_it);
            }
            next.back() = no_edge_index_;
            begin = middle;
            middle = std::lower_bound(middle, end, std::make_pair(next, empty));
            mergeStep(left, bag_idx, idx + 1, found_solution, cut_paths, paths, left_result, begin, middle);
            next.back() = invalid_index_;
            begin = middle;
            middle = std::lower_bound(middle, end, std::make_pair(next, empty));
            mergeStep(left, bag_idx, idx + 1, found_solution, cut_paths, paths, left_result, begin, middle);
            mergeStep(left, bag_idx, idx + 1, found_solution, cut_paths, paths, left_result, middle, end);
        } else if(right[idx + 1] == two_edge_index_) {
            // open interval constructor
            Frontier next(right.begin(), right.begin() + idx + 1);
            std::vector<Edge_weight> empty;
            next.push_back(no_edge_index_);
            auto middle = std::lower_bound(begin, end, std::make_pair(next, empty));
            mergeStep(left, bag_idx, idx + 1, found_solution, cut_paths, paths, left_result, begin, middle);
            next.back() = invalid_index_;
            begin = middle;
            middle = std::lower_bound(middle, end, std::make_pair(next, empty));
            mergeStep(left, bag_idx, idx + 1, found_solution, cut_paths, paths, left_result, begin, middle);
            mergeStep(left, bag_idx, idx + 1, found_solution, cut_paths, paths, left_result, middle, end);
        } else if(right[idx + 1] == no_edge_index_) {
            // open interval constructor
            Frontier next(right.begin(), right.begin() + idx + 1);
            std::vector<Edge_weight> empty;
            next.push_back(invalid_index_);
            auto middle = std::lower_bound(begin, end, std::make_pair(next, empty));
            mergeStep(left, bag_idx, idx + 1, found_solution, cut_paths, paths, left_result, begin, middle);
            mergeStep(left, bag_idx, idx + 1, found_solution, cut_paths, paths, left_result, middle, end);
        } else if(right[idx + 1] == invalid_index_) {
            mergeStep(left, bag_idx, idx + 1, found_solution, cut_paths, paths, left_result, begin, end);
        }
    } else {
        // check if we want to keep the completely merged frontier
        // and if so merge
        assert(begin + 1 == end);
        size_t thread_id = omp_get_thread_num();
        merges_[thread_id]++;
        if(!finalizeMerge(left, bag_idx, found_solution, cut_paths, paths, left_result, begin->second)) {
            unsuccessful_merges_[thread_id]++;
        }
    }
    restoreStep(left, cut_paths, paths, restore, cut_paths_size, paths_size);
}

bool NautyPathwidthSearch::finalizeMerge(
      Frontier left,
      size_t bag_idx,
      bool found_solution,
      std::vector<frontier_index_t> const& cut_paths,
      std::vector<frontier_index_t> const& paths,
      std::vector<Edge_weight> const& left_result,
      std::vector<Edge_weight> const& right_result
    ) {
    bool left_empty = left_result[0] == 0;
    bool right_empty = right_result[0] == 0;

    assert(!left_empty || !right_empty || !found_solution);

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
    new_result.resize(std::min(new_result.size(), max_length_ - new_result[0] + 2));

    // possibly include solutions
    size_t thread_id = omp_get_thread_num();
    if(found_solution) {
        if(left_empty || right_empty) {
            return false;
        }
        // we cannot continue this either way but if there are not other partial paths we have to include the solution
        if(paths.size() + cut_paths.size() == 0) {
            for(Edge_length res_length = 1; res_length < new_result.size() && new_result[0] + res_length - 1 <= max_length_; res_length++) {
                thread_local_result_[thread_id][new_result[0] + res_length - 1] += new_result[res_length];
            }
        }
        return false;
    }
    if(!left_empty && !right_empty && paths.size()/2 + cut_paths.size() == 1) {
        // if either is empty then we already counted the solutions
        // this way there is currently exactly one path and we have not seen it yet
        for(Edge_length res_length = 1; res_length < new_result.size() && new_result[0] + res_length - 1 <= max_length_; res_length++) {
            thread_local_result_[thread_id][new_result[0] + res_length - 1] += new_result[res_length];
        }
    }

    if(cut_paths.size() > 2) {
        return false;
    }

    // - 1 since have to take one edge less
    if(paths.size()/2 + cut_paths.size() > 1) {
        if(paths.size()/2 + cut_paths.size() + new_result[0] > max_length_ + 1) {
            return false;
        }
    } else if(paths.size()/2 + cut_paths.size() + new_result[0] > max_length_) {
        return false;
    }

    if(!distancePrune(left, paths, cut_paths, bag_idx, new_result[0])) {
        return false;
    }
    // incorporate into cache
    advance(left, bag_idx);
    size_t last_idx = bag_idx;
    size_t new_idx = decomposition_[bag_idx].parent;
    assert(new_idx != size_t(-1));
    bool takeable = false;
    bool skippable = false;
    if(decomposition_[new_idx].type != JOIN) {
        takeable = canTake(left, new_idx, new_result);
        skippable = canSkip(left, new_idx, new_result);
        if(!takeable) {
            includeSolutions(left, new_idx, new_result);
        }
    }
    propagateLoop(left, new_idx, last_idx, new_result, takeable, skippable, thread_id);
    return true;
}

bool NautyPathwidthSearch::merge(Frontier& left, Frontier const& right, size_t bag_idx, std::vector<Edge_weight>& left_result, std::vector<Edge_weight> const& right_result) {
    // merge the frontiers and mark out of scope edge ends and populate paths, cut_paths
    bool found_solution = false;
    std::vector<frontier_index_t> paths;
    std::vector<frontier_index_t> cut_paths;
    for(size_t idx = 0; idx < right.size(); idx++) {
        // start goal case
        if(!is_all_pair_ && (bag_local_idx_map_[bag_idx][terminals_[0]] == idx || bag_local_idx_map_[bag_idx][terminals_[1]] == idx)) {
            assert(left[idx] == invalid_index_ || left[idx] == two_edge_index_);
            assert(right[idx] == invalid_index_ || right[idx] == two_edge_index_);
            // merge
            if(left[idx] == invalid_index_) {
                left[idx] = right[idx];
            } else if(right[idx] == two_edge_index_) {
                return false;
            }
            // paths
            if(left[idx] == invalid_index_) {
                if(!remaining_edges_after_this_[bag_idx][idx]) {
                    return false;
                }
                cut_paths.push_back(idx);
            }
            continue;
        }
        // rest
        if(right[idx] == no_edge_index_) {
            // standard path tracking
        } else if(right[idx] == two_edge_index_) {
            if(left[idx] != no_edge_index_) {
                return false;
            }
            left[idx] = two_edge_index_;
            // no path tracking
        } else if(right[idx] == invalid_index_) {
            if(left[idx] == no_edge_index_) {
                left[idx] = invalid_index_;
                // standard path tracking
            } else if(left[idx] == two_edge_index_) {
                return false;
            } else if(left[idx] == invalid_index_) {
                if(found_solution) {
                    return false;
                }
                left[idx] = two_edge_index_;
                found_solution = true;
                // no path tracking
            } else {
                assert(left[left[idx]] == idx);
                if(left[idx] < idx) {
                    // time travel path tracking!
                    frontier_index_t aux_idx = left[idx]; 
                    left[left[idx]] = invalid_index_;
                    left[idx] = two_edge_index_;
                    if(!remaining_edges_after_this_[bag_idx][aux_idx]) {
                        if(found_solution) {
                            return false;
                        }
                        left[aux_idx] = two_edge_index_;
                        found_solution = true;
                    } else {
                        cut_paths.push_back(aux_idx);
                    }
                } else {
                    // otherwise standard path tracking
                    left[left[idx]] = invalid_index_;
                    left[idx] = two_edge_index_;
                }
            }
        } else if(right[idx] > idx) {
            assert(right[right[idx]] == idx);
            if(left[idx] == no_edge_index_) {
                if(left[right[idx]] == no_edge_index_) {
                    left[idx] = right[idx];
                    left[right[idx]] = idx;
                    // standard path tracking
                } else if(left[right[idx]] == two_edge_index_) {
                    return false;
                } else if(left[right[idx]] == invalid_index_) {
                    left[idx] = invalid_index_;
                    left[right[idx]] = two_edge_index_;
                    // standard path tracking
                } else {
                    left[idx] = left[right[idx]];
                    left[left[right[idx]]] = idx;
                    left[right[idx]] = two_edge_index_;
                    // standard path tracking
                }
            } else if(left[idx] == two_edge_index_) {
                return false;
            } else if(left[idx] == invalid_index_) {
                if(left[right[idx]] == no_edge_index_) {
                    left[idx] = two_edge_index_;
                    left[right[idx]] = invalid_index_;
                    // standard path tracking
                } else if(left[right[idx]] == two_edge_index_) {
                    return false;
                } else if(left[right[idx]] == invalid_index_) {
                    if(found_solution) {
                        return false;
                    }
                    left[idx] = two_edge_index_;
                    left[right[idx]] = two_edge_index_;
                    found_solution = true;
                    // no path tracking
                } else {
                    if(left[right[idx]] < idx) {
                        // time travel path tracking!
                        frontier_index_t aux_idx = left[right[idx]]; 
                        left[idx] = two_edge_index_;
                        left[left[right[idx]]] = invalid_index_;
                        left[right[idx]] = two_edge_index_;
                        if(!remaining_edges_after_this_[bag_idx][aux_idx]) {
                            if(found_solution) {
                                return false;
                            }
                            left[aux_idx] = two_edge_index_;
                            found_solution = true;
                        } else {
                            cut_paths.push_back(aux_idx);
                        }
                    } else {
                        // otherwise standard path tracking
                        left[idx] = two_edge_index_;
                        left[left[right[idx]]] = invalid_index_;
                        left[right[idx]] = two_edge_index_;
                    }
                }
            } else {
                if(left[right[idx]] == no_edge_index_) {
                    left[left[idx]] = right[idx];
                    left[right[idx]] = left[idx];
                    left[idx] = two_edge_index_;
                    // standard path tracking
                } else if(left[right[idx]] == two_edge_index_) {
                    return false;
                } else if(left[right[idx]] == invalid_index_) {
                    if(left[idx] < idx) {
                        // time travel path tracking!
                        frontier_index_t aux_idx = left[idx]; 
                        left[left[idx]] = invalid_index_;
                        left[idx] = two_edge_index_;
                        left[right[idx]] = two_edge_index_;
                        if(!remaining_edges_after_this_[bag_idx][aux_idx]) {
                            if(found_solution) {
                                return false;
                            }
                            left[aux_idx] = two_edge_index_;
                            found_solution = true;
                        } else {
                            cut_paths.push_back(aux_idx);
                        }
                    } else {
                        // otherwise standard path tracking
                        left[left[idx]] = invalid_index_;
                        left[idx] = two_edge_index_;
                        left[right[idx]] = two_edge_index_;
                    }
                } else {
                    if(left[right[idx]] == idx) {
                        assert(left[idx] == right[idx]);
                        assert(right[left[idx]] == idx);
                        // found a loop
                        return false;
                    }
                    if(left[right[idx]] < idx && left[idx] < idx) {
                        frontier_index_t aux_idx = left[idx];
                        left[left[idx]] = left[right[idx]];
                        left[left[right[idx]]] = left[idx];
                        left[idx] = two_edge_index_;
                        left[right[idx]] = two_edge_index_;
                        // time travel path tracking!
                        if(!remaining_edges_after_this_[bag_idx][aux_idx] && !remaining_edges_after_this_[bag_idx][left[aux_idx]]) {
                            if(found_solution) {
                                return false;
                            }
                            left[left[aux_idx]] = two_edge_index_;
                            left[aux_idx] = two_edge_index_;
                            found_solution = true;
                        } else if(!remaining_edges_after_this_[bag_idx][aux_idx] && remaining_edges_after_this_[bag_idx][left[aux_idx]]) {
                            left[left[aux_idx]] = invalid_index_;
                            cut_paths.push_back(left[aux_idx]);
                            left[aux_idx] = two_edge_index_;
                        } else if(remaining_edges_after_this_[bag_idx][aux_idx] && !remaining_edges_after_this_[bag_idx][left[aux_idx]]) {
                            left[left[aux_idx]] = two_edge_index_;
                            left[aux_idx] = invalid_index_;
                            cut_paths.push_back(aux_idx);
                        } else {
                            paths.push_back(aux_idx);
                            paths.push_back(left[aux_idx]);
                        }
                    } else {
                        left[left[idx]] = left[right[idx]];
                        left[left[right[idx]]] = left[idx];
                        left[idx] = two_edge_index_;
                        left[right[idx]] = two_edge_index_;
                        // standard path tracking
                    }
                }
            }
        }
        // standard path tracking
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
    return finalizeMerge(left, bag_idx, found_solution, cut_paths, paths, left_result, right_result);
}

void NautyPathwidthSearch::advance(Frontier& frontier, size_t bag_idx) {
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
}

sparsegraph NautyPathwidthSearch::construct_sparsegraph(Frontier const& frontier, size_t last_idx) {
    if(last_idx == size_t(-1)) {
        // LEAF
        // just take full graph
        return graph_.to_canon_nauty(false);
    }
    size_t bag_idx = decomposition_[last_idx].parent;
    sparsegraph sg;
    SG_INIT(sg);
    auto const& base_sg = sparsegraph_after_this_[last_idx];
    sg.v = (edge_t *)malloc(sizeof(edge_t)*(base_sg.vlen + base_sg.elen) + sizeof(degree_t)*base_sg.dlen);
    sg.d = (degree_t *)(sg.v + base_sg.vlen);
    sg.e = (edge_t *)(sg.d + base_sg.vlen);
    sg.nv = base_sg.nv;
    sg.nde = base_sg.nde;
    sg.vlen = base_sg.vlen;
    sg.dlen = base_sg.vlen;
    sg.elen = base_sg.elen;
    std::memcpy(sg.v, base_sg.v, sizeof(edge_t)*(base_sg.vlen + base_sg.elen) + sizeof(degree_t)*base_sg.dlen);
    // std::cout << sg.nv << " " << sg.nde << std::endl;
    // std::cout << sg.vlen << " " << sg.dlen << " " << sg.elen << std::endl;
    // size_t e_count = 0;
    // for(size_t v = 0; v < sg.nv; v++) {
    //     std::vector<int> neighs(sg.nv, 0);
    //     for(size_t i = sg.v[v]; i < sg.v[v] + sg.d[v]; i++) {
    //         neighs[sg.e[i]] = 1;
    //         e_count++;
    //     }
    //     for(size_t i = 0; i < sg.nv; i++) {
    //         std::cout << neighs[i];
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << std::endl;
    int *lab = (int *)malloc(sizeof(int)*sg.nv);
    int *ptn = (int *)malloc(sizeof(int)*sg.nv);
    int *orbits = (int *)malloc(sizeof(int)*sg.nv);
    std::vector<frontier_index_t> two_edge;
    std::vector<frontier_index_t> no_edge;
    std::vector<frontier_index_t> invalid;
    std::vector<frontier_index_t> active;
    size_t offset = decomposition_[last_idx].bag.size();
    // now sg is just a copy of base_sg
    // we add the edges of the frontier 
    for(frontier_index_t idx = 0; idx < frontier.size(); idx++) {
        auto v = bag_local_vertex_map_[bag_idx][idx];
        auto old_v_idx = bag_local_idx_map_[last_idx][v];
        if(old_v_idx == invalid_index_) {
            assert(frontier[idx] == no_edge_index_);
            continue;
        }
        if(frontier[idx] == no_edge_index_) {
            // dont do anything
            no_edge.push_back(old_v_idx);
        } else if(frontier[idx] == two_edge_index_) {
            // remove adjacent
            two_edge.push_back(old_v_idx);
        } else if(frontier[idx] == invalid_index_) {
            // single connect fake node to real node
            assert(sg.d[old_v_idx] == 0);
            assert(sg.e[sg.v[old_v_idx]] == 0);
            sg.e[sg.v[old_v_idx]] = old_v_idx + offset;
            sg.d[old_v_idx]++;
            assert(sg.d[old_v_idx + offset] > 0);
            assert(sg.e[sg.v[old_v_idx + offset] + sg.d[old_v_idx + offset]] == 0);
            sg.e[sg.v[old_v_idx + offset] + sg.d[old_v_idx + offset]] = old_v_idx;
            sg.d[old_v_idx + offset]++;
            sg.nde += 2;
            invalid.push_back(old_v_idx);
        } else if(idx < frontier[idx]) {
            // add the edges
            auto w = bag_local_vertex_map_[bag_idx][frontier[idx]];
            auto old_w_idx = bag_local_idx_map_[last_idx][w];
            assert(old_w_idx != invalid_index_);
            sg.nde += 4;
            assert(sg.e[sg.v[old_v_idx] + sg.d[old_v_idx]] == 0);
            sg.e[sg.v[old_v_idx] + sg.d[old_v_idx]++] = old_w_idx + offset;
            assert(sg.e[sg.v[old_v_idx] + sg.d[old_v_idx]] == 0);
            sg.e[sg.v[old_v_idx] + sg.d[old_v_idx]++] = old_v_idx + offset;
            assert(sg.e[sg.v[old_v_idx + offset] + sg.d[old_v_idx + offset]] == 0);
            assert(sg.e[sg.v[old_w_idx + offset] + sg.d[old_w_idx + offset]] == 0);
            sg.e[sg.v[old_v_idx + offset] + sg.d[old_v_idx + offset]++] = old_v_idx;
            sg.e[sg.v[old_w_idx + offset] + sg.d[old_w_idx + offset]++] = old_v_idx;
            active.push_back(old_v_idx);
            active.push_back(old_w_idx);
        }
    }
    size_t found = 0;
    for(auto idx : no_edge) {
        lab[found] = idx;
        ptn[found++] = 1;
    }
    if(found > 0) {
        ptn[found - 1] = 0;
    }
    for(auto idx : two_edge) {
        lab[found] = idx;
        ptn[found++] = 1;
    }
    if(found > 0) {
        ptn[found - 1] = 0;
    }
    for(auto idx : invalid) {
        lab[found] = idx;
        ptn[found++] = 1;
    }
    if(found > 0) {
        ptn[found - 1] = 0;
    }
    for(auto idx : active) {
        lab[found] = idx;
        ptn[found++] = 1;
    }
    if(found > 0) {
        ptn[found - 1] = 0;
    }
    for(size_t j = 0; j < offset; j++) {
        if(remaining_edges_after_this_[last_idx][j] == 0 
            && (is_all_pair_ 
                || (bag_local_idx_map_[last_idx][terminals_[0]] != j 
                    && bag_local_idx_map_[last_idx][terminals_[1]] != j))) {
            lab[found] = j;
            ptn[found++] = 0;
        }
    }
    assert(found == offset);
    for(size_t j = offset; j < sg.nv; j++) {
        lab[j] = j;
        ptn[j] = 1;
    }
    ptn[offset - 1] = 0;
    // remove the edges that are incident to two_edge vertices
    for(auto old_v_idx : two_edge) {
        auto actual = old_v_idx + offset;
        // for all neighbors remove old_v_idx
        for(size_t i = 0; i < sg.d[actual]; i++) {
            auto neigh = sg.e[sg.v[actual] + i];
            for(size_t j = 0; j < sg.d[neigh]; j++) {
                if(sg.e[sg.v[neigh] + j] == actual) {
                    std::swap(sg.e[sg.v[neigh] + j], sg.e[sg.v[neigh] + sg.d[neigh] - 1]);
                    sg.e[sg.v[neigh] + --sg.d[neigh]] = 0;
                    break;
                }
            }
            sg.e[sg.v[actual] + i] = 0;
        }
        sg.nde -= 2*sg.d[actual];
        sg.d[actual] = 0;
        assert(sg.d[old_v_idx] == 0);
        assert(sg.e[sg.v[old_v_idx]] == 0);
        sg.e[sg.v[old_v_idx]] = actual;
        sg.d[old_v_idx]++;
        assert(sg.e[sg.v[actual]] == 0);
        sg.e[sg.v[actual]] = old_v_idx;
        sg.d[actual]++;
        sg.nde += 2;
    }
    // std::cout << sg.nv << " " << sg.nde << std::endl;
    // std::cout << sg.vlen << " " << sg.dlen << " " << sg.elen << std::endl;
    // size_t e_count = 0;
    // for(size_t v = 0; v < sg.nv; v++) {
    //     std::vector<int> neighs(sg.nv, 0);
    //     for(size_t i = sg.v[v]; i < sg.v[v] + sg.d[v]; i++) {
    //         neighs[sg.e[i]] = 1;
    //         e_count++;
    //     }
    //     for(size_t i = 0; i < sg.nv; i++) {
    //         std::cout << neighs[i];
    //     }
    //     std::cout << std::endl;
    // }
    // for(size_t v = 0; v < sg.nv; v++) {
    //     std::cout << lab[v] << " ";
    // }
    // std::cout << std::endl;
    // for(size_t v = 0; v < sg.nv; v++) {
    //     std::cout << ptn[v] << " ";
    // }
    // std::cout << std::endl;
    // assert(e_count == sg.nde);
    DEFAULTOPTIONS_SPARSEGRAPH(options);
    options.getcanon = true;
    options.defaultptn = false;
    // options.schreier = true;
    // options.tc_level = 1000;
    statsblk stats;
    SG_DECL(canon_sg);
    canon_sg.v = (edge_t *)malloc(sizeof(edge_t)*(sg.nv + sg.nde) + sizeof(degree_t)*sg.nv);
    canon_sg.d = (degree_t *)(canon_sg.v + sg.nv);
    canon_sg.e = (edge_t *)(canon_sg.d + sg.nv);
    canon_sg.nv = sg.nv;
    canon_sg.nde = sg.nde;
    canon_sg.vlen = sg.nv;
    canon_sg.dlen = sg.nv;
    canon_sg.elen = sg.nde;
    sparsenauty(&sg,lab,ptn,orbits,&options,&stats,&canon_sg);
    sortlists_sg(&canon_sg);
    free(sg.v);
    free(lab);
    free(ptn);
    free(orbits);
    return canon_sg;
}


void NautyPathwidthSearch::print_stats() {
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
