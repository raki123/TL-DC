#include "treewidth_search.h"
#include <algorithm>
namespace fpc {

TreewidthSearch::TreewidthSearch(Graph& input, std::vector<std::pair<Edge, std::vector<vertex_t>>> path_decomposition, size_t nthreads) 
    :   nthreads_(nthreads),
        graph_(input),
        max_length_(graph_.max_length()),
        is_all_pair_(graph_.is_all_pair()),
        terminals_(graph_.terminals()),
        path_decomposition_(path_decomposition),
        remaining_edges_after_this_(path_decomposition_.size()),
        bag_local_idx_map_(
            path_decomposition_.size(),
            std::vector<frontier_index_t>(graph_.adjacency_.size(), invalid_index_)),
        bag_local_vertex_map_(path_decomposition_.size()),
        bag_local_distance_(path_decomposition_.size()),
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
    std::vector<vertex_t> degree;
    for(size_t v = 0; v < graph_.adjacency_.size(); v++) {
        degree.push_back(graph_.neighbors(v).size());
    }
    std::vector<vertex_t> last_remaining;
    std::vector<frontier_index_t> last_idx;
    if(!is_all_pair_) {
        for(auto &[_, bag] : path_decomposition_) {
            if(std::find(bag.begin(), bag.end(), terminals_[0]) == bag.end()) {
                bag.push_back(terminals_[0]);
            }
            if(std::find(bag.begin(), bag.end(), terminals_[1]) == bag.end()) {
                bag.push_back(terminals_[1]);
            }
        }
    }
    for(size_t bag_idx = 0; bag_idx < path_decomposition_.size(); bag_idx++) {
        auto edge = path_decomposition_[bag_idx].first;
        auto &bag = path_decomposition_[bag_idx].second;
        auto &idx = bag_local_idx_map_[bag_idx];
        auto &vertex = bag_local_vertex_map_[bag_idx];
        auto &remaining = remaining_edges_after_this_[bag_idx];
        remaining.resize(bag.size());
        for(size_t i = 0; i < bag.size(); i++) {
            if(degree[bag[i]] == 0) {
                std::swap(bag[i], bag.back());
                bag.pop_back();
                i--;
            } else if(is_all_pair_ || (bag[i] != terminals_[0] && bag[i] != terminals_[1])) {
                if(degree[bag[i]] == graph_.neighbors(bag[i]).size() && bag[i] != edge.first && bag[i] != edge.second) {
                    std::swap(bag[i], bag.back());
                    bag.pop_back();
                    i--;
                }
            }
        }
        std::sort(bag.begin(), bag.end());
        frontier_index_t cur_idx = 0;
        for(auto v : bag) {
            assert(cur_idx != invalid_index_);
            idx[v] = cur_idx++;
            vertex.push_back(v);
            if(seen[v]) {
                assert(v < idx.size());
                frontier_index_t new_idx = idx[v];
                assert(new_idx != invalid_index_);
                assert(v < last_idx.size());
                frontier_index_t old_idx = last_idx[v];
                assert(old_idx != invalid_index_);
                assert(new_idx < remaining.size());
                assert(old_idx < last_remaining.size());
                remaining[new_idx] = last_remaining[old_idx];
            } else {
                seen[v] = true;
                remaining[idx[v]] = graph_.neighbors(v).size();
            }
        }
        assert(degree[edge.first] != 0);
        assert(remaining[idx[edge.first]] != 0);
        remaining[idx[edge.first]]--;
        degree[edge.first]--;
        assert(degree[edge.second] != 0);
        assert(remaining[idx[edge.second]] != 0);
        remaining[idx[edge.second]]--;
        degree[edge.second]--;
        assert(degree[edge.first] == remaining[idx[edge.first]]);
        assert(degree[edge.second] == remaining[idx[edge.second]]);
        last_remaining = remaining;
        last_idx = idx;
    }
    for(vertex_t v = 0; v < graph_.adjacency_.size(); v++) {
        if(degree[v] != 0) {
            std::cerr << v << std::endl;
        }
        assert(degree[v] == 0);
    }
    std::set<vertex_t> empty;
    for(size_t bag_idx = 0; bag_idx < path_decomposition_.size(); bag_idx++) {
        auto edge = path_decomposition_[bag_idx].first;
        auto &bag = path_decomposition_[bag_idx].second;
        auto &idx = bag_local_idx_map_[bag_idx];
        bag_local_distance_[bag_idx].resize(bag.size());
        graph_.remove_edge(edge);
        for(auto v : bag) {
            std::vector<Edge_length> distance(graph_.adjacency_.size(), invalid_distance_);
            graph_.dijkstra(v, distance, empty);
            std::vector<Edge_length> local_distance(bag.size(), invalid_distance_);
            for(auto other : bag) {
                if(other != v) {
                    local_distance[idx[other]] = distance[other];
                }
            }
            bag_local_distance_[bag_idx][idx[v]] = local_distance;
        }
    }

    Frontier initial_frontier(path_decomposition_[0].second.size(), no_edge_index_);
    if(!is_all_pair_) {
        assert(bag_local_idx_map_[0][terminals_[0]] != invalid_index_);
        assert(bag_local_idx_map_[0][terminals_[1]] != invalid_index_);
        initial_frontier[bag_local_idx_map_[0][terminals_[0]]] = invalid_index_;
        initial_frontier[bag_local_idx_map_[0][terminals_[1]]] = invalid_index_;
    }
    // cached vectors are {offset, results}
    std::vector<Edge_weight> initial_result = {0, 1};
    cache_[0][initial_frontier] = initial_result;
    includeSolutions(initial_frontier, 0, initial_result);
}

std::vector<Edge_weight> TreewidthSearch::search() {
    for(size_t bag_idx = 0; bag_idx < path_decomposition_.size(); bag_idx++) {
        #pragma omp parallel for /*default(none)*/ shared(bag_idx) shared(max_length_) shared(cache_) shared(path_decomposition_) shared(thread_local_result_) shared(pos_hits_) shared(neg_hits_) shared(bag_local_idx_map_) shared(bag_local_vertex_map_)
        for(size_t bucket = 0; bucket < cache_[bag_idx].bucket_count(); bucket++) {
            size_t thread_id = omp_get_thread_num();
            for(auto task_it = cache_[bag_idx].begin(bucket); task_it != cache_[bag_idx].end(bucket); ++task_it) {
                auto const& old_frontier = task_it->first;
                auto const& result = task_it->second;
                for(auto pp : {std::make_pair(true, false), std::make_pair(false, true)}) {
                    bool takeable = pp.first;
                    bool skippable = pp.second;
                    Frontier new_frontier = old_frontier;
                    size_t new_idx = bag_idx;
                    std::vector<Edge_weight> new_result = result;
                    edges_[thread_id]++;
                    propagations_[thread_id]--;
                    // propagate while only one of the two is possible
                    while(takeable ^ skippable && new_idx + 1 < path_decomposition_.size()) {
                        propagations_[thread_id]++;
                        // Edge edge = path_decomposition_[new_idx].first;
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
                        new_idx++;
                        if(new_idx < path_decomposition_.size()) {
                            // Edge edge = path_decomposition_[new_idx].first;
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
                    if(takeable && skippable && new_idx < path_decomposition_.size()) {
                        #pragma omp critical 
                        {
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
        }
        cache_[bag_idx] = {};
        std::cerr << bag_idx << std::endl;
        // std::cerr << results() << std::endl;
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
    Edge edge = path_decomposition_[bag_idx].first;
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
            frontier[v_idx] = no_edge_index_;
        }
    }
    if(w_remaining == 0) {
        assert(frontier[w_idx] != invalid_index_);
        assert(frontier[w_idx] != no_edge_index_);
        if(frontier[w_idx] != two_edge_index_) {
            frontier[frontier[w_idx]] = invalid_index_;
            frontier[w_idx] = no_edge_index_;
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
    Edge edge = path_decomposition_[bag_idx].first;
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
    if(cut_paths.size() == 2) {
        if(paths.size() == 0) {
            // we have to connect the cut paths
            return distance[cut_paths[0]][cut_paths[1]] + offset <= max_length_;
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
        for(auto v : paths) {
            Edge_length cur_min_dist = std::min(distance[v][cut_paths[0]], distance[v][cut_paths[1]]);
            for(auto other : paths) {
                if(frontier[v] != other) {
                    cur_min_dist = std::min(cur_min_dist, distance[v][other]);
                }
            }
            path_max_dist = std::max(path_max_dist, cur_min_dist);
            path_min_dist += cur_min_dist;
        }
        path_min_dist /= 2;
        if(path_min_dist > path_max_dist) {
            path_min_dist -= path_max_dist;
        } else {
            path_min_dist = 0;
        }
        if(path_min_dist + min_dist + offset >  max_length_) {
            return false;
        }
    }
    if(cut_paths.size() == 1) {
        if(paths.size() == 0) {
            // + 1 because we want to take another edge later
            // we dont have to do anything
            return offset <= max_length_;
        }
        // there is at least one other path that we have to incorporate
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
        for(auto v : paths) {
            Edge_length cur_min_dist = distance[v][cut_paths[0]];
            for(auto other : paths) {
                if(frontier[v] != other) {
                    cur_min_dist = std::min(cur_min_dist, distance[v][other]);
                }
            }
            path_max_dist = std::max(path_max_dist, cur_min_dist);
            path_min_dist += cur_min_dist;
        }
        path_min_dist /= 2;
        if(path_min_dist > path_max_dist) {
            path_min_dist -= path_max_dist;
        } else {
            path_min_dist = 0;
        }
        if(path_min_dist + min_dist + offset >  max_length_) {
            return false;
        }
    }
    if(cut_paths.size() == 0) {
        if(paths.size() <= 1) {
            // + 1 because we want to take another edge later
            // we dont have to do anything
            return offset + 1 <= max_length_;
        }
        // there is at least one other path that we have to incorporate
        size_t path_min_dist = 0;
        Edge_length path_max_dist = 0;
        for(auto v : paths) {
            Edge_length cur_min_dist = invalid_distance_;
            for(auto other : paths) {
                if(frontier[v] != other) {
                    cur_min_dist = std::min(cur_min_dist, distance[v][other]);
                }
            }
            path_max_dist = std::max(path_max_dist, cur_min_dist);
            path_min_dist += cur_min_dist;
        }
        path_min_dist /= 2;
        if(path_min_dist > path_max_dist) {
            path_min_dist -= path_max_dist;
        } else {
            path_min_dist = 0;
        }
        if(path_min_dist + offset > max_length_) {
            return false;
        }
    }
    return true;
}

void TreewidthSearch::take(Frontier& frontier, size_t bag_idx) {
    Edge edge = path_decomposition_[bag_idx].first;
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

void TreewidthSearch::advance(Frontier& frontier, size_t bag_idx) {
    static __thread Frontier old;
    old = frontier;
    auto &bag = path_decomposition_[bag_idx + 1].second;
    frontier.resize(bag.size());
    std::fill(frontier.begin(), frontier.end(), no_edge_index_);
    auto &old_idx = bag_local_idx_map_[bag_idx];
    auto &new_idx = bag_local_idx_map_[bag_idx + 1];
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
    std::cerr << "Cache hit rate: " << 100*pos_hits/(double)(pos_hits + neg_hits) << "% (" << pos_hits << "/" << pos_hits + neg_hits << ")" << std::endl;
    std::cerr << "#Edges: " << edges << " #Propagations: " << propagations << std::endl;
}


}
