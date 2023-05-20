#include "search.h"
#include <queue>
#include <limits>

clhasher hasher__(UINT64_C(0x23a23cf5033c3c81),UINT64_C(0xb3816f6a2c68e530));

Search::Search(Graph& input) :  enable_dag_(true),
                                max_length_(input.max_length_),
                                terminals_(input.terminals_),
                                neighbors_(input.neighbors_.size()),
                                adjacency_(input.adjacency_.size(), std::vector<std::vector<std::pair<Edge_length, Edge_weight>>>(input.adjacency_.size())),
                                invalid_(std::numeric_limits<Edge_length>::max() - max_length_ - 1),
                                distance_to_goal_(adjacency_.size(), invalid_),
                                distance_(adjacency_.size(), std::vector<Edge_length>(adjacency_.size(), invalid_)),
                                exclude_(adjacency_.size()),
                                visited_(adjacency_.size(), false),
                                cache_( 
                                    adjacency_.size(), 
                                    std::vector<std::unordered_map<CacheKey, std::pair<Edge_length, std::vector<Edge_weight>>, vector_hash>>(
                                        adjacency_.size(), 
                                        std::unordered_map<CacheKey, std::pair<Edge_length, std::vector<Edge_weight>>, vector_hash>()
                                    )
                                ),
                                ap_disc_(adjacency_.size()),
                                ap_low_(adjacency_.size()),
                                ap_visited_(adjacency_.size())  {
    assert(terminals_.size() == 2);
    for(Vertex v = 0; v < adjacency_.size();v++) {
        // fill neighbors
        neighbors_[v] = std::vector<Vertex>(input.neighbors_[v].begin(), input.neighbors_[v].end());
        for(Vertex w : neighbors_[v]) {
            adjacency_[v][w] = std::vector<std::pair<Edge_length, Edge_weight>>(input.adjacency_[v][w].begin(), input.adjacency_[v][w].end());
        }
        for(Vertex w : input.exclusion_classes_[input.exclude_[v]]) {
            if(w != v) {
                enable_dag_ = false;
                exclude_[v].push_back(w);
            }
        }
    }
    dijkstra(terminals_[1], distance_to_goal_, max_length_);
    for(Vertex v = 0; v < adjacency_.size(); v++) {
        dijkstra(v, distance_[v], max_length_);
    }
}

std::vector<Edge_weight> Search::search(Vertex start, Edge_length budget) {
    edges++;
    if(edges % 1000000 == 0) {
        print_stats();
    }
    if(start == terminals_[1]) {
        // we have one path of length zero.
        return {1};
    }
    // visit the current start
    assert(!visited_[start]);
    visited_[start] = true;
    for(Vertex excluded : exclude_[start]) {
        assert(!visited_[excluded]);
        visited_[excluded] = true;
    }
    // update distance based on new budged
    distance_to_goal_ = std::vector<Edge_length>(adjacency_.size(), invalid_);
    pruning_dijkstra(terminals_[1], start, distance_to_goal_, budget);

    // first check if there is something simple we can do
    // only process articulation points and cache otherwise
    // simple means: 
    //      * there is no edge (can only happen when we have exclusion constraints)
    //      * there are only dag edges
    //      * there is only one edge
    std::vector<Vertex> poss_non_dag, poss_dag;
    for(auto v : neighbors(start)) {
        if(budget > distance_to_goal_[v] + adjacency_[start][v].begin()->first) {
            poss_non_dag.push_back(v);
        }
        if(budget == distance_to_goal_[v] + adjacency_[start][v].begin()->first) {
            poss_dag.push_back(v);
        }
    }
    // there is no edge
    if(poss_non_dag.size() + poss_dag.size() == 0) {
        visited_[start] = false;
        for(Vertex excluded : exclude_[start]) {
            visited_[excluded] = false;
        }
        return {};
    }
    // there are only dag edges
    if(poss_non_dag.size() == 0 && enable_dag_) {
        // dont cache and do everything here
        std::vector<Edge_weight> ret(budget + 1, 0);
        for(auto v : poss_dag) {
            dags++;
            std::vector<Edge_weight> tmp = dag_search(v, budget - adjacency_[start][v].begin()->first);
            for(size_t i = 0; i < tmp.size(); i++) {
                for(auto &[length, weight] : adjacency_[start][v]) {
                    if(length + i > budget) {
                        break;
                    }
                    ret[length + i] += weight*tmp[i];
                }
            }
        }
        visited_[start] = false;
        for(Vertex excluded : exclude_[start]) {
            visited_[excluded] = false;
        }
        return ret;
    }
    // there is only one edge
    if(poss_dag.size() + poss_non_dag.size() == 1) {
        propagations++;
        std::vector<Edge_weight> tmp;
        Vertex v;
        if(poss_dag.size() > 0) {
            v = poss_dag[0];
        } else {
            v = poss_non_dag[0];
        }
        tmp = search(v, budget - adjacency_[start][v].begin()->first);
        std::vector<Edge_weight> ret(budget + 1, 0);
        for(size_t i = 0; i < tmp.size(); i++) {
            for(auto &[length, weight] : adjacency_[start][v]) {
                if(length + i > budget) {
                    break;
                }
                ret[length + i] += weight*tmp[i];
            }
        }
        visited_[start] = false;
        for(Vertex excluded : exclude_[start]) {
            visited_[excluded] = false;
        }
        return ret;
    }
    // now do more complicated stuff 
    // based on articulaiton points we:
    //      * prune unreachable parts of the graph
    //      * can split "start - G_1 - ap - G_2 - t_1"
    //        into "start - G_1 - ap" and "ap - G_2 - t_1"
    //        and solve independently
    prune_articulation(start);
    if(!ap_components_.empty() && enable_dag_) {
        splits += ap_components_.size() - 1;
        assert(ap_components_.size() >= 2);
        // we found a split and should process the ap components separately.
        
        // remember the old configuration
        Vertex prev_goal = terminals_[1];
        auto prev_distance_to_goal = distance_to_goal_;
        auto prev_visited = visited_;
        auto cur_ap_components = ap_components_;
        auto cur_ap_start_goal = ap_start_goal_;

        for(Vertex v : neighbors(start)) {
            if(!visited_[v]) {
                prev_distance_to_goal[start] = std::min(prev_distance_to_goal[start], Edge_length(prev_distance_to_goal[v] + adjacency_[start][v].begin()->first));
            }
        }
        // we start with one path of length 0
        std::vector<Edge_weight> result = {1};
        Edge_weight used_budget = 0;
        for(int i = cur_ap_components.size() - 1; i >= 0; i--) {
            terminals_[1] = cur_ap_start_goal[i].second;
            std::fill(visited_.begin(), visited_.end(), 1);
            for(Vertex v : cur_ap_components[i]) {
                visited_[v] = false;
            }
            assert(budget >= used_budget + prev_distance_to_goal[cur_ap_start_goal[i].second]);
            auto partial_res = search(cur_ap_start_goal[i].first, budget - used_budget - prev_distance_to_goal[cur_ap_start_goal[i].second]);
            // incorporate the partial result
            std::vector<Edge_weight> new_result(budget + 1 - prev_distance_to_goal[cur_ap_start_goal[i].second], 0);
            for(Edge_length prev_length = 0; prev_length < result.size(); prev_length++) {
                for(Edge_length cur_length = 0; cur_length + prev_length < new_result.size(); cur_length++) {
                    new_result[prev_length + cur_length] += result[prev_length]*partial_res[cur_length];
                }
            }
            result = new_result;
            used_budget += (prev_distance_to_goal[cur_ap_start_goal[i].first] - prev_distance_to_goal[cur_ap_start_goal[i].second]);
        }

        terminals_[1] = prev_goal;
        visited_ = prev_visited;
        visited_[start] = false;
        for(Vertex excluded : exclude_[start]) {
            visited_[excluded] = false;
        }
        // FIXME: Do we want to cache here?
        // only makes sense if we also check the cache before going in here
        return result;
    }
    
    auto cached_result = cache_[start][terminals_[1]].find(ap_visited_);
    if(cached_result != cache_[start][terminals_[1]].end()) {
        if(cached_result->second.first >= budget) {
            pos_hits++;
            visited_[start] = false;
            for(Vertex excluded : exclude_[start]) {
                visited_[excluded] = false;
            }
            std::vector<Edge_weight> ret(cached_result->second.second.begin(), cached_result->second.second.begin() + budget + 1);
            return ret;
        }
    }
    neg_hits++;
    auto cached_position = cache_[start][terminals_[1]].insert(std::make_pair(ap_visited_, std::make_pair(budget, std::vector<Edge_weight>())));
    std::vector<Edge_weight> ret(budget + 1, 0);
    for(auto v : poss_dag) {
        std::vector<Edge_weight> tmp;
        if(enable_dag_) {
            dags++;
            tmp = dag_search(v, budget - adjacency_[start][v].begin()->first);
        } else {
            tmp = search(v, budget - adjacency_[start][v].begin()->first);
        }
        for(size_t i = 0; i < tmp.size(); i++) {
            for(auto &[length, weight] : adjacency_[start][v]) {
                if(length + i > budget) {
                    break;
                }
                ret[length + i] += weight*tmp[i];
            }
        }
    }for(auto v : poss_non_dag) {
        std::vector<Edge_weight> tmp = search(v, budget - adjacency_[start][v].begin()->first);
        for(size_t i = 0; i < tmp.size(); i++) {
            for(auto &[length, weight] : adjacency_[start][v]) {
                if(length + i > budget) {
                    break;
                }
                ret[length + i] += weight*tmp[i];
            }
        }
    }
    cached_position.first->second.second = ret;
    visited_[start] = false;
    for(Vertex excluded : exclude_[start]) {
        visited_[excluded] = false;
    }
    return ret;
}

std::vector<Edge_weight> Search::dag_search(Vertex start, Edge_length budget) {
    std::priority_queue<std::pair<Edge_length, Vertex>> biggest_first_queue;
    std::vector<Edge_weight> dp(adjacency_.size(), 0);
    std::vector<char> in_queue(adjacency_.size(), false);
    dp[start] = 1;
    biggest_first_queue.push(std::make_pair(distance_to_goal_[start], start));
    while(!biggest_first_queue.empty()) {
        auto [cur_cost, cur_vertex] = biggest_first_queue.top();
        biggest_first_queue.pop();
        for(auto &w : neighbors(cur_vertex)) {
            if(distance_to_goal_[w] == invalid_) {
                continue;
            }
            auto edge_cost = adjacency_[w][cur_vertex].begin()->first;
            // this edge can be part of a shorted path in direction w -> cur_vertex
            if(distance_to_goal_[w] == distance_to_goal_[cur_vertex] + edge_cost) {
                auto factor = adjacency_[w][cur_vertex].begin()->second;
                dp[cur_vertex] += factor*dp[w];
            // this edge can be part of a shorted path in direction cur_vertex -> w
            } else if(distance_to_goal_[w] + edge_cost == distance_to_goal_[cur_vertex]) {
                if(!in_queue[w]) {
                    biggest_first_queue.push(std::make_pair(distance_to_goal_[w], w));
                    in_queue[w] = true;
                }
            }
        }
    }
    std::vector<Edge_weight> result(budget + 1, 0);
    result[budget] = dp[terminals_[1]];
    assert(result[budget] >= 1);
    return result;
}

void Search::prune_articulation(Vertex start) {
    ap_start_goal_.clear();
    ap_components_.clear();
    std::fill(ap_visited_.data_, ap_visited_.data_ + ap_visited_.chunks_, 0);
    std::fill(ap_disc_.begin(), ap_disc_.end(), 0);
    std::fill(ap_low_.begin(), ap_low_.end(), 0);
    last_ap_ = terminals_[1];
    int time = 0;
    ap_util(start, ap_visited_, ap_disc_, ap_low_, time, -1, start);
    if(last_ap_ != terminals_[1]) {
        ap_start_goal_.push_back(std::make_pair(start, last_ap_));
        ap_components_.push_back({start});
        component_util(start, ap_disc_);
    }
}

bool Search::ap_util(Vertex u, Bitset& visited, std::vector<Vertex>& disc, std::vector<Vertex>& low, int& time, int parent, Vertex start) {
    // We do not need to check whether the root is an articulation point
    // since if it is, then the goal can only be in one of the induced components
    // for all the other components we cannot enter them since we cannot reach the goal from them
    // but this means that we have already pruned them using dijkstra
    
    // Mark the current node as visited
    visited.SetTrue(u);

    bool found_elsewhere = false;
 
    // Initialize discovery time and low value
    disc[u] = low[u] = ++time;
 
    // Go through all vertices adjacent to this
    for (auto v : neighbors(u)) {
        // If v is not visited yet, then make it a child of u
        // in DFS tree and recur for it            
        if(v != start && (distance_to_goal_[v] == invalid_ || visited_[v])) {
            continue;
        }
        if (!visited.Get(v)) {
            bool found_here = ap_util(v, visited, disc, low, time, u, start);
            found_elsewhere |= found_here;
 
            // Check if the subtree rooted with v has
            // a connection to one of the ancestors of u
            low[u] = std::min(low[u], low[v]);
 
            // If u is not root and low value of one of
            // its child is more than discovery value of u.
            if (parent != -1 && low[v] >= disc[u]) {
                // AP
                if(!found_here) {
                    auto tmp = distance_to_goal_[u];
                    distance_to_goal_[u] = invalid_;
                    // prune the rest
                    distance_to_goal_[v] = invalid_;
                    visited.SetFalse(v);
                    prune_util(v);
                    distance_to_goal_[u] = tmp;
                } else {
                    auto tmp = disc[u];
                    disc[u] = invalid_;
                    ap_start_goal_.push_back(std::make_pair(u, last_ap_));
                    ap_components_.push_back({u, v});
                    component_util(v, disc);
                    disc[u] = tmp;
                    disc[v] = invalid_;
                    last_ap_ = u;
                }
            }
        } else if (v != parent) {
            low[u] = std::min(low[u], disc[v]);
        }
    }
    return found_elsewhere || u == terminals_[1];
}

void Search::prune_util(Vertex u) {
    for (auto v : neighbors(u)) {
        // If v is not visited yet, then make it a child of u
        // in DFS tree and recur for it
        if (distance_to_goal_[v] != invalid_ && !visited_[v]) {
            distance_to_goal_[v] = invalid_;
            ap_visited_.SetFalse(v);
            prune_util(v);
        }
    }
}
void Search::component_util(Vertex u, std::vector<Vertex>& disc) {
    for (auto v : neighbors(u)) {
        // If v is not visited yet, then make it a child of u
        // in DFS tree and recur for it
        if (distance_to_goal_[v] != invalid_ && !visited_[v] && disc[v] != invalid_) {
            disc[v] = invalid_;
            ap_components_.back().push_back(v);
            component_util(v, disc);
        }
    }
}
 

void Search::dijkstra(Vertex start, std::vector<Edge_length>& distance, Edge_length budget) {
    DijkstraQueue queue;
    queue.push(std::make_pair(0, start));
    distance[start] = 0;
    while(!queue.empty()) {
        auto [cur_cost, cur_vertex] = queue.top();
        queue.pop();
        if(cur_cost > distance[cur_vertex]) {
            continue;
        }
        for(auto &w : neighbors(cur_vertex)) {
            if(cur_cost + 1 >= distance[w] || visited_[w]) {
                continue;
            }
            Edge_length min_cost = adjacency_[cur_vertex][w].begin()->first;
            if(cur_cost + min_cost < distance[w] && cur_cost + min_cost <= budget) {
                distance[w] = min_cost + cur_cost;
                queue.push(std::make_pair(cur_cost + min_cost, w));
            }
        }
    }
}
void Search::pruning_dijkstra(Vertex start, Vertex prune, std::vector<Edge_length>& distance, Edge_length budget) {
    std::deque<Vertex> queue;
    queue.push_back(start);
    distance[start] = 0;
    while(!queue.empty()) {
        auto cur_vertex = queue.front();
        auto cur_cost = distance[cur_vertex];
        queue.pop_front();
        for(auto &w : neighbors(cur_vertex)) {
            if(cur_cost + 1 >= distance[w] || visited_[w]) {
                continue;
            }
            Edge_length min_cost = adjacency_[cur_vertex][w].begin()->first;
            if(cur_cost + min_cost < distance[w] && cur_cost + min_cost + distance_[prune][w] <= budget) {
                distance[w] = min_cost + cur_cost;
                queue.push_back(w);
            }
        }
    }
}

void Search::print_stats() {
    std::cerr << "Cache hit rate: " << 100*pos_hits/(double)(pos_hits + neg_hits) << "% (" << pos_hits << "/" << pos_hits + neg_hits << ")" << std::endl;
    std::cerr << "#DAG searches: " << dags << " #Splits: " << splits << std::endl;
    std::cerr << "#Edges: " << edges << " #Propagations: " << propagations << std::endl;
}
