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
                                cache_(adjacency_.size(), std::unordered_map<CacheKey, std::pair<Edge_length, std::vector<Edge_weight>>>())  {
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
    if(edges % 100000 == 0) {
        print_stats();
    }
    if(start == terminals_[1]) {
        // we have one path of length zero.
        return {1};
    }
    assert(!visited_[start]);
    visited_[start] = true;
    std::vector<Vertex> restore;
    for(Vertex excluded : exclude_[start]) {
        if(!visited_[excluded]) {
            visited_[excluded] = true;
            restore.push_back(excluded);
        }
    }
    distance_to_goal_ = std::vector<Edge_length>(adjacency_.size(), invalid_);
    pruning_dijkstra(terminals_[1], start, distance_to_goal_, budget);
    prune_articulation(start);
    // only cache if there is more than one edge we can take
    std::vector<Vertex> poss_non_dag, poss_dag;
    for(auto v : neighbors(start)) {
        if(budget > distance_to_goal_[v] + adjacency_[start][v].begin()->first) {
            poss_non_dag.push_back(v);
        }
        if(budget == distance_to_goal_[v] + adjacency_[start][v].begin()->first) {
            poss_dag.push_back(v);
        }
    }
    if(poss_non_dag.size() + poss_dag.size() == 0) {
        visited_[start] = false;
        for(Vertex excluded : restore) {
            visited_[excluded] = false;
        }
        return {};
    }
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
        for(Vertex excluded : restore) {
            visited_[excluded] = false;
        }
        return ret;
    }
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
        for(Vertex excluded : restore) {
            visited_[excluded] = false;
        }
        return ret;
    }
    auto cached_result = cache_[start].find(distance_to_goal_);
    if(cached_result != cache_[start].end()) {
        if(cached_result->second.first >= budget) {
            pos_hits++;
            visited_[start] = false;
            for(Vertex excluded : restore) {
                visited_[excluded] = false;
            }
            std::vector<Edge_weight> ret(cached_result->second.second.begin(), cached_result->second.second.begin() + budget + 1);
            return ret;
        }
    }
    neg_hits++;
    auto cached_position = cache_[start].insert(std::make_pair(distance_to_goal_, std::make_pair(budget, std::vector<Edge_weight>())));
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
    for(Vertex excluded : restore) {
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
    std::vector<Vertex> disc(adjacency_.size(), 0);
    std::vector<Vertex> low(adjacency_.size(), 0);
    std::vector<char> visited(adjacency_.size(), false);
    int time = 0;
    ap_util(start, visited, disc, low, time, -1, start);
}

bool Search::ap_util(Vertex u, std::vector<char>& visited, std::vector<Vertex>& disc, std::vector<Vertex>& low, int& time, int parent, Vertex start) {
    // Count of children in DFS Tree
    int children = 0;
 
    // Mark the current node as visited
    visited[u] = true;

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
        if (!visited[v]) {
            children++;
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
                    prune_util(v);
                    distance_to_goal_[u] = tmp;
                }
            }
        } else if (v != parent) {
            low[u] = std::min(low[u], disc[v]);
        }
    }
 
    // If u is root of DFS tree and has two or more children.
    if (parent == -1 && children > 1) {
        // AP
    }
    return found_elsewhere || u == terminals_[1];
}

void Search::prune_util(Vertex u) {
    for (auto v : neighbors(u)) {
        // If v is not visited yet, then make it a child of u
        // in DFS tree and recur for it
        if (distance_to_goal_[v] != invalid_ && !visited_[v]) {
            distance_to_goal_[v] = invalid_;
            prune_util(v);
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
            if(cur_cost + min_cost < distance[w] && cur_cost + min_cost + distance_[prune][w] <= budget) {
                distance[w] = min_cost + cur_cost;
                queue.push(std::make_pair(cur_cost + min_cost, w));
            }
        }
    }
}

void Search::print_stats() {
    std::cerr << "Cache hit rate: " << 100*pos_hits/(double)(pos_hits + neg_hits) << "% (" << pos_hits << "/" << pos_hits + neg_hits << ")" << std::endl;
    std::cerr << "#DAG searches: " << dags << std::endl;
    std::cerr << "#Edges: " << edges << " #Propagations: " << propagations << std::endl;
}
