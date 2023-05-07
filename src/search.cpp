#include "search.h"
#include <queue>
#include <limits>

clhasher hasher__(UINT64_C(0x23a23cf5033c3c81),UINT64_C(0xb3816f6a2c68e530));

Search::Search(Graph& input) :  max_length_(input.max_length_),
                                terminals_(input.terminals_),
                                neighbors_(input.neighbors_),
                                adjacency_(input.adjacency_),
                                invalid_(std::numeric_limits<Edge_length>::max() - max_length_ - 1),
                                distance_to_goal_(adjacency_.size(), invalid_),
                                distance_(adjacency_.size(), std::vector<Edge_length>(adjacency_.size(), invalid_)),
                                visited_(adjacency_.size(), false),
                                cache_(adjacency_.size(), std::unordered_map<CacheKey, std::pair<Edge_length, std::vector<Edge_weight>>>())  {
    assert(terminals_.size() == 2);
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
    visited_[start] = true;
    distance_to_goal_ = std::vector<Edge_length>(adjacency_.size(), invalid_);
    pruning_dijkstra(terminals_[1], start, distance_to_goal_, budget);
    prune_singleout(start);
    // only cache if there is more than one edge we can take
    std::vector<std::pair<Vertex, bool>> poss;
    for(auto v : neighbors(start)) {
        if(budget >= distance_to_goal_[v] + adjacency_[start][v].begin()->first) {
            poss.push_back(std::make_pair(v, budget == distance_to_goal_[v] + adjacency_[start][v].begin()->first));
        }
    }
    assert(poss.size() >= 1);
    if(poss.size() == 1) {
        propagations++;
        Vertex v = poss[0].first;
        assert(!poss[0].second);
        auto tmp = search(v, budget - adjacency_[start][v].begin()->first);
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
        return ret;
    }
    auto cached_result = cache_[start].find(distance_to_goal_);
    if(cached_result != cache_[start].end()) {
        if(cached_result->second.first >= budget) {
            pos_hits++;
            visited_[start] = false;
            std::vector<Edge_weight> ret(cached_result->second.second.begin(), cached_result->second.second.begin() + budget + 1);
            return ret;
        }
    }
    neg_hits++;
    std::vector<Edge_weight> ret(budget + 1, 0);
    for(auto [v, dag] : poss) {
        std::vector<Edge_weight> tmp;
        if(dag) {
            dags++;
            distance_to_goal_ = std::vector<Edge_length>(adjacency_.size(), invalid_);
            pruning_dijkstra(terminals_[1], start, distance_to_goal_, budget);
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
    }
    distance_to_goal_ = std::vector<Edge_length>(adjacency_.size(), invalid_);
    pruning_dijkstra(terminals_[1], start, distance_to_goal_, budget);
    prune_singleout(start);
    cache_[start][distance_to_goal_] = std::make_pair(budget, ret);
    visited_[start] = false;
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

void Search::prune_singleout(Vertex start) {
    std::vector<Vertex> to_check;
    for(Vertex v = 0; v < adjacency_.size(); v++) {
        if(distance_to_goal_[v] == invalid_ || v == terminals_[1]) {
            continue;
        }
        int neighbor = -1;
        for(auto neigh : neighbors(v)) {
            if(distance_to_goal_[neigh] != invalid_ || neigh == start) {
                if(neighbor == -1) {
                    neighbor = neigh;
                } else {
                    neighbor = -1;
                    break;
                }
            }
        }
        if(neighbor != -1) {
            distance_to_goal_[v] = invalid_;
            if(neighbor < v) {
                to_check.push_back(neighbor);
            }
        }
    }
    while(!to_check.empty()) {
        auto v = to_check.back();
        to_check.pop_back();
        if(distance_to_goal_[v] == invalid_ || v == terminals_[1]) {
            continue;
        }
        int neighbor = -1;
        for(auto neigh : neighbors(v)) {
            if(distance_to_goal_[neigh] != invalid_ || neigh == start) {
                if(neighbor == -1) {
                    neighbor = neigh;
                } else {
                    neighbor = -1;
                    break;
                }
            }
        }
        if(neighbor != -1) {
            distance_to_goal_[v] = invalid_;
            if(neighbor < v) {
                to_check.push_back(neighbor);
            }
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
