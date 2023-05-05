#include "search.h"
#include <queue>

Search::Search(Graph& input) :  max_length_(input.max_length_),
                                terminals_(input.terminals_),
                                neighbors_(input.neighbors_),
                                adjacency_(input.adjacency_),
                                invalid_(std::numeric_limits<Edge_length>::max() - max_length_ - 1),
                                distance_to_goal_(adjacency_.size(), invalid_),
                                visited_(adjacency_.size(), false),
                                cache_(adjacency_.size(), std::map<CacheKey, std::pair<Edge_length, std::vector<Edge_weight>>>())  {
    assert(terminals_.size() == 2);
    dijkstra(terminals_[1], distance_to_goal_, max_length_);
}

std::vector<Edge_weight> Search::search(Vertex start, Edge_length budget) {
    edges++;
    if(start == terminals_[1]) {
        // we have one path of length zero.
        return {1};
    }
    auto cached_result = cache_[start].find(distance_to_goal_);
    if(cached_result != cache_[start].end()) {
        if(cached_result->second.first >= budget) {
            pos_hits++;
            std::vector<Edge_weight> ret(cached_result->second.second.begin(), cached_result->second.second.begin() + budget + 1);
            return ret;
        }
    }
    neg_hits++;
    visited_[start] = true;
    // Edge_length visited = distance_to_goal_[start];
    // distance_to_goal_[start] = invalid_;
    distance_to_goal_ = std::vector<Edge_length>(adjacency_.size(), invalid_);
    dijkstra(terminals_[1], distance_to_goal_, budget);
    std::vector<Edge_weight> ret(budget + 1, 0);
    for(auto v : neighbors(start)) {
        if(budget >= distance_to_goal_[v] + adjacency_[start][v].begin()->first) {
            std::vector<Edge_weight> tmp;
            if(budget == distance_to_goal_[v] + adjacency_[start][v].begin()->first) {
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
        }
    }
    // distance_to_goal_[start] = visited;
    visited_[start] = false;
    distance_to_goal_ = std::vector<Edge_length>(adjacency_.size(), invalid_);
    dijkstra(terminals_[1], distance_to_goal_, budget);
    cache_[start][distance_to_goal_] = std::make_pair(budget, ret);
    return ret;
}

std::vector<Edge_weight> Search::dag_search(Vertex start, Edge_length budget) {
    std::vector<Edge_length> actual_distance(adjacency_.size(), invalid_);
    DijkstraQueue queue;
    queue.push(std::make_pair(0, terminals_[1]));
    actual_distance[terminals_[1]] = 0;
    while(!queue.empty()) {
        auto [cur_cost, cur_vertex] = queue.top();
        queue.pop();
        if(cur_cost > actual_distance[cur_vertex]) {
            continue;
        }
        for(auto &w : neighbors(cur_vertex)) {
            if(cur_cost + 1 >= actual_distance[w] || distance_to_goal_[w] == invalid_) {
                continue;
            }
            Edge_length min_cost = std::numeric_limits<Edge_length>::max();
            auto weights = adjacency_[cur_vertex][w];
            assert(weights.size() > 0);
            for(auto &weight : weights) {
                min_cost = std::min(weight.first, min_cost);
            }
            if(cur_cost + min_cost < actual_distance[w] && cur_cost + min_cost <= budget) {
                actual_distance[w] = min_cost + cur_cost;
                queue.push(std::make_pair(cur_cost + min_cost, w));
            }
        }
    }
    if(actual_distance[start] > budget) {
        fails++;
        return {};
    }
    assert(actual_distance[start] == budget);
    std::priority_queue<std::pair<Edge_length, Vertex>> biggest_first_queue;
    std::vector<Edge_weight> dp(adjacency_.size(), 0);
    std::vector<char> in_queue(adjacency_.size(), false);
    dp[start] = 1;
    biggest_first_queue.push(std::make_pair(actual_distance[start], start));
    while(!biggest_first_queue.empty()) {
        auto [cur_cost, cur_vertex] = biggest_first_queue.top();
        biggest_first_queue.pop();
        for(auto &w : neighbors(cur_vertex)) {
            if(actual_distance[w] == invalid_) {
                continue;
            }
            auto edge_cost = adjacency_[w][cur_vertex].begin()->first;
            // this edge can be part of a shorted path in direction w -> cur_vertex
            if(actual_distance[w] == actual_distance[cur_vertex] + edge_cost) {
                auto factor = adjacency_[w][cur_vertex].begin()->second;
                dp[cur_vertex] += factor*dp[w];
            // this edge can be part of a shorted path in direction cur_vertex -> w
            } else if(actual_distance[w] + edge_cost == actual_distance[cur_vertex]) {
                if(!in_queue[w]) {
                    biggest_first_queue.push(std::make_pair(actual_distance[w], w));
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
            Edge_length min_cost = std::numeric_limits<Edge_length>::max();
            auto weights = adjacency_[cur_vertex][w];
            assert(weights.size() > 0);
            for(auto &weight : weights) {
                min_cost = std::min(weight.first, min_cost);
            }
            if(cur_cost + min_cost < distance[w] && cur_cost + min_cost <= budget) {
                distance[w] = min_cost + cur_cost;
                queue.push(std::make_pair(cur_cost + min_cost, w));
            }
        }
    }
}

void Search::print_stats() {
    std::cerr << "Cache hits: " << pos_hits << " out of " << pos_hits + neg_hits << " tries." << std::endl;
    std::cerr << "Fails: " << fails << " out of " << dags << " DAG searches." << std::endl;
    std::cerr << "Edges: " << edges << std::endl;
}
