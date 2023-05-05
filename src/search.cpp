#include "search.h"
#include <queue>

Search::Search(Graph& input) :  max_length_(input.max_length_),
                                terminals_(input.terminals_),
                                neighbors_(input.neighbors_),
                                adjacency_(input.adjacency_),
                                invalid_(std::numeric_limits<Edge_length>::max() - max_length_ - 1),
                                distance_to_goal_(adjacency_.size(), invalid_)  {
    assert(terminals_.size() == 2);
    dijkstra(terminals_[1], distance_to_goal_);
}

std::vector<Edge_weight> Search::search(Vertex start, Edge_length budget) {
    if(start == terminals_[1]) {
        // we have one path of length zero.
        return {1};
    }
    Edge_length visited = distance_to_goal_[start];
    distance_to_goal_[start] = invalid_;
    std::vector<Edge_weight> ret(budget + 1, 0);
    for(auto v : neighbors(start)) {
        if(budget >= distance_to_goal_[v] + adjacency_[start][v].begin()->first) {
            std::vector<Edge_weight> tmp;
            if(budget == distance_to_goal_[v] + adjacency_[start][v].begin()->first) {
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
    distance_to_goal_[start] = visited;
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

void Search::dijkstra(Vertex start, std::vector<Edge_length>& distance) {
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
            if(cur_cost + 1 >= distance[w]) {
                continue;
            }
            Edge_length min_cost = std::numeric_limits<Edge_length>::max();
            auto weights = adjacency_[cur_vertex][w];
            assert(weights.size() > 0);
            for(auto &weight : weights) {
                min_cost = std::min(weight.first, min_cost);
            }
            if(cur_cost + min_cost < distance[w]) {
                distance[w] = min_cost + cur_cost;
                queue.push(std::make_pair(cur_cost + min_cost, w));
            }
        }
    }
}
