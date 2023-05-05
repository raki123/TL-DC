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
    std::vector<Edge_weight> ret(budget + 1, 0);
    for(auto v : neighbors(start)) {
        if(budget >= distance_to_goal_[v] + (adjacency_[start][v].begin())->first) {
            Edge_length visited = distance_to_goal_[v];
            distance_to_goal_[v] = invalid_;
            auto tmp = search(v, budget - (adjacency_[start][v].begin())->first);
            for(size_t i = 0; i < tmp.size(); i++) {
                for(auto &[length, weight] : adjacency_[start][v]) {
                    if(length + i > budget) {
                        break;
                    }
                    ret[length + i] += weight*tmp[i];
                }
            }
            distance_to_goal_[v] = visited;
        }
    }
    return ret;
}

void Search::dijkstra(Vertex start, std::vector<Edge_length>& distance) {
    std::priority_queue<std::pair<Edge_length, Vertex>> queue;
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
