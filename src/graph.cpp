#include "graph.h"
#include <assert.h>
#include <queue>
#include <algorithm>


Graph::Graph(std::istream &input) {
    char dec;
    input >> dec;
    Vertex nr_vertices;
    Vertex nr_edges;
    while(!input.eof()) {
        switch (dec)
        {
        case 'c':
            break;
        case 'p':
        {
            std::string str;
            input >> str >> nr_vertices >> nr_edges;
            assert(str == "edge");
            adjacency_ = std::vector<std::unordered_map<Vertex, std::vector<Weight>>>(nr_vertices, std::unordered_map<Vertex, std::vector<Weight>>());
            break;
        }
        case 'e':
            Vertex v,w;
            input >> v >> w;
            add_edge(Edge(v - 1,w - 1), Weight(1,1));
            break;
        case 'l':
            input >> max_length_;
            break;
        case 't':
            terminals_.resize(2);
            input >> terminals_[0] >> terminals_[1];
            terminals_[0]--;
            terminals_[1]--;
            break;
        default:
            std::cerr << "Invalid character " << dec << " at beginning of line." << std::endl;
            break;
        }
        input >> dec;
    }
}

void Graph::preprocess() {
    bool found = true;
    Vertex isolated_removed = 0;
    Vertex forwarder_removed = 0;
    Vertex unreachable_removed = 0;
    Vertex unusable_edge_removed = 0;
    while(found) {
        found = false;
        Vertex cur_isolated_removed = preprocess_isolated();
        found |= cur_isolated_removed > 0;
        isolated_removed += cur_isolated_removed;
        Vertex cur_forwarder_removed = preprocess_forwarder();
        found |= cur_forwarder_removed > 0;
        forwarder_removed += cur_forwarder_removed;
        if(terminals_.size() > 0) {
            Vertex cur_unreachable_removed = preprocess_unreachable();
            found |= cur_unreachable_removed > 0;
            unreachable_removed += cur_unreachable_removed;
            if(!found) {
                Vertex cur_unusable_edge_removed = preprocess_unusable_edge();
                found |= cur_unusable_edge_removed > 0;
                unusable_edge_removed += cur_unusable_edge_removed;
            }
        }
    }
    std::cerr << "Removed isolated: " << isolated_removed << std::endl;
    std::cerr << "Removed forwarder: " << forwarder_removed << std::endl;
    std::cerr << "Removed unreachable: " << unreachable_removed << std::endl;
    std::cerr << "Removed unusable edge: " << unusable_edge_removed << std::endl;
}

void Graph::add_edge(Edge edge, Weight weight) {
    assert(edge.first >= 0 && edge.first < adjacency_.size());
    assert(edge.second >= 0 && edge.second < adjacency_.size());
    adjacency_[edge.first][edge.second].push_back(weight);
    adjacency_[edge.second][edge.first].push_back(weight);
}

void Graph::remove_vertex(Vertex v) {
    assert(v >= 0 && v < adjacency_.size());
    for(auto neighbor : neighbors(v)) {
        adjacency_[neighbor].erase(v);
        adjacency_[v].clear();
    }
}

std::vector<Vertex> Graph::neighbors(Vertex v) {
    assert(v >= 0 && v < adjacency_.size());
    std::vector<Vertex> ret;
    for(auto & [w, weight] : adjacency_[v]) {
        ret.push_back(w);
    }
    return ret;
}

void Graph::remove_edge(Edge edge) {
    assert(edge.first >= 0 && edge.first < adjacency_.size());
    assert(edge.second >= 0 && edge.second < adjacency_.size());
    adjacency_[edge.first].erase(edge.second);
    adjacency_[edge.second].erase(edge.first);
}

void Graph::dijkstra(Vertex start, std::vector<Edge_length>& distance, const std::set<Vertex>& forbidden) {
    std::priority_queue<std::pair<Edge_length, Vertex>> queue;
    queue.push(std::make_pair(0, start));
    distance[start] = 0;
    while(!queue.empty()) {
        auto [cur_cost, cur_vertex] = queue.top();
        queue.pop();
        if(cur_cost > distance[cur_vertex]) {
            continue;
        }
        for(auto &[w, weights] : adjacency_[cur_vertex]) {
            if(cur_cost + 1 >= distance[w]) {
                continue;
            }
            if(forbidden.count(w) > 0) {
                continue;
            }
            Edge_length min_cost = std::numeric_limits<Edge_length>::max();
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

Vertex Graph::preprocess_isolated() {
    Vertex found = 0;
    for(size_t v = 0; v < adjacency_.size(); v++) {
        std::vector<Vertex> neighbors_ = neighbors(v);
        if(neighbors_.size() == 1) {
            Vertex neighbor = neighbors_[0];
            if(terminals_.size() > 0 
                && (terminals_[0] == v || terminals_[1] == v)
                && (terminals_[0] == neighbor || terminals_[1] == neighbor)) {
                    // if the graph consists only of v -- neighbor, we cannot preprocess this node away
                continue;
            }
            found++;
            if(terminals_.size() > 0 && (terminals_[0] == v || terminals_[1] == v)) {
                // we will remove a terminal -> use a different terminal instead
                if(terminals_[0] == v) {
                    terminals_[0] = neighbor;
                } else {
                    terminals_[1] = neighbor; 
                }
                // in this case we use the edge(s) of the removed vertex and need to adapt the edges of the neighbor accordingly.
                for(auto &[w, old_weights] : adjacency_[neighbor]) {
                    if(w == v) {
                        continue;
                    }
                    std::unordered_map<Edge_length, Edge_weight> new_weights;
                    for(auto old_weight : old_weights) {
                        for(auto weight : adjacency_[v][neighbor]) {
                            new_weights[old_weight.first + weight.first] += old_weight.second * weight.second;
                        }
                    }
                    old_weights = std::vector<Weight>(new_weights.begin(), new_weights.end());
                }
            }
            remove_vertex(v);
        }
    }
    return found;
}

Vertex Graph::preprocess_forwarder() {
    Vertex found = 0;
    for(size_t v = 0; v < adjacency_.size(); v++) {
        if(terminals_.size() > 0 && (terminals_[0] == v || terminals_[1] == v)) {
            // we cannot remove terminals this way
            continue;
        }
        std::vector<Vertex> neighbors_ = neighbors(v);
        if(neighbors_.size() == 2) {
            found++;
            // if the neighbors are w1 and w2, then we can use any combination of weighted edges (w1,v),(v,w2) as a single edge (w1,w2)
            Vertex w1 = neighbors_[0], w2 = neighbors_[1];
            std::unordered_map<Edge_length, Edge_weight> new_weights;
            for(auto old_weight : adjacency_[w1][w2]) {
                new_weights[old_weight.first] += old_weight.second;
            }
            for(auto w1_weight : adjacency_[v][w1]) {
                for(auto w2_weight : adjacency_[v][w2]) {
                    new_weights[w1_weight.first + w2_weight.first] += w1_weight.second * w2_weight.second;
                }
            }
            adjacency_[w1][w2] = std::vector<Weight>(new_weights.begin(), new_weights.end());
            adjacency_[w2][w1] = std::vector<Weight>(new_weights.begin(), new_weights.end());
            remove_vertex(v);
        }
    }
    return found;
}

Vertex Graph::preprocess_unreachable() {
    assert(terminals_.size() > 0);
    std::vector<Edge_length> distance_from_start(adjacency_.size(), std::numeric_limits<Edge_length>::max());
    dijkstra(terminals_[0], distance_from_start, {});

    std::vector<Edge_length> distance_to_goal(adjacency_.size(), std::numeric_limits<Edge_length>::max());
    dijkstra(terminals_[1], distance_to_goal, {});

    Vertex found = 0;
    for(Vertex v = 0; v < adjacency_.size(); v++) {
        if(adjacency_[v].size() == 0) {
            continue;
        }
        if(distance_from_start[v] + distance_to_goal[v] > max_length_) {
            found++;
            remove_vertex(v);
        }
    }
    return found;
}

Vertex Graph::preprocess_unusable_edge() {
    assert(terminals_.size() > 0);
    std::vector<std::vector<Edge_length>> distances_from_start(adjacency_.size(), std::vector<Edge_length>(adjacency_.size(), std::numeric_limits<Edge_length>::max()));
    std::vector<std::vector<Edge_length>> distances_to_goal(adjacency_.size(), std::vector<Edge_length>(adjacency_.size(), std::numeric_limits<Edge_length>::max()));
    Vertex found = 0;
    for(Vertex v = 0; v < adjacency_.size(); v++) {
        if(adjacency_[v].size() == 0) {
            continue;
        }
        dijkstra(terminals_[0], distances_from_start[v], {v});
        dijkstra(terminals_[1], distances_to_goal[v], {v});
    }
    for(Vertex v = 0; v < adjacency_.size(); v++) {
        std::vector<Vertex> remove_completely;
        for(auto &[w, weights] : adjacency_[v]) {
            Edge_length min_without_edge = std::min(distances_from_start[v][w] + distances_to_goal[w][v], distances_from_start[w][v] + distances_to_goal[v][w]);
            std::vector<Weight> new_weights;
            for(auto &weight : weights) {
                if(weight.first + min_without_edge <= max_length_) {
                    new_weights.push_back(weight);
                } else {
                    found++;
                }
            }
            if(new_weights.size() == 0) {
                remove_completely.push_back(w);
            } else {
                weights = new_weights;
            }
        }
        for(auto w : remove_completely) {
            remove_edge(Edge(v,w));
        }
    }
    return found;
}