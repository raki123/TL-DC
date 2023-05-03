#include "graph.h"
#include <queue>
#include <algorithm>
#include <map>

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
            adjacency_ = std::vector<std::vector<std::map<Edge_length,Edge_weight>>>(nr_vertices, 
                                    std::vector<std::map<Edge_length,Edge_weight>>(nr_vertices, 
                                                            std::map<Edge_length,Edge_weight>()));
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
    Vertex position_determined_removed = 0;
    Vertex twin_edges_removed = 0;
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
        Vertex cur_twin_edges_removed = preprocess_twins();
        found |= cur_twin_edges_removed > 0;
        twin_edges_removed += cur_twin_edges_removed;
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
        if(!found) {
            Vertex cur_position_determined_removed = preprocess_position_determined();
            found |= cur_position_determined_removed > 0;
            position_determined_removed += cur_position_determined_removed;
        }
    }
    std::cerr << "Removed isolated: " << isolated_removed << std::endl;
    std::cerr << "Removed forwarder: " << forwarder_removed << std::endl;
    std::cerr << "Removed position determined: " << position_determined_removed << std::endl;
    std::cerr << "Removed twin edges: " << twin_edges_removed << std::endl;
    std::cerr << "Removed unreachable: " << unreachable_removed << std::endl;
    std::cerr << "Removed unusable edge: " << unusable_edge_removed << std::endl;
}

void Graph::encode_unary(std::ostream& output) {
    output << "reach(X, 0) :- start(X).\n";
    output << ":- goal(X)";
    for(int i = 0; i <= max_length_; i++) {
        output << ", not reach(X, " << i << ")";
    }
    output << ".\n";
    output << "start(" << terminals_[0] << ").\n";
    output << "goal(" << terminals_[1] << ").\n";
    output << ":- reach(X, L), reach(X, L'), L != L'.\n";
    std::vector<Edge_length> distance_from_start(adjacency_.size(), std::numeric_limits<Edge_length>::max());
    dijkstra(terminals_[0], distance_from_start, {});

    std::vector<Edge_length> distance_to_goal(adjacency_.size(), std::numeric_limits<Edge_length>::max());
    dijkstra(terminals_[1], distance_to_goal, {});

    for(Vertex v = 0; v < adjacency_.size(); v++) {
        if(v == terminals_[1] || adjacency_[v].empty()) {
            continue;
        }
        assert(distance_from_start[v] + distance_to_goal[v] <= max_length_);
        std::vector<std::pair<Vertex,Edge_length>> edges;
        for(auto &w : neighbors(v)) {
            if(w == terminals_[0]) {
                continue;
            }
            for(auto &weight : adjacency_[v][w]) {
                assert(std::find(edges.begin(), edges.end(), std::make_pair(w, weight.first)) == edges.end());
                edges.push_back(std::make_pair(w, weight.first));
                // if(weight.second != 1) {
                //     std::cerr << "Found edge with weight " << weight.second << std::endl;
                // }
            }
        }
        for(auto &[w, length] : edges) {
            assert(distance_to_goal[w] <= max_length_);
            output << "reach(" << w << ",L + L') :- reach(" << v << ", L), edge(" << v << "," << w << ",L'), L <= " << max_length_ - distance_to_goal[w] - length;
            output << ", L>= " << distance_from_start[v] << ".\n";
            output << "edge(" << v << "," << w << "," << length << ") :- reach(" << v << ", L), L <= " << max_length_ - distance_to_goal[w] - length;
            output << ", L>= " << distance_from_start[v];
            for(auto &[wp, lengthp] : edges) {
                if(w != wp || length != lengthp) {
                    output << ", not edge(" << v << "," << wp << "," << lengthp << ")";
                }
            }
            output << ".\n";
        }
    }
}

void Graph::add_edge(Edge edge, Weight weight) {
    assert(edge.first != edge.second);
    assert(edge.first >= 0 && edge.first < adjacency_.size());
    assert(edge.second >= 0 && edge.second < adjacency_.size());
    assert(!adjacency_[edge.first].empty());
    assert(!adjacency_[edge.second].empty());
    if(weight.first > max_length_) {
        return;
    }
    adjacency_[edge.first][edge.second][weight.first] += weight.second;
    adjacency_[edge.second][edge.first][weight.first] += weight.second;
}

void Graph::remove_vertex(Vertex v) {
    assert(v >= 0 && v < adjacency_.size());
    assert(!adjacency_[v].empty());
    for(auto neighbor : neighbors(v)) {
        adjacency_[neighbor][v].clear();
    }
    adjacency_[v].clear();
}

std::vector<Vertex> Graph::neighbors(Vertex v) {
    assert(v >= 0 && v < adjacency_.size());
    std::vector<Vertex> ret;
    if(adjacency_[v].empty()) {
        return ret;
    }
    for(Vertex w = 0; w < adjacency_.size(); w++) {
        if(adjacency_[v][w].size() > 0) {
            ret.push_back(w);
        }
    }
    return ret;
}

void Graph::remove_edge(Edge edge) {
    assert(edge.first >= 0 && edge.first < adjacency_.size());
    assert(edge.second >= 0 && edge.second < adjacency_.size());
    assert(!adjacency_[edge.first].empty());
    assert(!adjacency_[edge.second].empty());
    adjacency_[edge.first][edge.second].clear();
    adjacency_[edge.second][edge.first].clear();
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
        for(auto &w : neighbors(cur_vertex)) {
            if(cur_cost + 1 >= distance[w]) {
                continue;
            }
            if(forbidden.count(w) > 0) {
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

Vertex Graph::preprocess_start_goal_edges() {
    assert(terminals_.size() > 0);
    Vertex found = 0;
    if(adjacency_[terminals_[0]].size() > 0) {
        for(auto &[length, weight] : adjacency_[terminals_[0]][terminals_[1]]) {
            if(length <= max_length_) {
                extra_paths_ += weight;
            }
        }
        found = adjacency_[terminals_[0]][terminals_[1]].size();
        remove_edge(Edge(terminals_[0], terminals_[1]));
    }
    return found;
}

Vertex Graph::preprocess_isolated() {
    Vertex found = 0;
    for(size_t v = 0; v < adjacency_.size(); v++) {
        std::vector<Vertex> cur_neighbors = neighbors(v);
        if(cur_neighbors.size() == 1) {
            Vertex neighbor = cur_neighbors[0];
            found++;
            if(terminals_.size() > 0 && (terminals_[0] == v || terminals_[1] == v)) {
                if(terminals_.size() > 0 && (terminals_[0] == neighbor || terminals_[1] == neighbor)) {
                    // if the graph consists only of v -- neighbor, we have a solution
                    return preprocess_start_goal_edges() + found;
                }
                // we will remove a terminal -> use a different terminal instead
                if(terminals_[0] == v) {
                    terminals_[0] = neighbor;
                } else {
                    terminals_[1] = neighbor; 
                }
                // in this case we use the edge(s) of the removed vertex and need to adapt the edges of the neighbor accordingly.
                for(auto w : neighbors(neighbor)) {
                    if(w == v) {
                        continue;
                    }
                    std::map<Edge_length, Edge_weight> new_weights;
                    for(auto old_weight : adjacency_[neighbor][w]) {
                        for(auto weight : adjacency_[v][neighbor]) {
                            new_weights[old_weight.first + weight.first] += old_weight.second * weight.second;
                        }
                    }
                    adjacency_[neighbor][w] = new_weights;
                    adjacency_[w][neighbor] = new_weights;
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
        std::vector<Vertex> cur_neighbors = neighbors(v);
        if(cur_neighbors.size() == 2) {
            found++;
            // if the neighbors are w1 and w2, then we can use any combination of weighted edges (w1,v),(v,w2) as a single edge (w1,w2)
            Vertex w1 = cur_neighbors[0], w2 = cur_neighbors[1];
            std::map<Edge_length, Edge_weight> new_weights;
            for(auto old_weight : adjacency_[w1][w2]) {
                new_weights[old_weight.first] += old_weight.second;
            }
            for(auto w1_weight : adjacency_[v][w1]) {
                for(auto w2_weight : adjacency_[v][w2]) {
                    new_weights[w1_weight.first + w2_weight.first] += w1_weight.second * w2_weight.second;
                }
            }
            adjacency_[w1][w2] = new_weights;
            adjacency_[w2][w1] = new_weights;
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

Vertex Graph::preprocess_twins() {
    Vertex found = 0;
    assert(terminals_.size() > 0);
    for(int i = 0; i < 2; i++) {
        Vertex v = terminals_[i];
        preprocess_start_goal_edges();
        auto neighs = neighbors(v);
        assert(std::find(neighs.begin(), neighs.end(), terminals_[1 - i]) == neighs.end());
        std::map<Vertex, std::vector<Vertex>> eq;
        for(size_t j = 0; j < neighs.size(); j++) {
            eq[neighs[j]] = {neighs[j]};
        }
        for(size_t j = 0; j < neighs.size(); j++) {
            for(size_t k = j+1; k < neighs.size(); k++) {
                Vertex w1 = neighs[j];
                Vertex w2 = neighs[k];
                // compare the neighborhoods of w1 and w2 to see if they are equal (up to w1/w2)
                bool problem = false;
                auto w1_neighs = neighbors(w1);
                auto w2_neighs = neighbors(w2);
                if(w1_neighs.size() != w2_neighs.size()) {
                    continue;
                }
                auto it = std::find(w1_neighs.begin(), w1_neighs.end(), w2);
                if(it != w1_neighs.end()) {
                    w1_neighs.erase(it);
                }
                it = std::find(w2_neighs.begin(), w2_neighs.end(), w1);
                if(it != w2_neighs.end()) {
                    w2_neighs.erase(it);
                }
                if(w1_neighs != w2_neighs) {
                    continue;
                }
                // we have the same neighbors. 
                // make sure that also the weights are the same
                for(auto w : w1_neighs) {
                    if(adjacency_[w2][w].size() != adjacency_[w1][w].size()) {
                        problem = true;
                        break;
                    }
                    if(adjacency_[w2][w] != adjacency_[w1][w]) {
                        problem = true;
                        break;
                    }
                }
                if(!problem) {
                    eq[neighs[j]].push_back(neighs[k]);
                    std::swap(neighs[k],neighs.back());
                    neighs.pop_back();
                    k--;
                }
            }
        }
        for(auto &[representative, others] : eq) {
            if(others.size() == 1) {
                continue;
            }
            // multiply the weight of the edges between terminal and representative by the number of twins
            Edge_weight factor = others.size();
            for(auto &[length, weight] : adjacency_[v][representative]) {
                weight *= factor;
            }
            for(auto &[length, weight] : adjacency_[representative][v]) {
                weight *= factor;
            }
            // remove the edges between the terminal and the other twins
            for(auto other : others) {
                if(other != representative) {
                    found += adjacency_[v][other].size();
                    remove_edge(Edge(v, other));
                }
            }
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
        for(auto &w : neighbors(v)) {
            Edge_length min_without_edge = std::min(distances_from_start[v][w] + distances_to_goal[w][v], distances_from_start[w][v] + distances_to_goal[v][w]);
            std::map<Edge_length, Edge_weight> new_weights;
            for(auto &weight : adjacency_[v][w]) {
                if(weight.first + min_without_edge <= max_length_) {
                    new_weights.insert(weight);
                } else {
                    found++;
                }
            }
            if(new_weights.size() == 0) {
                remove_completely.push_back(w);
            } else {
                adjacency_[v][w] = new_weights;
            }
        }
        for(auto w : remove_completely) {
            remove_edge(Edge(v,w));
        }
    }
    return found;
}

Vertex Graph::preprocess_position_determined() {
    assert(terminals_.size() > 0);
    std::vector<Edge_length> distance_from_start(adjacency_.size(), std::numeric_limits<Edge_length>::max());
    dijkstra(terminals_[0], distance_from_start, {});

    std::vector<Edge_length> distance_to_goal(adjacency_.size(), std::numeric_limits<Edge_length>::max());
    dijkstra(terminals_[1], distance_to_goal, {});

    Vertex found = 0;
    for(Vertex v = 0; v < adjacency_.size(); v++) {
        if(adjacency_[v].empty() || neighbors(v).empty()) {
            continue;
        }
        if(terminals_[0] == v || terminals_[1] == v) {
            continue;
        }
        if(distance_from_start[v] + distance_to_goal[v] == max_length_) {
            found++;
            std::vector<std::tuple<Edge, Edge_length, Edge_weight>> to_add;
            for(auto &w : neighbors(v)) {
                for(auto &wp : neighbors(v)) {
                    if(w >= wp) {
                        continue;
                    }
                    for(auto &weight : adjacency_[v][w]) {
                        for(auto &weightp : adjacency_[v][wp]) {
                            to_add.push_back(std::make_tuple(Edge(w,wp), weight.first + weightp.first, weight.second * weightp.second));
                        }
                    }
                }
            }
            remove_vertex(v);
            for(auto &[edge, length, weight] : to_add) {
                add_edge(edge, Weight(length, weight));
            }
        }
    }
    return found;
}