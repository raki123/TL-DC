#include "graph.h"
#include "search.h"
#include <clingo.hh>
#include <algorithm>
#include <map>
#include <limits>
#include <sstream>

Graph::Graph(std::istream &input) {
    char dec;
    input >> dec;
    Vertex nr_vertices;
    Vertex nr_edges;
    std::string line;
    while(!input.eof()) {
        switch (dec)
        {
        case 'c':
            std::getline(input, line);  
            break;
        case 'p':
        {
            std::string str;
            input >> str >> nr_vertices >> nr_edges;
            assert(str == "edge");
            adjacency_ = std::vector<std::vector<std::map<Edge_length,Edge_weight>>>(nr_vertices + 2, 
                                    std::vector<std::map<Edge_length,Edge_weight>>(nr_vertices + 2, 
                                                            std::map<Edge_length,Edge_weight>()));
            neighbors_ = std::vector<std::set<Vertex>>(nr_vertices + 2, std::set<Vertex>());
            max_length_ = nr_vertices + 2;
            break;
        }
        case 'e':
            Vertex v,w;
            input >> v >> w;
            add_edge(Edge(v - 1,w - 1), Weight(1,1));
            break;
        case 'l':
            input >> max_length_;
            extra_paths_ = std::vector<Edge_weight>(max_length_ + 1, 0);
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
    assert(max_length_ + 1 == extra_paths_.size());
    assert(adjacency_.size() > 0);
    if(terminals_.empty()) {
        terminals_ = {nr_vertices, Vertex(nr_vertices + 1)};
        max_length_ += 2;
        extra_paths_.push_back(0);
        extra_paths_.push_back(0);
        for(Vertex v = 0; v < nr_vertices; v++) {
            add_edge(Edge(v,terminals_[0]), Weight(1,1));
            add_edge(Edge(v,terminals_[1]), Weight(1,1));
        }
        all_pair_ = true;
    }
}

void Graph::preprocess() {
    std::cerr << terminals_[0] << " " << terminals_[1] << std::endl;
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
    std::cerr << preprocess_tiny_separator() << std::endl;
    std::cerr << "Removed isolated: " << isolated_removed << std::endl;
    std::cerr << "Removed forwarder: " << forwarder_removed << std::endl;
    std::cerr << "Removed position determined: " << position_determined_removed << std::endl;
    std::cerr << "Removed twin edges: " << twin_edges_removed << std::endl;
    std::cerr << "Removed unreachable: " << unreachable_removed << std::endl;
    std::cerr << "Removed unusable edge: " << unusable_edge_removed << std::endl;
}

void Graph::print_stats() {
    Vertex nr_edges = 0;
    for(Vertex v = 0; v < adjacency_.size(); v++) {
        nr_edges += neighbors(v).size();
    }
    nr_edges /= 2;
    std::vector<Edge_length> distance_to_goal(adjacency_.size(), std::numeric_limits<Edge_length>::max());
    dijkstra(terminals_[1], distance_to_goal, {});
    std::cerr << "#vertices " << adjacency_.size() << " #edges " << nr_edges;
    std::cerr << " max. length " << max_length_ << " min. length " << distance_to_goal[terminals_[0]] << std::endl;
}

void Graph::normalize() {
    Vertex unnamed = std::numeric_limits<Vertex>::max();
    std::vector<Vertex> new_name(adjacency_.size(), unnamed);
    Vertex cur_name = 0;
    for(Vertex v = 0; v < adjacency_.size(); v++) {
        if(adjacency_[v].empty() && v != terminals_[0] && v != terminals_[1]) {
            continue;
        }
        if(new_name[v] == unnamed) {
            new_name[v] = cur_name++;
        }
    }
    auto new_adjacency = std::vector<std::vector<std::map<Edge_length,Edge_weight>>>(cur_name, 
                            std::vector<std::map<Edge_length,Edge_weight>>(cur_name, 
                                                    std::map<Edge_length,Edge_weight>()));
    auto new_neighbors = std::vector<std::set<Vertex>>(cur_name, std::set<Vertex>());
    for(Vertex v = 0; v < adjacency_.size(); v++) {
        if(adjacency_[v].empty()) {
            continue;
        }
        for(auto neigh : neighbors(v)) {
            assert(new_name[neigh] != unnamed);
            new_neighbors[new_name[v]].insert(new_name[neigh]);
            new_adjacency[new_name[v]][new_name[neigh]] = adjacency_[v][neigh];
        }
    }
    adjacency_ = new_adjacency;
    neighbors_ = new_neighbors;
    terminals_[0] = new_name[terminals_[0]];
    terminals_[1] = new_name[terminals_[1]];
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
    neighbors_[edge.first].insert(edge.second);
    neighbors_[edge.second].insert(edge.first);
}

void Graph::remove_vertex(Vertex v) {
    assert(v >= 0 && v < adjacency_.size());
    assert(!adjacency_[v].empty());
    for(auto neighbor : neighbors(v)) {
        adjacency_[neighbor][v].clear();
        neighbors_[neighbor].erase(v);
    }
    adjacency_[v].clear();
    neighbors_[v].clear();
}

void Graph::remove_edge(Edge edge) {
    assert(edge.first >= 0 && edge.first < adjacency_.size());
    assert(edge.second >= 0 && edge.second < adjacency_.size());
    assert(!adjacency_[edge.first].empty());
    assert(!adjacency_[edge.second].empty());
    adjacency_[edge.first][edge.second].clear();
    adjacency_[edge.second][edge.first].clear();
    neighbors_[edge.first].erase(edge.second);
    neighbors_[edge.second].erase(edge.first);
}

std::set<Vertex> Graph::neighbors(Vertex v) {
    assert(v >= 0 && v < adjacency_.size());
    return neighbors_[v];
}

Graph::Graph(Vertex n) :    max_length_(-1),
                            neighbors_(n, std::set<Vertex>()),
                            adjacency_(n, std::vector<std::map<Edge_length, Edge_weight>>(n, std::map<Edge_length, Edge_weight>())) {

}

Graph Graph::subgraph(std::vector<Vertex> restrict_to) {
    Graph ret(restrict_to.size());
    ret.max_length_ = max_length_;
    std::vector<Vertex> new_name(adjacency_.size(), std::numeric_limits<Vertex>::max());
    Vertex cur_name = 0;
    for(Vertex v : restrict_to) {
        new_name[v] = cur_name++;
    }
    for(Vertex v : restrict_to) {
        for(Vertex neigh : neighbors(v)) {
            if(new_name[neigh] != std::numeric_limits<Vertex>::max()) {
                ret.neighbors_[new_name[v]].insert(new_name[neigh]);
                ret.adjacency_[new_name[v]][new_name[neigh]] = adjacency_[v][neigh];
            }
        }
    }
    ret.extra_paths_ = std::vector<Edge_weight>(max_length_ + 1, 0);
    return ret;
}

void Graph::dijkstra(Vertex start, std::vector<Edge_length>& distance, const std::set<Vertex>& forbidden) {
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
            if(forbidden.count(w) > 0) {
                continue;
            }
            Edge_length min_cost = adjacency_[cur_vertex][w].begin()->first;
            if(cur_cost + min_cost < distance[w]) {
                distance[w] = min_cost + cur_cost;
                queue.push(std::make_pair(cur_cost + min_cost, w));
            }
        }
    }
}

std::vector<Vertex> Graph::find_separator(size_t size) {
    // build the program
    std::stringstream prog_str;
    prog_str << "{sep(X) : v(X)}" << size << ".\n\
    {r(X)}:- v(X), not sep(X).\n\
    :- e(X,Y), r(X), not r(Y), not sep(Y).\n\
    :- e(Y,X), r(X), not r(Y), not sep(Y).\n\
    ok_nr(X) :- v(X), not sep(X), not r(X).\n";
    for(Vertex v = 0; v < adjacency_.size(); v++) {
        if(!adjacency_[v].empty()) {
            prog_str << "v(" << v << ").\n";
        }
    }
    prog_str << "r(" << terminals_[0] << ").\n";
    prog_str << "r(" << terminals_[1] << ").\n";
    prog_str << ":- sep(" << terminals_[0] << "), sep(" << terminals_[1] << ").\n";
    prog_str << ":- ";
    bool first = true;
    for(Vertex v = 0; v < adjacency_.size(); v++) {
        if(!adjacency_[v].empty()) {
            if(first) {
                first = false;
            } else {
                prog_str << ", ";
            }
            prog_str << "not r(" << v << ")";
        }
    }
    prog_str << ".\n";
    prog_str << ":- ";
    first = true;
    for(Vertex v = 0; v < adjacency_.size(); v++) {
        if(!adjacency_[v].empty()) {
            if(first) {
                first = false;
            } else {
                prog_str << ", ";
            }
            prog_str << "not ok_nr(" << v << ")";
        }
    }
    prog_str << ".\n";
    for(Vertex v = 0; v < adjacency_.size(); v++) {
        for(Vertex w : neighbors(v)) {
            prog_str << "e(" << v << "," << w << ").\n";
        }
    }
    prog_str << "#show sep/1.\n";
    // initialize clingo
    Clingo::Logger logger = [](Clingo::WarningCode, char const *) {
        // std::cerr << message << std::endl;
    };
    Clingo::Control ctl{{}, logger, 20};
    ctl.add("base", {}, prog_str.str().c_str());
    ctl.ground({{"base", {}}});
    auto handle = ctl.solve();
    std::vector<Vertex> ret;
    if(handle.get().is_satisfiable()) {
        for(auto symbol : handle.model().symbols()) {
            auto substr = symbol.to_string().substr(4,symbol.to_string().length() - 4 - 1);
            ret.push_back(std::stoi(substr));
        }
    }
    return ret;
}

Vertex Graph::preprocess_start_goal_edges() {
    assert(terminals_.size() > 0);
    Vertex found = 0;
    if(adjacency_[terminals_[0]].size() > 0) {
        for(auto &[length, weight] : adjacency_[terminals_[0]][terminals_[1]]) {
            if(length <= max_length_) {
                extra_paths_[length] += weight;
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
        auto cur_neighbors = neighbors(v);
        if(cur_neighbors.size() == 1) {
            Vertex neighbor = *cur_neighbors.begin();
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
        auto cur_neighbors = neighbors(v);
        if(cur_neighbors.size() == 2) {
            found++;
            // if the neighbors are w1 and w2, then we can use any combination of weighted edges (w1,v),(v,w2) as a single edge (w1,w2)
            Vertex w1 = *cur_neighbors.begin(), w2 = *(++cur_neighbors.begin());
            for(auto w1_weight : adjacency_[v][w1]) {
                for(auto w2_weight : adjacency_[v][w2]) {
                    add_edge(Edge(w1,w2), Weight(w1_weight.first + w2_weight.first, w1_weight.second * w2_weight.second));
                }
            }
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
        auto tmp = neighbors(v);
        std::vector<Vertex> neighs(tmp.begin(), tmp.end());
        assert(std::find(neighs.begin(), neighs.end(), terminals_[1 - i]) == neighs.end());
        std::map<Vertex, std::vector<Vertex>> eq;
        for(auto neigh : neighs) {
            eq[neigh] = {neigh};
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
                auto it = w1_neighs.find(w2);;
                if(it != w1_neighs.end()) {
                    w1_neighs.erase(it);
                }
                it = w2_neighs.find(w1);
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

Vertex Graph::preprocess_tiny_separator() {
    std::vector<Vertex> separator = find_separator(2);
    if(separator.size() == 0) {
        return 0;
    }
    assert(separator.size() == 2);
    std::vector<char> right(adjacency_.size(), false);
    std::vector<Vertex> queue;
    if(terminals_[0] != separator[0] && terminals_[0] != separator[1]) {
        right[terminals_[0]] = true;
        queue.push_back(terminals_[0]);   
    }
    if(terminals_[1] != separator[0] && terminals_[1] != separator[1]) {
        right[terminals_[1]] = true;
        queue.push_back(terminals_[1]);   
    }
    while(!queue.empty()) {
        Vertex cur = queue.back();
        queue.pop_back();
        for(auto neigh : neighbors(cur)) {
            if(right[neigh] || neigh == separator[0] || neigh == separator[1]) {
                continue;
            }
            right[neigh] = true;
            queue.push_back(neigh);
        }
    }
    std::vector<Vertex> left = separator;
    for(Vertex v = 0; v < adjacency_.size(); v++) {
        if(right[v] || adjacency_[v].empty() || v == separator[0] || v == separator[1]) {
            continue;
        }
        left.push_back(v);
    }
    Graph left_graph = subgraph(left);
    left_graph.terminals_ = {0, 1};
    left_graph.print_stats();
    left_graph.preprocess();
    left_graph.print_stats();
    Search search(left_graph);
    auto res = search.search();

    return left.size() - 2;
}