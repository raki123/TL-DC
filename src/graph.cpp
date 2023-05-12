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
            exclude_ = std::vector<std::vector<Vertex>>(nr_vertices + 2, std::vector<Vertex>());
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
    bool found = true;
    Vertex isolated_removed = 0;
    Vertex forwarder_removed = 0;
    Vertex position_determined_removed = 0;
    Vertex twin_edges_removed = 0;
    Vertex unreachable_removed = 0;
    Vertex unusable_edge_removed = 0;
    Vertex two_sep_removed = 0;
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
        if(!found) {
            Vertex cur_two_sep_removed = preprocess_two_separator();
            found |= cur_two_sep_removed > 0;
            two_sep_removed += cur_two_sep_removed;
        }
    }
    std::cerr << "Removed isolated: " << isolated_removed << std::endl;
    std::cerr << "Removed forwarder: " << forwarder_removed << std::endl;
    std::cerr << "Removed position determined: " << position_determined_removed << std::endl;
    std::cerr << "Removed twin edges: " << twin_edges_removed << std::endl;
    std::cerr << "Removed unreachable: " << unreachable_removed << std::endl;
    std::cerr << "Removed unusable edge: " << unusable_edge_removed << std::endl;
    std::cerr << "Removed due to 2-separation: " << two_sep_removed << std::endl;
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
    auto new_exclude = std::vector<std::vector<Vertex>>(cur_name, std::vector<Vertex>());
    for(Vertex v = 0; v < adjacency_.size(); v++) {
        if(adjacency_[v].empty()) {
            continue;
        }
        for(auto neigh : neighbors(v)) {
            assert(new_name[neigh] != unnamed);
            new_neighbors[new_name[v]].insert(new_name[neigh]);
            new_adjacency[new_name[v]][new_name[neigh]] = adjacency_[v][neigh];
        }
        for(auto excluded : exclude_[v]) {
            assert(new_name[excluded] != unnamed);
            new_exclude[new_name[v]].push_back(new_name[excluded]);
        }
    }
    adjacency_ = new_adjacency;
    neighbors_ = new_neighbors;
    exclude_ = new_exclude;
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

void Graph::add_exclude(Vertex v, Vertex w) {
    assert(v != w);
    assert(v >= 0 && v < adjacency_.size());
    assert(w >= 0 && w < adjacency_.size());
    assert(!adjacency_[v].empty());
    assert(!adjacency_[w].empty());
    exclude_[v].push_back(w);
    exclude_[w].push_back(v);
}

void Graph::remove_vertex(Vertex v) {
    assert(v >= 0 && v < adjacency_.size());
    assert(!adjacency_[v].empty());
    assert(exclude_[v].empty());
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
                            adjacency_(n, std::vector<std::map<Edge_length, Edge_weight>>(n, std::map<Edge_length, Edge_weight>())),
                            exclude_(n, std::vector<Vertex>()) {

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
        for(Vertex excluded : exclude_[v]) {
            if(new_name[excluded] != std::numeric_limits<Vertex>::max()) {
                // FIXME, why does this happen??
                assert(new_name[excluded] != std::numeric_limits<Vertex>::max());
                ret.exclude_[new_name[v]].push_back(new_name[excluded]);
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

std::vector<std::vector<Vertex>> Graph::components(const std::set<Vertex>& forbidden) {
    std::vector<std::vector<Vertex>> ret;
    std::vector<char> visited(adjacency_.size(), false);
    std::vector<Vertex> queue;
    for(Vertex v = 0; v < adjacency_.size(); v++) {
        if(adjacency_[v].empty() || visited[v] || forbidden.count(v) > 0) {
            continue;
        }
        std::vector<Vertex> component = {v};
        visited[v] = true;
        queue.push_back(v);
        while(!queue.empty()) {
            Vertex cur = queue.back();
            queue.pop_back();
            for(auto w : neighbors(cur)) {
                if(visited[w] || forbidden.count(w) > 0) {
                    continue;
                }
                visited[w] = true;
                queue.push_back(w);
                component.push_back(w);
            }
        }
        ret.push_back(component);
    }
    return ret;
}

std::vector<Vertex> Graph::find_separator(size_t size, size_t min_component_size, bool terminals_in_same) {
    // build the program
    std::stringstream prog_str;
    prog_str << size << "{sep(X) : v(X), not fixed(X)}" << size << ".\n\
    {r(X)}:- v(X), not sep(X).\n\
    :- e(X,Y), r(X), not r(Y), not sep(Y).\n\
    ok_nr(X) :- v(X), not sep(X), not r(X).\n\
    :- #count{ X : r(X)} < " << min_component_size << ".\n\
    :- #count{ X : ok_nr(X)} < " << min_component_size << ".\n";
    for(Vertex v = 0; v < adjacency_.size(); v++) {
        if(!adjacency_[v].empty()) {
            prog_str << "v(" << v << ").\n";
        }
        if(fixed(v)) {
            // dont use fixed vertices as separators
            prog_str << "fixed(" << v << ").\n";
            // and have fixed vertices in the same component as the others that are in an exclusion constraint with them
            for(auto excluded : exclude_[v]) {
                prog_str << "e(" << v << "," << excluded << ").\n";
            }
        }
    }
    if(terminals_in_same) {
        prog_str << ":- r(" << terminals_[0] << "), not sep(" << terminals_[0] << ").\n";
        prog_str << ":- r(" << terminals_[1] << "), not sep(" << terminals_[1] << ").\n";
    }
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
            // std::cout << substr << " " << symbol << std::endl;
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
        if(fixed(v)) {
            continue;
        }
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
        if(fixed(v)) {
            // we cannot remove fixed vertices this way
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
        if(fixed(v)) {
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
        // first make sure that none of the neighbors are fixed because then we cannot touch them here
        bool any_fixed = false;
        for(auto neigh : neighs) {
            if(fixed(neigh)) {
                any_fixed = true;
                break;
            }
        }
        if(any_fixed) {
            continue;
        }
        // gather equivalent vertices (i.e. twins) and merge them into one vertex each
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
        // merge the equivalent vertices
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
        if(fixed(v)) {
            continue;
        }
        std::vector<Vertex> remove_completely;
        for(auto &w : neighbors(v)) {
            if(fixed(w)) {
                continue;
            }
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
        if(fixed(v)) {
            // we cannot touch fixed vertices here
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

Vertex Graph::preprocess_two_separator() {
    Vertex found = 0;
    // if the separator can split start and goal it should have a minimum size
    // otherwise we can separate if the start or the goal have less than two neighbors
    std::vector<Vertex> separator = find_separator(2, 6, false);
    if(separator.size() == 0) {
        separator = find_separator(2, 1, true);
    }
    if(separator.size() == 0) {
        return 0;
    }
    assert(separator.size() == 2);

    std::vector<Edge_length> distance_from_start(adjacency_.size(), std::numeric_limits<Edge_length>::max());
    dijkstra(terminals_[0], distance_from_start, {});

    std::vector<Edge_length> distance_to_goal(adjacency_.size(), std::numeric_limits<Edge_length>::max());
    dijkstra(terminals_[1], distance_to_goal, {});

    // compute the components induced by the separator
    std::set<Vertex> forbidden(separator.begin(), separator.end());
    std::vector<std::vector<Vertex>> comps = components(forbidden);
    assert(comps.size() > 1);
    // try to reduce each of the induced components
    for(auto &comp : comps) {
        // different cases depending on whether there 0/1/2 of the terminals in the component
        bool found_start = std::find(comp.begin(), comp.end(), terminals_[0]) != comp.end();
        bool found_goal = std::find(comp.begin(), comp.end(), terminals_[1]) != comp.end();
        if(found_goal && found_start) {
            // nothing we can do
            continue;
        } 
        if(found_goal || found_start) {
            // let t be the terminal found and s_1, s_2 the separators
            // solve four subqueries:
            // C(1,Y) = number of paths from t to s_1 that may use s_2
            // C(1,N) = number of paths from t to s_1 that may not use s_2
            // C(2,Y) = number of paths from t to s_2 that may use s_1
            // C(2,N) = number of paths from t to s_2 that may not use s_1
            // then we modify the graph by replacing the component with 
            //      x---t---y
            //      |  / \  |
            //      | /   \ |
            //      |/     \|
            //     s_1(---)s_2
            // where:
            // {x,t}, {y,t} have (0,1)
            // {t,s_1} has C(1,N)
            // {t,s_2} has C(2,N)
            // {x,s_1} has C(1,Y) - C(1,N) (i.e. all paths from t to s_1 that use s_2)
            // {y,s_2} has C(2,Y) - C(2,N) (i.e. all paths from t to s_2 that use s_1)
            // and 
            // x and s_2 exclude each other
            // y and s_1 exclude each other

            // minus three because we also add two vertices and keep t
            found += comp.size() - 3;
            Vertex terminal = found_goal?terminals_[1]:terminals_[0];
            // compute C(1,Y)
            std::vector<Vertex> subset = comp;
            subset.push_back(separator[0]);
            std::swap(subset[0], subset.back());
            subset.push_back(separator[1]);
            std::swap(subset[1], subset.back());
            auto term_it = std::find(subset.begin(), subset.end(), terminal);
            size_t term_index = std::distance(subset.begin(), term_it);
            assert(term_index != subset.size());
            Graph comp_graph = subgraph(subset);
            comp_graph.max_length_ = max_length_;
            comp_graph.terminals_ = {0, (Vertex)term_index};
            comp_graph.preprocess();
            comp_graph.normalize();
            Search search(comp_graph);
            auto res = search.search();
            auto res_extra = comp_graph.extra_paths();
            res.resize(max_length_ + 1);
            res_extra.resize(max_length_ + 1);
            std::vector<Edge_weight> c_one_y;
            for(size_t length = 0; length < res.size(); length++) {
                c_one_y.push_back(res[length] + res_extra[length]);
            }
            // compute C(1,N)
            subset = comp;
            subset.push_back(separator[0]);
            std::swap(subset[0], subset.back());
            term_it = std::find(subset.begin(), subset.end(), terminal);
            term_index = std::distance(subset.begin(), term_it);
            assert(term_index != subset.size());
            comp_graph = subgraph(subset);
            comp_graph.max_length_ = max_length_;
            comp_graph.terminals_ = {0, (Vertex)term_index};
            comp_graph.preprocess();
            comp_graph.normalize();
            search = Search(comp_graph);
            res = search.search();
            res_extra = comp_graph.extra_paths();
            res.resize(max_length_ + 1);
            res_extra.resize(max_length_ + 1);
            std::vector<Edge_weight> c_one_n;
            for(size_t length = 0; length < res.size(); length++) {
                c_one_n.push_back(res[length] + res_extra[length]);
            }
            // compute C(2,Y)
            subset = comp;
            subset.push_back(separator[0]);
            std::swap(subset[0], subset.back());
            subset.push_back(separator[1]);
            std::swap(subset[1], subset.back());
            term_it = std::find(subset.begin(), subset.end(), terminal);
            term_index = std::distance(subset.begin(), term_it);
            assert(term_index != subset.size());
            comp_graph = subgraph(subset);
            comp_graph.max_length_ = max_length_;
            comp_graph.terminals_ = {1, (Vertex)term_index};
            comp_graph.preprocess();
            comp_graph.normalize();
            search = Search(comp_graph);
            res = search.search();
            res_extra = comp_graph.extra_paths();
            res.resize(max_length_ + 1);
            res_extra.resize(max_length_ + 1);
            std::vector<Edge_weight> c_two_y;
            for(size_t length = 0; length < res.size(); length++) {
                c_two_y.push_back(res[length] + res_extra[length]);
            }
            // compute C(2,N)
            subset = comp;
            subset.push_back(separator[1]);
            std::swap(subset[0], subset.back());
            term_it = std::find(subset.begin(), subset.end(), terminal);
            term_index = std::distance(subset.begin(), term_it);
            assert(term_index != subset.size());
            comp_graph = subgraph(subset);
            comp_graph.max_length_ = max_length_;
            comp_graph.terminals_ = {0, (Vertex)term_index};
            comp_graph.preprocess();
            comp_graph.normalize();
            search = Search(comp_graph);
            res = search.search();
            res_extra = comp_graph.extra_paths();
            res.resize(max_length_ + 1);
            res_extra.resize(max_length_ + 1);
            std::vector<Edge_weight> c_two_n;
            for(size_t length = 0; length < res.size(); length++) {
                c_two_n.push_back(res[length] + res_extra[length]);
            }
            // modify the component accordingly
            // remove all vertices in the component (apart from the separator, the terminal, and two vertices that we will reuse)
            term_it = std::find(comp.begin(), comp.end(), terminal);
            term_index = std::distance(comp.begin(), term_it);
            assert(term_index != comp.size());
            std::swap(comp[term_index], comp[2]);
            // we keep the first three: two vertices to reuse and the terminal (in the third position)
            for(size_t i = 3; i < comp.size(); i++) {
                remove_vertex(comp[i]);
            }
            // for the kept vertices we need to remove the edges though
            for(size_t i = 0; i < 3; i++) {
                remove_edge(Edge(comp[i], separator[0]));
                remove_edge(Edge(comp[i], separator[1]));
                for(size_t j = i + 1; j < 3; j++) {
                    remove_edge(Edge(comp[i], comp[j]));
                }
            }
            // now readd appropriate edges
            assert(c_one_y[0] == 0);
            assert(c_one_n[0] == 0);
            assert(c_two_y[0] == 0);
            assert(c_two_n[0] == 0);
            assert(c_one_y[1] - c_one_n[1] == 0);
            assert(c_two_y[1] - c_two_n[1] == 0);
            add_edge(Edge(comp[0], terminal), Weight(1,1));
            add_edge(Edge(comp[1], terminal), Weight(1,1));
            for(Edge_length length = 0; length <= max_length_; length++) {
                if(c_one_n[length] > 0) {
                    add_edge(Edge(terminal, separator[0]), Weight(length, c_one_n[length]));
                }
                if(c_two_n[length] > 0) {
                    add_edge(Edge(terminal, separator[1]), Weight(length, c_two_n[length]));
                }
                if(c_one_y[length] - c_one_n[length] > 0) {
                    add_edge(Edge(comp[0], separator[0]), Weight(length - 1, c_one_y[length] - c_one_n[length]));
                }
                if(c_two_y[length] - c_two_n[length] > 0) {
                    add_edge(Edge(comp[1], separator[1]), Weight(length - 1, c_two_y[length] - c_two_n[length]));
                }
            }
            add_exclude(comp[0], separator[1]);
            add_exclude(comp[1], separator[0]);
            continue;
        } 
        // we have to enter AND leave the component, meaning we to traverse s_1 -> G[comp] -> s_2 (or the other way around)
        // thus we can reduce to {s_1, s_2}
        found += comp.size();
        comp.push_back(separator[0]);
        std::swap(comp[0], comp.back());
        comp.push_back(separator[1]);
        std::swap(comp[1], comp.back());
        Graph comp_graph = subgraph(comp);
        Edge_length d1 = 0, d2 = 0;
        if(max_length_ > distance_from_start[separator[0]] + distance_to_goal[separator[1]]) {
            d1 = max_length_ - distance_from_start[separator[0]] - distance_to_goal[separator[1]];
        }
        if(max_length_ > distance_from_start[separator[1]] + distance_to_goal[separator[0]]) {
            d2 = max_length_ - distance_from_start[separator[1]] - distance_to_goal[separator[0]];
        }
        comp_graph.max_length_ = std::max(d1, d2);
        comp_graph.terminals_ = {0, 1};
        comp_graph.print_stats();
        comp_graph.preprocess();
        comp_graph.normalize();
        // comp_graph.print_stats();
        Search search(comp_graph);
        auto res = search.search();
        // search.print_stats();
        auto res_extra = comp_graph.extra_paths();
        res.resize(std::max(res.size(), res_extra.size()));
        res_extra.resize(std::max(res.size(), res_extra.size()));
        // remove all vertices in the component (apart from the separator)
        for(size_t i = 2; i < comp.size(); i++) {
            remove_vertex(comp[i]);
        }
        // remove all previous edges between the separator vertices
        // we covered them in the result of the search
        remove_edge(Edge(separator[0], separator[1]));
        assert(res[0] == 0);
        assert(res_extra[0] == 0);
        for(Edge_length length = 1; length < res.size(); length++) {
            Edge_weight weight = res[length] + res_extra[length];
            if(weight > 0) {
                add_edge(Edge(separator[0], separator[1]), Weight(length, weight));
            }
        }
    }
    return found;
}