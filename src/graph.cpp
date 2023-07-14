// TLDC "Too long; Didn't Count" A length limited path counter.
// Copyright (C) 2023 Rafael Kiesel, Markus Hecher

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

#include "graph.h"
#include "parallel_search.h"
#include <clingo.hh>
#include <algorithm>
#include <map>
#include <limits>
#include <sstream>
#include "nauty2_8_6/gtools.h"
#include "math.h"

namespace fpc {

Graph::Graph(std::istream &input) {
    char dec;
    input >> dec;
    Vertex nr_vertices;
    Vertex nr_edges;
    std::string line;
    all_pair_ = true;
    max_length_ = std::numeric_limits<Edge_length>::max();
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
            adjacency_ = std::vector<std::vector<std::map<Edge_length,Edge_weight>>>(nr_vertices, 
                                    std::vector<std::map<Edge_length,Edge_weight>>(nr_vertices, 
                                                            std::map<Edge_length,Edge_weight>()));
            neighbors_ = std::vector<std::set<Vertex>>(nr_vertices, std::set<Vertex>());
            exclusion_classes_ = std::vector<std::set<Vertex>>(nr_vertices, std::set<Vertex>());
            exclude_ = std::vector<size_t>(nr_vertices);
            for(size_t i = 0; i < nr_vertices; i++) {
                exclusion_classes_[i].insert(i);
                exclude_[i] = i;
            }
            break;
        }
        case 'e':
            Vertex v,w;
            input >> v >> w;
            add_edge(Edge(v - 1,w - 1), Weight(1,1));
            break;
        case 'l':
            size_t tmp;
            input >> tmp;
            assert(tmp <= 127);
            max_length_ = Edge_length(tmp);
            extra_paths_ = std::vector<Edge_weight>(max_length_ + 1, 0);
            break;
        case 't':
            terminals_.resize(2);
            input >> terminals_[0] >> terminals_[1];
            terminals_[0]--;
            terminals_[1]--;
            all_pair_ = false;
            break;
        default:
            std::cerr << "Invalid character " << dec << " at beginning of line." << std::endl;
            break;
        }
        input >> dec;
    }
    assert(max_length_ + 1 == extra_paths_.size());
    assert(adjacency_.size() > 0);
    assert(all_pair_ || terminals_.size() == 2);
    // max_length_ = Edge_length(std::min(nr_vertices - 1, int(max_length_)));
    extra_paths_ = std::vector<Edge_weight>(max_length_ + 1, 0);
}

void Graph::preprocess() {
    bool found = true;
    while(found) {
        found = false;
        preprocess_start_goal_edges();
        Vertex cur_isolated_removed = preprocess_isolated();
        found |= cur_isolated_removed > 0;
        isolated_removed += cur_isolated_removed;
        preprocess_start_goal_edges();
        Vertex cur_unreachable_removed = preprocess_unreachable();
        found |= cur_unreachable_removed > 0;
        unreachable_removed += cur_unreachable_removed;
        if(!found) {
            preprocess_start_goal_edges();
            Vertex cur_unusable_edge_removed = preprocess_unusable_edge();
            found |= cur_unusable_edge_removed > 0;
            unusable_edge_removed += cur_unusable_edge_removed;
        }
        if(!found && max_length_decrease == 0) {
            preprocess_start_goal_edges();
            Vertex cur_max_length_decrease = limit_max_length();
            found |= cur_max_length_decrease > 0;
            max_length_decrease += cur_max_length_decrease;
        }
    }
    preprocess_start_goal_edges();
    assert(all_pair_ || neighbors(terminals_[0]).count(terminals_[1]) == 0);
    assert(all_pair_ || neighbors(terminals_[1]).count(terminals_[0]) == 0);
}

void Graph::print_stats() {
    Vertex nr_edges = 0;
    for(Vertex v = 0; v < adjacency_.size(); v++) {
        nr_edges += neighbors(v).size();
    }
    nr_edges /= 2;
    std::vector<Edge_length> distance_to_goal(adjacency_.size(), std::numeric_limits<Edge_length>::max());
    if(!all_pair_) {
        dijkstra(terminals_[1], distance_to_goal, {});
    }
    if(isolated_removed)            std::cerr << "Removed isolated: " << isolated_removed << std::endl;
    if(forwarder_removed)           std::cerr << "Removed forwarder: " << forwarder_removed << std::endl;
    if(position_determined_removed) std::cerr << "Removed position determined: " << position_determined_removed << std::endl;
    if(twin_edges_removed)          std::cerr << "Removed twin edges: " << twin_edges_removed << std::endl;
    if(unreachable_removed)         std::cerr << "Removed unreachable: " << unreachable_removed << std::endl;
    if(unusable_edge_removed)       std::cerr << "Removed unusable edge: " << unusable_edge_removed << std::endl;
    if(two_sep_removed)             std::cerr << "Removed due to 2-separation: " << two_sep_removed << std::endl;
    if(three_sep_removed)           std::cerr << "Removed due to 3-separation: " << three_sep_removed << std::endl;
    if(max_length_decrease)         std::cerr << "Max length decreased by: " << max_length_decrease << std::endl;
    std::cerr << "#vertices " << adjacency_.size() << " #edges " << nr_edges;
    if(!all_pair_) {
        std::cerr << " max. length " << static_cast<size_t>(max_length_) << " min. length " << static_cast<size_t>(distance_to_goal[terminals_[0]]) << std::endl;
        std::cerr << "terminals: " << terminals_[0] << "," << terminals_[1] << std::endl;
    } else {
        std::cerr << std::endl;
    }
}

Edge_length Graph::min_length() {
    if(all_pair_) {
        return 1;
    }
    std::vector<Edge_length> distance_to_goal(adjacency_.size(), std::numeric_limits<Edge_length>::max());
    dijkstra(terminals_[1], distance_to_goal, {});
    return distance_to_goal[terminals_[0]];
}

void Graph::normalize(bool reorder) {
    Vertex unnamed = std::numeric_limits<Vertex>::max();
    std::vector<Vertex> new_name(adjacency_.size(), unnamed);
    Vertex cur_name = 0;
    if(reorder && !all_pair_) {
        new_name[terminals_[0]] = cur_name++;
        new_name[terminals_[1]] = cur_name++;
    }
    for(Vertex v = 0; v < adjacency_.size(); v++) {
        if(neighbors(v).empty()) {
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
        if(neighbors(v).empty()) {
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
    if(!all_pair_) {
        assert(new_name[terminals_[0]] != unnamed);
        assert(new_name[terminals_[1]] != unnamed);
        terminals_[0] = new_name[terminals_[0]];
        terminals_[1] = new_name[terminals_[1]];
    }
}

size_t Graph::nr_vertices() {
    size_t ret = 0;
    for(Vertex v = 0; v < adjacency_.size(); v++) {
        if(!neighbors(v).empty()) {
            ret++;
        }
    }
    return ret;
}

size_t Graph::nr_edges() {
    size_t nr_edges = 0;
    for(Vertex v = 0; v < adjacency_.size(); v++) {
        nr_edges += neighbors(v).size();
    }
    assert(nr_edges % 2 == 0);
    return nr_edges/2;
}

sparsegraph Graph::to_canon_nauty(bool reorder) {
    normalize(reorder);
    DEFAULTOPTIONS_SPARSEGRAPH(options);
    options.getcanon = true;
    options.defaultptn = false;
    SG_DECL(sg);
    int m = SETWORDSNEEDED(adjacency_.size());
    nauty_check(WORDSIZE,m,adjacency_.size(),NAUTYVERSIONID);
    size_t nr_edges = 0;
    for(Vertex v = 0; v < adjacency_.size(); v++) {
        nr_edges += neighbors(v).size();
    }
    sg.v = (edge_t *)malloc(sizeof(edge_t)*(adjacency_.size() + nr_edges) + sizeof(degree_t)*adjacency_.size());
    sg.d = (degree_t *)(sg.v + adjacency_.size());
    sg.e = (edge_t *)(sg.d + adjacency_.size());
    sg.nv = adjacency_.size();
    sg.nde = nr_edges;
    sg.vlen = adjacency_.size();
    sg.dlen = adjacency_.size();
    sg.elen = nr_edges;
    nr_edges = 0;
    for(Vertex v = 0; v < adjacency_.size(); v++) {
        sg.v[v] = nr_edges;
        for(Vertex w : neighbors(v)) {
            sg.e[nr_edges++] = w;
        }
        sg.d[v] = neighbors(v).size();
    }
    int *lab = (int *)malloc(sg.nv*sizeof(int));
    int *ptn = (int *)malloc(sg.nv*sizeof(int));
    int *orbits = (int *)malloc(sg.nv*sizeof(int));
    ptn[0] = 0;
    ptn[1] = 0;
    lab[0] = 0;
    lab[1] = 1;
    for(Vertex v = 2; v < adjacency_.size(); v++) {
        ptn[v] = 1;
        lab[v] = v;
    }
    statsblk stats;
    SG_DECL(canon_sg);
    canon_sg.v = (edge_t *)malloc(sizeof(edge_t)*(adjacency_.size() + nr_edges) + sizeof(degree_t)*adjacency_.size());
    canon_sg.d = (degree_t *)(canon_sg.v + adjacency_.size());
    canon_sg.e = (edge_t *)(canon_sg.d + adjacency_.size());
    canon_sg.nv = adjacency_.size();
    canon_sg.nde = nr_edges;
    canon_sg.vlen = sg.vlen;
    canon_sg.dlen = sg.dlen;
    canon_sg.elen = sg.elen;
    sparsenauty(&sg,lab,ptn,orbits,&options,&stats,&canon_sg);
    sortlists_sg(&canon_sg);
    free(sg.v);
    free(lab);
    free(ptn);
    free(orbits);
    return canon_sg;
}


double Graph::nr_automorphisms() {
    normalize(false);
    DEFAULTOPTIONS_SPARSEGRAPH(options);
    options.getcanon = true;
    options.defaultptn = false;
    SG_DECL(sg);
    int m = SETWORDSNEEDED(adjacency_.size());
    nauty_check(WORDSIZE,m,adjacency_.size(),NAUTYVERSIONID);
    size_t nr_edges = 0;
    for(Vertex v = 0; v < adjacency_.size(); v++) {
        nr_edges += neighbors(v).size();
    }
    sg.v = (edge_t *)malloc(sizeof(edge_t)*(adjacency_.size() + nr_edges) + sizeof(degree_t)*adjacency_.size());
    sg.d = (degree_t *)(sg.v + adjacency_.size());
    sg.e = (edge_t *)(sg.d + adjacency_.size());
    sg.nv = adjacency_.size();
    sg.nde = nr_edges;
    sg.vlen = adjacency_.size();
    sg.dlen = adjacency_.size();
    sg.elen = nr_edges;
    nr_edges = 0;
    for(Vertex v = 0; v < adjacency_.size(); v++) {
        sg.v[v] = nr_edges;
        for(Vertex w : neighbors(v)) {
            sg.e[nr_edges++] = w;
        }
        sg.d[v] = neighbors(v).size();
    }
    int *lab = (int *)malloc(sg.nv*sizeof(int));
    int *ptn = (int *)malloc(sg.nv*sizeof(int));
    int *orbits = (int *)malloc(sg.nv*sizeof(int));
    ptn[0] = 0;
    ptn[1] = 0;
    lab[0] = 0;
    lab[1] = 1;
    for(Vertex v = 2; v < adjacency_.size(); v++) {
        ptn[v] = 1;
        lab[v] = v;
    }
    statsblk stats;
    SG_DECL(canon_sg);
    canon_sg.v = (edge_t *)malloc(sizeof(edge_t)*(adjacency_.size() + nr_edges) + sizeof(degree_t)*adjacency_.size());
    canon_sg.d = (degree_t *)(canon_sg.v + adjacency_.size());
    canon_sg.e = (edge_t *)(canon_sg.d + adjacency_.size());
    canon_sg.nv = adjacency_.size();
    canon_sg.nde = nr_edges;
    canon_sg.vlen = sg.vlen;
    canon_sg.dlen = sg.dlen;
    canon_sg.elen = sg.elen;
    sparsenauty(&sg,lab,ptn,orbits,&options,&stats,&canon_sg);
    sortlists_sg(&canon_sg);
    free(sg.v);
    free(lab);
    free(ptn);
    free(orbits);
    free(canon_sg.v);
    return stats.grpsize1*std::pow(10,stats.grpsize2);
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
    if(exclude_[v] == exclude_[w]) {
        assert(exclusion_classes_[exclude_[v]].count(w) > 0);
        assert(exclusion_classes_[exclude_[v]].count(v) > 0);
        return;
    }
    exclusion_classes_[exclude_[v]].insert(exclusion_classes_[exclude_[w]].begin(), exclusion_classes_[exclude_[w]].end());
    auto prev = exclude_[w];
    for(auto other : exclusion_classes_[exclude_[w]]) {
        exclude_[other] = exclude_[v];
    }
    exclusion_classes_[prev].clear();
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
    exclusion_classes_[exclude_[v]].erase(v);
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
                            exclusion_classes_(n),
                            exclude_(n) {
    for(size_t i = 0; i < n; i++) {
        exclusion_classes_[i].insert(i);
        exclude_[i] = i;
    }
}

Graph Graph::subgraph(std::vector<Vertex> restrict_to) {
    Graph ret(restrict_to.size());
    ret.max_length_ = max_length_;
    auto unnamed = std::numeric_limits<Vertex>::max();
    std::vector<Vertex> new_name(adjacency_.size(), unnamed);
    Vertex cur_name = 0;
    for(Vertex v : restrict_to) {
        new_name[v] = cur_name++;
    }
    ret.exclusion_classes_ = std::vector<std::set<Vertex>>();
    ret.exclude_ = std::vector<size_t>(cur_name, -1);
    Vertex cur_exclude = 0;
    std::vector<Vertex> new_exclude_name(exclusion_classes_.size(), unnamed);
    for(Vertex v : restrict_to) {
        for(Vertex neigh : neighbors(v)) {
            if(new_name[neigh] != unnamed) {
                ret.neighbors_[new_name[v]].insert(new_name[neigh]);
                ret.adjacency_[new_name[v]][new_name[neigh]] = adjacency_[v][neigh];
            }
        }
        if(new_exclude_name[exclude_[v]] == unnamed) {
            ret.exclusion_classes_.push_back({});
            new_exclude_name[exclude_[v]] = cur_exclude++;
            for(auto excluded : exclusion_classes_[exclude_[v]]) {
                assert(new_name[excluded] != unnamed);
                ret.exclusion_classes_[new_exclude_name[exclude_[v]]].insert(new_name[excluded]);
            }
        }
        ret.exclude_[new_name[v]] = new_exclude_name[exclude_[v]];
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
            assert(!adjacency_[v].empty() || !fixed(v));
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
            for(auto w : exclusion_classes_[exclude_[cur]]) {
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
            for(auto excluded : exclusion_classes_[exclude_[v]]) {
                prog_str << "e(" << v << "," << excluded << ").\n";
            }
        }
    }
    if(terminals_in_same) {
        prog_str << ":- r(" << terminals_[0] << "), not sep(" << terminals_[0] << ").\n";
        prog_str << ":- r(" << terminals_[1] << "), not sep(" << terminals_[1] << ").\n";
    }
    for(Vertex v = 0; v < adjacency_.size(); v++) {
        for(Vertex w : neighbors(v)) {
            prog_str << "e(" << v << "," << w << ").\n";
        }
    }
    prog_str << "#show sep/1.\n";
    // initialize clingo
    Clingo::Logger logger = [](Clingo::WarningCode, char const *) {
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
    if(all_pair_) {
        return 0;
    }
    assert(terminals_.size() == 2);
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
    if(all_pair_) {
        return 0;
    }
    assert(terminals_.size() == 2);
    Vertex found = 0;
    for(size_t v = 0; v < adjacency_.size(); v++) {
        auto cur_neighbors = neighbors(v);
        if(cur_neighbors.size() == 1) {
            Vertex neighbor = *cur_neighbors.begin();
            if(terminals_[0] == v || terminals_[1] == v) {
                assert(neighbor != terminals_[0] && neighbor != terminals_[1]);
                assert(max_length_ > 0);
                found++;
                remove_vertex(v);
                if(terminals_[0] == v) {
                    terminals_[0] = neighbor;
                } else {
                    terminals_[1] = neighbor;
                }
                max_length_--;
            } else {
                found++;
                remove_vertex(v);
            }
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
    if(all_pair_) {
        return 0;
    }
    assert(terminals_.size() == 2);
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
    if(all_pair_) {
        return 0;
    }
    assert(terminals_.size() == 2);
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
    std::vector<Vertex> separator = find_separator(2, 1, true);
    // if(separator.size() == 0) {
    //     separator = find_separator(2, 4, false);
    // }
    if(separator.size() == 0) {
        return 0;
    }
    assert(separator.size() == 2);
    assert(!fixed(separator[0]));
    assert(!fixed(separator[1]));
    std::map<Vertex, std::set<Vertex>> to_exclude;
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
        if(found_goal || found_start) {
            // nothing we can do
            continue;
        } 
        if(found_goal || found_start) {
            if(comp.size() < 4) {
                // it could happen that we have multiple components that are smaller 
                // but whose size adds up to more than 4
                // we could in principle look at them together
                continue;
            }
            // let t be the terminal found, and s_1 and s_2 be the separators
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
            // first bound the length of the path within the component
            Edge_length c_one_length, c_two_length;
            if(found_goal) {
                // if we found the goal in this component
                // we first need to get from the start to the separator
                assert(distance_from_start[separator[0]] <= max_length_);
                assert(distance_from_start[separator[1]] <= max_length_);
                c_one_length = max_length_ - distance_from_start[separator[0]];
                c_two_length = max_length_ - distance_from_start[separator[1]];
            } else {
                // if we found the start in this component
                // we first need to get from the goal to the separator
                assert(distance_to_goal[separator[0]] <= max_length_);
                assert(distance_to_goal[separator[1]] <= max_length_);
                c_one_length = max_length_ - distance_to_goal[separator[0]];
                c_two_length = max_length_ - distance_to_goal[separator[1]];
            }
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
            if(neighbors(separator[0]).count(separator[1]) > 0) {
                comp_graph.remove_edge(Edge(0,1));
            }
            comp_graph.max_length_ = c_one_length;
            comp_graph.terminals_ = {0, (Vertex)term_index};
            comp_graph.preprocess();
            comp_graph.normalize();
            ParallelSearch search(comp_graph.to_canon_nauty(true), c_one_length, 12);
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
            comp_graph.max_length_ = c_one_length;
            comp_graph.terminals_ = {0, (Vertex)term_index};
            comp_graph.preprocess();
            comp_graph.normalize();
            search = ParallelSearch(comp_graph.to_canon_nauty(true), c_one_length, 12);
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
            if(neighbors(separator[0]).count(separator[1]) > 0) {
                comp_graph.remove_edge(Edge(0,1));
            }
            comp_graph.max_length_ = c_two_length;
            comp_graph.terminals_ = {1, (Vertex)term_index};
            comp_graph.preprocess();
            comp_graph.normalize();
            search = ParallelSearch(comp_graph.to_canon_nauty(true), c_two_length, 12);
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
            comp_graph.max_length_ = c_two_length;
            comp_graph.terminals_ = {0, (Vertex)term_index};
            comp_graph.preprocess();
            comp_graph.normalize();
            search = ParallelSearch(comp_graph.to_canon_nauty(true), c_two_length, 12);
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
                exclusion_classes_[exclude_[comp[i]]].erase(comp[i]);
                exclusion_classes_.push_back({comp[i]});
                exclude_[comp[i]] = exclusion_classes_.size() - 1;
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
            if(comp[0] < separator[1]) {
                to_exclude[comp[0]].insert(separator[1]);
            } else {
                to_exclude[separator[1]].insert(comp[0]);
            }
            if(comp[1] < separator[0]) {
                to_exclude[comp[1]].insert(separator[0]);
            } else {
                to_exclude[separator[0]].insert(comp[1]);
            }
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
        comp_graph.preprocess();
        comp_graph.normalize();
        // comp_graph.print_stats();
        ParallelSearch search(comp_graph.to_canon_nauty(true), comp_graph.max_length_, 12);
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
        // there cannot be any exclusion constraints on the separator
        // so we do not need to remove those here
        assert(res[0] == 0);
        assert(res_extra[0] == 0);
        for(Edge_length length = 1; length < res.size(); length++) {
            Edge_weight weight = res[length] + res_extra[length];
            if(weight > 0) {
                add_edge(Edge(separator[0], separator[1]), Weight(length, weight));
            }
        }
    }
    for(auto &[ex_1, ex_set] : to_exclude) {
        for(auto ex_2 : ex_set) {
            add_exclude(ex_1, ex_2);
        }
    }
    return found;
}


Vertex Graph::preprocess_three_separator() {
    Vertex found = 0;
    // here we can only do something if neither start nor goal are in the component
    // furthermore, we at least need 
    std::vector<Vertex> separator = find_separator(3, 7, true);
    if(separator.size() == 0) {
        return 0;
    }
    assert(separator.size() == 3);
    assert(!fixed(separator[0]));
    assert(!fixed(separator[1]));
    assert(!fixed(separator[2]));
    std::map<Vertex, std::set<Vertex>> to_exclude;
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
        if(found_goal || found_start) {
            // nothing we can do
            continue;
        } 
        found += comp.size() - 6;
        // let s_0, s_1, and s_2 be the separators
        // solve four subqueries:
        // counts[i][j] for i = 0, 1, 2, j = 0, 1
        // if j == 0, counts[i][j] contains the number of paths 
        // between the two other separators
        // in the component that may use separator[i]
        // if j == 1, counts[i][j] contains the number of paths 
        // between the two other separators
        // in the component that may *not* use separator[i]
        // then we modify the graph by replacing the component with 
        //     e_0-----s_0-----e_5
        //      |      / \      |
        //      |     /   \     |
        //      |    /     \    |
        //      |  e_1     e_4  |
        //      |  /         \  |
        //      | / /--e_3--\ \ |
        //      |/ /         \ \|
        //      s_1           s_2
        //         \         /
        //          \--e_2--/
        // where:
        // for i = 0, 1, 2:
        //      {s_i,e_2*i}, {s_i,e_2*1 + 1}    have (0,1)
        //      {s_i,e_2*(i-1)}                 has counts[i+1][1] 
        //      {s_i,e_2*(i-1) + 1}             has counts[i+1][0] - counts[i+1][1]
        // and 
        // e_1 and s_2 exclude each other (implies that e_1 and all other e_i's exclude each other)
        // e_3 and s_0 exclude each other (implies that e_3 and all other e_i's exclude each other)
        // e_5 and s_1 exclude each other (implies that e_5 and all other e_i's exclude each other)
        // e_0 and all other e_i's exclude each other (e_0 in principle already excludes e_1,e_3,e_5 so excluding e_2,e_4 is enough)
        // e_2 and all other e_i's exclude each other (e_2 in principle already excludes e_1,e_3,e_5 so excluding e_0,e_4 is enough)
        // e_4 and all other e_i's exclude each other (e_4 in principle already excludes e_1,e_3,e_5 so excluding e_0,e_2 is enough)
        std::vector<Edge_weight> counts[3][2];
        for(size_t i = 0; i < 3; i++) {
            // separator[i] is not a terminal
            for(size_t j = 0; j < 2; j++) {
                // we may use separator[i] iff j == 0
                // compute counts[i][j]
                // construct the subgraph
                std::vector<Vertex> subset = comp;
                size_t start_idx = (i + 1) % 3;
                size_t goal_idx = (i + 2) % 3;
                subset.push_back(separator[start_idx]);
                std::swap(subset[0], subset.back());
                subset.push_back(separator[goal_idx]);
                std::swap(subset[1], subset.back());
                if(j == 0) {
                    subset.push_back(separator[i]);
                }
                Graph comp_graph = subgraph(subset);
                if(neighbors(separator[start_idx]).count(separator[goal_idx])) {
                    comp_graph.remove_edge(Edge(0,1));
                }
                if(j == 0) {
                    if(neighbors(separator[start_idx]).count(separator[i])) {
                        comp_graph.remove_edge(Edge(0,subset.size() - 1));
                    }
                    if(neighbors(separator[goal_idx]).count(separator[i])) {
                        comp_graph.remove_edge(Edge(1,subset.size() - 1));
                    }
                }
                // restrict local maximum length
                Edge_length d1 = 0, d2 = 0;
                if(max_length_ > distance_from_start[separator[start_idx]] + distance_to_goal[separator[goal_idx]]) {
                    d1 = max_length_ - distance_from_start[separator[start_idx]] - distance_to_goal[separator[goal_idx]];
                }
                if(max_length_ > distance_from_start[separator[goal_idx]] + distance_to_goal[separator[start_idx]]) {
                    d2 = max_length_ - distance_from_start[separator[goal_idx]] - distance_to_goal[separator[start_idx]];
                }
                comp_graph.max_length_ = std::max(d1, d2);
                comp_graph.terminals_ = {0, 1};
                comp_graph.preprocess();
                comp_graph.normalize();
                ParallelSearch search(comp_graph.to_canon_nauty(true), comp_graph.max_length_, 12);
                auto res = search.search();
                auto res_extra = comp_graph.extra_paths();
                res.resize(max_length_ + 1);
                res_extra.resize(max_length_ + 1);
                for(size_t length = 0; length < res.size(); length++) {
                    counts[i][j].push_back(res[length] + res_extra[length]);
                }
                assert(counts[i][j][0] == 0);
                if(j == 1) {
                    assert(counts[i][0][1] - counts[i][1][1] == 0);
                }
            }
        }
        // modify the graph accordingly
        // remove all vertices in the component (apart from the separator, the terminal, and two vertices that we will reuse)
        // we keep the first six to reuse
        for(size_t i = 6; i < comp.size(); i++) {
            remove_vertex(comp[i]);
        }
        // for the kept vertices we need to remove the edges though
        for(size_t i = 0; i < 6; i++) {
            exclusion_classes_[exclude_[comp[i]]].erase(comp[i]);
            exclusion_classes_.push_back({comp[i]});
            exclude_[comp[i]] = exclusion_classes_.size() - 1;
            remove_edge(Edge(comp[i], separator[0]));
            remove_edge(Edge(comp[i], separator[1]));
            remove_edge(Edge(comp[i], separator[2]));
            for(size_t j = i + 1; j < 6; j++) {
                remove_edge(Edge(comp[i], comp[j]));
            }
        }
        // now readd appropriate edges
        // first half, unweighted
        for(size_t i = 0; i < 3; i++) {
            add_edge(Edge(comp[2*i], separator[i]), Weight(1,1));
            add_edge(Edge(comp[2*i + 1], separator[i]), Weight(1,1));
        }
        // second half, weighted
        for(size_t i = 0; i < 3; i++) {
            assert(!counts[i][1][1]);
        }
        for(Edge_length length = 2; length <= max_length_; length++) {
            for(size_t i = 0; i < 3; i++) {
                if(counts[(i + 1) % 3][1][length] > 0) {
                    add_edge(
                                Edge(separator[i], comp[(2*(i + 2)) % 6]),
                                Weight(length - 1, counts[(i + 1) % 3][1][length])
                            );
                }
                if(counts[(i + 1) % 3][0][length] - counts[(i + 1) % 3][1][length] > 0) {
                    add_edge(
                                Edge(separator[i], comp[(2*(i + 2) + 1) % 6]),
                                Weight(length - 1, counts[(i + 1) % 3][0][length] - counts[(i + 1) % 3][1][length])
                            );
                }
            }
        }
        // add exclusion constraints
        for(size_t i = 0; i < 3; i++) {
            auto comp_idx = (2*i + 3) % 6;
            if(comp[comp_idx] < separator[i]) {
                to_exclude[comp[comp_idx]].insert(separator[i]);
            } else {
                to_exclude[separator[i]].insert(comp[comp_idx]);
            }
        }
        for(size_t i = 1; i < 2; i++) {
            auto comp_idx = 2*i;
            if(comp[comp_idx] < comp[0]) {
                to_exclude[comp[comp_idx]].insert(comp[0]);
            } else {
                to_exclude[comp[0]].insert(comp[comp_idx]);
            }
        }
        if(comp[2] < comp[4]) {
            to_exclude[comp[2]].insert(comp[4]);
        } else {
            to_exclude[comp[4]].insert(comp[2]);
        }
    }
    for(auto &[ex_1, ex_set] : to_exclude) {
        for(auto ex_2 : ex_set) {
            add_exclude(ex_1, ex_2);
        }
    }
    return found;
}

Vertex Graph::limit_max_length() {
    if(all_pair_) {
        return 0;
    }
    assert(terminals_.size() == 2);
    // build the program
    std::stringstream prog_str;
    prog_str << "reach(X, 0) :- start(X).\n";
    // prog_str << ":~ goal(X), reach(X,Y). [-Y]\n";
    prog_str << "start(" << terminals_[0] << ").\n";
    prog_str << "goal(" << terminals_[1] << ").\n";
    prog_str << ":- reach(X, L), reach(X, L'), L != L'.\n";
    std::vector<Edge_length> distance_from_start(adjacency_.size(), std::numeric_limits<Edge_length>::max());
    dijkstra(terminals_[0], distance_from_start, {});

    std::vector<Edge_length> distance_to_goal(adjacency_.size(), std::numeric_limits<Edge_length>::max());
    dijkstra(terminals_[1], distance_to_goal, {});
    if(distance_to_goal[terminals_[0]] > max_length_) {
        return 0;
    }

    for(Vertex v = 0; v < adjacency_.size(); v++) {
        if(v == terminals_[1] || adjacency_[v].empty()) {
            continue;
        }
        std::vector<std::pair<Vertex,Edge_length>> edges;
        for(auto &w : neighbors(v)) {
            if(w == terminals_[0]) {
                continue;
            }
            for(auto &weight : adjacency_[v][w]) {
                edges.push_back(std::make_pair(w, weight.first));
            }
        }
        for(auto &[w, length] : edges) {
            prog_str << "reach(" << w << ",L + L') :- reach(" << v << ", L), edge(" << v << "," << w << ",L'), L <= " << static_cast<size_t>(max_length_ - distance_to_goal[w] - length);
            prog_str << ", L>= " << static_cast<size_t>(distance_from_start[v]) << ".\n";
            prog_str << "edge(" << v << "," << w << "," << static_cast<size_t>(length) << ") :- reach(" << v << ", L), L <= " << static_cast<size_t>(max_length_ - distance_to_goal[w] - length);
            prog_str << ", L>= " << static_cast<size_t>(distance_from_start[v]);
            for(auto &[wp, lengthp] : edges) {
                if(w != wp || length != lengthp) {
                    prog_str << ", not edge(" << v << "," << wp << "," << static_cast<size_t>(lengthp) << ")";
                }
            }
            prog_str << ".\n";
        }
    }
    // std::cerr << prog_str.str();
    // return 0;
    // prog_str << "#show sep/1.\n";
    // initialize clingo
    Clingo::Logger logger = [](Clingo::WarningCode, char const *) {
    };
    Clingo::Control ctl{{}, logger, 20};
    ctl.add("base", {}, prog_str.str().c_str());
    ctl.ground({{"base", {}}});
    ctl.configuration()["solve"]["solve_limit"] = "1000";
    int64_t best_cost = max_length_;
    bool impossible = true;
    while(impossible) {
        impossible = false;
        auto atom = Clingo::Function("reach", {Clingo::Number(terminals_[1]), Clingo::Number(best_cost)});
        auto atom_it = ctl.symbolic_atoms().find(atom);
        if(atom_it == ctl.symbolic_atoms().end()) {
            impossible = true;
            assert(best_cost > 0);
            best_cost--;
            continue;
        }
        auto handle = ctl.solve(Clingo::LiteralSpan{(*atom_it).literal()});
        if(handle.get().is_unsatisfiable()) {
            impossible = true;
            assert(best_cost > 0);
            best_cost--;
        }
    }
    Edge_length prev_max = max_length_;
    max_length_ = best_cost;
    return prev_max - max_length_;
}

} // namespace fpc
