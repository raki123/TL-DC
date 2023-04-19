#include "graph.h"
#include <assert.h>


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
    while(found) {
        found = false;
        Vertex cur_isolated_removed = preprocess_isolated();
        found |= cur_isolated_removed > 0;
        isolated_removed += cur_isolated_removed;
    }
    std::cerr << "Removed isolated: " << isolated_removed << std::endl;
}

void Graph::add_edge(Edge edge, Weight weight) {
    if(adjacency_[edge.first].find(edge.second) == adjacency_[edge.first].end()) {
        assert(adjacency_[edge.second].find(edge.first) == adjacency_[edge.second].end());
        adjacency_[edge.first][edge.second] = {};
        adjacency_[edge.second][edge.first] = {};
    }
    adjacency_[edge.first][edge.second].push_back(weight);
    adjacency_[edge.second][edge.first].push_back(weight);
}

void Graph::remove_vertex(Vertex v) {
    for(auto neighbor : neighbors(v)) {
        adjacency_[neighbor].erase(v);
        adjacency_[v].clear();
    }
}

std::vector<Vertex> Graph::neighbors(Vertex v) {
    std::vector<Vertex> ret;
    for(auto & [w, weight] : adjacency_[v]) {
        ret.push_back(w);
    }
    return ret;
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
                for(auto & [w, old_weights] : adjacency_[neighbor]) {
                    if(w == v) {
                        continue;
                    }
                    std::vector<Weight> new_weights(old_weights.size() * adjacency_[v][neighbor].size());
                    for(auto old_weight : old_weights) {
                        for(auto weight : adjacency_[v][neighbor]) {
                            new_weights.push_back(Weight(old_weight.first + weight.first, old_weight.second * weight.second));
                        }
                    }
                    old_weights = new_weights;
                }
            }
            remove_vertex(v);
        }
    }
    return found;
}