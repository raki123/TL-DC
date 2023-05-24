#include "parallel_search.h"
#include <queue>
#include <limits>


ParallelSearch::ParallelSearch(Graph& input) :  
                                nthreads_(OMP_NUM_THREADS),
                                enable_dag_(true),
                                max_length_(input.max_length_),
                                terminals_(input.terminals_),
                                neighbors_(input.neighbors_.size()),
                                adjacency_(input.adjacency_.size(), std::vector<std::vector<std::pair<Edge_length, Edge_weight>>>(input.adjacency_.size())),
                                invalid_(std::numeric_limits<Edge_length>::max() - max_length_ - 1),
                                distance_(adjacency_.size(), std::vector<Edge_length>(adjacency_.size(), invalid_)),
                                visited_(adjacency_.size(), false),
                                cache_( 
                                    max_length_, 
                                    std::vector<std::unordered_map<CacheKey, std::vector<Edge_weight>, vector_hash>>(
                                        adjacency_.size(), 
                                        std::unordered_map<CacheKey, std::vector<Edge_weight>, vector_hash>()
                                    )
                                ),
                                result_(max_length_ + 1, 0),
                                thread_local_result_(
                                    nthreads_,
                                    std::vector<Edge_weight>(max_length_ + 1, 0)
                                ),
                                pos_hits_(OMP_NUM_THREADS, 0),
                                neg_hits_(OMP_NUM_THREADS, 0),
                                edges_(OMP_NUM_THREADS, 0),
                                propagations_(OMP_NUM_THREADS, 0),
                                dags_(OMP_NUM_THREADS, 0)  {
    assert(terminals_.size() == 2);
    for(Vertex v = 0; v < adjacency_.size();v++) {
        // fill neighbors
        neighbors_[v] = std::vector<Vertex>(input.neighbors_[v].begin(), input.neighbors_[v].end());
        for(Vertex w : neighbors_[v]) {
            adjacency_[v][w] = std::vector<std::pair<Edge_length, Edge_weight>>(input.adjacency_[v][w].begin(), input.adjacency_[v][w].end());
        }
        // for(Vertex w : input.exclusion_classes_[input.exclude_[v]]) {
        //     if(w != v) {
        //         enable_dag_ = false;
        //         exclude_[v].push_back(w);
        //     }
        // }
    }
    for(Vertex v = 0; v < adjacency_.size(); v++) {
        dijkstra(v, distance_[v]);
    }
}

std::vector<Edge_weight> ParallelSearch::search() {
    std::vector<char> first_key(adjacency_.size(), false);
    first_key[terminals_[0]] = true;
    cache_[0][terminals_[0]][first_key] = {1};
    Vertex nr_vertices = adjacency_.size();
    // omp_set_num_threads(1);
    for(Edge_length length = 0; length < max_length_; length++) {
        for(Vertex start = 0; start < nr_vertices; start++) {
            #pragma omp parallel for default(none) shared(std::cerr) shared(length) shared(start) shared(adjacency_) shared(cache_) shared(distance_) shared(invalid_) shared(thread_local_result_) shared(terminals_)
            for(size_t bucket = 0; bucket < cache_[length][start].bucket_count(); bucket++) {
                Edge_length budget = max_length_ - length;
                size_t thread_id = omp_get_thread_num();
                for(auto task_it = cache_[length][start].begin(bucket); task_it != cache_[length][start].end(bucket); ++task_it) {
                    auto const& visited = task_it->first;
                    auto const& result = task_it->second;
                    for(auto v : neighbors(start)) {
                        if(visited[v]) {
                            continue;
                        }
                        edges_[thread_id]++;
                        // update the partial result
                        std::vector<Edge_weight> new_result(result.size() + adjacency_[start][v].back().first, 0);
                        for(Edge_length res_length = 0; res_length < result.size(); res_length++) {
                            for(auto [e_length, e_weight] : adjacency_[start][v]) {
                                new_result[res_length + e_length] += result[res_length] * e_weight;
                            }
                        }
                        new_result.resize(std::min(Edge_length(new_result.size()), Edge_length(max_length_ + 1)));
                        if(v == terminals_[1]) {
                            for(Edge_length res_length = 0; res_length < new_result.size(); res_length++) {
                                thread_local_result_[thread_id][res_length] += new_result[res_length];
                            }
                            continue;
                        }
                        Edge_length v_budget = budget - adjacency_[start][v][0].first;
                        std::vector<Edge_length> distance_to_goal(adjacency_.size(), invalid_);
                        // make sure we disable paths going through v
                        distance_to_goal[v] = 0;
                        pruning_dijkstra(terminals_[1], v, distance_to_goal, visited, v_budget);
                        // unset the hack value for v
                        distance_to_goal[v] = invalid_;
                        std::vector<Vertex> poss_non_dag, poss_dag;
                        for(auto w : neighbors(v)) {
                            if(v_budget > distance_to_goal[w] + adjacency_[v][w].begin()->first) {
                                poss_non_dag.push_back(w);
                            }
                            if(v_budget == distance_to_goal[w] + adjacency_[v][w].begin()->first) {
                                poss_dag.push_back(w);
                            }
                        }
                        // // there is no edge
                        // if(poss_non_dag.size() + poss_dag.size() == 0) {
                        //     visited_[start] = false;
                        //     for(Vertex excluded : exclude_[start]) {
                        //         visited_[excluded] = false;
                        //     }
                        //     return {};
                        // }
                        // // there are only dag edges
                        // if(poss_non_dag.size() == 0 && enable_dag_) {
                        //     // dont cache and do everything here
                        //     std::vector<Edge_weight> ret(budget + 1, 0);
                        //     for(auto v : poss_dag) {
                        //         dags++;
                        //         std::vector<Edge_weight> tmp = dag_search(v, budget - adjacency_[start][v].begin()->first);
                        //         for(size_t i = 0; i < tmp.size(); i++) {
                        //             for(auto &[length, weight] : adjacency_[start][v]) {
                        //                 if(length + i > budget) {
                        //                     break;
                        //                 }
                        //                 ret[length + i] += weight*tmp[i];
                        //             }
                        //         }
                        //     }
                        //     visited_[start] = false;
                        //     for(Vertex excluded : exclude_[start]) {
                        //         visited_[excluded] = false;
                        //     }
                        //     return ret;
                        // }
                        // // there is only one edge
                        // if(poss_dag.size() + poss_non_dag.size() == 1) {
                        //     propagations++;
                        //     std::vector<Edge_weight> tmp;
                        //     Vertex v;
                        //     if(poss_dag.size() > 0) {
                        //         v = poss_dag[0];
                        //     } else {
                        //         v = poss_non_dag[0];
                        //     }
                        //     tmp = search(v, budget - adjacency_[start][v].begin()->first);
                        //     std::vector<Edge_weight> ret(budget + 1, 0);
                        //     for(size_t i = 0; i < tmp.size(); i++) {
                        //         for(auto &[length, weight] : adjacency_[start][v]) {
                        //             if(length + i > budget) {
                        //                 break;
                        //             }
                        //             ret[length + i] += weight*tmp[i];
                        //         }
                        //     }
                        //     visited_[start] = false;
                        //     for(Vertex excluded : exclude_[start]) {
                        //         visited_[excluded] = false;
                        //     }
                        //     return ret;
                        // }
                        std::vector<char> new_visited(adjacency_.size(), true);
                        // for(size_t i = 0; i < adjacency_.size(); i++) {
                        //     new_visited[i] = distance_to_goal[i] == invalid_;
                        // }
                        prune_articulation(v, new_visited, distance_to_goal);
                        new_visited[v] = true;
                        #pragma omp critical
                        {
                            auto ins = cache_[length + adjacency_[start][v][0].first][v].insert(
                                std::make_pair(new_visited, new_result)
                            );
                            if(!ins.second) {
                                pos_hits_[thread_id]++;
                                // there is already an element with that key
                                // instead increase the partial result for that key
                                if(ins.first->second.size() < new_result.size()) {
                                    ins.first->second.resize(new_result.size());
                                }
                                for(Edge_length res_length = 0; res_length < new_result.size(); res_length++) {
                                    ins.first->second[res_length] += new_result[res_length];
                                }
                            } else {
                                neg_hits_[thread_id]++;
                            }
                        }
                    }
                }
            }
        }
        cache_[length].resize(0);
        // std::cerr << length << std::endl;
        // print_stats();
    }
    for(Edge_length length = 0; length <= max_length_; length++) {
        for(size_t id = 0; id < nthreads_; id++) {
            result_[length] += thread_local_result_[id][length];
        }
    }
    return result_;
}


void ParallelSearch::prune_articulation(Vertex start, std::vector<char>& visited, std::vector<Edge_length>& distance) {
    std::vector<Edge_length> ap_disc(adjacency_.size(), 0);
    std::vector<Edge_length> ap_low(adjacency_.size(), 0);
    int time = 0;
    ap_util(start, visited, ap_disc, ap_low, time, -1, start, distance);
}

bool ParallelSearch::ap_util(Vertex u, std::vector<char>& unvisited, std::vector<Vertex>& disc, std::vector<Vertex>& low, int& time, int parent, Vertex start, std::vector<Edge_length>& distance) {
    // We do not need to check whether the root is an articulation point
    // since if it is, then the goal can only be in one of the induced components
    // for all the other components we cannot enter them since we cannot reach the goal from them
    // but this means that we have already pruned them using dijkstra
    
    // Mark the current node as visited
    unvisited[u] = false;

    bool found_elsewhere = false;
 
    // Initialize discovery time and low value
    disc[u] = low[u] = ++time;
 
    // Go through all vertices adjacent to this
    for (auto v : neighbors(u)) {
        // If v is not visited yet, then make it a child of u
        // in DFS tree and recur for it            
        if(v != start && distance[v] == invalid_) {
            continue;
        }
        if (unvisited[v]) {
            bool found_here = ap_util(v, unvisited, disc, low, time, u, start, distance);
            found_elsewhere |= found_here;
 
            // Check if the subtree rooted with v has
            // a connection to one of the ancestors of u
            low[u] = std::min(low[u], low[v]);
 
            // If u is not root and low value of one of
            // its child is more than discovery value of u.
            if (parent != -1 && low[v] >= disc[u]) {
                // AP
                if(!found_here) {
                    unvisited[u] = true;
                    prune_util(v, unvisited);
                    unvisited[u] = false;
                }
            }
        } else if (v != parent) {
            low[u] = std::min(low[u], disc[v]);
        }
    }
    return found_elsewhere || u == terminals_[1];
}

void ParallelSearch::prune_util(Vertex u, std::vector<char>& unvisited) {
    for (auto v : neighbors(u)) {
        if (!unvisited[v]) {
            unvisited[v] = true;
            prune_util(v, unvisited);
        }
    }
}

void ParallelSearch::dijkstra(Vertex start, std::vector<Edge_length>& distance) {
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
            Edge_length min_cost = adjacency_[cur_vertex][w].begin()->first;
            if(cur_cost + min_cost < distance[w]) {
                distance[w] = min_cost + cur_cost;
                queue.push(std::make_pair(cur_cost + min_cost, w));
            }
        }
    }
}
void ParallelSearch::pruning_dijkstra(Vertex start, Vertex prune, std::vector<Edge_length>& distance, std::vector<char> const& visited, Edge_length budget) {
    std::deque<Vertex> queue;
    queue.push_back(start);
    distance[start] = 0;
    while(!queue.empty()) {
        auto cur_vertex = queue.front();
        auto cur_cost = distance[cur_vertex];
        queue.pop_front();
        for(auto &w : neighbors(cur_vertex)) {
            if(cur_cost + 1 >= distance[w] || visited[w]) {
                continue;
            }
            Edge_length min_cost = adjacency_[cur_vertex][w].begin()->first;
            if(cur_cost + min_cost < distance[w] && cur_cost + min_cost + distance_[prune][w] <= budget) {
                distance[w] = min_cost + cur_cost;
                queue.push_back(w);
            }
        }
    }
}

void ParallelSearch::print_stats() {
    size_t pos_hits = 0, neg_hits = 0;
    for(size_t i = 0; i < OMP_NUM_THREADS; i++) {
        pos_hits += pos_hits_[i];
        neg_hits += neg_hits_[i];
    }
    size_t dags = 0, splits = 0;
    for(size_t i = 0; i < OMP_NUM_THREADS; i++) {
        dags += dags_[i];
        // splits += neg_hits_[i];
    }
    size_t edges = 0, propagations = 0;
    for(size_t i = 0; i < OMP_NUM_THREADS; i++) {
        edges += edges_[i];
        propagations += propagations_[i];
    }
    std::cerr << "Cache hit rate: " << 100*pos_hits/(double)(pos_hits + neg_hits) << "% (" << pos_hits << "/" << pos_hits + neg_hits << ")" << std::endl;
    std::cerr << "#DAG searches: " << dags << " #Splits: " << splits << std::endl;
    std::cerr << "#Edges: " << edges << " #Propagations: " << propagations << std::endl;
}