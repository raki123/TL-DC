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
                                    adjacency_.size() + 1,
                                    std::vector<std::unordered_map<PCacheKey, std::vector<Edge_weight>, pvector_hash>>(
                                        adjacency_.size(), 
                                        std::unordered_map<PCacheKey, std::vector<Edge_weight>, pvector_hash>()
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
    omp_set_num_threads(nthreads_);
}

std::vector<Edge_weight> ParallelSearch::search() {
    std::vector<Edge_length> first_key(adjacency_.size(), invalid_);
    dijkstra(terminals_[1], first_key);
    first_key[terminals_[0]] = invalid_;
    cache_[adjacency_.size() - 1][terminals_[0]][first_key] = {1};
    Vertex nr_vertices = adjacency_.size();
    // omp_set_num_threads(1);
    for(Vertex remaining_size = adjacency_.size(); remaining_size-- > 0; ) {
        for(Vertex start = 0; start < nr_vertices; start++) {
            #pragma omp parallel for default(none) shared(max_length_) shared(remaining_size) shared(start) shared(adjacency_) shared(cache_) shared(distance_) shared(invalid_) shared(thread_local_result_) shared(terminals_)
            for(size_t bucket = 0; bucket < cache_[remaining_size][start].bucket_count(); bucket++) {

                size_t thread_id = omp_get_thread_num();
                for(auto task_it = cache_[remaining_size][start].begin(bucket); task_it != cache_[remaining_size][start].end(bucket); ++task_it) {
                    auto const& old_distance_to_goal = task_it->first;
                    auto const& result = task_it->second;
                    Edge_length budget = max_length_;
                    while(result[max_length_ - budget] == 0) {
                        budget--;
                    }
                    std::vector<Vertex> poss;
                    for(auto v : neighbors(start)) {
                        if(v == terminals_[1]) {
                            edges_[thread_id]++;
                            // update the partial result
                            for(Edge_length res_length = 0; res_length < result.size(); res_length++) {
                                for(auto [e_length, e_weight] : adjacency_[start][v]) {
                                    if(res_length + e_length > max_length_) {
                                        break;
                                    }
                                    thread_local_result_[thread_id][res_length + e_length] += result[res_length] * e_weight;
                                }
                            }
                            continue;
                        }
                        if(budget > old_distance_to_goal[v] + adjacency_[start][v][0].first) {
                            edges_[thread_id]++;
                            poss.push_back(v);
                        }
                        if(budget == old_distance_to_goal[v] + adjacency_[start][v][0].first) {
                            edges_[thread_id]++;
                            dags_[thread_id]++;
                            Edge_weight result_until = result[max_length_ - budget] * adjacency_[start][v][0].second;
                            Edge_weight remaining_result = dag_search(v, old_distance_to_goal);
                            thread_local_result_[thread_id][max_length_] += result_until*remaining_result;
                            continue;
                        }
                    }
                    for(auto v : poss) {
                        // update the partial result
                        std::vector<Edge_weight> new_result(result.size() + adjacency_[start][v].back().first);
                        for(Edge_length res_length = 0; res_length < result.size(); res_length++) {
                            for(auto [e_length, e_weight] : adjacency_[start][v]) {
                                new_result[res_length + e_length] += result[res_length] * e_weight;
                            }
                        }
                        new_result.resize(std::min(Edge_length(new_result.size()), Edge_length(max_length_)));
                        Edge_length v_budget = budget - adjacency_[start][v][0].first;
                        Edge_length extra_length = adjacency_[start][v][0].first;
                        std::vector<Edge_length> distance_to_goal(adjacency_.size(), invalid_);
                        // make sure we disable paths going through v
                        distance_to_goal[v] = 0;
                        pruning_dijkstra(terminals_[1], v, distance_to_goal, old_distance_to_goal, v_budget);
                        // unset the hack value for v
                        distance_to_goal[v] = invalid_;
                        std::vector<Vertex> poss_non_dag, poss_dag;
                        for(auto w : neighbors(v)) {
                            if(w == terminals_[1]) {
                                continue;
                            }
                            if(v_budget == distance_to_goal[w] + adjacency_[v][w][0].first) {
                                poss_dag.push_back(w);
                            } else if(v_budget > distance_to_goal[w] + adjacency_[v][w][0].first) {
                                poss_non_dag.push_back(w);
                            }
                        }
                        // there is only one non_dag edge
                        std::vector<Vertex> extra = {v};
                        Vertex last = v;
                        while(poss_non_dag.size() == 1) {
                            edges_[thread_id]++;
                            propagations_[thread_id]++;
                            // first check if the goal is a neighbor
                            if(adjacency_[last][terminals_[1]].size() > 0) {
                                edges_[thread_id]++;
                                // update the partial result
                                for(Edge_length res_length = 0; res_length < new_result.size(); res_length++) {
                                    for(auto [e_length, e_weight] : adjacency_[last][terminals_[1]]) {
                                        if(res_length + e_length > max_length_) {
                                            break;
                                        }
                                        thread_local_result_[thread_id][res_length + e_length] += new_result[res_length] * e_weight;
                                    }
                                }
                            }
                            for(Vertex dag_v : poss_dag) {
                                edges_[thread_id]++;
                                dags_[thread_id]++;
                                Edge_weight result_until = new_result[max_length_ - v_budget] * adjacency_[last][dag_v][0].second;
                                Edge_weight remaining_result = dag_search(dag_v, distance_to_goal);
                                thread_local_result_[thread_id][max_length_] += result_until*remaining_result;
                            }
                            Vertex cur = poss_non_dag[0];
                            extra.push_back(cur);
                            new_result.resize(new_result.size() + adjacency_[last][cur].back().first);
                            for(Edge_length res_length = new_result.size() - 1; res_length >= adjacency_[last][cur][0].first; res_length--) {
                                new_result[res_length] = 0;
                                for(auto [e_length, e_weight] : adjacency_[last][cur]) {
                                    if(e_length > res_length) {
                                        break;
                                    }
                                    new_result[res_length] += new_result[res_length - e_length] * e_weight;
                                }
                            }
                            for(Edge_length res_length = 0 ; res_length < adjacency_[last][cur][0].first; res_length++) {
                                new_result[res_length] = 0;
                            }
                            new_result.resize(std::min(Edge_length(new_result.size()), Edge_length(max_length_)));
                            v_budget -= adjacency_[last][cur][0].first;
                            extra_length += adjacency_[last][cur][0].first;
                            std::fill(distance_to_goal.begin(), distance_to_goal.end(), invalid_);
                            for(Vertex extra_v : extra) {
                                distance_to_goal[extra_v] = 0;
                            }
                            pruning_dijkstra(terminals_[1], cur, distance_to_goal, old_distance_to_goal, v_budget);
                            for(Vertex extra_v : extra) {
                                distance_to_goal[extra_v] = invalid_;
                            }
                            poss_dag.clear();
                            poss_non_dag.clear();
                            for(auto w : neighbors(cur)) {
                                if(w == terminals_[1]) {
                                    continue;
                                }
                                if(v_budget == distance_to_goal[w] + adjacency_[cur][w][0].first) {
                                    poss_dag.push_back(w);
                                } else if(v_budget > distance_to_goal[w] + adjacency_[cur][w][0].first) {
                                    poss_non_dag.push_back(w);
                                }
                            }
                            last = cur;
                        }
                        if(poss_non_dag.size() == 0) {
                            if(adjacency_[last][terminals_[1]].size() > 0) {
                                edges_[thread_id]++;
                                // update the partial result
                                for(Edge_length res_length = 0; res_length < new_result.size(); res_length++) {
                                    for(auto [e_length, e_weight] : adjacency_[last][terminals_[1]]) {
                                        if(res_length + e_length > max_length_) {
                                            break;
                                        }
                                        thread_local_result_[thread_id][res_length + e_length] += new_result[res_length] * e_weight;
                                    }
                                }
                            }
                            for(Vertex dag_v : poss_dag) {
                                edges_[thread_id]++;
                                dags_[thread_id]++;
                                Edge_weight result_until = new_result[max_length_ - v_budget] * adjacency_[last][dag_v][0].second;
                                Edge_weight remaining_result = dag_search(dag_v, distance_to_goal);
                                thread_local_result_[thread_id][max_length_] += result_until*remaining_result;
                            }
                            continue;
                        }
                        prune_articulation(last, distance_to_goal);
                        assert(last != terminals_[1]);
                        assert(last < adjacency_.size());
                        assert(new_result.size() <= max_length_ + 1);
                        Vertex new_remaining_size = 0;
                        for(Vertex i = 0; i < adjacency_.size(); i++) {
                            if(distance_to_goal[i] != invalid_) {
                                new_remaining_size++;
                            }
                        }
                        #pragma omp critical
                        {
                            auto ins = cache_[new_remaining_size][last].insert(
                                std::make_pair(distance_to_goal, new_result)
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
        cache_[remaining_size].resize(0);
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

Edge_weight ParallelSearch::dag_search(Vertex start, std::vector<Edge_length> const& distance_to_goal) {
    std::priority_queue<std::pair<Edge_length, Vertex>> biggest_first_queue;
    std::vector<Edge_weight> dp(adjacency_.size(), 0);
    std::vector<char> in_queue(adjacency_.size(), false);
    dp[start] = 1;
    biggest_first_queue.push(std::make_pair(distance_to_goal[start], start));
    while(!biggest_first_queue.empty()) {
        auto [cur_cost, cur_vertex] = biggest_first_queue.top();
        biggest_first_queue.pop();
        for(auto w : neighbors(cur_vertex)) {
            if(distance_to_goal[w] == invalid_) {
                continue;
            }
            auto edge_cost = adjacency_[w][cur_vertex][0].first;
            // this edge can be part of a shorted path in direction w -> cur_vertex
            if(distance_to_goal[w] == distance_to_goal[cur_vertex] + edge_cost) {
                auto factor = adjacency_[w][cur_vertex][0].second;
                dp[cur_vertex] += factor*dp[w];
            // this edge can be part of a shorted path in direction cur_vertex -> w
            } else if(distance_to_goal[w] + edge_cost == distance_to_goal[cur_vertex]) {
                if(!in_queue[w]) {
                    biggest_first_queue.push(std::make_pair(distance_to_goal[w], w));
                    in_queue[w] = true;
                }
            }
        }
    }
    assert(dp[terminals_[1]] >= 1);
    return dp[terminals_[1]];
}

void ParallelSearch::prune_articulation(Vertex start, std::vector<Edge_length>& distance) {
    std::vector<Vertex> ap_disc(adjacency_.size(), 0);
    std::vector<Vertex> ap_low(adjacency_.size(), 0);
    std::vector<char> ap_visited(adjacency_.size(), false);
    int time = 0;
    ap_util(start, ap_visited, ap_disc, ap_low, time, -1, start, distance);
}

bool ParallelSearch::ap_util(Vertex u, std::vector<char>& visited, std::vector<Vertex>& disc, std::vector<Vertex>& low, int& time, int parent, Vertex start, std::vector<Edge_length>& distance) {
    // We do not need to check whether the root is an articulation point
    // since if it is, then the goal can only be in one of the induced components
    // for all the other components we cannot enter them since we cannot reach the goal from them
    // but this means that we have already pruned them using dijkstra
    
    // Mark the current node as visited
    visited[u] = true;

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
        if (!visited[v]) {
            bool found_here = ap_util(v, visited, disc, low, time, u, start, distance);
            found_elsewhere |= found_here;
 
            // Check if the subtree rooted with v has
            // a connection to one of the ancestors of u
            low[u] = std::min(low[u], low[v]);
 
            // If u is not root and low value of one of
            // its child is more than discovery value of u.
            if (parent != -1 && low[v] >= disc[u]) {
                // AP
                if(!found_here) {
                    auto tmp = distance[u];
                    distance[u] = invalid_;
                    distance[v] = invalid_;
                    prune_util(v, distance);
                    distance[u] = tmp;
                }
            }
        } else if (v != parent) {
            low[u] = std::min(low[u], disc[v]);
        }
    }
    return found_elsewhere || u == terminals_[1];
}

void ParallelSearch::prune_util(Vertex u, std::vector<Edge_length>& distance) {
    for (auto v : neighbors(u)) {
        if (distance[v] != invalid_) {
            distance[v] = invalid_;
            prune_util(v, distance);
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
void ParallelSearch::pruning_dijkstra(Vertex start, Vertex prune, std::vector<Edge_length>& distance, std::vector<Edge_length> const& old_distance, Edge_length budget) {
    std::deque<Vertex> queue;
    queue.push_back(start);
    distance[start] = 0;
    while(!queue.empty()) {
        auto cur_vertex = queue.front();
        auto cur_cost = distance[cur_vertex];
        queue.pop_front();
        for(auto &w : neighbors(cur_vertex)) {
            if(cur_cost + 1 >= distance[w] || old_distance[w] == invalid_) {
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