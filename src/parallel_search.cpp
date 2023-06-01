#include "parallel_search.h"
#include <queue>
#include <limits>

namespace fpc {

ParallelSearch::ParallelSearch(sparsegraph input, Edge_length max_length) :  
                                nthreads_(OMP_NUM_THREADS),
                                max_length_(max_length),
                                invalid_(std::numeric_limits<Edge_length>::max() - max_length_ - 1),
                                initial_(input),
                                cache_( 
                                    input.nv + 1,
                                    std::unordered_map<PCacheKey, std::vector<Edge_weight>, sg_hash, sg_equal>()
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
                                dags_(OMP_NUM_THREADS, 0) {
    omp_set_num_threads(nthreads_);
}

std::vector<Edge_weight> ParallelSearch::search() {
    cache_[initial_.nv][initial_] = {1};
    Vertex nr_vertices = initial_.nv;
    for(Vertex remaining_size = nr_vertices + 1; remaining_size-- > 0; ) {
        //#pragma omp parallel for default(none) shared(max_length_) shared(remaining_size) shared(cache_) shared(invalid_) shared(thread_local_result_) shared(dispatch_sparse)
        for(size_t bucket = 0; bucket < cache_[remaining_size].bucket_count(); bucket++) {
            size_t thread_id = 0; //omp_get_thread_num();
            for(auto task_it = cache_[remaining_size].begin(bucket); task_it != cache_[remaining_size].end(bucket); ++task_it) {
                auto const& old_sg = task_it->first;
                auto const& result = task_it->second;
                Edge_length budget = max_length_;
                while(result[max_length_ - budget] == 0) {
                    budget--;
                }
                std::vector<Vertex> poss;
                for(int i = 0; i < old_sg.d[0]; i++) {
                    Vertex v = old_sg.e[i];
                    if(v == 1) {
                        edges_[thread_id]++;
                        // update the partial result
                        for(Edge_length res_length = 0; res_length < result.size(); res_length++) {
                            thread_local_result_[thread_id][res_length + 1] += result[res_length];
                        }
                        continue;
                    }
                    // if(budget > old_distance_to_goal[v] + adjacency_[start][v][0].first) {
                    edges_[thread_id]++;
                    poss.push_back(v);
                    // }
                    // if(budget == old_distance_to_goal[v] + adjacency_[start][v][0].first) {
                    //     edges_[thread_id]++;
                    //     dags_[thread_id]++;
                    //     Edge_weight result_until = result[max_length_ - budget] * adjacency_[start][v][0].second;
                    //     Edge_weight remaining_result = dag_search(v, old_distance_to_goal);
                    //     thread_local_result_[thread_id][max_length_] += result_until*remaining_result;
                    //     continue;
                    // }
                }
                for(Vertex v : poss) {
                    // update the partial result
                    std::vector<Edge_weight> new_result(result.size() + 1);
                    for(Edge_length res_length = 0; res_length < result.size(); res_length++) {
                            new_result[res_length + 1] += result[res_length];
                    }
                    Edge_length v_budget = budget - 1;
                    Edge_length extra_length = 1;
                    std::vector<Edge_length> distance_to_goal(old_sg.nv, invalid_);
                    // make sure we disable paths going through v
                    distance_to_goal[0] = 0;
                    distance_to_goal[v] = 0;
                    pruning_dijkstra(old_sg, v, distance_to_goal, v_budget);
                    // unset the hack value for v
                    distance_to_goal[0] = invalid_;
                    distance_to_goal[v] = invalid_;
                    std::vector<Vertex> poss_non_dag, poss_dag;
                    for(int i = 0; i < old_sg.d[v]; i++) {
                        Vertex w = old_sg.e[old_sg.v[v] + i];
                        if(w == 1) {
                            continue;
                        }
                        if(v_budget == distance_to_goal[w] + 1) {
                            poss_dag.push_back(w);
                        } else if(v_budget > distance_to_goal[w] + 1) {
                            poss_non_dag.push_back(w);
                        }
                    }
                    // there is only one non_dag edge
                    std::vector<Vertex> extra = {0,v};
                    Vertex last = v;
                    while(poss_non_dag.size() == 1) {
                        edges_[thread_id]++;
                        propagations_[thread_id]++;
                        // first check if the goal is a neighbor
                        for(int i = 0; i < old_sg.d[last]; i++) {
                            Vertex w = old_sg.e[old_sg.v[last] + i];
                            if(w == 1) {
                                edges_[thread_id]++;
                                // update the partial result
                                for(Edge_length res_length = 0; res_length < new_result.size(); res_length++) {
                                    thread_local_result_[thread_id][res_length + 1] += new_result[res_length];
                                }
                            }
                        }

                        for(Vertex dag_v : poss_dag) {
                            edges_[thread_id]++;
                            dags_[thread_id]++;
                            Edge_weight result_until = new_result[max_length_ - v_budget];
                            Edge_weight remaining_result = dag_search(old_sg, dag_v, distance_to_goal);
                            thread_local_result_[thread_id][max_length_] += result_until*remaining_result;
                        }
                        Vertex cur = poss_non_dag[0];
                        extra.push_back(cur);
                        new_result.resize(new_result.size() + 1);
                        for(Edge_length res_length = new_result.size() - 1; res_length >= 1; res_length--) {
                            new_result[res_length] = new_result[res_length - 1];
                        }
                        new_result[0] = 0;
                        v_budget -= 1;
                        extra_length += 1;
                        std::fill(distance_to_goal.begin(), distance_to_goal.end(), invalid_);
                        for(Vertex extra_v : extra) {
                            distance_to_goal[extra_v] = 0;
                        }
                        pruning_dijkstra(old_sg, cur, distance_to_goal, v_budget);
                        for(Vertex extra_v : extra) {
                            distance_to_goal[extra_v] = invalid_;
                        }
                        poss_dag.clear();
                        poss_non_dag.clear();
                        for(int i = 0; i < old_sg.d[cur]; i++) {
                            Vertex w = old_sg.e[old_sg.v[cur] + i];
                            if(w == 1) {
                                continue;
                            }
                            if(v_budget == distance_to_goal[w] + 1) {
                                poss_dag.push_back(w);
                            } else if(v_budget > distance_to_goal[w] + 1) {
                                poss_non_dag.push_back(w);
                            }
                        }
                        last = cur;
                    }
                    if(poss_non_dag.size() == 0) {
                        // first check if the goal is a neighbor
                        for(int i = 0; i < old_sg.d[last]; i++) {
                            Vertex w = old_sg.e[old_sg.v[last] + i];
                            if(w == 1) {
                                edges_[thread_id]++;
                                // update the partial result
                                for(Edge_length res_length = 0; res_length < new_result.size(); res_length++) {
                                    thread_local_result_[thread_id][res_length + 1] += new_result[res_length];
                                }
                            }
                        }
                        for(Vertex dag_v : poss_dag) {
                            edges_[thread_id]++;
                            dags_[thread_id]++;
                            Edge_weight result_until = new_result[max_length_ - v_budget];
                            Edge_weight remaining_result = dag_search(old_sg, dag_v, distance_to_goal);
                            thread_local_result_[thread_id][max_length_] += result_until*remaining_result;
                        }
                        continue;
                    }
                    prune_articulation(old_sg, last, distance_to_goal);
                    assert(last != 1);
                    assert(last < old_sg.nv);
                    assert(new_result.size() <= max_length_ + 1);
                    // construct the new canonical graph
                    DEFAULTOPTIONS_SPARSEGRAPH(options);
                    options.getcanon = true;
                    options.defaultptn = false;
                    SG_DECL(sg);
                    sg.v = (size_t *)malloc(sizeof(size_t)*old_sg.nv + sizeof(int)*(old_sg.nv + old_sg.nde) + 100);
                    sg.d = (int *)(sg.v + old_sg.nv);
                    sg.e = sg.d + old_sg.nv;
                    size_t nr_edges = 0;
                    int skipped = -1;
                    std::vector<int> mapping(old_sg.nv, -1);
                    mapping[last] = 0;
                    mapping[1] = 1;
                    int ctr = 2;
                    for(size_t i = 0; i < old_sg.nv; i++) {
                        if(distance_to_goal[i] != invalid_ && mapping[i] == -1) {
                            mapping[i] = ctr++;
                        }
                    }
                    sg.v[0] = 0;
                    for(size_t j = 0; j < old_sg.d[last]; j++) {
                        assert(old_sg.v[last] + j < old_sg.nde);
                        Vertex w = old_sg.e[old_sg.v[last] + j];
                        assert(w < old_sg.nv);
                        if(distance_to_goal[w] != invalid_) {
                            sg.e[nr_edges++] = mapping[w];
                        }
                    }
                    sg.d[0] = nr_edges;
                    for(size_t i = 0; i < old_sg.nv; i++) {
                        if(distance_to_goal[i] == invalid_ ) {
                            skipped++;
                            continue;
                        }
                        assert(i != last);
                        sg.v[mapping[i]] = nr_edges;
                        for(size_t j = 0; j < old_sg.d[i]; j++) {
                            assert(old_sg.v[i] + j < old_sg.nde);
                            Vertex w = old_sg.e[old_sg.v[i] + j];
                            assert(w < old_sg.nv);
                            if(distance_to_goal[w] != invalid_) {
                                sg.e[nr_edges++] = mapping[w];
                            }
                            if(w == last) {
                                sg.e[nr_edges++] = 0;
                            }
                        }
                        sg.d[mapping[i]] = nr_edges - sg.v[mapping[i]];
                    }
                    sg.nv = old_sg.nv - skipped;
                    sg.nde = nr_edges;
                    assert(nr_edges < old_sg.nde);
                    int *lab = (int *)malloc(sg.nv*sizeof(int));
                    int *ptn = (int *)malloc(sg.nv*sizeof(int));
                    int *orbits = (int *)malloc(sg.nv*sizeof(int));
                    ptn[0] = 0;
                    ptn[1] = 0;
                    lab[0] = 0;
                    lab[1] = 1;
                    for(size_t i = 2; i < sg.nv; i++) {
                        ptn[i] = 1;
                        lab[i] = i;
                    }
                    statsblk stats;
                    SG_DECL(canon_sg);
                    canon_sg.v = (size_t *)malloc(sizeof(size_t)*sg.nv + sizeof(int)*(sg.nv + sg.nde) + 100);
                    canon_sg.d = (int *)(canon_sg.v + sg.nv);
                    canon_sg.e = canon_sg.d + sg.nv;
                    canon_sg.nv = sg.nv;
                    canon_sg.nde = sg.nde;
                    sparsenauty(&sg,lab,ptn,orbits,&options,&stats,&canon_sg);
                    sortlists_sg(&canon_sg);
                    free(lab);
                    free(ptn);
                    free(orbits);
                    #pragma omp critical
                    {
                        auto ins = cache_[canon_sg.nv].insert(
                            std::make_pair(canon_sg, new_result)
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
                            free(canon_sg.v);
                        } else {
                            neg_hits_[thread_id]++;
                        }
                    }
                }
                free(old_sg.v);
            }
        }
        cache_[remaining_size].clear();
        // std::cerr << remaining_size << std::endl;
        // print_stats();
    }
    for(Edge_length length = 0; length <= max_length_; length++) {
        for(size_t id = 0; id < nthreads_; id++) {
            result_[length] += thread_local_result_[id][length];
        }
    }
    return result_;
}

Edge_weight ParallelSearch::dag_search(sparsegraph const& sg, Vertex start, std::vector<Edge_length> const& distance_to_goal) {
    std::deque<std::pair<Edge_length, Vertex>> biggest_first_queue;
    std::vector<Edge_weight> dp(sg.nv, 0);
    dp[start] = 1;
    biggest_first_queue.push_back(std::make_pair(distance_to_goal[start], start));
    while(biggest_first_queue.front().first != 0) {
        auto [cur_cost, cur_vertex] = biggest_first_queue.front();
        biggest_first_queue.pop_front();
        for(int i = 0; i < sg.d[cur_vertex]; i++) {
            Vertex w = sg.e[sg.v[cur_vertex] + i];
            if(distance_to_goal[w] == invalid_ || dp[cur_vertex] == 0) {
                continue;
            }
            if(distance_to_goal[w] + 1 == distance_to_goal[cur_vertex]) {
                dp[w] += dp[cur_vertex];
                biggest_first_queue.push_back(std::make_pair(distance_to_goal[w], w));
            }
        }
        dp[cur_vertex] = 0;
    }
    assert(dp[1] >= 1);
    return dp[1];
}

void ParallelSearch::prune_articulation(sparsegraph const& sg, Vertex start, std::vector<Edge_length>& distance) {
    std::vector<Vertex> ap_disc(sg.nv, 0);
    std::vector<Vertex> ap_low(sg.nv, 0);
    std::vector<char> ap_visited(sg.nv, false);
    int time = 0;
    ap_util(sg, start, ap_visited, ap_disc, ap_low, time, -1, start, distance);
}

bool ParallelSearch::ap_util(sparsegraph const& sg, Vertex u, std::vector<char>& visited, std::vector<Vertex>& disc, std::vector<Vertex>& low, int& time, int parent, Vertex start, std::vector<Edge_length>& distance) {
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
    for(int i = 0; i < sg.d[u]; i++) {
        Vertex v = sg.e[sg.v[u] + i];
        // If v is not visited yet, then make it a child of u
        // in DFS tree and recur for it            
        if(v != start && distance[v] == invalid_) {
            continue;
        }
        if (!visited[v]) {
            bool found_here = ap_util(sg, v, visited, disc, low, time, u, start, distance);
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
                    prune_util(sg, v, distance);
                    distance[u] = tmp;
                }
            }
        } else if (v != parent) {
            low[u] = std::min(low[u], disc[v]);
        }
    }
    return found_elsewhere || u == 1;
}

void ParallelSearch::prune_util(sparsegraph const& sg, Vertex u, std::vector<Edge_length>& distance) {
    for(int i = 0; i < sg.d[u]; i++) {
        Vertex v = sg.e[sg.v[u] + i];
        if (distance[v] != invalid_) {
            distance[v] = invalid_;
            prune_util(sg, v, distance);
        }
    }
}

void ParallelSearch::pruning_dijkstra(sparsegraph const& sg, Vertex prune, std::vector<Edge_length>& distance, Edge_length budget) {
    std::deque<Vertex> queue;
    queue.push_back(1);
    distance[1] = 0;
    while(!queue.empty()) {
        auto cur_vertex = queue.front();
        auto cur_cost = distance[cur_vertex];
        queue.pop_front();
        assert(cur_vertex < sg.nv);
        for(int i = 0; i < sg.d[cur_vertex]; i++) {
            assert(sg.v[cur_vertex] + i < sg.nde);
            Vertex w = sg.e[sg.v[cur_vertex] + i];
            if(w >= sg.nv) {
                std::cerr << w << " " << sg.nv << std::endl;
            }
            assert(w < sg.nv);
            if(cur_cost + 1 >= distance[w] || cur_cost + 1 > budget) {
                continue;
            }
            distance[w] = 1 + cur_cost;
            if(cur_cost + 1 < budget) {
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

} // namespace fpc