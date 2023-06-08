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
                                thread_local_sg_(OMP_NUM_THREADS),
                                thread_local_lab_(OMP_NUM_THREADS),
                                thread_local_ptn_(OMP_NUM_THREADS),
                                thread_local_orbits_(OMP_NUM_THREADS),
                                pos_hits_(OMP_NUM_THREADS, 0),
                                neg_hits_(OMP_NUM_THREADS, 0),
                                edges_(OMP_NUM_THREADS, 0),
                                propagations_(OMP_NUM_THREADS, 0),
                                dags_(OMP_NUM_THREADS, 0) {
    omp_set_num_threads(nthreads_);
    for(size_t i = 0; i < OMP_NUM_THREADS; i++) {
        SG_INIT(thread_local_sg_[i]);
        thread_local_sg_[i].v = (edge_t *)malloc(sizeof(edge_t)*(input.nv + input.nde) + sizeof(degree_t)*input.nv);
        thread_local_lab_[i] = (int *)malloc(3*input.nv*sizeof(int));
        thread_local_ptn_[i] = thread_local_lab_[i] + input.nv;
        thread_local_orbits_[i] = thread_local_ptn_[i] + input.nv;
        thread_local_ptn_[i][0] = 0;
        thread_local_ptn_[i][1] = 0;
        thread_local_lab_[i][0] = 0;
        thread_local_lab_[i][1] = 1;
        for(size_t j = 2; j < input.nv; j++) {
            thread_local_ptn_[i][j] = 1;
            thread_local_lab_[i][j] = j;
        }
    }
}

std::vector<Edge_weight> ParallelSearch::search() {
    cache_[initial_.nv][initial_] = {1};
    Vertex nr_vertices = initial_.nv;
    for(Vertex remaining_size = nr_vertices + 1; remaining_size-- > 0; ) {
        #pragma omp parallel for default(none) shared(max_length_) shared(remaining_size) shared(cache_) shared(invalid_) shared(thread_local_result_) shared(thread_local_sg_) shared(thread_local_ptn_) shared(thread_local_lab_) shared(thread_local_orbits_) shared(dispatch_sparse)
        for(size_t bucket = 0; bucket < cache_[remaining_size].bucket_count(); bucket++) {
            size_t thread_id = omp_get_thread_num();
            for(auto task_it = cache_[remaining_size].begin(bucket); task_it != cache_[remaining_size].end(bucket); ++task_it) {
                auto const& old_sg = task_it->first;
                auto const& result = task_it->second;
                Edge_length budget = max_length_;
                while(result[max_length_ - budget] == 0) {
                    budget--;
                }
                for(int o = 0; o < old_sg.d[0]; o++) {
                    Vertex last = 0;
                    assert(old_sg.e[o] != 1);
                    edges_[thread_id]++;
                    std::vector<Vertex> poss_non_dag = {Vertex(old_sg.e[o])};
                    std::vector<Vertex> extra = {0};
                    Edge_length v_budget = budget;
                    std::vector<Edge_length> distance_to_goal(old_sg.nv, invalid_);
                    std::vector<Edge_weight> new_result = result;
                    // there is only one non_dag edge
                    while(poss_non_dag.size() == 1) {
                        edges_[thread_id]++;
                        propagations_[thread_id]++;
                        Vertex cur = poss_non_dag[0];
                        extra.push_back(cur);
                        new_result.resize(std::min(new_result.size() + 1, size_t(max_length_)));
                        for(Edge_length res_length = new_result.size() - 1; res_length >= 1; res_length--) {
                            new_result[res_length] = new_result[res_length - 1];
                        }
                        new_result[0] = 0;
                        v_budget -= 1;
                        std::fill(distance_to_goal.begin(), distance_to_goal.end(), invalid_);
                        for(Vertex extra_v : extra) {
                            distance_to_goal[extra_v] = 0;
                        }
                        pruning_dijkstra(old_sg, distance_to_goal, v_budget);
                        for(Vertex extra_v : extra) {
                            distance_to_goal[extra_v] = invalid_;
                        }
                        poss_non_dag.clear();
                        for(int i = 0; i < old_sg.d[cur]; i++) {
                            Vertex w = old_sg.e[old_sg.v[cur] + i];
                            if(w == 1) {
                                edges_[thread_id]++;
                                // update the partial result
                                for(Edge_length res_length = 0; res_length < new_result.size(); res_length++) {
                                    thread_local_result_[thread_id][res_length + 1] += new_result[res_length];
                                }
                            } else if(v_budget == distance_to_goal[w] + 1) {
                                extra.push_back(w);
                                edges_[thread_id]++;
                                dags_[thread_id]++;
                                Edge_weight result_until = new_result[max_length_ - v_budget];
                                Edge_weight remaining_result = dag_search(old_sg, w, distance_to_goal);
                                thread_local_result_[thread_id][max_length_] += result_until*remaining_result;
                            } else if(v_budget > distance_to_goal[w] + 1) {
                                poss_non_dag.push_back(w);
                            }
                        }
                        last = cur;
                    }
                    if(poss_non_dag.size() == 0) {
                        continue;
                    }
                    std::vector<Edge_length> distance_to_start(old_sg.nv, invalid_);
                    for(Vertex extra_v : extra) {
                        distance_to_start[extra_v] = 0;
                    }
                    reverse_pruning_dijkstra(old_sg, last, distance_to_start, distance_to_goal, v_budget);
                    for(Vertex extra_v : extra) {
                        distance_to_start[extra_v] = invalid_;
                    }
                    prune_articulation(old_sg, last, distance_to_start);
                    assert(last != 1);
                    assert(last < old_sg.nv);
                    // construct the new canonical graph
                    DEFAULTOPTIONS_SPARSEGRAPH(options);
                    options.getcanon = true;
                    options.defaultptn = false;
                    sparsegraph sg = thread_local_sg_[thread_id];
                    sg.d = (degree_t *)(sg.v + old_sg.nv);
                    sg.e = (edge_t *)(sg.d + old_sg.nv);
                    sg.vlen = old_sg.nv;
                    sg.dlen = old_sg.nv;
                    sg.elen = old_sg.nde;
                    size_t nr_edges = 0;
                    int skipped = 0;
                    std::vector<int> mapping(old_sg.nv, -1);
                    mapping[last] = 0;
                    mapping[1] = 1;
                    int ctr = 2;
                    for(size_t i = 0; i < old_sg.nv; i++) {
                        if(distance_to_start[i] != invalid_ && mapping[i] == -1) {
                            mapping[i] = ctr++;
                        }
                    }
                    sg.v[0] = 0;
                    for(size_t j = 0; j < old_sg.d[last]; j++) {
                        assert(old_sg.v[last] + j < old_sg.nde);
                        Vertex w = old_sg.e[old_sg.v[last] + j];
                        assert(w < old_sg.nv);
                        if(distance_to_start[w] != invalid_ && w != 1) {
                            sg.e[nr_edges++] = mapping[w];
                        }
                    }
                    sg.d[0] = nr_edges;
                    assert(distance_to_start[1] != invalid_);
                    sg.v[1] = nr_edges;
                    for(size_t j = 0; j < old_sg.d[1]; j++) {
                        assert(old_sg.v[1] + j < old_sg.nde);
                        Vertex w = old_sg.e[old_sg.v[1] + j];
                        assert(w < old_sg.nv);
                        if(distance_to_start[w] != invalid_ && w != last) {
                            sg.e[nr_edges++] = mapping[w];
                        }
                    }
                    sg.d[1] = nr_edges - sg.v[1];
                    for(size_t i = 2; i < old_sg.nv; i++) {
                        if(distance_to_start[i] == invalid_ ) {
                            skipped++;
                            continue;
                        }
                        assert(i != last);
                        sg.v[mapping[i]] = nr_edges;
                        for(size_t j = 0; j < old_sg.d[i]; j++) {
                            assert(old_sg.v[i] + j < old_sg.nde);
                            Vertex w = old_sg.e[old_sg.v[i] + j];
                            assert(w < old_sg.nv);
                            if(distance_to_start[w] != invalid_) {
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
                    statsblk stats;
                    SG_DECL(canon_sg);
                    canon_sg.v = (edge_t *)malloc(sizeof(edge_t)*(sg.nv + sg.nde) + sizeof(degree_t)*sg.nv);
                    canon_sg.d = (degree_t *)(canon_sg.v + sg.nv);
                    canon_sg.e = (edge_t *)(canon_sg.d + sg.nv);
                    canon_sg.nv = sg.nv;
                    canon_sg.nde = sg.nde;
                    canon_sg.vlen = sg.nv;
                    canon_sg.dlen = sg.nv;
                    canon_sg.elen = sg.nde;
                    auto old_v = canon_sg.v;
                    auto old_d = canon_sg.d;
                    auto old_e = canon_sg.e;
                    int *lab = thread_local_lab_[thread_id];
                    int *ptn = thread_local_ptn_[thread_id];
                    int *orbits = thread_local_orbits_[thread_id];
                    ptn[0] = 0;
                    ptn[1] = 0;
                    std::fill(ptn + 2, ptn + sg.nv, 1);
                    for(size_t j = 0; j < sg.nv; j++) {
                        lab[j] = j;
                    }
                    sparsenauty(&sg,lab,ptn,orbits,&options,&stats,&canon_sg);
                    // We could apply some additional check to see if every neighbor is in the same orbit to save some time
                    // int least = orbits[sg.e[0]];
                    // for(size_t i = 1; i < sg.d[0]; i++) {
                    //     if(least != orbits[sg.e[i]]) {
                    //         least = -1;
                    //         break;
                    //     }
                    // }
                    sortlists_sg(&canon_sg);
                    assert(canon_sg.v == old_v);
                    assert(canon_sg.d == old_d);
                    assert(canon_sg.e == old_e);
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
    std::deque<Vertex> biggest_first_queue;
    std::vector<Edge_weight> dp(sg.nv, 0);
    dp[start] = 1;
    biggest_first_queue.push_back(start);
    while(biggest_first_queue.front() != 1) {
        auto cur_vertex = biggest_first_queue.front();
        biggest_first_queue.pop_front();
        for(int i = 0; i < sg.d[cur_vertex]; i++) {
            Vertex w = sg.e[sg.v[cur_vertex] + i];
            if(distance_to_goal[w] == invalid_ || dp[cur_vertex] == 0) {
                continue;
            }
            if(distance_to_goal[w] + 1 == distance_to_goal[cur_vertex]) {
                dp[w] += dp[cur_vertex];
                biggest_first_queue.push_back(w);
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

void ParallelSearch::pruning_dijkstra(sparsegraph const& sg, std::vector<Edge_length>& distance, Edge_length budget) {
    std::deque<Vertex> queue;
    queue.push_back(1);
    distance[1] = 0;
    while(!queue.empty()) {
        auto cur_vertex = queue.front();
        auto cur_cost = distance[cur_vertex];
        queue.pop_front();
        for(int i = 0; i < sg.d[cur_vertex]; i++) {
            Vertex w = sg.e[sg.v[cur_vertex] + i];
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


void ParallelSearch::reverse_pruning_dijkstra(sparsegraph const& sg, Vertex prune, std::vector<Edge_length>& distance, std::vector<Edge_length> const& forward_distance, Edge_length budget) {
    std::deque<Vertex> queue;
    queue.push_back(prune);
    distance[prune] = 0;
    while(!queue.empty()) {
        auto cur_vertex = queue.front();
        auto cur_cost = distance[cur_vertex];
        queue.pop_front();
        for(int i = 0; i < sg.d[cur_vertex]; i++) {
            Vertex w = sg.e[sg.v[cur_vertex] + i];
            if(forward_distance[w] == invalid_ || cur_cost + 1 >= distance[w] || cur_cost + 1 + forward_distance[w] > budget) {
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