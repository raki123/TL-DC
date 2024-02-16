// TLDC "Too long; Didn't Count" A length limited path counter.
// Copyright (C) 2023-2024 Rafael Kiesel, Markus Hecher

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

#include "parallel_directed_search.h"
#include "logging.h"
#include <deque>
#include <limits>

namespace fpc {

template <template <typename> typename Count_structure, typename count_t>
ParallelDirectedSearch<Count_structure, count_t>::ParallelDirectedSearch(
    Digraph const &input, std::optional<Vertex> must_use, size_t nthreads)
    : nr_vertices_(input.nr_vertices()), nthreads_(nthreads),
      max_length_(input.max_length()), terminals_(input.terminals()),
      in_neighbors_(nr_vertices_), out_neighbors_(nr_vertices_),
      must_use_(must_use),
      invalid_(std::numeric_limits<Edge_length>::max() - max_length_ - 1),
      distance_(nr_vertices_, std::vector<Edge_length>(nr_vertices_, invalid_)),
      canonizer_(must_use_ ? Canonizer::noop_canonizer()
                           : Canonizer(input.equivalence_classes(), invalid_,
                                       terminals_[1])),
      cache_(nr_vertices_ + 1), thread_local_result_(nthreads_),
      pos_hits_(nthreads_, 0), neg_hits_(nthreads_, 0), edges_(nthreads_, 0),
      propagations_(nthreads_, 0), dags_(nthreads_, 0) {
  assert(terminals_.size() == 2);
  for (Vertex v = 0; v < nr_vertices_; v++) {
    // fill neighbors
    in_neighbors_[v] = input.in_neighbors(v);
    out_neighbors_[v] = input.out_neighbors(v);
  }
  for (Vertex v = 0; v < nr_vertices_; v++) {
    dijkstra(v, distance_[v], true);
  }
  omp_set_num_threads(nthreads_);
}

template <template <typename> typename Count_structure, typename count_t>
mpz_class ParallelDirectedSearch<Count_structure, count_t>::search() {
  CacheKey first_key;
  first_key.start = terminals_[0];
  first_key.distance_to_goal = std::vector<Edge_length>(nr_vertices_, invalid_);
  dijkstra(terminals_[1], first_key.distance_to_goal, false);
  first_key.distance_to_goal[terminals_[0]] = invalid_;
  cache_[nr_vertices_ - 1].emplace(first_key, Partial_result(max_length_));

  for (Vertex remaining_size = nr_vertices_; remaining_size-- > 0;) {
    std::size_t entries = cache_[remaining_size].size();
    LOG << "\r" << nr_vertices_ - remaining_size << " / " << nr_vertices_;
    LOG << ", with " << entries << " entries.";
#pragma omp parallel for default(shared)
    for (size_t bucket = 0; bucket < cache_[remaining_size].bucket_count();
         bucket++) {

      size_t thread_id = omp_get_thread_num();
      for (auto task_it = cache_[remaining_size].begin(bucket);
           task_it != cache_[remaining_size].end(bucket); ++task_it) {
        auto const &[start, old_distance_to_goal] = task_it->first;
        auto const &result = task_it->second;
        auto initially_satisfied =
            must_use_ ? old_distance_to_goal[*must_use_] == invalid_ : true;
        Edge_length budget = max_length_ - result.offset();
        std::vector<Vertex> poss;
        for (auto v : out_neighbors(start)) {
          if (v == terminals_[1]) {
            if (!initially_satisfied) {
              // we did not visit the must use vertex and cannot include the
              // found paths in the result
              continue;
            }
            edges_[thread_id]++;
            // update the partial result
            thread_local_result_[thread_id] +=
                result.total_count(max_length_ - 1);
            continue;
          }
          if (budget >= old_distance_to_goal[v] + 1) {
            edges_[thread_id]++;
            poss.push_back(v);
          }
        }
        for (auto v : poss) {
          // update the partial result
          auto new_result = result;
          new_result.increment_offset();
          auto currently_satisfied = initially_satisfied || v == *must_use_;
          Edge_length v_budget = budget - 1;
          std::vector<Edge_length> distance_to_goal(nr_vertices_, invalid_);
          // make sure we disable paths going through v
          distance_to_goal[v] = 0;
          pruning_dijkstra(terminals_[1], v, distance_to_goal,
                           old_distance_to_goal, v_budget);
          if (!currently_satisfied &&
              distance_to_goal[*must_use_] == invalid_) {
            // we cannot visit the must use vertex anymore if we continue
            // along v.
            // -> skip v
            continue;
          }
          // unset the hack value for v
          distance_to_goal[v] = invalid_;
          std::vector<Vertex> poss_non_dag;
          bool found_goal = false;
          for (auto w : out_neighbors(v)) {
            if (w == terminals_[1]) {
              found_goal = true;
              continue;
            }
            if (v_budget >= distance_to_goal[w] + 1) {
              poss_non_dag.push_back(w);
            }
          }
          // there is only one non_dag edge
          std::vector<Vertex> extra = {v};
          Vertex last = v;
          while (poss_non_dag.size() == 1) {
            edges_[thread_id]++;
            propagations_[thread_id]++;
            // first check if the goal is a neighbor
            if (found_goal && currently_satisfied) {
              edges_[thread_id]++;
              // update the partial result
              thread_local_result_[thread_id] +=
                  new_result.total_count(max_length_ - 1);
            }
            Vertex cur = poss_non_dag[0];
            extra.push_back(cur);
            new_result.increment_offset();
            currently_satisfied = currently_satisfied || cur == *must_use_;
            v_budget--;
            std::fill(distance_to_goal.begin(), distance_to_goal.end(),
                      invalid_);
            for (Vertex extra_v : extra) {
              distance_to_goal[extra_v] = 0;
            }
            pruning_dijkstra(terminals_[1], cur, distance_to_goal,
                             old_distance_to_goal, v_budget);
            for (Vertex extra_v : extra) {
              distance_to_goal[extra_v] = invalid_;
            }
            poss_non_dag.clear();
            found_goal = false;
            for (auto w : out_neighbors(cur)) {
              if (w == terminals_[1]) {
                found_goal = true;
                continue;
              }
              if (v_budget >= distance_to_goal[w] + 1) {
                poss_non_dag.push_back(w);
              }
            }
            last = cur;
          }
          if (poss_non_dag.size() == 0) {
            if (found_goal && currently_satisfied) {
              edges_[thread_id]++;
              // update the partial result
              thread_local_result_[thread_id] +=
                  new_result.total_count(max_length_ - 1);
            }
            continue;
          }
          prune_domination(last, distance_to_goal);
          if (!currently_satisfied &&
              distance_to_goal[*must_use_] == invalid_) {
            // we cannot visit the must use vertex anymore if we continue
            // here
            // -> skip
            continue;
          }
          assert(last != terminals_[1]);
          assert(last < nr_vertices_);
          assert(new_result.offset() < max_length_);
          Vertex new_remaining_size = 0;
          for (Vertex i = 0; i < nr_vertices_; i++) {
            if (distance_to_goal[i] != invalid_) {
              new_remaining_size++;
            }
          }
          CacheKey new_key{last, std::move(distance_to_goal)};
          canonizer_(new_key.start, new_key.distance_to_goal);
#pragma omp critical
          {
            auto ins = cache_[new_remaining_size].insert(
                std::make_pair(new_key, new_result));
            if (!ins.second) {
              pos_hits_[thread_id]++;
              // there is already an element with that key
              // instead increase the partial result for that key
              ins.first->second += new_result;
            } else {
              neg_hits_[thread_id]++;
            }
          }
        }
      }
    }
    cache_[remaining_size].clear();
  }

  mpz_class res = 0;
  for (size_t id = 0; id < nthreads_; id++) {
    res += to_mpz(thread_local_result_[id]);
  }
  return res;
}

template <template <typename> typename Count_structure, typename count_t>
void ParallelDirectedSearch<Count_structure, count_t>::prune_domination(
    Vertex start, std::vector<Edge_length> &distance) {
  std::vector<Vertex> ap_disc(in_neighbors_.size(), 0);
  std::vector<Vertex> ap_low(in_neighbors_.size(), 0);
  std::vector<char> ap_visited(in_neighbors_.size(), false);
  int time = 0;
  dom_util(start, ap_visited, ap_disc, ap_low, time, -1, start, distance);
}

template <template <typename> typename Count_structure, typename count_t>
bool ParallelDirectedSearch<Count_structure, count_t>::dom_util(
    Vertex u, std::vector<char> &visited, std::vector<Vertex> &disc,
    std::vector<Vertex> &low, int &time, int parent, Vertex start,
    std::vector<Edge_length> &distance) {
  // We do not need to check whether the root is an articulation point
  // since if it is, then the goal can only be in one of the induced
  // components for all the other components we cannot enter them since we
  // cannot reach the goal from them but this means that we have already
  // pruned them using dijkstra

  // Mark the current node as visited
  visited[u] = true;

  bool found_elsewhere = false;

  // Initialize discovery time and low value
  disc[u] = low[u] = ++time;

  // Go through all vertices adjacent to this
  for (auto v : out_neighbors(u)) {
    // If v is not visited yet, then make it a child of u
    // in DFS tree and recur for it
    if (v != start && distance[v] == invalid_) {
      continue;
    }
    if (!visited[v]) {
      bool found_here =
          dom_util(v, visited, disc, low, time, u, start, distance);
      found_elsewhere |= found_here;

      // Check if the subtree rooted with v has
      // a connection to one of the ancestors of u
      low[u] = std::min(low[u], low[v]);

      // If u is not root and low value of one of
      // its child is more than discovery value of u.
      if (parent != -1 && low[v] >= disc[u]) {
        // AP
        if (!found_here) {
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

  for (auto v : in_neighbors(u)) {
    // If v is not visited yet, then make it a child of u
    // in DFS tree and recur for it
    if (v != start && distance[v] == invalid_) {
      continue;
    }
    if (!visited[v]) {
      bool found_here =
          dom_util(v, visited, disc, low, time, u, start, distance);
      found_elsewhere |= found_here;

      // Check if the subtree rooted with v has
      // a connection to one of the ancestors of u
      low[u] = std::min(low[u], low[v]);

      // If u is not root and low value of one of
      // its child is more than discovery value of u.
      if (parent != -1 && low[v] >= disc[u]) {
        // AP
        if (!found_here) {
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

template <template <typename> typename Count_structure, typename count_t>
void ParallelDirectedSearch<Count_structure, count_t>::prune_util(
    Vertex u, std::vector<Edge_length> &distance) {
  for (auto v : out_neighbors(u)) {
    if (distance[v] != invalid_) {
      distance[v] = invalid_;
      prune_util(v, distance);
    }
  }
  for (auto v : in_neighbors(u)) {
    if (distance[v] != invalid_) {
      distance[v] = invalid_;
      prune_util(v, distance);
    }
  }
}

template <template <typename> typename Count_structure, typename count_t>
void ParallelDirectedSearch<Count_structure, count_t>::dijkstra(
    Vertex start, std::vector<Edge_length> &distance, bool outgoing) {
  std::deque<Vertex> queue;
  queue.push_back(start);
  distance[start] = 0;
  while (!queue.empty()) {
    auto cur_vertex = queue.front();
    auto cur_cost = distance[cur_vertex];
    queue.pop_front();
    for (auto w :
         outgoing ? out_neighbors(cur_vertex) : in_neighbors(cur_vertex)) {
      if (cur_cost + 1 >= distance[w]) {
        continue;
      }
      if (cur_cost + 1 < distance[w]) {
        distance[w] = 1 + cur_cost;
        queue.push_back(w);
      }
    }
  }
}

template <template <typename> typename Count_structure, typename count_t>
void ParallelDirectedSearch<Count_structure, count_t>::pruning_dijkstra(
    Vertex start, Vertex prune, std::vector<Edge_length> &distance,
    std::vector<Edge_length> const &old_distance, Edge_length budget) {
  std::deque<Vertex> queue;
  queue.push_back(start);
  distance[start] = 0;
  while (!queue.empty()) {
    auto cur_vertex = queue.front();
    auto cur_cost = distance[cur_vertex];
    queue.pop_front();
    for (auto w : in_neighbors(cur_vertex)) {
      if (cur_cost + 1 >= distance[w] || old_distance[w] == invalid_) {
        continue;
      }
      if (cur_cost + 1 < distance[w] &&
          cur_cost + 1 + distance_[prune][w] <= budget) {
        distance[w] = 1 + cur_cost;
        queue.push_back(w);
      }
    }
  }
}

template <template <typename> typename Count_structure, typename count_t>
void ParallelDirectedSearch<Count_structure, count_t>::print_stats() const {
  size_t pos_hits = 0, neg_hits = 0;
  for (size_t i = 0; i < nthreads_; i++) {
    pos_hits += pos_hits_[i];
    neg_hits += neg_hits_[i];
  }
  size_t dags = 0, splits = 0;
  for (size_t i = 0; i < nthreads_; i++) {
    dags += dags_[i];
  }
  size_t edges = 0, propagations = 0;
  for (size_t i = 0; i < nthreads_; i++) {
    edges += edges_[i];
    propagations += propagations_[i];
  }
  LOG << "Cache hit rate: " << 100 * pos_hits / (double)(pos_hits + neg_hits)
      << "% (" << pos_hits << "/" << pos_hits + neg_hits << ")" << std::endl;
  LOG << "#DAG searches: " << dags << " #Splits: " << splits << std::endl;
  LOG << "#Edges: " << edges << " #Propagations: " << propagations << std::endl;
}
} // namespace fpc
