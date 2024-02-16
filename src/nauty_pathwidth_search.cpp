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

#include "nauty_pathwidth_search.h"
#include "logging.h"
#include <algorithm>

namespace fpc {

template <template <typename> typename Count_structure, typename count_t>
NautyPathwidthSearch<Count_structure, count_t>::NautyPathwidthSearch(
    Graph const &input, AnnotatedDecomposition decomposition, size_t nthreads)
    : nthreads_(nthreads), graph_(input), max_length_(graph_.max_length()),
      is_all_pair_(graph_.is_all_pair()), terminals_(graph_.terminals()),
      decomposition_(decomposition),
      remaining_edges_after_this_(decomposition_.size()),
      bag_local_idx_map_(
          decomposition_.size(),
          std::vector<frontier_index_t>(graph_.nr_vertices(), invalid_index_)),
      bag_local_vertex_map_(decomposition_.size()),
      bag_local_distance_(decomposition_.size()),
      sparsegraph_after_this_(decomposition_.size()),
      thread_local_result_(nthreads_, 0), cache_(decomposition_.size()),
      pos_hits_(nthreads_, 0), neg_hits_(nthreads_, 0), edges_(nthreads_, 0),
      propagations_(nthreads_, 0) {

  omp_set_num_threads(nthreads_);

  if (!is_all_pair_) {
    decomposition_.foreach_post_order([this](auto bag_idx) {
      auto &bag = decomposition_[bag_idx].bag;
      if (std::find(bag.begin(), bag.end(), terminals_[0]) == bag.end()) {
        bag.push_back(terminals_[0]);
      }
      if (std::find(bag.begin(), bag.end(), terminals_[1]) == bag.end()) {
        bag.push_back(terminals_[1]);
      }
    });
  }

  std::set<Vertex> empty;
  std::vector<Vertex> all;
  for (Vertex i = 0; i < graph_.nr_vertices(); i++) {
    all.push_back(i);
  }
  Graph cur_graph = graph_.subgraph(all);
  decomposition_.foreach_post_order([&](auto bag_idx) {
    auto &node = decomposition_[bag_idx];
    auto &idx = bag_local_idx_map_[bag_idx];
    auto &vertex = bag_local_vertex_map_[bag_idx];
    auto &remaining = remaining_edges_after_this_[bag_idx];
    switch (node.type) {
    case NodeType::LEAF:
      cur_graph = graph_.subgraph(all);
      break;
    case NodeType::PATH_LIKE:
      break;
    case NodeType::JOIN:
      cur_graph = graph_.subgraph(all);
      decomposition_.foreach_post_order(
          [&](auto top) {
            auto &top_node = decomposition_[top];
            switch (top_node.type) {
            case NodeType::LEAF:
            case NodeType::PATH_LIKE:
              cur_graph.remove_edge(top_node.edge.first, top_node.edge.second);
              break;
            case NodeType::JOIN:
              break;
            }
          },
          bag_idx);
      break;
    }
    for (size_t i = 0; i < node.bag.size(); i++) {
      if (!is_all_pair_ &&
          (node.bag[i] == terminals_[0] || node.bag[i] == terminals_[1])) {
        continue;
      }
      if (node.type != NodeType::JOIN &&
          cur_graph.neighbors(node.bag[i]).size() == 0) {
        std::swap(node.bag[i], node.bag.back());
        node.bag.pop_back();
        i--;
      } else if (cur_graph.neighbors(node.bag[i]).size() ==
                     graph_.neighbors(node.bag[i]).size() &&
                 node.bag[i] != node.edge.first &&
                 node.bag[i] != node.edge.second) {
        std::swap(node.bag[i], node.bag.back());
        node.bag.pop_back();
        i--;
      } else if (node.type == NodeType::JOIN) {
        auto v = node.bag[i];
        auto left_child = node.children.first;
        auto right_child = node.children.second;
        auto left_idx = bag_local_idx_map_[left_child][v];
        auto right_idx = bag_local_idx_map_[right_child][v];
        if (right_idx != invalid_index_ &&
            remaining_edges_after_this_[right_child][right_idx] == 0) {
          assert(left_idx == invalid_index_);
          std::swap(node.bag[i], node.bag.back());
          node.bag.pop_back();
          i--;
        } else if (left_idx != invalid_index_ &&
                   remaining_edges_after_this_[left_child][left_idx] == 0) {
          assert(right_idx == invalid_index_);
          std::swap(node.bag[i], node.bag.back());
          node.bag.pop_back();
          i--;
        }
      }
    }
    std::sort(node.bag.begin(), node.bag.end());
    bag_local_distance_[bag_idx].resize(node.bag.size());
    remaining.resize(node.bag.size());
    if (node.type != NodeType::JOIN) {
      cur_graph.remove_edge(node.edge.first, node.edge.second);
    }
    frontier_index_t cur_idx = 0;
    for (auto v : node.bag) {
      assert(cur_idx != invalid_index_);
      idx[v] = cur_idx++;
      vertex.push_back(v);
      remaining[idx[v]] = cur_graph.neighbors(v).size();
    }
    for (auto v : node.bag) {
      std::vector<Edge_length> distance(graph_.nr_vertices(),
                                        invalid_distance_);
      cur_graph.dijkstra(v, distance, empty);
      std::vector<Edge_length> local_distance(node.bag.size(),
                                              invalid_distance_);
      for (auto other : node.bag) {
        if (other != v) {
          local_distance[idx[other]] = distance[other];
        }
      }
      bag_local_distance_[bag_idx][idx[v]] = local_distance;
    }
    auto &sg = sparsegraph_after_this_[bag_idx];
    SG_INIT(sg);
    size_t nr_vertices = cur_graph.nr_vertices();
    for (auto v : node.bag) {
      if (cur_graph.neighbors(v).empty()) {
        nr_vertices++;
      }
    }
    nr_vertices += node.bag.size();
    size_t nr_edges = 2 * cur_graph.nr_edges();
    nr_edges += 3 * node.bag.size();
    std::vector<size_t> new_name(graph_.nr_vertices(), size_t(-1));
    std::vector<size_t> reverse(nr_vertices, size_t(-1));
    for (auto v : node.bag) {
      new_name[v] = bag_local_idx_map_[bag_idx][v] + node.bag.size();
      reverse[new_name[v]] = v;
    }
    size_t cur_name = 2 * node.bag.size();
    for (Vertex v = 0; v < graph_.nr_vertices(); v++) {
      if (cur_graph.neighbors(v)
              .empty()) { // && v != node.edge.first && v != node.edge.second) {
        continue;
      }
      if (new_name[v] == size_t(-1)) {
        new_name[v] = cur_name++;
        reverse[new_name[v]] = v;
      }
    }
    sg.v = (edge_t *)calloc(sizeof(edge_t) * (nr_vertices + nr_edges) +
                                sizeof(degree_t) * nr_vertices,
                            1);
    sg.d = (degree_t *)(sg.v + nr_vertices);
    sg.e = (edge_t *)(sg.d + nr_vertices);
    sg.nv = nr_vertices;
    sg.nde = nr_edges - 3 * node.bag.size();
    sg.vlen = nr_vertices;
    sg.dlen = nr_vertices;
    sg.elen = nr_edges;
    size_t cur_e_idx = 0;
    for (size_t i = 0; i < node.bag.size(); i++) {
      sg.v[i] = cur_e_idx;
      sg.d[i] = 0;
      cur_e_idx += 2;
    }
    for (size_t i = node.bag.size(); i < nr_vertices; i++) {
      size_t v = reverse[i];
      assert(v != size_t(-1));
      sg.v[i] = cur_e_idx;
      for (auto w : cur_graph.neighbors(v)) {
        sg.e[cur_e_idx++] = new_name[w];
      }
      sg.d[i] = cur_e_idx - sg.v[i];
      if (i < 2 * node.bag.size()) {
        cur_e_idx++;
      }
    }
    assert(cur_e_idx == nr_edges);
  });
  for (Vertex v = 0; v < graph_.nr_vertices(); v++) {
    assert(cur_graph.neighbors(v).size() == 0);
  }
  decomposition_.foreach_post_order([&](auto bag_idx) {
    auto &node = decomposition_[bag_idx];
    if (node.type == NodeType::LEAF) {
      Frontier initial_frontier(decomposition_[bag_idx].bag.size(),
                                no_edge_index_);
      if (!is_all_pair_) {
        assert(bag_local_idx_map_[bag_idx][terminals_[0]] != invalid_index_);
        assert(bag_local_idx_map_[bag_idx][terminals_[1]] != invalid_index_);
        initial_frontier[bag_local_idx_map_[bag_idx][terminals_[0]]] =
            invalid_index_;
        initial_frontier[bag_local_idx_map_[bag_idx][terminals_[1]]] =
            invalid_index_;
      }
      auto initial_result = Partial_result(max_length_);
      auto sg = construct_sparsegraph(initial_frontier, size_t(-1));
      cache_[bag_idx].first.emplace(std::make_pair(sg, initial_frontier),
                                    initial_result);
    }
  });
}

template <template <typename> typename Count_structure, typename count_t>
mpz_class NautyPathwidthSearch<Count_structure, count_t>::search() {
  for (size_t bag_idx = 0; bag_idx < decomposition_.size(); bag_idx++) {
    LOG << "\r" << bag_idx << " / " << decomposition_.size() - 1 << " of type ";
    if (decomposition_[bag_idx].type == NodeType::LEAF) {
      LOG << "leaf, with " << cache_[bag_idx].first.size() << " entries.";
    } else if (decomposition_[bag_idx].type == NodeType::PATH_LIKE) {
      LOG << "path, with " << cache_[bag_idx].first.size() << " entries.";
    } else if (decomposition_[bag_idx].type == NodeType::JOIN) {
      LOG << "join, with (" << cache_[bag_idx].first.size() << ","
          << cache_[bag_idx].second.size() << ") entries.";
    }
    if (decomposition_[bag_idx].type == LEAF ||
        decomposition_[bag_idx].type == PATH_LIKE) {
// PATH_LIKE/LEAF
#pragma omp parallel for default(shared)
      for (size_t bucket = 0; bucket < cache_[bag_idx].first.bucket_count();
           bucket++) {
        size_t thread_id = omp_get_thread_num();
        for (auto task_it = cache_[bag_idx].first.begin(bucket);
             task_it != cache_[bag_idx].first.end(bucket); ++task_it) {
          auto copy_frontier = task_it->first.second;
          auto copy_result = task_it->second;
          propagateLoop(copy_frontier, bag_idx, -1, copy_result, true, false,
                        thread_id);
          copy_frontier = task_it->first.second;
          copy_result = task_it->second;
          propagateLoop(copy_frontier, bag_idx, -1, copy_result, false, true,
                        thread_id);
          free(task_it->first.first.v);
        }
      }
    } else {
      assert(false);
    }

    cache_[bag_idx] = {};
  }

  LOG << std::endl;
  print_stats();

  mpz_class rv = 0;
  for (size_t id = 0; id < nthreads_; id++) {
    rv += to_mpz(thread_local_result_[id]);
  }
  for (size_t bag_idx = 0; bag_idx < decomposition_.size(); bag_idx++) {
    free(sparsegraph_after_this_[bag_idx].v);
  }
  return rv;
}

template <template <typename> typename Count_structure, typename count_t>
void NautyPathwidthSearch<Count_structure, count_t>::propagateLoop(
    Frontier &frontier, size_t bag_idx, size_t last_idx,
    Partial_result &partial_results, bool takeable, bool skippable,
    size_t thread_id) {
  size_t new_idx = bag_idx;
  if (decomposition_[new_idx].type != JOIN && (takeable ^ skippable)) {
    edges_[thread_id]++;
    propagations_[thread_id]--;
  }
  // propagate while only one of the two is possible
  while ((takeable ^ skippable) && new_idx + 1 < decomposition_.size() &&
         decomposition_[new_idx].type != JOIN) {
    propagations_[thread_id]++;
    if (takeable) {
      take(frontier, new_idx);
      partial_results.increment_offset();
      size_t paths = 0;
      for (frontier_index_t idx : frontier) {
        if (idx <= 252) {
          paths++;
        } else if (idx == invalid_index_) {
          paths += 2;
        }
      }
      if (paths / 2 == 1) {
        thread_local_result_[thread_id] +=
            partial_results.total_count(max_length_);
      }
    } else {
      skip(frontier, new_idx);
    }
    last_idx = new_idx;
    new_idx = decomposition_[new_idx].parent;
    if (new_idx < decomposition_.size() &&
        decomposition_[new_idx].type != JOIN) {
      takeable = canTake(frontier, new_idx, partial_results);
      skippable = canSkip(frontier, new_idx, partial_results);
      if (!takeable) {
        includeSolutions(frontier, new_idx, partial_results);
      }
    }
  }
  // both are possible, so we have a new decision edge
  // put it into the cache
  if (new_idx < decomposition_.size() &&
      ((takeable && skippable) ||
       (decomposition_[new_idx].type == JOIN &&
        last_idx == decomposition_[new_idx].children.first))) {
    auto sg = construct_sparsegraph(frontier, last_idx);
#pragma omp critical
    {
      auto ins = cache_[new_idx].first.insert(
          std::make_pair(std::make_pair(sg, frontier), partial_results));
      if (!ins.second) {
        free(sg.v);
        pos_hits_[thread_id]++;
        // there is already an element with that key
        // instead increase the partial result for that key
        ins.first->second += partial_results;
      } else {
        neg_hits_[thread_id]++;
      }
    }
  } else if (new_idx < decomposition_.size() &&
             decomposition_[new_idx].type == JOIN) {
    assert(last_idx == decomposition_[new_idx].children.second);
    auto sg = construct_sparsegraph(frontier, last_idx);
#pragma omp critical
    {
      auto ins = cache_[new_idx].second.insert(
          std::make_pair(std::make_pair(sg, frontier), partial_results));
      if (!ins.second) {
        free(sg.v);
        pos_hits_[thread_id]++;
        // there is already an element with that key
        // instead increase the partial result for that key
        ins.first->second += partial_results;
      } else {
        neg_hits_[thread_id]++;
      }
    }
  }
}

template <template <typename> typename Count_structure, typename count_t>
void NautyPathwidthSearch<Count_structure, count_t>::includeSolutions(
    Frontier const &frontier, size_t bag_idx,
    Partial_result const &partial_result) {
  assert(partial_result.offset() + 1 <= max_length_);
  Edge edge = decomposition_[bag_idx].edge;
  auto v_idx = bag_local_idx_map_[bag_idx][edge.first];
  auto w_idx = bag_local_idx_map_[bag_idx][edge.second];
  // must not be closing a path to a loop
  if (frontier[v_idx] == w_idx) {
    assert(frontier[w_idx] == v_idx);
    return;
  }
  assert(frontier[w_idx] != v_idx);

  if (frontier[w_idx] == two_edge_index_ ||
      frontier[v_idx] == two_edge_index_) {
    return;
  }

  size_t thread_id = omp_get_thread_num();
  size_t number_paths = 0;
  size_t number_cut_paths = 0;
  for (frontier_index_t idx = 0; idx < frontier.size(); idx++) {
    if (frontier[idx] != no_edge_index_ && frontier[idx] != two_edge_index_) {
      number_paths++;
    }
    if (frontier[idx] == invalid_index_) {
      number_cut_paths++;
      number_paths--;
    }
  }
  assert(number_paths % 2 == 0);
  number_paths /= 2;
  number_paths += number_cut_paths;
  // no way this leads to a complete path
  if (number_paths > 2) {
    return;
  }
  // this can lead to a path if the first path and the second path connect
  if (number_paths == 2) {
    // cannot connect if they arent ends of paths
    if (frontier[v_idx] == no_edge_index_ ||
        frontier[v_idx] == two_edge_index_ ||
        frontier[w_idx] == no_edge_index_ ||
        frontier[w_idx] == two_edge_index_) {
      return;
    }
    // we connected them
    // add the partial result
    // current offset + 1 for the additional edge
    thread_local_result_[thread_id] +=
        partial_result.total_count(max_length_ - 1);
    return;
  }
  // this can lead to a path if we continue the existing path
  if (number_paths == 1) {
    // cannot connect if neither is the end of the path
    if ((frontier[v_idx] == no_edge_index_ ||
         frontier[v_idx] == two_edge_index_) &&
        (frontier[w_idx] == no_edge_index_ ||
         frontier[w_idx] == two_edge_index_)) {
      return;
    }
    // we connected them
    // add the partial result
    // current offset + 1 for the additional edge
    thread_local_result_[thread_id] +=
        partial_result.total_count(max_length_ - 1);
    return;
  }
  assert(number_paths == 0);
  assert(partial_result.offset() == 0);
  thread_local_result_[thread_id] += 1;
}

template <template <typename> typename Count_structure, typename count_t>
bool NautyPathwidthSearch<Count_structure, count_t>::canTake(
    Frontier &frontier, size_t bag_idx, Partial_result const &partial_result) {
  assert(partial_result.offset() + 1 <= max_length_);
  Edge edge = decomposition_[bag_idx].edge;
  auto v_idx = bag_local_idx_map_[bag_idx][edge.first];
  auto w_idx = bag_local_idx_map_[bag_idx][edge.second];

  // must not be closing a path to a loop
  if (frontier[v_idx] == w_idx) {
    assert(frontier[w_idx] == v_idx);
    return false;
  }
  assert(frontier[w_idx] != v_idx);

  // must not take more than two edges
  if (frontier[v_idx] == two_edge_index_ ||
      frontier[w_idx] == two_edge_index_) {
    return false;
  }

  // can connect two invalid paths but cannot continue after
  if (frontier[v_idx] == invalid_index_ && frontier[w_idx] == invalid_index_) {
    return false;
  }

  auto v_remaining = remaining_edges_after_this_[bag_idx][v_idx];
  auto w_remaining = remaining_edges_after_this_[bag_idx][w_idx];

  // can connect two invalid paths but cannot continue after
  if (!v_remaining && frontier[v_idx] == no_edge_index_ &&
      frontier[w_idx] == invalid_index_) {
    return false;
  }
  // can connect two invalid paths but cannot continue after
  if (!w_remaining && frontier[w_idx] == no_edge_index_ &&
      frontier[v_idx] == invalid_index_) {
    return false;
  }
  // can connect two invalid paths but cannot continue after
  if (!v_remaining && !w_remaining && frontier[v_idx] == no_edge_index_ &&
      frontier[w_idx] == no_edge_index_) {
    return false;
  }

  assert(w_idx < frontier.size());
  assert(v_idx < frontier.size());
  assert(frontier[w_idx] != v_idx);
  // length based pruning
  auto old_v = frontier[v_idx];
  auto old_w = frontier[w_idx];
  std::vector<std::pair<frontier_index_t, frontier_index_t>> restore(
      {std::make_pair(v_idx, old_v), std::make_pair(w_idx, old_w)});
  if (old_v == no_edge_index_) {
    if (old_w == no_edge_index_) {
      frontier[v_idx] = w_idx;
      frontier[w_idx] = v_idx;
    } else if (old_w == invalid_index_) {
      frontier[w_idx] = two_edge_index_;
      frontier[v_idx] = invalid_index_;
    } else {
      restore.push_back(std::make_pair(old_w, frontier[old_w]));
      frontier[old_w] = v_idx;
      frontier[w_idx] = two_edge_index_;
      frontier[v_idx] = old_w;
    }
  } else if (old_v == invalid_index_) {
    if (old_w == no_edge_index_) {
      frontier[v_idx] = two_edge_index_;
      frontier[w_idx] = invalid_index_;
    } else if (old_w == invalid_index_) {
      assert(false);
    } else {
      restore.push_back(std::make_pair(old_w, frontier[old_w]));
      frontier[old_w] = invalid_index_;
      frontier[w_idx] = two_edge_index_;
      frontier[v_idx] = two_edge_index_;
    }
  } else {
    restore.push_back(std::make_pair(old_v, frontier[old_v]));
    if (old_w == no_edge_index_) {
      frontier[old_v] = w_idx;
      frontier[v_idx] = two_edge_index_;
      frontier[w_idx] = old_v;
    } else if (old_w == invalid_index_) {
      frontier[old_v] = invalid_index_;
      frontier[v_idx] = two_edge_index_;
      frontier[w_idx] = two_edge_index_;
    } else {
      restore.push_back(std::make_pair(old_w, frontier[old_w]));
      frontier[old_w] = old_v;
      frontier[old_v] = old_w;
      frontier[w_idx] = two_edge_index_;
      frontier[v_idx] = two_edge_index_;
    }
  }
  if (v_remaining == 0) {
    assert(frontier[v_idx] != invalid_index_);
    assert(frontier[v_idx] != no_edge_index_);
    if (frontier[v_idx] != two_edge_index_) {
      frontier[frontier[v_idx]] = invalid_index_;
      frontier[v_idx] = two_edge_index_;
    }
  }
  if (w_remaining == 0) {
    assert(frontier[w_idx] != invalid_index_);
    assert(frontier[w_idx] != no_edge_index_);
    if (frontier[w_idx] != two_edge_index_) {
      frontier[frontier[w_idx]] = invalid_index_;
      frontier[w_idx] = two_edge_index_;
    }
  }
  std::vector<frontier_index_t> paths;
  std::vector<frontier_index_t> cut_paths;
  for (frontier_index_t idx = 0; idx < frontier.size(); idx++) {
    if (frontier[idx] != no_edge_index_ && frontier[idx] != two_edge_index_) {
      if (frontier[idx] == invalid_index_) {
        cut_paths.push_back(idx);
      } else {
        paths.push_back(idx);
      }
    }
  }
  if (cut_paths.size() > 2) {
    for (auto [idx, rest] : restore) {
      frontier[idx] = rest;
    }
    return false;
  }
  assert(paths.size() % 2 == 0);

  for (auto [idx, rest] : restore) {
    frontier[idx] = rest;
  }
  // + 1 - 1 since were taking the edge
  if (paths.size() / 2 + cut_paths.size() > 1) {
    if (paths.size() / 2 + cut_paths.size() + partial_result.offset() >
        max_length_) {
      return false;
    }
  } else if (paths.size() / 2 + cut_paths.size() + partial_result.offset() + 1 >
             max_length_) {
    return false;
  }

  return distancePrune(frontier, paths, cut_paths, bag_idx,
                       partial_result.offset() + 1);
}

template <template <typename> typename Count_structure, typename count_t>
bool NautyPathwidthSearch<Count_structure, count_t>::canSkip(
    Frontier &frontier, size_t bag_idx, Partial_result const &partial_result) {
  assert(partial_result.offset() + 1 <= max_length_);
  Edge edge = decomposition_[bag_idx].edge;
  auto v_idx = bag_local_idx_map_[bag_idx][edge.first];
  auto w_idx = bag_local_idx_map_[bag_idx][edge.second];
  auto v_remaining = remaining_edges_after_this_[bag_idx][v_idx];
  auto w_remaining = remaining_edges_after_this_[bag_idx][w_idx];

  // cannot let an outside path die
  if ((!v_remaining && frontier[v_idx] == invalid_index_) ||
      (!w_remaining && frontier[w_idx] == invalid_index_)) {
    return false;
  }
  // both ends will be outside the outside path dies
  if (v_remaining == 0 && w_remaining == 0 && frontier[v_idx] == w_idx) {
    return false;
  }
  // length based pruning
  std::vector<frontier_index_t> paths;
  std::vector<frontier_index_t> cut_paths;
  for (frontier_index_t idx = 0; idx < frontier.size(); idx++) {
    if (frontier[idx] != no_edge_index_ && frontier[idx] != two_edge_index_) {
      if (frontier[idx] == invalid_index_) {
        cut_paths.push_back(idx);
      } else {
        paths.push_back(idx);
      }
    }
  }
  assert(cut_paths.size() <= 2);
  assert(paths.size() % 2 == 0);
  bool last_and_incomplete_v = v_remaining == 0 &&
                               frontier[v_idx] != two_edge_index_ &&
                               frontier[v_idx] != no_edge_index_;
  bool last_and_incomplete_w = w_remaining == 0 &&
                               frontier[w_idx] != two_edge_index_ &&
                               frontier[w_idx] != no_edge_index_;
  if (last_and_incomplete_v && last_and_incomplete_w) {
    if (cut_paths.size() > 0) {
      return false;
    }
  } else if (last_and_incomplete_v || last_and_incomplete_w) {
    if (cut_paths.size() > 1) {
      return false;
    }
  }

  return distancePrune(frontier, paths, cut_paths, bag_idx,
                       partial_result.offset());
}

template <template <typename> typename Count_structure, typename count_t>
bool NautyPathwidthSearch<Count_structure, count_t>::distancePrune(
    Frontier &frontier, std::vector<frontier_index_t> const &paths,
    std::vector<frontier_index_t> const &cut_paths, size_t bag_idx,
    size_t offset) {
  // advanced length based pruning
  auto &distance = bag_local_distance_[bag_idx];
  if (cut_paths.size() == 2 && paths.size() == 0) {
    // we have to connect the cut paths
    return distance[cut_paths[0]][cut_paths[1]] + offset <= max_length_;
  }
  if (cut_paths.size() + paths.size() <= 1) {
    return offset + 1 <= max_length_;
  }
  // there is at least one other path that we have to incorporate
  // connecting the cut paths is not an option
  size_t min_dist = 0;
  for (auto v : cut_paths) {
    Edge_length cur_min_dist = invalid_distance_;
    for (auto other : paths) {
      cur_min_dist = std::min(cur_min_dist, distance[v][other]);
    }
    if (cur_min_dist == invalid_distance_) {
      // we cannot connect this cut path
      return false;
    }
    min_dist += cur_min_dist;
  }
  if (min_dist + offset > max_length_) {
    return false;
  }
  size_t path_min_dist = 0;
  Edge_length path_max_dist = 0;
  Edge_length path_snd_max_dist = 0;
  for (auto v : paths) {
    Edge_length cur_min_dist = invalid_distance_;
    for (auto cut_end : cut_paths) {
      cur_min_dist = std::min(cur_min_dist, distance[v][cut_end]);
    }
    for (auto other : paths) {
      if (frontier[v] != other) {
        cur_min_dist = std::min(cur_min_dist, distance[v][other]);
      }
    }
    if (cur_min_dist >= path_max_dist) {
      path_snd_max_dist = path_max_dist;
      path_max_dist = cur_min_dist;
    } else if (cur_min_dist > path_snd_max_dist) {
      path_snd_max_dist = cur_min_dist;
    }
    path_min_dist += cur_min_dist;
  }
  assert(path_min_dist >= path_max_dist + path_snd_max_dist);
  path_min_dist -= path_max_dist;
  path_min_dist -= path_snd_max_dist;
  path_min_dist /= 2;
  if (path_min_dist + min_dist + offset > max_length_) {
    return false;
  }
  return true;
}

template <template <typename> typename Count_structure, typename count_t>
void NautyPathwidthSearch<Count_structure, count_t>::take(Frontier &frontier,
                                                          size_t bag_idx) {
  Edge edge = decomposition_[bag_idx].edge;
  auto v_idx = bag_local_idx_map_[bag_idx][edge.first];
  auto w_idx = bag_local_idx_map_[bag_idx][edge.second];
  auto old_v = frontier[v_idx];
  auto old_w = frontier[w_idx];

  if (old_v == no_edge_index_) {
    if (old_w == no_edge_index_) {
      frontier[v_idx] = w_idx;
      frontier[w_idx] = v_idx;
    } else if (old_w == invalid_index_) {
      frontier[w_idx] = two_edge_index_;
      frontier[v_idx] = invalid_index_;
    } else {
      frontier[old_w] = v_idx;
      frontier[w_idx] = two_edge_index_;
      frontier[v_idx] = old_w;
    }
  } else if (old_v == invalid_index_) {
    if (old_w == no_edge_index_) {
      frontier[v_idx] = two_edge_index_;
      frontier[w_idx] = invalid_index_;
    } else if (old_w == invalid_index_) {
      assert(false);
    } else {
      frontier[old_w] = invalid_index_;
      frontier[w_idx] = two_edge_index_;
      frontier[v_idx] = two_edge_index_;
    }
  } else {
    if (old_w == no_edge_index_) {
      frontier[old_v] = w_idx;
      frontier[v_idx] = two_edge_index_;
      frontier[w_idx] = old_v;
    } else if (old_w == invalid_index_) {
      frontier[old_v] = invalid_index_;
      frontier[v_idx] = two_edge_index_;
      frontier[w_idx] = two_edge_index_;
    } else {
      frontier[old_w] = old_v;
      frontier[old_v] = old_w;
      frontier[w_idx] = two_edge_index_;
      frontier[v_idx] = two_edge_index_;
    }
  }

  advance(frontier, bag_idx);
}

template <template <typename> typename Count_structure, typename count_t>
void NautyPathwidthSearch<Count_structure, count_t>::skip(Frontier &frontier,
                                                          size_t bag_idx) {
  advance(frontier, bag_idx);
}

template <template <typename> typename Count_structure, typename count_t>
void NautyPathwidthSearch<Count_structure, count_t>::advance(Frontier &frontier,
                                                             size_t bag_idx) {
  Frontier old = frontier;
  size_t next_idx = decomposition_[bag_idx].parent;
  auto &bag = decomposition_[next_idx].bag;
  frontier.resize(bag.size());
  std::fill(frontier.begin(), frontier.end(), no_edge_index_);
  auto &old_idx = bag_local_idx_map_[bag_idx];
  auto &new_idx = bag_local_idx_map_[next_idx];
  auto &old_vertex = bag_local_vertex_map_[bag_idx];
  int found_two = 0;
  int found_invalid = 0;
  int found_path = 0;
  for (auto v : bag) {
    if (old_idx[v] != invalid_index_) {
      if (old[old_idx[v]] == no_edge_index_) {
        frontier[new_idx[v]] = no_edge_index_;
      } else if (old[old_idx[v]] == two_edge_index_) {
        frontier[new_idx[v]] = two_edge_index_;
        found_two += 1;
      } else if (old[old_idx[v]] == invalid_index_) {
        frontier[new_idx[v]] = invalid_index_;
        found_invalid += 1;
      } else {
        frontier[new_idx[v]] = new_idx[old_vertex[old[old_idx[v]]]];
        if (frontier[new_idx[v]] != invalid_index_) {
          found_path += 1;
        } else {
          found_path += 2;
        }
      }
      assert(frontier[new_idx[v]] == no_edge_index_ ||
             frontier[new_idx[v]] == two_edge_index_ ||
             remaining_edges_after_this_[bag_idx][old_idx[v]] > 0);
    }
  }
  assert(found_path % 2 == 0);
  assert(!found_two || found_invalid || found_path);
  assert(found_invalid <= 2);
}

template <template <typename> typename Count_structure, typename count_t>
sparsegraph
NautyPathwidthSearch<Count_structure, count_t>::construct_sparsegraph(
    Frontier const &frontier, size_t last_idx) {
  if (last_idx == size_t(-1)) {
    // LEAF
    // just take full graph
    return graph_.to_canon_nauty(false);
  }
  size_t bag_idx = decomposition_[last_idx].parent;
  sparsegraph sg;
  SG_INIT(sg);
  auto const &base_sg = sparsegraph_after_this_[last_idx];
  sg.v = (edge_t *)malloc(sizeof(edge_t) * (base_sg.vlen + base_sg.elen) +
                          sizeof(degree_t) * base_sg.dlen);
  sg.d = (degree_t *)(sg.v + base_sg.vlen);
  sg.e = (edge_t *)(sg.d + base_sg.vlen);
  sg.nv = base_sg.nv;
  sg.nde = base_sg.nde;
  sg.vlen = base_sg.vlen;
  sg.dlen = base_sg.vlen;
  sg.elen = base_sg.elen;
  std::memcpy(sg.v, base_sg.v,
              sizeof(edge_t) * (base_sg.vlen + base_sg.elen) +
                  sizeof(degree_t) * base_sg.dlen);
  // std::cout << sg.nv << " " << sg.nde << std::endl;
  // std::cout << sg.vlen << " " << sg.dlen << " " << sg.elen << std::endl;
  // size_t e_count = 0;
  // for(size_t v = 0; v < sg.nv; v++) {
  //     std::vector<int> neighs(sg.nv, 0);
  //     for(size_t i = sg.v[v]; i < sg.v[v] + sg.d[v]; i++) {
  //         neighs[sg.e[i]] = 1;
  //         e_count++;
  //     }
  //     for(size_t i = 0; i < sg.nv; i++) {
  //         std::cout << neighs[i];
  //     }
  //     std::cout << std::endl;
  // }
  // std::cout << std::endl;
  int *lab = (int *)malloc(sizeof(int) * sg.nv);
  int *ptn = (int *)malloc(sizeof(int) * sg.nv);
  int *orbits = (int *)malloc(sizeof(int) * sg.nv);
  std::vector<frontier_index_t> two_edge;
  std::vector<frontier_index_t> no_edge;
  std::vector<frontier_index_t> invalid;
  std::vector<frontier_index_t> active;
  size_t offset = decomposition_[last_idx].bag.size();
  // now sg is just a copy of base_sg
  // we add the edges of the frontier
  for (frontier_index_t idx = 0; idx < frontier.size(); idx++) {
    auto v = bag_local_vertex_map_[bag_idx][idx];
    auto old_v_idx = bag_local_idx_map_[last_idx][v];
    if (old_v_idx == invalid_index_) {
      assert(frontier[idx] == no_edge_index_);
      continue;
    }
    if (frontier[idx] == no_edge_index_) {
      // dont do anything
      no_edge.push_back(old_v_idx);
    } else if (frontier[idx] == two_edge_index_) {
      // remove adjacent
      two_edge.push_back(old_v_idx);
    } else if (frontier[idx] == invalid_index_) {
      // single connect fake node to real node
      assert(sg.d[old_v_idx] == 0);
      assert(sg.e[sg.v[old_v_idx]] == 0);
      sg.e[sg.v[old_v_idx]] = old_v_idx + offset;
      sg.d[old_v_idx]++;
      assert(sg.d[old_v_idx + offset] > 0);
      assert(sg.e[sg.v[old_v_idx + offset] + sg.d[old_v_idx + offset]] == 0);
      sg.e[sg.v[old_v_idx + offset] + sg.d[old_v_idx + offset]] = old_v_idx;
      sg.d[old_v_idx + offset]++;
      sg.nde += 2;
      invalid.push_back(old_v_idx);
    } else if (idx < frontier[idx]) {
      // add the edges
      auto w = bag_local_vertex_map_[bag_idx][frontier[idx]];
      auto old_w_idx = bag_local_idx_map_[last_idx][w];
      assert(old_w_idx != invalid_index_);
      sg.nde += 4;
      assert(sg.e[sg.v[old_v_idx] + sg.d[old_v_idx]] == 0);
      sg.e[sg.v[old_v_idx] + sg.d[old_v_idx]++] = old_w_idx + offset;
      assert(sg.e[sg.v[old_v_idx] + sg.d[old_v_idx]] == 0);
      sg.e[sg.v[old_v_idx] + sg.d[old_v_idx]++] = old_v_idx + offset;
      assert(sg.e[sg.v[old_v_idx + offset] + sg.d[old_v_idx + offset]] == 0);
      assert(sg.e[sg.v[old_w_idx + offset] + sg.d[old_w_idx + offset]] == 0);
      sg.e[sg.v[old_v_idx + offset] + sg.d[old_v_idx + offset]++] = old_v_idx;
      sg.e[sg.v[old_w_idx + offset] + sg.d[old_w_idx + offset]++] = old_v_idx;
      active.push_back(old_v_idx);
      active.push_back(old_w_idx);
    }
  }
  size_t found = 0;
  for (auto idx : no_edge) {
    lab[found] = idx;
    ptn[found++] = 1;
  }
  if (found > 0) {
    ptn[found - 1] = 0;
  }
  for (auto idx : two_edge) {
    lab[found] = idx;
    ptn[found++] = 1;
  }
  if (found > 0) {
    ptn[found - 1] = 0;
  }
  for (auto idx : invalid) {
    lab[found] = idx;
    ptn[found++] = 1;
  }
  if (found > 0) {
    ptn[found - 1] = 0;
  }
  for (auto idx : active) {
    lab[found] = idx;
    ptn[found++] = 1;
  }
  if (found > 0) {
    ptn[found - 1] = 0;
  }
  for (size_t j = 0; j < offset; j++) {
    if (remaining_edges_after_this_[last_idx][j] == 0 &&
        (is_all_pair_ || (bag_local_idx_map_[last_idx][terminals_[0]] != j &&
                          bag_local_idx_map_[last_idx][terminals_[1]] != j))) {
      lab[found] = j;
      ptn[found++] = 0;
    }
  }
  assert(found == offset);
  for (size_t j = offset; j < sg.nv; j++) {
    lab[j] = j;
    ptn[j] = 1;
  }
  ptn[offset - 1] = 0;
  // remove the edges that are incident to two_edge vertices
  for (auto old_v_idx : two_edge) {
    auto actual = old_v_idx + offset;
    // for all neighbors remove old_v_idx
    for (size_t i = 0; i < sg.d[actual]; i++) {
      auto neigh = sg.e[sg.v[actual] + i];
      for (size_t j = 0; j < sg.d[neigh]; j++) {
        if (sg.e[sg.v[neigh] + j] == actual) {
          std::swap(sg.e[sg.v[neigh] + j], sg.e[sg.v[neigh] + sg.d[neigh] - 1]);
          sg.e[sg.v[neigh] + --sg.d[neigh]] = 0;
          break;
        }
      }
      sg.e[sg.v[actual] + i] = 0;
    }
    sg.nde -= 2 * sg.d[actual];
    sg.d[actual] = 0;
    assert(sg.d[old_v_idx] == 0);
    assert(sg.e[sg.v[old_v_idx]] == 0);
    sg.e[sg.v[old_v_idx]] = actual;
    sg.d[old_v_idx]++;
    assert(sg.e[sg.v[actual]] == 0);
    sg.e[sg.v[actual]] = old_v_idx;
    sg.d[actual]++;
    sg.nde += 2;
  }
  // std::cout << sg.nv << " " << sg.nde << std::endl;
  // std::cout << sg.vlen << " " << sg.dlen << " " << sg.elen << std::endl;
  // size_t e_count = 0;
  // for(size_t v = 0; v < sg.nv; v++) {
  //     std::vector<int> neighs(sg.nv, 0);
  //     for(size_t i = sg.v[v]; i < sg.v[v] + sg.d[v]; i++) {
  //         neighs[sg.e[i]] = 1;
  //         e_count++;
  //     }
  //     for(size_t i = 0; i < sg.nv; i++) {
  //         std::cout << neighs[i];
  //     }
  //     std::cout << std::endl;
  // }
  // for(size_t v = 0; v < sg.nv; v++) {
  //     std::cout << lab[v] << " ";
  // }
  // std::cout << std::endl;
  // for(size_t v = 0; v < sg.nv; v++) {
  //     std::cout << ptn[v] << " ";
  // }
  // std::cout << std::endl;
  // assert(e_count == sg.nde);
  DEFAULTOPTIONS_SPARSEGRAPH(options);
  options.getcanon = true;
  options.defaultptn = false;
  // options.schreier = true;
  // options.tc_level = 1000;
  statsblk stats;
  SG_DECL(canon_sg);
  canon_sg.v = (edge_t *)malloc(sizeof(edge_t) * (sg.nv + sg.nde) +
                                sizeof(degree_t) * sg.nv);
  canon_sg.d = (degree_t *)(canon_sg.v + sg.nv);
  canon_sg.e = (edge_t *)(canon_sg.d + sg.nv);
  canon_sg.nv = sg.nv;
  canon_sg.nde = sg.nde;
  canon_sg.vlen = sg.nv;
  canon_sg.dlen = sg.nv;
  canon_sg.elen = sg.nde;
  sparsenauty(&sg, lab, ptn, orbits, &options, &stats, &canon_sg);
  sortlists_sg(&canon_sg);
  free(sg.v);
  free(lab);
  free(ptn);
  free(orbits);
  return canon_sg;
}

template <template <typename> typename Count_structure, typename count_t>
void NautyPathwidthSearch<Count_structure, count_t>::print_stats() const {
  size_t pos_hits = 0, neg_hits = 0;
  for (size_t i = 0; i < nthreads_; i++) {
    pos_hits += pos_hits_[i];
    neg_hits += neg_hits_[i];
  }
  size_t edges = 0, propagations = 0;
  for (size_t i = 0; i < nthreads_; i++) {
    edges += edges_[i];
    propagations += propagations_[i];
  }
  LOG << "Cache hit rate: " << 100 * pos_hits / (double)(pos_hits + neg_hits)
      << "% (" << pos_hits << "/" << pos_hits + neg_hits << ")" << std::endl;
  LOG << "#Edges: " << edges << " #Propagations: " << propagations << std::endl;
}

} // namespace fpc
