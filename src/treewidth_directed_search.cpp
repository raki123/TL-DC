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

#include "treewidth_directed_search.h"
#include "logging.h"
#include <algorithm>

namespace fpc {

namespace {

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
void prettyPrint(DirectedFrontier frontier) {
  for (auto idx : frontier) {
    if (idx == no_edge_index) {
      LOG << "*"
          << " ";
    } else if (idx == outgoing_index) {
      LOG << "<"
          << " ";
    } else if (idx == incoming_index) {
      LOG << ">"
          << " ";
    } else if (idx == two_edge_index) {
      LOG << "#"
          << " ";
    } else if (idx & incoming_flag) {
      LOG << "-" << std::size_t(idx ^ incoming_flag) << " ";
    } else {
      LOG << std::size_t(idx) << " ";
    }
  }
  LOG << std::endl;
}
#pragma GCC diagnostic pop

} // namespace

template <template <typename> typename Count_structure, typename count_t>
TreewidthDirectedSearch<Count_structure, count_t>::TreewidthDirectedSearch(
    Digraph const &input, AnnotatedDecomposition decomposition, size_t nthreads)
    : nthreads_(nthreads), graph_(input), max_length_(graph_.max_length()),
      terminals_(graph_.terminals()), decomposition_(decomposition),
      remaining_in_arcs_after_this_(decomposition_.size()),
      remaining_out_arcs_after_this_(decomposition_.size()),
      remaining_vertices_after_this_(decomposition_.size()),
      eliminated_after_this_(decomposition_.size()),
      bag_local_distance_(decomposition_.size()),
      bag_lock_(decomposition_.size()), thread_local_result_(nthreads_, 0),
      mod_val_(nthreads_), mod_lock_(mod_val_), cur_cache_(mod_val_),
      next_cache_(mod_val_), join_cache_(decomposition_.size()),
      pos_hits_(nthreads_, 0), neg_hits_(nthreads_, 0), edges_(nthreads_, 0),
      propagations_(nthreads_, 0), merges_(nthreads_, 0),
      unsuccessful_merges_(nthreads_, 0) {

  omp_set_num_threads(nthreads_);
  for (auto i = 0; i < mod_val_; i++) {
    omp_init_lock(&mod_lock_[i]);
  }
  decomposition_.foreach_post_order([&](auto bag_idx) {
    if (decomposition_[bag_idx].type == JOIN) {
      join_cache_[bag_idx].resize(mod_val_);
    }
  });

  max_width_ = decomposition_.width();

  decomposition_.foreach_post_order([this](auto bag_idx) {
    auto &bag = decomposition_[bag_idx].bag;
    if (std::find(bag.begin(), bag.end(), terminals_[0]) == bag.end()) {
      bag.push_back(terminals_[0]);
    }
    if (std::find(bag.begin(), bag.end(), terminals_[1]) == bag.end()) {
      bag.push_back(terminals_[1]);
    }
    max_width_ = std::max(max_width_, bag.size());
  });

  if (max_width_ >= 0b0'1111110) {
    throw std::runtime_error("width of decomposition too high");
  }

  decomposition_.foreach_post_order([&](auto bag_idx) {
    assert(bag_idx == arc_weights_.size());
    auto type = decomposition_[bag_idx].type;
    if (type != LEAF && type != PATH_LIKE) {
      arc_weights_.push_back(std::nullopt);
      return;
    }
    auto [from, to] = decomposition_[bag_idx].edge;
    if (auto val = graph_.arc_weight(from, to)) {
      arc_weights_.push_back(
          convert_to<count_t, Count_structure>(*val, max_length_));
    } else {
      arc_weights_.push_back(std::nullopt);
    }
  });

  auto invalid_index = two_edge_index;
  std::vector<directed_frontier_index_t> vertex_to_idx(graph_.nr_vertices(),
                                                       invalid_index);

  std::vector<std::set<directed_frontier_index_t>> free_slots(
      decomposition_.size());
  decomposition_.foreach_pre_order([&](auto bag_idx) {
    auto &bag = decomposition_[bag_idx].bag;
    if (decomposition_.get_root() == bag_idx) {
      for (directed_frontier_index_t idx{0}; idx < max_width_;
           idx = directed_frontier_index_t(idx + 1)) {
        free_slots[bag_idx].insert(idx);
      }
    } else {
      auto parent = decomposition_.parent(bag_idx);
      free_slots[bag_idx] = free_slots[parent];
      for (auto prev : decomposition_[parent].bag) {
        if (auto it = srg::find(bag, prev); it == bag.end()) {
          free_slots[bag_idx].insert(vertex_to_idx[prev]);
        }
      }
    }
    for (auto v : bag) {
      if (vertex_to_idx[v] == invalid_index) {
        assert(!free_slots[bag_idx].empty());
        auto idx = *free_slots[bag_idx].begin();
        free_slots[bag_idx].erase(free_slots[bag_idx].begin());
        vertex_to_idx[v] = idx;
      }
    }
  });

  terminals_idx_[0] = vertex_to_idx[terminals_[0]];
  terminals_idx_[1] = vertex_to_idx[terminals_[1]];

  std::set<Vertex> empty;
  std::vector<Vertex> all;
  for (Vertex i = 0; i < graph_.nr_vertices(); i++) {
    all.push_back(i);
  }
  Digraph cur_graph = graph_.subgraph(all);
  auto cur_arcs = [&](Vertex v) {
    return cur_graph.in_neighbors(v).size() + cur_graph.out_neighbors(v).size();
  };
  auto overall_arcs = [&](Vertex v) {
    return graph_.in_neighbors(v).size() + graph_.out_neighbors(v).size();
  };
  decomposition_.foreach_post_order([&](auto bag_idx) {
    omp_init_lock(&bag_lock_[bag_idx]);
    auto &node = decomposition_[bag_idx];
    auto &remaining_in = remaining_in_arcs_after_this_[bag_idx];
    auto &remaining_out = remaining_out_arcs_after_this_[bag_idx];
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
              cur_graph.remove_arc(top_node.edge.first, top_node.edge.second);
              break;
            case NodeType::JOIN:
              break;
            }
          },
          bag_idx);
      break;
    }
    for (size_t i = 0; i < node.bag.size(); i++) {
      if (node.bag[i] == terminals_[0] || node.bag[i] == terminals_[1]) {
        continue;
      }
      if (node.type != NodeType::JOIN && cur_arcs(node.bag[i]) == 0) {
        std::swap(node.bag[i], node.bag.back());
        node.bag.pop_back();
        i--;
      } else if (cur_arcs(node.bag[i]) == overall_arcs(node.bag[i]) &&
                 node.bag[i] != node.edge.first &&
                 node.bag[i] != node.edge.second) {
        std::swap(node.bag[i], node.bag.back());
        node.bag.pop_back();
        i--;
      } else if (node.type == NodeType::JOIN) {
        auto v = node.bag[i];

        auto left_child = node.children.first;
        auto it = srg::find(decomposition_[left_child].bag, v);
        auto in_left = it != decomposition_[left_child].bag.end();

        auto right_child = node.children.second;
        it = srg::find(decomposition_[right_child].bag, v);
        auto in_right = it != decomposition_[right_child].bag.end();

        auto idx = vertex_to_idx[v];
        if (in_right &&
            remaining_in_arcs_after_this_[right_child][idx] +
                    remaining_out_arcs_after_this_[right_child][idx] ==
                0) {
          assert(!in_left);
          std::swap(node.bag[i], node.bag.back());
          node.bag.pop_back();
          i--;
        } else if (in_left &&
                   remaining_in_arcs_after_this_[left_child][idx] +
                           remaining_out_arcs_after_this_[left_child][idx] ==
                       0) {
          assert(!in_right);
          std::swap(node.bag[i], node.bag.back());
          node.bag.pop_back();
          i--;
        }
      }
    }
    bag_local_distance_[bag_idx].resize(max_width_);
    remaining_in.resize(max_width_);
    remaining_out.resize(max_width_);
    if (node.type != NodeType::JOIN) {
      cur_graph.remove_arc(node.edge.first, node.edge.second);
    }
    for (auto v : node.bag) {
      assert(vertex_to_idx[v] != invalid_index);
      remaining_in[vertex_to_idx[v]] = cur_graph.in_neighbors(v).size();
      remaining_out[vertex_to_idx[v]] = cur_graph.out_neighbors(v).size();
    }
    remaining_vertices_after_this_[bag_idx] = cur_graph.nr_vertices();
    for (auto v : node.bag) {
      std::vector<Edge_length> distance(graph_.nr_vertices(),
                                        invalid_distance_);
      cur_graph.dijkstra(v, distance, true);
      std::vector<Edge_length> local_distance(max_width_, invalid_distance_);
      for (auto other : node.bag) {
        if (other != v) {
          local_distance[vertex_to_idx[other]] = distance[other];
        }
      }
      bag_local_distance_[bag_idx][vertex_to_idx[v]] = local_distance;
    }
  });
  std::vector<char> seen(graph_.nr_vertices(), false);
  decomposition_.foreach_pre_order([&](auto bag_idx) {
    auto &bag = decomposition_[bag_idx].bag;
    for (auto v : bag) {
      if (!seen[v]) {
        eliminated_after_this_[bag_idx].push_back(vertex_to_idx[v]);
        seen[v] = true;
      }
    }
  });

  for (Vertex v = 0; v < graph_.nr_vertices(); v++) {
    assert(cur_arcs(v) == 0);
  }

  decomposition_.foreach_post_order([&](auto bag_idx) {
    auto &node = decomposition_[bag_idx];
    if (node.type != NodeType::JOIN) {
      node.edge = {vertex_to_idx[node.edge.first],
                   vertex_to_idx[node.edge.second]};
    }
  });
}

template <template <typename> typename Count_structure, typename count_t>
mpz_class TreewidthDirectedSearch<Count_structure, count_t>::search() {
  for (size_t bag_idx = 0; bag_idx < decomposition_.size(); bag_idx++) {

    if (decomposition_[bag_idx].type == LEAF) {
      DirectedFrontier initial_frontier(max_width_, no_edge_index);
      initial_frontier[terminals_idx_[0]] = incoming_index;
      initial_frontier[terminals_idx_[1]] = outgoing_index;
      initial_frontier.rehash();
      auto initial_result = Partial_result(max_length_);
      cur_cache_[0].first.emplace(initial_frontier, initial_result);
    } else if (decomposition_[bag_idx].type == JOIN) {
#pragma omp parallel for default(shared)
      for (size_t mod = 0; mod < mod_val_; mod++) {
        std::swap(join_cache_[bag_idx][mod], cur_cache_[mod]);
      }
    }

    LOG << "\r" << bag_idx << " / " << decomposition_.size() - 1 << " of type ";
    if (decomposition_[bag_idx].type == LEAF) {
      auto entries = 0;
      for (auto mod = 0; mod < mod_val_; mod++) {
        entries += cur_cache_[mod].first.size();
      }
      LOG << "leaf, with " << entries << " entries.";
    } else if (decomposition_[bag_idx].type == PATH_LIKE) {
      auto entries = 0;
      for (auto mod = 0; mod < mod_val_; mod++) {
        entries += cur_cache_[mod].first.size();
      }
      LOG << "path, with " << entries << " entries.";
    } else if (decomposition_[bag_idx].type == JOIN) {
      auto entries_first = 0;
      auto entries_second = 0;
      for (auto mod = 0; mod < mod_val_; mod++) {
        entries_first += cur_cache_[mod].first.size();
        entries_second += cur_cache_[mod].second.size();
      }
      LOG << "join, with (" << entries_first << "," << entries_second
          << ") entries.";
    }
    if (decomposition_[bag_idx].type == LEAF ||
        decomposition_[bag_idx].type == PATH_LIKE) {
// PATH_LIKE/LEAF
#pragma omp parallel for default(shared)
      for (size_t mod = 0; mod < mod_val_; mod++) {
        size_t thread_id = omp_get_thread_num();
        for (auto &[orig_frontier, orig_result] : cur_cache_[mod].first) {
          auto copy_frontier = orig_frontier;
          assert(copy_frontier.size() > decomposition_[bag_idx].edge.first);
          assert(copy_frontier.size() > decomposition_[bag_idx].edge.second);
          auto copy_result = orig_result;
          size_t new_idx = decomposition_[bag_idx].parent;
          assert(new_idx != size_t(-1));
          if (canTake(copy_frontier, bag_idx, copy_result)) {
            if (arc_weights_[bag_idx]) {
              copy_result *= *arc_weights_[bag_idx];
            } else {
              copy_result.increment_offset();
            }
            advance(copy_frontier, bag_idx);
            copy_frontier.rehash();
            cache(copy_frontier, new_idx, copy_result, thread_id);
            copy_frontier = orig_frontier;
            copy_result = orig_result;
          }
          if (canSkip(copy_frontier, bag_idx, copy_result)) {
            advance(copy_frontier, bag_idx);
            copy_frontier.rehash();
            cache(copy_frontier, new_idx, copy_result, thread_id);
          }
        }
      }
    } else {
      // JOIN
      // FIXME: check if its better to take either bag as the outer one
#pragma omp parallel for default(shared)
      for (size_t mod = 0; mod < mod_val_; mod++) {
        size_t thread_id = omp_get_thread_num();
        for (auto &[orig_frontier, orig_result] : cur_cache_[mod].first) {
          // FIXME: can we just modify the object somehow?
          for (size_t other_mod = 0; other_mod < mod_val_; other_mod++) {
            for (auto const &[right_frontier, right_result] :
                 cur_cache_[other_mod].second) {
              auto left_frontier = orig_frontier;
              auto left_result = orig_result;
              merges_[thread_id]++;
              if (!merge(left_frontier, right_frontier, bag_idx, left_result,
                         right_result)) {
                unsuccessful_merges_[thread_id]++;
                continue;
              }
              size_t new_idx = decomposition_[bag_idx].parent;
              assert(new_idx != size_t(-1));
              left_frontier.rehash();
              cache(left_frontier, new_idx, left_result, thread_id);
            }
          }
        }
      }
    }

#pragma omp parallel for default(shared)
    for (size_t mod = 0; mod < mod_val_; mod++) {
      cur_cache_[mod] = {};
      size_t new_idx = decomposition_[bag_idx].parent;
      if (new_idx >= decomposition_.size()) {
        continue;
      }
      if (decomposition_[new_idx].type == JOIN &&
          bag_idx == decomposition_[new_idx].children.first) {
        std::swap(join_cache_[new_idx][mod].first, next_cache_[mod]);
      } else if (decomposition_[new_idx].type == JOIN &&
                 bag_idx == decomposition_[new_idx].children.second) {
        std::swap(join_cache_[new_idx][mod].second, next_cache_[mod]);
      } else {
        assert(decomposition_[new_idx].type == PATH_LIKE);
        std::swap(cur_cache_[mod].first, next_cache_[mod]);
        next_cache_[mod].reserve(cur_cache_[mod].first.size() * 2);
      }
    }
  }

  LOG << std::endl;
  print_stats();

  mpz_class rv = 0;
  for (size_t id = 0; id < nthreads_; id++) {
    rv += to_mpz(thread_local_result_[id]);
  }
  return rv;
}

template <template <typename> typename Count_structure, typename count_t>
void TreewidthDirectedSearch<Count_structure, count_t>::cache(
    DirectedFrontier &frontier, size_t new_idx, Partial_result &partial_results,
    size_t thread_id) {
  edges_[thread_id]++;
  // put it into the cache
  auto mod = DirectedFrontierHash{}(frontier) % mod_val_;
  if (new_idx < decomposition_.size()) {
    omp_set_lock(&mod_lock_[mod]);
    auto ins = next_cache_[mod].try_emplace(frontier, partial_results);
    if (!ins.second) {
      pos_hits_[thread_id]++;
      // there is already an element with that key
      // instead increase the partial result for that key
      ins.first->second += partial_results;
    } else {
      neg_hits_[thread_id]++;
    }
    omp_unset_lock(&mod_lock_[mod]);
  }
}

template <template <typename> typename Count_structure, typename count_t>
void TreewidthDirectedSearch<Count_structure, count_t>::includeSolutions(
    DirectedFrontier const &frontier, Partial_result const &partial_result,
    size_t bag_idx) {
  assert(partial_result.offset() < max_length_);

  if (arc_weights_[bag_idx] &&
      partial_result.offset() + arc_weights_[bag_idx]->offset() > max_length_) {
    // no paths short enough
    return;
  }

  if (srg::any_of(frontier, [this](auto v) {
        return v != no_edge_index && v != two_edge_index;
      })) {
    return;
  }

  size_t thread_id = omp_get_thread_num();
  if (!arc_weights_[bag_idx]) {
    thread_local_result_[thread_id] +=
        partial_result.total_count(max_length_ - 1);
  } else {
    auto cp = partial_result;
    cp *= *arc_weights_[bag_idx];
    thread_local_result_[thread_id] += cp.total_count(max_length_);
  }
}

template <template <typename> typename Count_structure, typename count_t>
bool TreewidthDirectedSearch<Count_structure, count_t>::canTake(
    DirectedFrontier &frontier, size_t bag_idx,
    Partial_result &partial_result) {
  assert(partial_result.offset() + 1 <= max_length_);
  Edge edge = decomposition_[bag_idx].edge;
  auto v_idx = directed_frontier_index_t(edge.first) | incoming_flag;
  auto w_idx = directed_frontier_index_t(edge.second);
  assert(w_idx < frontier.size());
  assert((v_idx ^ incoming_flag) < frontier.size());

  // must not take more than two edges
  if (frontier[v_idx] == two_edge_index || frontier[w_idx] == two_edge_index) {
    return false;
  }

  // two incoming arcs
  if (frontier[w_idx] & incoming_flag) {
    return false;
  }

  // two outgoing arcs
  if (frontier[v_idx] != no_edge_index &&
      (frontier[v_idx] & incoming_flag) == 0) {
    return false;
  }

  // if there are already partial paths, they are now in the correct direction

  // must not be closing a path to a loop
  if (frontier[v_idx] == (w_idx ^ incoming_flag)) {
    assert(frontier[w_idx] == (v_idx ^ incoming_flag));
    return false;
  }
  assert(frontier[w_idx] != (v_idx ^ incoming_flag));
  assert(frontier[w_idx] != incoming_index);
  assert(frontier[v_idx] != outgoing_index);

  // can connect two invalid paths but cannot continue after
  if (frontier[v_idx] == incoming_index && frontier[w_idx] == outgoing_index) {
    frontier[v_idx] = two_edge_index;
    frontier[w_idx] = two_edge_index;
    includeSolutions(frontier, partial_result, bag_idx);
    frontier[v_idx] = incoming_index;
    frontier[w_idx] = outgoing_index;
    return false;
  }

  auto v_remaining =
      remaining_in_arcs_after_this_[bag_idx][v_idx ^ incoming_flag];
  auto w_remaining = remaining_out_arcs_after_this_[bag_idx][w_idx];

  // can connect two invalid paths but cannot continue after
  if (!v_remaining && frontier[v_idx] == no_edge_index &&
      frontier[w_idx] == outgoing_index) {
    return false;
  }
  // can connect two invalid paths but cannot continue after
  if (!w_remaining && frontier[w_idx] == no_edge_index &&
      frontier[v_idx] == incoming_index) {
    return false;
  }
  // can connect two invalid paths but cannot continue after
  if (!v_remaining && !w_remaining && frontier[v_idx] == no_edge_index &&
      frontier[w_idx] == no_edge_index) {
    return false;
  }

  // length based pruning
  auto old_v = frontier[v_idx];
  auto old_w = frontier[w_idx];
  std::vector<std::pair<directed_frontier_index_t, directed_frontier_index_t>>
      restore({std::make_pair(v_idx, old_v), std::make_pair(w_idx, old_w)});

  if (old_v == no_edge_index) {
    if (old_w == no_edge_index) {
      frontier[v_idx] = w_idx;
      frontier[w_idx] = v_idx;
    } else if (old_w == outgoing_index) {
      frontier[w_idx] = two_edge_index;
      frontier[v_idx] = outgoing_index;
    } else if (old_w == incoming_index) {
      assert(false);
    } else if (old_w == two_edge_index) {
      assert(false);
    } else {
      assert((old_w & incoming_flag) == 0);
      restore.push_back(std::make_pair(old_w, frontier[old_w]));
      frontier[old_w] = v_idx;
      frontier[w_idx] = two_edge_index;
      frontier[v_idx] = old_w;
    }
  } else if (old_v == outgoing_index) {
    assert(false);
  } else if (old_v == incoming_index) {
    if (old_w == no_edge_index) {
      frontier[v_idx] = two_edge_index;
      frontier[w_idx] = incoming_index;
    } else if (old_w == outgoing_index) {
      assert(false);
    } else if (old_w == incoming_index) {
      assert(false);
    } else if (old_w == two_edge_index) {
      assert(false);
    } else {
      restore.push_back(std::make_pair(old_w, frontier[old_w]));
      frontier[old_w] = incoming_index;
      frontier[w_idx] = two_edge_index;
      frontier[v_idx] = two_edge_index;
    }
  } else if (old_v == two_edge_index) {
    assert(false);
  } else {
    assert(old_v & incoming_flag);
    restore.push_back(std::make_pair(old_v, frontier[old_v]));
    if (old_w == no_edge_index) {
      frontier[old_v] = w_idx;
      frontier[v_idx] = two_edge_index;
      frontier[w_idx] = old_v;
    } else if (old_w == outgoing_index) {
      frontier[old_v] = outgoing_index;
      frontier[v_idx] = two_edge_index;
      frontier[w_idx] = two_edge_index;
    } else if (old_w == incoming_index) {
      assert(false);
    } else if (old_w == two_edge_index) {
      assert(false);
    } else {
      restore.push_back(std::make_pair(old_w, frontier[old_w]));
      frontier[old_w] = old_v;
      frontier[old_v] = old_w;
      frontier[w_idx] = two_edge_index;
      frontier[v_idx] = two_edge_index;
    }
  }
  if (v_remaining == 0) {
    assert(frontier[v_idx] != incoming_index);
    assert(frontier[v_idx] != outgoing_index);
    assert(frontier[v_idx] != no_edge_index);
    if (frontier[v_idx] != two_edge_index) {
      frontier[frontier[v_idx]] = incoming_index;
      frontier[v_idx] = two_edge_index;
    }
  }
  if (w_remaining == 0) {
    assert(frontier[w_idx] != incoming_index);
    assert(frontier[w_idx] != outgoing_index);
    assert(frontier[w_idx] != no_edge_index);
    if (frontier[w_idx] != two_edge_index) {
      frontier[frontier[w_idx]] = outgoing_index;
      frontier[w_idx] = two_edge_index;
    }
  }
  std::vector<directed_frontier_index_t> paths;
  std::vector<directed_frontier_index_t> cut_paths;
  for (directed_frontier_index_t idx{0}; idx < frontier.size();
       idx = directed_frontier_index_t(idx + 1)) {
    if (frontier[idx] != no_edge_index && frontier[idx] != two_edge_index) {
      if (frontier[idx] == outgoing_index || frontier[idx] == incoming_index) {
        cut_paths.push_back(idx);
      } else {
        paths.push_back(idx);
      }
    }
  }

  auto revert = [&]() {
    for (auto [idx, rest] : restore) {
      frontier[idx] = rest;
    }
  };

  if (cut_paths.size() > 2) {
    revert();
    return false;
  }
  assert(paths.size() % 2 == 0);

  Edge_length arc_length =
      arc_weights_[bag_idx] ? arc_weights_[bag_idx]->offset() : 1;
  // + 1 - 1 since were taking the edge
  if (paths.size() / 2 + cut_paths.size() > 1) {
    if (paths.size() / 2 + cut_paths.size() + partial_result.offset() +
            arc_length - 1 >
        max_length_) {
      revert();
      return false;
    }
  } else if (paths.size() / 2 + cut_paths.size() + partial_result.offset() +
                 arc_length >
             max_length_) {
    revert();
    return false;
  }

  auto rv = distancePrune(frontier, paths, cut_paths, bag_idx,
                          partial_result.offset() + arc_length);
  if (!rv) {
    revert();
  }
  return rv;
}

template <template <typename> typename Count_structure, typename count_t>
bool TreewidthDirectedSearch<Count_structure, count_t>::canSkip(
    DirectedFrontier &frontier, size_t bag_idx,
    Partial_result const &partial_result) {
  assert(partial_result.offset() + 1 <= max_length_);
  Edge edge = decomposition_[bag_idx].edge;
  auto v_idx = directed_frontier_index_t(edge.first);
  auto w_idx = directed_frontier_index_t(edge.second);

  auto v_remaining = remaining_out_arcs_after_this_[bag_idx][v_idx];
  auto w_remaining = remaining_in_arcs_after_this_[bag_idx][w_idx];

  // cannot let an outside path die
  if ((!v_remaining && frontier[v_idx] == incoming_index) ||
      (!w_remaining && frontier[w_idx] == outgoing_index)) {
    return false;
  }
  // both ends will be outside the outside path dies
  if (!v_remaining && !w_remaining &&
      frontier[v_idx] == (w_idx | incoming_flag)) {
    return false;
  }
  // length based pruning
  std::vector<directed_frontier_index_t> paths;
  std::vector<directed_frontier_index_t> cut_paths;
  for (directed_frontier_index_t idx{0}; idx < frontier.size();
       idx = directed_frontier_index_t(idx + 1)) {
    if (frontier[idx] != no_edge_index && frontier[idx] != two_edge_index) {
      if (frontier[idx] == incoming_index) {
        cut_paths.push_back(idx);
      } else if (frontier[idx] == outgoing_index) {
        cut_paths.push_back(idx);
      } else {
        paths.push_back(idx);
      }
    }
  }
  assert(cut_paths.size() == 2);
  assert(paths.size() % 2 == 0);
  bool last_and_incomplete_v =
      !v_remaining && frontier[v_idx] != two_edge_index &&
      frontier[v_idx] != no_edge_index && (frontier[v_idx] & incoming_flag);
  bool last_and_incomplete_w = !w_remaining &&
                               frontier[w_idx] != two_edge_index &&
                               frontier[w_idx] != no_edge_index &&
                               ((frontier[w_idx] & incoming_flag) == 0);
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
bool TreewidthDirectedSearch<Count_structure, count_t>::distancePrune(
    DirectedFrontier const &frontier,
    std::vector<directed_frontier_index_t> const &paths,
    std::vector<directed_frontier_index_t> &cut_paths, size_t bag_idx,
    size_t offset) {
  // early exit, does not make sense to check here
  if (offset + remaining_vertices_after_this_[bag_idx] <= max_length_) {
    return true;
  }
  // advanced length based pruning

  // make sure that the outside path that is incoming is first
  if (cut_paths.size() == 2) {
    if (frontier[cut_paths.back()] & incoming_flag) {
      std::swap(cut_paths.front(), cut_paths.back());
    }
    assert(frontier[cut_paths.front()] & incoming_flag);
    assert(!(frontier[cut_paths.back()] & incoming_flag));
  }
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
    if (frontier[v] & incoming_flag) {
      for (auto other : paths) {
        if (frontier[other] & incoming_flag) {
          continue;
        }
        cur_min_dist = std::min(cur_min_dist, distance[v][other]);
      }
    } else {
      for (auto other : paths) {
        if (frontier[other] & incoming_flag) {
          cur_min_dist = std::min(cur_min_dist, distance[other][v]);
        }
      }
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
      if ((frontier[v] & incoming_flag) &&
          !(frontier[cut_end] & incoming_flag)) {
        cur_min_dist = std::min(cur_min_dist, distance[v][cut_end]);
      } else if (!(frontier[v] & incoming_flag) &&
                 (frontier[cut_end] & incoming_flag)) {
        cur_min_dist = std::min(cur_min_dist, distance[cut_end][v]);
      }
    }
    for (auto other : paths) {
      if (frontier[v] != other) {
        if ((frontier[v] & incoming_flag) &&
            !(frontier[other] & incoming_flag)) {
          cur_min_dist = std::min(cur_min_dist, distance[v][other]);
        } else if (!(frontier[v] & incoming_flag) &&
                   (frontier[other] & incoming_flag)) {
          cur_min_dist = std::min(cur_min_dist, distance[other][v]);
        }
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
bool TreewidthDirectedSearch<Count_structure, count_t>::merge(
    DirectedFrontier &left, DirectedFrontier const &right, size_t bag_idx,
    Partial_result &left_result, Partial_result const &right_result) {
  if (left_result.offset() + right_result.offset() > max_length_) {
    return false;
  }
  bool left_empty = left_result.offset() == 0;
  bool right_empty = right_result.offset() == 0;
  // merge the frontiers
  bool found_solution = false;
  for (directed_frontier_index_t idx{0}; idx < right.size();
       idx = directed_frontier_index_t(idx + 1)) {
    // start goal case
    if ((terminals_idx_[0] == idx || terminals_idx_[1] == idx)) {
      assert(left[idx] == incoming_index || left[idx] == outgoing_index ||
             left[idx] == two_edge_index);
      assert(right[idx] == incoming_index || right[idx] == outgoing_index ||
             right[idx] == two_edge_index);
      if (left[idx] == incoming_index) {
        left[idx] = right[idx];
      } else if (left[idx] == outgoing_index) {
        left[idx] = right[idx];
      } else if (right[idx] == two_edge_index) {
        return false;
      }
      continue;
    }
    // rest
    if (right[idx] == no_edge_index) {
      continue;
    } else if (right[idx] == two_edge_index) {
      if (left[idx] != no_edge_index) {
        return false;
      }
      left[idx] = two_edge_index;
    } else if (right[idx] == incoming_index) {
      if (left[idx] == no_edge_index) {
        left[idx] = incoming_index;
      } else if (left[idx] == two_edge_index) {
        return false;
      } else if (left[idx] == incoming_index) {
        return false;
      } else if (left[idx] == outgoing_index) {
        if (found_solution) {
          return false;
        }
        left[idx] = two_edge_index;
        found_solution = true;
      } else if (left[idx] & incoming_flag) {
        return false;
      } else {
        assert(left[left[idx]] == (idx ^ incoming_flag));
        left[left[idx]] = incoming_index;
        left[idx] = two_edge_index;
      }
    } else if (right[idx] == outgoing_index) {
      if (left[idx] == no_edge_index) {
        left[idx] = outgoing_index;
      } else if (left[idx] == two_edge_index) {
        return false;
      } else if (left[idx] == incoming_index) {
        if (found_solution) {
          return false;
        }
        left[idx] = two_edge_index;
        found_solution = true;
      } else if (left[idx] == outgoing_index) {
        return false;
      } else if (left[idx] & incoming_flag) {
        assert(left[left[idx]] == idx);
        left[left[idx]] = outgoing_index;
        left[idx] = two_edge_index;
      } else {
        return false;
      }
    } else if ((right[idx] & incoming_flag) &&
               (right[idx] ^ incoming_flag) > idx) {
      assert(right[right[idx]] == idx);
      if (left[idx] == no_edge_index) {
        if (left[right[idx]] == no_edge_index) {
          left[idx] = right[idx];
          left[right[idx]] = idx;
        } else if (left[right[idx]] == two_edge_index) {
          return false;
        } else if (left[right[idx]] == incoming_index) {
          left[idx] = incoming_index;
          left[right[idx]] = two_edge_index;
        } else if (left[right[idx]] == outgoing_index) {
          return false;
        } else if (left[right[idx]] & incoming_flag) {
          left[idx] = left[right[idx]];
          left[left[right[idx]]] = idx;
          left[right[idx]] = two_edge_index;
        } else {
          return false;
        }
      } else if (left[idx] == two_edge_index) {
        return false;
      } else if (left[idx] == incoming_index) {
        return false;
      } else if (left[idx] == outgoing_index) {
        if (left[right[idx]] == no_edge_index) {
          left[idx] = two_edge_index;
          left[right[idx]] = outgoing_index;
        } else if (left[right[idx]] == two_edge_index) {
          return false;
        } else if (left[right[idx]] == incoming_index) {
          if (found_solution) {
            return false;
          }
          left[idx] = two_edge_index;
          left[right[idx]] = two_edge_index;
          found_solution = true;
        } else if (left[right[idx]] == outgoing_index) {
          return false;
        } else if (left[right[idx]] & incoming_flag) {
          left[idx] = two_edge_index;
          left[left[right[idx]]] = outgoing_index;
          left[right[idx]] = two_edge_index;
        } else {
          return false;
        }
      } else if (left[idx] & incoming_flag) {
        return false;
      } else {
        if (left[right[idx]] == no_edge_index) {
          left[left[idx]] = right[idx];
          left[right[idx]] = left[idx];
          left[idx] = two_edge_index;
        } else if (left[right[idx]] == two_edge_index) {
          return false;
        } else if (left[right[idx]] == incoming_index) {
          left[left[idx]] = incoming_index;
          left[idx] = two_edge_index;
          left[right[idx]] = two_edge_index;
        } else if (left[right[idx]] == outgoing_index) {
          return false;
        } else if (left[right[idx]] & incoming_flag) {
          if ((left[right[idx]] ^ incoming_flag) == idx) {
            assert((left[idx] ^ incoming_flag) == right[idx]);
            assert(right[left[idx]] == idx);
            // found a loop
            return false;
          }
          left[left[idx]] = left[right[idx]];
          left[left[right[idx]]] = left[idx];
          left[idx] = two_edge_index;
          left[right[idx]] = two_edge_index;
        } else {
          return false;
        }
      }
    } else if (!(right[idx] & incoming_flag) && right[idx] > idx) {
      assert((right[right[idx]] ^ incoming_flag) == idx);
      if (left[idx] == no_edge_index) {
        if (left[right[idx]] == no_edge_index) {
          left[idx] = right[idx];
          left[right[idx]] = (idx ^ incoming_flag);
        } else if (left[right[idx]] == two_edge_index) {
          return false;
        } else if (left[right[idx]] == incoming_index) {
          return false;
        } else if (left[right[idx]] == outgoing_index) {
          left[idx] = outgoing_index;
          left[right[idx]] = two_edge_index;
        } else if (left[right[idx]] & incoming_flag) {
          return false;
        } else {
          left[idx] = left[right[idx]];
          left[left[right[idx]]] = (idx ^ incoming_flag);
          left[right[idx]] = two_edge_index;
        }
      } else if (left[idx] == two_edge_index) {
        return false;
      } else if (left[idx] == incoming_index) {
        if (left[right[idx]] == no_edge_index) {
          left[idx] = two_edge_index;
          left[right[idx]] = incoming_index;
        } else if (left[right[idx]] == two_edge_index) {
          return false;
        } else if (left[right[idx]] == incoming_index) {
          return false;
        } else if (left[right[idx]] == outgoing_index) {
          if (found_solution) {
            return false;
          }
          left[idx] = two_edge_index;
          left[right[idx]] = two_edge_index;
          found_solution = true;
        } else if (left[right[idx]] & incoming_flag) {
          return false;
        } else {
          left[idx] = two_edge_index;
          left[left[right[idx]]] = incoming_index;
          left[right[idx]] = two_edge_index;
        }
      } else if (left[idx] == outgoing_index) {
        return false;
      } else if (left[idx] & incoming_flag) {
        if (left[right[idx]] == no_edge_index) {
          left[left[idx]] = right[idx];
          left[right[idx]] = left[idx];
          left[idx] = two_edge_index;
        } else if (left[right[idx]] == two_edge_index) {
          return false;
        } else if (left[right[idx]] == incoming_index) {
          return false;
        } else if (left[right[idx]] == outgoing_index) {
          left[left[idx]] = outgoing_index;
          left[idx] = two_edge_index;
          left[right[idx]] = two_edge_index;
        } else if (left[right[idx]] & incoming_flag) {
          return false;
        } else {
          if (left[right[idx]] == idx) {
            assert((left[idx] ^ incoming_flag) == right[idx]);
            assert((right[left[idx]] ^ incoming_flag) == idx);
            // found a loop
            return false;
          }
          left[left[idx]] = left[right[idx]];
          left[left[right[idx]]] = left[idx];
          left[idx] = two_edge_index;
          left[right[idx]] = two_edge_index;
        }
      } else {
        return false;
      }
    }
  }
  assert(!left_empty || !right_empty || !found_solution);
  // mark out of scope edge ends and populate paths, cut_paths
  std::vector<directed_frontier_index_t> paths;
  std::vector<directed_frontier_index_t> cut_paths;
  for (directed_frontier_index_t idx{0}; idx < right.size();
       idx = directed_frontier_index_t(idx + 1)) {
    // start goal case
    if (terminals_idx_[0] == idx || terminals_idx_[1] == idx) {
      if (left[idx] == outgoing_index) {
        if (!remaining_in_arcs_after_this_[bag_idx][idx]) {
          return false;
        }
        cut_paths.push_back(idx);
      }
      if (left[idx] == incoming_index) {
        if (!remaining_out_arcs_after_this_[bag_idx][idx]) {
          return false;
        }
        cut_paths.push_back(idx);
      }
      continue;
    }
    // rest
    if (left[idx] == incoming_index) {
      if (!remaining_out_arcs_after_this_[bag_idx][idx]) {
        if (found_solution) {
          return false;
        }
        left[idx] = two_edge_index;
        found_solution = true;
      } else {
        cut_paths.push_back(idx);
      }
    } else if (left[idx] == outgoing_index) {
      if (!remaining_in_arcs_after_this_[bag_idx][idx]) {
        if (found_solution) {
          return false;
        }
        left[idx] = two_edge_index;
        found_solution = true;
      } else {
        cut_paths.push_back(idx);
      }
    } else if (left[idx] != no_edge_index && left[idx] != two_edge_index &&
               (left[idx] & incoming_flag) &&
               idx > (left[idx] ^ incoming_flag)) {
      assert(left[left[idx]] == idx);
      if (!remaining_out_arcs_after_this_[bag_idx][idx] &&
          !remaining_in_arcs_after_this_[bag_idx][left[idx] ^ incoming_flag]) {
        if (found_solution) {
          return false;
        }
        left[left[idx]] = two_edge_index;
        left[idx] = two_edge_index;
        found_solution = true;
      } else if (!remaining_out_arcs_after_this_[bag_idx][idx] &&
                 remaining_in_arcs_after_this_[bag_idx]
                                              [left[idx] ^ incoming_flag]) {
        left[left[idx]] = outgoing_index;
        cut_paths.push_back(left[idx] ^ incoming_flag);
        left[idx] = two_edge_index;
      } else if (remaining_out_arcs_after_this_[bag_idx][idx] &&
                 !remaining_in_arcs_after_this_[bag_idx]
                                               [left[idx] ^ incoming_flag]) {
        left[left[idx]] = two_edge_index;
        left[idx] = incoming_index;
        cut_paths.push_back(idx);
      } else {
        paths.push_back(idx);
        paths.push_back(left[idx] ^ incoming_flag);
      }
    } else if (left[idx] != no_edge_index && left[idx] != two_edge_index &&
               !(left[idx] & incoming_flag) && idx > left[idx]) {
      assert(left[left[idx]] == (idx ^ incoming_flag));
      if (!remaining_in_arcs_after_this_[bag_idx][idx] &&
          !remaining_out_arcs_after_this_[bag_idx][left[idx]]) {
        if (found_solution) {
          return false;
        }
        left[left[idx]] = two_edge_index;
        left[idx] = two_edge_index;
        found_solution = true;
      } else if (!remaining_in_arcs_after_this_[bag_idx][idx] &&
                 remaining_out_arcs_after_this_[bag_idx][left[idx]]) {
        left[left[idx]] = incoming_index;
        cut_paths.push_back(left[idx]);
        left[idx] = two_edge_index;
      } else if (remaining_in_arcs_after_this_[bag_idx][idx] &&
                 !remaining_out_arcs_after_this_[bag_idx][left[idx]]) {
        left[left[idx]] = two_edge_index;
        left[idx] = outgoing_index;
        cut_paths.push_back(idx);
      } else {
        paths.push_back(idx);
        paths.push_back(left[idx]);
      }
    }
  }

  assert(paths.size() % 2 == 0);

  // possibly include solutions
  size_t thread_id = omp_get_thread_num();
  if (found_solution) {
    if (left_empty || right_empty) {
      return false;
    }
    // we cannot continue this either way but if there are not other partial
    // paths we have to include the solution
    if (paths.size() + cut_paths.size() == 0) {
      // adapt result
      left_result *= right_result;
      thread_local_result_[thread_id] += left_result.total_count(max_length_);
    }
    return false;
  }
  if (!left_empty && !right_empty && paths.size() / 2 + cut_paths.size() == 1) {
    // if either is empty then we already counted the solutions
    // this way there is currently exactly one path and we have not seen it
    // yet
    // adapt result
    left_result *= right_result;
    thread_local_result_[thread_id] += left_result.total_count(max_length_);
  }

  if (cut_paths.size() > 2) {
    return false;
  }

  auto offset = left_result.offset() + right_result.offset();

  // - 1 since have to take one edge less
  if (paths.size() / 2 + cut_paths.size() > 1) {
    if (paths.size() / 2 + cut_paths.size() + offset > max_length_ + 1) {
      return false;
    }
  } else if (paths.size() / 2 + cut_paths.size() + offset > max_length_) {
    return false;
  }

  if (!distancePrune(left, paths, cut_paths, bag_idx, offset)) {
    return false;
  }
  // adapt result
  left_result *= right_result;
  advance(left, bag_idx);
  return true;
}

template <template <typename> typename Count_structure, typename count_t>
void TreewidthDirectedSearch<Count_structure, count_t>::advance(
    DirectedFrontier &frontier, size_t bag_idx) {
  int found_two = 0;
  int found_invalid = 0;
  int found_path = 0;
  for (auto v : frontier) {
    if (v == two_edge_index) {
      found_two += 1;
    } else if (v == outgoing_index) {
      found_invalid += 1;
      found_path += 2;
    } else if (v == incoming_index) {
      found_invalid += 1;
      found_path += 2;
    } else if (v != no_edge_index) {
      found_path += 1;
    }
  }

  assert(found_path % 2 == 0);
  assert(!found_two || found_invalid || found_path);
  assert(found_invalid == 2);
  for (auto idx : eliminated_after_this_[bag_idx]) {
    assert(remaining_in_arcs_after_this_[bag_idx][idx] == 0);
    assert(remaining_out_arcs_after_this_[bag_idx][idx] == 0);
    assert(frontier[idx] != outgoing_index);
    assert(frontier[idx] != incoming_index);
    if (frontier[idx] == no_edge_index || frontier[idx] == two_edge_index) {
      frontier[idx] = no_edge_index;
    } else {
      if (frontier[idx] ^ incoming_flag) {
        frontier[frontier[idx]] = outgoing_index;
      } else {
        frontier[frontier[idx]] = incoming_index;
      }
      frontier[idx] = no_edge_index;
    }
  }
}

template <template <typename> typename Count_structure, typename count_t>
void TreewidthDirectedSearch<Count_structure, count_t>::print_stats() const {
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
  size_t merges = 0, unsuccessful_merges = 0;
  for (size_t i = 0; i < nthreads_; i++) {
    merges += merges_[i];
    unsuccessful_merges += unsuccessful_merges_[i];
  }
  LOG << "Cache hit rate: " << 100 * pos_hits / (double)(pos_hits + neg_hits)
      << "% (" << pos_hits << "/" << pos_hits + neg_hits << ")" << std::endl;
  LOG << "#Edges: " << edges << " #Propagations: " << propagations << std::endl;
  LOG << "#Merges: " << merges
      << " #Unsuccessful merges: " << unsuccessful_merges << std::endl;
}

} // namespace fpc
