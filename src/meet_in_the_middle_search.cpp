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

#include "meet_in_the_middle_search.h"

namespace fpc {

template <template <typename> typename Count_structure, typename count_t>
MITMSearch<Count_structure, count_t>::MITMSearch(Graph const &input,
                                                 size_t nthreads)
    : nr_vertices_(input.nr_vertices()), nthreads_(nthreads),
      max_length_(
          std::min(input.max_length(), Edge_length(input.nr_vertices()))),
      terminals_(input.terminals()), graph_(input), neighbors_(nr_vertices_),
      invalid_(std::numeric_limits<Edge_length>::max() - max_length_ - 1),
      distance_from_start_(nr_vertices_, invalid_),
      distance_to_goal_(nr_vertices_, invalid_), cache_(nr_vertices_ + 1),
      local_result_(nr_vertices_, 0) {
  assert(terminals_.size() == 2);
  for (Vertex v = 0; v < nr_vertices_; v++) {
    // fill neighbors
    neighbors_[v] = input.neighbors(v);
  }
  input.dijkstra(terminals_[0], distance_from_start_, {});
  input.dijkstra(terminals_[1], distance_to_goal_, {});
  omp_set_num_threads(nthreads_);
}

template <template <typename> typename Count_structure, typename count_t>
mpz_class MITMSearch<Count_structure, count_t>::search() {
  // strongly inspired by diodrm implementation

  auto eq_classes = graph_.equivalence_classes();
  std::vector<Vertex> representatives;
  std::vector<Vertex> size;
  for (auto &eq_class : eq_classes) {
    auto found_start = srg::find(eq_class, terminals_[0]) != eq_class.end();
    auto found_goal = srg::find(eq_class, terminals_[1]) != eq_class.end();
    for (auto v : eq_class) {
      if (v == terminals_[0]) {
        continue;
      }
      if (v == terminals_[1]) {
        continue;
      }
      representatives.push_back(v);
      size.push_back(eq_class.size());
      if (found_goal) {
        size.back()--;
      }
      if (found_start) {
        size.back()--;
      }
      break;
    }
  }
  assert(representatives.size() == size.size());

  // maximum m-s path length (m: midpoint of paths)
  const int ms_maxlen = (max_length_ + 1) / 2;
  // maximum m-t path length (m: midpoint of paths)
  const int mt_maxlen = max_length_ / 2;

#pragma omp parallel for schedule(dynamic, 1) default(shared)
  for (auto idx = 0; idx < representatives.size(); ++idx) {
    // midpoint
    auto mid = representatives[idx];

    assert(mid != terminals_[0]);
    assert(mid != terminals_[1]);

    if (distance_from_start_[mid] > ms_maxlen ||
        distance_to_goal_[mid] > mt_maxlen) {
      // wont find relevant paths here
      continue;
    }

    enumerate_ms_paths(mid, ms_maxlen, mt_maxlen);
    enumerate_mt_paths(mid, ms_maxlen, mt_maxlen);

    cache_[mid] = {};
  }
  mpz_class res = 0;
  for (auto idx = 0; idx < representatives.size(); ++idx) {
    auto mid = representatives[idx];
    res += to_mpz(local_result_[mid] * size[idx]);
  }
  return res;
}

template <template <typename> typename Count_structure, typename count_t>
std::pair<double, double>
MITMSearch<Count_structure, count_t>::percent_done_after(
    std::size_t paths) const {
  assert(nr_vertices_ > 2);
  assert(paths > 0);

  const int ms_max_len = (max_length_ + 1) / 2;
  const int mt_max_len = max_length_ / 2;

  std::size_t steps_taken = 0;
  std::size_t emergency_break = std::min(paths * 100, std::size_t(100'000'000));
  std::size_t path_size_sum = 0;
  std::size_t found_paths = 0;
  auto compute_for = [&](auto mid) {
    std::vector<Edge_length> distance_from_mid(nr_vertices_, invalid_);
    graph_.dijkstra(mid, distance_from_mid, {});
    std::vector<Vertex> path;
    std::vector<char> visited(nr_vertices_, false);
    visited[mid] = true;

    Edge_length length = 0;
    // run dfs for enumerating mid-start paths
    struct Data {
      Vertex neigh_idx;
      Vertex u;
    };
    std::vector<Data> stack;
    stack.emplace_back(0, mid);
    auto percent_done = [&]() {
      double rv = 0.0;
      visited[terminals_[0]] = false;
      --length;
      for (auto [neigh_idx, u] :
           stack | std::views::reverse | std::views::drop(1)) {
        double nr_real_neighbors = 0.0;
        double nr_real_neighbors_done = 0.0;
        std::size_t cur_idx = 0;
        for (auto v : neighbors(u)) {
          cur_idx++;
          if (visited[v]) {
            continue;
          }
          if (length + 1 + distance_from_start_[v] > ms_max_len) {
            continue;
          }
          if (v == terminals_[1]) {
            continue;
          }
          if (neigh_idx > cur_idx) {
            nr_real_neighbors_done += 1.0;
          }
          nr_real_neighbors += 1.0;
        }
        if (nr_real_neighbors > 0.0) {
          assert(nr_real_neighbors_done >= 0.0);
          // rv is currently how much the neighbor we took was processed.
          // add the number of other neighbors we completely processed.
          rv += nr_real_neighbors_done;
          // get the percentage of neighbors processed here
          rv /= nr_real_neighbors;
        }
        visited[u] = false;
        assert(length != 0 || u == mid);
        --length;
      }
      return rv;
    };
    while (!stack.empty()) {
      auto [neigh_idx, u] = stack.back();
      if (u == terminals_[0]) {
        if (u != mid &&
            distance_from_mid[u] + distance_to_goal_[u] <= mt_max_len) {
          path.pop_back();
        }
        path_size_sum += path.size();
        if (++found_paths == paths) {
          // calculate percent done
          return percent_done();
        }
        if (steps_taken > emergency_break) {
          if (found_paths == 0) {
            found_paths = 1;
            path_size_sum = 100;
          }
          paths = found_paths;
          return percent_done();
        }
        --length;
        visited[u] = false;
        stack.pop_back();
        continue;
      }
      if (neigh_idx < neighbors(u).size()) {
        ++steps_taken;
        if (steps_taken > emergency_break && length >= 10) {
          if (found_paths == 0) {
            found_paths = 1;
            path_size_sum = 100;
          }
          paths = found_paths;
          return percent_done();
        }
        stack.back().neigh_idx++;
        auto v = neighbors(u)[neigh_idx];
        if (visited[v]) {
          continue;
        }
        if (length + 1 + distance_from_start_[v] > ms_max_len) {
          continue;
        }
        if (v == terminals_[1]) {
          continue;
        }
        visited[v] = true;
        ++length;
        if (distance_from_mid[v] + distance_to_goal_[v] <= mt_max_len) {
          path.push_back(v);
        }
        stack.emplace_back(0, v);
      } else {
        if (u != mid &&
            distance_from_mid[u] + distance_to_goal_[u] <= mt_max_len) {
          path.pop_back();
        }
        --length;
        visited[u] = false;
        stack.pop_back();
      }
    }

    return 1.0;
  };

  // choose a representative middle vertex
  Vertex mid = -1;
  Edge_length min_sum = Edge_length(-1);
  Edge_length min_diff = Edge_length(-1);
  for (Vertex v = 0; v < graph_.nr_vertices(); v++) {
    auto this_sum = distance_from_start_[v] + distance_to_goal_[v];
    auto this_diff =
        this_sum - 2 * std::min(distance_from_start_[v], distance_to_goal_[v]);
    if (this_sum < min_sum) {
      mid = v;
      min_sum = this_sum;
      min_diff = this_diff;
    } else if (this_sum == min_sum && this_diff < min_diff) {
      mid = v;
      min_sum = this_sum;
      min_diff = this_diff;
    }
  }
  auto rv = compute_for(mid);
  while (found_paths < paths) {
    auto remaining = paths - found_paths;
    mid++;
    mid %= nr_vertices_;
    if (distance_from_start_[mid] == 0) {
      mid++;
      mid %= nr_vertices_;
    }
    if (distance_to_goal_[mid] == 0) {
      mid++;
      mid %= nr_vertices_;
    }
    rv = compute_for(mid);
    rv = std::min(1.0, double(paths) / remaining * rv);
  }
  return {rv, double(path_size_sum) / paths};
}

template <template <typename> typename Count_structure, typename count_t>
void MITMSearch<Count_structure, count_t>::insert_all_subsets(
    std::vector<Vertex> const &path, Vertex mid, Edge_length length) {
  auto &count = cache_[mid];

  auto sorted_path = path;
  if (!sorted_path.empty() && sorted_path.back() == terminals_[0]) {
    sorted_path.pop_back(); // remove start
  }
  auto const p_size = sorted_path.size();

  // order is irrelevant for what were doing
  std::sort(sorted_path.begin(), sorted_path.end());

  std::vector<Vertex> subset;
  subset.reserve(p_size);
  for (int bits = 0; bits < 1 << p_size; ++bits) {
    subset.clear();
    // elements are added in ascending order
    for (int i = 0; i < p_size; ++i) {
      if ((bits >> i) & 1) {
        subset.push_back(sorted_path[i]);
      }
    }
    auto [it, _] =
        count.try_emplace(subset, Limited_count<count_t>::zero(max_length_));
    it->second.add_one(length);
  }
}

template <template <typename> typename Count_structure, typename count_t>
void MITMSearch<Count_structure, count_t>::enumerate_ms_paths(
    Vertex mid, Edge_length ms_max_len, Edge_length mt_max_len) {
  std::vector<Edge_length> distance_from_mid(nr_vertices_, invalid_);
  graph_.dijkstra(mid, distance_from_mid, {});

  std::vector<char> visited(nr_vertices_, false);
  visited[mid] = true;

  std::vector<Vertex> path;
  Edge_length length = 0;
  // run dfs for enumerating mid-start paths
  struct Data {
    Vertex neigh_idx;
    Vertex u;
  };
  std::vector<Data> stack;
  stack.emplace_back(0, mid);
  while (!stack.empty()) {
    auto [neigh_idx, u] = stack.back();
    if (u == terminals_[0]) {
      insert_all_subsets(path, mid, length);
      if (u != mid &&
          distance_from_mid[u] + distance_to_goal_[u] <= mt_max_len) {
        path.pop_back();
      }
      --length;
      visited[u] = false;
      stack.pop_back();
      continue;
    }
    if (neigh_idx < neighbors(u).size()) {
      stack.back().neigh_idx++;
      auto v = neighbors(u)[neigh_idx];
      if (visited[v]) {
        continue;
      }
      if (length + 1 + distance_from_start_[v] > ms_max_len) {
        continue;
      }
      if (v == terminals_[1]) {
        continue;
      }
      visited[v] = true;
      ++length;
      if (distance_from_mid[v] + distance_to_goal_[v] <= mt_max_len) {
        path.push_back(v);
      }
      stack.emplace_back(0, v);
    } else {
      if (u != mid &&
          distance_from_mid[u] + distance_to_goal_[u] <= mt_max_len) {
        path.pop_back();
      }
      --length;
      visited[u] = false;
      stack.pop_back();
    }
  }
}

template <template <typename> typename Count_structure, typename count_t>
void MITMSearch<Count_structure, count_t>::inclusion_exclusion(
    std::vector<Vertex> const &path, Vertex mid, Edge_length length) {
  auto &count = cache_[mid];
  auto sorted_path = path;
  if (!sorted_path.empty() && sorted_path.back() == terminals_[1]) {
    sorted_path.pop_back(); // remove goal
  }
  auto const p_size = sorted_path.size();
  std::sort(sorted_path.begin(), sorted_path.end());

  std::vector<Vertex> subset;
  // Here we run dfs for enumerating subsets since the number of valid
  // subsets is expected to be small.
  auto all_subsets = [&](auto rec, Edge_length nxt_begin) -> void {
    auto it = count.find(subset);
    // pruning
    if (it == count.end()) {
      return;
    }

    const auto &c = it->second;
    int sgn = (subset.size() & 1) ? -1 : +1;
    if (c.size() >= length) {
      local_result_[mid] += sgn * c[length];
      if (c.size() >= length + 1) {
        local_result_[mid] += sgn * c[length + 1];
      }
    }

    // enumerating larger subsets by recursive calling
    for (Edge_length nxt = nxt_begin; nxt < p_size; ++nxt) {
      subset.push_back(sorted_path[nxt]);
      rec(rec, nxt + 1);
      subset.pop_back();
    }
  };
  all_subsets(all_subsets, 0);
}

template <template <typename> typename Count_structure, typename count_t>
void MITMSearch<Count_structure, count_t>::enumerate_mt_paths(
    Vertex mid, Edge_length ms_max_len, Edge_length mt_max_len) {
  std::vector<Edge_length> distance_to_mid(nr_vertices_, invalid_);
  graph_.dijkstra(mid, distance_to_mid, {});

  std::vector<char> visited(nr_vertices_, false);
  visited[mid] = true;

  std::vector<Vertex> path;
  Edge_length length = 0;
  // run dfs for enumerating mid-goal paths
  struct Data {
    Vertex neigh_idx;
    Vertex u;
  };
  std::vector<Data> stack;
  stack.emplace_back(0, mid);
  while (!stack.empty()) {
    auto [neigh_idx, u] = stack.back();
    if (u == terminals_[1]) {
      inclusion_exclusion(path, mid, length);
      if (u != mid &&
          distance_to_mid[u] + distance_from_start_[u] <= ms_max_len) {
        path.pop_back();
      }
      --length;
      visited[u] = false;
      stack.pop_back();
      continue;
    }
    if (neigh_idx < neighbors(u).size()) {
      stack.back().neigh_idx++;
      auto v = neighbors(u)[neigh_idx];
      if (visited[v]) {
        continue;
      }
      if (length + 1 + distance_to_goal_[v] > mt_max_len) {
        continue;
      }
      if (v == terminals_[0]) {
        continue;
      }
      visited[v] = true;
      ++length;
      if (distance_to_mid[v] + distance_from_start_[v] <= ms_max_len) {
        path.push_back(v);
      }
      stack.emplace_back(0, v);
    } else {
      if (u != mid &&
          distance_to_mid[u] + distance_from_start_[u] <= ms_max_len) {
        path.pop_back();
      }
      --length;
      visited[u] = false;
      stack.pop_back();
    }
  }
}

template <template <typename> typename Count_structure, typename count_t>
void MITMSearch<Count_structure, count_t>::print_stats() const {}

} // namespace fpc
