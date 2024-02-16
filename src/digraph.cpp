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

#include "digraph.h"
#include "dgraph.h"
#include "logging.h"
#include "math.h"
#include "parallel_directed_search.h"
#include <algorithm>
#include <deque>
#include <limits>
#include <sstream>
#include <unordered_set>

namespace fpc {

Digraph::Digraph(std::istream &input)
    : nr_vertices_(0), nr_arcs_(0), max_length_(0) {
  char dec;
  input >> dec;
  Vertex nr_arcs;
  Vertex nr_vertices;
  std::string line;
  max_length_ = std::numeric_limits<Edge_length>::max();
  while (!input.eof()) {
    switch (dec) {
    case 'c':
      std::getline(input, line);
      break;
    case 'p': {
      std::string str;
      input >> str >> nr_vertices >> nr_arcs;
      assert(str == "edge");
      in_neighbors_ =
          std::vector<std::vector<Vertex>>(nr_vertices, std::vector<Vertex>());
      out_neighbors_ =
          std::vector<std::vector<Vertex>>(nr_vertices, std::vector<Vertex>());
      break;
    }
    case 'a':
      Vertex v, w;
      input >> v >> w;
      --v;
      --w;
      if (v != w && !has_arc(v, w)) {
        add_arc(v, w);
      } else {
        nr_arcs--;
      }
      break;
    case 'l':
      size_t tmp;
      input >> tmp;
      assert(tmp < std::numeric_limits<Edge_length>::max());
      max_length_ = Edge_length(tmp);
      extra_paths_ = 0;
      break;
    case 't':
      terminals_.resize(2);
      input >> terminals_[0] >> terminals_[1];
      terminals_[0]--;
      terminals_[1]--;
      break;
    default:
      LOG << "Invalid character " << dec << " at beginning of line."
          << std::endl;
      break;
    }
    input >> dec;
  }
  assert(nr_arcs_ == nr_arcs);
  assert(nr_vertices_ == nr_vertices);
}

void Digraph::preprocess() {
  bool found = true;
  while (found) {
    found = false;
    preprocess_start_goal_edges();
    Vertex cur_isolated_removed = preprocess_isolated();
    found |= cur_isolated_removed > 0;
    isolated_removed += cur_isolated_removed;
    preprocess_start_goal_edges();
    Vertex cur_unreachable_removed = preprocess_unreachable();
    found |= cur_unreachable_removed > 0;
    unreachable_removed += cur_unreachable_removed;
    if (!is_length_limited()) {
      Vertex cur_forwarder_removed = preprocess_forwarder();
      found |= cur_forwarder_removed > 0;
      forwarder_removed += cur_forwarder_removed;
    }
    preprocess_start_goal_edges();
    Vertex cur_dominator_arcs_removed = preprocess_dominator_arcs();
    found |= cur_dominator_arcs_removed > 0;
    dominator_arcs_removed += cur_dominator_arcs_removed;
    if (!found) {
      Vertex cur_position_determined = preprocess_position_determined(1);
      found |= cur_position_determined > 0;
      position_determined_removed += cur_position_determined;
    }
    if (!found) {
      preprocess_start_goal_edges();
      Vertex cur_unusable_arcs_removed = preprocess_unusable_arcs();
      found |= cur_unusable_arcs_removed > 0;
      unusable_arcs_removed += cur_unusable_arcs_removed;
    }
  }
  preprocess_start_goal_edges();
  assert(!has_arc(terminals_[0], terminals_[1]));
}

void Digraph::weighted_preprocess() {
  is_length_limited_ = is_length_limited();
  bit_upper_bound_ = bit_upper_bound();
  bool found = true;
  while (found) {
    found = false;
    preprocess_start_goal_edges();
    Vertex cur_isolated_removed = preprocess_isolated();
    found |= cur_isolated_removed > 0;
    isolated_removed += cur_isolated_removed;
    preprocess_start_goal_edges();
    Vertex cur_unreachable_removed = preprocess_unreachable();
    found |= cur_unreachable_removed > 0;
    unreachable_removed += cur_unreachable_removed;
    preprocess_start_goal_edges();
    Vertex cur_forwarder_removed = weighted_preprocess_forwarder();
    found |= cur_forwarder_removed > 0;
    forwarder_removed += cur_forwarder_removed;
    preprocess_start_goal_edges();
    Vertex cur_dominator_arcs_removed = preprocess_dominator_arcs();
    found |= cur_dominator_arcs_removed > 0;
    dominator_arcs_removed += cur_dominator_arcs_removed;
    if (!found) {
      preprocess_start_goal_edges();
      Vertex cur_unusable_arcs_removed = preprocess_unusable_arcs();
      found |= cur_unusable_arcs_removed > 0;
      unusable_arcs_removed += cur_unusable_arcs_removed;
    }
  }
  preprocess_start_goal_edges();
  assert(!has_arc(terminals_[0], terminals_[1]));
}

void Digraph::reduce_modulo_equivalence() {
  auto classes = equivalence_classes();
  for (auto &eq : classes) {
    for (auto other : eq | srg::views::drop(1)) {
      if (other == terminals_[0]) {
        terminals_[0] = eq.front();
      }
      if (other == terminals_[1]) {
        terminals_[1] = eq.front();
      }
      remove_vertex(other);
    }
  }
}

void Digraph::print_stats() const {
  if (isolated_removed)
    LOG << "Removed isolated: " << isolated_removed << std::endl;
  if (forwarder_removed)
    LOG << "Removed forwarder: " << forwarder_removed << std::endl;
  if (position_determined_removed)
    LOG << "Removed position determined: " << position_determined_removed
        << std::endl;
  if (unreachable_removed)
    LOG << "Removed unreachable: " << unreachable_removed << std::endl;
  if (dominator_arcs_removed)
    LOG << "Removed dominator arcs: " << dominator_arcs_removed << std::endl;
  if (unusable_arcs_removed)
    LOG << "Removed unusable arcs: " << unusable_arcs_removed << std::endl;
  if (max_length_decrease)
    LOG << "Max length decreased by: " << max_length_decrease << std::endl;
  if (arc_weights_.empty()) {
    LOG << "#vertices " << nr_vertices() << " (" << equivalence_classes().size()
        << ")"
        << " #arcs " << nr_arcs();
  } else {
    LOG << "#vertices " << nr_vertices() << " (?)"
        << " #arcs " << nr_arcs();
  }
  LOG << " #paths upper bound 2^" << bit_upper_bound();
  LOG << " max. length " << static_cast<size_t>(max_length_) << " min. length "
      << static_cast<size_t>(min_length()) << std::endl;
  LOG << "terminals: " << terminals_[0] << "," << terminals_[1] << std::endl;
}

void Digraph::print_ordered(std::ostream &os,
                            std::vector<std::size_t> const &order) const {
  os << "p edge " << nr_vertices() << " " << nr_arcs() << std::endl;
  std::set<Vertex> active;
  std::vector<char> taken(nr_vertices(), false);
  for (auto v : order) {
    active.insert(v);
    for (auto neigh : out_neighbors(v)) {
      if (taken[neigh]) {
        continue;
      }
      active.insert(neigh);
      os << "a " << v + 1 << " " << neigh + 1 << std::endl;
    }
    for (auto neigh : in_neighbors(v)) {
      if (taken[neigh]) {
        continue;
      }
      active.insert(neigh);
      os << "a " << neigh + 1 << " " << v + 1 << std::endl;
    }
    taken[v] = true;
    active.erase(v);
  }
  os << "l " << max_length() << std::endl;
  os << "t " << terminals().front() + 1 << " " << terminals().back() + 1
     << std::endl;
}

Edge_length Digraph::min_length() const {
  std::vector<Edge_length> distance_to_goal(
      in_neighbors_.size(), std::numeric_limits<Edge_length>::max());
  dijkstra(terminals_[1], distance_to_goal, false);
  return distance_to_goal[terminals_[0]];
}

std::size_t Digraph::bit_upper_bound() const {
  if (bit_upper_bound_.has_value()) {
    return *bit_upper_bound_;
  }
  assert(arc_weights_.empty());
  std::vector<mpz_class> degrees;
  for (auto v = 0; v < nr_vertices(); ++v) {
    if (terminals_[0] != v && terminals_[1] != v) {
      degrees.push_back(out_neighbors(v).size());
    }
  }
  srg::sort(degrees);
  auto limit = std::min(size_t(max_length()), nr_vertices());
  if (limit >= 2) {
    limit -= 2;
  } else {
    limit = 0;
  }
  mpz_class cur = 1;
  mpz_class res = 1;
  for (auto l = 0; l < limit; ++l) {
    cur *= degrees[degrees.size() - l - 1];
    res += cur;
  }
  res *= out_neighbors(terminals_[0]).size();
  res *= in_neighbors(terminals_[1]).size();
  res += extra_paths_;

  std::size_t bits = 1;
  mpz_class values = 1;
  while (values < res) {
    values *= 2;
    bits += 1;
  }
  return bits;
}

Digraph::Digraph(Vertex n)
    : nr_vertices_(0), nr_arcs_(0), max_length_(-1), in_neighbors_(n),
      out_neighbors_(n) {}

Digraph Digraph::subgraph(std::vector<Vertex> const &restrict_to) const {
  Digraph ret(restrict_to.size());
  ret.max_length_ = max_length_;
  auto unnamed = std::numeric_limits<Vertex>::max();
  std::vector<Vertex> new_name(in_neighbors_.size(), unnamed);
  Vertex cur_name = 0;
  for (Vertex v : restrict_to) {
    new_name[v] = cur_name++;
  }
  for (Vertex v : restrict_to) {
    for (Vertex neigh : in_neighbors(v)) {
      if (new_name[neigh] != unnamed) {
        ret.add_arc(new_name[neigh], new_name[v]);
        if (auto it = arc_weights_.find(Edge(neigh, v));
            it != arc_weights_.end()) {
          ret.arc_weights_.try_emplace(Edge(new_name[neigh], new_name[v]),
                                       it->second);
        }
      }
    }
  }
  ret.extra_paths_ = 0;
  if (!terminals_.empty()) {
    ret.terminals_ = {Vertex(-1), Vertex(-1)};
    if (new_name[terminals_[0]] != unnamed) {
      ret.terminals_[0] = new_name[terminals_[0]];
    }
    if (new_name[terminals_[1]] != unnamed) {
      ret.terminals_[1] = new_name[terminals_[1]];
    }
  }
  assert(ret.count_arcs() == ret.nr_arcs_);
  assert(ret.count_vertices() == ret.nr_vertices_);
  assert(ret.in_neighbors_.size() > 0);
  assert(ret.out_neighbors_.size() > 0);
  return ret;
}

Digraph Digraph::copy() const {
  std::vector<Vertex> all;
  all.reserve(in_neighbors_.size());
  for (Vertex v = 0; v < in_neighbors_.size(); v++) {
    all.push_back(v);
  }
  auto rv = subgraph(all);
  rv.extra_paths_ = extra_paths_;
  return rv;
}

void Digraph::normalize(bool reorder) {
  Vertex unnamed = std::numeric_limits<Vertex>::max();
  std::vector<Vertex> new_name(in_neighbors_.size(), unnamed);
  Vertex cur_name = 0;
  if (reorder) {
    new_name[terminals_[0]] = cur_name++;
    new_name[terminals_[1]] = cur_name++;
  }
  for (Vertex v = 0; v < in_neighbors_.size(); v++) {
    if (in_neighbors(v).empty() && out_neighbors(v).empty() &&
        v != terminals_[0] && v != terminals_[1]) {
      continue;
    }
    if (new_name[v] == unnamed) {
      new_name[v] = cur_name++;
    }
  }
  auto new_in_neighbors =
      std::vector<std::vector<Vertex>>(cur_name, std::vector<Vertex>());
  auto new_out_neighbors =
      std::vector<std::vector<Vertex>>(cur_name, std::vector<Vertex>());
  decltype(arc_weights_) new_arc_weights;
  for (Vertex v = 0; v < in_neighbors_.size(); v++) {
    if (in_neighbors(v).empty() && out_neighbors(v).empty()) {
      continue;
    }
    for (auto neigh : in_neighbors(v)) {
      assert(new_name[neigh] != unnamed);
      new_in_neighbors[new_name[v]].push_back(new_name[neigh]);
      new_out_neighbors[new_name[neigh]].push_back(new_name[v]);
      if (auto it = arc_weights_.find(Edge(neigh, v));
          it != arc_weights_.end()) {
        new_arc_weights.emplace(Edge(new_name[neigh], new_name[v]), it->second);
      }
    }
  }
  in_neighbors_ = new_in_neighbors;
  out_neighbors_ = new_out_neighbors;
  arc_weights_ = new_arc_weights;
  assert(new_name[terminals_[0]] != unnamed);
  assert(new_name[terminals_[1]] != unnamed);
  terminals_[0] = new_name[terminals_[0]];
  terminals_[1] = new_name[terminals_[1]];
}

bool Digraph::is_length_limited() const {
  if (is_length_limited_.has_value()) {
    return *is_length_limited_;
  }
  assert(arc_weights_.empty());
  return max_length_ + 1 < nr_vertices();
}

void Digraph::dijkstra(Vertex start, std::vector<Edge_length> &distance,
                       bool outgoing, std::optional<Vertex> forbidden) const {
  std::deque<Vertex> queue;
  queue.push_back(start);
  distance[start] = 0;
  while (!queue.empty()) {
    auto cur_vertex = queue.front();
    auto cur_cost = distance[cur_vertex];
    queue.pop_front();
    auto const &neighbors =
        outgoing ? out_neighbors(cur_vertex) : in_neighbors(cur_vertex);
    for (auto w : neighbors) {
      if (forbidden && w == *forbidden) {
        continue;
      }
      auto extra =
          arc_length(outgoing ? cur_vertex : w, outgoing ? w : cur_vertex);
      if (cur_cost + extra >= distance[w]) {
        continue;
      }
      if (cur_cost + extra < distance[w]) {
        distance[w] = extra + cur_cost;
        queue.push_back(w);
      }
    }
  }
}

std::vector<std::vector<Vertex>> Digraph::strong_connected_components() const {
  std::vector<Vertex> disc(nr_vertices_, 0);
  std::vector<Vertex> low(nr_vertices_, 0);
  std::vector<char> on_stack(nr_vertices_, false);
  std::vector<Vertex> stack;

  std::vector<std::vector<Vertex>> sccs;

  Vertex time = 0;
  scc_dfs(terminals_[0], disc, low, on_stack, stack, time, sccs);
  return sccs;
}

void Digraph::scc_dfs(Vertex start, std::vector<Vertex> &disc,
                      std::vector<Vertex> &low, std::vector<char> &on_stack,
                      std::vector<Vertex> &stack, Vertex &time,
                      std::vector<std::vector<Vertex>> &sccs) const {
  disc[start] = low[start] = ++time;
  on_stack[start] = true;
  stack.push_back(start);
  for (auto child : out_neighbors(start)) {
    if (disc[child] == 0) {
      scc_dfs(child, disc, low, on_stack, stack, time, sccs);
      low[start] = std::min(low[start], low[child]);
    } else if (on_stack[child]) {
      low[start] = std::min(low[start], disc[child]);
    }
  }
  if (low[start] == disc[start]) {
    // on top of the scc
    std::vector<Vertex> scc;
    while (stack.back() != start) {
      on_stack[stack.back()] = false;
      scc.push_back(stack.back());
      stack.pop_back();
    }
    scc.push_back(start);
    on_stack[start] = false;
    stack.pop_back();
    // assert that if the graph is not trivial then we find the scc containing
    // the goal first
    assert(nr_arcs_ == 0 || !sccs.empty() || scc.size() != 1 ||
           scc.front() == terminals_[1]);
    sccs.push_back(std::move(scc));
  }
}

std::unordered_set<Vertex>
Digraph::must_use_subgraph(Vertex must_use,
                           std::vector<Edge_length> const &distance_from_start,
                           std::vector<Edge_length> const &distance_to_goal,
                           Edge_length backward_allowance) const {
  assert(arc_weights_.empty());
  std::unordered_set<Vertex> restrict_to;
  auto find_usable = [&restrict_to, must_use, backward_allowance,
                      this](auto const &distance, bool outgoing) {
    std::vector<Edge_length> min_allowance_used(in_neighbors_.size(),
                                                max_length_);
    min_allowance_used[must_use] = 0;
    min_allowance_used[terminals_[0]] = 0;
    min_allowance_used[terminals_[1]] = 0;
    std::deque<Vertex> fifo;
    fifo.push_back(must_use);
    while (!fifo.empty()) {
      auto cur_vertex = fifo.front();
      auto cur_allowance = min_allowance_used[cur_vertex];
      assert(cur_allowance <= backward_allowance);
      fifo.pop_front();
      for (auto next :
           outgoing ? out_neighbors(cur_vertex) : in_neighbors(cur_vertex)) {
        if (distance[next] < distance[cur_vertex] &&
            cur_allowance < min_allowance_used[next]) {
          min_allowance_used[next] = cur_allowance;
          fifo.push_back(next);
        } else if (int new_allowance = cur_allowance + distance[next] -
                                       distance[cur_vertex] + 1;
                   new_allowance < min_allowance_used[next] &&
                   new_allowance <= backward_allowance) {
          assert(new_allowance >= 0);
          assert(new_allowance <= backward_allowance);
          min_allowance_used[next] = new_allowance;
          fifo.push_back(next);
        }
      }
    }
    for (Vertex v = 0; v < in_neighbors_.size(); v++) {
      if (min_allowance_used[v] <= backward_allowance) {
        restrict_to.insert(v);
      }
    }
  };
  find_usable(distance_from_start, false);
  find_usable(distance_to_goal, true);
  return restrict_to;
}

std::vector<std::vector<Vertex>> Digraph::equivalence_classes() const {
  assert(arc_weights_.empty());
  std::unordered_map<std::vector<int>, std::vector<Vertex>, vector_hash<int>>
      false_twins;
  std::unordered_map<std::vector<int>, std::vector<Vertex>, vector_hash<int>>
      true_twins;
  for (Vertex v = 0; v < in_neighbors_.size(); v++) {
    std::vector<int> sorted_neighbors;
    for (auto out : out_neighbors(v)) {
      sorted_neighbors.push_back(-int(out + 1));
    }
    for (auto in : in_neighbors(v)) {
      sorted_neighbors.push_back(int(in + 1));
    }
    srg::sort(sorted_neighbors);
    auto [it, inserted_false] =
        false_twins.try_emplace(sorted_neighbors, std::vector<Vertex>{v});
    if (!inserted_false) {
      it->second.push_back(v);
    }
    sorted_neighbors.push_back(int(v + 1));
    sorted_neighbors.push_back(-int(v + 1));
    srg::sort(sorted_neighbors);
    auto [other_it, inserted_true] =
        true_twins.try_emplace(sorted_neighbors, std::vector<Vertex>{v});
    if (!inserted_true) {
      other_it->second.push_back(v);
    }
    assert(inserted_false || inserted_true);
  }
  std::vector<std::vector<Vertex>> rv;
  std::vector<char> seen(in_neighbors_.size(), false);
  for (auto &[_, eq_class] : false_twins) {
    if (eq_class.size() > 1) {
      rv.push_back(eq_class);
      for (auto v : eq_class) {
        assert(!seen[v]);
        seen[v] = true;
      }
    }
  }

  for (auto &[_, eq_class] : true_twins) {
    if (eq_class.size() > 1) {
      rv.push_back(eq_class);
      for (auto v : eq_class) {
        assert(!seen[v]);
        seen[v] = true;
      }
    }
  }
  for (Vertex v = 0; v < in_neighbors_.size(); v++) {
    if (!seen[v]) {
      rv.push_back({v});
    }
  }

  return rv;
}

void Digraph::print_dfvs_instance() const {
  std::cout << nr_vertices() << " " << nr_arcs() << std::endl;
  for (auto v = 0; v < nr_vertices(); v++) {
    for (auto neigh : out_neighbors(v)) {
      std::cout << neigh << " ";
    }
    std::cout << std::endl;
  }
}

mpz_class Digraph::naive_dfs() const {
  assert(arc_weights_.empty());
  std::vector<char> on_stack(in_neighbors_.size(), false);
  std::vector<std::pair<std::size_t, std::size_t>> stack;
  std::vector<Edge_length> distance_to_goal(in_neighbors_.size(), max_length_);
  dijkstra(terminals_[1], distance_to_goal, false);
  stack.emplace_back(0, terminals_[0]);
  on_stack[terminals_[0]] = true;
  uint128_t res = 0;
  while (!stack.empty()) {
    auto [neigh_idx, cur] = stack.back();
    if (neigh_idx < out_neighbors_[cur].size()) {
      stack.back().first++;
      // add the child
      auto neigh = out_neighbors_[cur][neigh_idx];
      if (stack.size() + distance_to_goal[neigh] > max_length_) {
        continue;
      }
      if (neigh == terminals_[1]) {
        ++res;
        continue;
      }
      if (on_stack[neigh]) {
        continue;
      }
      on_stack[neigh] = true;
      stack.emplace_back(0, neigh);
    } else {
      // nothing left. self has been handled before being added to stack
      stack.pop_back();
      on_stack[cur] = false;
    }
  }
  return to_mpz(res);
}

void Digraph::contract_arcs(Vertex v, bool unsafe) {
  for (auto in : in_neighbors(v)) {
    for (auto out : out_neighbors(v)) {
      if (in == out) {
        assert(unsafe);
        continue;
      }
      if (!has_arc(in, out)) {
        arc_weights_.try_emplace(
            Edge(in, out),
            Limited_count<mpz_class>(max_length_, max_length_, {}));
        add_arc(in, out);
      } else {
        arc_weights_.try_emplace(Edge(in, out),
                                 Limited_count<mpz_class>(max_length_, 1));
      }
      Limited_count<mpz_class> in_to_v(max_length_, 1);
      if (auto it = arc_weights_.find(Edge(in, v)); it != arc_weights_.end()) {
        in_to_v = it->second;
      }
      Limited_count<mpz_class> v_to_out(max_length_, 1);
      if (auto it = arc_weights_.find(Edge(v, out)); it != arc_weights_.end()) {
        v_to_out = it->second;
      }
      in_to_v *= v_to_out;
      arc_weights_.at(Edge(in, out)) += in_to_v;
    }
  }
  remove_vertex(v);
}

Vertex Digraph::preprocess_start_goal_edges() {
  assert(terminals_.size() == 2);
  Vertex found = 0;
  if (has_arc(terminals_[0], terminals_[1])) {
    found = 1;

    if (auto it = arc_weights_.find(Edge(terminals_[0], terminals_[1]));
        it != arc_weights_.end()) {
      extra_paths_ += it->second.total_count(max_length_);
    } else {
      extra_paths_ += 1;
    }

    remove_arc(terminals_[0], terminals_[1]);
  }
  auto out_end = out_neighbors(terminals_[1]);
  for (auto out : out_end) {
    remove_arc(terminals_[1], out);
  }
  auto in_start = in_neighbors(terminals_[0]);
  for (auto in : in_start) {
    remove_arc(in, terminals_[0]);
  }
  return found;
}

Vertex Digraph::preprocess_isolated() {
  assert(terminals_.size() == 2);
  Vertex found = 0;
  for (size_t v = 0; v < out_neighbors_.size(); v++) {
    if (in_neighbors(v).empty() && out_neighbors(v).empty()) {
      continue;
    }
    if (v == terminals_[0] || v == terminals_[1]) {
      // handled separately
      continue;
    }
    if (in_neighbors(v).size() == 0) {
      found++;
      remove_vertex(v);
    } else if (out_neighbors(v).size() == 0) {
      found++;
      remove_vertex(v);
    }
  }

  auto do_start = [&]() {
    if (out_neighbors(terminals_[0]).size() != 1) {
      return 0;
    }
    Vertex neighbor = out_neighbors(terminals_[0]).front();
    if (neighbor == terminals_[1]) {
      preprocess_start_goal_edges();
      return 1;
    }
    assert(max_length_ > 0);
    found++;
    if (auto it = arc_weights_.find(Edge(terminals_[0], neighbor));
        it != arc_weights_.end()) {
      while (it->second.offset() > 0) {
        max_length_--;
        it->second.decrement_offset();
      }
      // cannot come from anywhere to the new terminal
      auto neighs = in_neighbors(neighbor);
      for (auto other : neighs) {
        if (other == terminals_[0]) {
          continue;
        }
        remove_arc(other, neighbor);
      }
      contract_arcs(neighbor);
    } else {
      max_length_--;
      remove_vertex(terminals_[0]);
      terminals_[0] = neighbor;
    }
    return 1;
  };

  auto do_end = [&]() {
    if (in_neighbors(terminals_[1]).size() != 1) {
      return 0;
    }
    Vertex neighbor = in_neighbors(terminals_[1]).front();
    if (neighbor == terminals_[0]) {
      preprocess_start_goal_edges();
      return 1;
    }
    assert(max_length_ > 0);
    found++;
    if (auto it = arc_weights_.find(Edge(neighbor, terminals_[1]));
        it != arc_weights_.end()) {
      while (it->second.offset() > 0) {
        max_length_--;
        it->second.decrement_offset();
      }
      // cannot go anywhere from the new terminal
      auto neighs = out_neighbors(neighbor);
      for (auto other : neighs) {
        if (other == terminals_[1]) {
          continue;
        }
        remove_arc(neighbor, other);
      }
      contract_arcs(neighbor);
    } else {
      max_length_--;
      remove_vertex(terminals_[1]);
      terminals_[1] = neighbor;
    }
    return 1;
  };
  found += do_start();
  found += do_end();

  return found;
}

Vertex Digraph::preprocess_unreachable() {
  assert(terminals_.size() == 2);
  std::vector<Edge_length> distance_from_start(in_neighbors_.size(),
                                               max_length_ + 1);
  dijkstra(terminals_[0], distance_from_start, true);

  std::vector<Edge_length> distance_to_goal(in_neighbors_.size(),
                                            max_length_ + 1);
  dijkstra(terminals_[1], distance_to_goal, false);

  Vertex found = 0;
  for (Vertex v = 0; v < in_neighbors_.size(); v++) {
    if (in_neighbors(v).empty() && out_neighbors(v).empty()) {
      continue;
    }
    if (distance_from_start[v] + distance_to_goal[v] > max_length_) {
      found++;
      remove_vertex(v);
    }
  }
  return found;
}

Vertex Digraph::preprocess_forwarder() {
  assert(arc_weights_.empty());
  Vertex found = 0;
  for (size_t v = 0; v < in_neighbors_.size(); v++) {
    if (terminals_[0] == v || terminals_[1] == v) {
      // we cannot remove terminals this way
      continue;
    }
    if (in_neighbors(v).size() == 1) {
      Vertex in = in_neighbors(v).front();
      // currently not counting edge multiplicity
      if (srg::any_of(out_neighbors(v),
                      [&](auto out) { return in == out || has_arc(in, out); }))
        continue;

      found++;
      for (auto out : out_neighbors(v)) {
        add_arc(in, out);
      }
      remove_vertex(v);
    } else if (out_neighbors(v).size() == 1) {
      Vertex out = out_neighbors(v).front();
      // currently not counting edge multiplicity
      if (srg::any_of(in_neighbors(v),
                      [&](auto in) { return in == out || has_arc(in, out); }))
        continue;

      found++;
      for (auto in : in_neighbors(v)) {
        add_arc(in, out);
      }
      remove_vertex(v);
    } else if (in_neighbors(v).size() == 2 && out_neighbors(v).size() == 2) {
      auto in_1 = in_neighbors(v).front();
      auto in_2 = in_neighbors(v).back();
      auto out_1 = out_neighbors(v).front();
      auto out_2 = out_neighbors(v).back();
      if (in_1 != out_1 && in_1 != out_2) {
        continue;
      }
      if (in_2 != out_1 && in_2 != out_2) {
        continue;
      }
      if (has_arc(in_1, in_2) || has_arc(in_2, in_1)) {
        continue;
      }
      assert(in_1 != in_2);
      found++;
      add_arc(in_1, in_2);
      add_arc(in_2, in_1);

      remove_vertex(v);
    }
  }
  return found;
}

Vertex Digraph::weighted_preprocess_forwarder() {
  Vertex found = 0;
  for (size_t v = 0; v < in_neighbors_.size(); v++) {
    if (terminals_[0] == v || terminals_[1] == v) {
      // we cannot remove terminals this way
      continue;
    }
    if (in_neighbors(v).size() == 1) {
      found++;
      contract_arcs(v, true);
    } else if (out_neighbors(v).size() == 1) {
      found++;
      contract_arcs(v, true);
    } else if (in_neighbors(v).size() == 2 && out_neighbors(v).size() == 2) {
      auto in_1 = in_neighbors(v).front();
      auto in_2 = in_neighbors(v).back();
      auto out_1 = out_neighbors(v).front();
      auto out_2 = out_neighbors(v).back();
      if (in_1 != out_1 && in_1 != out_2) {
        continue;
      }
      if (in_2 != out_1 && in_2 != out_2) {
        continue;
      }
      assert(in_1 != in_2);
      found++;
      contract_arcs(v, true);
    }
  }
  return found;
}

Vertex Digraph::preprocess_dominator_arcs() {
  assert(terminals_.size() == 2);

  DominatorGraph d_graph;
  int *arclist = new int[2 * nr_arcs()];
  auto cur = 0;
  for (Vertex v = 0; v < in_neighbors_.size(); v++) {
    for (auto w : out_neighbors(v)) {
      arclist[cur++] = v + 1;
      arclist[cur++] = w + 1;
    }
  }
  d_graph.buildGraph(in_neighbors_.size(), nr_arcs(), terminals_[0] + 1,
                     arclist, false);
  int *idom = new int[in_neighbors_.size() + 1];
  d_graph.slt(d_graph.getSource(), idom);

  std::unordered_set<Vertex> to_check;
  Vertex found_forward = 0;
  assert(idom[0] == 0);
  for (Vertex v = 1; v < in_neighbors_.size() + 1; v++) {
    auto w = idom[v] - 1;
    if (w != -1) {
      to_check.emplace(w);
      if (has_arc(v - 1, w)) {
        remove_arc(v - 1, w);
        found_forward++;
      }
    }
  }

  d_graph.destroy();
  cur = 0;
  for (Vertex v = 0; v < in_neighbors_.size(); v++) {
    for (auto w : out_neighbors(v)) {
      arclist[cur++] = w + 1;
      arclist[cur++] = v + 1;
    }
  }
  d_graph.buildGraph(in_neighbors_.size(), nr_arcs(), terminals_[1] + 1,
                     arclist, false);
  srg::fill(idom, idom + in_neighbors_.size() + 1, 0);
  d_graph.slt(d_graph.getSource(), idom);

  to_check.clear();
  Vertex found_backward = 0;
  for (Vertex v = 1; v < in_neighbors_.size() + 1; v++) {
    auto w = idom[v] - 1;
    if (w != -1) {
      to_check.emplace(w);
      if (has_arc(w, v - 1)) {
        remove_arc(w, v - 1);
        found_backward++;
      }
    }
  }
  delete[] arclist;
  delete[] idom;
  return found_forward + found_backward;
}

Vertex Digraph::preprocess_unusable_arcs() {
  assert(terminals_.size() == 2);

  if (nr_arcs_ > 10000) {
    return 0;
  }

  Vertex found = 0;
  // forward edges
  std::vector<Edge_length> normal_distance_from_start(in_neighbors_.size(),
                                                      max_length_ + 1);
  dijkstra(terminals_[0], normal_distance_from_start, true);
  std::vector<Edge_length> local_distance_to_goal(in_neighbors_.size(),
                                                  max_length_ + 1);
  for (Vertex v = 0; v < in_neighbors_.size(); v++) {
    if (in_neighbors(v).empty() && out_neighbors(v).empty()) {
      continue;
    }
    if (v == terminals_[0])
      continue;
    if (v == terminals_[1]) {
      continue;
    }
    srg::fill(local_distance_to_goal, max_length_ + 1);
    dijkstra(terminals_[1], local_distance_to_goal, false, v);

    std::vector<Vertex> remove_completely;
    for (auto w : out_neighbors(v)) {
      Edge_length min_with_arc = normal_distance_from_start[v] +
                                 local_distance_to_goal[w] + arc_length(v, w);
      if (min_with_arc > max_length_) {
        found++;
        remove_completely.push_back(w);
      }
    }
    for (auto w : remove_completely) {
      remove_arc(v, w);
    }
  }

  // backward edges
  std::vector<Edge_length> normal_distance_to_goal(in_neighbors_.size(),
                                                   max_length_ + 1);
  std::vector<Edge_length> local_distance_from_start(in_neighbors_.size(),
                                                     max_length_ + 1);
  dijkstra(terminals_[1], normal_distance_to_goal, false);
  for (Vertex v = 0; v < in_neighbors_.size(); v++) {
    if (in_neighbors(v).empty() && out_neighbors(v).empty()) {
      continue;
    }
    if (v == terminals_[0])
      continue;
    if (v == terminals_[1]) {
      continue;
    }

    srg::fill(local_distance_from_start, max_length_ + 1);
    dijkstra(terminals_[0], local_distance_from_start, true, v);
    std::vector<Vertex> remove_completely;
    for (auto w : in_neighbors(v)) {
      Edge_length min_with_arc = local_distance_from_start[w] +
                                 normal_distance_to_goal[v] + arc_length(w, v);
      if (min_with_arc > max_length_) {
        found++;
        remove_completely.push_back(w);
      }
    }
    for (auto w : remove_completely) {
      remove_arc(w, v);
    }
  }
  return found;
}

Vertex Digraph::preprocess_position_determined(Edge_length budget) {
  assert(!terminals_.empty());
  assert(arc_weights_.empty());
  if (max_length_ <= budget) {
    return 0;
  }

  if (nr_arcs_ > 10000) {
    return 0;
  }

  if (nr_arcs_ == 0) {
    return 0;
  }
  std::vector<Edge_length> distance_from_start(in_neighbors_.size(),
                                               max_length_ + 1);
  dijkstra(terminals_[0], distance_from_start, true);

  std::vector<Edge_length> distance_to_goal(in_neighbors_.size(),
                                            max_length_ + 1);
  dijkstra(terminals_[1], distance_to_goal, false);

  assert(!out_neighbors(terminals_[0]).empty());
  assert(!in_neighbors(terminals_[1]).empty());
  Vertex found = 0;
  for (Vertex v = 0; v < in_neighbors_.size(); v++) {
    if (in_neighbors(v).empty() && out_neighbors(v).empty()) {
      continue;
    }
    if (v == terminals_[0])
      continue;
    if (v == terminals_[1]) {
      continue;
    }

    if (distance_from_start[v] + distance_to_goal[v] >= max_length_ - budget) {
      found++;
      auto relevant =
          must_use_subgraph(v, distance_from_start, distance_to_goal, budget);
      assert(relevant.contains(v) && relevant.contains(terminals_[0]) &&
             relevant.contains(terminals_[1]));
      // put v, start, and goal to the back so we know where they are
      relevant.erase(v);
      relevant.erase(terminals_[0]);
      relevant.erase(terminals_[1]);
      std::vector<Vertex> restrict_to(relevant.begin(), relevant.end());
      restrict_to.push_back(v);
      restrict_to.push_back(terminals_[0]);
      restrict_to.push_back(terminals_[1]);
      auto search_graph = subgraph(restrict_to);
      assert(
          !search_graph.out_neighbors(Vertex(restrict_to.size() - 2)).empty());
      assert(
          !search_graph.in_neighbors(Vertex(restrict_to.size() - 1)).empty());
      search_graph.terminals_ = {Vertex(restrict_to.size() - 2),
                                 Vertex(restrict_to.size() - 1)};

      // TODO: use more threads
      auto solver = make_solver<ParallelDirectedSearch>(
          search_graph, Vertex(restrict_to.size() - 3), 1);
      extra_paths_ += solver->search();
      remove_vertex(v);
      break;
    }
  }
  return found;
}
} // namespace fpc
