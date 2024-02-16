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

#include "graph.h"
#include "logging.h"
#include "math.h"
#include "nauty2_8_6/gtools.h"
#include "parallel_search.h"
#include <algorithm>
#include <deque>
#include <limits>
#include <map>
#include <sstream>

namespace fpc {

Graph::Graph(std::istream &input)
    : nr_vertices_(0), nr_edges_(0), max_length_(0) {
  char dec;
  input >> dec;
  Vertex nr_edges;
  Vertex nr_vertices;
  std::string line;
  all_pair_ = true;
  max_length_ = std::numeric_limits<Edge_length>::max();
  while (!input.eof()) {
    switch (dec) {
    case 'c':
      std::getline(input, line);
      break;
    case 'p': {
      std::string str;
      input >> str >> nr_vertices >> nr_edges;
      assert(str == "edge");
      adjacency_ = std::vector<std::vector<char>>(
          nr_vertices, std::vector<char>(nr_vertices, false));
      neighbors_ =
          std::vector<std::vector<Vertex>>(nr_vertices, std::vector<Vertex>());
      break;
    }
    case 'e':
      Vertex v, w;
      input >> v >> w;
      --v;
      --w;
      if (v != w && !has_edge(v, w)) {
        add_edge(v, w);
      } else {
        nr_edges--;
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
      all_pair_ = false;
      break;
    default:
      LOG << "Invalid character " << dec << " at beginning of line."
          << std::endl;
      break;
    }
    input >> dec;
  }
  assert(nr_edges_ == nr_edges);
  assert(nr_vertices_ == nr_vertices);
  assert(adjacency_.size() > 0);
  assert(all_pair_ || terminals_.size() == 2);
}

void Graph::preprocess() {
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
    Vertex cur_articulation_point_removed = preprocess_articulation_points();
    found |= cur_articulation_point_removed > 0;
    articulation_point_removed += cur_articulation_point_removed;
    if (!found) {
      Vertex cur_position_determined = preprocess_position_determined(1);
      found |= cur_position_determined > 0;
      position_determined_removed += cur_position_determined;
    }
    if (!found) {
      preprocess_start_goal_edges();
      Vertex cur_unusable_edge_removed = preprocess_unusable_edge();
      found |= cur_unusable_edge_removed > 0;
      unusable_edge_removed += cur_unusable_edge_removed;
    }
  }
  preprocess_start_goal_edges();
  assert(all_pair_ || !has_edge(terminals_[0], terminals_[1]));
}

void Graph::weighted_preprocess() {
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
    Vertex cur_articulation_point_removed = preprocess_articulation_points();
    found |= cur_articulation_point_removed > 0;
    articulation_point_removed += cur_articulation_point_removed;
    if (!found) {
      preprocess_start_goal_edges();
      Vertex cur_unusable_edge_removed = preprocess_unusable_edge();
      found |= cur_unusable_edge_removed > 0;
      unusable_edge_removed += cur_unusable_edge_removed;
    }
  }
  preprocess_start_goal_edges();
  assert(!has_edge(terminals_[0], terminals_[1]));
}

void Graph::reduce_modulo_equivalence() {
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

void Graph::print_stats() const {
  if (isolated_removed)
    LOG << "Removed isolated: " << isolated_removed << std::endl;
  if (forwarder_removed)
    LOG << "Removed forwarder: " << forwarder_removed << std::endl;
  if (position_determined_removed)
    LOG << "Removed position determined: " << position_determined_removed
        << std::endl;
  if (articulation_point_removed)
    LOG << "Removed behind articulation point: " << articulation_point_removed
        << std::endl;
  if (unreachable_removed)
    LOG << "Removed unreachable: " << unreachable_removed << std::endl;
  if (unusable_edge_removed)
    LOG << "Removed unusable edge: " << unusable_edge_removed << std::endl;
  if (max_length_decrease)
    LOG << "Max length decreased by: " << max_length_decrease << std::endl;
  if (edge_weights_.empty()) {
    LOG << "#vertices " << nr_vertices() << " (" << equivalence_classes().size()
        << ")"
        << " #edges " << nr_edges();
  } else {

    LOG << "#vertices " << nr_vertices() << " (?)"
        << " #edges " << nr_edges();
  }
  LOG << " #paths upper bound 2^" << bit_upper_bound();
  LOG << " max. length " << static_cast<size_t>(max_length_);
  if (!all_pair_) {
    LOG << " min. length " << static_cast<size_t>(min_length()) << std::endl;
    LOG << "terminals: " << terminals_[0] << "," << terminals_[1] << std::endl;
  } else {
    LOG << std::endl;
  }
}

void Graph::print_ordered(std::ostream &os,
                          std::vector<std::size_t> const &order) const {
  os << "p edge " << nr_vertices() << " " << nr_edges() << std::endl;
  std::set<Vertex> active;
  std::vector<char> taken(nr_vertices(), false);
  std::vector<Vertex> degree(nr_vertices(), 0);
  for (Vertex v = 0; v < nr_vertices(); v++) {
    degree[v] = neighbors(v).size();
  }
  for (auto v : order) {
    active.insert(v);
    for (auto neigh : neighbors(v)) {
      if (taken[neigh]) {
        continue;
      }
      active.insert(neigh);
      os << "e " << v + 1 << " " << neigh + 1 << std::endl;
    }
    taken[v] = true;
    active.erase(v);
  }
  os << "l " << max_length() << std::endl;
  if (!is_all_pair()) {
    os << "t " << terminals().front() + 1 << " " << terminals().back() + 1
       << std::endl;
  }
}

Edge_length Graph::min_length() const {
  if (all_pair_) {
    return 1;
  }
  std::vector<Edge_length> distance_to_goal(
      adjacency_.size(), std::numeric_limits<Edge_length>::max());
  dijkstra(terminals_[1], distance_to_goal, {});
  return distance_to_goal[terminals_[0]];
}

std::size_t Graph::bit_upper_bound() const {
  if (bit_upper_bound_.has_value()) {
    return *bit_upper_bound_;
  }
  assert(edge_weights_.empty());
  std::vector<mpz_class> degrees;
  for (auto v = 0; v < nr_vertices(); ++v) {
    if (is_all_pair() || (terminals_[0] != v && terminals_[1] != v)) {
      degrees.push_back(neighbors(v).size() - 1);
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
  if (is_all_pair()) {
    mpz_class sum = 0;
    for (auto v = 0; v < nr_vertices(); ++v) {
      for (auto other = v + 1; other < nr_vertices(); ++other) {
        sum += neighbors(v).size() * neighbors(other).size();
      }
    }
    res *= sum;
  } else {
    res *= neighbors(terminals_[0]).size();
    res *= neighbors(terminals_[1]).size();
  }
  res += extra_paths_;

  std::size_t bits = 1;
  mpz_class values = 1;
  while (values < res) {
    values *= 2;
    bits += 1;
  }
  return bits;
}

void Graph::normalize(bool reorder) {
  Vertex unnamed = std::numeric_limits<Vertex>::max();
  std::vector<Vertex> new_name(adjacency_.size(), unnamed);
  Vertex cur_name = 0;
  if (reorder && !all_pair_) {
    new_name[terminals_[0]] = cur_name++;
    new_name[terminals_[1]] = cur_name++;
  }
  for (Vertex v = 0; v < adjacency_.size(); v++) {
    if (neighbors(v).empty() &&
        (all_pair_ || (v != terminals_[0] && v != terminals_[1]))) {
      continue;
    }
    if (new_name[v] == unnamed) {
      new_name[v] = cur_name++;
    }
  }
  auto new_adjacency = std::vector<std::vector<char>>(
      cur_name, std::vector<char>(cur_name, false));
  auto new_neighbors =
      std::vector<std::vector<Vertex>>(cur_name, std::vector<Vertex>());
  decltype(edge_weights_) new_arc_weights;
  for (Vertex v = 0; v < adjacency_.size(); v++) {
    if (neighbors(v).empty()) {
      continue;
    }
    for (auto neigh : neighbors(v)) {
      assert(new_name[neigh] != unnamed);
      new_neighbors[new_name[v]].push_back(new_name[neigh]);
      new_adjacency[new_name[v]][new_name[neigh]] = true;
      if (auto it = find_weight(neigh, v); it != edge_weights_.end()) {
        new_arc_weights.emplace(Edge(std::min(new_name[neigh], new_name[v]),
                                     std::max(new_name[neigh], new_name[v])),
                                it->second);
      }
    }
  }
  adjacency_ = new_adjacency;
  neighbors_ = new_neighbors;
  assert(edge_weights_.size() == new_arc_weights.size());
  edge_weights_ = new_arc_weights;
  if (!all_pair_) {
    assert(new_name[terminals_[0]] != unnamed);
    assert(new_name[terminals_[1]] != unnamed);
    terminals_[0] = new_name[terminals_[0]];
    terminals_[1] = new_name[terminals_[1]];
  }
}

sparsegraph Graph::to_canon_nauty(bool reorder) {
  assert(edge_weights_.empty());
  normalize(reorder);
  DEFAULTOPTIONS_SPARSEGRAPH(options);
  options.getcanon = true;
  options.defaultptn = false;
  SG_DECL(sg);
  int m = SETWORDSNEEDED(adjacency_.size());
  nauty_check(WORDSIZE, m, adjacency_.size(), NAUTYVERSIONID);
  size_t nr_edges = 0;
  for (Vertex v = 0; v < adjacency_.size(); v++) {
    nr_edges += neighbors(v).size();
  }
  sg.v = (edge_t *)malloc(sizeof(edge_t) * (adjacency_.size() + nr_edges) +
                          sizeof(degree_t) * adjacency_.size());
  sg.d = (degree_t *)(sg.v + adjacency_.size());
  sg.e = (edge_t *)(sg.d + adjacency_.size());
  sg.nv = adjacency_.size();
  sg.nde = nr_edges;
  sg.vlen = adjacency_.size();
  sg.dlen = adjacency_.size();
  sg.elen = nr_edges;
  nr_edges = 0;
  for (Vertex v = 0; v < adjacency_.size(); v++) {
    sg.v[v] = nr_edges;
    for (Vertex w : neighbors(v)) {
      sg.e[nr_edges++] = w;
    }
    sg.d[v] = neighbors(v).size();
  }
  int *lab = (int *)malloc(sg.nv * sizeof(int));
  int *ptn = (int *)malloc(sg.nv * sizeof(int));
  int *orbits = (int *)malloc(sg.nv * sizeof(int));
  ptn[0] = 0;
  ptn[1] = 0;
  lab[0] = 0;
  lab[1] = 1;
  for (Vertex v = 2; v < adjacency_.size(); v++) {
    ptn[v] = 1;
    lab[v] = v;
  }
  statsblk stats;
  SG_DECL(canon_sg);
  canon_sg.v =
      (edge_t *)malloc(sizeof(edge_t) * (adjacency_.size() + nr_edges) +
                       sizeof(degree_t) * adjacency_.size());
  canon_sg.d = (degree_t *)(canon_sg.v + adjacency_.size());
  canon_sg.e = (edge_t *)(canon_sg.d + adjacency_.size());
  canon_sg.nv = adjacency_.size();
  canon_sg.nde = nr_edges;
  canon_sg.vlen = sg.vlen;
  canon_sg.dlen = sg.dlen;
  canon_sg.elen = sg.elen;
  sparsenauty(&sg, lab, ptn, orbits, &options, &stats, &canon_sg);
  sortlists_sg(&canon_sg);
  free(sg.v);
  free(lab);
  free(ptn);
  free(orbits);
  return canon_sg;
}

std::vector<sparsegraph> Graph::all_pair_nauty(bool reorder) {
  assert(edge_weights_.empty());
  assert(!is_all_pair());
  normalize(reorder);
  std::vector<sparsegraph> result;
  std::vector<Vertex> all;
  for (Vertex v = 0; v < adjacency_.size(); v++) {
    all.push_back(v);
  }
  for (Vertex t1 = 0; t1 < adjacency_.size(); t1++) {
    for (Vertex t2 = t1 + 1; t2 < adjacency_.size(); t2++) {
      Graph copy = subgraph(all);
      copy.all_pair_ = false;
      copy.terminals_ = {t1, t2};
      // copy.preprocess();
      result.push_back(copy.to_canon_nauty(reorder));
    }
  }
  return result;
}

double Graph::nr_automorphisms() {
  assert(edge_weights_.empty());
  normalize(false);
  DEFAULTOPTIONS_SPARSEGRAPH(options);
  options.getcanon = true;
  options.defaultptn = false;
  SG_DECL(sg);
  int m = SETWORDSNEEDED(adjacency_.size());
  nauty_check(WORDSIZE, m, adjacency_.size(), NAUTYVERSIONID);
  size_t nr_edges = 0;
  for (Vertex v = 0; v < adjacency_.size(); v++) {
    nr_edges += neighbors(v).size();
  }
  sg.v = (edge_t *)malloc(sizeof(edge_t) * (adjacency_.size() + nr_edges) +
                          sizeof(degree_t) * adjacency_.size());
  sg.d = (degree_t *)(sg.v + adjacency_.size());
  sg.e = (edge_t *)(sg.d + adjacency_.size());
  sg.nv = adjacency_.size();
  sg.nde = nr_edges;
  sg.vlen = adjacency_.size();
  sg.dlen = adjacency_.size();
  sg.elen = nr_edges;
  nr_edges = 0;
  for (Vertex v = 0; v < adjacency_.size(); v++) {
    sg.v[v] = nr_edges;
    for (Vertex w : neighbors(v)) {
      sg.e[nr_edges++] = w;
    }
    sg.d[v] = neighbors(v).size();
  }
  int *lab = (int *)malloc(sg.nv * sizeof(int));
  int *ptn = (int *)malloc(sg.nv * sizeof(int));
  int *orbits = (int *)malloc(sg.nv * sizeof(int));
  ptn[0] = 0;
  ptn[1] = 0;
  lab[0] = 0;
  lab[1] = 1;
  for (Vertex v = 2; v < adjacency_.size(); v++) {
    ptn[v] = 1;
    lab[v] = v;
  }
  statsblk stats;
  SG_DECL(canon_sg);
  canon_sg.v =
      (edge_t *)malloc(sizeof(edge_t) * (adjacency_.size() + nr_edges) +
                       sizeof(degree_t) * adjacency_.size());
  canon_sg.d = (degree_t *)(canon_sg.v + adjacency_.size());
  canon_sg.e = (edge_t *)(canon_sg.d + adjacency_.size());
  canon_sg.nv = adjacency_.size();
  canon_sg.nde = nr_edges;
  canon_sg.vlen = sg.vlen;
  canon_sg.dlen = sg.dlen;
  canon_sg.elen = sg.elen;
  sparsenauty(&sg, lab, ptn, orbits, &options, &stats, &canon_sg);
  sortlists_sg(&canon_sg);
  free(sg.v);
  free(lab);
  free(ptn);
  free(orbits);
  free(canon_sg.v);
  return stats.grpsize1 * std::pow(10, stats.grpsize2);
}

std::vector<std::vector<Vertex>> Graph::equivalence_classes() const {
  assert(edge_weights_.empty());
  std::unordered_map<std::vector<Vertex>, std::vector<Vertex>,
                     vector_hash<Vertex>>
      false_twins;
  std::unordered_map<std::vector<Vertex>, std::vector<Vertex>,
                     vector_hash<Vertex>>
      true_twins;
  for (Vertex v = 0; v < adjacency_.size(); v++) {
    std::vector<Vertex> sorted_neighbors = neighbors(v);
    srg::sort(sorted_neighbors);
    auto [it, inserted_false] =
        false_twins.try_emplace(sorted_neighbors, std::vector<Vertex>{v});
    if (!inserted_false) {
      it->second.push_back(v);
    }
    sorted_neighbors.push_back(v);
    srg::sort(sorted_neighbors);
    auto [other_it, inserted_true] =
        true_twins.try_emplace(sorted_neighbors, std::vector<Vertex>{v});
    if (!inserted_true) {
      other_it->second.push_back(v);
    }
    assert(inserted_false || inserted_true);
  }
  std::vector<std::vector<Vertex>> rv;
  std::vector<char> seen(adjacency_.size(), false);
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
  for (Vertex v = 0; v < adjacency_.size(); v++) {
    if (!seen[v]) {
      rv.push_back({v});
    }
  }

  return rv;
}

Graph::Graph(Vertex n)
    : nr_vertices_(0), nr_edges_(0), max_length_(-1),
      neighbors_(n, std::vector<Vertex>()),
      adjacency_(n, std::vector<char>(n, false)) {}

Graph Graph::subgraph(std::vector<Vertex> restrict_to) const {
  Graph ret(restrict_to.size());
  ret.max_length_ = max_length_;
  auto unnamed = std::numeric_limits<Vertex>::max();
  std::vector<Vertex> new_name(adjacency_.size(), unnamed);
  Vertex cur_name = 0;
  for (Vertex v : restrict_to) {
    new_name[v] = cur_name++;
  }
  for (Vertex v : restrict_to) {
    for (Vertex neigh : neighbors(v)) {
      if (new_name[neigh] == unnamed) {
        continue;
      }
      if (neigh < v) {
        continue;
      }
      ret.add_edge(new_name[v], new_name[neigh]);

      if (auto it = find_weight(v, neigh); it != edge_weights_.end()) {
        ret.edge_weights_.try_emplace(
            Edge(std::min(new_name[v], new_name[neigh]),
                 std::max(new_name[v], new_name[neigh])),
            it->second);
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
  assert(ret.count_edges() == ret.nr_edges_);
  assert(ret.adjacency_.size() > 0);
  return ret;
}

Graph Graph::copy() const {
  std::vector<Vertex> all;
  all.reserve(adjacency_.size());
  for (Vertex v = 0; v < adjacency_.size(); v++) {
    all.push_back(v);
  }
  auto rv = subgraph(all);
  rv.extra_paths_ = extra_paths_;
  return rv;
}

std::unordered_set<Vertex>
Graph::must_use_subgraph(Vertex must_use,
                         std::vector<Edge_length> const &distance_from_start,
                         std::vector<Edge_length> const &distance_to_goal,
                         Edge_length backward_allowance) const {
  assert(edge_weights_.empty());
  std::unordered_set<Vertex> restrict_to;
  auto find_usable = [&restrict_to, must_use, backward_allowance,
                      this](auto const &distance) {
    std::vector<Edge_length> min_allowance_used(
        adjacency_.size(), std::numeric_limits<Edge_length>::max());
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
      for (auto next : neighbors(cur_vertex)) {
        if (distance[next] < distance[cur_vertex] &&
            cur_allowance < min_allowance_used[next]) {
          min_allowance_used[next] = cur_allowance;
          fifo.push_back(next);
        } else if (auto new_allowance = cur_allowance + distance[next] -
                                        distance[cur_vertex] + 1;
                   new_allowance < min_allowance_used[next] &&
                   new_allowance <= backward_allowance) {
          min_allowance_used[next] = new_allowance;
          fifo.push_back(next);
        }
      }
    }
    for (Vertex v = 0; v < adjacency_.size(); v++) {
      if (min_allowance_used[v] <= backward_allowance) {
        restrict_to.insert(v);
      }
    }
  };
  find_usable(distance_from_start);
  find_usable(distance_to_goal);
  return restrict_to;
}

bool Graph::is_length_limited() const {
  if (is_length_limited_.has_value()) {
    return *is_length_limited_;
  }
  assert(edge_weights_.empty());
  return max_length_ + 1 < nr_vertices();
}

void Graph::dijkstra(Vertex start, std::vector<Edge_length> &distance,
                     const std::set<Vertex> &forbidden) const {
  std::deque<Vertex> queue;
  queue.push_back(start);
  distance[start] = 0;
  while (!queue.empty()) {
    auto cur_vertex = queue.front();
    auto cur_cost = distance[cur_vertex];
    queue.pop_front();
    for (auto w : neighbors(cur_vertex)) {
      auto extra = edge_length(cur_vertex, w);
      if (cur_cost + extra >= distance[w]) {
        continue;
      }
      if (forbidden.contains(w)) {
        continue;
      }
      if (cur_cost + extra < distance[w]) {
        distance[w] = extra + cur_cost;
        queue.push_back(w);
      }
    }
  }
}

void Graph::contract_edges(Vertex v) {
  for (auto in : neighbors(v)) {
    for (auto out : neighbors(v)) {
      if (in <= out) {
        continue;
      }
      if (!has_edge(in, out)) {
        edge_weights_.try_emplace(
            Edge(std::min(in, out), std::max(in, out)),
            Limited_count<mpz_class>(max_length_, max_length_, {}));
        add_edge(in, out);
      } else {
        edge_weights_.try_emplace(Edge(std::min(in, out), std::max(in, out)),
                                  Limited_count<mpz_class>(max_length_, 1));
      }
      Limited_count<mpz_class> in_to_v(max_length_, 1);
      if (auto it = find_weight(in, v); it != edge_weights_.end()) {
        in_to_v = it->second;
      }
      Limited_count<mpz_class> v_to_out(max_length_, 1);
      if (auto it = find_weight(v, out); it != edge_weights_.end()) {
        v_to_out = it->second;
      }
      in_to_v *= v_to_out;
      edge_weights_.at(Edge(std::min(in, out), std::max(in, out))) += in_to_v;
    }
  }
  remove_vertex(v);
}

void Graph::contract_terminal_edges(Vertex v, Vertex terminal) {
  for (auto out : neighbors(v)) {
    if (terminal == out) {
      continue;
    }
    if (!has_edge(terminal, out)) {
      edge_weights_.try_emplace(
          Edge(std::min(terminal, out), std::max(terminal, out)),
          Limited_count<mpz_class>(max_length_, max_length_, {}));
      add_edge(terminal, out);
    } else {
      edge_weights_.try_emplace(
          Edge(std::min(terminal, out), std::max(terminal, out)),
          Limited_count<mpz_class>(max_length_, 1));
    }
    Limited_count<mpz_class> in_to_v(max_length_, 1);
    if (auto it = find_weight(terminal, v); it != edge_weights_.end()) {
      in_to_v = it->second;
    }
    Limited_count<mpz_class> v_to_out(max_length_, 1);
    if (auto it = find_weight(v, out); it != edge_weights_.end()) {
      v_to_out = it->second;
    }
    in_to_v *= v_to_out;
    edge_weights_.at(Edge(std::min(terminal, out), std::max(terminal, out))) +=
        in_to_v;
  }
  remove_vertex(v);
}

Vertex Graph::preprocess_start_goal_edges() {
  if (all_pair_) {
    return 0;
  }
  assert(terminals_.size() == 2);
  Vertex found = 0;
  if (has_edge(terminals_[0], terminals_[1])) {
    found = 1;
    if (auto it = find_weight(terminals_[0], terminals_[1]);
        it != edge_weights_.end()) {
      extra_paths_ += it->second.total_count(max_length_);
    } else {
      extra_paths_ += 1;
    }
    remove_edge(terminals_[0], terminals_[1]);
  }
  return found;
}

Vertex Graph::preprocess_start_goal_forwarders() {
  assert(edge_weights_.empty());
  if (all_pair_) {
    return 0;
  }
  Vertex found = 0;
  std::vector<char> visited(adjacency_.size(), false);
  visited[terminals_[0]] = true;
  std::vector<char> on_stack(adjacency_.size(), false);
  std::vector<Vertex> stack = neighbors(terminals_[0]);
  for (auto neigh : neighbors(terminals_[0])) {
    on_stack[neigh] = true;
  }
  while (!stack.empty()) {
    auto cur = stack.back();
    stack.pop_back();
    on_stack[cur] = false;
    std::optional<Vertex> non_visited;
    bool found_two = false;
    bool found_goal = false;
    for (auto neigh : neighbors(cur)) {
      if (neigh == terminals_[1]) {
        // TODO: count the paths
        found_goal = true;
        continue;
      }
      if (!visited[neigh]) {
        found_two = !!non_visited;
        non_visited = neigh;
      }
    }
    if (found_two) {
      continue;
    }
    if (found_goal && !non_visited) {
      found++;
      remove_vertex(cur);
      continue;
    }
    assert(non_visited);
    visited[cur] = true;
    if (!on_stack[*non_visited]) {
      stack.push_back(*non_visited);
      on_stack[*non_visited] = true;
    }
  }
  return found;
}

Vertex Graph::preprocess_isolated() {
  if (all_pair_) {
    return 0;
  }
  assert(terminals_.size() == 2);
  Vertex found = 0;
  for (size_t v = 0; v < adjacency_.size(); v++) {
    if (terminals_[0] == v || terminals_[1] == v) {
      continue;
    }
    if (neighbors(v).size() == 1) {
      found++;
      remove_vertex(v);
    }
  }

  auto do_terminal = [&](Vertex &t) {
    if (neighbors(t).size() != 1) {
      return 0;
    }
    Vertex neighbor = neighbors(t).front();
    if (neighbor == terminals_[0] || neighbor == terminals_[1]) {
      preprocess_start_goal_edges();
      return 1;
    }
    assert(max_length_ > 0);
    found++;
    if (auto it = find_weight(t, neighbor); it != edge_weights_.end()) {
      while (it->second.offset() > 0) {
        max_length_--;
        it->second.decrement_offset();
      }
      contract_terminal_edges(neighbor, t);
    } else {
      max_length_--;
      remove_vertex(t);
      t = neighbor;
    }
    return 1;
  };

  found += do_terminal(terminals_[0]);
  found += do_terminal(terminals_[1]);
  return found;
}

Vertex Graph::preprocess_unreachable() {
  if (all_pair_) {
    return 0;
  }
  assert(terminals_.size() == 2);
  std::vector<Edge_length> distance_from_start(
      adjacency_.size(), std::numeric_limits<Edge_length>::max());
  dijkstra(terminals_[0], distance_from_start, {});

  std::vector<Edge_length> distance_to_goal(
      adjacency_.size(), std::numeric_limits<Edge_length>::max());
  dijkstra(terminals_[1], distance_to_goal, {});

  Vertex found = 0;
  for (Vertex v = 0; v < adjacency_.size(); v++) {
    if (neighbors(v).empty()) {
      continue;
    }
    if (distance_from_start[v] + distance_to_goal[v] > max_length_) {
      found++;
      remove_vertex(v);
    }
  }
  return found;
}

Vertex Graph::preprocess_forwarder() {
  assert(edge_weights_.empty());
  if (all_pair_) {
    return 0;
  }
  Vertex found = 0;
  for (size_t v = 0; v < adjacency_.size(); v++) {
    if (terminals_[0] == v || terminals_[1] == v) {
      // we cannot remove terminals this way
      continue;
    }
    if (neighbors(v).size() == 2) {
      Vertex w1 = neighbors(v).front(), w2 = neighbors(v).back();
      // currently not counting edge multiplicity
      if (has_edge(w1, w2))
        continue;

      found++;
      add_edge(w1, w2);
      remove_vertex(v);
    }
  }
  return found;
}

Vertex Graph::weighted_preprocess_forwarder() {
  if (all_pair_) {
    return 0;
  }
  Vertex found = 0;
  for (size_t v = 0; v < adjacency_.size(); v++) {
    if (terminals_[0] == v || terminals_[1] == v) {
      // we cannot remove terminals this way
      continue;
    }
    if (neighbors(v).size() == 2) {
      found++;
      contract_edges(v);
    }
  }
  return found;
}

Vertex Graph::preprocess_unusable_edge() {
  if (all_pair_) {
    return 0;
  }
  if (nr_edges_ > 10000) {
    return 0;
  }
  assert(terminals_.size() == 2);
  std::vector<std::vector<Edge_length>> distances_from_start(
      adjacency_.size(),
      std::vector<Edge_length>(adjacency_.size(),
                               std::numeric_limits<Edge_length>::max()));
  std::vector<std::vector<Edge_length>> distances_to_goal(
      adjacency_.size(),
      std::vector<Edge_length>(adjacency_.size(),
                               std::numeric_limits<Edge_length>::max()));
  Vertex found = 0;
  for (Vertex v = 0; v < adjacency_.size(); v++) {
    if (neighbors(v).empty()) {
      continue;
    }
    dijkstra(terminals_[0], distances_from_start[v], {v});
    dijkstra(terminals_[1], distances_to_goal[v], {v});
  }
  for (Vertex v = 0; v < adjacency_.size(); v++) {
    std::vector<Vertex> remove_completely;
    for (auto w : neighbors(v)) {
      Edge_length min_without_edge =
          std::min(distances_from_start[v][w] + distances_to_goal[w][v],
                   distances_from_start[w][v] + distances_to_goal[v][w]);
      if (edge_length(v, w) + min_without_edge > max_length_) {
        found++;
        remove_completely.push_back(w);
      }
    }
    for (auto w : remove_completely) {
      remove_edge(v, w);
    }
  }
  return found;
}

Vertex Graph::preprocess_articulation_points() {
  if (all_pair_) {
    return 0;
  }
  std::vector<Vertex> ap_disc(neighbors_.size(), 0);
  std::vector<Vertex> ap_low(neighbors_.size(), 0);
  std::vector<char> ap_visited(neighbors_.size(), false);
  int time = 0;
  return ap_util(terminals_[0], ap_visited, ap_disc, ap_low, time, -1).second;
}

std::pair<bool, Vertex> Graph::ap_util(Vertex u, std::vector<char> &visited,
                                       std::vector<Vertex> &disc,
                                       std::vector<Vertex> &low, int &time,
                                       int parent) {

  // Count of children in DFS Tree
  std::vector<std::pair<Vertex, bool>> children;

  // Mark the current node as visited
  visited[u] = true;

  bool found_elsewhere = false;
  Vertex eliminated_elsewhere = 0;

  // Initialize discovery time and low value
  disc[u] = low[u] = ++time;
  // Go through all vertices adjacent to this
  auto neighs = neighbors(u); // copy because we invalidate the vector
  for (auto v : neighs) {
    // If v is not visited yet, then make it a child of u
    // in DFS tree and recur for it
    if (!visited[v]) {
      auto [found_here, eliminated_here] =
          ap_util(v, visited, disc, low, time, u);
      children.emplace_back(v, found_here);
      found_elsewhere |= found_here;
      eliminated_elsewhere += eliminated_here;

      // Check if the subtree rooted with v has
      // a connection to one of the ancestors of u
      low[u] = std::min(low[u], low[v]);

      // If u is not root and low value of one of
      // its child is more than discovery value of u.
      if (parent != -1 && low[v] >= disc[u]) {
        // AP
        if (!found_here) {
          eliminated_elsewhere += prune_util(u, v);
        }
      }
    } else if (v != parent) {
      low[u] = std::min(low[u], disc[v]);
    }
  }

  // If u is root of DFS tree and has two or more children.
  if (parent == -1 && children.size() > 1) {
    for (auto [child, found_goal] : children) {
      if (found_goal) {
        continue;
      }
      eliminated_elsewhere += prune_util(u, child);
    }
  }
  return std::make_pair(found_elsewhere || u == terminals_[1],
                        eliminated_elsewhere);
}

Vertex Graph::prune_util(Vertex ap, Vertex prune) {
  Vertex found = 1;
  auto neighs = neighbors(prune);
  assert(prune != terminals_[0]);
  assert(prune != terminals_[1]);
  remove_vertex(prune);
  for (auto v : neighs) {
    if (v == ap) {
      continue;
    }
    found += prune_util(ap, v);
  }
  return found;
}

Vertex Graph::preprocess_position_determined(Edge_length budget) {
  assert(edge_weights_.empty());
  if (terminals_.empty()) {
    return 0;
  }
  std::vector<Edge_length> distance_from_start(
      adjacency_.size(), std::numeric_limits<Edge_length>::max());
  dijkstra(terminals_[0], distance_from_start, {});

  std::vector<Edge_length> distance_to_goal(
      adjacency_.size(), std::numeric_limits<Edge_length>::max());
  dijkstra(terminals_[1], distance_to_goal, {});
  if (max_length_ <= budget) {
    return 0;
  }
  Vertex found = 0;
  for (Vertex v = 0; v < adjacency_.size(); v++) {
    if (neighbors(v).empty()) {
      continue;
    }
    if (terminals_[0] == v || terminals_[1] == v) {
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
      search_graph.terminals_ = {Vertex(restrict_to.size() - 2),
                                 Vertex(restrict_to.size() - 1)};
      // TODO: use more threads
      auto solver = make_solver<ParallelSearch>(
          search_graph, Vertex(restrict_to.size() - 3), 1);
      extra_paths_ += solver->search();
      remove_vertex(v);
      break;
    }
  }
  return found;
}
} // namespace fpc
