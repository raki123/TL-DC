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

#pragma once
#include "clhash/clhash.h"
extern "C" {
#include "nauty2_8_6/nausparse.h"
}
#include <vector>

namespace fpc {

extern clhasher hasher__;

// vectors
template <typename T> struct vector_hash {
  using is_avalanching = void;
  size_t operator()(const std::vector<T> &key) const { return hasher__(key); }
};

// nauty sparsegraphs
typedef sparsegraph NPCacheKey;
struct sg_hash {
  size_t operator()(NPCacheKey const &key) const {
    return hasher__((const char *)key.v, sizeof(edge_t) * (key.nv + key.nde) +
                                             sizeof(degree_t) * key.nv);
  }
};
struct sg_equal {
  bool operator()(NPCacheKey const &lhs, NPCacheKey const &rhs) const {
    return lhs.nv == rhs.nv && lhs.nde == rhs.nde &&
           std::memcmp(lhs.v, rhs.v,
                       sizeof(edge_t) * (lhs.nv + lhs.nde) +
                           sizeof(degree_t) * lhs.nv) == 0;
  }
};

// undirected fbs
typedef uint8_t frontier_index_t;
typedef std::vector<frontier_index_t> Frontier;
struct nFrontier {
  nFrontier(std::size_t size, frontier_index_t idx) : data_(size, idx){};

  frontier_index_t &operator[](frontier_index_t idx) { return data_[idx]; }
  frontier_index_t const &operator[](frontier_index_t idx) const {
    return data_[idx];
  }

  using Iter = std::vector<frontier_index_t>::iterator;
  using ConstIter = std::vector<frontier_index_t>::const_iterator;

  Iter begin() { return data_.begin(); }
  ConstIter begin() const { return data_.begin(); }
  Iter end() { return data_.end(); }
  ConstIter end() const { return data_.end(); }
  std::size_t size() const { return data_.size(); }

  void rehash() { hash_ = vector_hash<frontier_index_t>{}(data_); }

private:
  friend struct nFrontierHash;
  friend bool operator==(nFrontier const &, nFrontier const &);

  std::size_t hash_ = 0;
  std::vector<frontier_index_t> data_;
};

struct nFrontierHash {
  using is_avalanching = void;
  size_t operator()(nFrontier const &key) const { return key.hash_; }
};

inline bool operator==(nFrontier const &lhs, nFrontier const &rhs) {
  return lhs.data_ == rhs.data_;
}

// directed fbs
enum directed_frontier_index_t : uint8_t {
  no_edge_index = 0b0'1111110,
  two_edge_index = 0b1'1111110,
  outgoing_index = 0b0'1111111,
  incoming_index = 0b1'1111111,
  incoming_flag = 0b1'0000000
};

inline directed_frontier_index_t operator&(directed_frontier_index_t lhs,
                                           directed_frontier_index_t rhs) {
  return (directed_frontier_index_t)(static_cast<std::underlying_type_t<
                                         directed_frontier_index_t>>(lhs) &
                                     static_cast<std::underlying_type_t<
                                         directed_frontier_index_t>>(rhs));
}
inline directed_frontier_index_t operator+(directed_frontier_index_t lhs,
                                           directed_frontier_index_t rhs) {
  return (directed_frontier_index_t)(static_cast<std::underlying_type_t<
                                         directed_frontier_index_t>>(lhs) +
                                     static_cast<std::underlying_type_t<
                                         directed_frontier_index_t>>(rhs));
}
inline directed_frontier_index_t operator|(directed_frontier_index_t lhs,
                                           directed_frontier_index_t rhs) {
  return (directed_frontier_index_t)(static_cast<std::underlying_type_t<
                                         directed_frontier_index_t>>(lhs) |
                                     static_cast<std::underlying_type_t<
                                         directed_frontier_index_t>>(rhs));
}
inline directed_frontier_index_t operator^(directed_frontier_index_t lhs,
                                           directed_frontier_index_t rhs) {
  return (directed_frontier_index_t)(static_cast<std::underlying_type_t<
                                         directed_frontier_index_t>>(lhs) ^
                                     static_cast<std::underlying_type_t<
                                         directed_frontier_index_t>>(rhs));
}

static_assert(sizeof(directed_frontier_index_t) == 1);
struct DirectedFrontier {
  DirectedFrontier(std::size_t size, directed_frontier_index_t idx)
      : data_(size, idx){};

  directed_frontier_index_t &operator[](directed_frontier_index_t idx) {
    return (idx & incoming_flag ? data_[idx ^ incoming_flag] : data_[idx]);
  }
  directed_frontier_index_t const &
  operator[](directed_frontier_index_t idx) const {
    return (idx & incoming_flag ? data_[idx ^ incoming_flag] : data_[idx]);
  }

  using Iter = std::vector<directed_frontier_index_t>::iterator;
  using ConstIter = std::vector<directed_frontier_index_t>::const_iterator;

  Iter begin() { return data_.begin(); }
  ConstIter begin() const { return data_.begin(); }
  Iter end() { return data_.end(); }
  ConstIter end() const { return data_.end(); }
  std::size_t size() const { return data_.size(); }

  void rehash() { hash_ = vector_hash<directed_frontier_index_t>{}(data_); }

private:
  friend struct DirectedFrontierHash;
  friend bool operator==(DirectedFrontier const &, DirectedFrontier const &);

  std::size_t hash_;
  std::vector<directed_frontier_index_t> data_;
};

struct DirectedFrontierHash {
  using is_avalanching = void;
  size_t operator()(DirectedFrontier const &key) const { return key.hash_; }
};

inline bool operator==(DirectedFrontier const &lhs,
                       DirectedFrontier const &rhs) {
  return lhs.data_ == rhs.data_;
}

// nauty_pathwidth_search
typedef std::pair<sparsegraph, Frontier> TWCacheKey;
struct twc_hash {
  size_t operator()(TWCacheKey const &key) const {
    return sg_hash{}(key.first);
  }
};
struct twc_equal {
  bool operator()(TWCacheKey const &lhs, TWCacheKey const &rhs) const {
    return sg_equal{}(lhs.first, rhs.first);
  }
};
} // namespace fpc
