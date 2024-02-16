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

#include <gmpxx.h>
#include <boost/multiprecision/cpp_int.hpp>
#include <sstream>
#include <utility>
#include <vector>

template <typename T>
inline void hash_combine(std::size_t &seed, T const &key) {
  std::hash<T> hasher;
  seed ^= hasher(key) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

namespace std {
template <typename T1, typename T2> struct hash<std::pair<T1, T2>> {
  std::size_t operator()(std::pair<T1, T2> const &p) const {
    std::size_t seed(0);
    hash_combine(seed, p.first);
    hash_combine(seed, p.second);
    return seed;
  }
};
} // namespace std

namespace fpc {

typedef uint32_t Vertex;
typedef std::pair<Vertex, Vertex> Edge;
typedef uint16_t Edge_length;
typedef uint64_t uint64_t;
typedef boost::multiprecision::uint128_t uint128_t;
typedef boost::multiprecision::uint256_t uint256_t;
typedef boost::multiprecision::uint512_t uint512_t;
typedef boost::multiprecision::uint1024_t uint1024_t;
typedef mpz_class Edge_weight;

size_t get_offset(std::vector<uint64_t> const &result);
size_t get_offset(std::vector<mpz_class> const &result);

template <typename count_t> inline mpz_class to_mpz(count_t const &count) {
  std::stringstream ss;
  ss << count;
  mpz_class val(ss.str());
  return val;
};
} // namespace fpc