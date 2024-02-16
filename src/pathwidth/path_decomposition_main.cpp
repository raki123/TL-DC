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

#include "heuristics.hpp"
#include "path_decomposer.hpp"
#include <iostream>
#include <optional>
#include <string.h>
#include <unordered_map>

namespace {
std::string strip(std::string in) {
  std::string res;
  for (int i = 0; i < in.length(); i++) {
    if (!isspace(in[i]))
      res += in[i];
  }
  return res;
}

} // namespace

int main(int argc, char **argv) {

  std::unordered_map<size_t, size_t> vertex_remap;
  std::vector<size_t> reverse_map;
  std::vector<std::vector<std::size_t>> neighbors;
  size_t cur_v = 0;
  auto get_vertex = [&](auto orig) {
    if (auto it = vertex_remap.find(orig); it != vertex_remap.end()) {
      return it->second;
    }
    neighbors.push_back({});
    reverse_map.push_back(orig);
    return vertex_remap[orig] = cur_v++;
  };
  std::optional<std::size_t> first;
  std::optional<std::size_t> last;

  size_t arg_idx = 0;
  while (++arg_idx < argc) {
    if (strcmp(argv[arg_idx], "-f") == 0 ||
        strcmp(argv[arg_idx], "--first") == 0) {
      arg_idx++;
      if (arg_idx == argc) {
        std::cerr << "Not enough arguments!" << std::endl;
        return -1;
      }
      first = get_vertex(std::stoi(argv[arg_idx]));
    } else if (strcmp(argv[arg_idx], "-l") == 0 ||
               strcmp(argv[arg_idx], "--last") == 0) {
      arg_idx++;
      if (arg_idx == argc) {
        std::cerr << "Not enough arguments!" << std::endl;
        return -1;
      }
      last = get_vertex(std::stoi(argv[arg_idx]));
    } else {
      std::cerr << "Unknown argument " << argv[arg_idx] << std::endl;
      return -1;
    }
  }

  std::string line;
  while (!std::cin.eof()) {
    std::getline(std::cin, line);
    if (line.size() == 0) {
      continue;
    }
    line = strip(line);
    auto idx = line.find(',');
    auto v_str = line.substr(2, idx - 2);
    auto v = get_vertex(stoi(v_str));
    auto w_str = line.substr(idx + 1, line.size() - idx - 2);
    auto w = get_vertex(stoi(w_str));
    neighbors[v].push_back(w);
    neighbors[w].push_back(v);
  }

  fpt::PathDecomposer<fpt::H1> decomposer(neighbors, first, last);

  auto [min_width, best_order] = [&]() {
    auto width = decomposer.make_order();
    return std::make_pair(width, decomposer.get_order());
  }();
  for (auto i = 0; i < 100; i++) {
    auto width = decomposer.make_probabilistic_order();
    if (width < min_width) {
      best_order = decomposer.get_order();
    }
  }
  for (auto v : best_order) {
    std::cout << reverse_map[v] << " ";
  }
  std::cout << std::endl;
}
