// TLDC "Too long; Didn't Count" A length limited path counter.
// Copyright (C) 2023 Rafael Kiesel, Markus Hecher

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
#include "nauty_pathwidth_search.h"
#include "parallel_search.h"
#include "treewidth_search.h"

#include <sys/resource.h>

#include "Decomposer.hpp"
#include <iostream>

void print_help() {
  std::cerr
      << "TLDC \"Too long; Didn't Count\" A length limited path counter.\n"
      << "Copyright (C) 2023 Rafael Kiesel, Markus Hecher\n"
      << "License GPLv3+: GNU GPL version 3 or later "
         "<http://gnu.org/licenses/gpl.html>\n"
      << "\n"
      << "Usage: tldc [-s .] [-t .] [-m .] [-h]\n"
      << "       -s         --strategy STRATEGY     use STRATEGY.\n"
      << "                                          Must be one of:\n"
      << "                                          * auto: choose "
         "automatically (default)\n"
      << "                                          * pathfbs: frontier-based "
         "search along a path decomposition\n"
      << "                                          * nautyfbs: frontier-based "
         "search along a path decomposition with modulo automorphisms\n"
      << "                                          * treefbs: frontier-based "
         "search along a tree decomposition\n"
      << "                                          * nautydfs: depth first "
         "search with modulo automorphisms\n"
      << "       -t         --threads NUM           use NUM threads.\n"
      << "       -m         --memory_limit LIMIT    set a hardlimit of LIMIT "
         "MB memory.\n"
      << "       -h         --help                  print this help.\n";
}

void details(size_t bags, size_t nr_bags, size_t joinch, size_t max_join,
             size_t nr_joins) {
  std::cerr << "max bag size " << bags << ", nr bags " << nr_bags
            << ", max join child bags " << joinch << ", max join bag size "
            << max_join << ", nr of joins " << nr_joins << std::endl;
}

std::ostream &operator<<(std::ostream &o, const unsigned __int128 &x) {
  if (x < 10)
    return o << (char)(x + '0');
  return o << x / 10 << (char)(x % 10 + '0');
}

enum SolverStrategy { AUTO, PATH_FBS, NAUTY_FBS, TREE_FBS, NAUTY_DFS };

int main(int argc, char **argv) {

  SolverStrategy strategy = SolverStrategy::AUTO;
  size_t number_threads = omp_get_max_threads();
  Decomposer decomposer;
  AnnotatedDecomposition decomposition;

  size_t max_bagsize, nr_bags;
  size_t t_max_bagsize, max_join_child, max_join, nr_joins, t_nr_bags;

  size_t arg_idx = 0;

  while (++arg_idx < argc) {
    if (strcmp(argv[arg_idx], "-s") == 0 ||
        strcmp(argv[arg_idx], "--strategy") == 0) {
      assert(strategy == SolverStrategy::AUTO);
      arg_idx++;
      if (arg_idx == argc) {
        std::cerr << "Not enough arguments!" << std::endl;
        return -1;
      }
      if (strcmp(argv[arg_idx], "pathfbs") == 0) {
        strategy = SolverStrategy::PATH_FBS;
      } else if (strcmp(argv[arg_idx], "nautyfbs") == 0) {
        strategy = SolverStrategy::NAUTY_FBS;
      } else if (strcmp(argv[arg_idx], "treefbs") == 0) {
        strategy = SolverStrategy::TREE_FBS;
      } else if (strcmp(argv[arg_idx], "nautydfs") == 0) {
        strategy = SolverStrategy::NAUTY_DFS;
      } else if (strcmp(argv[arg_idx], "auto") == 0) {
        strategy = SolverStrategy::AUTO;
      } else {
        std::cerr << "Unknown solver strategy " << argv[arg_idx] << std::endl;
        return -1;
      }
    } else if (strcmp(argv[arg_idx], "-m") == 0 ||
               strcmp(argv[arg_idx], "--memory_limit") == 0) {
      arg_idx++;
      if (arg_idx == argc) {
        std::cerr << "Not enough arguments!" << std::endl;
        return -1;
      }
      std::size_t mem_limit = atol(argv[arg_idx]) * 1024 * 1024;
      struct rlimit64 lim;
      getrlimit64(RLIMIT_AS, &lim);
      lim.rlim_cur = std::min(mem_limit, lim.rlim_cur);
      if (setrlimit64(RLIMIT_AS, &lim) == -1) {
        perror(strerror(errno));
        return -1;
      }
    } else if (strcmp(argv[arg_idx], "-t") == 0 ||
               strcmp(argv[arg_idx], "--threads") == 0) {
      arg_idx++;
      if (arg_idx == argc) {
        std::cerr << "Not enough arguments!" << std::endl;
        return -1;
      }
      number_threads = atol(argv[arg_idx]);
    } else if (strcmp(argv[arg_idx], "-h") == 0 ||
               strcmp(argv[arg_idx], "--help") == 0) {
      print_help();
      return 0;
    } else {
      std::cerr << "Unknown argument " << argv[arg_idx] << std::endl;
      print_help();
      return -1;
    }
  }

  fpc::Graph initial_graph(std::cin);
  initial_graph.preprocess();
  initial_graph.normalize();
  initial_graph.print_stats();

  if (strategy == SolverStrategy::AUTO) {
    // portfolio parameters
    double nr_automporphisms = initial_graph.nr_automorphisms();
    Edge_length max_length = initial_graph.max_length();

    // path
    auto rp =
        decomposer.tree_decompose(initial_graph, true, &max_bagsize, &nr_bags);

    std::cerr << "Path decomposition: ";
    details(max_bagsize, nr_bags, 0, 0, 0);

    // tree
    auto rt = decomposer.tree_decompose(initial_graph, false, &t_max_bagsize,
                                        &t_nr_bags, &max_join_child, &max_join,
                                        &nr_joins);

    std::cerr << "Tree decomposition: ";
    details(t_max_bagsize, t_nr_bags, max_join_child, max_join, nr_joins);

    bool is_all_pair = initial_graph.is_all_pair(), use_pw = false;

    // decide whether we prefer a path decomposition over a tree decomposition.
    // tree decompositions perform worse at the same width since joins are
    // quadratic in the number of entries so if we have low enough width the
    // path decomposition will be faster all pair instances have more entries
    // therefore we prefer path decompositions for higher widths here
    if ((!is_all_pair && max_bagsize <= 18) ||
        (is_all_pair && max_bagsize <= 22)) {
      use_pw = true;
    }

    // case one pair:
    // we can increase our performance on high (but doable) path width
    // if there are many automorphisms by doing caching modulo automorphism
    // if there are too many automorphisms though computing a canonical
    // representation becomes expensive case all pair: if we want to use path
    // width but it is high we have a better chance when we reduce the number of
    // entries by caching modulo automorphism
    if ((!is_all_pair && max_bagsize >= 16 && max_bagsize <= 25 &&
         nr_automporphisms >= 1e+11 && nr_automporphisms <= 1e+45) ||
        (is_all_pair && use_pw && max_bagsize >= 21)) {
      decomposition = rp;
      strategy = SolverStrategy::NAUTY_FBS;
      std::cerr << "Using strategy nautyfbs.\n";
    }
    // case all pair:
    // we have to use tree width search because parallel search only works for
    // one pair instances case use pw: we want to exploit low path width again
    // this is our last chance to do so case high length: when the length is
    // high we cannot prune based on it this is highly relevant for the last
    // solver therefore if the width and join width is reasonable use a tree
    // decomposition case low tree width: tree width and join size are low we
    // have a good chance using treewidth based search
    else if (is_all_pair || use_pw ||
             (max_length >= 60 &&
              std::max(max_join_child, t_max_bagsize) <= 22) ||
             std::max(max_join_child, t_max_bagsize) <= 20) {
      decomposition = use_pw ? rp : rt;
      strategy = use_pw ? SolverStrategy::PATH_FBS : SolverStrategy::TREE_FBS;
      std::cerr << "Using strategy " << (use_pw ? "pathfbs.\n" : "treefbs.\n");
    }
    // last resort, when none of the other cases we are good at trigger
    // do backtracking search with caching modulo automorphisms
    else {
      strategy = SolverStrategy::NAUTY_DFS;
      std::cerr << "Using strategy nautydfs.\n";
    }
  } else if (strategy == SolverStrategy::PATH_FBS ||
             strategy == SolverStrategy::NAUTY_FBS) {
    decomposition =
        decomposer.tree_decompose(initial_graph, true, &max_bagsize, &nr_bags);
    std::cerr << "Path decomposition: ";
    details(max_bagsize, nr_bags, 0, 0, 0);
  } else if (strategy == SolverStrategy::TREE_FBS) {
    decomposition = decomposer.tree_decompose(
        initial_graph, false, &t_max_bagsize, &t_nr_bags, &max_join_child,
        &max_join, &nr_joins);
    std::cerr << "Tree decomposition: ";
    details(t_max_bagsize, t_nr_bags, max_join_child, max_join, nr_joins);
  }

  std::vector<Edge_weight> res;
  switch (strategy) {
  case SolverStrategy::PATH_FBS:
  case SolverStrategy::TREE_FBS: {
    fpc::TreewidthSearch search(initial_graph, decomposition, number_threads);
    res = search.search();
    std::cerr << std::endl;
    search.print_stats();
  } break;

  case SolverStrategy::NAUTY_FBS: {
    fpc::NautyPathwidthSearch search(initial_graph, decomposition,
                                     number_threads);
    res = search.search();
    std::cerr << std::endl;
    search.print_stats();
  } break;

  case SolverStrategy::NAUTY_DFS: {
    std::vector<sparsegraph> sgs;
    if (initial_graph.is_all_pair()) {
      sgs = initial_graph.all_pair_nauty(true);
    } else {
      sgs.push_back(initial_graph.to_canon_nauty(true));
    }
    fpc::ParallelSearch search(initial_graph.to_canon_nauty(true),
                               initial_graph.max_length(), number_threads);
    search.add_to_initial(sgs);
    res = search.search();
    std::cerr << std::endl;
    search.print_stats();
  } break;

  default:
    assert(false);
  }

  fpc::Edge_weight final_result = 0;
  for (fpc::Edge_length l = 0; l < res.size(); l++) {
    final_result += res[l];
  }
  res = initial_graph.extra_paths();
  for (fpc::Edge_length l = 0; l < res.size(); l++) {
    final_result += res[l];
  }
  std::cout << final_result << std::endl;
  return 0;
}
