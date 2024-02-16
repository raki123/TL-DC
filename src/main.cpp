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

#include "annotated_decomposition.hpp"
#include "branch_decomposition.hpp"
#include "digraph.h"
#include "graph.h"
#include "logging.h"
#include "make_decomposition.h"
#include "math.h"
#include "meet_in_the_middle_directed_search.h"
#include "meet_in_the_middle_search.h"
#include "nauty_pathwidth_search.h"
#include "parallel_directed_search.h"
#include "parallel_nauty_search.h"
#include "parallel_search.h"
#include "tree_decomposition.hpp"
#include "treewidth_directed_search.h"
#include "treewidth_search.h"

#include <sys/resource.h>

#include <chrono>
#include <iostream>

typedef std::chrono::high_resolution_clock::time_point TimePoint;

using namespace fpc;

void print_help() {
  std::cerr
      << "TLDC \"Too long; Didn't Count\" A length limited path counter.\n"
      << "Copyright (C) 2023-2024 Rafael Kiesel, Markus Hecher\n"
      << "License GPLv3+: GNU GPL version 3 or later "
         "<http://gnu.org/licenses/gpl.html>\n"
      << "\n"
      << "Usage: tldc [-s .] [-t .] [-m .] [-d] [-h]\n"
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
      << "                                          * branchfbs: "
         "frontier-based search along a branch decomposition\n"
      << "                                          * nautydfs: depth first "
         "search with modulo automorphisms\n"
      << "                                          * dfs: depth first search\n"
      << "                                          * mitm: meet in the middle "
         "search\n"
      << "       -t         --threads NUM           use NUM threads.\n"
      << "       -m         --memory_limit LIMIT    set a hardlimit of LIMIT "
         "MB memory.\n"
      << "       -d         --directed              expect a directed "
         "instance\n"
      << "       -h         --help                  print this help.\n";
}

template <typename count_t>
void print_result(std::vector<mpz_class> const &extra_paths,
                  std::vector<count_t> search_res, TimePoint start) {
  mpz_class final_result = 0;

  for (auto const &count : search_res) {
    final_result += to_mpz(count);
  }
  for (auto const &count : extra_paths) {
    final_result += count;
  }
  auto end = std::chrono::high_resolution_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
  std::cout << final_result << std::endl;
  LOG << "Time taken: " << duration << " seconds" << std::endl;
}

void print_result(mpz_class const &extra_paths, mpz_class const &search_res,
                  TimePoint start) {
  print_result<mpz_class>({extra_paths}, {search_res}, start);
}

enum SolverStrategy {
  AUTO,
  PATH_FBS,
  NAUTY_FBS,
  TREE_FBS,
  BRANCH_FBS,
  NAUTY_DFS,
  DFS,
  MITM
};

void print_stategy(SolverStrategy strategy) {
  switch (strategy) {
  case SolverStrategy::TREE_FBS:
    LOG << "tw ";
    break;
  case SolverStrategy::PATH_FBS:
    LOG << "pw ";
    break;
  case SolverStrategy::DFS:
    LOG << "dfs ";
    break;
  case SolverStrategy::MITM:
    LOG << "mitm ";
    break;
  case SolverStrategy::NAUTY_FBS:
    LOG << "pw(nauty) ";
    break;
  case SolverStrategy::BRANCH_FBS:
    LOG << "bw ";
    break;
  case SolverStrategy::NAUTY_DFS:
    LOG << "dfs(nauty) ";
    break;
  default:
    assert(false);
  }
}

void set_memory_limit(size_t mem_limit) {
  if (mem_limit == size_t(-1)) {
    return;
  }
  struct rlimit64 lim;
  getrlimit64(RLIMIT_AS, &lim);
  lim.rlim_cur = std::min(mem_limit, lim.rlim_cur);
  if (setrlimit64(RLIMIT_AS, &lim) == -1) {
    perror(strerror(errno));
    exit(-1);
  }
}

int solve_digraph(SolverStrategy strategy, size_t number_threads,
                  size_t mem_limit) {
  auto start_time = std::chrono::high_resolution_clock::now();
  AnnotatedDecomposition decomposition = AnnotatedDecomposition::invalid();
  Digraph initial_graph(std::cin);
  if (initial_graph.terminals().front() == initial_graph.terminals().back()) {
    LOG << "Terminals are equal, exiting with zero count" << std::endl;
    print_result(0, 0, start_time);
    return 0;
  }
  if (initial_graph.min_length() > initial_graph.max_length()) {
    LOG << "Minimum length less than max length, exiting with zero count"
        << std::endl;
    print_result(0, 0, start_time);
    return 0;
  }
  initial_graph.preprocess();
  initial_graph.normalize();
  initial_graph.print_stats();

  if (initial_graph.nr_arcs() == 0) {
    // solved by preprocessing
    LOG << "pp ";
    print_result(initial_graph.extra_paths(), 0, start_time);
    return 0;
  }

  auto weighted_graph = initial_graph.copy();
  weighted_graph.weighted_preprocess();
  weighted_graph.normalize();
  weighted_graph.print_stats();

  if (weighted_graph.nr_arcs() == 0) {
    // solved by weighted preprocessing
    LOG << "wpp ";
    print_result(weighted_graph.extra_paths(), 0, start_time);
    return 0;
  }

  if (strategy == SolverStrategy::AUTO) {
    auto pd_decomp =
        make_via_path_decomposition(number_threads, weighted_graph);
    LOG << "Path decomposition: ";
    pd_decomp.print_stats();
    auto td_decomp =
        make_via_tree_decomposition(number_threads, weighted_graph);
    LOG << "Tree decomposition: ";
    td_decomp.print_stats();
    auto join_adapted_tw =
        td_decomp.width() +
        (td_decomp.join_width() > td_decomp.width()
             ? (td_decomp.join_width() - td_decomp.width()) / 3
             : 0);
    if (join_adapted_tw < pd_decomp.width()) {
      decomposition = std::move(td_decomp);
      strategy = SolverStrategy::TREE_FBS;
    } else {
      decomposition = std::move(pd_decomp);
      strategy = SolverStrategy::PATH_FBS;
    }
    // DFS check
    std::vector<Vertex> all;
    all.reserve(initial_graph.nr_vertices());
    for (Vertex v = 0; v < initial_graph.nr_vertices(); v++) {
      all.push_back(v);
    }
    auto copy_graph = initial_graph.subgraph(all);
    copy_graph.reduce_modulo_equivalence();
    auto delta_length = copy_graph.max_length() - copy_graph.min_length();
    double avg_degree = copy_graph.nr_arcs() / double(copy_graph.nr_vertices());
    auto weighted_delta_length = delta_length * (avg_degree / 6.0);
    if (weighted_delta_length < std::min(join_adapted_tw, pd_decomp.width())) {
      strategy = SolverStrategy::DFS;
    }
    // MITM check
    MITMDirectedSearch<Unlimited_count, std::size_t> check_param(initial_graph,
                                                                 1);
    auto [perc, avg_path_size] = check_param.percent_done_after(1'000'000);
    double mitm_width = 37 * std::max(avg_path_size, 0.1) /
                        log((1'000'000 * std::max(perc, 0.000002)));
    LOG << "MITM width: " << mitm_width << std::endl;
    LOG << "DFS width: " << weighted_delta_length << std::endl;
    if (mitm_width < std::min({weighted_delta_length, double(join_adapted_tw),
                               double(pd_decomp.width())})) {
      strategy = SolverStrategy::MITM;
    }
  } else if (strategy == SolverStrategy::PATH_FBS) {

    decomposition = make_via_path_decomposition(number_threads, weighted_graph);
    LOG << "Path decomposition: ";
    decomposition.print_stats();
  } else if (strategy == SolverStrategy::TREE_FBS) {

    decomposition = make_via_tree_decomposition(number_threads, weighted_graph);
    LOG << "Tree decomposition: ";
    decomposition.print_stats();
  }

  print_stategy(strategy);
  set_memory_limit(mem_limit);

  switch (strategy) {
  case SolverStrategy::DFS: {
    auto solver = make_solver<ParallelDirectedSearch>(
        initial_graph, std::nullopt, number_threads);
    auto res = solver->search();
    solver->print_stats();
    print_result(initial_graph.extra_paths(), res, start_time);
  } break;
  case SolverStrategy::MITM: {
    auto solver =
        make_solver<MITMDirectedSearch>(initial_graph, number_threads);
    auto res = solver->search();
    solver->print_stats();
    print_result(initial_graph.extra_paths(), res, start_time);
  } break;
  case SolverStrategy::PATH_FBS:
  case SolverStrategy::TREE_FBS:
  case SolverStrategy::BRANCH_FBS: {
    auto solver = make_solver<TreewidthDirectedSearch>(
        weighted_graph, decomposition, number_threads);
    print_result(weighted_graph.extra_paths(), solver->search(), start_time);
    break;
  }
  default:
    assert(false);
  }
  return 0;
}

int solve_graph(SolverStrategy strategy, size_t number_threads,
                size_t mem_limit) {
  auto start_time = std::chrono::high_resolution_clock::now();
  AnnotatedDecomposition decomposition = AnnotatedDecomposition::invalid();
  Graph initial_graph(std::cin);
  if (!initial_graph.is_all_pair() &&
      initial_graph.terminals().front() == initial_graph.terminals().back()) {
    LOG << "Terminals are equal, exiting with zero count" << std::endl;
    print_result(0, 0, start_time);
    return 0;
  }
  if (initial_graph.min_length() > initial_graph.max_length()) {
    LOG << "Minimum length less than max length, exiting with zero count"
        << std::endl;
    print_result(0, 0, start_time);
    return 0;
  }
  initial_graph.preprocess();
  initial_graph.normalize();
  initial_graph.print_stats();

  if (initial_graph.nr_edges() == 0) {
    // solved by preprocessing
    LOG << "pp ";
    print_result(initial_graph.extra_paths(), 0, start_time);
    return 0;
  }

  auto weighted_graph = initial_graph.copy();
  weighted_graph.weighted_preprocess();
  weighted_graph.normalize();
  weighted_graph.print_stats();

  if (weighted_graph.nr_edges() == 0) {
    // solved by weighted preprocessing
    LOG << "wpp ";
    print_result(weighted_graph.extra_paths(), 0, start_time);
    return 0;
  }

  if (strategy == SolverStrategy::AUTO) {
    auto pd_decomp =
        make_via_path_decomposition(number_threads, weighted_graph);
    LOG << "Path decomposition: ";
    pd_decomp.print_stats();
    auto td_decomp =
        make_via_tree_decomposition(number_threads, weighted_graph);
    LOG << "Tree decomposition: ";
    td_decomp.print_stats();
    auto join_adapted_tw =
        td_decomp.width() +
        (td_decomp.join_width() > td_decomp.width()
             ? (td_decomp.join_width() - td_decomp.width()) / 3
             : 0);
    if (join_adapted_tw < pd_decomp.width()) {
      decomposition = std::move(td_decomp);
      strategy = SolverStrategy::TREE_FBS;
    } else {
      decomposition = std::move(pd_decomp);
      strategy = SolverStrategy::PATH_FBS;
    }
    // DFS check
    std::vector<Vertex> all;
    all.reserve(initial_graph.nr_vertices());
    for (Vertex v = 0; v < initial_graph.nr_vertices(); v++) {
      all.push_back(v);
    }
    auto copy_graph = initial_graph.subgraph(all);
    copy_graph.reduce_modulo_equivalence();
    auto delta_length = copy_graph.max_length() - copy_graph.min_length();
    // times 2, because every edge is like two arcs
    double avg_degree =
        2.0 * copy_graph.nr_edges() / double(copy_graph.nr_vertices());
    auto weighted_delta_length = delta_length * (avg_degree / 6.0);
    if (weighted_delta_length < std::min(join_adapted_tw, pd_decomp.width())) {
      strategy = SolverStrategy::DFS;
    }
    // MITM check
    MITMSearch<Unlimited_count, std::size_t> check_param(initial_graph, 1);
    auto [perc, avg_path_size] = check_param.percent_done_after(1'000'000);
    double mitm_width = 37 * std::max(avg_path_size, 0.1) /
                        log((1'000'000 * std::max(perc, 0.000002)));
    LOG << "MITM width: " << mitm_width << std::endl;
    LOG << "DFS width: " << weighted_delta_length << std::endl;
    if (mitm_width < std::min({weighted_delta_length, double(join_adapted_tw),
                               double(pd_decomp.width())})) {
      strategy = SolverStrategy::MITM;
    }
  } else if (strategy == SolverStrategy::NAUTY_FBS) {

    decomposition = make_via_path_decomposition(number_threads, initial_graph);
    LOG << "Path decomposition: ";
    decomposition.print_stats();

  } else if (strategy == SolverStrategy::PATH_FBS) {

    decomposition = make_via_path_decomposition(number_threads, weighted_graph);
    LOG << "Path decomposition: ";
    decomposition.print_stats();

  } else if (strategy == SolverStrategy::TREE_FBS) {
    decomposition = make_via_tree_decomposition(number_threads, weighted_graph);
    LOG << "Tree decomposition: ";
    decomposition.print_stats();
  } else if (strategy == SolverStrategy::BRANCH_FBS) {
    auto bd = BranchDecomposition::from_graph(
        weighted_graph, BranchDecomposition::Strategy::HICKS);
    decomposition =
        AnnotatedDecomposition::from_branch_decomposition(bd, weighted_graph);
    LOG << "Branch decomposition: ";
    decomposition.print_stats();
  }

  print_stategy(strategy);
  set_memory_limit(mem_limit);

  switch (strategy) {
  case SolverStrategy::PATH_FBS:
  case SolverStrategy::TREE_FBS:
  case SolverStrategy::BRANCH_FBS: {
    auto solver = make_solver<TreewidthSearch>(weighted_graph, decomposition,
                                               number_threads);
    print_result(weighted_graph.extra_paths(), solver->search(), start_time);
  } break;

  case SolverStrategy::NAUTY_FBS: {
    auto solver = make_solver<NautyPathwidthSearch>(
        initial_graph, decomposition, number_threads);
    print_result(initial_graph.extra_paths(), solver->search(), start_time);
  } break;

  case SolverStrategy::NAUTY_DFS: {
    std::vector<sparsegraph> sgs;
    if (initial_graph.is_all_pair()) {
      sgs = initial_graph.all_pair_nauty(true);
    } else {
      sgs.push_back(initial_graph.to_canon_nauty(true));
    }
    ParallelNautySearch search(initial_graph.to_canon_nauty(true),
                               initial_graph.max_length(), number_threads);
    search.add_to_initial(sgs);
    print_result({initial_graph.extra_paths()}, search.search(), start_time);
  } break;

  case SolverStrategy::DFS: {
    auto solver = make_solver<ParallelSearch>(initial_graph, std::nullopt,
                                              number_threads);
    print_result({initial_graph.extra_paths()}, solver->search(), start_time);
  } break;

  case SolverStrategy::MITM: {
    auto solver = make_solver<MITMSearch>(initial_graph, number_threads);
    auto res = solver->search();
    solver->print_stats();
    print_result(initial_graph.extra_paths(), res, start_time);
  } break;

  default:
    assert(false);
  }
  return 0;
}

int main(int argc, char **argv) {

  SolverStrategy strategy = SolverStrategy::AUTO;
  size_t number_threads = omp_get_max_threads();
  bool directed = false;
  size_t mem_limit = -1;

  size_t arg_idx = 0;
  while (++arg_idx < argc) {
    if (strcmp(argv[arg_idx], "-s") == 0 ||
        strcmp(argv[arg_idx], "--strategy") == 0) {
      assert(strategy == SolverStrategy::AUTO);
      arg_idx++;
      if (arg_idx == argc) {
        LOG << "Not enough arguments!" << std::endl;
        return -1;
      }
      if (strcmp(argv[arg_idx], "pathfbs") == 0) {
        strategy = SolverStrategy::PATH_FBS;
      } else if (strcmp(argv[arg_idx], "nautyfbs") == 0) {
        strategy = SolverStrategy::NAUTY_FBS;
      } else if (strcmp(argv[arg_idx], "treefbs") == 0) {
        strategy = SolverStrategy::TREE_FBS;
      } else if (strcmp(argv[arg_idx], "branchfbs") == 0) {
        strategy = SolverStrategy::BRANCH_FBS;
      } else if (strcmp(argv[arg_idx], "nautydfs") == 0) {
        strategy = SolverStrategy::NAUTY_DFS;
      } else if (strcmp(argv[arg_idx], "dfs") == 0) {
        strategy = SolverStrategy::DFS;
      } else if (strcmp(argv[arg_idx], "mitm") == 0) {
        strategy = SolverStrategy::MITM;
      } else if (strcmp(argv[arg_idx], "auto") == 0) {
        strategy = SolverStrategy::AUTO;
      } else {
        LOG << "Unknown solver strategy " << argv[arg_idx] << std::endl;
        return -1;
      }
    } else if (strcmp(argv[arg_idx], "-m") == 0 ||
               strcmp(argv[arg_idx], "--memory_limit") == 0) {
      arg_idx++;
      if (arg_idx == argc) {
        LOG << "Not enough arguments!" << std::endl;
        return -1;
      }
      mem_limit = atol(argv[arg_idx]) * 1024 * 1024;
    } else if (strcmp(argv[arg_idx], "-d") == 0 ||
               strcmp(argv[arg_idx], "--directed") == 0) {
      directed = true;
    } else if (strcmp(argv[arg_idx], "-t") == 0 ||
               strcmp(argv[arg_idx], "--threads") == 0) {
      arg_idx++;
      if (arg_idx == argc) {
        LOG << "Not enough arguments!" << std::endl;
        return -1;
      }
      number_threads = atol(argv[arg_idx]);
    } else if (strcmp(argv[arg_idx], "-h") == 0 ||
               strcmp(argv[arg_idx], "--help") == 0) {
      print_help();
      return 0;
    } else {
      LOG << "Unknown argument " << argv[arg_idx] << std::endl;
      print_help();
      return -1;
    }
  }
  if (directed) {
    return solve_digraph(strategy, number_threads, mem_limit);
  }
  return solve_graph(strategy, number_threads, mem_limit);
}
