#pragma once

#include "graph.h"
#include "search.h"
#include <limits>
#include <omp.h>


namespace fpc {

typedef uint8_t frontier_index_t;
typedef std::vector<frontier_index_t> Frontier;

struct vec_hash {
public:
    size_t operator()(Frontier const& key) const {
        return hasher__(key);
    }  
};

class TreewidthSearch {
  public:
    TreewidthSearch(Graph& input, std::vector<std::pair<Edge, std::vector<vertex_t>>> path_decomposition, size_t nthreads);

    std::vector<Edge_weight> search();
    void print_stats();
  private:
    size_t nthreads_;
    Graph graph_;
    Edge_length max_length_;
    bool is_all_pair_;
    std::vector<vertex_t> terminals_;
    std::vector<std::pair<Edge, std::vector<vertex_t>>> path_decomposition_;
    std::vector<std::vector<vertex_t>> remaining_edges_after_this_;

    frontier_index_t invalid_index_ = std::numeric_limits<frontier_index_t>::max();
    frontier_index_t no_edge_index_ = std::numeric_limits<frontier_index_t>::max() - 1;
    frontier_index_t two_edge_index_ = std::numeric_limits<frontier_index_t>::max() - 2;
    std::vector<std::vector<frontier_index_t>> bag_local_idx_map_;
    std::vector<std::vector<vertex_t>> bag_local_vertex_map_;
    Edge_length invalid_distance_ = std::numeric_limits<Edge_length>::max();
    std::vector<std::vector<std::vector<Edge_length>>> bag_local_distance_;

    std::vector<Edge_weight> result_;
    std::vector<std::vector<Edge_weight>> thread_local_result_;

    std::vector<std::unordered_map<Frontier, std::vector<Edge_weight>, vec_hash>> cache_;

    std::vector<std::vector<char>> takeable_ = {
      {true, false, true, true, true, false, true, true},
      {false, false, false, false, false, false, false, false},
      {true, false, false, false, true, false, true, false},
      {true, false, false, false, true, false, true, false},
      {true, false, true, true, true, false, true, true},
      {false, false, false, false, false, false, false, false},
      {true, false, true, true, true, false, true, true},
      {true, false, false, false, true, false, true, false}
    };
    std::vector<std::vector<char>> skippable_ = {
      {true, true, true, false, true, true, true, true},
      {true, true, true, false, true, true, true, true},
      {true, true, true, false, true, true, true, true},
      {false, false, false, false, false, false, false, false},
      {true, true, true, false, true, true, true, true},
      {true, true, true, false, true, true, true, true},
      {true, true, true, false, true, true, true, true},
      {true, true, true, false, true, true, true, true}
    };

    void includeSolutions(Frontier const& frontier, size_t bag_idx, std::vector<Edge_weight> const& partial_result);

    bool canTake(Frontier& frontier, size_t bag_idx, std::vector<Edge_weight> const& partial_result);
    bool canSkip(Frontier& frontier, size_t bag_idx, std::vector<Edge_weight> const& partial_result);

    bool distancePrune(
      Frontier& frontier,
      std::vector<frontier_index_t> const& paths,
      std::vector<frontier_index_t> const& cut_paths,
      size_t bag_idx,
      size_t offset);

    void take(Frontier& frontier, size_t bag_idx);
    void skip(Frontier& frontier, size_t bag_idx);

    void advance(Frontier& frontier, size_t bag_idx);

        // stats
    std::vector<size_t> pos_hits_;
    std::vector<size_t> neg_hits_;

    // size_t splits = 0;

    std::vector<size_t> edges_;
    std::vector<size_t> propagations_;

};

} // namespace