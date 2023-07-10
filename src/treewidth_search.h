#pragma once

#include "annotated_decomposition.hpp"
#include "graph.h"
#include "search.h"
#include <limits>
#include <omp.h>


namespace fpc {

typedef uint8_t frontier_index_t;
typedef std::vector<frontier_index_t> Frontier;
typedef std::pair<sparsegraph, Frontier> TWCacheKey;

struct twc_hash {
public:
    size_t operator()(TWCacheKey const& key) const {
      auto const& sg = key.first;
      return hasher__((const char *)sg.v, sizeof(edge_t)*(sg.nv + sg.nde) + sizeof(degree_t)*sg.nv);
    }  
};
struct twc_equal {
public:
    bool operator()(TWCacheKey const& lhs, TWCacheKey const& rhs) const {
        return lhs.first.nv == rhs.first.nv && lhs.first.nde == rhs.first.nde && std::memcmp(lhs.first.v, rhs.first.v, sizeof(edge_t)*(lhs.first.nv + lhs.first.nde) + sizeof(degree_t)*lhs.first.nv) == 0;
    }
};


typedef std::unordered_map<TWCacheKey, std::vector<Edge_weight>, twc_hash, twc_equal> Cache;

class TreewidthSearch {
  public:
    TreewidthSearch(Graph& input, AnnotatedDecomposition decomposition, size_t nthreads);

    std::vector<Edge_weight> search();
    void print_stats();

  const frontier_index_t invalid_index_ = std::numeric_limits<frontier_index_t>::max();
  const frontier_index_t no_edge_index_ = std::numeric_limits<frontier_index_t>::max() - 1;
  const frontier_index_t two_edge_index_ = std::numeric_limits<frontier_index_t>::max() - 2;
  private:
    size_t nthreads_;
    Graph graph_;
    Edge_length max_length_;
    bool is_all_pair_;
    std::vector<vertex_t> terminals_;
    AnnotatedDecomposition decomposition_;
    std::vector<std::vector<vertex_t>> remaining_edges_after_this_;

    std::vector<std::vector<frontier_index_t>> bag_local_idx_map_;
    std::vector<std::vector<vertex_t>> bag_local_vertex_map_;
    Edge_length invalid_distance_ = std::numeric_limits<Edge_length>::max();
    std::vector<std::vector<std::vector<Edge_length>>> bag_local_distance_;
    std::vector<sparsegraph> sparsegraph_after_this_;

    std::vector<Edge_weight> result_;
    std::vector<std::vector<Edge_weight>> thread_local_result_;

    std::vector<std::pair<Cache, Cache>> cache_;

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

    void propagateLoop(Frontier &frontier, size_t bag_idx, size_t last_idx, std::vector<Edge_weight>& partial_results, bool takeable, bool skippable, size_t thread_id);

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

    using block_iter = std::vector<std::pair<Frontier, std::vector<Edge_weight>>>::const_iterator;

    void restoreStep(
      Frontier &left, 
      std::vector<frontier_index_t>& cut_paths,
      std::vector<frontier_index_t>& paths,
      std::vector<std::pair<frontier_index_t, frontier_index_t>> restore,
      size_t cut_paths_size,
      size_t paths_size
    );

    void mergeStep(
      Frontier &left,
      size_t bag_idx,
      frontier_index_t idx,
      bool found_solution,
      std::vector<frontier_index_t>& cut_paths,
      std::vector<frontier_index_t>& paths,
      std::vector<Edge_weight>& left_result,
      block_iter begin,
      block_iter end
    );

    bool finalizeMerge(
      Frontier left,
      size_t bag_idx,
      bool found_solution,
      std::vector<frontier_index_t> const& cut_paths,
      std::vector<frontier_index_t> const& paths,
      std::vector<Edge_weight> const& left_result,
      std::vector<Edge_weight> const& right_result
    );
    bool merge(Frontier& left, Frontier const& right, size_t bag_idx, std::vector<Edge_weight>& left_result, std::vector<Edge_weight> const& right_result);

    void advance(Frontier& frontier, size_t bag_idx);

    sparsegraph construct_sparsegraph(Frontier const& frontier, size_t last_idx);

        // stats
    std::vector<size_t> pos_hits_;
    std::vector<size_t> neg_hits_;

    // size_t splits = 0;

    std::vector<size_t> edges_;
    std::vector<size_t> propagations_;
    std::vector<size_t> merges_;
    std::vector<size_t> unsuccessful_merges_;

};

} // namespace