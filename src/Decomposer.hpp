#include <string>
#include "graph.h"
#include <vector>
#include <utility>
#include <limits>

using namespace fpc;

#ifndef _DECOMPOSER_HPP_
#define _DECOMPOSER_HPP_

#define DEF_DECOMPOSER "../utils/htd"
#define DEF_DECOMP_PARAMS "--child-limit 1"

class Decomposer {
	public:
		Decomposer(const char* const decomposer=DEF_DECOMPOSER, const char* const params=DEF_DECOMP_PARAMS);
		virtual ~Decomposer() {}

		std::vector<std::pair<Edge, std::vector<vertex_t>>> decompose(/*const*/ Graph& graph, int child_nodes=1);
	protected:
		const char* decomposer, *params;
};


#endif
