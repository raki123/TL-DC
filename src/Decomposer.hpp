#include <string>
#include "graph.h"
#include <vector>
#include <utility>
#include <limits>

using namespace fpc;

#ifndef _DECOMPOSER_HPP_
#define _DECOMPOSER_HPP_

#define DEF_DECOMPOSER "../utils/htd"
#define DEF_DECOMP_PARAMS "SEE popen2.cpp"

class Decomposer {
	public:
		Decomposer(const char* const decomposer=DEF_DECOMPOSER, const char* const params=DEF_DECOMP_PARAMS);
		virtual ~Decomposer() {}

		std::pair<Graph, std::map<int, std::pair<std::vector<Edge>, std::vector<vertex_t>>>>
		//std::vector<std::pair<Edge, std::vector<vertex_t>>> 
		decompose(/*const*/ Graph& graph);
	protected:
		const char* decomposer, *params;
};


#endif
