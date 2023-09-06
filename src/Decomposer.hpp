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

#include <string>
#include "graph.h"
#include <vector>
#include <utility>
#include <limits>
#include "annotated_decomposition.hpp"

using namespace fpc;

#ifndef _DECOMPOSER_HPP_
#define _DECOMPOSER_HPP_

#define DEF_DECOMPOSER "./utils/htd" 
#define DEF_PATH_DECOMPOSER "./utils/path_decomp.py"
#define DEF_DECOMP_PARAMS "SEE popen2.cpp"

typedef std::tuple<int, std::vector<int>, std::map<int, std::vector<int>>, std::map<int, std::pair<std::vector<Edge>, std::vector<vertex_t>>>> Td_t;


class Decomposer {
	public:
		Decomposer(const char* const decomposer=DEF_DECOMPOSER, const char* const path_decomposer=DEF_PATH_DECOMPOSER, const char* const params=DEF_DECOMP_PARAMS);
		virtual ~Decomposer() {}

		//std::pair<Graph, std::map<int, std::pair<std::vector<Edge>, std::vector<vertex_t>>>>
		//std::pair<int, std::pair<std::map<int, std::vector<int>>, std::map<int, std::pair<std::vector<Edge>, std::vector<vertex_t>>>>>
		//std::vector<std::pair<Edge, std::vector<vertex_t>>> 

		AnnotatedDecomposition tree_decompose(/*const*/ Graph& graph,  bool path=false, size_t* width=nullptr, size_t* nr_bags = nullptr, size_t* max_join_child_bags=nullptr, size_t* max_join_bag=nullptr, size_t* nr_joins=nullptr);

		std::vector<std::pair<Edge, std::vector<vertex_t>>> path_decompose(/*const*/ Graph& graph);



		void stats(const Td_t&);
		void stats(const AnnotatedDecomposition&);

		Td_t	
		decompose(/*const*/ Graph& graph, bool path=false);
	protected:
		const char* decomposer, *path_decomposer, *params;


		void update_Join_bag(AnnotatedNode& c, std::vector<vertex_t> &b1, std::vector<vertex_t> &b2);

		std::pair<int, int> insertEdges(AnnotatedDecomposition& r, std::vector<vertex_t>& bag, std::vector<Edge>& edges, std::set<Edge>&used, size_t child, NodeType type, bool join);
};


#endif
