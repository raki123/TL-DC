#include "Decomposer.hpp"
#include <cstdio>
#include "popen2.h"
#include <unistd.h>
#include <sstream>
#include <algorithm>

Decomposer::Decomposer(const char* const decomposer, const char* params)
{
	this->decomposer = decomposer;
	//this->params = params;
}


#define BUF_SIZE 1024

std::tuple<int, std::vector<int>, std::map<int, std::vector<int>>, std::map<int, std::pair<std::vector<Edge>, std::vector<vertex_t>>>>
//std::vector<std::pair<Edge, std::vector<vertex_t>>> 
Decomposer::decompose(/*const*/ Graph& graph)
{

	std::map<int, std::pair<std::vector<Edge>, std::vector<vertex_t>>> bgs; 
	//std::map<int, int> succ;
	std::map<int, std::vector<int>> succ;

	std::vector<int> leaves;

	//FIXME: add this to graph.cpp
	Vertex nr_edges = 0;
    	for (Vertex v = 0; v < graph.adjacency_.size(); v++) 
	{
        	nr_edges += graph.neighbors(v).size();
		//for (auto i=graph.neighbors(v).begin(); i!= graph.neighbors(v).end(); ++i)
		//	std::cerr << i << std::endl;
    	}
    	nr_edges /= 2;

	int in, out;
	int root = 0;
	if (popen2(this->decomposer, NULL /*this->params*/, &in, &out) > 0)
	{
		std::stringstream s;
		s << "p tw " << graph.neighbors_.size() << " " << //171 << std::endl; //
		nr_edges << std::endl;
		//std::cerr << s.str();
		write(in, s.str().c_str(), strlen(s.str().c_str()));
		int i = 0;
    		for (Vertex v = 0; v < graph.adjacency_.size(); v++) {
		//for (auto it = graph.neighbors_.begin(); it != graph.neighbors_.end(); ++it)
		//{
			//std::set<Vertex>& s = *it;
			auto ngbs = std::move(graph.neighbors(v));
			for (auto jt = ngbs.begin(); jt != ngbs.end(); ++jt) 
			{
				std::stringstream s;
				if (v < *jt) 
				{
					++i;
					s << (v+1) << " " << *(jt)+1 << std::endl;
					std::string ss = s.str();
					//std::cerr << ss;
					write(in, ss.c_str(), strlen(ss.c_str()));
				}
			}
		}
		//std::cerr << i << std::endl;
		//fsync(in);
		close(in);
	
		FILE* fout = fdopen(out, "r");
			
		//std::cerr << "done " << std::endl;

		char buf[BUF_SIZE] = {0};

		int bags, width, verts;


		//std::cerr << "preread " << fout << std::endl;
		//if (read(out, buf, BUF_SIZE) > 0)
		//std::cerr << fscanf(fout, "s td %d %d %d", &bags, &width, &verts) << std::endl;
		if (fscanf(fout, "s td %d %d %d\n", &bags, &width, &verts) == 3)
		{
			std::set<Edge> edges; //only cover an edge at most once
			//std::cerr << "read " << std::endl;
			//sscanf(buf, "s td %d %d %d\n", &bags, &width, &verts);
			std::cerr << bags << ",w: " << width << "," << verts << std::endl;
			int b = 1;
			for (; b <= bags && fgets(buf, BUF_SIZE-1, fout) != NULL; ++b)
			{
				if (buf[0] != 'c') 
				{
					assert(buf[0] == 'b');
					//std::cerr << buf << std::endl;
					char* pos = buf + 2; //don't read the 'b'
					int v1;

					int bid = -1;
					
					std::vector<vertex_t> bag;

					while ((sscanf(pos, "%d", &v1)) > 0) 
					{
						if (bid < 0)
							bid = v1;
						else 
						{
							bag.push_back(v1-1);
						}
						//std::cerr << v1 << std::endl;
							while (*pos != '\0' && *pos != ' ')
								++pos;
							if (*pos == ' ')
								++pos;
					}

					std::vector<Edge> td; //, std::vector<vertex_t>>> td;
					std::sort(bag.begin(), bag.end());
					
					for (auto jt = bag.begin(); jt != bag.end(); ++jt)	
					{
						std::set<Vertex> ngbs = std::move(
							graph.neighbors(*jt));
						
						//std::sort(ngbs.begin(), ngbs.end());

						std::set<Vertex> intersect;
						std::set_intersection(std::make_move_iterator(ngbs.begin()),
                    							std::make_move_iterator(ngbs.end()),
                    							bag.begin(), bag.end(),
									std::inserter(intersect, intersect.begin()));
						//ngbs.swap(intersect);
						//set<Vertex> intersect;
						//set_intersection(bag.begin(), bag.end(), ngbs.begin(), ngbs.end(), std::inserter(intersect, intersect.begin()));
						/*if ((*jt) == 18) {
							std::cerr << *jt << "," << bag.size() << "," << ngbs.size() << "," << intersect.size() << std::endl;
						for (auto it = ngbs.begin(); it != ngbs.end(); ++it) 
							std::cerr << "n " << *it << std::endl;
						for (auto it = bag.begin(); it != bag.end(); ++it) 
							std::cerr << "b " << *it << std::endl;}*/

						for (auto it = intersect.begin(); it != intersect.end(); ++it) 
						{
							/*if ((*jt) == 18)
								std::cerr << "inner " << *it << std::endl;*/
							if (*jt < *it) 
							{
								//vertex_t mi = std::min(*jt, *it), ma = std::max(*jt, *it);

								//Edge e(mi, ma);
								Edge e(*jt, *it);
								if (edges.find(e) == edges.end())
								{
									edges.insert(e);
									td.push_back(e); //, std::vector<vertex_t>> (std::move(e), bag));
								}
							}
						}
					}
					bgs.insert({bid, std::pair<std::vector<Edge>, std::vector<vertex_t>>(std::move(td), std::move(bag))});
				}
				else 
				{
					if (root == 0)
						sscanf(buf, "c r %d\n", &root);
					b--;
				}	
			}
			/*for (Vertex v = 0; v < graph.adjacency_.size(); v++) 
			{
				auto ngbs = std::move(graph.neighbors(v));
				for (auto jt = ngbs.begin(); jt != ngbs.end(); ++jt) 
				{
					if (v < *jt && edges.find(std::make_pair(v, (*jt))) == edges.end()) 
					{
						std::cerr << "MISSING EDGE " << v << "," << *jt << std::endl;
					}
				}
			}*/

			//std::cerr << "edges " << nr_edges << "; found edges " << edges.size() << std::endl;
			assert(edges.size() == nr_edges);
			int b1, b2;
			b = bags - 1;
			std::cerr << "root " << root << std::endl;
			std::map<int, std::vector<int>> neighs;
			for (; b > 0 && fgets(buf, BUF_SIZE-1, fout) != NULL && sscanf(buf, "%d %d\n", &b1, &b2); --b)
			//for (; b > 0 && !feof(fout) && fscanf(fout, "%d %d\n", &b1, &b2) == 2; --b)
			{
				neighs[b1].push_back(b2);
				neighs[b2].push_back(b1);
				std::cerr << b << "," << b1 << "," << b2 << std::endl;
			}
			for(root = 1; neighs[root].size() > 1; root++){}
			std::set<int> seen;
			std::vector<int> stack = {root};
			int cur = root;
			while (stack.size() > 0) {
			/*cur == root || neighs[cur].size() > 1) {
				if(seen.count(neighs[cur][0]) > 0) {
					succ[cur].push_back(neighs[cur][1]);
				} else {
					succ[cur].push_back(neighs[cur][0]);
				}*/
				cur = stack[stack.size() - 1];
				stack.pop_back();
				int sz = stack.size();
				for (auto it = neighs[cur].begin(); it != neighs[cur].end(); ++it)
					if (seen.count(*it) == 0) {
						//succ[cur].push_back(*it); //direction from the root
						succ[*it].push_back(cur); //direction towards the root
						stack.push_back(*it);
					}
				if (sz == stack.size())
					leaves.push_back(cur);
				seen.insert(cur);
				/*if (succ.find(cur) == succ.end())
					succ[cur]
				cur = succ[cur];*/
			}
			//std::cerr << seen.size() << "," << bags << std::endl;
			assert(seen.size() == bags);
		}
		fclose(fout);
		//std::cerr << "out" << std::endl;
		close(out);
	}
	//pair < root, pair<neighbors, map<int, {pair<edges, bag>}>>>
	/*return std::pair<int, std::pair<std::map<int, std::vector<int>>, std::map<int, std::pair<std::vector<Edge>, std::vector<vertex_t>>>>>(root, 
			std::pair<std::map<int, std::vector<int>>, std::map<int, std::pair<std::vector<Edge>, std::vector<vertex_t>>>>(succ, std::move(bgs)));	*/
	return std::make_tuple(root, std::move(leaves), std::move(succ), std::move(bgs));
}

