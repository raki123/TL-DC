#include "Decomposer.hpp"
#include <cstdio>
#include "popen2.h"
#include <unistd.h>
#include <sstream>
#include <algorithm>

Decomposer::Decomposer(const char* const decomposer, const char* params)
{
	this->decomposer = decomposer;
	this->params = params;
}


#define BUF_SIZE 1024


std::vector<std::pair<Edge, std::vector<vertex_t>>> Decomposer::decompose(/*const*/ Graph& graph, int child_nodes)
{
	std::vector<std::pair<Edge, std::vector<vertex_t>>> td;

	//FIXME: add this to graph.cpp
	Vertex nr_edges = 0;
    	for (Vertex v = 0; v < graph.adjacency_.size(); v++) 
	{
        	nr_edges += graph.neighbors(v).size();
		//for (auto i=graph.neighbors(v).begin(); i!= graph.neighbors(v).end(); ++i)
		//	std::cout << i << std::endl;
    	}
    	nr_edges /= 2;

	int in, out;
	if (popen2(this->decomposer, this->params, &in, &out) > 0)
	{
		std::stringstream s;
		s << "p tw " << graph.neighbors_.size() << " " << //171 << std::endl; //
		nr_edges << std::endl;
		//std::cout << s.str();
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
					//std::cout << ss;
					write(in, ss.c_str(), strlen(ss.c_str()));
				}
			}
		}
		//std::cout << i << std::endl;
		//fsync(in);
		close(in);
	
		FILE* fout = fdopen(out, "r");
			
		//std::cout << "done " << std::endl;

		char buf[BUF_SIZE] = {0};

		int bags, width, verts;

		//std::cout << "preread " << fout << std::endl;
		//if (read(out, buf, BUF_SIZE) > 0)
		//std::cout << fscanf(fout, "s td %d %d %d", &bags, &width, &verts) << std::endl;
		if (fscanf(fout, "s td %d %d %d\n", &bags, &width, &verts) == 3)
		{
			std::set<Edge> edges; //only cover an edge at most once
			//std::cout << "read " << std::endl;
			//sscanf(buf, "s td %d %d %d\n", &bags, &width, &verts);
			//std::cout << bags << "," << width << "," << verts << std::endl;
			for (; bags > 0 && fgets(buf, BUF_SIZE-1, fout) != NULL; --bags)
			{
				if (buf[0] != 'c') 
				{
					//std::cout << buf << std::endl;
					char* pos = buf + 2; //don't read the 'b'
					int v1;
					std::vector<vertex_t> bag;
					while ((sscanf(pos, "%d", &v1)) > 0) 
					{
						bag.push_back(v1-1);
						//std::cout << v1 << std::endl;
						while (*pos != '\0' && *pos != ' ')
							++pos;
						if (*pos == ' ')
							++pos;

					//todo: avoid copying (move?)

					}


					for (auto jt = bag.begin(); jt != bag.end(); ++jt)	
					{
						auto ngbs = std::move(graph.neighbors(*jt));
						std::set<Vertex> intersect;
						std::set_intersection(std::make_move_iterator(ngbs.begin()),
                    							std::make_move_iterator(ngbs.end()),
                    							bag.begin(), bag.end(),
									std::inserter(intersect, intersect.begin()));
						ngbs.swap(intersect);
						//set<Vertex> intersect;
						//set_intersection(bag.begin(), bag.end(), ngbs.begin(), ngbs.end(), std::inserter(intersect, intersect.begin()));
						for (auto it = ngbs.begin(); it != ngbs.end(); ++it) 
						{
							if (*jt < *it) 
							{
								Edge e(*jt, *it);
								if (edges.find(e) == edges.end())
								{
									edges.insert(e);
									td.push_back(std::pair<Edge, std::vector<vertex_t>> (std::move(e), bag));
								}
							}
						}
					}
				}
				else
				{
					bags++;
				}
			}
		}
		fclose(fout);
		//std::cout << "out" << std::endl;
		close(out);
	}
	return std::move(td);	
}
