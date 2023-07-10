#include "Decomposer.hpp"
#include <cstdio>
#include "popen2.h"
#include <unistd.h>
#include <sstream>
#include <algorithm>
#include <unordered_map>

Decomposer::Decomposer(const char* const decomposer, const char* const path_decomposer, const char* params)
{
	this->decomposer = decomposer;
	this->path_decomposer = path_decomposer;
	//this->params = params;
}


#define BUF_SIZE 4096

std::pair<int, int> Decomposer::insertEdges(AnnotatedDecomposition& r, std::vector<vertex_t>& bag, std::vector<Edge>& edges,  std::set<Edge>&used, size_t child, NodeType type, bool join)
{
	int first = -1;
	int idx = -1;
	for(auto it = edges.begin(); (join && it == edges.end() && first < 0) || it != edges.end();) 
	{
		AnnotatedNode c;
		
		if (it != edges.end()) 
		{
			if (used.count(*it) > 0) {
				// std::cerr << "SKIP " << it->first << "," << it->second << std::endl;
				continue;
			}
			used.insert(*it);
			c.edge = (*it);
			c.type = type;
		}
		else
		{
			// std::cerr << "CREATE EMPTY JOIN " << child << std::endl;
			c.type = JOIN;	
			c.edge = std::make_pair((size_t)-1, (size_t)-1);
		}
		//std::cerr << it->first << "," << it->second << std::endl;
		c.parent = (size_t)-1;
		c.bag = bag;
        	idx = r.size();
		if (first < 0)	//link to the outside world
		{
			c.children = std::make_pair(child, (size_t)-1);
			if (child != (size_t)-1)
				r[child].parent = idx;
			
			first = idx;
			type = PATH_LIKE;
		}
		else	//linking from within
		{
			c.children = std::make_pair(idx - 1, (size_t)-1);
			r[idx-1].parent = idx;
		}
		r.push_back(std::move(c));

		if (it != edges.end())
			++it;
	}
	//first and last
	return std::make_pair(first, idx);
}

/*int Decomposer::createNode(AnnotatedDecomposition& r, std::vector<vertex_t>& bag, std::vector<Edge>::iterator& it, NodeType type)
{
	AnnotatedNode c;
	c.type = type;
	c.bag = bag;
	
	if (it != edges.end())
	{
		c.edge = *it;
		it++;
	}
	r.push_back(std::move(c));
	return r.size()-1;
}*/

void Decomposer::update_Join_bag(AnnotatedNode& c, std::vector<vertex_t> &b1, std::vector<vertex_t> &b2)
{
	std::set<vertex_t> c1 = std::set<vertex_t>(b1.begin(), b1.end());
	std::set<vertex_t> c2 = std::set<vertex_t>(b2.begin(), b2.end());
						
	//parent bag, compute a better one; more suited for join
	for (auto it = c.bag.begin(); it != c.bag.end(); ) {
		if (c1.count(*it) == 0 && c2.count(*it) == 0)	//new one, don't need in join node
			it = c.bag.erase(it);
		else
			++it;
	}
}

AnnotatedDecomposition Decomposer::tree_decompose(/*const*/ Graph& graph, bool path, size_t* width, size_t* nr_bags, size_t* max_join_child_bags, size_t* max_join_bag, size_t* nr_joins)
{
    auto td = std::move(decompose(graph,path));
    //stats(td);
    
    AnnotatedDecomposition r;
    auto actual_td = std::get<3>(td);
    auto succ = std::get<2>(td);
    //int cur = std::get<0>(td);

    std::unordered_map<int, int> td2r;
    std::vector<std::pair<int, int>> stack;
    std::set<Edge> used;
    //leaves
    //stack.insert(stack.end(), std::get<1>(td).begin(), std::get<1>(td).end());

    std::unordered_map<int, int> joins;

    if (width != nullptr)
    	*width = 0;
    if (max_join_child_bags != nullptr)
        *max_join_child_bags = 0;
    if (max_join_bag != nullptr)
        *max_join_bag = 0;
    if (nr_joins != nullptr)
        *nr_joins = 0;
    if (nr_bags != nullptr)
        *nr_bags = actual_td.size();

    for (auto it = succ.begin(); it != succ.end(); ++it)
    {
    	auto f = joins.find(it->second[0]);
	if (f == joins.end())
		joins.emplace(it->second[0],1);
	else
		f->second = f->second + 1;
    	//if (joins.count(it->first) > 0)
		//joins[it->second]++;
	/*else
    		joins[it->first]*/
    }


	//leaves
    for (auto it = std::get<1>(td).begin(); it != std::get<1>(td).end(); ++it)
    {
    	stack.push_back(std::make_pair(-1, *it));

    	/*auto jt = actual_td[*it].first.begin();
	if (jt != actual_td[*it].first.end())
    	int chld = createNode(r, actual_td[*it].second, jt, LEAF);
	auto idx = insertEdges(r, actual_td[*it].second, jt);
	r[chld].parent = idx.first;
	//r idx, td idx
	stack.push_back(std::make_pair(idx.second, *it));*/
    }

    //int cur = std::get<1>(td)[0];	//FIXME: just take any leave for now
    //std::cerr << "leaf " << cur << std::endl;
    while (stack.size() > 0) { //actual_td.count(cur) != 0) {
	auto cur = stack.front();
	stack.erase(stack.begin());
    	
	auto edges = actual_td[cur.second].first;
	auto bag = actual_td[cur.second].second;
	
	std::pair<int,int> idx = std::make_pair(-1, -1);

	// std::cerr << cur.second << " pred: " << cur.first << " join: " << joins[cur.second] << std::endl;

	if (width != nullptr)
		*width = std::max(*width, bag.size());

	if (td2r.count(cur.second) == 0)
		idx = insertEdges(r, bag, edges, used, (size_t)cur.first, cur.first < 0 ? LEAF : PATH_LIKE, joins[cur.second] > 1);
	else {
		// std::cerr << "already inserted!" << std::endl;
	}


	/*if (idx.first == -1)	//no edges, skip it then
	{
		assert(idx.second == -1);
		if (succ.count(cur.second) > 0) {	//only non-empty tds
			assert(td2r.count(succ[cur.second][0])==0);
			stack.push_back(std::make_pair(cur.first, succ[cur.second][0]));	
		} else {
			std::cerr << cur.second << " DIES OUT, pred: " << cur.first  << std::endl;
		}
		//otherwise this branch dies out. maybe a good thing?
	}
	else*/
	{
		if (idx.first != -1) {
			// map tdidx to r idx
			td2r[cur.second] = idx.first;
		}
		else {
			idx.second = (size_t)cur.first;
			// std::cerr << "no insert " << cur.second << "," << idx.first << std::endl;
			//cur.second = cur.first;
			//assert(cur.first > 0);
		}
		if (succ.count(cur.second) > 0)
		{
			// std::cerr << "count ok" << std::endl;
			//node already exists if idx.first == -1 AND it actually had edges
			auto par = idx.first == -1 && td2r.count(cur.second) > 0 ? cur.second : succ[cur.second][0];
			if (td2r.count(par) > 0)	//parent already exists, can only be a join node, right?
			{
				// std::cerr << "map ok" << std::endl;
				auto ridx = td2r[par];
				assert(r[ridx].children.first != (size_t)-1);
				{	//add intermediate join node
					assert(r[ridx].children.second == (size_t)-1);


					AnnotatedNode* join = nullptr;

					auto& cn = r[r[ridx].children.first];	//join child 1 bag
					if (r[ridx].type != JOIN)	//create edge-empty JOIN node
					{
						int pos = r.size();

						if (idx.first == -1)
							// map tdidx to r idx
							td2r[cur.second] = pos;
						//	else	//update successor
							//	td2r[par] = pos;

						// std::cerr << "ADD EMPTY JOIN " << pos  << std::endl;
						AnnotatedNode c;
						c.type = JOIN;
						c.edge = std::make_pair((size_t)-1, (size_t)-1);
						c.parent = ridx;

						/*std::set<vertex_t> c1 = std::set<vertex_t>(cn.bag.begin(), cn.bag.end());
						std::set<vertex_t> c2 = std::set<vertex_t>(bag.begin(), bag.end());*/
						
						//parent bag, compute a better one; more suited for join
						c.bag = r[ridx].bag; //bag;
						update_Join_bag(c, cn.bag, bag);
						/*for (auto it = c.bag.begin(); it != c.bag.end(); ) {
							if (c1.count(*it) == 0 && c2.count(*it) == 0)	//new one, don't need in join node
								it = c.bag.erase(it);
							else
								++it;
						}*/

						c.children = std::make_pair(r[ridx].children.first, idx.second);
						r[r[ridx].children.first].parent = pos;
						r[idx.second].parent = pos;	

						r.push_back(std::move(c));

						r[ridx].children = std::make_pair(pos, (size_t)-1);

						join = &r[r.size()-1];						
					} else {	//already (edge-empty) JOIN node created
						r[ridx].children.second = idx.second;
						r[idx.second].parent = ridx;
						update_Join_bag(r[ridx], cn.bag, bag);

						join = &r[ridx];
					}
					if (max_join_child_bags != nullptr)
						*max_join_child_bags = std::max(*max_join_child_bags, r[join->children.first].bag.size() + 
													r[join->children.second].bag.size()); 
					if (max_join_bag != nullptr)
						*max_join_bag = std::max(*max_join_bag, join->bag.size());
					if (nr_joins != nullptr)
						(*nr_joins)++;

				}
				/*else //first child or not a join node
				{
					r[ridx].children.first = idx.second;
				}*/
			}
			else {//if (idx.first != -1) {	//parent requires building 
				stack.push_back(std::make_pair(idx.second, par));	
			}/*} else
				stack.push_back(std::make_pair(cur.first, par));*/	

		}	//otherwise: root
		else {
			// std::cerr << cur.second << " DIES OUT, pred: " << cur.first  << std::endl;
		}
	}

	
	/*for(; it != edges.end()-1; ++it) 
	{
        	int idx = r.size();
		AnnotatedNode c;
		c.type = PATH_LIKE;
		c.parent = idx;
		c.edge = (*it);
		c.bag = bag;
            	r.push_back(std::move(c));
	}*/

	    	//std::cerr << cur << std::endl;
	/*if (std::get<2>(td).count(cur) == 0)	//no successor
		break;
	else
		//FIXME: extend to TDs (first element / one successor sufficient for PDs)
        	cur = std::get<2>(td)[cur][0];*/
    }
    //stats(r);
    return std::move(r);
}

std::vector<std::pair<Edge, std::vector<vertex_t>>> Decomposer::path_decompose(/*const*/ Graph& graph)
{
    auto td = std::move(decompose(graph,true));
    // stats(td);
    std::vector<std::pair<Edge, std::vector<vertex_t>>> r;
    auto actual_td = std::get<3>(td);
    //int cur = std::get<0>(td);

    int cur = std::get<1>(td)[0];	
    std::cerr << "leaf " << cur << std::endl;
    while (true) { //actual_td.count(cur) != 0) {
        auto edges = actual_td[cur].first;
        auto bag = actual_td[cur].second;
        for(auto edge : edges) {
            r.push_back(std::make_pair(edge, bag));
        }

    	//std::cerr << cur << std::endl;
	if (std::get<2>(td).count(cur) == 0)	//no successor
		break;
	else
        	cur = std::get<2>(td)[cur][0];
    }
    return std::move(r);
}


void AnnotatedNode::stats() const
{
	std::cerr << "TYPE " << type << ", parent: " << parent << ", edge: (" << edge.first << "," << edge.second << ") children: (" << children.first << "," << children.second << ") bag {";
	for (auto it = bag.begin(); it != bag.end(); ++it)
		std::cerr << *it << ",";
	std::cerr << "}" << std::endl;
}

void Decomposer::stats(const AnnotatedDecomposition& td)
{
	int pos = 0;
	std::cerr << "decomposed & annotated" << std::endl;
	for (auto it = td.begin(); it != td.end(); ++it)
	{
		std::cerr << "<" << (pos) << "> ";
		it->stats();
		++pos;
	}
}

void Decomposer::stats(const Td_t& td)
{
	std::cerr << "decomposed, root " << std::get<0>(td) << std::endl;

	for (auto it = std::get<1>(td).begin(); it != std::get<1>(td).end(); ++it)
		std::cerr << "leaf " << *it << std::endl;

	for (auto it = std::get<2>(td).begin(); it != std::get<2>(td).end(); ++it)
		for (auto jt = it->second.begin(); jt != it->second.end(); ++jt)
			std::cerr << "td edge " << it->first << " -> " << *jt << std::endl;

	for (auto it = std::get<3>(td).begin(); it != std::get<3>(td).end(); ++it)
	{	
		std::cerr << "node " << it->first << std::endl;
		std::cerr << "edges (" << it->second.first.size() << ") { ";
		for (auto jt = it->second.first.begin(); jt != it->second.first.end(); ++jt)
			std::cerr << jt->first << "," << jt->second << "; ";
		std::cerr << "}" << std::endl;
	
		std::cerr << "bag (" << it->second.second.size() << ") { ";
		for (auto jt = it->second.second.begin(); jt != it->second.second.end(); ++jt)
			std::cerr << *jt << "; ";
		std::cerr << "}" << std::endl;
	}
}

Td_t
Decomposer::decompose(/*const*/ Graph& graph, bool path)
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
	if (popen2(path ? this->path_decomposer : this->decomposer, NULL /*this->params*/, &in, &out) > 0)
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
			std::cerr << "TD with " << bags << " bags, width " << width << " and " << verts << " vertices" << std::endl;
			std::cout << width << " ";
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
			// std::cerr << "root " << root << std::endl;
			std::map<int, std::vector<int>> neighs;
			for (; b > 0 && fgets(buf, BUF_SIZE-1, fout) != NULL && sscanf(buf, "%d %d\n", &b1, &b2); --b)
			//for (; b > 0 && !feof(fout) && fscanf(fout, "%d %d\n", &b1, &b2) == 2; --b)
			{
				neighs[b1].push_back(b2);
				neighs[b2].push_back(b1);
				// std::cerr << b << "," << b1 << "," << b2 << std::endl;
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

