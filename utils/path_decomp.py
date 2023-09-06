#!/usr/bin/env python3

#import random
#import networkx as nx
import sys

class Graph:
	def __init__(self):
		self._adj = {}
		self._nodes = set()

	def edges(self):
		return self._adj

	def nodes(self):
		return self._nodes

	def _add(self, x, y):
		if x not in self._adj:
			self._adj[x] = set()
		self._adj[x].add(y)

	def add_edge(self, x, y):
		self._add(x,y)
		self._add(y,x)
		self._nodes.add(x)
		self._nodes.add(y)

	def neighbors(self, x):
		return self._adj[x]

#x=Graph()
#x.add_edge(2,3)
#print(x.neighbors(2))
#print(x.neighbors(3))
#
#x.add_edge(4,2)
#
#print(x.neighbors(2))
#print(x.neighbors(4))
#print(x.neighbors(3))
#
#print(x.nodes())
#
#name = sys.argv[1]
#graph = nx.Graph()
graph = Graph()
length = -1
terminals = []
#with open(name, 'r') as in_file:
with sys.stdin as in_file:
    for line in in_file:
        if line[0] == 'c':
            continue
        elif line[0] == 'p':
            continue
        elif line[0] == 'l':
            line = [ int(v) for v in line.split(' ')[1:] ]
            length = line[0]
        elif line[0] == 't':
            line = [ int(v) for v in line.split(' ')[1:] ]
            terminals = line
        else: # line[0] == 'e'
            if line[0] == 'e':
                line = line[2:]
            line = [ int(v) for v in line.split(' ')[0:] ]
            graph.add_edge(line[0], line[1])

#print(graph.nodes())
#print(graph.edges())

## print(length, terminals)
#bef = len(graph.nodes())
#def preprocess(graph, terminals):
#    found = True
#    while found:
#        found = False
#        nodes = list(graph.nodes())
#        for node in nodes:
#            if node in terminals:
#                continue
#            if len(graph.edges(node)) <= 2:
#                if len(graph.edges(node)) <= 1:
#                    graph.remove_node(node)
#                # else:
#                #     neigh = list(graph.neighbors(node))
#                #     graph.remove_node(node)
#                #     graph.add_edge(neigh[0], neigh[1])
#                # found = True
#
#
def eliminate(graph):
	bags = []
	pbag = []
	bsize = 0
	removed = []

	ngbs_done  = {}

	for i in graph.nodes():
		ngbs_done[i] = set(graph.neighbors(i))


	order = sorted(graph.nodes())

#	for i in range(2):
#		x=random.randint(0, len(order)-1)
#		y=random.randint(0, len(order)-1)
#
#		if x != y:
#			t = order[x]
#			order[x] = order[y]
#			order[y] = t
#
	for i in order:
		if i not in ngbs_done:
			continue

		bag = set(graph.neighbors(i))
		bag.update(pbag) #previous bag

		bag.add(i)
		bag.difference_update(removed)

		removed.append(i)
		del ngbs_done[i]


		if not bag.issubset(pbag): #otherwise skip this bag
			bags.append(bag)
			bsize = max(bsize, len(bag))

			pbag = set(bag)	#copy
			
		pbag.remove(i)
	
		# update open edges and delete if possible
		for j in bag:
			if j != i:
				ngbs_done[j].difference_update(bag)
				if len(ngbs_done[j]) == 0:
					removed.append(j)
					pbag.remove(j)
					del ngbs_done[j]



	return bags, bsize

def other_eliminate(graph):
	bags = []
	bsize = 0
	active = set()
	degree = { i : len(list(graph.neighbors(i))) for i in graph.nodes() }
	order = sorted(graph.nodes())
	for v in order:
		rem_edges_from = active.intersection(graph.neighbors(v))
		for rem in rem_edges_from:
			degree[rem] -= 1
		degree[v] -= len(rem_edges_from)
		active.add(v)
		bags.append(active)
		bsize = max(bsize, len(active))
		active = { x for x in active if degree[x] > 0 }

	return bags, bsize

bags, bsize = other_eliminate(graph)
#bags, bsize = eliminate(graph)

print("s td {} {} {}".format(len(bags), bsize, len(graph.nodes())))

bn = 1
for b in bags:
	print("b {} {}".format(bn, " ".join(map(str,b))))
	bn += 1

for b in range(1,len(bags)+1):
	print("{} {}".format(b, b + 1))

