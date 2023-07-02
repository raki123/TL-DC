#!/home/hecher/minconda34/bin/python

import networkx as nx
import sys

#name = sys.argv[1]
graph = nx.Graph()
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


def separator(graph, sep, nbs):
	bags = []
	cmin = graph.nodes()
	sz = len(cmin)
	ngbs = cmin
	i = 10
	for c in nx.all_node_cuts(graph):
		ngb = set(c)
		for n in c:
			ngb.update(graph.neighbors(n))
		if sz > len(ngbs):
			sz = len(ngbs)
			cmin = c
			ngbs = ngb
		i += 1
		if i >= 10:
			break
	bags.append(sep + list(cmin))
	graph = nx.Graph(graph)
	graph.remove_nodes_from(cmin)
	if len(graph.nodes()) < 5:
		bags.append(list(graph.nodes()) + sep)
	else:
		for co in nx.connected_components(graph):
			bags.extend([list(s) + sep for s in separator(nx.induced_subgraph(graph, co), list(c))])	
	return bags

def eliminate(graph):
	bags = []
	pbag = []
	bsize = 0
	removed = [] 
	for i in sorted(graph.nodes()):
		bag = set(graph.neighbors(i))
		bag.update(pbag) #previous bag
		bag.add(i)
		bag.difference_update(removed)

		removed.append(i)

		if not bag.issubset(pbag): #otherwise skip this bag
			bags.append(bag)
			bsize = max(bsize, len(bag))

			pbag = set(bag)	#copy
			
		pbag.remove(i)


	return bags, bsize


b = separator(graph, [])
print(b, max([len(s) for s in b]))
exit(1)

bags, bsize = eliminate(graph)

print("s td {} {} {}".format(len(bags), bsize, len(graph.nodes())))

bn = 1
for b in bags:
	print("b {} {}".format(bn, " ".join(map(str,b))))
	bn += 1

for b in range(1,len(bags)+1):
	print("{} {}".format(b, b + 1))
