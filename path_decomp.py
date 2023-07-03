#!/home/rafael/miniconda3/bin/python
#!/home/hecher/minconda34/bin/python

import random
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


# print(length, terminals)
bef = len(graph.nodes())
def preprocess(graph, terminals):
    found = True
    while found:
        found = False
        nodes = list(graph.nodes())
        for node in nodes:
            if node in terminals:
                continue
            if len(graph.edges(node)) <= 2:
                if len(graph.edges(node)) <= 1:
                    graph.remove_node(node)
                # else:
                #     neigh = list(graph.neighbors(node))
                #     graph.remove_node(node)
                #     graph.add_edge(neigh[0], neigh[1])
                # found = True


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



bags, bsize = eliminate(graph)

print("s td {} {} {}".format(len(bags), bsize, len(graph.nodes())))

bn = 1
for b in bags:
	print("b {} {}".format(bn, " ".join(map(str,b))))
	bn += 1

for b in range(1,len(bags)+1):
	print("{} {}".format(b, b + 1))

