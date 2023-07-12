import networkx as nx
import sys

name = sys.argv[1]
graph = nx.Graph()
length = -1
terminals = []
with open(name, 'r') as in_file:
    for line in in_file:
        if line[0] == 'c':
            continue
        if line[0] == 'p':
            continue
        if line[0] == 'e':
            line = [ int(v) for v in line.split(' ')[1:] ]
            graph.add_edge(line[0], line[1])
        if line[0] == 'l':
            line = [ int(v) for v in line.split(' ')[1:] ]
            length = line[0]
        if line[0] == 't':
            line = [ int(v) for v in line.split(' ')[1:] ]
            terminals = line
sums = 0
for i in graph.nodes():
    for j in graph.nodes():
        if i < j:
            gen = nx.all_simple_paths(graph, i, j, length)
            sums += sum(1 for _ in gen)
print(sums)
