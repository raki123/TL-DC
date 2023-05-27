import networkx as nx
import sys

name = sys.argv[1]
graph = nx.Graph()
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

from graphillion import GraphSet
GraphSet.set_universe(graph.edges())

# Count the paths
n = len(GraphSet.paths(*terminals))

# Output
print(n)  # 2 paths
