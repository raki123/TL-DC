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


def find_2_sep(graph):
    prog_str = """
    {sep(X) : v(X)}2.
    {r(X)}:- v(X), not sep(X).
    :- e(X,Y), r(X), not r(Y), not sep(Y).
    :- e(Y,X), r(X), not r(Y), not sep(Y).
    ok_nr(X) :- v(X), not sep(X), not r(X).
    """
    for node in graph.nodes():
        prog_str += f"v({node}).\n"
    prog_str += ":- " + ','.join(f"not r({node})" for node in graph.nodes()) + ".\n"
    prog_str += ":- " + ','.join(f"not ok_nr({node})" for node in graph.nodes()) + ".\n"
    for edge in graph.edges():
        prog_str += f"e({edge[0]},{edge[1]}).\n"

def plot(graph):
    import matplotlib.pyplot as plt
    from networkx.drawing.nx_pydot import graphviz_layout
    labels = { node : str(node) for node in graph.nodes() }
    pos = graphviz_layout(graph, prog="dot")
    nx.draw_networkx(graph, pos)
    nx.draw_networkx_labels(graph, pos, labels)
    plt.axis("off")
    plt.show()

plot(graph)

# plot(graph)
# preprocess(graph, terminals)
# plot(graph)

def compute_fvs(graph):
    import subprocess
    enc = ""
    for v,w in graph.edges():
        enc += f"{v} {w}\n"
    q = subprocess.Popen(["/home/staff/rkiesel/projects/aspmc/aspmc/external/fvs/src/build/FeedbackVertexSet", "Appx"], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    output, err = q.communicate(input=enc.encode(), timeout = 30.0)

    res = [ int(v) for v in output.decode().split()[1:] ]
    return res

# print(bef, len(graph.nodes()), len(compute_fvs(graph)))
# from aspmc.graph.treedecomposition import from_graph

# td = from_graph(graph)
# print(td.width)

from aspmc.programs.program import Program
from aspmc.main import logger
logger.setLevel("DEBUG")

def std_enc(graph, terminals):
    prog_str = '''
    reach(X) :- start(X).
    reach(X) :- reach(Y), use(Y,X).
    :- not reach(X), goal(X). 
    '''
    prog_str += f"start({terminals[0]}).\n"
    prog_str += f"goal({terminals[1]}).\n"
    for n in graph.nodes():
        edges = graph.edges(n)
        for i, e in enumerate(edges):
            prog_str += f"use({e[0]},{e[1]}) :- reach({n})"
            for ep in edges:
                if e != ep:
                    prog_str += f", not use({ep[0]}, {ep[1]})"
            prog_str += ".\n"
    return prog_str

def unary_enc(graph, terminals, length):
    prog_str = "reach(X, 0) :- start(X).\n"
    prog_str += ":- goal(X)"
    for i in range(1,length+1):
        prog_str += f", not reach(X, {i})"
    prog_str += ".\n"
    prog_str += f"start({terminals[0]}).\n"
    prog_str += f"goal({terminals[1]}).\n"
    prog_str += f":- reach(X, L), reach(X, L'), L != L'.\n"
    shortest_from_start = nx.shortest_path(graph, source=terminals[0])
    shortest_to_goal = nx.shortest_path(graph, source=terminals[1])
    for n in graph.nodes():
        if n == terminals[1] \
            or (len(shortest_from_start[n]) + len(shortest_to_goal[n]) - 2) > length:
            continue
        edges = graph.edges(n)
        for i, e in enumerate(edges):
            if e[1] == terminals[0]:
                continue
            prog_str += f"reach({e[1]},L + 1) :- reach({e[0]}, L), L <= {length - len(shortest_to_goal[n]) + 1}, L>= {len(shortest_from_start[n]) - 1}"
            for ep in edges:
                if e != ep:
                    prog_str += f", not reach({ep[1]}, L+1)"
            prog_str += ".\n"
    return prog_str
print(unary_enc(graph,terminals,length))
# program = Program(program_str=unary_enc(graph, terminals, length))
# program.tpUnfold()
# program.choose_clark_completion()
# cnf = program._cnf
# cnf.evaluate()
