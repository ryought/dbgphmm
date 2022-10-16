from .dbgphmm import *
import networkx as nx

def hoge():
    print('hoge')

def to_networkx(dbg: PyDbg):
    G = nx.DiGraph()
    nodes = dbg.nodes()
    edges = dbg.edges()
    for (kmer, copy_num) in nodes:
        G.add_node(kmer, copy_num=copy_num)
    for (source, target, copy_num) in edges:
        G.add_edge(nodes[source][0], nodes[target][0], copy_num=copy_num)
    return G

def test():
    dbg = PyDbg(4, [
        StyledSequence("L:ATCGATTCGATTTAG"),
    ])
    print(dbg.k, dbg.n_nodes, dbg.n_edges)
    print(dbg.nodes())
    G = to_networkx(dbg)
    print(G.edges)
