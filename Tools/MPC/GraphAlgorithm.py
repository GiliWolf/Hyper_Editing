# region Imports
# region GGPS Imports
import math

from Outputs.Outputers.CSVOutputer import CSVOutputer
# endregion
# region External Python Libraries
import gmpy2
import networkx as nx
import networkx.algorithms.community as community
import matplotlib.pyplot as plt
import pydot
from networkx.drawing.nx_pydot import graphviz_layout
import sys

# sys.setdefaultencoding() does not exist, here!
reload(sys)  # Reload does the trick!
sys.setdefaultencoding('UTF8')
import os

os.environ["PATH"] += os.pathsep + 'C:/Program Files (x86)/Graphviz2.38/bin/'


# endregion

# endregion
def intersection(lst1, lst2):
    temp = set(lst2)
    lst3 = [value for value in lst1 if value in temp]
    return lst3


def get_semi_identical_relations(num_of_bits, max_bit_diff):
    top_limit = 2 ** num_of_bits
    relations = dict()
    for i in xrange(top_limit):
        for j in xrange(i, top_limit):
            semi_identical = gmpy2.popcount(i ^ j) <= max_bit_diff
            relations.setdefault(i, dict())[j] = semi_identical
            relations.setdefault(j, dict())[i] = semi_identical
    return relations


def get_semi_identical_graph(num_of_bits, max_bit_diff):
    top_limit = 2 ** num_of_bits
    identity_graph = nx.Graph()
    identity_graph.add_nodes_from(xrange(top_limit))
    for i in xrange(top_limit):
        for j in xrange(i, top_limit):
            semi_identical = gmpy2.popcount(i ^ j) <= max_bit_diff
            if semi_identical:
                identity_graph.add_edge(i, j)
    return identity_graph


def common_identity(relations, node1, node2):
    commons = list()
    for other in relations[node1]:
        if relations[node1][other] and relations[node2][other]:
            commons.append(other)
    return commons


def create_graph_clique_minimisers(identity_graph, max_bit_diff):
    # type: (nx.Graph, int) -> None
    cliques = [c for c in nx.algorithms.clique.find_cliques(identity_graph) if len(c) > 2 ** max_bit_diff]
    identity_graph.add_node("RootOfAllEvil")
    for clique in cliques:
        clique_id = str(hash(tuple(clique)))
        identity_graph.add_node(clique_id, color="red")
        identity_graph.add_edge("RootOfAllEvil", clique_id)
        for node_id in clique:
            identity_graph.add_edge(node_id, clique_id)
            for other_node in clique:
                if other_node != node_id and identity_graph.has_edge(node_id, other_node):
                    identity_graph.remove_edge(node_id, other_node)


def get_graph_data(num_of_bits, max_bit_diff, output_path):
    identity_graph = get_semi_identical_graph(num_of_bits=num_of_bits, max_bit_diff=max_bit_diff)

    create_graph_clique_minimisers(identity_graph, max_bit_diff)

    with open(output_path, "wb") as o:
        # nx.write_gml(identity_graph, o)
        nx.write_edgelist(identity_graph, o)

    # pos = graphviz_layout(identity_graph, prog='twopi', args='')#nx.draw_shell(identity_graph, with_labels=True)
    # plt.figure(figsize=(20, 20))
    # nx.draw(identity_graph, pos, node_size=100, alpha=0.5, node_color="blue", with_labels=True)
    # plt.axis('equal')
    # plt.show()

    # nx.draw_circular(identity_graph, with_labels=True)
    # plt.show()

