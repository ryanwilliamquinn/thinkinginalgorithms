"""
Utility functions for manipulating graphs
"""

import random
import copy
import DPATrial
import UPATrial
import time
from bfs import networkimport

EX_GRAPH0 = {0: {1, 2}, 1: set([]), 2: set([])}
EX_GRAPH1 = {0: {1, 4, 5}, 1: {2, 6}, 2: {3}, 3: {0}, 4: {1}, 5: {2}, 6: set([])}
EX_GRAPH2 = {0: {1, 4, 5}, 1: {2, 6}, 2: {3, 7}, 3: {7}, 4: {1}, 5: {2}, 6: set([]), 7: {3}, 8: {1, 2},
             9: {0, 3, 4, 5, 6, 7}}

# represent directed graphs with dictionary
# one entry for each node.  value is a dictionary of in and out nodes

# takes the number of nodes num_nodes and returns a dictionary corresponding to a complete directed graph with the specified number of nodes
# i.e. 4 = {0:{1,2,3}, 1: {0,2,3}, 2: {0,1,3}, 3: {0,1,2}}
# if num_nodes is not positive, return the empty graph {}
def make_complete_graph(num_nodes):
    """ makes a directed graph with all edges filled in """
    if num_nodes < 1:
        return {}
    nodelist = list(range(0, num_nodes))
    graph = dict.fromkeys(nodelist)
    for node in range(0, num_nodes):
        graph[node] = set([item for item in nodelist if item != node])
    return graph


# takes a directed graph digraph (represented as a dictionary) and computes the in-degrees for the nodes in the graph
# the function should return a dictionary with the same set of keys as digraph whose corresponding values are the number of edges whose head matches a particular node
# for graph 0, this should be: {0: 0, 1: 1, 2: 1}
def compute_in_degrees(digraph):
    """ computes in-degrees for given graph """
    degree_dict = dict.fromkeys(digraph.keys(), 0)
    for value in digraph.values():
        for val in value:
            degree_dict[val] += 1
    return degree_dict


# takes a directed graph digraph (represented as a dictionary) and computes the unnormalized ditribution of the in-degrees of the graph
# the function should return a dictionary whose keys correspond to in-degrees of nodes in the graph.
# the value associated with each particular in-degree is the number of nodes tht that in-degree.
# in-degrees with no corresponding nodes in the graph are not included in the dictionary
def in_degree_distribution(digraph):
    """ computes in-degrees distribution for given graph """
    in_degrees = compute_in_degrees(digraph)
    distribution = {}
    for value in in_degrees.values():
        if (value in distribution):
            distribution[value] += 1
        else:
            distribution[value] = 1
    return distribution


def add(x, y): return x + y


def sum_in_degrees(in_degree_graph):
    return reduce(add, in_degree_graph.values())


def avg_in_degree(in_degree_graph):
    numerator = 0;
    denominator = 0;
    for degree, frequency in in_degree_graph.iteritems():
        numerator += degree * frequency
        denominator += frequency
    return numerator / float(denominator)


def make_undirected_graph(num_nodes, probability):
    if num_nodes < 1:
        return {}
    nodelist = list(range(0, num_nodes))
    graph = dict.fromkeys(nodelist)
    for node in graph:
        graph[node] = set([])
    for node in range(0, num_nodes):
        for other_node in range(0, num_nodes):
            rand = random.random()
            if node != other_node and rand < probability:
                graph[node].add(other_node)
                graph[other_node].add(node)
    return graph


def er(n, p):
    outgraph = dict(enumerate(range(n)))
    for key in outgraph.keys():
        outgraph[key] = set([])
        nodes = set([item for item in outgraph.keys() if item != key])
        for node in nodes:
            rand = random.random()
            if rand < p:
                outgraph[key].add(node)
    return outgraph


def dpa(n, m):
    outgraph = make_complete_graph(m)
    dpat = DPATrial.DPATrial(m)
    i = m
    while i < n:
        outgraph[i] = dpat.run_trial(m)
        i += 1
    return outgraph


def upa(num_nodes, num_edges):
    outgraph = make_complete_graph(num_edges)
    upat = UPATrial.UPATrial(num_edges)
    node_index = num_edges
    while node_index < num_nodes:
        nodes_to_add = upat.run_trial(num_edges)
        outgraph[node_index] = nodes_to_add
        for node in nodes_to_add:
            if (outgraph.has_key(node)):
                outgraph[node].add(node_index)
            else:
                outgraph[node] = {node_index}

        node_index += 1
    return outgraph


def random_order(ugraph):
    """ takes a graph and returns a list of the nodes in the graph in some random order """
    cugraph = copy.copy(ugraph)
    node_list = list()
    while (len(cugraph) > 0):
        node = random.choice(cugraph.keys())
        node_list.append(node)
        del cugraph[node]
    return node_list

def fast_targeted_order(ugraph):
    """ Computes a node of the maximum degree in ugraph. If multiple nodes have the maximum degree, it chooses any of them (arbitrarily). Removes that node (and its incident edges) from ugraph. """
    graph = networkimport.copy_graph(ugraph)
    degree_sets = dict()
    removed_nodes = []
    max_in_degree = -1

    for node_index in graph.keys():
        degree = len(graph[node_index])

        if degree > max_in_degree:
            max_in_degree = degree
        if degree not in degree_sets:
            degree_sets[degree] = set()
        degree_sets[degree].add(node_index)

    for degree in range(max_in_degree):
        if degree not in degree_sets:
            degree_sets[degree] = set()

    for in_degree in reversed(range(max_in_degree + 1)):

        degree_set = degree_sets[in_degree]
        while len(degree_set) > 0:
            rand_node = degree_set.pop()
            for neighbor in graph[rand_node]:
                neighbor_degree = len(graph[neighbor])
                degree_sets[neighbor_degree].remove(neighbor)
                degree_sets[neighbor_degree - 1].add(neighbor)
                graph[neighbor].remove(rand_node)
            removed_nodes.append(rand_node)
            del graph[rand_node]

    return removed_nodes



# 1 degree_sets = dict()
# 1 removed_nodes = []
# 1 max_in_degree = -1
#
# n for node_index in graph.keys():
#  1   degree = len(graph[node_index])
#
#  1   if degree > max_in_degree:
#  1       max_in_degree = degree
#  1   if not degree_sets.has_key(degree):
#  1       degree_sets[degree] = set()
#  1   degree_sets[degree].add(node_index)
#
# n for degree in range(max_in_degree):
#   1  if not degree_sets.has_key(degree):
#   1      degree_sets[degree] = set()
#
# n for in_degree in reversed(range(max_in_degree)):
#  1   degree_set = degree_sets[in_degree]
#  n   while len(degree_set) > 0:
#   1      rand_node = degree_set.pop()
#   5      for neighbor in graph[rand_node]:
#    1         neighbor_degree = len(graph[neighbor])
#    1         degree_sets[neighbor_degree].remove(neighbor)
#    1         degree_sets[neighbor_degree - 1].add(neighbor)
#    1         graph[neighbor].remove(rand_node)
#   1      removed_nodes.append(rand_node)
#   1      del graph[rand_node]