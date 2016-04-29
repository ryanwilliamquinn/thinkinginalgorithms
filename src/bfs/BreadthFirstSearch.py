"""
Do breadth-first search of a graph
"""

from collections import deque
import random

def bfs_visited(ugraph, start_node):
    """ Takes the undirected graph ugraph and the node start_node and returns the set consisting of all nodes that are visited by a breadth-first search that starts at start_node. """
    queue = deque()
    visited = set([start_node])
    queue.append(ugraph[start_node])
    while (len(queue) > 0):
        node = queue.pop()
        for neighbor in node:
            if (neighbor not in visited):
                visited.add(neighbor)
                queue.append(ugraph[neighbor])
    return visited

def cc_visited(ugraph):
    """ Takes the undirected graph ugraph and returns a list of sets, where each set consists of all the nodes (and nothing else) in a connected component, and there is exactly one set in the list for each connected component in ugraph and nothing else. """
    remaining_nodes = set(ugraph.keys())
    connected = list([])
    while(len(remaining_nodes) > 0):
        node = remaining_nodes.pop()
        connected_nodes = bfs_visited(ugraph, node)
        connected.append(connected_nodes)
        remaining_nodes = remaining_nodes.difference(connected_nodes)
    return connected

def largest_cc_size(ugraph):
    """ Takes the undirected graph ugraph and returns the size (an integer) of the largest connected component in ugraph. """
    connected_sets = cc_visited(ugraph)
    max_connected_length = 0
    for connected in connected_sets:
        connected_length = len(connected)
        if (connected_length > max_connected_length):
            max_connected_length = connected_length
    return max_connected_length

def compute_resilience(ugraph, attack_order):
    """ Takes the undirected graph ugraph, a list of nodes attack_order and iterates through the nodes in attack_order. For each node in the list, the function removes the given node and its edges from the graph and then computes the size of the largest connected component for the resulting graph. """
    longest_cc_sizes = list([largest_cc_size(ugraph)])
    for node_to_remove in attack_order:
        ugraph.pop(node_to_remove, None)
        for edges in ugraph.values():
            edges.discard(node_to_remove)
        longest_cc_sizes.append(largest_cc_size(ugraph))
    return longest_cc_sizes
