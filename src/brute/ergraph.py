import random

from matplotlib import pyplot

import graph


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

percentage = 0.01
num_nodes = 27770

in_degrees = graph.in_degree_distribution(er(num_nodes, percentage))


def divide(x):
    return x / float(num_nodes)

normalized = map(divide, in_degrees.values())

pyplot.loglog(in_degrees.keys(), normalized, 'bo')
pyplot.xlabel('degree')
pyplot.ylabel('normalized frequency')
pyplot.title('ER graph - normalized degree distribution')
pyplot.show()