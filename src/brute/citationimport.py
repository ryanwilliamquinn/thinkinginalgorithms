"""
Provided code for Application portion of Module 1

Imports physics citation graph
"""

# general imports
import urllib2


# Set timeout for CodeSkulptor if necessary
#import codeskulptor
#codeskulptor.set_timeout(20)


###################################
# Code for loading citation graph
import graph

CITATION_URL = "http://storage.googleapis.com/codeskulptor-alg/alg_phys-cite.txt"

def load_graph(graph_url):
    """
    Function that loads a graph given the URL
    for a text representation of the graph

    Returns a dictionary that models a graph
    """
    graph_file = urllib2.urlopen(graph_url)
    graph_text = graph_file.read()
    graph_lines = graph_text.split('\n')
    graph_lines = graph_lines[ : -1]

    print "Loaded graph with", len(graph_lines), "nodes"

    answer_graph = {}
    for line in graph_lines:
        neighbors = line.split(' ')
        node = int(neighbors[0])
        answer_graph[node] = set([])
        for neighbor in neighbors[1 : -1]:
            answer_graph[node].add(int(neighbor))

    return answer_graph

citation_graph = load_graph(CITATION_URL)

total = float(27770)

def divide(x):
    """ division method """
    return x / total;

in_degrees = graph.in_degree_distribution(citation_graph)

print(in_degrees)
print graph.avg_in_degree(in_degrees)

# normalized = map(divide, in_degrees.values())
# print(len(normalized), "length")
#
# def add(x,y): return x + y
#
# print(reduce(add, normalized))
#
# print normalized
#
# pyplot.loglog(in_degrees.keys(), normalized, 'bo')
# pyplot.xlabel('degree')
# pyplot.ylabel('normalized frequency')
# pyplot.title('citation graph - normalized degree distribution')
# pyplot.show()



# how to determine the average out degree of this graph


