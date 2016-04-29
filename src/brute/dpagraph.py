import DPATrial
import graph


def dpa(n, m):
    outgraph = graph.make_complete_graph(m)
    dpat = DPATrial.DPATrial(m)
    i = m
    while i < n:
        outgraph[i] = dpat.run_trial(m)
        i += 1
    return outgraph




num_nodes = 100
avg_out_degree = 13

dpaGraph= dpa(num_nodes, avg_out_degree)

in_degrees = graph.in_degree_distribution(dpaGraph)

# def divide(x):
#     return x / float(num_nodes)
#
# normalized = map(divide, in_degrees.values())
#
#
#
# pyplot.loglog(in_degrees.keys(), normalized, 'bo')
# pyplot.xlabel('degree')
# pyplot.ylabel('normalized frequency')
# pyplot.title('DPA graph - normalized degree distribution')
# pyplot.show()