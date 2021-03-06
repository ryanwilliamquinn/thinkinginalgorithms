import graph
import networkimport
import BreadthFirstSearch
from matplotlib import pyplot
import time

undirected_er_graph = graph.make_undirected_graph(1239, .002)   # make the undirected graph with 1239 nodes, ~3049 edges
network_graph = networkimport.load_graph(networkimport.NETWORK_URL) # make the undirected network graph with 1239 nodes, ~3100 edges
upa_graph = graph.upa(1239, 3)    # make the undirected upa graph with 1239 nodes, ~3700 edges

undirected_er_graph_attack_order = graph.random_order(undirected_er_graph)
undirected_er_graph_resilience = BreadthFirstSearch.compute_resilience(undirected_er_graph, undirected_er_graph_attack_order)

network_graph_attack_order = graph.random_order(network_graph)
network_graph_resilience = BreadthFirstSearch.compute_resilience(network_graph, network_graph_attack_order)

upa_graph_attack_order = graph.random_order(upa_graph)
upa_graph_resilience = BreadthFirstSearch.compute_resilience(upa_graph, upa_graph_attack_order)

print(undirected_er_graph_resilience)
print(network_graph_resilience)
print(upa_graph_resilience)

num_nodes_removed = range(1240)
print len(num_nodes_removed)
print len(undirected_er_graph_resilience)
print len(upa_graph_resilience)
print len(network_graph_resilience)

pyplot.plot(num_nodes_removed, undirected_er_graph_resilience, '-b', label='ER Graph, p: .0019')
pyplot.plot(num_nodes_removed, upa_graph_resilience, '-r', label='UPA Graph, m: 3')
pyplot.plot(num_nodes_removed, network_graph_resilience, '-g', label='Network Graph')
pyplot.legend(loc='upper right')
pyplot.ylabel('size of the largest connected component')
pyplot.xlabel('number of nodes removed')
pyplot.title('comparison of graph resiliency')

# Examine the shape of the three curves from your plot in Question 1.
# Which of the three graphs are resilient under random attacks as the first 20% of their nodes are removed?
# Note that there is no need to compare the three curves against each other in your answer to this question.

print("order is ER, UPA, network")
for giraffe in [undirected_er_graph_resilience, upa_graph_resilience, network_graph_resilience]:
    print("at index 248, must be above 744 to be resilient:", giraffe[248])
pyplot.show()

# node_count_list = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, 370, 380, 390, 400, 410, 420, 430, 440, 450, 460, 470, 480, 490, 500, 510, 520, 530, 540, 550, 560, 570, 580, 590, 600, 610, 620, 630, 640, 650, 660, 670, 680, 690, 700, 710, 720, 730, 740, 750, 760, 770, 780, 790, 800, 810, 820, 830, 840, 850, 860, 870, 880, 890, 900, 910, 920, 930, 940, 950, 960, 970, 980, 990]
# fast_target_times = [6.40000000000085e-05, 0.00011099999999999999, 0.00016000000000000736, 0.00020300000000000873, 0.0002550000000000052, 0.00031199999999999284, 0.00036599999999999133, 0.0004200000000000037, 0.0004640000000000061, 0.0005240000000000106, 0.0005690000000000001, 0.0006159999999999916, 0.0006509999999999988, 0.0007440000000000085, 0.000793000000000002, 0.0008390000000000064, 0.0009440000000000004, 0.0009450000000000014, 0.001022000000000009, 0.0012849999999999945, 0.0011040000000000078, 0.0011759999999999965, 0.0011929999999999996, 0.0012179999999999969, 0.0012550000000000061, 0.0014139999999999986, 0.001376000000000016, 0.0013960000000000083, 0.001487000000000016, 0.0015550000000000008, 0.0016329999999999956, 0.001696000000000003, 0.0017159999999999953, 0.0018139999999999823, 0.0018479999999999885, 0.0018400000000000083, 0.0019829999999999848, 0.002017999999999992, 0.0020449999999999913, 0.002029000000000003, 0.0020840000000000025, 0.0022039999999999837, 0.0021599999999999953, 0.0022370000000000168, 0.0026810000000000445, 0.002350999999999992, 0.002385000000000026, 0.002581999999999973, 0.0024399999999999977, 0.002572999999999992, 0.002681999999999962, 0.0026710000000000345, 0.00270999999999999, 0.0027870000000000394, 0.0029000000000000137, 0.0030140000000000167, 0.0028229999999999644, 0.003593000000000013, 0.0030180000000000207, 0.0030789999999999984, 0.0032300000000000106, 0.003292000000000017, 0.003279000000000032, 0.003295000000000048, 0.0040060000000000096, 0.003383999999999998, 0.004074000000000022, 0.004288999999999987, 0.004173000000000038, 0.004315000000000013, 0.004446999999999979, 0.004322999999999966, 0.0045319999999999805, 0.004553999999999947, 0.0045840000000000325, 0.004502000000000006, 0.004841000000000095, 0.004713000000000078, 0.004989000000000021, 0.004369000000000067, 0.004804999999999948, 0.004434000000000049, 0.004497000000000084, 0.00468500000000005, 0.005259999999999931, 0.0048520000000000785, 0.005083999999999977, 0.00504899999999997, 0.0049599999999999644, 0.005029999999999979, 0.005057000000000089, 0.005224000000000006, 0.005557999999999952, 0.005265999999999993, 0.005299000000000054, 0.005519999999999969, 0.00541999999999998, 0.00558599999999998, 0.005877000000000021]
# target_times = [3.999999999999837e-05, 8.600000000000274e-05, 0.00015099999999999836, 0.00021200000000000385, 0.0002979999999999927, 0.0004009999999999986, 0.0005090000000000094, 0.0006339999999999957, 0.000792000000000001, 0.0009730000000000016, 0.001117999999999994, 0.001362000000000002, 0.0014850000000000002, 0.001689999999999997, 0.001902000000000001, 0.002120999999999998, 0.0024160000000000015, 0.002607999999999999, 0.0029570000000000013, 0.003253999999999993, 0.0034460000000000046, 0.003719, 0.004034999999999983, 0.0045439999999999925, 0.00464500000000001, 0.005120999999999987, 0.005535000000000012, 0.005841000000000013, 0.006336000000000008, 0.006733000000000017, 0.007197999999999982, 0.007440999999999975, 0.007697999999999983, 0.008383000000000002, 0.009309999999999985, 0.009628999999999999, 0.00989699999999999, 0.01079399999999997, 0.011093999999999993, 0.011807000000000012, 0.01660600000000001, 0.013554999999999984, 0.013141999999999987, 0.014130000000000031, 0.014597000000000027, 0.015049999999999952, 0.01605200000000001, 0.016687000000000007, 0.01728200000000002, 0.018513, 0.018727999999999967, 0.01990000000000003, 0.020485000000000086, 0.021255000000000024, 0.021509, 0.02262299999999995, 0.022764999999999924, 0.023615000000000053, 0.024442999999999993, 0.02527900000000005, 0.025712999999999986, 0.02635899999999991, 0.029897000000000062, 0.02869200000000005, 0.028939999999999855, 0.029554999999999998, 0.030442999999999998, 0.03572500000000001, 0.0320600000000002, 0.033552000000000026, 0.03465700000000016, 0.035230000000000095, 0.037016999999999856, 0.03751899999999986, 0.03830499999999981, 0.03839999999999999, 0.04049499999999995, 0.041828999999999894, 0.04231499999999988, 0.042327000000000004, 0.044459000000000026, 0.04469699999999999, 0.046466999999999814, 0.0468059999999999, 0.04801000000000011, 0.049745999999999846, 0.050120000000000164, 0.05202899999999966, 0.05197699999999994, 0.05349899999999996, 0.05518500000000026, 0.05459499999999995, 0.057141999999999804, 0.05859900000000007, 0.05919499999999989, 0.06054400000000015, 0.06207000000000029, 0.06294999999999984, 0.06374700000000022]

# def timed_fast_target(range):
#     upa = graph.upa(range, 5)
#     start = time.time()
#     graph.fast_targeted_order(upa)
#     return time.time() - start
#
# def target_times(range):
#     upa = graph.upa(range, 5)
#     start = time.time()
#     networkimport.targeted_order(upa)
#     return time.time() - start
#
#
# node_count_list = range(10, 1000, 10)
# print node_count_list
# fast_target_times = map(timed_fast_target, node_count_list)
# target_times = map(target_times, node_count_list)
#
# pyplot.plot(node_count_list, fast_target_times, '-b', label='fast targeted order')
# pyplot.plot(node_count_list, target_times, '-r', label='targeted order')
# pyplot.legend(loc='upper right')
# pyplot.ylabel('running time')
# pyplot.xlabel('number of nodes in graph')
# pyplot.title('targeted removal order run time (desktop Python)')
# pyplot.show()

#
# undirected_er_graph = graph.make_undirected_graph(1239, .002)   # make the undirected graph with 1239 nodes, ~3100 edges
# network_graph= networkimport.load_graph(networkimport.NETWORK_URL) # make the undirected network graph with 1239 nodes, ~3100 edges
# upa_graph = graph.upa(1239, 3)    # make the undirected upa graph with 1239 nodes, ~3700 edges
#
#
# undirected_er_graph_attack_order = graph.fast_targeted_order(undirected_er_graph)
# targeted_er_graph_resilience = BreadthFirstSearch.compute_resilience(undirected_er_graph, undirected_er_graph_attack_order)
#
# network_graph_attack_order = graph.fast_targeted_order(network_graph)
# targeted_network_graph_resilience = BreadthFirstSearch.compute_resilience(network_graph, network_graph_attack_order)
#
# upa_graph_attack_order = graph.fast_targeted_order(upa_graph)
# targeted_upa_graph_resilience = BreadthFirstSearch.compute_resilience(upa_graph, upa_graph_attack_order)
#
# for giraffe in [targeted_er_graph_resilience, targeted_upa_graph_resilience, targeted_network_graph_resilience]:
#     print("at index 248, must be above 744 to be resilient:", giraffe[248])
# # pyplot.show()
# #
# num_nodes_removed = range(1240)
#
# # for i in num_nodes_removed:
# #     print("er graph: ", targeted_er_graph_resilience[i], i)
# #     print("upa graph: ", targeted_upa_graph_resilience[i])
# #     print("network graph: ", targeted_network_graph_resilience[i])
# #     i += 100
#
# pyplot.plot(num_nodes_removed, targeted_er_graph_resilience, '-b', label='ER Graph, p: .0019')
# pyplot.plot(num_nodes_removed, targeted_upa_graph_resilience, '-r', label='UPA Graph, m: 3')
# pyplot.plot(num_nodes_removed, targeted_network_graph_resilience, '-g', label='Network Graph')
# pyplot.legend(loc='upper right')
# pyplot.ylabel('size of the largest connected component')
# pyplot.xlabel('number of nodes removed')
# pyplot.title('comparison of graph resiliency to targeted attack')
#
# pyplot.show()
