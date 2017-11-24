import graph
import copy
import time
from bfs import networkimport


UNDIRECTED_EX_GRAPH0 = {0: {1,2,3,4}, 1: {0,3,4}, 2: {0,4}, 3: {0,1}, 4: {0,1,2}}
EX_GRAPH0 = {0: {1,2}, 1: {}, 2: {}}
EX_GRAPH1 = {0: {1,4,5}, 1: {2,6}, 2: {3}, 3: {0}, 4: {1}, 5: {2}, 6: {}}
EX_GRAPH2 = {0: {1,4,5}, 1: {2,6}, 2: {3,7}, 3: {7}, 4: {1}, 5: {2}, 6: {}, 7: {3}, 8: {1,2}, 9: {0,3,4,5,6,7}}
import unittest

class GraphTest(unittest.TestCase):
    # def test_make_complete_graph(self):
    #     giraffe = graph.make_complete_graph(10)
    #     self.assertEqual(len(giraffe), 10)
    #     self.assertEqual(len(giraffe[0]), 9)
    #     self.assertFalse(0 in giraffe[0])
    #
    # def test_make_complete_graph2(self):
    #     giraffe = graph.make_complete_graph(0)
    #     self.assertEqual(len(giraffe), 0)
    #
    # def test_compute_in_degrees(self):
    #     giraffe = graph.compute_in_degrees(EX_GRAPH1)
    #     self.assertEqual(len(giraffe), 7)
    #     self.assertEqual(giraffe[0], 1)
    #     self.assertEqual(giraffe[2], 2)
    #
    # def test_compute_in_degrees2(self):
    #     giraffe = graph.compute_in_degrees(graph.make_complete_graph(0))
    #     self.assertEqual(len(giraffe), 0)
    #
    # def test_in_degree_distribution(self):
    #     giraffe = graph.in_degree_distribution(EX_GRAPH0)
    #     self.assertEqual(len(giraffe), 2)
    #     self.assertEqual(giraffe[1], 2)
    #     self.assertEqual(giraffe[0], 1)
    #
    # def test_sum_in_degrees(self):
    #     sum = graph.sum_in_degrees(graph.compute_in_degrees(EX_GRAPH1))
    #     self.assertEqual(sum, 9)
    #
    # def test_avg_in_degree(self):
    #     avg = graph.avg_in_degree(graph.in_degree_distribution(EX_GRAPH0))
    #     self.assertEqual(avg, 2/float(3))
    #
    # def test_make_undirected_graph(self):
    #     # giraffe = graph.make_undirected_graph(5,.1)
    #     giraffe2 = graph.make_undirected_graph(1347, .001715)   # make the undirected graph with 1327 nodes, ~6200 edges
    #     sum = 0
    #     for value in giraffe2.values():
    #         sum += len(value)
    #     print(sum)
    #
    # def test_network_import(self):
    #     giraffe = networkimport.load_graph(networkimport.NETWORK_URL) # make the undirected network graph with 1327 nodes, ~6200 edges
    #     sum = 0
    #     for value in giraffe.values():
    #         sum += len(value)
    #     print(sum)
    #
    # def test_make_upa_graph(self):
    #     giraffe = graph.upa(1347, 2)    # make the undirected upa graph with 1327 nodes, ~5300 edges
    #     sum = 0
    #     for value in giraffe.values():
    #         sum += len(value)
    #     print(sum)
    #
    # def test_random_order(self):
    #     list = graph.random_order(copy.copy(EX_GRAPH0))
    #     self.assertEqual(len(list), len(EX_GRAPH0.keys()))

    def test_fast_targeted_order(self):
        sum = 0;
        trials = 100;
        node_count_list = []
        time_list = []
        for i in range(10, 1000, 10):
            giraffe = graph.upa(i, 5)
            start = time.clock()
            targeted = graph.fast_targeted_order(giraffe)
            end = time.clock()
            node_count_list.append(i)
            time_list.append(end-start)
        print(time_list)
        print(node_count_list)

    def test_targeted_order(self):
        sum = 0;
        trials = 100;
        node_count_list = []
        time_list = []
        for i in range(10, 1000, 10):
            giraffe = graph.upa(i, 5)
            start = time.clock()
            targeted = networkimport.targeted_order(giraffe)
            end = time.clock()
            node_count_list.append(i)
            time_list.append(end-start)
        print(time_list)
        print(node_count_list)

