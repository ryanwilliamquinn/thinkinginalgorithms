from bfs import BreadthFirstSearch
import copy

EX_GRAPH0 = {0: {1, 2}, 1: {0}, 2: {0}, 3: set([]), 4: {5, 6}, 5: {4}, 6: {4, 7, 8}, 7: {6, 8}, 7: {7, 8}, 8: {6, 7}}
EX_GRAPH1 = {3: set([2]),
             0: set([1]),
             2: set([1, 3]),
             1: set([0, 2])}

import unittest


class BfsTest(unittest.TestCase):
    def test_bfs_visited(self):
        giraffe = BreadthFirstSearch.bfs_visited(EX_GRAPH0, 0)
        self.assertEqual(len(giraffe), 3)
        self.assertTrue(0 in giraffe)
        self.assertTrue(1 in giraffe)
        self.assertTrue(0 in giraffe)

    def test_cc_visited(self):
        giraffe = BreadthFirstSearch.cc_visited(EX_GRAPH0)
        self.assertEqual(len(giraffe), 3)
        self.assertTrue(set([4, 5, 6, 7, 8]) in giraffe)

    def test_largest_cc_size(self):
        self.assertTrue(BreadthFirstSearch.largest_cc_size(EX_GRAPH0), 5)
        self.assertTrue(BreadthFirstSearch.largest_cc_size(EX_GRAPH1), 4)


    def test_compute_resilience(self):
        longest_ccc_sizes = BreadthFirstSearch.compute_resilience(EX_GRAPH1, [1,2])
        self.assertEqual(longest_ccc_sizes, [4,2,1])


