from fractions import Fraction as Q
import unittest

from make_winding_certificate import Box, Edge, coarsen_edges


class CoarseningTests(unittest.TestCase):
    def test_merges_straight_segments_and_keeps_corners(self):
        box = Box(Q(1), Q(2), Q(-1), Q(1))
        points = [(0, 0), (1, 0), (2, 0), (2, 1), (0, 1), (0, 0)]
        edges = [Edge(box, box, box, 0, a, b) for a, b in zip(points, points[1:])]
        result = coarsen_edges(edges)
        self.assertEqual(len(result), 4)
        self.assertEqual(result[0].start_point, (0, 0))
        self.assertEqual(result[0].end_point, (2, 0))
        self.assertEqual(sum(e.increment for e in result), 0)
        for edge in result:
            self.assertTrue(edge.box.contains(edge.start))
            self.assertTrue(edge.box.contains(edge.end))

    def test_backtracking_is_not_replaced_by_chord(self):
        box = Box(Q(1), Q(2), Q(-1), Q(1))
        points = [(0, 0), (2, 0), (1, 0), (1, 1), (0, 0)]
        edges = [Edge(box, box, box, 0, a, b) for a, b in zip(points, points[1:])]
        self.assertEqual(coarsen_edges(edges), edges)


if __name__ == "__main__":
    unittest.main()
