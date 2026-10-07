from fractions import Fraction as Q
import unittest
from sharp_guarded_bridge import schur_bound


class SharpBridgeTests(unittest.TestCase):
    def test_asymmetric_coupling_and_rejection(self):
        self.assertEqual(schur_bound(Q(10), Q(1,1000), Q(1,2), Q(1,2)), Q(101,200))
        self.assertGreater(schur_bound(Q(10), Q(1), Q(1,2), Q(1,2)), 1)
        with self.assertRaises(ValueError):
            schur_bound(Q(-1), Q(1), Q(1), Q(0))
