from fractions import Fraction as Q
import unittest
from guarded_exterior import guarded_fields, outward_bound, ETA
from sharp_propagation import SOURCE


class GuardedExteriorTests(unittest.TestCase):
    def test_support_threshold(self):
        fields = [{(2, 4, 0): Q(3), (3, 4, 1): Q(5), (0, 5, 0): Q(7)}]
        self.assertEqual(guarded_fields(fields, (4, 8), (6, 12)),
                         [{(3, 4, 1): Q(5), (0, 5, 0): Q(7)}])
        with self.assertRaises(ValueError):
            guarded_fields(fields, (4, 8), (3, 12))

    @unittest.skipUnless((SOURCE/'GAMMA2_macro.c').exists(), 'fonte RT ausente')
    def test_low_support_vanishes_but_uncertainty_remains(self):
        fields = [{(0, 0, 0): Q(1)} for _ in range(4)]
        self.assertEqual(outward_bound(fields, (6, 18), (12, 36), Q(6), Q(0)), 0)
        error = Q(1, 10000)
        self.assertEqual(outward_bound(fields, (6, 18), (12, 36), Q(6), error),
                         36*max(ETA)/min(ETA)*error)


if __name__ == '__main__':
    unittest.main()
