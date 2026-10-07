from fractions import Fraction as Q
import unittest

from constraint_witness import violation_lower_bound


class ConstraintWitnessTests(unittest.TestCase):
    def test_strict_positive_witness(self):
        result = violation_lower_bound(Q(1, 10), Q(1, 1000), 2, Q(1, 1000), Q(1, 1000))
        self.assertTrue(result['excludes_constraint_kernel'])
        self.assertGreater(Q(result['lower']), Q(9, 100))
        self.assertFalse(result['physical_certificate'])

    def test_approximate_residual_alone_is_insufficient(self):
        result = violation_lower_bound(Q(1, 10), Q(1, 10), 1, 0, 0)
        self.assertFalse(result['excludes_constraint_kernel'])
        with self.assertRaises(ValueError):
            violation_lower_bound(-1, 0, 0, 0, 0)


if __name__ == '__main__':
    unittest.main()
