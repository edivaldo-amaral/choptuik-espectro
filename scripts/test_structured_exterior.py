from fractions import Fraction as Q
import unittest

from build_tile import exact_inverse
from check_free_preconditioner import free_matrix, input_tail
from structured_exterior import structured_input_tail, report


class StructuredExteriorTests(unittest.TestCase):
    def test_improvement_and_threshold(self):
        self.assertEqual(structured_input_tail(64, 192), Q(81, 193))
        self.assertLess(structured_input_tail(64, 192), input_tail(64, 192, Q(0)))
        self.assertFalse(report()['exterior_contraction'])
        result = report(outer=(364, 1964))
        self.assertTrue(result['exterior_contraction'])
        self.assertFalse(result['physical_certificate'])

    def test_exact_finite_columns_both_scalar_blocks(self):
        for d in (0, 1):
            degrees = range(1, 10, 2) if d == 0 else range(10)
            for m in (0, 4, 16):
                addresses = [(d, p, m, n) for n in degrees
                             for p in ((0,) if m == 0 else (0, 1))]
                matrix = free_matrix(addresses, Q(1, 6))
                inverse = exact_inverse(matrix)
                weights = [(2-(n == 0))*Q(5, 4)**n for _, _, _, n in addresses]
                for col, (_, _, _, n) in enumerate(addresses):
                    if m < 4 and n < 4:
                        continue
                    norm = sum(abs(row[col])*weights[i]/weights[col]
                               for i, row in enumerate(inverse))
                    self.assertLessEqual(norm, structured_input_tail(4, 4))

    def test_real_scalar_adjacent_cancellation(self):
        mu, base = Q(1, 6), Q(3, 5)
        for d, degrees, k in ((0, [1, 3, 5], 2), (1, [0, 1, 2], 1)):
            matrix = free_matrix([(d, 0, 0, n) for n in degrees], mu)
            for i in range(len(matrix)):
                matrix[i][i] += base
            inv = exact_inverse(matrix)
            for col in range(1, len(matrix)):
                n = degrees[col]
                self.assertEqual(inv[col-1][col],
                    -2*mu*n/((mu*(n+1)+base)*(mu*(n-k+1)+base)))

    def test_invalid_parameters(self):
        with self.assertRaises(ValueError):
            structured_input_tail(0, 10)
        with self.assertRaises(ValueError):
            structured_input_tail(4, 1)


if __name__ == '__main__':
    unittest.main()
