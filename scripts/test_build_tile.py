"""Verificacao exata da identidade que liga tiles ao pencil espectral."""
from fractions import Fraction as Q
import unittest

from build_tile import coefficient_matrix, exact_inverse, round_to_grid, one_norm


class ResolventTests(unittest.TestCase):
    def test_weighted_resolvent_uses_weighted_column_norm(self):
        matrix = [[2, 1], [1, 2]]
        numerators, denominator = coefficient_matrix(matrix, 0, Q(0),
                                                    weights=[Q(1), Q(10)])
        self.assertEqual([[Q(x, denominator) for x in row] for row in numerators],
                         [[Q(2, 3), Q(-1, 30)], [Q(-10, 3), Q(2, 3)]])
        self.assertEqual(one_norm(numerators, denominator), Q(4))

    def test_relative_pencil_identity(self):
        matrix = [[2, 1], [3, 4]]
        for exponent in (-3, 0, 2):
            base, delta = Q(1, 3), Q(2, 7)
            numerators, denominator = coefficient_matrix(matrix, exponent, base)
            resolvent = [[Q(value, denominator) for value in row] for row in numerators]
            a = [[value * Q(2) ** exponent for value in row] for row in matrix]
            for i in range(2):
                for j in range(2):
                    product = sum((a[i][k] + (base + delta if i == k else 0))
                                  * resolvent[k][j] for k in range(2))
                    self.assertEqual(product, Q(i == j) + delta * resolvent[i][j])

    def test_singular_base_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "singular"):
            coefficient_matrix([[-1, 0], [0, 2]], 0, Q(1))

    def test_compression_and_inversion_do_not_commute(self):
        full = exact_inverse([[Q(2), Q(1)], [Q(1), Q(2)]])
        compressed = exact_inverse([[Q(2)]])
        self.assertEqual(full[0][0], Q(2, 3))
        self.assertEqual(compressed[0][0], Q(1, 2))
        self.assertNotEqual(full[0][0], compressed[0][0])

    def test_rational_center_rounding_residual(self):
        for value in (Q(1, 3), Q(-5, 7), Q(9, 2)):
            grid = 16
            rounded = Q(round_to_grid(value.numerator * grid, value.denominator), grid)
            self.assertLessEqual(abs(value - rounded), Q(1, 2 * grid))


if __name__ == "__main__":
    unittest.main()
