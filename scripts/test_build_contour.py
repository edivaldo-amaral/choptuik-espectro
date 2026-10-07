"""Regressoes exatas da reconstrucao CRT; executar por unittest."""
import math
import unittest
from fractions import Fraction as Q

from build_contour import _coefficient_bound, bareiss_determinant, charpoly_exact
from build_contour import rectangle_nodes, integer_lattice_nodes


class ExactPolynomialTests(unittest.TestCase):
    def test_sector_b_contour_is_exact_translation(self):
        original = rectangle_nodes(Q(1, 8), Q(1), Q(1, 4), 4, 3)
        shifted = rectangle_nodes(Q(1, 8), Q(1), Q(1, 4), 4, 3, Q(1, 2))
        self.assertEqual(shifted, [(x, y + Q(1, 2)) for x, y in original])
        self.assertEqual(min(y for _, y in shifted), Q(1, 4))
        self.assertEqual(max(y for _, y in shifted), Q(3, 4))
        for exponent in (-71, 0, 3):
            scale, integer_nodes = integer_lattice_nodes(shifted, exponent)
            for (x, y), (ix, iy) in zip(shifted, integer_nodes):
                self.assertEqual(Q(ix, scale) * Q(2) ** exponent, x)
                self.assertEqual(Q(iy, scale) * Q(2) ** exponent, y)

    def test_hadamard_bound_rounds_outward(self):
        # k=3 tem raiz irracional; B grande amplifica o arredondamento antigo.
        for dimension in (1, 2, 3, 5, 8):
            for largest in (0, 1, 10, 1000):
                matrix = [[largest] * dimension for _ in range(dimension)]
                bound = _coefficient_bound(matrix, dimension)
                for size in range(1, dimension + 1):
                    squared = (math.comb(dimension, size) ** 2
                               * size ** size * largest ** (2 * size))
                    self.assertGreaterEqual(bound * bound, squared)

    def test_crt_agrees_with_independent_determinants(self):
        matrices = [
            [[0]],
            [[0, 2], [3, -1]],
            [[0, 2, -3], [0, 4, 5], [7, -8, 9]],
            [[1000, -2000, 3000], [4000, 5000, -6000],
             [-7000, 8000, 9000]],
        ]
        for matrix in matrices:
            dimension = len(matrix)
            coefficients = charpoly_exact(matrix, dimension, verbose=False)
            # n+1 pontos distintos determinam um polinomio de grau <= n.
            for point in range(-2, dimension - 1):
                shifted = [[value + (point if row == column else 0)
                            for column, value in enumerate(values)]
                           for row, values in enumerate(matrix)]
                evaluated = sum(value * point ** degree
                                for degree, value in enumerate(coefficients))
                self.assertEqual(evaluated,
                                 bareiss_determinant(shifted, dimension))


if __name__ == "__main__":
    unittest.main()
