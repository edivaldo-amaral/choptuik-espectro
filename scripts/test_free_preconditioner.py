from fractions import Fraction as Q
import unittest

from build_tile import exact_inverse
from check_free_preconditioner import free_matrix, input_tail, multiply


class FreePreconditionerTests(unittest.TestCase):
    def test_j_matches_differentiation_of_chebyshev_polynomial(self):
        # 2*T3=8*x^3-6*x. ((1+x)*d/dx+1)(2*T3)
        # =32*x^3+24*x^2-12*x-6 =6*T0+12*T1+12*T2+8*T3.
        addresses = [(1, 0, 0, n) for n in range(4)]
        matrix = free_matrix(addresses, Q(1))
        self.assertEqual([row[3] for row in matrix], [6, 6, 6, 4])

    def test_k_matches_parity_projection(self):
        # (x*d/dx+1)(2*T3)=32*x^3-12*x=8*T3+12*T1.
        matrix = free_matrix([(0, 0, 0, 1), (0, 0, 0, 3)], Q(1))
        self.assertEqual(matrix, [[2, 6], [0, 4]])

    def test_fourier_derivative_is_half_index_with_correct_sign(self):
        addresses = [(3, 0, 3, 0), (3, 1, 3, 0)]
        self.assertEqual(free_matrix(addresses, Q(1, 6)),
                         [[Q(1, 6), Q(-3, 2)], [Q(3, 2), Q(1, 6)]])

    def test_free_inverse_restricts_exactly(self):
        full = free_matrix([(1, 0, 0, n) for n in range(5)], Q(1, 6))
        for i in range(5):
            full[i][i] += 3
        small = [row[:3] for row in full[:3]]
        self.assertEqual([row[:3] for row in exact_inverse(full)[:3]], exact_inverse(small))

    def test_preconditioning_preserves_determinant_up_to_constant(self):
        principal = free_matrix([(1, 0, 0, n) for n in range(2)], Q(1, 6))
        principal[0][0] += 1
        principal[1][1] += 1
        inverse = exact_inverse(principal)
        det = lambda a: a[0][0]*a[1][1] - a[0][1]*a[1][0]
        for s in (-1, 0, 2):
            pencil = [[Q(1+s), Q(2)], [Q(3), Q(4+s)]]
            self.assertEqual(det(multiply(pencil, inverse)), det(pencil)/det(principal))

    def test_current_meshes_fail_even_with_zero_inverse_bound(self):
        for m, n in ((4, 12), (6, 18), (8, 24), (10, 30), (12, 36)):
            self.assertGreater(24 * input_tail(m, n, Q(0)), 1)
        with self.assertRaises(ValueError):
            input_tail(4, 1, Q(0))


if __name__ == "__main__":
    unittest.main()
