from fractions import Fraction as Q
import unittest
from sharp_exterior import free_tail, sharp_majorant, finite_background_error
from build_tile import exact_inverse
from check_free_preconditioner import free_matrix
from sharp_propagation import SOURCE


class SharpExteriorTests(unittest.TestCase):
    def test_weak_tail_and_exact_columns(self):
        self.assertEqual(free_tail(64, 192), Q(153, 193))
        for d in (0, 1):
            for m in (0, 4, 16):
                degrees = range(1, 12, 2) if d == 0 else range(12)
                addresses = [(d, p, m, n) for n in degrees
                             for p in ((0,) if m == 0 else (0, 1))]
                inv = exact_inverse(free_matrix(addresses, Q(1, 6)))
                weights = [(2-(n == 0))*Q(9, 8)**n for _, _, _, n in addresses]
                for j, (_, _, _, n) in enumerate(addresses):
                    if m >= 4 or n >= 4:
                        value = sum(abs(row[j])*weights[i]/weights[j] for i, row in enumerate(inv))
                        self.assertLessEqual(value, free_tail(4, 4))

    @unittest.skipUnless((SOURCE/'GAMMA2_SHARP_macro.c').exists(), 'fonte RT ausente')
    def test_constant_background_factor_one(self):
        fields = [{}, {(0, 0, 0): Q(1, 4)}, {}, {}]
        self.assertEqual(sharp_majorant(fields, mu_max=Q(0)),
                         [[Q(1, 4), Q(1, 4)], [Q(0), Q(1, 2)]])

    def test_finite_error_removes_temporal_part(self):
        addresses = [(1, p, 100, 0) for p in (0, 1)]
        result = finite_background_error(addresses, [Q(1), Q(1)], Q(1), Q(0))
        self.assertEqual(result, 1+Q(144, 17))

    def test_invalid_parameters(self):
        with self.assertRaises(ValueError):
            free_tail(0, 2)
        with self.assertRaises(ValueError):
            free_tail(2, 2, radial=Q(1))


if __name__ == '__main__':
    unittest.main()
