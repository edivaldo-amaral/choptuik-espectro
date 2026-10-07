from fractions import Fraction as Q
from pathlib import Path
import tempfile
import unittest

from background_couplings import (guard_tail, tail, reverse_coupling,
                                  weighted_coefficients, RT_BACKGROUND_ERROR)
from check_frequency_blocks import audit_blocks, two_block_test, transport_defect


class BlockBoundsTests(unittest.TestCase):
    def test_transport_component_norm(self):
        bounds = [[Q(1, 10**9), Q(2, 10**9)], [Q(3, 10**9), Q(1, 10**9)]]
        new = transport_defect(bounds, [Q(1)]*4, list(map(Q, [916, 4096, 243, 1233])))
        self.assertEqual(new[0][0], Q(4096, 243)*bounds[0][0])
        self.assertTrue(two_block_test(new[0][0], new[0][1], new[1][0], new[1][1])['accepted'])

    def test_common_rescaling_of_norm_costs_nothing(self):
        bounds = [[Q(1, 10)]]
        self.assertEqual(transport_defect(bounds, [Q(2), Q(3)], [Q(8), Q(12)]), bounds)
        with self.assertRaises(ValueError):
            transport_defect(bounds, [Q(0)], [Q(1)])

    def test_asymmetric_couplings_pass_after_scaling(self):
        result = two_block_test(Q(1,10), Q(4), Q(1,100), Q(1,5))
        self.assertTrue(result["accepted"])
        r = Q(result["weight"])
        self.assertLess(Q(1,10)+r/100, 1)
        self.assertLess(Q(1,5)+4/r, 1)

    def test_diagonal_margins_alone_are_not_sufficient(self):
        self.assertFalse(two_block_test(Q(0), Q(2), Q(1), Q(0))["accepted"])
        self.assertFalse(two_block_test(Q(2), Q(0), Q(0), Q(2))["accepted"])

    def test_large_coupling_does_not_make_matrix_singular(self):
        matrix = [[Q(1), Q(2)], [Q(-2), Q(1)]]
        separate = audit_blocks(matrix, [0], [1], 32, True)
        coupled = audit_blocks(matrix, [0], [1], 32, True, coupled=True)
        self.assertFalse(separate["accepted"])
        self.assertTrue(coupled["accepted"])

    def test_guard_uses_exact_first_omitted_index(self):
        coeffs = [(2,0,Q(1)), (3,0,Q(2)), (0,7,Q(3))]
        self.assertEqual(guard_tail(coeffs, (4,12), (6,18)), Q(5))
        self.assertEqual(tail(coeffs, 4,8), 0)
        self.assertGreater(6*(tail(coeffs,4,8)+RT_BACKGROUND_ERROR), 0)

    def test_reverse_passage_remains_when_background_tail_vanishes(self):
        result = reverse_coupling([], (4,12), (12,36), (24,72))
        self.assertGreater(Q(result["q_middle_from_far"]), 0)
        self.assertGreater(Q(result["reverse_H_bound"]), 0)
        with self.assertRaises(ValueError):
            reverse_coupling([], (4,12), (3,36), (24,72))

    def test_background_parser_keeps_multiplicities_and_dyadic_scale(self):
        block = "Field\nTwoExp -1\noff_m 1\noff_n 1\nnum_m 1\nnum_n 1\n((3,4))\n"
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory)/"background.dat"
            source.write_text("mu_MultiField\nMultiField\n"+block*4, encoding="ascii")
            coefficients = weighted_coefficients(source)
        self.assertEqual(coefficients, [(1,1,Q(2275,128))]*4)


if __name__ == "__main__":
    unittest.main()
