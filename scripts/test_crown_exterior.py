import unittest
from fractions import Fraction as Q

from background_couplings import crown_exterior


class CrownExteriorTests(unittest.TestCase):
    def test_exact_uniform_bounds(self):
        result = crown_exterior((12, 36), (64, 192))
        self.assertEqual(Q(result['bounds']['T_from_S']), Q(2484, 35))
        self.assertEqual(Q(result['bounds']['S_from_T']), Q(2619, 191))
        self.assertFalse(result['exterior_contraction'])

    def test_delta_cancels_only_outward(self):
        a = crown_exterior((12, 36), (64, 192), distance=Q(0))
        b = crown_exterior((12, 36), (64, 192), distance=Q(2))
        self.assertEqual(a['bounds']['T_from_S'], b['bounds']['T_from_S'])
        self.assertLess(Q(a['bounds']['TT_defect']), Q(b['bounds']['TT_defect']))

    def test_exterior_alone_is_not_physical_certificate(self):
        result = crown_exterior((12, 36), (874, 2621))
        self.assertTrue(result['exterior_contraction'])
        self.assertFalse(result['physical_certificate'])

    def test_invalid_boxes(self):
        for core, outer in [((12, 36), (12, 36)), ((12, 36), (10, 40))]:
            with self.assertRaises(ValueError):
                crown_exterior(core, outer)


if __name__ == '__main__':
    unittest.main()
