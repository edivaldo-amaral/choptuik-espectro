from fractions import Fraction as Q
import unittest
from radius_transfer import transfer_bounds


class RadiusTransferTests(unittest.TestCase):
    def test_large_box_and_failure_at_small_box(self):
        result = transfer_bounds()
        self.assertTrue(result['both_contract'])
        self.assertLess(Q(result['strong_defect']), Q(23, 100))
        self.assertLess(Q(result['weak_defect']), Q(786, 1000))
        self.assertFalse(transfer_bounds(6, 18)['both_contract'])

    def test_division_norm_conversion(self):
        self.assertEqual(Q(144, 17)/Q(40, 9), Q(162, 85))


if __name__ == '__main__':
    unittest.main()
