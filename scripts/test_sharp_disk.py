from fractions import Fraction as Q
import tempfile
from pathlib import Path
import unittest
from certify_sharp_disk import certify, check_compression
from build_contour import read_spectral_matrix


class SharpDiskTests(unittest.TestCase):
    def test_preconditioned_bound_uses_right_residual(self):
        result = certify([[2]], 0, [Q(1)], Q(0), Q(1), principal=[[Q(3)]])
        self.assertLess(Q(result['right_uniform_defect']), 1)
        self.assertLess(Q(result['preconditioned_inverse_bound']), Q(31,10))
        self.assertGreaterEqual(Q(result['preconditioned_inverse_bound']), 3)

    def test_nested_compression_exact_and_mismatch(self):
        addresses = [(1, 0, 0, 0), (0, 0, 0, 1)]
        result = check_compression([[8, 3], [4, 5]], -2, addresses,
                                   [[2]], 0, addresses[:1])
        self.assertEqual(result['shell_dimension'], 1)
        with self.assertRaises(ValueError):
            check_compression([[8, 3], [4, 5]], -2, addresses,
                              [[3]], 0, addresses[:1])

    def test_scalar_accepted_and_rejected(self):
        self.assertTrue(certify([[2]], 0, [Q(1)], Q(0), Q(1))['accepted'])
        self.assertFalse(certify([[2]], 0, [Q(1)], Q(0), Q(3))['accepted'])

    def test_background_error_is_not_dropped(self):
        result = certify([[2]], 0, [Q(1)], Q(0), Q(0), background_error=Q(3))
        self.assertTrue(result['background_error_included'])
        self.assertFalse(result['accepted'])
        with self.assertRaises(ValueError):
            certify([[2]], 0, [Q(1)], background_error=Q(-1))

    def test_sharp_magic_is_explicit(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'matrix.dat'
            path.write_text('CHOPTUIK_RT_SHARP_MATRIX_V1\ndimension 1\nBEGIN_ENTRIES\n0 0 2 0\n')
            with self.assertRaises(SystemExit):
                read_spectral_matrix(path)
            self.assertEqual(read_spectral_matrix(path, 'CHOPTUIK_RT_SHARP_MATRIX_V1'),
                             (1, [[2]], 0))


if __name__ == '__main__':
    unittest.main()
