from fractions import Fraction as Q
from pathlib import Path
import tempfile
import unittest

from reference_error_budget import read_mu, error_budget


class ReferenceErrorTests(unittest.TestCase):
    def test_mu_error_is_not_discarded_when_omega_error_zero(self):
        result = error_budget(Q(1, 6), Q(1, 6)+Q(1, 10**10), Q(0))
        self.assertGreater(Q(result['total_H_error']), Q(0))
        self.assertGreater(Q(result['total_H_error']), Q(result['principal_relative_error']))
        self.assertFalse(result['physical_certificate'])

    def test_reader_uses_only_header(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'mu.dat'
            path.write_text('mu_MultiField\nDyadicQ\nTwoExp -3\ncoeff 1\nignored tail\n')
            self.assertEqual(read_mu(path), Q(1, 8))
            path.write_text('wrong\nDyadicQ\nTwoExp -3\ncoeff 1\n')
            with self.assertRaises(ValueError):
                read_mu(path)


if __name__ == '__main__':
    unittest.main()
