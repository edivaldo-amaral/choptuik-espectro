from pathlib import Path
import subprocess
import tempfile
import unittest

from analyze_spectrum import read_matrix
from analyze_sharp import constraint_addresses

ROOT = Path(__file__).resolve().parents[1]
BINARY = ROOT/'build/sharp-matrix'


@unittest.skipUnless(BINARY.exists(), 'compile primeiro com bash scripts/build_sharp.sh')
class SharpExportTests(unittest.TestCase):
    def test_zero_background_principal_and_constraint_coupling(self):
        field = 'Field\nTwoExp 0\noff_m 0\noff_n 0\nnum_m 1\nnum_n 1\n(((0,0)))\n'
        with tempfile.TemporaryDirectory() as directory:
            background = Path(directory)/'zero.dat'
            output = Path(directory)/'sharp.dat'
            background.write_text('mu_MultiField\nDyadicQ\nTwoExp -3\ncoeff 1\nMultiField\nnum_comp 4\n'+field*4)
            command = [str(BINARY), str(background), '1', '4', str(output)]
            subprocess.run(command, check=True, capture_output=True)
            with self.assertRaises(ValueError):
                read_matrix(output)  # nao deixar sharp entrar no analisador selecionado
            matrix, addresses = read_matrix(output, 'CHOPTUIK_RT_SHARP_MATRIX_V1')
            lookup = {a: i for i, a in enumerate(addresses)}
            self.assertEqual(len(matrix), 6)
            for i, (_, _, _, n) in enumerate(addresses):
                self.assertEqual(matrix[i, i], (n+1)/8)
            self.assertEqual(matrix[lookup[1, 0, 0, 0], lookup[0, 0, 0, 1]], 1/4)
            self.assertEqual(matrix[lookup[1, 0, 0, 0], lookup[0, 0, 0, 3]], -1/4)
            self.assertEqual(matrix[lookup[1, 0, 0, 2], lookup[0, 0, 0, 3]], 1/4)
            previous = output.read_bytes()
            result = subprocess.run(command, capture_output=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertEqual(output.read_bytes(), previous)

    def test_reject_invalid_cutoff(self):
        result = subprocess.run([str(BINARY), 'unused', '-1', '4', 'unused'], capture_output=True)
        self.assertEqual(result.returncode, 2)

    def test_sharp_background_uses_factor_one(self):
        field = 'Field\nTwoExp -2\noff_m 0\noff_n 0\nnum_m 1\nnum_n 1\n(((0,0)))\n'
        fields = [field, field.replace('(((0,0)))', '(((1,0)))'), field, field]
        with tempfile.TemporaryDirectory() as directory:
            background, output = Path(directory)/'constant.dat', Path(directory)/'sharp.dat'
            background.write_text('mu_MultiField\nDyadicQ\nTwoExp -3\ncoeff 1\nMultiField\nnum_comp 4\n'+''.join(fields))
            subprocess.run([str(BINARY), str(background), '1', '4', str(output)],
                           check=True, capture_output=True)
            matrix, addresses = read_matrix(output, 'CHOPTUIK_RT_SHARP_MATRIX_V1')
            for i, (component, _, _, n) in enumerate(addresses):
                # w=(0,1/4,0,0): Gamma_sharp diagonals 1/4 and 1/2.
                self.assertEqual(matrix[i, i], (n+1)/8+(1/4 if component == 0 else 1/2))


class SharpAddressTests(unittest.TestCase):
    def test_bad_address_order_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'constraints'
            path.write_text('BEGIN_ADDRESSES\n1 0 0 0 1\nBEGIN_ENTRIES\n')
            with self.assertRaises(ValueError):
                constraint_addresses(path)


if __name__ == '__main__':
    unittest.main()
