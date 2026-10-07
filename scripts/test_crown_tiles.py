from fractions import Fraction as Q
import unittest
from unittest.mock import patch

import numpy as np

from certify_crown_tiles import (tile, nearest_integer, ceil_grid,
                                 sparse_right_product, cover_upper_boundary, exact_product)


class CrownTilesTests(unittest.TestCase):
    def test_rounding_bounds(self):
        for q in (Q(1, 3), Q(-1, 2), Q(-5, 2), Q(7, 2)):
            self.assertLessEqual(abs(nearest_integer(q)-q), Q(1, 2))
            self.assertGreaterEqual(ceil_grid(q, 8), q)

    def test_sparse_product_exact(self):
        a, b = [[2, 0], [3, 7]], [[0, 5], [11, 0]]
        self.assertEqual(sparse_right_product(a, b).tolist(), [[0, 10], [77, 15]])

    def test_limb_product_large_signed_integers(self):
        a = np.array([[2**130+7, -2**90, 0], [-31, 2**61+3, 2**75]], dtype=object)
        b = np.array([[-2**85, 3], [2**93+1, -2**111], [0, 2**63]], dtype=object)
        self.assertEqual(exact_product(a, b).tolist(), (a@b).tolist())
        self.assertEqual(exact_product(a*0, b).tolist(), [[0, 0], [0, 0]])

    def test_sparse_left_and_right_limb_products(self):
        a = np.eye(16, dtype=object)*(2**105+7)
        a[0, 3] = -2**75
        b = np.array([[(-1)**(i+j)*(2**80+3*i+j) for j in range(16)]
                      for i in range(16)], dtype=object)
        self.assertEqual(exact_product(a, b).tolist(), (a@b).tolist())
        self.assertEqual(exact_product(b, a).tolist(), (b@a).tolist())

    def test_scalar_tile_and_singular_rejection(self):
        grid = 2**48
        data = dict(h=np.array([[2*grid]], dtype=object),
                    q=np.array([[grid]], dtype=object), grid=grid,
                    qerror=Q(0), anorm=Q(2), rounding=Q(0), size=1)
        result = tile(data, Q(1), Q(1, 4))
        self.assertTrue(result['accepted'])
        self.assertLess(Q(result['uniform_defect']), 1)
        self.assertFalse(tile(data, Q(-2), Q(0))['accepted'])

    def test_coverage_has_no_gaps(self):
        def fake_tile(data, real, imag):
            return dict(center=[str(real), str(imag)], accepted=True, radius='1/8')
        with patch('certify_crown_tiles.tile', fake_tile), patch('builtins.print'):
            tiles, segments, covered = cover_upper_boundary({})
        self.assertTrue(covered)
        self.assertEqual(segments[0]['start'], ['1/8', '0'])
        self.assertEqual(segments[-1]['end'], ['1', '0'])
        for a, b in zip(segments, segments[1:]):
            self.assertEqual(a['end'], b['start'])
        for segment in segments:
            distance = sum(abs(Q(a)-Q(b)) for a, b in zip(segment['start'], segment['end']))
            self.assertLessEqual(distance, Q(tiles[segment['tile']]['radius']))


if __name__ == '__main__':
    unittest.main()
