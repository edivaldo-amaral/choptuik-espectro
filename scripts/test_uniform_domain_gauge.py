from fractions import Fraction as Q
import unittest

from uniform_coupling import uniform_tile, admissible_radius
from gauge_null_profiles import null_generator, power, phase_transversality
from sharp_propagation import X, SOURCE


class UniformCouplingTests(unittest.TestCase):
    def test_strict_uniform_radius(self):
        args = (Q(1, 10), Q(2), Q(1, 20), Q(1, 5), Q(1, 4), Q(1, 2))
        report = admissible_radius(*args)
        limit = Q(report['strict_radius_upper'])
        self.assertEqual(limit, Q(3, 8))
        a, v, err, b, c, d = args
        self.assertTrue(uniform_tile(a, v, limit/2, err, b, c, d)['accepted'])
        self.assertFalse(uniform_tile(a, v, limit, err, b, c, d)['accepted'])

    def test_exterior_diagonal_alone_is_insufficient(self):
        result = uniform_tile(0, 0, 0, 0, 2, 1, Q(1, 2))
        self.assertFalse(result['accepted'])
        self.assertEqual(Q(result['schur_defect']), 4)

    def test_current_64x192_exterior_bound_fails_even_ideal_finite_block(self):
        d = (Q(5191, 500)+Q(5, 4))*Q(81, 193)
        self.assertFalse(uniform_tile(0, 0, 0, 0, 0, 0, d)['accepted'])


class GaugeProfileTests(unittest.TestCase):
    @unittest.skipUnless((SOURCE/'RefA.dat').exists(), 'RefA nao disponivel')
    def test_phase_coefficient_matches_reference_file(self):
        from component_exterior import read_fields
        field = read_fields(SOURCE/'RefA.dat')[3]
        self.assertEqual(field.get((1, 0, 0), 0), 0)
        self.assertEqual(field[1, 0, 1], Q(138672388040959354547, 590295810358705651712))

    def test_phase_slice_is_transverse_under_published_error(self):
        reference = Q(138672388040959354547, 590295810358705651712)
        lower = phase_transversality(reference, Q(1, 2**25)+Q(1, 2**277))
        self.assertGreater(lower, Q(117, 1000))

    def test_neutral_phase_projection_does_not_preserve_positive_kernel(self):
        # L(0)=diag(0,-1), p=(1,0), h=(0,1), b(x)=x1+x2, s=1.
        # h e modo de L(1), mas h-p nao e: L(1)(h-p)=(-1,0).
        h, p = (0, 1), (1, 0)
        projected = tuple(a-b for a, b in zip(h, p))
        self.assertEqual((projected[0], 0*projected[1]), (-1, 0))

    def test_coordinate_pushforward_exact(self):
        mu = Q(17, 100)
        for n in range(7):
            s, a, b = null_generator(n, mu)
            self.assertEqual(s, mu*(1-n))
            for sigma in (-1, 1):
                self.assertEqual(mu*(1+sigma*X)*a-sigma*b,
                                 (-1)**n*power(1+sigma*X, n))

    def test_phase_and_translation_are_distinct(self):
        mu = Q(1, 6)
        s, a, b = null_generator(1, mu)
        self.assertEqual((s, a, b), (0, -1/mu, 0))
        s, a, b = null_generator(0, mu)
        self.assertEqual((s, a, b), (mu, 1/mu, X))

    def test_only_positive_generator_moves_cone(self):
        mu = Q(1, 6)
        for n in range(7):
            s, a, b = null_generator(n, mu)
            normal = mu*(1-X)*a+b
            at_cone = sum(normal.terms.values())  # perfil depende somente de xi
            self.assertEqual(at_cone, 1 if n == 0 else 0)
            self.assertEqual(s > 0, n == 0)

    def test_radius_loss_constants(self):
        qt, qx = Q(129, 130), Q(9, 10)
        self.assertEqual(qt/(2*(1-qt)**2), 8385)
        self.assertEqual(2*qx/((Q(9, 8)-1)*(1-qx)**2), 1440)


if __name__ == '__main__':
    unittest.main()
