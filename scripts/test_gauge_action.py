from fractions import Fraction as Q
import unittest
from gauge_null_profiles import translation_profile
from sharp_propagation import (SOURCE, TAU, X, Poly, omega, linearized_omega,
                               split_plus)


@unittest.skipUnless((SOURCE/'GAMMA2_macro.c').exists(), 'fonte RT ausente')
class GaugeActionTests(unittest.TestCase):
    def test_translation_full_offshell_covariance(self):
        for mu in (Q(1, 6), Q(17, 100)):
            w = [X*(1+TAU+X*X), TAU*TAU+X, 2+TAU*X, 3+X*X+TAU]
            h = [translation_profile(p, mu) for p in w]
            expected = [translation_profile(p, mu) for p in omega(w, mu)]
            self.assertEqual(linearized_omega(w, h, mu, mu), expected)
            self.assertTrue(any(p != 0 for p in expected))

    def test_selected_and_constraint_covariance(self):
        mu = Q(1, 6)
        w = [X*(1+TAU), TAU+X, 2+X*X, TAU*X+1]
        f, c = split_plus(omega(w, mu))
        h = [translation_profile(p, mu) for p in w]
        df, dc = split_plus(linearized_omega(w, h, mu, mu))
        self.assertEqual(df, [translation_profile(p, mu)+p for p in f])
        self.assertEqual(dc, [translation_profile(p, mu) for p in c])

    def test_phase_offshell_covariance(self):
        w = [TAU*X, TAU*TAU+X, 1-X, TAU+X*X]
        h = [p.derivative(0) for p in w]
        self.assertEqual(linearized_omega(w, h, Q(1, 6), 0),
                         [p.derivative(0) for p in omega(w, Q(1, 6))])


if __name__ == '__main__':
    unittest.main()
