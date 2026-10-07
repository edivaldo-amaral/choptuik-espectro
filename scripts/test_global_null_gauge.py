from fractions import Fraction as Q
import unittest
from sharp_propagation import Poly, TAU, X
from global_null_gauge import (transport, inverse_polynomial, coordinate_profiles,
                               inverse_periodic_mode)


class GlobalNullGaugeTests(unittest.TestCase):
    def test_periodic_modes_solve_differential_equation_exactly(self):
        fr, fi = 1+X+X*X*X, 2-X*X
        mu = Q(1,6)
        for sr in (Q(1,12), mu, Q(3,4)):
            for m in (-3, 0, 4):
                for center in (-1, 1):
                    si = Q(1,7)
                    yr, yi = inverse_periodic_mode(fr, fi, mu, sr, m, center, si)
                    self.assertEqual(transport(yr, mu, sr, center)-(si+Q(m,2))*yi, fr)
                    self.assertEqual(transport(yi, mu, sr, center)+(si+Q(m,2))*yr, fi)

    def test_periodic_resonance_and_alias_obstruction(self):
        mu = Q(1,6)
        for m in (-2, 0, 3):
            with self.assertRaisesRegex(ValueError, 'obstrucao Floquet'):
                inverse_periodic_mode(1, 0, mu, mu, m, 1, -Q(m,2))
            yr, yi = inverse_periodic_mode((X-1)*(1+X*X), 0,
                                           mu, mu, m, 1, -Q(m,2))
            self.assertEqual(transport(yr, mu, mu, 1), (X-1)*(1+X*X))
            self.assertEqual(yi, 0)
            # A cadeia T y1=-y0 e obstruida para o modo constante do nucleo.
            with self.assertRaises(ValueError):
                inverse_periodic_mode(-1, 0, mu, mu, m, 1, -Q(m,2))

    def test_fixed_cone_trace_and_periodic_reflection(self):
        mu, sr = Q(1,6), Q(3,4)
        fr, fi = (X-1)*(1+X*X), (X-1)*X
        ym = inverse_periodic_mode(fr, fi, mu, sr, 2, 1)
        yp = inverse_periodic_mode(fr.parity(), fi.parity(), mu, sr, 2, -1)
        for minus, plus in zip(ym, yp):
            self.assertEqual(plus, minus.parity())
            self.assertEqual(sum(minus.terms.values()), 0)  # xi=1
            a, b = coordinate_profiles(minus, plus, mu)
            self.assertEqual(a.parity(), a)
            self.assertEqual(b.parity(), -b)

    def test_periodic_input_rejects_nonperiodic_polynomial(self):
        with self.assertRaises(ValueError):
            inverse_periodic_mode(TAU, 0, Q(1,6), Q(3,4), 0, 1)

    def test_transport_inverse_on_both_sides_of_mu(self):
        source = 1+TAU*TAU+X*X*X+TAU*X
        for s in (Q(3,4), Q(1,12)):
            for center in (-1,1):
                y = inverse_polynomial(source, Q(1,6), s, center)
                self.assertEqual(transport(y, Q(1,6), s, center), source)

    def test_resonance_obstruction_and_kernel(self):
        mu = Q(1,6)
        self.assertEqual(transport(Poly(1), mu, mu, 1), 0)
        with self.assertRaises(ValueError):
            inverse_polynomial(Poly(1), mu, mu, 1)
        source = (X-1)*(1+TAU*TAU+X*X)
        y = inverse_polynomial(source, mu, mu, 1)
        self.assertEqual(transport(y, mu, mu, 1), source)

    def test_reflection_preserves_center(self):
        mu, s = Q(1,6), Q(3,4)
        source = 1+X+TAU*X*X
        ym = inverse_polynomial(source, mu, s, 1)
        yp = inverse_polynomial(source.parity(), mu, s, -1)
        self.assertEqual(yp, ym.parity())
        a, b = coordinate_profiles(ym, yp, mu)
        self.assertEqual(a.parity(), a)
        self.assertEqual(b.parity(), -b)


if __name__ == '__main__':
    unittest.main()
