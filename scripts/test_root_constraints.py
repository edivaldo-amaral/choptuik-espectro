from fractions import Fraction as Q
import unittest
from sharp_propagation import Poly, X
from gauge_null_profiles import power


class RootConstraintTests(unittest.TestCase):
    def test_chain_needs_derivative_of_constraint(self):
        mu = Q(1, 6)
        # A=[[-mu,1],[0,-mu]], C1=[1,0], C0=C1 A; D=0.
        a = lambda h: [-mu*h[0]+h[1], -mu*h[1]]
        c0 = lambda h: -mu*h[0]+h[1]
        c1 = lambda h: h[0]
        h0, h1 = [Q(1), Q(0)], [Q(0), Q(-1)]
        self.assertEqual([x+mu*y for x, y in zip(a(h1), h1)], [-1, 0])
        self.assertNotEqual(c0(h1)+mu*c1(h1), 0)
        self.assertEqual(c0(h1)+mu*c1(h1)+c1(h0), 0)
        for h in (h0, h1):
            self.assertEqual(c0(h)-c1(a(h)), 0)

    def test_intertwining_nonzero_constraint_map(self):
        mu = Q(1, 6)
        a = lambda h: [-mu*h[0]+h[1], -mu*h[1]]
        d = lambda h: h[1]
        for h in ([Q(1), Q(0)], [Q(0), Q(1)], [Q(2), Q(3)]):
            self.assertEqual(-mu*d(h), d(a(h)))

    def test_coordinate_gauge_has_no_polynomial_jordan_extension(self):
        mu = Q(1, 6)
        shifted = lambda f: -mu*X*f.derivative(1)
        for n in range(7):
            f = power(X, n)
            value = f
            for k in range(1, 5):
                value = shifted(value)
                self.assertEqual(value, (-mu*n)**k*f)
                self.assertEqual(value == Poly(0), n == 0)


if __name__ == '__main__':
    unittest.main()
