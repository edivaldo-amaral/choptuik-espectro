from fractions import Fraction as Q
import unittest

from sharp_propagation import (SOURCE, Poly, TAU, X, SHIFT, omega,
    linearized_omega, split_plus, inject_constraints, lift_selected,
    complete, gamma_sharp_full, propagation, forcing, j)


@unittest.skipUnless((SOURCE/'GAMMA2_macro.c').exists(), 'fonte RT nao disponivel')
class SharpPropagationTests(unittest.TestCase):
    def setUp(self):
        self.mu = Q(1, 6)
        self.w = [X*(1+TAU), 1+X+TAU*X, 2+X*X, TAU+X]
        self.h = [X*(2+TAU*TAU+X*X), TAU+2*X*X, 1-X+TAU, 3+TAU*X]

    def test_full_nonlinear_identity_off_shell(self):
        q = omega(self.w, self.mu)
        self.assertTrue(any(p != 0 for p in q))
        self.assertEqual(complete(self.w, [-p.parity() for p in q], self.mu), [0, 0])

    def test_selected_constraint_reconstruction(self):
        for q in (omega(self.w, self.mu),
                  linearized_omega(self.w, self.h, self.mu)):
            f, c = split_plus(q)
            self.assertEqual([a+b for a, b in zip(inject_constraints(c), lift_selected(f))],
                             [-p.parity() for p in q])

    def test_full_floquet_linearized_identity_with_background_defect(self):
        q = omega(self.w, self.mu)
        dq = linearized_omega(self.w, self.h, self.mu)
        e = [-p.parity() for p in q]
        de = [-p.parity() for p in dq]
        correction = gamma_sharp_full(self.h, e)
        # s e^{s tau} surge em TODAS as derivadas do residuo linearizado.
        result = complete(self.w, de, self.mu, SHIFT)
        self.assertEqual([a+X*b for a, b in zip(result, correction)], [0, 0])
        self.assertTrue(any(p != 0 for p in correction))

    def test_explicit_forcing_sign(self):
        f, c = split_plus(linearized_omega(self.w, self.h, self.mu))
        e = [-p.parity() for p in omega(self.w, self.mu)]
        correction = gamma_sharp_full(self.h, e)
        k = propagation(self.w, c, self.mu)
        m = forcing(self.w, f, self.mu)
        self.assertEqual([a-b+d for a, b, d in zip(k, m, correction)], [0, 0])

    def test_exact_zero_background_has_no_correction(self):
        w = [Poly() for _ in range(4)]
        f, c = split_plus(linearized_omega(w, self.h, self.mu))
        self.assertEqual(propagation(w, c, self.mu), forcing(w, f, self.mu))

    def test_expanded_source_operator(self):
        f, _ = split_plus(linearized_omega(self.w, self.h, self.mu))
        g = gamma_sharp_full(self.w, lift_selected(f))
        pf = [p.parity() for p in f]
        expanded = [-self.mu*(X*pf[0].derivative(1)+2*pf[0])-g[0],
                    j(X*pf[1], self.mu, SHIFT).parity()
                    -self.mu*(pf[0]+pf[1]-f[1]+2*f[2])-g[1]]
        self.assertEqual(forcing(self.w, f, self.mu), expanded)

    def test_identity_for_additional_rational_backgrounds(self):
        for scale in (Q(-2, 3), Q(0), Q(7, 5)):
            w = [scale*p for p in self.h]
            mu = Q(17, 100)
            q = omega(w, mu)
            self.assertEqual(complete(w, [-p.parity() for p in q], mu), [0, 0])

    def test_division_does_not_invert_arbitrary_residual(self):
        p = 1+X+TAU*X*X
        self.assertEqual(X*p.div_x(), p-1)
        self.assertNotEqual(X*p.div_x(), p)


if __name__ == '__main__':
    unittest.main()
