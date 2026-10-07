from fractions import Fraction as Q
import unittest
from sharp_propagation import SOURCE, Poly, TAU, X, omega, linearized_omega
from geometric_reconstruction import (W_of, null_derivative, temporal_form,
    radial_form, reconstruct_polynomial, inverse_tau_polynomial,
    reconstruct_polynomial_chain)


@unittest.skipUnless((SOURCE/'GAMMA2_macro.c').exists(), 'fonte RT ausente')
class GeometricReconstructionTests(unittest.TestCase):
    def test_W_residual_identity_and_linearization_offshell(self):
        mu, s = Q(1, 6), Q(3, 4)
        w = [X*(1+TAU), TAU*X+2, 1+X*X, TAU+X]
        h = [X*(2+TAU*TAU), X+TAU, 2+TAU*X, 1+X*X]
        W = W_of(w, mu)
        jump = -w[1].parity()+w[2].parity()
        residual = null_derivative(W, mu, 1)-W*jump
        om = omega(w, mu)
        self.assertEqual(residual, (om[0]-om[1]-om[2])/2)
        dw = W_of(h, Q(0))
        djump = -h[1].parity()+h[2].parity()
        dr = null_derivative(dw, mu, 1, s)-dw*jump-W*djump
        dom = linearized_omega(w, h, mu, s)
        self.assertEqual(dr, (dom[0]-dom[1]-dom[2])/2)

    def test_closed_form_residuals_without_on_shell_assumption(self):
        mu, s = Q(1, 6), Q(3, 4)
        w = [X*(1+TAU), TAU*X+2, 1+X*X, TAU+X]
        h = [X*(2+TAU*TAU), X+TAU, 2+TAU*X, 1+X*X]
        dom = linearized_omega(w, h, mu, s)
        for component, index in ((1, 1), (2, 3), (3, 4)):
            a, b = temporal_form(h[component]), radial_form(h[component], mu)
            closure = b.derivative(0)+s*b-a.derivative(1)
            self.assertEqual(-2*mu*X*closure, -dom[index].parity()-dom[index])

    def test_reconstruct_definitions_at_constant_W(self):
        mu, s = Q(1, 6), Q(3, 4)
        q = 1+TAU*X*X
        z = 2+TAU*TAU+X*X
        p = TAU+3*X*X
        ell = X*X*q/mu
        h = [2*X*q+2*mu*(z+ell).derivative(1),
             null_derivative(z+ell, mu, -1, s),
             null_derivative(z, mu, -1, s),
             null_derivative(p, mu, -1, s)]
        self.assertEqual(reconstruct_polynomial(h, mu, s), (q, z, p))
        for scalar in (z, p):
            self.assertEqual(inverse_tau_polynomial(scalar.derivative(0)+s*scalar, s), scalar)

    def test_reconstruction_of_jets_includes_previous_potential(self):
        mu, s = Q(1, 6), Q(3, 4)
        profiles = [(1+X*X, 2+TAU+X*X, 3+TAU*TAU),
                    (TAU+X*X, 1+TAU*TAU, 2+TAU*X*X)]
        chain, previous = [], (Poly(), Poly(), Poly())
        for q, z, p in profiles:
            old_q, old_z, old_p = previous
            f, old_f = z+X*X*q/mu, old_z+X*X*old_q/mu
            chain.append([2*X*q+2*mu*f.derivative(1),
                          null_derivative(f, mu, -1, s)-old_f,
                          null_derivative(z, mu, -1, s)-old_z,
                          null_derivative(p, mu, -1, s)-old_p])
            previous = q, z, p
        self.assertEqual(reconstruct_polynomial_chain(chain, mu, s), profiles)
        self.assertNotEqual(reconstruct_polynomial(chain[1], mu, s)[1], profiles[1][1])

    def test_geometric_residual_combinations_are_invertible(self):
        from build_tile import exact_inverse
        # Colunas: Omega1-, Omega1+, Omega2, Omega3, Omega4.
        matrix = [[Q(x) for x in row] for row in
                  [[0,0,4,0,0], [2,2,-4,4,0], [1,1,-4,-2,0],
                   [-1,1,0,0,0], [0,0,0,0,2]]]
        inverse = exact_inverse(matrix)
        product = [[sum(inverse[i][k]*matrix[k][j] for k in range(5))
                    for j in range(5)] for i in range(5)]
        self.assertEqual(product, [[Q(i == j) for j in range(5)] for i in range(5)])


if __name__ == '__main__':
    unittest.main()
