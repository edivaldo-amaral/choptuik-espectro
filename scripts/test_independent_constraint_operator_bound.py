"""Testes rapidos (<5 s) de independent_constraint_operator_bound.py (Fisico, rodada 3)."""
from fractions import Fraction as Fr
import unittest

import numpy as np

import independent_constraint_operator_bound as M
import independent_numeric_audit as NA
import independent_rt_symbolic as S


def D_formula(h):
    h1, h2, h3, _ = h
    xdx = lambda f: S.X * f.dx()
    return [S.MU * (xdx(h1) + 2 * h1),
            -1 * S.X * (h2.dt() + S.S * h2 + S.MU * h2) + S.MU * (xdx(h2) - S.X * xdx(h2)) + 2 * S.MU * (h2 - h3)]


def A_formula(w, h):
    g = S.gamma2_tex(w, h)
    return [2 * S.odd(S.X * g[0]), -2 * (S.X * g[2]).P()]


class Decomposition(unittest.TestCase):
    def test_C_equals_D_plus_A_universally(self):
        def chk(order):
            src = S.SymbolSource()
            h, w = src.multifield(order), src.multifield(order)
            _, C = S.select(S.lin_omega_fc(w, h, S.S))
            return [c - d - a for c, d, a in zip(C, D_formula(h), A_formula(w, h))]
        for deg, order in ((1, 2), (2, 1)):
            S.CAP[0] = deg
            self.assertTrue(S.all_zero(chk(order)))

    def test_dropping_xi_squared_term_fails(self):
        src = S.SymbolSource()
        S.CAP[0] = 1
        h = src.multifield(2)
        w = [S.Poly() for _ in range(4)]
        _, C = S.select(S.lin_omega_fc(w, h, S.S))
        wrong = D_formula(h)
        wrong[1] = wrong[1] + S.MU * S.X * (S.X * h[1].dx())
        self.assertFalse(S.all_zero([c - d for c, d in zip(C, wrong)]))


class ChebyshevExact(unittest.TestCase):
    def test_xi_and_xidxi_match_float_engine(self):
        rng = np.random.default_rng(3)
        coeffs = {n: Fr(int(rng.integers(-9, 10)), int(rng.integers(1, 9))) for n in range(12)}
        f = NA.Fld.zeros(0, 11)
        for n, c in coeffs.items():
            f.a[0, n] = float(c)
        for exact, flt in ((M.xi_apply(coeffs), f.times_xi()), (M.xidxi_apply(coeffs), f.xi_dxi())):
            for n in range(flt.N + 1):
                self.assertAlmostEqual(float(exact.get(n, 0)), flt.a[0, n].real, places=12)

    def test_C1_norm_is_nine_eighths(self):
        e0 = {0: Fr(1)}
        self.assertEqual(M.nnorm(M.xi_apply(e0), M.B2) / M.nnorm(e0, M.B), Fr(9, 8))
        for n in range(1, 30):
            e = {n: Fr(1)}
            self.assertLessEqual(M.nnorm(M.xi_apply(e), M.B2) / M.nnorm(e, M.B), Fr(9, 8))

    def test_tail_bounds_dominate_and_decrease(self):
        for n in range(15, 60, 2):
            e = {n: Fr(1)}
            den = M.nnorm(e, M.B)
            T1, X, Y = M.tail_bounds(n)
            self.assertGreaterEqual(T1, M.nnorm(M.add((1, M.xidxi_apply(e)), (2, e)), M.B2) / den)
            xd = M.xidxi_apply(e)
            y = M.nnorm(M.add((1, xd), (-1, M.xi_apply(xd)), (2, e), (-1, M.xi_apply(e))), M.B2) / den
            self.assertGreaterEqual(Y, y)
            self.assertGreaterEqual(X, M.nnorm(M.xi_apply(e), M.B2) / den)
            T1n, Xn, Yn = M.tail_bounds(n + 1)
            self.assertLessEqual(T1n, T1)
            self.assertLessEqual(Yn, Y)

    def test_fourier_sup_monotone_beyond_130(self):
        v130 = Fr(130, 2) * M.Q_A ** 130
        self.assertGreaterEqual(M.fourier_sup(), v130)
        self.assertLess(Fr(131, 130) * M.Q_A, 1)


if __name__ == '__main__':
    unittest.main()
