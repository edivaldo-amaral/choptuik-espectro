#!/usr/bin/env python3
"""G5: confere scripts/gauge_null_profiles.null_generator contra o pushforward calculado aqui com sympy
diretamente da troca de coordenadas u_sigma = -(1 + sigma xi) e^{-mu tau} (sem usar as formulas A_n, B_n).

Para f(u) = u^n, X = f(u_-) d_- + f(u_+) d_+; em (tau, xi), X = e^{s tau}(A(xi) d_tau + B(xi) d_xi) com
s = mu(1 - n).  Tambem confere o multiplicador de Floquet de (Theta^2)^* X: e^{4 pi s}."""
import sys
from fractions import Fraction

import sympy as sp

sys.path.insert(0, "<repo>/scripts")
from gauge_null_profiles import null_generator  # noqa: E402

tau, xi = sp.symbols("tau xi", real=True)
um, up = sp.symbols("u_m u_p", real=True)
falhas = 0
for mu_q in (Fraction(1, 6), Fraction(17, 100), Fraction(90359175, 536870912)):
    mu = sp.Rational(mu_q.numerator, mu_q.denominator)
    U_p = -(1 + xi) * sp.exp(-mu * tau)
    U_m = -(1 - xi) * sp.exp(-mu * tau)
    TAU_u = -sp.log(-(up + um) / 2) / mu
    XI_u = (up - um) / (up + um)
    for n in range(0, 7):
        Xt = (um**n * sp.diff(TAU_u, um) + up**n * sp.diff(TAU_u, up)).subs({um: U_m, up: U_p})
        Xx = (um**n * sp.diff(XI_u, um) + up**n * sp.diff(XI_u, up)).subs({um: U_m, up: U_p})
        s_n, Ptau, Pxi = null_generator(n, mu_q)
        # Poly -> sympy (variaveis (tau, xi, s); aqui so xi aparece)
        conv = lambda P: sum(sp.Rational(v.numerator, v.denominator) * tau**a * xi**b for (a, b, c), v in P.terms.items())
        A = conv(Ptau)
        B = conv(Pxi)
        s = sp.Rational(s_n.numerator, s_n.denominator)
        r1 = sp.simplify(sp.expand(Xt - sp.exp(s * tau) * A))
        r2 = sp.simplify(sp.expand(Xx - sp.exp(s * tau) * B))
        r3 = sp.simplify(s - mu * (1 - n))
        # multiplicador: X(tau + 4 pi) = e^{4 pi s} X
        r4 = sp.simplify((Xt.subs(tau, tau + 4 * sp.pi) - sp.exp(4 * sp.pi * s) * Xt))
        ok = (r1 == 0 and r2 == 0 and r3 == 0 and r4 == 0)
        falhas += (not ok)
        print(f"mu={mu_q} n={n}: s_n={s_n} (= mu(1-n): {r3 == 0}); tau ok={r1 == 0}; xi ok={r2 == 0}; "
              f"multiplicador e^(4 pi s) ok={r4 == 0}")
print("FALHAS:", falhas)
