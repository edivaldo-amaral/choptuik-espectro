#!/usr/bin/env python3
"""Revisao R1: teste proprio da formula (*) de GLOBAL_NULL_GAUGE_TRANSPORT.md com fonte ANALITICA NAO
POLINOMIAL, resolvida numericamente por quadratura e substituida na equacao por diferencas finitas.

  T_c(s) y = d_tau y + mu (xi - c) d_xi y + (s - mu) y = f,
  y = R_c(s) f_c + int_0^inf e^{-(s-mu)t} g(tau - t, c + e^{-mu t}(xi - c)) dt,  g = f - f_c, f_c = f(., c).

Confere tambem: periodicidade 4pi, reflexao y_+ = P y_- (c = -1 com f_+ = P f_-), paridades A par / B impar,
convergencia para 0 < Re s < mu, e a ressonancia s0 = mu - i m0/2 (solucao com fonte compativel e a
identidade secular T_1(s0)(tau e^{i m0 tau/2}) = e^{i m0 tau/2}).
Tambem conferido aqui: a derivacao de T_{+-1} a partir de d_{u_+} X^{u_-} = -dg_{++}/(2 g_{+-}) (sympy, se
disponivel no PYTHONPATH).
"""
import cmath
import math
import sys
import numpy as np
from scipy.integrate import quad

MU = 0.16830707896344996951


def f_minus(tau, xi):
    # analitica nao polinomial, 4pi-periodica, com polo em xi = 3 e xi = -2.5 (fora de [-1, 1])
    return (cmath.exp(0.7*xi)*(1 + 0.3*cmath.cos(tau/2)) + cmath.sin(tau)/(3 - xi)
            + 0.2j*cmath.exp(1j*tau/2)/(2.5 + xi))


def f_plus(tau, xi):
    return f_minus(tau, -xi)


def g_minus(tau, xi, c):
    """f_minus(tau, xi) - f_minus(tau, c), calculada sem cancelamento (d = xi - c pode ser ~1e-30)."""
    d = xi - c
    return (cmath.exp(0.7*c)*math.expm1(0.7*d)*(1 + 0.3*cmath.cos(tau/2))
            + cmath.sin(tau)*d/((3 - xi)*(3 - c))
            - 0.2j*cmath.exp(1j*tau/2)*d/((2.5 + xi)*(2.5 + c)))


def g_plus(tau, xi, c):
    return g_minus(tau, -xi, -c)


GDIFF = {}


def fourier_trace(f, c, mmax=40, npts=512):
    taus = np.arange(npts)*4*np.pi/npts
    vals = np.array([f(t, c) for t in taus])
    coef = np.fft.fft(vals)/npts           # coef[k] <-> e^{i k tau/2}
    ms = np.fft.fftfreq(npts, d=1.0/npts).astype(int)
    return {int(m): coef[k] for k, m in enumerate(ms) if abs(m) <= mmax}


def solve(f, c, s, skip_m=None, tmax=None):
    """y pela formula (*); skip_m: na ressonancia, inverte so no complemento do modo m0 (coef livre = 0)."""
    tr = fourier_trace(f, c)
    if tmax is None:
        tmax = 45.0/s.real

    def y(tau, xi):
        rc = 0j
        for m, a in tr.items():
            if skip_m is not None and m == skip_m:
                assert abs(a) < 1e-12, 'fonte incompativel na ressonancia'
                continue
            rc += a/(s - MU + 1j*m/2)*cmath.exp(1j*m*tau/2)

        def integrand(t, part):
            tt, xx = tau - t, c + math.exp(-MU*t)*(xi - c)
            gd = GDIFF.get(f)
            v = cmath.exp(-(s - MU)*t)*(gd(tt, xx, c) if gd else (f(tt, xx) - f(tt, c)))
            return v.real if part == 0 else v.imag
        # quadratura por blocos de comprimento 4pi (evita perda de precisao com decaimento lento)
        re = im = 0.0
        edges = np.arange(0.0, tmax + 4*math.pi, 4*math.pi)
        for a, b in zip(edges[:-1], edges[1:]):
            re += quad(integrand, a, b, args=(0,), limit=200, epsabs=1e-15, epsrel=1e-14)[0]
            im += quad(integrand, a, b, args=(1,), limit=200, epsabs=1e-15, epsrel=1e-14)[0]
        return rc + re + 1j*im
    return y


def residual(y, f, c, s, tau, xi, h=4e-3):
    # diferencas centrais de 4a ordem
    def d(fun, a, b, dir_):
        if dir_ == 0:
            return (-fun(a + 2*h, b) + 8*fun(a + h, b) - 8*fun(a - h, b) + fun(a - 2*h, b))/(12*h)
        return (-fun(a, b + 2*h) + 8*fun(a, b + h) - 8*fun(a, b - h) + fun(a, b - 2*h))/(12*h)
    lhs = d(y, tau, xi, 0) + MU*(xi - c)*d(y, tau, xi, 1) + (s - MU)*y(tau, xi)
    return abs(lhs - f(tau, xi)), abs(f(tau, xi))


def main():
    GDIFF[f_minus] = g_minus
    GDIFF[f_plus] = g_plus
    ok = True
    pts = [(0.3, -0.9), (1.7, 0.0), (5.0, 0.5), (11.0, 1.0), (2.2, 1.02)]
    for s in (0.09 + 0.3j, 0.7332 + 0.0j, 0.15 - 1.1j):
        ym = solve(f_minus, 1.0, s)
        yp = solve(f_plus, -1.0, s)
        print(f's = {s}:  (0 < Re s < mu)' if s.real < MU else f's = {s}:')
        for (t, x) in pts:
            r, n = residual(ym, f_minus, 1.0, s, t, x)
            per = abs(ym(t + 4*math.pi, x) - ym(t, x))
            refl = abs(yp(t, -x) - ym(t, x))
            print(f'  (tau, xi) = ({t}, {x}): |T_1 y - f| = {r:.2e} (|f| = {n:.2f}), '
                  f'|y(tau+4pi) - y| = {per:.1e}, |y_+(-xi) - y_-(xi)| = {refl:.1e}')
            ok &= r < 1e-7 and per < 1e-9 and refl < 1e-9
        # paridades de A = (y_+ + y_-)/(2 mu) e B = ((1+xi) y_- - (1-xi) y_+)/2
        t, x = 0.8, 0.37
        A = lambda xx: (yp(t, xx) + ym(t, xx))/(2*MU)
        B = lambda xx: ((1 + xx)*ym(t, xx) - (1 - xx)*yp(t, xx))/2
        pa, pb = abs(A(x) - A(-x)), abs(B(x) + B(-x))
        print(f'  A par: |A(xi) - A(-xi)| = {pa:.1e};  B impar: |B(xi) + B(-xi)| = {pb:.1e}')
        ok &= pa < 1e-9 and pb < 1e-9
    # ressonancia s0 = mu - i m0/2, fonte compativel: f - (traco no cone, modo m0) e^{i m0 tau/2}
    for m0 in (0, 1):
        s0 = MU - 1j*m0/2
        kappa = fourier_trace(f_minus, 1.0)[m0]
        fc = lambda tau, xi, k=kappa, m=m0: f_minus(tau, xi) - k*cmath.exp(1j*m*tau/2)
        GDIFF[fc] = g_minus
        y = solve(fc, 1.0, s0, skip_m=m0)
        r, n = residual(y, fc, 1.0, s0, 1.3, 0.4)
        sec = lambda tau, xi, m=m0: tau*cmath.exp(1j*m*tau/2)
        rs, _ = residual(sec, lambda tau, xi, m=m0: cmath.exp(1j*m*tau/2), 1.0, s0, 1.3, 0.4)
        print(f'ressonancia m0 = {m0}: kappa = {kappa:.4f}; |T_1(s0) y_reg - (f - kappa e)| = {r:.2e}; '
              f'|T_1(s0)(tau e) - e| = {rs:.1e}')
        ok &= r < 1e-7 and rs < 1e-9
    try:
        import sympy as sp
        tau, xi, mu, s = sp.symbols('tau xi mu s')
        up, um = -(1 + xi)*sp.exp(-mu*tau), -(1 - xi)*sp.exp(-mu*tau)
        Jm = sp.Matrix([[sp.diff(up, tau), sp.diff(up, xi)], [sp.diff(um, tau), sp.diff(um, xi)]])
        inv = Jm.inv()  # colunas: d/du_+ e d/du_- em (d_tau, d_xi)
        dplus = (inv[0, 0], inv[1, 0])
        y = sp.Function('y')(tau, xi)
        Xm = sp.exp((s - mu)*tau)*y
        lhs = sp.simplify(dplus[0]*sp.diff(Xm, tau) + dplus[1]*sp.diff(Xm, xi))
        T1 = sp.diff(y, tau) + mu*(xi - 1)*sp.diff(y, xi) + (s - mu)*y
        ratio = sp.simplify(lhs/(sp.exp(s*tau)*T1))
        print(f'sympy: d_{{u+}}(e^{{(s-mu)tau}} y) = {ratio} * e^{{s tau}} T_1(s) y  (esperado 1/(2 mu))')
        ok &= sp.simplify(ratio - 1/(2*mu)) == 0
    except ImportError:
        print('sympy ausente: derivacao simbolica pulada')
    print('R1 NUMERICO: CONFERE' if ok else 'R1 NUMERICO: FALHA')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
