#!/usr/bin/env python3
"""Revisao V (03/10/2026), V5/V6: confere as formulas de curvatura (dfkhdjhsdshkfd) de RT com o
Ricci 4D calculado DIRETO da metrica, com zeta, Q, phi funcoes GENERICAS de (tau, xi).

Metrica do ansatz (deduzida aqui do frame de RT, ekhfkjfhkfhfkhf):
    g^{-1} = e^{2 zeta}(-e0 e0 + sum_i e_i e_i),  e0 = d_tau + mu xi d_xi,
    radial: x^k e_k/xi = mu d_xi,  tangencial unitario t: t.e = W t.d,  W = mu + xi^2 Q,
    => g = e^{-2 zeta}[-(1-xi^2) dtau^2 - (2 xi/mu) dtau dxi + dxi^2/mu^2] + e^{-2 zeta} xi^2 W^{-2} g_S2.
Frame em (tau, xi, theta, varphi): e0 = (1, mu xi, 0, 0), e_r = (0, mu, 0, 0), e_th = (0, 0, W/xi, 0).

RT: Ric(e_i,e_j) - 2 e_i(phi) e_j(phi) = delta_ij A/(4xi) + x^i x^j B/(4 xi^3), logo
    Ric(e_r,e_r) - 2 e_r(phi)^2 = (A+B)/(4xi),  Ric(e_th,e_th) = A/(4xi),
    Ric(e0,e0) - 2 e0(phi)^2 = (-O2- -O2+ +O2#- +O2#+ -O3- -O3+)/(2xi),
    Ric(e0,e_r) - 2 e0(phi) e_r(phi) = (-O1- +O1+ -O3- +O3+)/(2xi),
    e^{-2 zeta} Box phi = (O4- + O4+)/(2xi).
Os omega vem de (ekjehee), os Omega de (eezuuirzr) (reuso de v123_identidades.Omegas).
Tambem: (a) o frame e ortonormal a menos de e^{-2zeta}; (b) as cinco identidades de definicao
(Omega1-Omega2-Omega2# = 0, Omega_j^- = Omega_j^+, j = 2,3,4) valem para omega de (ekjehee);
(c) controle negativo: alterar um coeficiente de Omega faz a comparacao falhar.
"""
import sympy as sp

import pathlib

# Importar v123 rodaria todas as verificacoes; copiamos so o necessario.
src = pathlib.Path(__file__).with_name('v123_identidades.py').read_text()
ns = {}
exec(src.split('resultados = {}')[0], ns)
Omegas, D = ns['Omegas'], ns['D']
tau, xi, mu = ns['tau'], ns['xi'], ns['mu']

th, ph = sp.symbols('theta varphi')
X = [tau, xi, th, ph]
zeta = sp.Function('zeta')(tau, xi)
Qf = sp.Function('Q')(tau, xi)
phi = sp.Function('phi')(tau, xi)
W = mu + xi**2*Qf
E = sp.exp(-2*zeta)

g = sp.zeros(4, 4)
g[0, 0] = -E*(1 - xi**2)
g[0, 1] = g[1, 0] = -E*xi/mu
g[1, 1] = E/mu**2
g[2, 2] = E*xi**2/W**2
g[3, 3] = E*xi**2/W**2*sp.sin(th)**2
ginv = sp.simplify(g.inv())

Gam = [[[sp.simplify(sum(ginv[a, d]*(sp.diff(g[d, b], X[c]) + sp.diff(g[d, c], X[b]) - sp.diff(g[b, c], X[d]))
                         for d in range(4))/2) for c in range(4)] for b in range(4)] for a in range(4)]


def ricci(b, d):
    r = 0
    for a in range(4):
        r += sp.diff(Gam[a][b][d], X[a]) - sp.diff(Gam[a][b][a], X[d])
        for e in range(4):
            r += Gam[a][a][e]*Gam[e][b][d] - Gam[a][d][e]*Gam[e][b][a]
    return r


Ric = sp.Matrix(4, 4, lambda b, d: ricci(b, d) if b <= d else 0)
for b in range(4):
    for d in range(b):
        Ric[b, d] = Ric[d, b]

e0 = sp.Matrix([1, mu*xi, 0, 0])
er = sp.Matrix([0, mu, 0, 0])
et = sp.Matrix([0, 0, W/xi, 0])
dphi = sp.Matrix([sp.diff(phi, x) for x in X])


def T(u, v, M):
    return (u.T*M*v)[0, 0]


def EE(u, v):
    return T(u, v, Ric) - 2*(u.T*dphi)[0, 0]*(v.T*dphi)[0, 0]


sqrtg = sp.sqrt(-g.det())
box = sum(sp.diff(sqrtg*sum(ginv[a, b]*sp.diff(phi, X[b]) for b in range(4)), X[a]) for a in range(4))/sqrtg

# omega de (ekjehee) com D^sigma de RT
L = zeta + sp.log(W)
w = {'w1': 2*xi*Qf + D(-1, L) + D(1, L)}
for sig in (1, -1):
    w[('w2', sig)] = D(sig, L) - sig*mu
    w[('w3', sig)] = D(sig, zeta) - sig*mu
    w[('w4', sig)] = D(sig, phi)
Om = {sig: Omegas(w, sig, lambda S, f: D(S, f)) for sig in (1, -1)}
A = Om[-1]['1'] + Om[1]['1'] + Om[-1]['2'] + Om[1]['2'] - Om[-1]['2s'] - Om[1]['2s']
B = Om[-1]['1'] + Om[1]['1'] - Om[-1]['2'] - Om[1]['2'] + Om[-1]['2s'] + Om[1]['2s'] + 2*Om[-1]['3'] + 2*Om[1]['3']

# substituicao de derivadas por simbolos, para simplificar como funcoes racionais
subs = []
for f, nm in ((zeta, 'z'), (Qf, 'q'), (phi, 'p')):
    for i in range(3, -1, -1):
        for j in range(3 - i, -1, -1):
            if i + j == 0:
                continue
            der = sp.Derivative(f, *([tau]*i + [xi]*j))
            subs.append((der, sp.Symbol(f'{nm}_{i}{j}')))
for f, nm in ((zeta, 'z'), (Qf, 'q'), (phi, 'p')):
    subs.append((f, sp.Symbol(f'{nm}_00')))


def red(expr):
    expr = expr.doit()
    for a, b in subs:
        expr = expr.subs(a, b)
    return sp.cancel(sp.together(sp.expand(expr)))


checks = {
    'frame g(e0,e0) = -e^{-2z}': red(T(e0, e0, g) + E),
    'frame g(er,er) = e^{-2z}': red(T(er, er, g) - E),
    'frame g(et,et) = e^{-2z}': red(T(et, et, g) - E),
    'frame g(e0,er) = 0': red(T(e0, er, g)),
    'Ric(er,er)-2(er phi)^2 = (A+B)/(4xi)': red(EE(er, er) - (A + B)/(4*xi)),
    'Ric(et,et) = A/(4xi)': red(EE(et, et) - A/(4*xi)),
    'Ric(e0,e0)-2(e0 phi)^2': red(EE(e0, e0) - (-Om[-1]['2'] - Om[1]['2'] + Om[-1]['2s'] + Om[1]['2s']
                                                - Om[-1]['3'] - Om[1]['3'])/(2*xi)),
    'Ric(e0,er)-2 e0phi er phi': red(EE(e0, er) - (-Om[-1]['1'] + Om[1]['1'] - Om[-1]['3'] + Om[1]['3'])/(2*xi)),
    'Ric(e0,et) = 0': red(EE(e0, et)),
    'Ric(er,et) = 0': red(EE(er, et)),
    'e^{-2z} Box phi = (O4- + O4+)/(2xi)': red(E*box - (Om[-1]['4'] + Om[1]['4'])/(2*xi)),
    'def: Omega1^+ - Omega2^+ - Omega2#^+ = 0': red(Om[1]['1'] - Om[1]['2'] - Om[1]['2s']),
    'def: Omega1^- - Omega2^- - Omega2#^- = 0': red(Om[-1]['1'] - Om[-1]['2'] - Om[-1]['2s']),
    'def: Omega2^- = Omega2^+': red(Om[-1]['2'] - Om[1]['2']),
    'def: Omega3^- = Omega3^+': red(Om[-1]['3'] - Om[1]['3']),
    'def: Omega4^- = Omega4^+': red(Om[-1]['4'] - Om[1]['4']),
}
ok = True
for k, v in checks.items():
    good = (v == 0)
    ok &= good
    print(f"{'OK ' if good else 'FALHA'} {k}: {'0' if good else v}")

# controle negativo: mudar o coeficiente 6 -> 5 em Omega1 (termo omega2^sigma omega1)
cn = red(EE(er, er) - (A + B + sp.Rational(1, 4)*xi*(w[('w2', 1)]*w['w1'] + w[('w2', -1)]*w['w1']))/(4*xi))
print(f"{'OK ' if cn != 0 else 'FALHA'} controle negativo (Omega1 alterado deve falhar): residuo {'!= 0' if cn != 0 else '= 0'}")
ok &= (cn != 0)


print('TODAS AS VERIFICACOES CONFEREM' if ok else 'HA DISCREPANCIAS')
