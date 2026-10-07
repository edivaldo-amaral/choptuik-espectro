#!/usr/bin/env python3
"""Revisao R4: confere, de forma independente, as formulas de curvatura de RT (eq. dfkhdjhsdshkfd)
e as cinco identidades de definicao (Omega1 - Omega2 - Omega2sharp = 0, Omega2/3/4 iguais nos dois sinais),
para (zeta, Q, phi) GENERICOS (funcoes simbolicas de (tau, xi)), sem supor as equacoes de campo.

Metodo: metrica 4D em (tau, xi, theta, varphi) construida a partir do frame de RT (ekhfkjfhkfhfkhf):
  g_2D = e^{-2 zeta}[ -(1 - xi^2) dtau^2 - (2 xi/mu) dtau dxi + dxi^2/mu^2 ],  r = e^{-zeta} xi / W,  W = mu + xi^2 Q,
  e0 = d_tau + mu xi d_xi,  e_r = mu d_xi,  e_theta = (W/xi) d_theta.
Ricci calculado diretamente pelos simbolos de Christoffel (sem formulas de produto torcido). Os dois lados sao
avaliados em jatos racionais aleatorios (zeta = 0 no ponto, pela invariancia conforme constante), de forma exata.
Requer sympy (PYTHONPATH do scratchpad)."""
import random
import sys
import sympy as sp

t, x, th, ph = sp.symbols('tau xi theta varphi', real=True)
mu = sp.Rational(1, 6) + sp.Rational(1, 997)   # mu generico racional
Z = sp.Function('zeta')(t, x)
Qf = sp.Function('Q')(t, x)
F = sp.Function('phi')(t, x)
W = mu + x**2*Qf
e2 = sp.exp(-2*Z)
r = sp.exp(-Z)*x/W

coords = [t, x, th, ph]
g = sp.zeros(4, 4)
g[0, 0] = -e2*(1 - x**2)
g[0, 1] = g[1, 0] = -e2*x/mu
g[1, 1] = e2/mu**2
g[2, 2] = r**2
g[3, 3] = r**2*sp.sin(th)**2
ginv = sp.zeros(4, 4)
g2inv = g[:2, :2].inv()
for i in range(2):
    for j in range(2):
        ginv[i, j] = g2inv[i, j]
ginv[2, 2] = 1/g[2, 2]
ginv[3, 3] = 1/g[3, 3]

Gam = [[[sum(ginv[a, d]*(sp.diff(g[d, b], coords[c]) + sp.diff(g[d, c], coords[b]) - sp.diff(g[b, c], coords[d]))
             for d in range(4))/2 for c in range(4)] for b in range(4)] for a in range(4)]


def ricci(b, c):
    s = 0
    for a in range(4):
        s += sp.diff(Gam[a][b][c], coords[a]) - sp.diff(Gam[a][b][a], coords[c])
        for d in range(4):
            s += Gam[a][a][d]*Gam[d][b][c] - Gam[a][c][d]*Gam[d][b][a]
    return s


Ric = {(b, c): ricci(b, c) for b in range(3) for c in range(b, 3)}


def R(b, c):
    return Ric[(min(b, c), max(b, c))]


# frame (componentes em (tau, xi, theta))
e0 = [1, mu*x, 0]
er = [0, mu, 0]
eth = [0, 0, W/x]


def ric_frame(u, v):
    return sum(u[i]*v[j]*R(i, j) for i in range(3) for j in range(3))


def dvec(u, f):
    return u[0]*sp.diff(f, t) + u[1]*sp.diff(f, x) + u[2]*sp.diff(f, th)


# box phi = (1/sqrt|g|) d_a (sqrt|g| g^{ab} d_b phi); sqrt|g| = sqrt(-det g2) r^2 sin(theta)
sq = sp.sqrt(-g[:2, :2].det())*r**2*sp.sin(th)
box = sum(sp.diff(sq*ginv[a, b]*sp.diff(F, coords[b]), coords[a]) for a in range(2) for b in range(2))/sq

LHS = {
    'theta': ric_frame(eth, eth) - 2*dvec(eth, F)**2,
    'rr': ric_frame(er, er) - 2*dvec(er, F)**2,
    '00': ric_frame(e0, e0) - 2*dvec(e0, F)**2,
    '0r': ric_frame(e0, er) - 2*dvec(e0, F)*dvec(er, F),
    'box': e2*box,
}


# lado de RT
def D(sig, f):
    return sig*sp.diff(f, t) + mu*(1 + sig*x)*sp.diff(f, x)


L = Z + sp.log(W)
w1 = 2*x*Qf + D(-1, L) + D(1, L)
w2 = {s: D(s, L) - s*mu for s in (-1, 1)}
w3 = {s: D(s, Z) - s*mu for s in (-1, 1)}
w4 = {s: D(s, F) for s in (-1, 1)}


def Om(s):
    m = -s
    O1 = x*(D(s, w1) + s*mu*w1) + mu*2*w1 + x/4*(w1*w1 - 2*w2[m]*w1 - 6*w2[s]*w1 + 4*w3[s]*w1 + w2[m]**2
                                                  - 2*w2[s]*w2[m] - 4*w3[s]*w2[m] + w2[s]**2 + 4*w3[s]*w2[s]
                                                  - 4*w4[s]**2)
    O2 = x*(D(s, w2[m]) + s*mu*w2[m]) + mu*(w1 + w2[m] + w2[s]) + x/4*(w1*w1 - 2*w2[-1]*w1 - 2*w2[1]*w1
                                                                      + w2[-1]**2 - 6*w2[1]*w2[-1] + w2[1]**2)
    O2s = x*(D(s, w2[s]) + s*mu*w2[s]) + mu*(2*w2[s] - 2*w3[s]) + x/4*(-4*w2[s]**2 + 8*w3[s]*w2[s] - 4*w4[s]**2)
    O3 = x*(D(s, w3[m]) + s*mu*w3[m]) + mu*(-w1) + x/4*(-w1*w1 + 2*w2[-1]*w1 + 2*w2[1]*w1 - w2[-1]**2
                                                       + 2*w2[1]*w2[-1] - w2[1]**2 - 4*w4[1]*w4[-1])
    O4 = x*(D(s, w4[m]) + s*mu*w4[m]) + mu*(w4[m] + w4[s]) + x/4*(-4*w4[1]*w2[-1] - 4*w4[-1]*w2[1])
    return dict(O1=O1, O2=O2, O2s=O2s, O3=O3, O4=O4)


Om_m, Om_p = Om(-1), Om(1)
A = Om_m['O1'] + Om_p['O1'] + Om_m['O2'] + Om_p['O2'] - Om_m['O2s'] - Om_p['O2s']
B = (Om_m['O1'] + Om_p['O1'] - Om_m['O2'] - Om_p['O2'] + Om_m['O2s'] + Om_p['O2s']
     + 2*Om_m['O3'] + 2*Om_p['O3'])
RHS = {
    'theta': A/(4*x),
    'rr': (A + B)/(4*x),
    '00': (-Om_m['O2'] - Om_p['O2'] + Om_m['O2s'] + Om_p['O2s'] - Om_m['O3'] - Om_p['O3'])/(2*x),
    '0r': (-Om_m['O1'] + Om_p['O1'] - Om_m['O3'] + Om_p['O3'])/(2*x),
    'box': (Om_m['O4'] + Om_p['O4'])/(2*x),
}
IDENT = {
    'O1-O2-O2s (-)': Om_m['O1'] - Om_m['O2'] - Om_m['O2s'],
    'O1-O2-O2s (+)': Om_p['O1'] - Om_p['O2'] - Om_p['O2s'],
    'O2- - O2+': Om_m['O2'] - Om_p['O2'],
    'O3- - O3+': Om_m['O3'] - Om_p['O3'],
    'O4- - O4+': Om_m['O4'] - Om_p['O4'],
}


def jet_subs(seed):
    rnd = random.Random(seed)
    subs = {}
    for fn in (Z, Qf, F):
        for a in range(4):
            for b in range(4 - a):
                if a == b == 0:
                    val = sp.Integer(0) if fn is Z else sp.Rational(rnd.randint(-9, 9), rnd.randint(1, 7))
                    subs[fn] = val
                else:
                    subs[sp.Derivative(fn, *([t]*a + [x]*b))] = sp.Rational(rnd.randint(-9, 9), rnd.randint(1, 7))
    pt = {t: sp.Rational(rnd.randint(-5, 5), 3), x: sp.Rational(rnd.randint(1, 9), 10), th: sp.pi/3}
    return subs, pt


def evaluate(expr, subs, pt):
    # derivadas de ordem alta primeiro
    keys = sorted(subs, key=lambda k: -(len(k.variables) if isinstance(k, sp.Derivative) else 0))
    e = expr.doit()
    for k in keys:
        e = e.subs(k, subs[k])
    e = e.subs(pt)
    return sp.nsimplify(sp.simplify(e))


def main():
    ok = True
    for seed in (1, 2, 3):
        subs, pt = jet_subs(seed)
        print(f'--- jato aleatorio {seed}: xi = {pt[x]}, tau = {pt[t]}')
        for k in LHS:
            d = evaluate(LHS[k] - RHS[k], subs, pt)
            lv = evaluate(LHS[k], subs, pt)
            print(f'  {k:6s}: LHS = {sp.N(lv, 12)},  LHS - RHS = {d}')
            ok &= (d == 0)
        for k, v in IDENT.items():
            d = evaluate(v, subs, pt)
            print(f'  identidade {k}: {d}')
            ok &= (d == 0)
    print('TUDO CONFERE' if ok else 'DIVERGENCIA')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
