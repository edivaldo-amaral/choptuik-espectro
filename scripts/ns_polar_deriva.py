#!/usr/bin/env python3
"""Perturbacoes polares (pares), l >= 2, campo escalar: derivacao simbolica (sympy) das equacoes de Gerlach-Sengupta
(Martin-Garcia & Gundlach 1999, eqs. inveqAB, inveqAeven, inveqwave, inveq3; tracos para T^2 = 0) no fundo de RT,
coordenadas (tau, xi). Gera os coeficientes das equacoes como funcoes do fundo para scripts/ns_polar.py.

Fundo: g_AB = e^{-2 zeta} [[-(1 - xi^2), -xi/mu], [-xi/mu, 1/mu^2]], r = e^{-zeta} xi/W, phi (normalizacao de RT,
Ric = 2 dphi dphi; nas formulas de MG, com 8 pi, entra phi/(2 sqrt(pi)) e Phi/(2 sqrt(pi))).
Perturbacao: k_AB = e^{-2 zeta} (A th+ th+ + B th- th-), th^s = mu(1 + s xi) dtau - s dxi (nulos), k, Phi.
Saida: build/naoesf/polar_coef.py com, para cada equacao, os coeficientes de cada derivada de cada incognita."""
from __future__ import annotations

import os
import sys
from pathlib import Path

sys.path.insert(0, os.environ.get('SYMPY_PATH', ''))
import sympy as sp

t, x, mu, l = sp.symbols('tau xi mu l', real=True)
zeta = sp.Function('zeta')(t, x); W = sp.Function('W')(t, x); ph = sp.Function('phi')(t, x)
X = (t, x)
Af, Bf, kf, Pf = [sp.Function(n)(t, x) for n in ('A', 'B', 'k', 'Phi')]
fourpi = sp.Integer(4)*sp.pi
s2pi = 2*sp.sqrt(sp.pi)
phM = ph/s2pi                              # phi de MG
PhM = Pf/s2pi                              # Phi de MG

ef = sp.exp(-2*zeta)
g = ef*sp.Matrix([[-(1 - x**2), -x/mu], [-x/mu, 1/mu**2]])
gi = sp.simplify(g.inv())
r = sp.exp(-zeta)*x/W
th = {+1: sp.Matrix([mu*(1 + x), -1]), -1: sp.Matrix([mu*(1 - x), 1])}
kAB = ef*(Af*th[1]*th[1].T + Bf*th[-1]*th[-1].T)

Gam = [[[sum(gi[c, d]*(sp.diff(g[d, a], X[b]) + sp.diff(g[d, b], X[a]) - sp.diff(g[a, b], X[d])) for d in range(2))/2
         for b in range(2)] for a in range(2)] for c in range(2)]


def d1(f):
    return [sp.diff(f, X[a]) for a in range(2)]


def hess(f):                                   # f_{|AB}
    return sp.Matrix(2, 2, lambda a, b: sp.diff(f, X[a], X[b]) - sum(Gam[c][a][b]*sp.diff(f, X[c]) for c in range(2)))


def cov_vec(V):                                # V_{A|B} (covetor)
    return sp.Matrix(2, 2, lambda a, b: sp.diff(V[a], X[b]) - sum(Gam[c][a][b]*V[c] for c in range(2)))


def cov_t2(T):                                 # T_{AB|C}
    return [[[sp.diff(T[a, b], X[c]) - sum(Gam[d][c][a]*T[d, b] + Gam[d][c][b]*T[a, d] for d in range(2))
              for c in range(2)] for b in range(2)] for a in range(2)]


def raise1(V):
    return [sum(gi[a, b]*V[b] for b in range(2)) for a in range(2)]


def tr(M):
    return sum(gi[a, b]*M[a, b] for a in range(2) for b in range(2))


v = [sp.diff(sp.log(r), X[a]) for a in range(2)]
vu = raise1(v)
vv = sum(v[a]*vu[a] for a in range(2))
V0 = r**-2 - vv                                # campo escalar (MG): V0 = r^-2 - v_A v^A
dk = d1(kf); dkU = raise1(dk)
dph = d1(phM); dphU = raise1(dph)
dPh = d1(PhM); dPhU = raise1(dPh)
kU = gi*kAB*gi                                 # k^{AB}
DkAB = cov_t2(kAB)                             # k_{AB|C} = DkAB[a][b][c]
divk = [sum(gi[b, c]*DkAB[a][b][c] for b in range(2) for c in range(2)) for a in range(2)]   # k_{AB}^{|B}
L = l*(l + 1)

# (inveqAeven), eta = 0, k^B_B = 0:  1/2 (k_AB^{|B} - k_{|A}) = 8 pi Phi phi_{|A}
EA = [sp.Rational(1, 2)*(divk[a] - dk[a]) - 8*sp.pi*PhM*dph[a] for a in range(2)]
# componente nula: contracao com D^- (vetor) -> transporta B? (escolhemos as duas e decidimos numericamente)
Dvec = {+1: sp.Matrix([1, mu*(1 + x)]), -1: sp.Matrix([-1, mu*(1 - x)])}
EAn = {s: sum(Dvec[s][a]*EA[a] for a in range(2)) for s in (1, -1)}

# (inveqAB) com eta = 0, k^D_D = 0; T_AB (invmatterAB)
Hk = hess(kf)
vD = cov_vec(v)                                 # v_{D|F}
vDF = sum(kU[a, b]*vD[a, b] for a in range(2) for b in range(2))   # k_DF v^{D|F}  (= k^{DF} v_{D|F})
kvv = sum(kAB[a, b]*vu[a]*vu[b] for a in range(2) for b in range(2))
lapk = tr(Hk)
divkU_v = sum(divk[a]*vu[a] for a in range(2))  # k_{CD}^{|D} v^C
term1 = sp.Matrix(2, 2, lambda a, b: sum((DkAB[c][a][b] + DkAB[c][b][a] - DkAB[a][b][c])*vu[c] for c in range(2)))
EAB_L = (term1 - g*(2*divkU_v)
         - sp.Matrix(2, 2, lambda a, b: dk[a]*v[b] + dk[b]*v[a] + Hk[a, b])
         + (V0 + L/(2*r**2))*kAB
         - g*(2*vDF + 3*kvv - lapk - 3*sum(dk[a]*vu[a] for a in range(2)) + (l - 1)*(l + 2)/(2*r**2)*kf))
phF_PhF = sum(dph[a]*dPhU[a] for a in range(2))
kphph = sum(kU[a, b]*dph[a]*dph[b] for a in range(2) for b in range(2))
phph = sum(dph[a]*dphU[a] for a in range(2))
TAB = sp.Matrix(2, 2, lambda a, b: dPh[a]*dph[b] + dPh[b]*dph[a] - g[a, b]*phF_PhF + g[a, b]*kphph/2 - kAB[a, b]*phph/2)
EAB = EAB_L - 8*sp.pi*TAB
Etr = tr(EAB)

# (inveqwave), eta = 0
lapPh = tr(hess(PhM))
rhs = sum(sp.diff(r**2*sum(kU[a, b]*dph[a] for a in range(2)), X[b]) for b in range(2))   # nao covariante! ver abaixo
# divergencia covariante de V^B = r^2 k^{AB} phi_A:  (1/sqrt|g|) d_B (sqrt|g| V^B)
sqg = sp.sqrt(-g.det())
VB = [r**2*sum(kU[a, b]*dph[a] for a in range(2)) for b in range(2)]
divV = sum(sp.diff(sqg*VB[b], X[b]) for b in range(2))/sqg
Ewave = (lapPh + 2*sum(vu[a]*dPh[a] for a in range(2)) - L*PhM/r**2
         - (divV/r**2 - sum(dk[b]*dphU[b] for b in range(2))))

EQ = {'wave': Ewave, 'trace': Etr, 'EAp': EAn[1], 'EAm': EAn[-1],
      'cpp': sum(Dvec[1][a]*Dvec[1][b]*EAB[a, b] for a in range(2) for b in range(2)),
      'cmm': sum(Dvec[-1][a]*Dvec[-1][b]*EAB[a, b] for a in range(2) for b in range(2))}

# homogeneidade em e^zeta: E(zeta + c) = e^{n c} E(zeta)
c = sp.symbols('c')
inc = {'A': Af, 'B': Bf, 'k': kf, 'Phi': Pf}
ders = [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2)]
out = ['# gerado por scripts/ns_polar_deriva.py; nao editar', 'import numpy as np', '', 'COEF = {']
bg = {}
names = {}
for nome, E in EQ.items():
    E = sp.expand(E)
    shifted = E.subs(zeta, zeta + c).doit()
    expo = None
    for n_try in (2, 0, 4, -2):
        if sp.simplify(sp.expand(shifted*sp.exp(-n_try*c) - E)) == 0:
            expo = n_try; break
    print(nome, 'homogeneo com expoente', expo, flush=True)
    assert expo is not None
    En = sp.expand(E*sp.exp(2*expo*0))
    out.append(f'  {nome!r}: {{')
    for un, f in inc.items():
        for (i, j) in ders:
            dd = f
            if i: dd = sp.diff(dd, t, i)
            if j: dd = sp.diff(dd, x, j)
            coef = sp.expand(En).coeff(dd)
            # remove dependencia de derivadas maiores (coeff pega so o termo linear exato)
            if coef == 0:
                continue
            # troca as derivadas do fundo por simbolos; DEPOIS zera zeta nao derivado (homogeneo)
            reps = {}
            for fn, nm in ((zeta, 'z'), (W, 'W'), (ph, 'p')):  # inclui zeta(t,x) -> 0 so depois
                for (a, b) in [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2), (3, 0), (2, 1), (1, 2), (0, 3)]:
                    if (a, b) == (0, 0):
                        if fn is zeta: continue
                        reps[fn] = sp.Symbol(nm)
                        continue
                    dfn = fn
                    if a: dfn = sp.diff(dfn, t, a)
                    if b: dfn = sp.diff(dfn, x, b)
                    reps[dfn] = sp.Symbol(f'{nm}_{"t"*a}{"x"*b}')
            # substitui das derivadas mais altas para as mais baixas
            coef = coef.xreplace(reps)            # simultaneo e estrutural: d2 zeta nao vira d(z_t) = 0
            coef = coef.subs(zeta, 0)
            assert not coef.has(zeta), coef
            coef = sp.simplify(coef)
            out.append(f'    ({un!r}, {i}, {j}): "{sp.sstr(coef)}",')
    out.append('  },')
out.append('}')
Path('build/naoesf').mkdir(parents=True, exist_ok=True)
Path('build/naoesf/polar_coef.py').write_text('\n'.join(out) + '\n')
print('gravado build/naoesf/polar_coef.py')
