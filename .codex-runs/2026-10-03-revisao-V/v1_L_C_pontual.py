#!/usr/bin/env python3
"""Revisao V (03/10/2026), V1: o L(s) e o C(s) dos certificados (signed_operator_L, signed_constraints)
coincidem com S dOmega^+ e S# dOmega^+ calculados PONTUALMENTE a partir do TeX de RT (eezuuirzr),
no fundo RefA, para h aleatorio do setor A e s complexo.

  S dOmega^+  = ( Xi^-1 (1/2)(1+P) dOmega_1^+, Xi^-1 dOmega_2^+, Xi^-1 dOmega_3^+, Xi^-1 dOmega_4^+ ),
  S# dOmega^+ = ( (1/2)(1-P) dOmega_1^+, -P dOmega_2sharp^+ ),     Xi^-1 f = (f - f|_{xi=0})/xi.

dOmega^+ e obtido por sympy da transcricao em v123_identidades.Omegas (linearizacao exata em
simbolos de valor), e avaliado com os valores do fundo e de h, sem nenhuma tabela C de RT.
Rodar da raiz do repositorio.
"""
import sys
from pathlib import Path
import numpy as np
import sympy as sp
from numpy.polynomial import chebyshev as Ch

sys.path.insert(0, str(Path('scripts').resolve()))
import signed_operator_L as sl  # noqa: E402
import signed_constraints as sc  # noqa: E402

rng = np.random.default_rng(20261003)
MU = sl.MU
s = 0.7 + 0.3j

# ---------- dOmega^+ simbolico em valores ----------
src = Path(__file__).with_name('v123_identidades.py').read_text()
ns = {}
exec(src.split('resultados = {}')[0], ns)
Omegas = ns['Omegas']
xi_s, mu_s, eps = ns['xi'], ns['mu'], ns['eps']
nomes = ['1', '2m', '2p', '3m', '3p', '4m', '4p']
Wv = {k: sp.Symbol('W' + k) for k in nomes}
Hv = {k: sp.Symbol('H' + k) for k in nomes}
DH = {k: sp.Symbol('DH' + k) for k in nomes}       # D_s^+ h_k (a derivada da perturbacao)
DW = {k: sp.Symbol('DW' + k) for k in nomes}
chave = {'1': 'w1', '2m': ('w2', -1), '2p': ('w2', 1), '3m': ('w3', -1), '3p': ('w3', 1),
         '4m': ('w4', -1), '4p': ('w4', 1)}
w = {chave[k]: Wv[k] + eps*Hv[k] for k in nomes}
mapa = {Wv[k]: DW[k] for k in nomes}
mapa.update({Hv[k]: DH[k] for k in nomes})


def Dop(S, f):
    assert S == 1
    return f.xreplace(mapa)


Om = Omegas(w, 1, Dop)
args = [xi_s, mu_s] + [Wv[k] for k in nomes] + [Hv[k] for k in nomes] + [DH[k] for k in nomes]
dOm = {k: sp.lambdify(args, sp.expand(sp.diff(v, eps).subs(eps, 0)), 'numpy') for k, v in Om.items()}

# ---------- fundo (RefA) e h em pontos ----------
g = sl.campos()


def eval_signed(coefs, tau, xi, dtau=False, dxi=False):
    """coefs: dict m -> array c_n (simetricos). f = sum_m e^{i m tau/2} (c_0 + 2 sum c_n T_n)."""
    out = np.zeros(np.broadcast(tau, xi).shape, complex)
    for m, c in coefs.items():
        wgt = np.full(len(c), 2.0); wgt[0] = 1.0
        cc = c*wgt
        if dxi:
            cc = Ch.chebder(cc)
        f = Ch.chebval(xi, cc)
        e = np.exp(0.5j*m*tau)
        if dtau:
            e = e*0.5j*m
        out += e*f
    return out


def fundo(j):
    d = {}
    for m, a in g[j].items():
        d[m] = a
        if m > 0:
            d[-m] = a.conj()
    return d


bg = [fundo(j) for j in range(4)]

# h aleatorio, setor A, |m| <= 2, n <= 5
Nin = 6
h = {c: {} for c in range(4)}
vec = {}
for m in range(-2, 3):
    r = sl.Radial(m, Nin)
    x = np.zeros(r.dim, complex)
    for c, nn in r.comps:
        o0, o1 = r.off[c]
        x[o0:o1] = rng.normal(size=len(nn)) + 1j*rng.normal(size=len(nn))
        full = np.zeros(Nin, complex)
        full[nn] = x[o0:o1]
        h[c][m] = full
    vec[m] = (r, x)


def valores(tau, xi):
    """valores codificados (omega^- = f(xi), omega^+ = -f(-xi)) do fundo e de h, e D_s^+ h."""
    Wd = {'1': eval_signed(bg[0], tau, xi)}
    Hd = {'1': eval_signed(h[0], tau, xi)}
    D = {}
    t1 = eval_signed(h[0], tau, xi, dtau=True); x1 = eval_signed(h[0], tau, xi, dxi=True)
    D['1'] = t1 + s*Hd['1'] + MU*(1 + xi)*x1
    for i, j in ((2, 1), (3, 2), (4, 3)):
        Wd[f'{i}m'] = eval_signed(bg[j], tau, xi)
        Wd[f'{i}p'] = -eval_signed(bg[j], tau, -xi)
        Hd[f'{i}m'] = eval_signed(h[j], tau, xi)
        Hd[f'{i}p'] = -eval_signed(h[j], tau, -xi)
        tm = eval_signed(h[j], tau, xi, dtau=True); xm = eval_signed(h[j], tau, xi, dxi=True)
        tp = -eval_signed(h[j], tau, -xi, dtau=True); xp = eval_signed(h[j], tau, -xi, dxi=True)
        D[f'{i}m'] = tm + s*Hd[f'{i}m'] + MU*(1 + xi)*xm
        D[f'{i}p'] = tp + s*Hd[f'{i}p'] + MU*(1 + xi)*xp
    return Wd, Hd, D


def dOmega_plus(tau, xi):
    Wd, Hd, D = valores(tau, xi)
    a = [xi, MU] + [Wd[k] for k in nomes] + [Hd[k] for k in nomes] + [D[k] for k in nomes]
    return {k: f(*a)*np.ones_like(xi, dtype=complex) for k, f in dOm.items()}


tau = rng.uniform(0, 4*np.pi, 40)
xi = rng.uniform(-1.02, 1.02, 40)
Op, Om_, O0 = dOmega_plus(tau, xi), dOmega_plus(tau, -xi), dOmega_plus(tau, 0*xi)
Xi = lambda f, f0: (f - f0)/xi
sel = [Xi((Op['1'] + Om_['1'])/2, O0['1']), Xi(Op['2'], O0['2']), Xi(Op['3'], O0['3']), Xi(Op['4'], O0['4'])]
con = [(Op['1'] - Om_['1'])/2, -Om_['2s']]

# ---------- L(s) h e C(s) h pelo codigo dos certificados ----------
Nout = 112
Lh = {c: {} for c in range(4)}
Ch_ = {c: {} for c in range(2)}
for mo in range(-44, 45):
    ro = sl.Radial(mo, Nout)
    y = np.zeros(ro.dim, complex)
    for mi, (ri, x) in vec.items():
        y += sl.bloco_B(g, ro, ri) @ x
        if mi == mo:
            J = sl.J_livre(ri, s)
            Jx = J @ x
            for c, nn in ri.comps:
                i0, i1 = ri.off[c]
                o0, _ = ro.off[c]
                pos = {int(n): k for k, n in enumerate(dict(ro.comps)[c])}
                for k, n in enumerate(nn):
                    y[o0 + pos[int(n)]] += Jx[i0 + k]
    for c, nn in ro.comps:
        o0, o1 = ro.off[c]
        full = np.zeros(Nout, complex); full[nn] = y[o0:o1]
        Lh[c][mo] = full
    if mo % 2 == 0:
        so = sc.Saida(mo, Nout)
        z = np.zeros(so.dim, complex)
        for mi, (ri, x) in vec.items():
            z += sc.bloco_C(g, so, ri, s) @ x
        for c, nn in so.comps:
            o0, o1 = so.off[c]
            full = np.zeros(Nout, complex); full[nn] = z[o0:o1]
            Ch_[c][mo] = full

ok = True
for c in range(4):
    v = eval_signed(Lh[c], tau, xi)
    err = np.max(np.abs(v - sel[c]))/max(1, np.max(np.abs(sel[c])))
    print(f'L componente {c}: max |TF(L h) - S dOmega^+| relativo = {err:.2e}   (max |.| = {np.max(np.abs(sel[c])):.3e})')
    ok &= err < 1e-10
for c in range(2):
    v = eval_signed(Ch_[c], tau, xi)
    err = np.max(np.abs(v - con[c]))/max(1, np.max(np.abs(con[c])))
    print(f'C componente {c}: max |TF(C h) - S# dOmega^+| relativo = {err:.2e}   (max |.| = {np.max(np.abs(con[c])):.3e})')
    ok &= err < 1e-10
# controle: trocar s muda o resultado (o teste e sensivel a s)
print('CONFERE' if ok else 'DISCREPANCIA')

# ---------- controles negativos ----------
# (a) s errado na avaliacao pontual; (b) coeficiente quadratico de Omega3 alterado.
s_ok = s
s = s_ok + 0.01
Op2, O02 = dOmega_plus(tau, xi), dOmega_plus(tau, 0*xi)
e_a = np.max(np.abs(eval_signed(Lh[2], tau, xi) - Xi(Op2['3'], O02['3'])))
s = s_ok
Om_alt = dict(Om)
Om_alt['3'] = Om['3'] + sp.Rational(1, 4)*xi_s*w[('w2', 1)]*w['w1']
f_alt = sp.lambdify(args, sp.expand(sp.diff(Om_alt['3'], eps).subs(eps, 0)), 'numpy')


def alt(tau_, xi_):
    Wd, Hd, D = valores(tau_, xi_)
    a = [xi_, MU] + [Wd[k] for k in nomes] + [Hd[k] for k in nomes] + [D[k] for k in nomes]
    return f_alt(*a)*np.ones_like(xi_, dtype=complex)


e_b = np.max(np.abs(eval_signed(Lh[2], tau, xi) - Xi(alt(tau, xi), alt(tau, 0*xi))))
print(f'controle (a) s+0,01: discrepancia {e_a:.2e} (deve ser >> 1e-10)')
print(f'controle (b) Omega3 alterado: discrepancia {e_b:.2e} (deve ser >> 1e-10)')
print('CONTROLES OK' if min(e_a, e_b) > 1e-6 else 'CONTROLES FALHARAM')
