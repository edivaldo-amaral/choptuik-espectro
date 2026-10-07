#!/usr/bin/env python3
"""O operador de constraints C(s) = C0 + s C1 (do espaco de L para o de K) na base de
Fourier COM SINAL, para o testemunho de S3b. Decomposicao (CONSTRAINT_WITNESS_NORMS.md):

    C(s)h = D(s)h + A_w h,
    D(s)h = ( mu(xi dxi + 2)h1 ,  -xi(dtau + s + mu)h2 + mu(xi dxi - xi^2 dxi)h2 + 2mu(h2 - h3) ),
    A_w h = ( 2 O(xi Gamma2_1(w, h)) ,  -2 P(xi Gamma2_3(w, h)) ),

com os coeficientes de Gamma2_1 e Gamma2_3 da tabela de H3
(independent_constraint_operator_bound.background_combo_norms). Entrada: layout de L
(signed_operator_L); saida: layout de K (componente 0 com n impar, 1 com todo n; m par).
Validado contra o exportador de RT (rt-A-*.dat.constraints) entrada a entrada (--validar).
"""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import signed_operator_L as sl
from fourier_tail_bound import mult, refl

MU = sl.MU

# (coeficiente, campo 0..3, refletido) -- da tabela de H3, campos 1..4 -> 0..3
COMBOS = {
    'k1': [(3, 1, 1), (-2, 2, 1), (1, 0, 0), (-1, 1, 0)],
    'k2e_r1': [(1, 1, 1), (2, 2, 1), (1, 0, 0), (1, 1, 0)],
    'k2e_r2': [(-2, 1, 1), (2, 2, 1)],
    'k2o_r1': [(-2, 0, 0)],
    'k2o_r2': [(2, 1, 1), (-2, 2, 1)],
    'k3e_r1': [(1, 1, 1), (-1, 0, 0), (1, 1, 0)],
    'k3e_r2': [(2, 1, 1)],
    'k3o_r1': [(1, 0, 0), (-1, 1, 1), (-1, 1, 0)],
    'k3o_r2': [(-2, 1, 1)],
    'k4_r1': [(2, 3, 1)],
    'k4_r2': [(2, 3, 1)],
}
# (componente de entrada, paridade de n ou None) -> [(saida 0 = Gamma2_1 ou 1 = Gamma2_3, combo, fator)]
MAPA = {
    (0, None): [(0, 'k1', 0.25)],
    (1, 'e'): [(0, 'k2e_r1', 0.5), (1, 'k2e_r2', 0.5)],
    (1, 'o'): [(0, 'k2o_r1', 0.5), (1, 'k2o_r2', 0.5)],
    (2, 'e'): [(0, 'k3e_r1', 0.5), (1, 'k3e_r2', 0.5)],
    (2, 'o'): [(0, 'k3o_r1', 0.5), (1, 'k3o_r2', 0.5)],
    # h4: a combinacao k4 age sobre -P h4 (sinal - na parte par, + na impar); a tabela de H3
    # registra so o modulo, que basta para normas. Fixado contra o exportador de RT.
    (3, 'e'): [(0, 'k4_r1', -0.5), (1, 'k4_r2', -0.5)],
    (3, 'o'): [(0, 'k4_r1', 0.5), (1, 'k4_r2', 0.5)],
}


def comb(g, termos, d):
    out = np.zeros(101, complex)
    for c, j, p in termos:
        a = sl.modo(g, j, d)
        out += c*(refl(a) if p else a)
    return out


def xi_vezes(a):
    """coeficientes simetricos de xi*f, f = a_0 + 2 sum a_n T_n: (xi f)_n = (a_{n-1} + a_{n+1})/2."""
    L = len(a)
    b = np.zeros(L + 1, complex)
    for n in range(L + 1):
        s = 0
        if n - 1 >= 0 and n - 1 < L:
            s += a[n - 1]
        elif n - 1 == -1:
            s += a[1] if L > 1 else 0
        if n + 1 < L:
            s += a[n + 1]
        b[n] = s/2
    return b


def xidxi(nin, nout):
    """matriz de xi d/dxi em coeficientes simetricos (independent_constraint_operator_bound.xidxi_apply)."""
    pos = {int(n): r for r, n in enumerate(nout)}
    X = np.zeros((len(nout), len(nin)))
    for c, n in enumerate(nin):
        n = int(n)
        if n in pos:
            X[pos[n], c] += n
        for j in range(n - 2, -1, -2):
            if j in pos:
                X[pos[j], c] += 2*n
    return X


class Saida:
    """Layout de K num modo par m: (0, n impar), (1, todo n)."""
    def __init__(self, m, N):
        self.m = m
        self.comps = [(0, np.arange(1, N, 2)), (1, np.arange(N))]
        self.off = {0: (0, len(self.comps[0][1])), 1: (len(self.comps[0][1]), len(self.comps[0][1]) + N)}
        self.dim = self.off[1][1]


def bloco_C(g, so: Saida, ri: sl.Radial, s=0.0, mu=MU):
    """Bloco (modo de saida so <- modo de entrada ri) de C(s), sem pesos."""
    d = so.m - ri.m
    X = np.zeros((so.dim, ri.dim), complex)
    nout0, nout1 = so.comps[0][1], so.comps[1][1]
    xi = np.zeros(2, complex); xi[1] = 0.5            # xi = T_1: coeficiente simetrico 1/2
    if d == 0:                                       # D(s)
        for cin, nin in ri.comps:
            i0, i1 = ri.off[cin]
            if cin == 0:
                o0, o1 = so.off[0]
                X[o0:o1, i0:i1] += mu*(xidxi(nin, nout0) + 2*mult(np.array([1.0]), nin, nout0))
            if cin == 1:
                o0, o1 = so.off[1]
                sig = 1j*so.m/2 + s + mu
                Xi = mult(xi, nin, nout1)
                nmid = np.arange(len(nout1) + 2)
                X[o0:o1, i0:i1] += (-sig*Xi + mu*xidxi(nin, nout1)
                                    - mu*mult(xi, nmid, nout1) @ xidxi(nin, nmid)
                                    + 2*mu*mult(np.array([1.0]), nin, nout1))
            if cin == 2:
                o0, o1 = so.off[1]
                X[o0:o1, i0:i1] += -2*mu*mult(np.array([1.0]), nin, nout1)
    for cin, nin in ri.comps:                        # A_w
        i0, i1 = ri.off[cin]
        for par in ((None,) if (cin, None) in MAPA else ('e', 'o')):
            sel = np.ones(len(nin), bool) if par is None else ((nin % 2 == 0) if par == 'e' else (nin % 2 == 1))
            for saida, nome, fator in MAPA[(cin, par)]:
                a = comb(g, COMBOS[nome], d)
                if not np.any(a):
                    continue
                c = xi_vezes(a)                       # xi * Gamma2_(...)
                if saida == 0:
                    o0, o1 = so.off[0]
                    M = 2*fator*mult(c, nin, nout0)   # 2 O(...): linhas impares ja sao a projecao
                else:
                    o0, o1 = so.off[1]
                    sinais = (-1.0)**nout1
                    M = -2*fator*sinais[:, None]*mult(c, nin, nout1)
                X[o0:o1, i0:i1] += M*sel[None, :]
    return X


def operador_C(M, N, s=0.0, g=None, mu=MU):
    g = sl.campos() if g is None else g
    ris = [sl.Radial(m, N) for m in sl.modos(M)]
    sos = [Saida(m, N) for m in sl.modos(M) if m % 2 == 0]
    oi = np.cumsum([0] + [r.dim for r in ris]); oo = np.cumsum([0] + [r.dim for r in sos])
    C = np.zeros((oo[-1], oi[-1]), complex)
    for i, so in enumerate(sos):
        for j, ri in enumerate(ris):
            if abs(so.m - ri.m) <= sl.BAND_M:
                C[oo[i]:oo[i+1], oi[j]:oi[j+1]] = bloco_C(g, so, ri, s, mu)
    return C, sos, ris, oo, oi


def validar(malha):
    import independent_numeric_audit as NA
    L, addr = NA.read_rt_matrix(NA.SPEC/f'rt-A-{malha}.dat', 'CHOPTUIK_RT_SPECTRAL_MATRIX_V1')
    C0, C1, caddr = NA.read_constraints(NA.SPEC/f'rt-A-{malha}.dat.constraints', len(L))
    M = max(a[2] for a in addr) + 1; N = max(a[3] for a in addr) + 1
    res = {}
    for s, Cr in ((0.0, C0), (1.0, C0 + C1)):
        C, sos, ris, oo, oi = operador_C(M, N, s)
        pin = {a: i for i, a in enumerate(addr)}; pout = {a: i for i, a in enumerate(caddr)}
        Tin = np.zeros((oi[-1], len(addr)), complex); Tout = np.zeros((oo[-1], len(caddr)), complex)
        for T, blocos, offs, pos in ((Tin, ris, oi, pin), (Tout, sos, oo, pout)):
            for i, r in enumerate(blocos):
                for comp, ns in r.comps:
                    o0, _ = r.off[comp]
                    for k, n in enumerate(ns):
                        lin = offs[i] + o0 + k
                        if r.m == 0:
                            T[lin, pos[(comp, 0, 0, int(n))]] = 1
                        else:
                            T[lin, pos[(comp, 0, abs(r.m), int(n))]] = 1
                            T[lin, pos[(comp, 1, abs(r.m), int(n))]] = 1j if r.m > 0 else -1j
        Cs = Tout @ Cr @ np.linalg.inv(Tin)
        dif = np.abs(C - Cs)
        res[s] = dif.max()
        print(f'C({s}) {malha}: max |C_sinal - T C_RT T^-1| = {dif.max():.2e}  (max |C| {np.abs(Cs).max():.2e})')
        if dif.max() > 1e-10:
            i, j = np.unravel_index(np.argmax(dif), dif.shape)
            print(f'   pior: linha {i}, coluna {j}: sinal {C[i, j]:.6f}  RT {Cs[i, j]:.6f}')
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--validar', nargs='+', default=['6x18'])
    a = ap.parse_args()
    for m in a.validar:
        validar(m)


if __name__ == '__main__':
    main()
