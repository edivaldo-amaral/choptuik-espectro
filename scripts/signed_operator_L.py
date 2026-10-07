#!/usr/bin/env python3
"""O operador selecionado L (setor A) na base de Fourier COM SINAL, para S3b:

    L = O_mu + mu S Gamma1 + 2 S Xi Gamma2(RefA, .).

Componentes (setor A): 0 = v1 (n impar, m par), 1 = v2 e 2 = v3 (todo n, m par),
3 = v4 (todo n, m IMPAR). Modo m com sinal, |m| < M: se m e par, carrega as componentes
0, 1, 2; se e impar, so a 3. Norma do toro de raios (65/64, 5/4): peso kappa1^(2m)
omega2(n), m com sinal.

  * O_mu: diagonal mu(n+1), d_tau = i m/2, 2 mu n acima da diagonal (passo 2 na
    componente 0, 1 nas outras) -- a mesma forma de K.
  * mu S Gamma1 (d = 0): componente 0 -> saidas 1 (+mu R) e 2 (-mu R); parte impar de 1
    -> saida 1 (2 mu R); parte impar de 3 -> saida 3 (2 mu R). R: T_n/xi -> 2(-1)^L em
    n - 1 - 2L (desce de qualquer n).
  * 2 S Xi Gamma2(RefA, .): tabela (23) de RT (independent_component_beta.selected_coefficients):
    entrada w1 (fator 1/4) e partes par/impar de w2, w3, w4 (fator 1/2), vezes 2.

Validado contra o exportador de RT (rt-A-*.dat) entrada a entrada por semelhanca (--validar).
"""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from component_exterior import read_fields
from fourier_tail_bound import mult, refl

REF = Path('.cache/rt-1203.3766v1/sourcecode/RefA.dat')
MU = 722873400/2**32
K1, K2 = 65/64, 5/4
BAND_M, BAND_N = 40, 100


def campos():
    """g[j][m] (m >= 0) = coeficientes simetricos n = 0..100 do modo m do campo j (4 campos)."""
    F = read_fields(REF)
    g = []
    for f in F[:4]:
        dd = {}
        for (m, n, part), v in f.items():
            dd.setdefault(m, np.zeros(101, complex))[n] += float(v)*(1j if part else 1)
        g.append(dd)
    return g


def modo(g, j, d):
    a = g[j].get(abs(d))
    if a is None:
        return np.zeros(101, complex)
    return a.conj() if d < 0 else a


# tabela (23): chave de entrada -> (fator, [lista por saida i de termos (coef, campo, refletido)])
TABELA = {
    ('w1', None): (0.25, [[(1, 1, 1), (1, 1, 0), (-1, 2, 1), (-1, 2, 0)],
                          [(1, 0, 0), (-1, 1, 0), (1, 1, 1)],
                          [(-1, 0, 0), (1, 1, 0), (-1, 1, 1)], []]),
    ('w2', 'e'): (0.5, [[(1, 0, 0), (1, 2, 1), (-1, 2, 0)],
                        [(2, 1, 0), (2, 1, 1)],
                        [(-1, 1, 0), (-1, 1, 1)],
                        [(1, 3, 1), (1, 3, 0)]]),
    ('w2', 'o'): (0.5, [[], [(-1, 0, 0), (1, 1, 1), (-1, 1, 0)],
                        [(1, 0, 0)], [(1, 3, 1), (-1, 3, 0)]]),
    ('w3', 'e'): (0.5, [[(-1, 0, 0)], [], [], []]),
    ('w3', 'o'): (0.5, [[(-1, 1, 0), (-1, 1, 1)], [], [], []]),
    ('w4', 'e'): (0.5, [[(1, 3, 0), (-1, 3, 1)], [],
                        [(1, 3, 1), (1, 3, 0)], [(1, 1, 1), (1, 1, 0)]]),
    ('w4', 'o'): (0.5, [[(1, 3, 0), (1, 3, 1)], [],
                        [(1, 3, 1), (-1, 3, 0)], [(1, 1, 1), (-1, 1, 0)]]),
}


def coeficiente(g, termos, d):
    """Modo de Fourier d da funcao-coeficiente sum coef (P?) v_campo."""
    out = np.zeros(101, complex)
    for c, j, p in termos:
        a = modo(g, j, d)
        out += c*(refl(a) if p else a)
    return out


class Radial:
    """Indices radiais de um modo: lista de (componente, array de n)."""
    def __init__(self, m, N):
        self.m = m
        if m % 2 == 0:
            self.comps = [(0, np.arange(1, N, 2)), (1, np.arange(N)), (2, np.arange(N))]
        else:
            self.comps = [(3, np.arange(N))]
        self.off = {}
        o = 0
        for c, ns in self.comps:
            self.off[c] = (o, o + len(ns)); o += len(ns)
        self.dim = o


def bloco_B(g, mo: Radial, mi: Radial, mu=MU):
    """Bloco (modo mo <- modo mi) do fundo mu S Gamma1 + 2 S Xi Gamma2, sem pesos."""
    d = mo.m - mi.m
    X = np.zeros((mo.dim, mi.dim), complex)
    if abs(d) > BAND_M:
        return X
    for cin, nin in mi.comps:
        for paridade in ((None,) if cin == 0 else ('e', 'o')):
            chave = ('w1', None) if cin == 0 else (f'w{cin+1}', paridade)
            fator, saidas = TABELA[chave]
            sel = np.ones(len(nin), bool) if cin == 0 else ((nin % 2 == 0) if paridade == 'e' else (nin % 2 == 1))
            if not sel.any():
                continue
            for cout, nout in mo.comps:
                termos = saidas[cout]
                if not termos:
                    continue
                a = coeficiente(g, termos, d)
                if not np.any(a):
                    continue
                M = 2*fator*mult(a, nin, nout)
                o0, o1 = mo.off[cout]; i0, i1 = mi.off[cin]
                X[o0:o1, i0:i1] += M*sel[None, :]
    if d == 0:                                   # mu S Gamma1
        fontes = {0: {1: 1.0, 2: -1.0}, 1: {1: 2.0}, 3: {3: 2.0}}
        for cin, nin in mi.comps:
            if cin not in fontes:
                continue
            i0, _ = mi.off[cin]
            for col, n in enumerate(nin):
                if cin in (1, 3) and n % 2 == 0:
                    continue
                for cout, fac in fontes[cin].items():
                    if cout not in mo.off:
                        continue
                    o0, _ = mo.off[cout]
                    nout = dict(mo.comps)[cout]
                    pos = {int(x): r for r, x in enumerate(nout)}
                    L = 0
                    while n - 1 - 2*L >= 0:
                        j = n - 1 - 2*L
                        if j in pos:
                            X[o0 + pos[j], i0 + col] += mu*fac*2*(-1)**L
                        L += 1
    return X


def J_livre(r: Radial, s, mu=MU):
    sig = s + 1j*r.m/2
    X = np.zeros((r.dim, r.dim), complex)
    for c, ns in r.comps:
        o0, o1 = r.off[c]
        J = np.diag(mu*(ns + 1.0)).astype(complex) + sig*np.eye(len(ns))
        for b in range(len(ns)):
            J[:b, b] += 2*mu*ns[b]
        X[o0:o1, o0:o1] = J
    return X


def modos(M):
    return list(range(-(M - 1), M))


def operador(M, N, s=0.0, mu=MU, g=None):
    g = campos() if g is None else g
    rs = [Radial(m, N) for m in modos(M)]
    offs = np.cumsum([0] + [r.dim for r in rs])
    L = np.zeros((offs[-1], offs[-1]), complex)
    for j, ri in enumerate(rs):
        L[offs[j]:offs[j+1], offs[j]:offs[j+1]] += J_livre(ri, s, mu)
        for i, ro in enumerate(rs):
            if abs(ro.m - ri.m) <= BAND_M:
                L[offs[i]:offs[i+1], offs[j]:offs[j+1]] += bloco_B(g, ro, ri, mu)
    return L, rs, offs


def validar(malha):
    from independent_crown_tiles import read_matrix
    from check_free_preconditioner import read_spectral_matrix
    p = Path(f'build/spectrum/rt-A-{malha}.dat')
    addr, _ = read_matrix(p)
    dim, raw, expo = read_spectral_matrix(p)
    A = np.array([[float(x) for x in row] for row in raw])*2.0**expo
    M = max(a[2] for a in addr) + 1; N = max(a[3] for a in addr) + 1
    L, rs, offs = operador(M, N)
    pos_real = {a: i for i, a in enumerate(addr)}
    T = np.zeros((len(L), len(addr)), complex)
    for i, r in enumerate(rs):
        for c, ns in r.comps:
            o0, _ = r.off[c]
            for k, n in enumerate(ns):
                linha = offs[i] + o0 + k
                if r.m == 0:
                    T[linha, pos_real[(c, 0, 0, int(n))]] = 1
                else:
                    T[linha, pos_real[(c, 0, abs(r.m), int(n))]] = 1
                    T[linha, pos_real[(c, 1, abs(r.m), int(n))]] = 1j if r.m > 0 else -1j
    dif = L - T @ A @ np.linalg.inv(T)
    print(f'L {malha} (M={M}, N={N}, dim {len(L)} vs {len(addr)}): max |L_sinal - T L_RT T^-1| = '
          f'{np.abs(dif).max():.2e}  (max |L| = {np.abs(L).max():.2e})')
    i, j = np.unravel_index(np.argmax(np.abs(dif)), dif.shape)
    return float(np.abs(dif).max())


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--validar', nargs='+', default=['6x18', '8x24'])
    a = ap.parse_args()
    for malha in a.validar:
        validar(malha)


if __name__ == '__main__':
    main()
