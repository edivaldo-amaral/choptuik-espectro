#!/usr/bin/env python3
"""Certificado de um tile para o operador de constraints K (S3a): K(s) injetivo
para |s - s0| <= r. Norma: L^2 do toro de raios (129/128, 9/8).

    H(s) = K(s) Q0 = I + (B + ds) Q0,   Q0 = (K_livre + s0)^-1,   B = mu R_x + Gamma2sharp.

Tres niveis: G (bloco interno, inverso aproximado V), near = G+ \\ G (faixa), far
(o resto, infinito). Inverso aproximado A = diag(V, I, I). Se a norma espectral da
matriz 3x3 de cotas dos blocos de E = I - A H(s) for < 1, A H e invertivel e K(s)
e injetivo em todo o tile.

  * blocos entre G e near: FINITOS e calculados exatamente, porque Q0 preserva
    caixas (K_livre so desce em n e nao muda m);
  * G <-> far: pequenos por construcao (o fundo truncado nao alcanca; o inverso
    livre so desce com peso (9/8)^-dn; o resto do fundo e ~1e-4);
  * far -> near: finito (o fundo so sobe 20 em m e 40 em n);
  * beta_f = ||(B + ds) Q0 P_far||: composicoes radiais 1D modo a modo
    (fourier_tail_bound).

ESTADO: prototipo. As normas finitas sao calculadas em ponto flutuante com margem
explicita; ver a secao de rigor no relatorio.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys
import time

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from sharp_operator_builder import operador, enderecos, MU

K1, K2 = 129/128, 9/8


def omega(k, j):
    return 1.0 if j == 0 else k**(2*j) + k**(-2*j)


def livre(addr, mu=MU):
    idx = {x: i for i, x in enumerate(addr)}
    n = len(addr)
    Kf = np.zeros((n, n))
    for col, (d, part, m, nn) in enumerate(addr):
        Kf[col, col] += mu*(nn + 1)
        if m:
            Kf[idx[(d, 1-part, m, nn)], col] += (m/2)*(1 if part == 0 else -1)
        passo = 2 if d == 0 else 1
        for low in range(nn - passo, -1, -passo):
            Kf[idx[(d, part, m, low)], col] += 2*mu*nn
    return Kf


def inverso_livre(Kf, addr, s0):
    """(K_livre + s0)^-1 por blocos (d, m): bloco-diagonal, exato a menos de arredondamento."""
    n = len(addr)
    Q = np.zeros((n, n), complex)
    grupos = {}
    for i, (d, part, m, nn) in enumerate(addr):
        grupos.setdefault((d, m), []).append(i)
    for idx in grupos.values():
        ix = np.ix_(idx, idx)
        Q[ix] = np.linalg.inv(Kf[ix] + s0*np.eye(len(idx)))
    return Q


def norma2(X, iters=80, seed=0):
    """Estimativa por potencia (cota INFERIOR); a cota superior vem de norma2_sup."""
    if X.size == 0:
        return 0.0
    rng = np.random.default_rng(seed)
    v = rng.standard_normal(X.shape[1]) + 0j
    v /= np.linalg.norm(v)
    for _ in range(iters):
        u = X.conj().T @ (X @ v)
        nu = np.linalg.norm(u)
        if nu == 0:
            return 0.0
        v = u/nu
    return float(np.linalg.norm(X @ v))


def norma2_sup(X):
    """Cota SUPERIOR de ||X||_2: t com t^2 I - X*X definida positiva (Cholesky), com
    margem para o arredondamento do produto e da fatoracao (n u ||X||_F^2)."""
    if X.size == 0:
        return 0.0
    G = X.conj().T @ X
    fro2 = float(np.sum(np.abs(X)**2))
    n = G.shape[0]
    folga = 4*(n + 2)*2.0**-53*fro2
    t = norma2(X)*1.0005 + 1e-14
    for _ in range(40):
        try:
            np.linalg.cholesky(t*t*np.eye(n) - G - folga*np.eye(n))
            return t
        except np.linalg.LinAlgError:
            t *= 1.01
    return float(np.linalg.norm(X, 'fro'))


def parte_finita(s0, G, Gp, r):
    MG, NG = G; MP, NP = Gp
    t0 = time.time()
    K, addr = operador(MP, NP)
    Kf = livre(addr)
    B = K - Kf
    del K
    Q0 = inverso_livre(Kf, addr, s0)
    dv = np.sqrt(np.array([omega(K1, m)*omega(K2, nn) for _, _, m, nn in addr]))
    W = lambda X: (dv[:, None]*X)/dv[None, :]
    H = W(np.eye(len(addr)) + B @ Q0)
    Qw = W(Q0)
    iG = np.array([i for i, (_, _, m, nn) in enumerate(addr) if m < MG and nn < NG])
    iN = np.array([i for i, (_, _, m, nn) in enumerate(addr) if not (m < MG and nn < NG)])
    HGG = H[np.ix_(iG, iG)]
    V = np.linalg.inv(HGG)
    blk = lambda X, a, b: X[np.ix_(a, b)]
    q = {}
    q['rho'] = norma2_sup(np.eye(len(iG)) - V @ HGG)
    q['V'] = norma2_sup(V)
    q['e_Gn'] = norma2_sup(V @ blk(H, iG, iN))
    q['e_nG'] = norma2_sup(blk(H, iN, iG))
    q['e_nn'] = norma2_sup(blk(H, iN, iN) - np.eye(len(iN)))
    q['VQ_GG'] = norma2_sup(V @ blk(Qw, iG, iG))
    q['VQ_Gn'] = norma2_sup(V @ blk(Qw, iG, iN))
    q['Q_nn'] = norma2_sup(blk(Qw, iN, iN))
    q['Q_GG'] = norma2_sup(blk(Qw, iG, iG))
    q['dim'] = dict(G=len(iG), near=len(iN))
    q['tempo'] = time.time() - t0
    return q


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--s0', type=complex, default=1.0)
    ap.add_argument('--r', type=float, default=0.0)
    ap.add_argument('--G', default='12x36')
    ap.add_argument('--Gp', default='20x60')
    ap.add_argument('--out', type=Path)
    a = ap.parse_args()
    G = tuple(map(int, a.G.split('x'))); Gp = tuple(map(int, a.Gp.split('x')))
    q = parte_finita(a.s0, G, Gp, a.r)
    for k, v in q.items():
        print(f'  {k:>8} = {v}' if not isinstance(v, float) else f'  {k:>8} = {v:.5f}')
    if a.out:
        a.out.write_text(json.dumps(dict(s0=str(a.s0), r=a.r, G=G, Gp=Gp, finita=q), indent=2) + '\n')


if __name__ == '__main__':
    main()
