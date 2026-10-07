#!/usr/bin/env python3
"""Conferencia independente das cotas de blocos de tile_perron.

Para linhas e colunas dentro de F, o bloco P_i E P_j do operador infinito coincide com o
da caixa F: Q0 e diagonal em m e triangular superior em n, logo Q0 P_j fica em F, e as
linhas pedidas estao em F. Monta-se entao, na caixa F e com o fundo COMPLETO (sem truncar
a banda), H(s) = K(s) P com P = Q0(sZ) P_Z + Q0(s1) P_T1 + Q0(s2) P_T2, e
E = I - diag(V, I, I) H(s), com o mesmo V do certificado. Para varios s no tile,
||P_i E P_j||_2 (verdadeiro, por SVD) tem de ficar abaixo da cota N_ij do certificado
avaliada com as distancias |s - s_l|.

Niveis testados: Z, T1, T2 (far nao cabe na caixa; suas entradas sao as cotas da forma
fechada e do simbolo).
"""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import tile_perron as tp
from signed_operator import operador, Jrad, modos, K1


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--Z', default='8x24'); ap.add_argument('--Gp', default='12x36')
    ap.add_argument('--F', default='32x77'); ap.add_argument('--D', type=int, default=4)
    ap.add_argument('--sZ', type=complex, default=0.32); ap.add_argument('--rZ', type=float, default=0.02)
    ap.add_argument('--s1', type=complex, default=0.3); ap.add_argument('--r1', type=float, default=0.05)
    ap.add_argument('--s2', type=complex, default=0.25); ap.add_argument('--r2', type=float, default=0.12)
    ap.add_argument('--rigoroso', action='store_true')
    a = ap.parse_args()
    tp.RIG = a.rigoroso
    cx = lambda t: tuple(map(int, t.split('x')))
    Z, Gp, F = cx(a.Z), cx(a.Gp), cx(a.F)
    (MZ, NZ), (MP, NP), (Mc, Nc) = Z, Gp, F
    silencio = lambda *x, **y: None
    ctx = tp.Ctx(Z, Gp, F, a.D)
    qZ, V = tp.parte_Z(ctx, a.sZ, silencio)
    q1, HZ1, QZ1 = tp.parte_T1(ctx, a.s1, silencio)
    q2, SZ, SQZ = tp.parte_T2(ctx, a.s2, silencio)
    N0, NZm, N1m, N2m, qa = tp.montar(ctx, qZ, V, q1, HZ1, QZ1, q2, SZ, SQZ, a.rZ, a.r1, a.r2, silencio)

    # ---- a verdade na caixa F, fundo completo, pesos do toro com sinal
    ms = modos(Mc)
    md = tp.Modo(Nc); k = len(md.n)
    peso = np.concatenate([K1**m*md.w for m in ms])
    niv = np.empty(len(ms)*k, int)
    for i, m in enumerate(ms):
        loc = np.full(k, 2)                                   # T2
        if abs(m) < MP:
            loc[md.n < NP] = 1                                # T1
        if abs(m) < MZ:
            loc[md.n < NZ] = 0                                # Z
        niv[i*k:(i+1)*k] = loc
    idx = [np.where(niv == l)[0] for l in range(3)]
    K0, _ = operador(Mc, Nc, s=0.0)                           # K(0) na caixa, fundo completo
    Kf = np.zeros_like(K0)                                    # K_livre (sem s)
    for i, m in enumerate(ms):
        Kf[i*k:(i+1)*k, i*k:(i+1)*k] = Jrad(m, 0.0, Nc)       # J_rad + i m/2
    W = lambda X: (peso[:, None]*X)/peso[None, :]
    I = np.eye(len(K0))
    Q = {c: np.linalg.inv(Kf + c*I) for c in {a.sZ, a.s1, a.s2}}
    P = np.zeros_like(K0)
    for l, c in enumerate((a.sZ, a.s1, a.s2)):
        P[:, idx[l]] = Q[c][:, idx[l]]
    # a ordem dos indices de Z na caixa F e a mesma do certificado (modos crescentes, c1 depois c2)
    pior = 0.0
    rng = np.random.default_rng(0)
    pontos = [a.sZ] + [a.sZ + a.rZ*np.exp(2j*np.pi*t) for t in rng.random(5)]
    for s in pontos:
        H = W((K0 + s*I) @ P)
        E = -H.copy()
        E[np.diag_indices(len(E))] += 1
        EZ = E[idx[0]]
        E[idx[0]] = EZ - (V - np.eye(len(V))) @ H[idx[0]]    # linhas de Z: I - V H
        dist = (abs(s - a.sZ), abs(s - a.s1), abs(s - a.s2))
        N = N0 + dist[0]*NZm + dist[1]*N1m + dist[2]*N2m
        linhas = []
        for i_ in range(3):
            for j_ in range(3):
                verdade = np.linalg.norm(E[np.ix_(idx[i_], idx[j_])], 2)
                cota = N[i_, j_]
                razao = verdade/cota if cota > 0 else (np.inf if verdade > 1e-12 else 0.0)
                pior = max(pior, razao)
                linhas.append(f'{"ZT1T2"[[0, 1, 3][i_]:[1, 3, 5][i_]]:>2}<-{"ZT1T2"[[0, 1, 3][j_]:[1, 3, 5][j_]]:<2}'
                              f' {verdade:.4f}/{cota:.4f}')
        print(f's = {s:.4f}: ' + '  '.join(linhas))
    print(f'pior razao verdade/cota = {pior:.4f}  ->  ' + ('OK' if pior <= 1 else 'COTA VIOLADA'))


if __name__ == '__main__':
    main()
