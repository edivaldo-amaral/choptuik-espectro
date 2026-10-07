#!/usr/bin/env python3
"""Monta o operador de constraints K (setor A) numa caixa M x N, DIRETO na memoria,
a partir dos campos de RefA -- sem passar pelo exportador de texto, inviavel em
dezenas de milhares de DOFs.

    K = K_livre + mu R_x + Gamma2sharp(omega, .)

K_livre: diag mu(n+1), d_tau = i m/2 (acopla as partes real/imaginaria), e 2 mu n
para todo n' < n (c2) ou n' < n de mesma paridade (c1).
mu R_x: saida c2, entrada c1: divisao regularizada por xi em coeficientes simetricos,
T_{2k+1}/xi -> 2(-1)^(k-j) em n = 2j, 0 <= j <= k.
Gamma2sharp (RT): (1/2) F1 x1 + (1/2) F2 E(x2) + (1/2) F3 O(x2), com
F1 = (v2+Pv2-Pv3-v3, v1-v2+Pv2), F2 = (-v1+Pv3-v3, -v1+v2+3Pv2),
F3 = (v2+Pv2+Pv3+v3, -v1+v2+3Pv2); P: coeficiente n -> (-1)^n.

Coordenadas reais do exportador: (d, parte, m, n), m par (setor A), m = 0 so parte 0,
c1 (d=0) com n impar, c2 (d=1) com todo n. Coeficientes simetricos em n e
conjugados em m. Validado contra rt-sharp-12x36.dat (ver --validar).
"""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from component_exterior import read_fields

REF = Path('.cache/rt-1203.3766v1/sourcecode/RefA.dat')
MU = 722873400/2**32
MB, NB = 40, 100                     # suporte do fundo em m e n


def enderecos(M: int, N: int):
    out = []
    for m in range(0, M, 2):
        for part in ((0,) if m == 0 else (0, 1)):
            for n in range(N):
                if n % 2 == 1:
                    out.append((0, part, m, n))
                out.append((1, part, m, n))
    return out


def campos_z2():
    """Arrays complexos Z^2 dos quatro campos: G[j][m+MB, n+NB], m em [-MB,MB], n em [-NB,NB]."""
    F = read_fields(REF)
    out = []
    for f in F:
        G = np.zeros((2*MB+1, 2*NB+1), complex)
        for (m, n, part), v in f.items():
            val = float(v)*(1j if part else 1)
            for nn in {n, -n}:
                G[MB+m, NB+nn] += val
                if m:
                    G[MB-m, NB+nn] += np.conj(val)
        out.append(G)
    return out


def P(G):
    sinais = (-1.0)**np.abs(np.arange(-NB, NB+1))
    return G*sinais[None, :]


def mult(G, ent, sai):
    """Matriz real (len(sai) x len(ent)) da multiplicacao por G, entradas/saidas em DOFs reais."""
    sai = np.asarray(sai); ent = np.asarray(ent)
    mo, no, po = sai[:, 2], sai[:, 3], sai[:, 1]
    R = np.zeros((len(sai), len(ent)))
    for col, (_, part, m, n) in enumerate(ent):
        u = 1.0 if part == 0 else 1j
        fontes = [(m, n, u)]
        if n:
            fontes.append((m, -n, u))
        if m:
            fontes += [(-m, s_n, np.conj(u)) for (_, s_n, _) in list(fontes)]
        acc = np.zeros(len(sai), complex)
        for ms, ns, us in fontes:
            dm, dn = mo - ms, no - ns
            ok = (np.abs(dm) <= MB) & (np.abs(dn) <= NB)
            acc[ok] += G[MB+dm[ok], NB+dn[ok]]*us
        R[:, col] = np.where(po == 0, acc.real, acc.imag)
    return R


def operador(M: int, N: int, mu: float = MU, s: complex = 0.0):
    addr = enderecos(M, N)
    idx = {a: i for i, a in enumerate(addr)}
    n_ = len(addr)
    K = np.zeros((n_, n_))
    # livre
    for col, (d, part, m, n) in enumerate(addr):
        K[col, col] += mu*(n + 1)
        if m:
            K[idx[(d, 1-part, m, n)], col] += (m/2)*(1 if part == 0 else -1)
        passo = 2 if d == 0 else 1
        for low in range(n - passo, -1, -passo):
            K[idx[(d, part, m, low)], col] += 2*mu*n
    # mu R_x: c2 (n par) <- c1 (n impar)
    for col, (d, part, m, n) in enumerate(addr):
        if d != 0:
            continue
        k = (n - 1)//2
        for j in range(k + 1):
            K[idx[(1, part, m, 2*j)], col] += mu*2*(-1)**(k - j)
    # Gamma2sharp
    v1, v2, v3, _ = campos_z2()
    P1, P2, P3 = P(v1), P(v2), P(v3)
    F1 = (v2+P2-P3-v3, v1-v2+P2)
    F2 = (-v1+P3-v3, -v1+v2+3*P2)
    F3 = (v2+P2+P3+v3, -v1+v2+3*P2)
    sai = [[i for i, a in enumerate(addr) if a[0] == dd] for dd in (0, 1)]
    c1 = [i for i, a in enumerate(addr) if a[0] == 0]
    c2p = [i for i, a in enumerate(addr) if a[0] == 1 and a[3] % 2 == 0]
    c2i = [i for i, a in enumerate(addr) if a[0] == 1 and a[3] % 2 == 1]
    for dd in (0, 1):
        rows = sai[dd]; ra = [addr[i] for i in rows]
        for cols, G in ((c1, F1[dd]), (c2p, F2[dd]), (c2i, F3[dd])):
            K[np.ix_(rows, cols)] += 0.5*mult(G, [addr[i] for i in cols], ra)
    if s:
        K = K + s*np.eye(n_)
    return K, addr


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--validar', type=Path, default=Path('build/spectrum/rt-sharp-12x36.dat'))
    a = ap.parse_args()
    from independent_sharp_matrix import read_sharp
    ea, E = read_sharp(a.validar)
    M = 2*(max(x[2] for x in ea)//2 + 1); N = max(x[3] for x in ea) + 1
    K, addr = operador(M, N)
    assert sorted(addr) == sorted(ea), 'enderecos diferentes do exportador'
    pos = {x: i for i, x in enumerate(addr)}
    perm = [pos[x] for x in ea]
    Km = K[np.ix_(perm, perm)]
    Ef = np.array([[float(x) for x in row] for row in E])
    dif = np.abs(Km - Ef)
    print(f'caixa {M}x{N}, dim {len(addr)}: max |mine - exportado| = {dif.max():.3e}'
          f'  (max |exportado| = {np.abs(Ef).max():.3e})')
    i, j = np.unravel_index(np.argmax(dif), dif.shape)
    print(f'  pior entrada: linha {ea[i]}, coluna {ea[j]}: mine {Km[i, j]:.6e}, exportado {Ef[i, j]:.6e}')


if __name__ == '__main__':
    main()
