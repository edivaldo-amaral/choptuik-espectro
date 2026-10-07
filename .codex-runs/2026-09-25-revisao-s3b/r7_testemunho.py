#!/usr/bin/env python3
"""R7 (revisao independente de S3b): o testemunho de constraints (A6).

  (1) c_ref = ||Pi_F C(lambda0) h0|| com montagem propria (laco, pesos e Pi_F escritos aqui;
      so os blocos sem peso vem de signed_constraints.bloco_C, validado contra o exportador de RT);
  (2) c_J = ||Pi_F C Q0 P_J|| DIRETO (SVD densa) para as janelas 32..47 e -15..0, e o porque de
      serem pequenos (modos de Fourier do fundo em |d| >= 21 e descida radial n >= 128 -> n < 48);
  (3) a perda de lambda: ||Pi_F C1|| (C1 h = C(1)h - C(0)h) medido numa caixa, contra nC1 = 1,5.
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path('scripts').resolve()))
import signed_operator_L as sl
import signed_constraints as sc

K1L, K2L = 65/64, 5/4
K1K, K2K = 129/128, 9/8
LAM0 = 0.4010247330
MZ, NZ, NC = 32, 128, 400
MF, NF = 12, 48


def om(n, k):
    n = np.asarray(n, float)
    return np.where(n == 0, 1.0, k**(2*n) + k**(-2*n))


def wL(m, N):
    r = sl.Radial(m, N)
    return K1L**m*np.sqrt(np.concatenate([om(ns, K2L) for _, ns in r.comps]))


def wK(m, N):
    so = sc.Saida(m, N)
    return K1K**m*np.sqrt(np.concatenate([om(ns, K2K) for _, ns in so.comps]))


G = sl.campos()
saidas = [m for m in range(-(MF - 1), MF) if m % 2 == 0]


def Cw(mo, mi, Nin, s):
    X = sc.bloco_C(G, sc.Saida(mo, NF), sl.Radial(mi, Nin), s)
    return (wK(mo, NF)[:, None]*X)/wL(mi, Nin)[None, :]


def main():
    rig = Path(sys.argv[1]) if len(sys.argv) > 1 else Path('build/s3b/rig')
    h0w = np.load(rig/'h0w.npy')
    msZ = list(range(-(MZ - 1), MZ))
    offs = np.cumsum([0] + [sl.Radial(m, NZ).dim for m in msZ])
    assert offs[-1] == len(h0w)
    # (1) c_ref
    partes = []
    for mo in saidas:
        acc = np.zeros(sc.Saida(mo, NF).dim, complex)
        for j, mi in enumerate(msZ):
            if abs(mo - mi) <= 40:
                acc += Cw(mo, mi, NZ, LAM0) @ h0w[offs[j]:offs[j+1]]
        partes.append(acc)
    c = np.concatenate(partes)
    print(f'(1) ||Pi_F C(lambda0) h0|| (montagem propria, F = {MF}x{NF}) = {np.linalg.norm(c):.8f}'
          f'   (D.json: c_ref >= 0,171739, ja descontados erros)')
    # sensibilidade: caixa de saida maior nao e necessaria (Pi_F e a projecao); conferir |h0| = 1
    print(f'    ||h0|| = {np.linalg.norm(h0w):.12f};  maior bloco de saida: modo {saidas[int(np.argmax([np.linalg.norm(p) for p in partes]))]}')
    # (2) c_J direto
    D = json.loads((rig/'D.json').read_text()); T3 = json.loads((rig/'T3.json').read_text())
    for nome in ('32..47', '-15..0', '16..31'):
        a, b = nome.split('..'); ms = list(range(int(a), int(b) + 1))
        blocos = []
        for mi in ms:
            r = sl.Radial(mi, NC)
            nn = np.concatenate([ns for _, ns in r.comps])
            cols = np.where(nn >= NZ)[0] if abs(mi) < MZ else np.arange(r.dim)
            w = wL(mi, NC)/K1L**mi
            Q = np.linalg.solve(sl.J_livre(r, LAM0), np.eye(r.dim))
            Qw = (w[:, None]*Q)/w[None, :]
            linhas = []
            for mo in saidas:
                if abs(mo - mi) <= 40:
                    linhas.append(Cw(mo, mi, NC, LAM0) @ Qw[:, cols])
                else:
                    linhas.append(np.zeros((sc.Saida(mo, NF).dim, len(cols)), complex))
            blocos.append(np.concatenate(linhas))
        X = np.concatenate(blocos, axis=1)
        v = np.linalg.norm(X, 2)
        i = T3['nome'].index(nome)
        print(f'(2) c_J[{nome}] direto (SVD, {X.shape}) = {v:.4e};   D.json: {D["cJ"][i]:.4e}')
    # por que pequenos: coeficientes do fundo em |d| = 21 (sem peso) e em d = 0
    for d in (0, 10, 21, 30):
        mx = max(float(np.max(np.abs(sl.modo(G, j, d)))) for j in range(4))
        print(f'    max |coeficiente do fundo| no modo de Fourier d = {d}: {mx:.2e}')
    # (3) ||Pi_F C1|| numa caixa (entrada |m| < 24, n < 120)
    ms = list(range(-23, 24))
    rows = []
    for mo in saidas:
        linha = [Cw(mo, mi, 120, 1.0) - Cw(mo, mi, 120, 0.0) if abs(mo - mi) <= 40 else
                 np.zeros((sc.Saida(mo, NF).dim, sl.Radial(mi, 120).dim)) for mi in ms]
        rows.append(np.concatenate(linha, axis=1))
    C1 = np.concatenate(rows)
    print(f'(3) ||Pi_F C1|| (secao finita, entrada 24x120) = {np.linalg.norm(C1, 2):.4f}   (usado: nC1 = 1,5)')


if __name__ == '__main__':
    main()
