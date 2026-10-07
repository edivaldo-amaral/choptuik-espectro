#!/usr/bin/env python3
"""Prototipo em FLOAT (DIAGNOSTICO, nao e prova) das arquiteturas de parametriz.

Caixa reduzida: matriz exportada rt-A-12x36.dat (1422 DOFs), H(s0)=W(A+s0)Q W^-1 na
norma eta x pesos RT (soma por coluna, |Re|+|Im}). Nucleo C=6x18 (333), casca O=resto.
Parametrizes P (defeito ||I-HP|| por colunas e raio espectral):
  P0  = H_CC^-1 (+) I                             (bloco, casca = identidade)
  P1  = H_CC^-1 (+) (I - X_OO)                    (estilo RT: 1a ordem na casca)
  P2k = Schur: casca V_O = omega sum_{j<k}(I-omega H_OO)^j, nucleo S^-1, S=H_CC-H_CO V_O H_OC
  P3k = refinamento iterativo de P1: P1 sum_{j<k}(I-H P1)^j
Mede tambem o tempo de aplicar H a uma coluna (matriz densa) como referencia.
"""
from __future__ import annotations

import json
import math
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
MU = 722873400/2**32
ETA = [916., 4096., 243., 1233.]


def read(path):
    with open(path) as fh:
        assert fh.readline().strip() == 'CHOPTUIK_RT_SPECTRAL_MATRIX_V1'
        dim = int(fh.readline().split()[1])
        line = fh.readline()
        while line.strip() != 'BEGIN_ADDRESSES':
            line = fh.readline()
        addr = [tuple(map(int, fh.readline().split()))[1:] for _ in range(dim)]
        assert fh.readline().strip() == 'BEGIN_ENTRIES'
        A = np.zeros((dim, dim))
        for line in fh:
            f = line.split()
            if f:
                A[int(f[0]), int(f[1])] += math.ldexp(float(int(f[2])), int(f[3]) if len(f) > 3 else 0)
    return A, addr


def free(addr):
    look = {a: i for i, a in enumerate(addr)}
    J = np.zeros((len(addr), len(addr)))
    for c, (d, p, m, n) in enumerate(addr):
        J[c, c] += MU*(n+1)
        if m:
            J[look[(d, 1-p, m, n)], c] += (m/2 if p == 0 else -m/2)
        step = 2 if d == 0 else 1
        for j in range(n-step, -1, -step):
            J[look[(d, p, m, j)], c] += 2*MU*n
    return J


def colnorm(M):
    return float(np.max(np.sum(np.abs(M.real)+np.abs(M.imag), axis=0)))


def main():
    A, addr = read(ROOT/'build/spectrum/rt-A-12x36.dat')
    n = len(addr)
    Q = np.linalg.inv(free(addr))
    w = np.array([ETA[d]*(2-(m == 0))*(2-(nn == 0))*(65/64)**m*(5/4)**nn for d, p, m, nn in addr])
    Wm, Wi = w[:, None], 1/w[None, :]
    core = np.array([m < 6 and nn < 18 for d, p, m, nn in addr])
    C, O = np.where(core)[0], np.where(~core)[0]
    out = dict(note='DIAGNOSTICO float; nao e prova', dim=n, core=int(len(C)), shell=int(len(O)), centers={})
    for s0 in (1.0, 0.125+0.25j, 0.733):
        H = (Wm*((A+s0*np.eye(n))@Q)*Wi).astype(complex)
        I = np.eye(n)
        X = H-I
        res = {}
        # P0, P1
        for name, shell in (('P0_casca_identidade', np.eye(len(O))), ('P1_RT_primeira_ordem', np.eye(len(O))-X[np.ix_(O, O)])):
            P = np.zeros((n, n), complex)
            P[np.ix_(C, C)] = np.linalg.inv(H[np.ix_(C, C)])
            P[np.ix_(O, O)] = shell
            D = I-H@P
            res[name] = dict(defect=colnorm(D), rho=float(np.max(np.abs(np.linalg.eigvals(D)))))
            if name.startswith('P1'):
                P1, D1 = P, D
        # P3k: refinamento de P1
        Dk = np.eye(n, dtype=complex)
        for k in range(1, 9):
            Dk = Dk@D1
            if k in (2, 4, 8):
                res[f'P3_refino_k{k}'] = dict(defect=colnorm(Dk))
        # P2k: Schur com Neumann amortecida na casca
        HOO = H[np.ix_(O, O)]
        for omega in (1.0, 0.5):
            E = np.eye(len(O))-omega*HOO
            res[f'casca_rho_I-omegaH_OO_w{omega}'] = float(np.max(np.abs(np.linalg.eigvals(E))))
            V = np.zeros_like(E)
            Pw = np.eye(len(O), dtype=complex)
            for k in range(1, 65):
                V = V+omega*Pw
                Pw = Pw@E
                if k in (8, 16, 32, 64):
                    S = H[np.ix_(C, C)]-H[np.ix_(C, O)]@V@H[np.ix_(O, C)]
                    Si = np.linalg.inv(S)
                    # P = LDU^-1 aproximada: [[Si, -Si H_CO V],[-V H_OC Si, V + V H_OC Si H_CO V]]
                    P = np.zeros((n, n), complex)
                    P[np.ix_(C, C)] = Si
                    P[np.ix_(C, O)] = -Si@H[np.ix_(C, O)]@V
                    P[np.ix_(O, C)] = -V@H[np.ix_(O, C)]@Si
                    P[np.ix_(O, O)] = V+V@H[np.ix_(O, C)]@Si@H[np.ix_(C, O)]@V
                    D = I-H@P
                    res[f'P2_schur_w{omega}_k{k}'] = dict(defect=colnorm(D), shell_defect=colnorm(Pw),
                                                          norm_P=colnorm(P))
        res['norm_H_inv'] = colnorm(np.linalg.inv(H))
        res['max_col_X'] = colnorm(X)
        out['centers'][str(s0)] = res
        print(str(s0), json.dumps(res), flush=True)
    # tempo: aplicar H denso a uma coluna nao representa o custo matriz-livre; so referencia
    t = time.time()
    for _ in range(200):
        H[:, 5]@H[5, :]
    out['dense_ref_seconds_per_200_dotproducts'] = time.time()-t
    (ROOT/'build/spectrum/g-parametrix-prototype.json').write_text(json.dumps(out, indent=2)+'\n')


if __name__ == '__main__':
    main()
