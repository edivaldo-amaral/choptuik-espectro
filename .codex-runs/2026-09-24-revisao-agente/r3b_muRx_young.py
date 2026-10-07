#!/usr/bin/env python3
"""R3b -- (1) a descida ilimitada de mu R_x nos blocos (Z, far) e (T1, far); (2) Young por modo.

(1) tile_perron.montar cota P_Z E P_far por nV (nB descida(NZ + BAND_N + 1, Nc) + dB Qf) e
P_T1 E P_far por nB descida(NP + BAND_N + 1, Nc) + dB Qf. Isso pressupoe que o fundo truncado B_t
nao leva n >= NZ + BAND_N + 1 para n < NZ ("o fundo truncado nao ... desce de far ate G+ sem a
folga radial", linha 475). Mas B_rad(d = 0) contem mu R_x: entrada c1 impar n -> c2 par n' para
TODO n' < n (fourier_tail_bound.py, linhas 93-99). Logo falta o termo
    T = P_Z mu R_x P_{c1, n >= J} Q0(s2) P_far        (e o analogo em T1, J = NP + 41).
Aqui: (a) exibe a entrada nao nula; (b) mede T numa truncagem grande; (c) da uma cota fechada
rigorosa de ||T||; (d) refaz o Perron dos 154 tiles gravados com o termo somado.

(2) norma_Bd_young (d != 0) contra a norma verdadeira pesada de B_d numa truncagem grande.
"""
from __future__ import annotations

from fractions import Fraction as Fr
import glob
import json
import math
from pathlib import Path
import sys

import numpy as np

AQUI = Path(__file__).resolve().parent
RAIZ = AQUI.parents[1]
sys.path.insert(0, str(RAIZ/'scripts')); sys.path.insert(0, str(AQUI))
import fourier_tail_bound as ft
import tile_perron as tp
from r1_forma_fechada import forma_fechada, lw

MU, K1, K2 = ft.MU, ft.K1, ft.K2
CW = math.sqrt(1 + K2**-4)


def parte_a():
    print('== (a) entrada de B_rad(d = 0) de c1 n = 401 para c2 n\' = 0')
    g = tp.truncar(ft.modos_campos(), 20, 40)          # o fundo TRUNCADO do certificado
    n1i = np.array([401]); n2i = np.array([], int)
    n1o = np.array([], int); n2o = np.array([0, 2, 78, 400])
    B = ft.B_rad(g, 0, n1i, n2i, n1o, n2o, MU)
    print(f'   B_t[(c2, n\'), (c1, 401)] para n\' = {n2o.tolist()}: {np.real(B[:, 0]).round(6).tolist()}')
    print(f'   (2 mu = {2*MU:.6f}); pesado por sqrt(omega(n\')/omega(401)): '
          f'{(np.abs(B[:, 0])*np.exp(lw(n2o) - lw(np.array([401])))).tolist()}')


def termo_verdade(NZ, J, Nc, s, m, Ntr):
    """||P_{c2, n'<NZ} muR_x P_{c1, J<=n} R^-1 P_{c1, n_c>=Nc}||, pesado, colunas n_c < Ntr."""
    ns = np.arange(1, Ntr, 2)                          # c1: n impar
    R = forma_fechada(ns, s + 1j*m/2)
    Rw = R*np.exp(lw(ns)[:, None] - lw(ns)[None, :])
    lin = ns >= J; col = ns >= Nc
    Rsub = Rw[np.ix_(lin, col)]
    nout = np.arange(0, NZ, 2)                         # c2 par < NZ
    nin = ns[lin]
    k = (nin - 1)//2
    jj = nout//2
    M = np.where(jj[:, None] <= k[None, :], 2*MU*(-1.0)**(k[None, :] - jj[:, None]), 0.0)
    Mw = M*np.exp(lw(nout)[:, None] - lw(nin)[None, :])
    return float(np.linalg.norm(Mw @ Rsub, 2))


def termo_cota(NZ, J, Nc):
    """Cota fechada (Schur) de ||P_{n'<NZ} muR_x P_{n>=J} Q0 P_{n_c>=Nc}||, Re sigma >= 0:
    |muR_x|_w(n', n) <= 2 mu c_w K2^(n'-n);  |Q0|_w(n, n_c) <= 2 c_w K2^(n-n_c)/(mu (n+1)) (n < n_c),
    <= 1/(mu (n_c+1)) (n = n_c). Entrada do produto <= 4 c_w^2 K2^(n'-n_c) (H(n_c+1) + 1).
    Linhas n' par < NZ, colunas n_c impar >= Nc; com h(n_c) = H(n_c+1) + 1 <= ln(n_c+1) + 2:
      coluna n_c: 4 c_w^2 h(n_c) K2^(-n_c) sum_{n'<NZ} K2^n' <= 4 c_w^2 h(n_c) K2^(NZ-n_c)/(K2-1)
      linha n'  : 4 c_w^2 K2^(n') sum_{n_c>=Nc} h(n_c) K2^(-n_c)."""
    c2 = 4*CW*CW
    # colunas: h(n) K2^-n e decrescente para n >= Nc >= 20, max em n_c = Nc
    hc = math.log(Nc + 1) + 2
    col = c2*hc*K2**(NZ - Nc)/(K2 - 1)
    # linhas: soma de h(n) K2^-n para n >= Nc <= h(Nc) K2^-Nc / (1 - K2^-1 (h(Nc+1)/h(Nc)))
    q = K2**-1*(math.log(Nc + 2) + 2)/hc
    lin = c2*K2**(NZ)*hc*K2**(-Nc)/(1 - q)
    return math.sqrt(col*lin)


def parte_b_c():
    print('\n== (b) termo que falta, medido (truncagem n_c < Nc + 700), e (c) cota fechada')
    out = {}
    for rot, (NZ, NP, Nc) in (('caixa pequena', (24, 36, 77)), ('laboratorio', (80, 160, 400))):
        for nome, (N0_, J) in (('Z <- far', (NZ, NZ + 41)), ('T1 <- far', (NP, NP + 41))):
            v = 0.0
            for s in (0.0, 0.125 + 0.125j, 0.25j, 1.0):
                for m in (0, 2, -6, 18):
                    v = max(v, termo_verdade(N0_, J, Nc, s, m, Nc + 700))
            c = termo_cota(N0_, J, Nc)
            out[rot, nome] = (v, c)
            print(f'   {rot:14s} {nome:9s}: medido {v:.3e}   cota fechada {c:.3e}')
    return out


def parte_d(eps):
    """Refaz o Perron dos tiles gravados somando o termo que falta: N0[0][3] += nV eps_Z,
    N0[1][3] += eps_1 (em racionais, com o v gravado)."""
    print('\n== (d) os 154 tiles com o termo somado (racionais exatos, v gravado)')
    eZ, e1 = Fr(eps['laboratorio', 'Z <- far'][1])*2, Fr(eps['laboratorio', 'T1 <- far'][1])*2
    ARRED = 1 + Fr(1, 10**9)
    pior = Fr(0); n = 0; falhas = 0
    for f in sorted(glob.glob(str(RAIZ/'build/cobertura_lab/grande-*.json'))):
        for t in json.load(open(f))['tiles']:
            if not t['ok']:
                continue
            rZ, r1, r2 = Fr(t['rZ']), Fr(t['r1']), Fr(t['r2'])
            N = [[ARRED*(Fr(t['N0'][i][j]) + rZ*Fr(t['NZ'][i][j]) + r1*Fr(t['N1'][i][j]) + r2*Fr(t['N2'][i][j]))
                  for j in range(4)] for i in range(4)]
            N[0][3] += Fr(t['V'])*eZ
            N[1][3] += e1
            v = [Fr(x) for x in t['v']]
            th = max(sum(N[i][j]*v[j] for j in range(4))/v[i] for i in range(4))
            pior = max(pior, th); n += 1
            falhas += th >= 1
    print(f'   {n} tiles; pior theta com o termo (x2 de folga): {float(pior):.10f}; falhas: {falhas}')
    return float(pior)


def parte_young():
    print('\n== (2) norma_Bd_young (d != 0) contra a norma verdadeira pesada de B_d (n < 400)')
    g = ft.modos_campos()
    N = 400
    n1 = np.arange(1, N, 2); n2 = np.arange(N)
    w = np.exp(np.concatenate([lw(n1), lw(n2)]))
    pior = 0.0
    for d in (2, -2, 4, 8, -12, 14, 20, 30, 40):
        X = ft.B_rad(g, d, n1, n2, n1, n2, MU)
        Xw = (w[:, None]*X)/w[None, :]
        v = float(np.linalg.norm(Xw, 2)); c = ft.norma_Bd_young(g, d, MU)
        pior = max(pior, v/c)
        print(f'   d = {d:+d}: verdade {v:.5e}  Young {c:.5e}  razao {v/c:.4f}')
    # o termo de mu R_x em d = 0 (retirado do certificado): conferir so por curiosidade
    X = ft.B_rad(g, 0, n1, n2, n1, n2, MU); Xw = (w[:, None]*X)/w[None, :]
    print(f'   d =  0 (com mu R_x): verdade {np.linalg.norm(Xw, 2):.5f}  Young {ft.norma_Bd_young(g, 0, MU):.5f}'
          f'  (nao usado no certificado atual)')
    return pior


if __name__ == '__main__':
    parte_a()
    eps = parte_b_c()
    parte_d(eps)
    parte_young()
