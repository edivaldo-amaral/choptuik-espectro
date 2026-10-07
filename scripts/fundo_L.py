#!/usr/bin/env python3
"""Perturbacao global do fundo de L (S3b), para o modo rigoroso do certificado de L.

L(s) = J_mu + B + s, B = mu S Gamma1 + 2 S Xi Gamma2(omega, .), na norma do toro de L
((65/64, 5/4), m com sinal; peso radial omega(n) = k2^2n + k2^-2n). O bloco (m + d <- m) de B e,
por componente de saida cout e entrada (cin, paridade de n), a multiplicacao radial
2 fator mult(a) pela funcao-coeficiente a = sum coef (P?) campo_d (tabela (23),
signed_operator_L.TABELA), mais mu S Gamma1 em d = 0. Cotas:

  Young radial: ||mult(a)|| <= l1(a) = |a_0| + 2 sum_n |a_n| k2^n; o bloco d tem norma
  <= ||M(d)||_2, M(d)[cout, (cin, par)] = 2 fator l1(a); mu S Gamma1 (d = 0) soma
  mu |fac| norma_divxi(k2) (Schur, closed_form). Fourier: ||B|| <= sum_d k1^|d| ||M(d)||.

  dB_L = (a) arredondamento da montagem: gamma_12 sum_d k1^|d| ||M_abs(d)||, com M_abs dos
             coeficientes em valor absoluto (combinacao <= 4 termos, mult, mu S Gamma1, pesos);
         (b) resto da banda: sum_{BANDA < |d| <= 40} k1^|d| ||M(d)|| (os dados param em |d| = 40
             e n = 100: nao ha outro truncamento);
         (c) erro do fundo de RT: 40 eps_omega (eps_omega = 2^-25 + 2^-277 por componente, bola
             publicada de RT nos raios (65/64, 5/4); RT dao a correcao <= 6 eps_omega; 40 como em K,
             auditoria C4 de S3a).
  eps_mu_L(s): ||delta_mu (dJ/dmu + S Gamma1) Q0(s)||, |delta_mu| <= DELTA_MU (tile_perron).
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import closed_form as cf
import signed_operator_L as sl
from fourier_tail_bound import refl
from tile_perron import DELTA_MU
from verified_norms import gamma_real

K1, K2, MU = sl.K1, sl.K2, sl.MU
EPS_OMEGA = 2.0**-25 + 2.0**-277
EPS_FUNDO_L = 40*EPS_OMEGA
# S4: com RefB explicito (RefAplusB.dat), o erro do fundo verdadeiro em B e ||B_L(RefB)|| + Corr,
# calculado por fundo_refB.py; a variavel de ambiente o substitui (e fica registrada nas saidas)
import os as _os
if _os.environ.get('CHOPTUIK_EPS_FUNDO_L'):
    EPS_FUNDO_L = float(_os.environ['CHOPTUIK_EPS_FUNDO_L'])
# entradas (cin, paridade) na ordem das colunas de M
ENTRADAS = [(0, None), (1, 'e'), (1, 'o'), (2, 'e'), (2, 'o'), (3, 'e'), (3, 'o')]


def l1(a):
    n = np.arange(len(a))
    return float(np.sum(np.abs(a)*np.where(n == 0, 1.0, 2*K2**n)))


def combo(g, termos, d, absoluto=False):
    out = np.zeros(101, complex)
    for c, j, p in termos:
        a = sl.modo(g, j, d)
        a = refl(a) if p else a
        out += (abs(c)*np.abs(a)) if absoluto else c*a
    return out


def matriz_M(g, d, absoluto=False, mu=MU):
    """M(d)[cout, (cin, par)]: cotas de Young dos blocos radiais de B_d."""
    M = np.zeros((4, len(ENTRADAS)))
    for k, (cin, par) in enumerate(ENTRADAS):
        chave = ('w1', None) if cin == 0 else (f'w{cin+1}', par)
        fator, saidas = sl.TABELA[chave]
        for cout in range(4):
            if saidas[cout]:
                M[cout, k] = 2*fator*l1(combo(g, saidas[cout], d, absoluto))
    if d == 0 and not absoluto:                       # mu S Gamma1 (so desce, Schur)
        nd = cf.norma_divxi(K2)
        M[1, 0] += mu*1*nd; M[2, 0] += mu*1*nd        # cin 0 -> cout 1 (+1), 2 (-1)
        M[1, 2] += mu*2*nd                            # cin 1, n impar -> cout 1
        M[3, 6] += mu*2*nd                            # cin 3, n impar -> cout 3
    return M


def norma_Bd(g, d, absoluto=False):
    return float(np.linalg.norm(matriz_M(g, d, absoluto), 2))*1.0001


def dB_L(g, banda):
    arred = gamma_real(12)*sum(K1**abs(d)*norma_Bd(g, d, True) for d in range(-sl.BAND_M, sl.BAND_M + 1))
    resto = sum(K1**abs(d)*norma_Bd(g, d) for d in range(-sl.BAND_M, sl.BAND_M + 1) if abs(d) > banda)
    return dict(arred=arred, resto=resto, rt=EPS_FUNDO_L, total=(arred + resto + EPS_FUNDO_L)*1.0001)


def eps_mu_L(s):
    """||delta_mu (J1 + S Gamma1) Q0(s)||: modo a modo, J1 R_m^-1 = (I - sigma_m R_m^-1)/MU e
    ||R_m^-1|| <= c/max(MU, |m|/2 - |Im s|) (forma fechada, Re s >= 0), c = 1 + 2 c_w/(k2 - 1);
    S Gamma1 = (mu S Gamma1)/mu: coeficientes |fac| <= 2, norma <= 2 norma_divxi (Schur, somando
    as duas saidas de cin 0)."""
    cw = math.sqrt(1 + K2**-4)
    c = 1 + 2*cw/(K2 - 1)
    ng1 = 2*cf.norma_divxi(K2)
    pior = 0.0
    for m in range(0, 4001):
        h = max(MU, m/2 - abs(s.imag))
        nR = c/h
        sig = abs(s.real) + m/2 + abs(s.imag)
        pior = max(pior, (1 + sig*nR)/MU + ng1*nR)
    pior = max(pior, (1 + c*(1 + 1e-3))/MU)
    return DELTA_MU*pior*1.0001


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--banda', type=int, default=16)
    ap.add_argument('--s', type=complex, default=0.4010247330)
    a = ap.parse_args()
    g = sl.campos()
    tot = sum(K1**abs(d)*norma_Bd(g, d) for d in range(-sl.BAND_M, sl.BAND_M + 1))
    print(f'||B_L|| por Young (Fourier x radial): {tot:.4f}   (simbolo em Arb: 3,964043)')
    db = dB_L(g, a.banda)
    print(f'dB_L (banda {a.banda}): arredondamento {db["arred"]:.2e}, resto da banda {db["resto"]:.2e},'
          f' erro de RT {db["rt"]:.2e};  total {db["total"]:.3e}')
    print(f'eps_mu_L({a.s}) = {eps_mu_L(a.s):.3e}')


if __name__ == '__main__':
    main()
