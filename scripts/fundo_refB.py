#!/usr/bin/env python3
"""Cota de Young do operador B_L(RefB) = 2 S Xi Gamma2(RefB, .) (a parte do fundo verdadeiro que
RefA nao tem), a partir de RefB = (RefA + RefB) - RefA EXATO, na norma de L (65/64, 5/4).

O fundo verdadeiro e w* = RefA + RefB + Corr, com RefB explicito (RefAplusB.dat) e ||Corr|| ~
2^-277 na bola de RT. Os certificados de S3b usaram o erro de RT como 40 eps_omega = 1,2e-6
(RefB tratado como desconhecido). Aqui: ||B_L(RefB)|| <= sum_d k1^|d| ||M_RefB(d)||_2, com
M_RefB(d)[cout, (cin, par)] = 2 fator l1(sum |coef| |RefB_campo,d|) (a mesma estrutura de
fundo_L.matriz_M, termos em valor absoluto), l1(a) = |a_0| + 2 sum |a_n| k2^n. O termo de Corr
entra como 1e-60 (generoso para 2^-277 vezes a estrutura). O mu verdadeiro continua em eps_mu_L."""
from __future__ import annotations

import argparse
import math
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fundo_L
import signed_operator_L as sl
from gauge_refB import ler
from fourier_tail_bound import refl


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--refA', default='.cache/rt-1203.3766v1/sourcecode/RefA.dat')
    ap.add_argument('--refAB', default='.cache/rt-1203.3766v1/RefAplusB.dat')
    a = ap.parse_args()
    A, _ = ler(a.refA); AB, _ = ler(a.refAB)
    K1, K2 = sl.K1, sl.K2
    # campos de RefB no formato de campos(): g[j][m] = array de coeficientes (n = 0..num_n-1), m >= 0
    gB = []
    for (hA, pA), (hB, pB) in zip(A, AB):
        eA, eB = hA['TwoExp'], hB['TwoExp']; esc = 2**(eA - eB)
        dd = {}
        for m in range(hB['num_m']):
            v = np.zeros(hB['num_n'], complex)
            for n in range(hB['num_n']):
                xr, xi = int(pB[m*hB['num_n'] + n][0]), int(pB[m*hB['num_n'] + n][1])
                if m < hA['num_m'] and n < hA['num_n']:
                    ar, ai = pA[m*hA['num_n'] + n]; xr -= int(ar)*esc; xi -= int(ai)*esc
                v[n] = complex(math.ldexp(xr, eB) if xr else 0.0, math.ldexp(xi, eB) if xi else 0.0)
            dd[m] = v
        gB.append(dd)
    nn = len(gB[0][0]); ns = np.arange(nn)
    peso = np.where(ns == 0, 1.0, 2*K2**ns)
    def modo(j, d):
        a_ = gB[j].get(abs(d))
        if a_ is None:
            return np.zeros(nn, complex)
        return a_.conj() if d < 0 else a_
    def l1_abs_combo(termos, d):
        out = np.zeros(nn)
        for c, j, p in termos:
            v = np.abs(modo(j, d))                      # |refl(a)| = |a|
            out += abs(c)*v
        return float(np.sum(out*peso))
    Mmax = max(len(gB[j]) for j in range(4)) - 1
    tot = 0.0
    for d in range(-Mmax, Mmax + 1):
        M = np.zeros((4, len(fundo_L.ENTRADAS)))
        for k, (cin, par) in enumerate(fundo_L.ENTRADAS):
            chave = ('w1', None) if cin == 0 else (f'w{cin+1}', par)
            fator, saidas = sl.TABELA[chave]
            for cout in range(4):
                if saidas[cout]:
                    M[cout, k] = 2*fator*l1_abs_combo(saidas[cout], d)
        nd = float(np.linalg.norm(M, 2))
        tot += K1**abs(d)*nd
    rt = tot*1.001 + 1e-60
    print(f'||B_L(RefB)|| <= {rt:.4e}  (Young, |d| <= {Mmax}); contra 40 eps_omega = {fundo_L.EPS_FUNDO_L:.4e}')


if __name__ == '__main__':
    main()
