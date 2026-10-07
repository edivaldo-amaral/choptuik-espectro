#!/usr/bin/env python3
"""S4, identificacao: ||(Z + 1) RefB|| e ||RefB|| nos pesos de L, com RefB = (RefA + RefB) - RefA
EXATO (diferenca de diadicos em inteiros), dos arquivos de RT. Z = mu^-1 d_tau + xi d_xi.

Norma de L (m com sinal): sum_m kappa1^(2m) omega2(n) |x_mn|^2; para campos reais (x_-m = conj x_m)
cada m > 0 entra com (kappa1^(2m) + kappa1^(-2m)). Coeficientes simetricos de Chebyshev; xi d_xi
por signed_constraints.xidxi (n na diagonal, 2n nas linhas n - 2, n - 4, ...), aplicado por soma
acumulada (O(n) por modo). Resultado em ponto flutuante com o valor calculado e uma margem de
arredondamento relativa 1e-10 (somas de ~10^6 termos positivos)."""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import math
from pathlib import Path
import re
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import signed_operator_L as sl

PAR = re.compile(r'\((-?\d+),(-?\d+)\)')


def ler(path):
    """[(TwoExp, off_m, off_n, num_m, num_n, lista de pares (re, im) inteiros)] por campo, e o mu."""
    txt = Path(path).read_text()
    mu = None
    mm = re.search(r'mu_MultiField\s+DyadicQ\s+TwoExp\s+(-?\d+)\s+coeff\s+(-?\d+)', txt)
    if mm:
        mu = Q(int(mm[2]))*Q(2)**int(mm[1])
    blocos = re.split(r'(?m)^Field\s*$', txt)[1:]
    out = []
    for b in blocos:
        h = {k: int(re.search(rf'(?m)^{k}\s+(-?\d+)\s*$', b)[1]) for k in ('TwoExp', 'off_m', 'off_n', 'num_m', 'num_n')}
        corpo = b[b.index('('):]
        pares = PAR.findall(corpo)
        assert len(pares) == h['num_m']*h['num_n'], (len(pares), h)
        out.append((h, pares))
    return out, mu


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--refA', default='.cache/rt-1203.3766v1/sourcecode/RefA.dat')
    ap.add_argument('--refAB', default='.cache/rt-1203.3766v1/RefAplusB.dat')
    a = ap.parse_args()
    A, muA = ler(a.refA)
    AB, muAB = ler(a.refAB)
    print(f'mu_RefA = {float(muA):.17f}; mu_RefA + mu_RefB = {float(muAB):.17f}; |mu_RefB| = {float(abs(muAB - muA)):.3e}', flush=True)
    mu = float(muA)
    K1, K2 = sl.K1, sl.K2
    tot_B = 0.0; tot_ZB = 0.0; tot_tB = 0.0
    for f, ((hA, pA), (hB, pB)) in enumerate(zip(A, AB)):
        eA, eB = hA['TwoExp'], hB['TwoExp']
        assert hA['off_m'] == hB['off_m'] == 0 and hA['off_n'] == hB['off_n'] == 0 and eB <= eA
        nmB, nnB = hB['num_m'], hB['num_n']
        escala = 2**(eA - eB)                            # RefA em unidades de 2^eB
        for m in range(nmB):
            re_ = np.zeros(nnB); im_ = np.zeros(nnB)
            for n in range(nnB):
                xr, xi = int(pB[m*nnB + n][0]), int(pB[m*nnB + n][1])
                if m < hA['num_m'] and n < hA['num_n']:
                    ar, ai = pA[m*hA['num_n'] + n]
                    xr -= int(ar)*escala; xi -= int(ai)*escala
                re_[n] = math.ldexp(xr, eB) if xr else 0.0
                im_[n] = math.ldexp(xi, eB) if xi else 0.0
            v = re_ + 1j*im_                              # coeficiente do modo m (convencao de campos())
            ns = np.arange(nnB)
            # (Z + 1) v = (1 + i m/(2 mu)) v + xi d_xi v; xi d_xi: n v_n + 2 sum_{k > n, k = n mod 2} k v_k
            s_par = np.zeros(nnB, complex)
            for par in (0, 1):
                idx = ns[ns % 2 == par]
                kv = (idx*v[idx])[::-1]
                acum = np.concatenate([[0], np.cumsum(kv)[:-1]])[::-1]   # soma dos k > n
                s_par[idx] = 2*acum
            zv = (1 + 1j*m/(2*mu))*v + ns*v + s_par
            w2 = np.where(ns == 0, 1.0, K2**(2.0*ns) + K2**(-2.0*ns))
            wm = 1.0 if m == 0 else (K1**(2.0*m) + K1**(-2.0*m))
            tot_B += wm*float(np.sum(w2*np.abs(v)**2))
            tot_ZB += wm*float(np.sum(w2*np.abs(zv)**2))
            tot_tB += wm*(m/2)**2*float(np.sum(w2*np.abs(v)**2))     # d_tau -> i m/2
        print(f'  campo {f + 1}: acumulado ||RefB|| {math.sqrt(tot_B):.3e}, ||(Z+1) RefB|| {math.sqrt(tot_ZB):.3e}', flush=True)
    nB = math.sqrt(tot_B)*(1 + 1e-10); nZB = math.sqrt(tot_ZB)*(1 + 1e-10); ntB = math.sqrt(tot_tB)*(1 + 1e-10)
    print(f'||RefB||_L = {nB:.4e};  ||(Z + 1) RefB||_L = {nZB:.4e};  ||d_tau RefB||_L = {ntB:.4e}  (pesos de L (65/64, 5/4), m com sinal)')


if __name__ == '__main__':
    main()
