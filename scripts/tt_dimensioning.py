#!/usr/bin/env python3
"""Dimensionamento fechado: majorante de tres blocos com os dados reais.

Fecha os dois buracos que TT_CEILING deixou: os acoplamentos P<->T passam a
usar a cauda medida do fundo (background-couplings.json) em vez de serem
desprezados, e ||V|| vira parametro explicito em vez de 1.

Blocos P (nucleo), S (coroa), T (exterior). O criterio de FREQUENCY_BLOCKS
(existe r>0 com sum_i r_i e_ij < r_j por coluna) equivale a rho(E)<1 para E
nao negativa, entao o teste e o raio espectral da matriz 3x3.

Duas escolhas livres que ninguem tinha feito:
  - a coroa precisa de largura ADITIVA ~(42,102), nao de um fator sobre o
    nucleo: RefA tem suporte finito e alem dessa separacao so resta o
    epsilon publicado;
  - a razao de aspecto otima para t=max(15/M,81/(N+1)) e N+1=5,4M, nao a
    3,59 de 107x384.

DIAGNOSTICO de dimensionamento. Nada aqui e certificado.
"""
from __future__ import annotations

import argparse
import json
from fractions import Fraction as Q

import numpy as np

RHO = Q(5, 4)
MU_MIN = Q(1, 6)
Q0_NORM = Q(18)/MU_MIN                       # 108, em z0=0
SEP_M, SEP_N = 42, 102                       # separacao de cauda nula de RefA
E_TP = Q(1931, 10**8)                        # 6(tail+eps)*||Q0|| em g=(42,102)
E_PT = Q(809, 10**6)                         # ||P D Q0 T|| medido em separacao comparavel
A_FIN = Q(228, 10**10)                       # parametriz acoplada transportada (so em s=1)
BASE = 107*384


def t_struct(M: int, N: int) -> Q:
    return max(Q(15, M), Q(81, N+1))


def majorante(beta: Q, P: tuple[int, int], G: tuple[int, int], v: Q, a_fin: Q):
    """Matriz majorante 3x3 nas ordens (P,S,T); e[i][j] = bloco de j para i."""
    tP, tG = t_struct(*P), t_struct(*G)
    e_TS = beta*tP
    e_SHT = (beta + RHO)*tG
    e_TT = e_SHT
    e_iT = v*(E_PT + e_SHT)                  # ||E_iT|| <= ||V||(||PHT||+||SHT||)
    return [[a_fin, a_fin, e_iT],
            [a_fin, a_fin, e_iT],
            [E_TP,  e_TS,  e_TT]]


def rho(E) -> float:
    return float(max(abs(np.linalg.eigvals(np.array([[float(x) for x in r] for r in E])))))


def menor_caixa(beta: Q, v: Q, a_fin: Q = A_FIN, aspecto: Q = Q(27, 5), cap: int = 200000):
    """Menor nucleo (logo menor caixa) com rho(E)<1, na razao de aspecto dada."""
    lo, hi = 8, 256
    while hi < cap:
        N = int(aspecto*hi)
        P, G = (hi, N), (hi + SEP_M, N + SEP_N)
        if rho(majorante(beta, P, G, v, a_fin)) < 1:
            break
        lo, hi = hi, hi*2
    else:
        return None
    while lo + 1 < hi:
        mid = (lo + hi)//2
        N = int(aspecto*mid)
        P, G = (mid, N), (mid + SEP_M, N + SEP_N)
        if rho(majorante(beta, P, G, v, a_fin)) < 1:
            hi = mid
        else:
            lo = mid
    N = int(aspecto*hi)
    P, G = (hi, N), (hi + SEP_M, N + SEP_N)
    E = majorante(beta, P, G, v, a_fin)
    return dict(core=P, box=G, rho=rho(E), dofs=G[0]*G[1],
                fator=G[0]*G[1]/BASE, E=[[str(x) for x in r] for r in E])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', default=None)
    args = ap.parse_args()

    betas = {'23 (uniforme)': Q(23),
             '10,372 (piso da reponderacao, hoje)': Q(103717, 10000),
             '8,378 (2b completo)': Q(8378, 1000),
             '5,000': Q(5),
             '3,000': Q(3)}
    vs = [Q(1), Q(2), Q(5), Q(10)]

    print(f"{'beta':<36}" + ''.join(f"{'||V||='+str(float(v)):>14}" for v in vs))
    print('-'*(36+14*len(vs)))
    linhas = []
    for nome, beta in betas.items():
        cels = []
        for v in vs:
            r = menor_caixa(beta, v)
            if r is None:
                cels.append('       —')
                continue
            cels.append(f"{r['fator']:>13.1f}x")
            linhas.append(dict(beta=str(beta), beta_nome=nome, V=str(v),
                               nucleo=f"{r['core'][0]}x{r['core'][1]}",
                               caixa=f"{r['box'][0]}x{r['box'][1]}",
                               rho=round(r['rho'], 6), fator=round(r['fator'], 2)))
        print(f"{nome:<36}" + ''.join(cels))

    print('\ndetalhe da coluna ||V||=2 (caixas):')
    for l in linhas:
        if l['V'] == '2':
            print(f"  beta={l['beta_nome']:<36} nucleo {l['nucleo']:>12}"
                  f"  caixa {l['caixa']:>13}  {l['fator']:>7.1f}x  rho={l['rho']}")

    rep = dict(nota='DIAGNOSTICO; e_TP, e_PT e a_fin vem de medidas existentes, V e parametro',
               data='2026-09-20', separacao=[SEP_M, SEP_N], e_TP=str(E_TP),
               e_PT=str(E_PT), a_fin=str(A_FIN), aspecto='27/5', linhas=linhas)
    if args.out:
        with open(args.out, 'w') as fh:
            json.dump(rep, fh, indent=2)
        print('\nescrito:', args.out)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
