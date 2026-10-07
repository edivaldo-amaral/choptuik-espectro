#!/usr/bin/env python3
"""Recertificacao INDEPENDENTE (auditoria 16/09) do disco sharp finito.

Disco |s-733/1000|<=1/32, pesos (129/128, 9/8), sem pesos de componentes.
Nao importa o certificador auditado. Constroi K_hat(s0)=W(K+s0 I)W^-1 com
erro de arredondamento de coluna calculado exatamente, V=inv em float (so
proposta), residuos inteiros exatos e orcamento de fundo derivado das
constantes RT:
   ||K_true,N - K_ref,N|| <= eps_mu (||J'_N|| + ||R||) + 3 eps_omega,
||R||=2 k^-1/(1-k^-2)=144/17 (Xi^-1_reg), 3 = Lipschitz de Gamma2^sharp (RT, secao 3),
eps_mu=43*2^-40+2^-277, eps_omega=2^-25+2^-277.
Modo --left : ||K_N(s)^-1|| <= ||V||/(1-q),  q=||I-V K_hat||+||V||(E+bg+r).
Modo --right: ||J_ref,N K_N(s)^-1|| <= ||J V||/(1-q'), q'=||I-K_hat V||+||V||(E+bg+r)
(bound da inversa PRECONDICIONADA H_FF^-1=J_ref,F K_FF^-1).
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import json
from operator import mul
from pathlib import Path

import numpy as np

from independent_crown_tiles import free_matrix, MU, rnd
from independent_sharp_matrix import read_sharp

ROOT = Path(__file__).resolve().parents[1]
G_BITS, V_BITS = 64, 56


def matmul(X, Y):
    cols = list(zip(*Y))
    return [[sum(map(mul, row, col)) for col in cols] for row in X]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('matrix', type=Path)
    parser.add_argument('--mode', choices=('left', 'right'), default='left')
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    s0, radius = Q(733, 1000), Q(1, 32)
    addresses, K = read_sharp(args.matrix)
    n = len(addresses)
    w = [(2-(m == 0))*(2-(nn == 0))*Q(129, 128)**m*Q(9, 8)**nn for d, p, m, nn in addresses]
    g = 1 << G_BITS
    Kw = [[(K[i][j]+(s0 if i == j else 0))*w[i]/w[j] for j in range(n)] for i in range(n)]
    Kr = [[rnd(g*x) for x in row] for row in Kw]
    EK = max(sum(abs(Kw[i][j]-Q(Kr[i][j], g)) for i in range(n)) for j in range(n))
    J1 = free_matrix(addresses, Q(1))
    J0 = free_matrix(addresses, Q(0))
    jprime = max(sum(abs(J1[i][j]-J0[i][j])*w[i] for i in range(n))/w[j] for j in range(n))
    eps_mu = Q(43, 2**40)+Q(1, 2**277)
    eps_om = Q(1, 2**25)+Q(1, 2**277)
    bg = eps_mu*(jprime+Q(144, 17))+3*eps_om
    V = np.linalg.inv(np.array(Kr, dtype=float)/g)
    gv = 1 << V_BITS
    Vr = [[round(float(x)*gv) for x in row] for row in V]
    vnorm = Q(max(sum(abs(Vr[i][j]) for i in range(n)) for j in range(n)), gv)
    prod = matmul(Vr, Kr) if args.mode == 'left' else matmul(Kr, Vr)
    den = gv*g
    res = Q(max(sum(abs(prod[i][j]-(den if i == j else 0)) for i in range(n)) for j in range(n)), den)
    q = res+vnorm*(EK+bg+radius)
    out = dict(matrix=str(args.matrix), mode=args.mode, dimension=n, residual=float(res),
               rounding=float(EK), background=str(bg), background_decimal=float(bg),
               jprime_norm=float(jprime), V_norm=float(vnorm),
               uniform_defect=str(q), uniform_defect_decimal=float(q), accepted=q < 1)
    if q < 1 and args.mode == 'left':
        out['inverse_bound_decimal'] = float(vnorm/(1-q))
        out['inverse_bound'] = str(vnorm/(1-q))
    if q < 1 and args.mode == 'right':
        Jw = free_matrix(addresses, MU)
        jcols = [[(i, Jw[i][j]*w[i]/w[j]) for i in range(n) if Jw[i][j]] for j in range(n)]
        # ||J V||: coluna k de J V = sum_j J[:,j] V[j,k]
        best = Q(0)
        for k in range(n):
            acc = {}
            for j in range(n):
                v = Vr[j][k]
                if v:
                    for i, x in jcols[j]:
                        acc[i] = acc.get(i, Q(0))+x*v
            best = max(best, sum((abs(x) for x in acc.values()), Q(0))/gv)
        out['JV_norm'] = float(best)
        out['preconditioned_inverse_bound'] = str(best/(1-q))
        out['preconditioned_inverse_bound_decimal'] = float(best/(1-q))
    text = json.dumps(out, indent=2)
    if args.out:
        args.out.write_text(text+'\n')
    print(text)


if __name__ == '__main__':
    main()
