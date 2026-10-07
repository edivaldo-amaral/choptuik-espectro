#!/usr/bin/env python3
"""Checagem INDEPENDENTE (auditoria de 16/09) da matriz exportada rt-A-NxM.dat.

Reconstroi, a partir SOMENTE das formulas do TeX de RT (arXiv:1203.3766v1),
o operador selecionado linearizado no fundo RefA,
    L(0) = O_mu + mu S Gamma1 + 2 S Xi Gamma2(RefA, .),
nas coordenadas reais (d, parte, m, n), e compara TODAS as entradas com a
matriz exportada (igualdade racional exata). Nao importa modulos do
repositorio; usa a leitura de RefA de independent_component_beta.py
(tambem independente) e o O_mu de independent_crown_tiles.py.

Convencoes (RT secao 3): campo v_{mn}, m,n em Z, v_{m,-n}=v_{mn},
v_{-m,-n}=conj(v_{mn}); convolucao em Z^2; (Pv)_{mn}=(-1)^n v_{mn};
Xi^-1_reg = 2 Band_{1,1} BandI_{2,1}, isto e (R v)_j = 2 sum_L (-1)^L v_{j+1+2L}.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import json
from pathlib import Path

from independent_component_beta import parse_refa, selected_coefficients
from independent_crown_tiles import read_matrix, free_matrix, MU

ROOT = Path(__file__).resolve().parents[1]


def expand(field_int, scale):
    """dict (m,n)>=0 -> (re,im) inteiros  ==>  dict Z^2 -> complex Fraction."""
    out = {}
    for (m, n), (a, b) in field_int.items():
        z = (Q(a)*scale, Q(b)*scale)
        for nn in {n, -n}:
            out[m, nn] = z
            out[-m, -nn] = (z[0], -z[1])
    return out


def unit_z2(m, n, part):
    z = (Q(1), Q(0)) if part == 0 else (Q(0), Q(1))
    out = {}
    for nn in {n, -n}:
        out[m, nn] = z
        out[-m, -nn] = (z[0], -z[1])
    return out


def conv(c, e):
    out = {}
    for (m2, n2), (a2, b2) in e.items():
        for (m1, n1), (a1, b1) in c.items():
            k = (m1+m2, n1+n2)
            x = out.get(k, (Q(0), Q(0)))
            out[k] = (x[0]+a1*a2-b1*b2, x[1]+a1*b2+b1*a2)
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('matrix', type=Path, nargs='?', default=ROOT/'build/spectrum/rt-A-6x18.dat')
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    addresses, A = read_matrix(args.matrix)
    look = {a: i for i, a in enumerate(addresses)}
    dim = len(addresses)
    J = free_matrix(addresses, MU)
    fields, scale = parse_refa()
    tab = selected_coefficients(fields)
    # Campos coeficiente expandidos em Z^2, com o fator 2 de 2 S Xi Gamma2.
    coef = {}
    for key, (c, vec) in tab.items():
        coef[key] = [({k: (2*c*a, 2*c*b) for k, (a, b) in expand(f, scale).items()} if f else {})
                     for f in vec]
    mismatches, max_abs = 0, Q(0)
    worst = None
    for col, (d, part, m, n) in enumerate(addresses):
        e = unit_z2(m, n, part)
        col_out = {}   # (i, M, N) -> complexo (so M>=0, N>=0)
        keys = [('w1', None)] if d == 0 else [(f'w{d+1}', 'e' if n % 2 == 0 else 'o')]
        for key in keys:
            for i, cf in enumerate(coef[key]):
                if not cf:
                    continue
                # (cf * e)_{MN} = sum_{(m2,n2) in supp e} cf_{M-m2,N-n2} e_{m2 n2},
                # avaliado apenas nas linhas do box (equivale a conv(cf, e)).
                for (_, _, M, N) in addresses:
                    if (i, M, N) in col_out:
                        continue
                    re_, im_ = Q(0), Q(0)
                    for (m2, n2), (a2, b2) in e.items():
                        a1, b1 = cf.get((M-m2, N-n2), (0, 0))
                        re_ += a1*a2-b1*b2
                        im_ += a1*b2+b1*a2
                    col_out[i, M, N] = (re_, im_)
        # mu S Gamma1 e = mu (0, R(e1+(1-P)e2), -R e1, R(1-P)e4); R so baixa n.
        src = {}
        z = (Q(1), Q(0)) if part == 0 else (Q(0), Q(1))
        if d == 0:
            src = {1: 1, 2: -1}
        elif d == 1 and n % 2 == 1:
            src = {1: 2}
        elif d == 3 and n % 2 == 1:
            src = {3: 2}
        for i, factor in src.items():
            L = 0
            while n-1-2*L >= 0:
                j = n-1-2*L
                s = MU*factor*2*(-1)**L
                x = col_out.get((i, m, j), (Q(0), Q(0)))
                col_out[i, m, j] = (x[0]+s*z[0], x[1]+s*z[1])
                L += 1
        for row, (i, p, M, N) in enumerate(addresses):
            val = col_out.get((i, M, N), (Q(0), Q(0)))[p]
            exported = A[row][col]-J[row][col]
            if val != exported:
                mismatches += 1
                diff = abs(val-exported)
                if diff > max_abs:
                    max_abs, worst = diff, (row, col, str(val), str(exported))
    result = dict(matrix=str(args.matrix), dimension=dim, entries_compared=dim*dim,
                  mismatches=mismatches, max_abs_difference=str(max_abs), worst=worst,
                  exact_match=mismatches == 0)
    text = json.dumps(result, indent=2)
    if args.out:
        args.out.write_text(text+'\n')
    print(text)


if __name__ == '__main__':
    main()
