#!/usr/bin/env python3
"""Checagem INDEPENDENTE (auditoria de 16/09) da matriz sharp exportada.

O exportador C (scripts/rt_sharp_matrix.c) foi escrito na janela nao auditada.
Aqui reconstruimos, SO a partir do TeX de RT (arXiv:1203.3766v1, secao 3,
"The sharp system", eq. (26)), o operador linear
    K(0) = O^sharp_mu + mu Xi^-1_reg Gamma1^sharp + Gamma2^sharp(RefA, .),
com Gamma1^sharp w = (0, w1), O^sharp = diag(O^B, O^A), e comparamos TODAS as
entradas com rt-sharp-NxM.dat (igualdade racional exata).
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import json
from pathlib import Path

from independent_component_beta import parse_refa, sharp_coefficients
from independent_crown_tiles import free_matrix, MU
from independent_selected_matrix import expand, unit_z2

ROOT = Path(__file__).resolve().parents[1]


def read_sharp(path):
    lines = path.read_text().splitlines()
    assert lines[0] == 'CHOPTUIK_RT_SHARP_MATRIX_V1'
    dim = int(lines[1].split()[1])
    a0 = lines.index('BEGIN_ADDRESSES')+1
    addresses = []
    for k in range(dim):
        idx, d, part, m, n = map(int, lines[a0+k].split())
        assert idx == k and d in (0, 1)
        addresses.append((d, part, m, n))
    e0 = lines.index('BEGIN_ENTRIES')+1
    A = [[Q(0)]*dim for _ in range(dim)]
    for line in lines[e0:]:
        f = line.split()
        if f:
            A[int(f[0])][int(f[1])] += Q(int(f[2]))*Q(2)**(int(f[3]) if len(f) > 3 else 0)
    return addresses, A


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('matrix', type=Path, nargs='?', default=ROOT/'build/spectrum/rt-sharp-8x24.dat')
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    addresses, A = read_sharp(args.matrix)
    dim = len(addresses)
    J = free_matrix(addresses, MU)   # d=0 -> O^B (passo 2), d=1 -> O^A
    fields, scale = parse_refa()
    coef = {key: [{k: (c*a, c*b) for k, (a, b) in expand(f, scale).items()} for f in vec]
            for key, (c, vec) in sharp_coefficients(fields).items()}
    mismatches, worst = 0, None
    for col, (d, part, m, n) in enumerate(addresses):
        e = unit_z2(m, n, part)
        key = ('w1', None) if d == 0 else ('w2', 'e' if n % 2 == 0 else 'o')
        out = {}
        for i, cf in enumerate(coef[key]):
            for (_, _, M, N) in addresses:
                if (i, M, N) in out:
                    continue
                re_, im_ = Q(0), Q(0)
                for (m2, n2), (a2, b2) in e.items():
                    a1, b1 = cf.get((M-m2, N-n2), (0, 0))
                    re_ += a1*a2-b1*b2
                    im_ += a1*b2+b1*a2
                out[i, M, N] = (re_, im_)
        if d == 0:  # mu R w1 na segunda componente
            z = (Q(1), Q(0)) if part == 0 else (Q(0), Q(1))
            L = 0
            while n-1-2*L >= 0:
                j, s = n-1-2*L, MU*2*(-1)**L
                x = out.get((1, m, j), (Q(0), Q(0)))
                out[1, m, j] = (x[0]+s*z[0], x[1]+s*z[1])
                L += 1
        for row, (i, p, M, N) in enumerate(addresses):
            val = out.get((i, M, N), (Q(0), Q(0)))[p]
            if val != A[row][col]-J[row][col]:
                mismatches += 1
                if worst is None:
                    worst = (row, col, str(val), str(A[row][col]-J[row][col]))
    result = dict(matrix=str(args.matrix), dimension=dim, entries_compared=dim*dim,
                  mismatches=mismatches, first_mismatch=worst, exact_match=mismatches == 0)
    text = json.dumps(result, indent=2)
    if args.out:
        args.out.write_text(text+'\n')
    print(text)


if __name__ == '__main__':
    main()
