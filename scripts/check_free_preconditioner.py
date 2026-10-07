#!/usr/bin/env python3
"""Audita a precondicao livre e o majorante de transferencia em um centro real.

As contas matriciais sao racionais. A cauda depende das hipoteses publicadas
de docs/ALGEBRAIC_TAIL.md; um ponto nao certifica um contorno nem o fundo exato.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import json
import math
from pathlib import Path

from build_contour import read_spectral_matrix
from build_tile import (approximate_inverse, exact_inverse, matrix_product,
                        one_norm, read_rt_weights)


def read_addresses(path: Path, dimension: int) -> list[tuple[int, int, int, int]]:
    addresses = []
    with path.open(encoding="ascii") as source:
        for line in source:
            if line.strip() == "BEGIN_ADDRESSES":
                break
        for index in range(dimension):
            row, d, part, m, n = map(int, next(source).split())
            if row != index or d not in range(4) or part not in (0, 1) or min(m, n) < 0:
                raise ValueError("enderecos invalidos")
            addresses.append((d, part, m, n))
    if len(set(addresses)) != dimension:
        raise ValueError("enderecos duplicados")
    m_cut = max(a[2] for a in addresses) + 1
    n_cut = max(a[3] for a in addresses) + 1
    expected = {(d, part, m, n)
                for m in range(m_cut) for n in range(n_cut) for d in range(4)
                for part in (0, 1)
                if m % 2 == (1 if d == 3 else 0)
                and (d != 0 or n % 2 == 1) and (m != 0 or part == 0)}
    if set(addresses) != expected:
        raise ValueError("a estimativa exige o retangulo completo de DOFs do setor A")
    return addresses


def free_matrix(addresses: list[tuple[int, int, int, int]], mu: Q) -> list[list[Q]]:
    """J_mu=diag(K_mu,J_mu,J_mu,J_mu), nas coordenadas REAIS de RT.

    Field.c: __Field_J_aux / __Field_K_aux seguidos de Field_V_k_c.
    Diagonal mu(n+1), derivada temporal im/2, e coeficientes 2*mu*n
    abaixo do grau de entrada; K preserva a paridade de n.
    """
    size = len(addresses)
    lookup = {address: index for index, address in enumerate(addresses)}
    result = [[Q(0) for _ in addresses] for _ in addresses]
    for column, (d, part, m, n) in enumerate(addresses):
        result[column][column] = mu * (n + 1)
        if m:
            other = lookup.get((d, 1 - part, m, n))
            if other is None:
                raise ValueError("corte nao preserva o par Fourier real/imaginario")
            result[other][column] = Q(m, 2) * (1 if part == 0 else -1)
        for lower in range(n - (2 if d == 0 else 1), -1, -(2 if d == 0 else 1)):
            row = lookup.get((d, part, m, lower))
            if row is None:
                raise ValueError("corte nao e invariante pelo operador livre")
            result[row][column] = 2 * mu * n
    return result


def integer_matrix(matrix: list[list[Q]]) -> tuple[list[list[int]], int]:
    denominator = math.lcm(*(x.denominator for row in matrix for x in row))
    return [[x.numerator * (denominator // x.denominator) for x in row]
            for row in matrix], denominator


def multiply(left: list[list[Q]], right: list[list[Q]]) -> list[list[Q]]:
    ln, ld = integer_matrix(left)
    rn, rd = integer_matrix(right)
    return [[Q(x, ld * rd) for x in row] for row in matrix_product(ln, rn)]


def input_tail(m: int, n: int, base: Q, mu_min: Q = Q(1, 6)) -> Q:
    if m < 1 or n < 2 or base < 0 or mu_min <= 0:
        raise ValueError("cauda exige M>=1, N>=2, base>=0, mu_min>0")
    return 18 * max(Q(2, m), 1 / (base + mu_min * (n - 1)))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("matrix", type=Path)
    parser.add_argument("--base", type=Q, default=Q(0))
    parser.add_argument("--center", type=Q, default=Q(1))
    parser.add_argument("--mu", type=Q, default=Q(722873400, 2**32),
                        help="mu do fundo exportado; padrao RefA")
    parser.add_argument("--bits", type=int, default=40)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    if args.base < 0 or not Q(1, 6) <= args.mu <= Q(17, 100) or args.bits < 0:
        parser.error("base>=0, 1/6<=mu<=17/100 e bits>=0 sao obrigatorios")
    dimension, raw, exponent = read_spectral_matrix(args.matrix)
    addresses = read_addresses(args.matrix, dimension)
    weights = read_rt_weights(args.matrix, dimension)
    principal = free_matrix(addresses, args.mu)
    for i in range(dimension):
        principal[i][i] += args.base
    inverse_free = exact_inverse(principal)
    pencil = [[Q(x) * Q(2)**exponent for x in row] for row in raw]
    for i in range(dimension):
        pencil[i][i] += args.center
    relative = multiply(pencil, inverse_free)
    weighted = [[x * weights[i] / weights[j] for j, x in enumerate(row)]
                for i, row in enumerate(relative)]
    numerators, denominator = integer_matrix(weighted)
    guess = approximate_inverse(numerators, denominator, args.bits)
    product = matrix_product(guess, numerators)
    scale = denominator * 2**args.bits
    residual = [[(scale if i == j else 0) - x for j, x in enumerate(row)]
                for i, row in enumerate(product)]
    eta = one_norm(residual, scale)
    inverse_bound = one_norm(guess, 2**args.bits) / (1 - eta) if eta < 1 else None
    m = max(a[2] for a in addresses) + 1
    n = max(a[3] for a in addresses) + 1
    b = 23 + abs(args.center - args.base)
    tail = input_tail(m, n, args.base)
    q0 = 18 / (Q(1, 6) + args.base)
    majorant = b * tail * (1 + b * q0 * inverse_bound) if inverse_bound is not None else None
    report = dict(matrix=str(args.matrix), dimension=dimension,
                  base=str(args.base), center=str(args.center), mu=str(args.mu),
                  fourier_cutoff=m, chebyshev_cutoff=n,
                  inverse_defect=str(eta),
                  finite_inverse_bound=str(inverse_bound) if inverse_bound is not None else None,
                  input_tail_bound=str(tail), transfer_floor=str(b * tail),
                  transfer_majorant=str(majorant) if majorant is not None else None,
                  finite_center_accepted=eta < 1,
                  conditional_transfer_inequality=majorant is not None and majorant < 1,
                  physical_certificate=False,
                  limitation="centro real; fundo RefA; erro do fundo exato e contorno nao certificados")
    output = json.dumps(report, indent=2) + "\n"
    print(output, end="")
    if args.out:
        args.out.write_text(output, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
