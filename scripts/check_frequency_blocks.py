#!/usr/bin/env python3
"""Bounds racionais separados para nucleo e coroa de uma matriz de Galerkin.

Nao confunde a coroa finita com a cauda infinita. Mede o defeito de uma
parametriz inversa diagonal por blocos e testa contracao em norma ponderada.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import json
from pathlib import Path

from build_contour import read_spectral_matrix
from build_tile import approximate_inverse, exact_inverse, read_rt_weights
from check_free_preconditioner import (free_matrix, integer_matrix, multiply,
                                       read_addresses)


def rectangular_product(left: list[list[int]], right: list[list[int]]) -> list[list[int]]:
    if not left or not right or any(len(row) != len(right) for row in left):
        raise ValueError("dimensoes incompativeis")
    columns = len(right[0])
    if any(len(row) != columns for row in right):
        raise ValueError("matriz irregular")
    result = [[0] * columns for _ in left]
    for i, row in enumerate(left):
        for k, coefficient in enumerate(row):
            if coefficient:
                for j, value in enumerate(right[k]):
                    result[i][j] += coefficient * value
    return result


def rectangular_norm(matrix: list[list[int]], denominator: int) -> Q:
    if denominator <= 0 or not matrix or not matrix[0]:
        raise ValueError("norma exige matriz nao vazia e denominador positivo")
    return Q(max(sum(abs(row[j]) for row in matrix) for j in range(len(matrix[0]))), denominator)


def two_block_test(a: Q, b: Q, c: Q, d: Q) -> dict:
    """E=[[a,b],[c,d]]; norma |x_L|+r|x_H|, coluna e nao linha."""
    if min(a, b, c, d) < 0:
        raise ValueError("majorantes negativos")
    gap = (1-a)*(1-d) - b*c
    accepted = a < 1 and d < 1 and gap > 0
    result = dict(accepted=accepted, schur_gap=str(gap), weight=None, contraction=None)
    if not accepted:
        return result
    lower = b/(1-d)
    upper = (1-a)/c if c else None
    weight = (lower+upper)/2 if upper is not None else lower+1
    contraction = max(a+weight*c, d+b/weight)
    assert weight > 0 and contraction < 1
    result.update(weight=str(weight), contraction=str(contraction))
    return result


def transport_defect(bounds: list[list[Q]], old_weights: list[Q],
                     new_weights: list[Q]) -> list[list[Q]]:
    """Conjugacao diagonal de uma parametriz JA validada, com perda segura.

Os pesos sao de componentes; cada bloco de frequencias pode conter todas.
Nao transfere a validade do centro para um contorno ou para outro fundo.
"""
    if (not old_weights or len(old_weights) != len(new_weights)
            or min(old_weights+new_weights) <= 0
            or any(v < 0 for row in bounds for v in row)):
        raise ValueError('pesos positivos e bounds nao negativos exigidos')
    ratios = [b/a for a, b in zip(old_weights, new_weights)]
    condition = max(ratios)/min(ratios)
    return [[condition*v for v in row] for row in bounds]


def audit_blocks(matrix: list[list[Q]], core: list[int], shell: list[int],
                 bits: int, invert_shell: bool, coupled: bool = False) -> dict:
    numerators, denominator = integer_matrix(matrix)
    groups = (core, shell)
    bounds = [[Q(0), Q(0)], [Q(0), Q(0)]]
    inverse_norms = []
    whole_inverse = approximate_inverse(numerators, denominator, bits) if coupled else None
    for target, rows in enumerate(groups):
        diagonal = [[numerators[i][j] for j in rows] for i in rows]
        grid = 2**bits
        inverse = ([[whole_inverse[i][j] for j in range(len(matrix))] for i in rows]
                   if whole_inverse is not None else
                   approximate_inverse(diagonal, denominator, bits)
                   if target == 0 or invert_shell else
                   [[grid if i == j else 0 for j in range(len(rows))]
                    for i in range(len(rows))])
        inverse_norms.append(rectangular_norm(inverse, grid))
        for source, columns in enumerate(groups):
            input_rows = range(len(matrix)) if coupled else rows
            block = [[numerators[i][j] for j in columns] for i in input_rows]
            product = rectangular_product(inverse, block)
            scale = denominator * grid
            residual = [[(scale if target == source and i == j else 0) - value
                         for j, value in enumerate(row)] for i, row in enumerate(product)]
            bounds[target][source] = rectangular_norm(residual, scale)
    a, b = bounds[0]
    c, d = bounds[1]
    return dict(core_dimension=len(core), shell_dimension=len(shell),
                shell_preconditioner="coupled" if coupled else "inverse" if invert_shell else "identity",
                defect_bounds=[[str(x) for x in row] for row in bounds],
                defect_bounds_decimal=[[float(x) for x in row] for row in bounds],
                inverse_norms=[str(x) for x in inverse_norms],
                unscaled_block_majorant=str(max(a+c, b+d)),
                **two_block_test(a, b, c, d))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("matrix", type=Path)
    parser.add_argument("--core-m", type=int, required=True)
    parser.add_argument("--core-n", type=int, required=True)
    parser.add_argument("--base", type=Q, default=Q(0))
    parser.add_argument("--center", type=Q, default=Q(1))
    parser.add_argument("--bits", type=int, default=40)
    parser.add_argument("--component-weights", default="1,1,1,1",
                        help="quatro pesos racionais positivos, separados por virgulas")
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    if min(args.core_m, args.core_n) < 1 or args.base < 0 or args.bits < 0:
        parser.error("cortes positivos, base>=0, bits>=0")
    try:
        component_weights = [Q(x) for x in args.component_weights.split(',')]
    except (ValueError, ZeroDivisionError):
        parser.error("pesos de componentes devem ser racionais")
    if len(component_weights) != 4 or min(component_weights) <= 0:
        parser.error("exige quatro pesos de componentes positivos")
    dimension, raw, exponent = read_spectral_matrix(args.matrix)
    addresses = read_addresses(args.matrix, dimension)
    core = [i for i, (_, _, m, n) in enumerate(addresses)
            if m < args.core_m and n < args.core_n]
    core_set = set(core)
    shell = [i for i in range(dimension) if i not in core_set]
    if not core or not shell:
        parser.error("nucleo e coroa devem ser nao vazios")
    weights = read_rt_weights(args.matrix, dimension)
    weights = [w*component_weights[address[0]] for w, address in zip(weights, addresses)]
    principal = free_matrix(addresses, Q(722873400, 2**32))
    for i in range(dimension):
        principal[i][i] += args.base
    inverse_free = exact_inverse(principal)
    pencil = [[Q(x)*Q(2)**exponent for x in row] for row in raw]
    for i in range(dimension):
        pencil[i][i] += args.center
    relative = multiply(pencil, inverse_free)
    weighted = [[x*weights[i]/weights[j] for j, x in enumerate(row)]
                for i, row in enumerate(relative)]
    results = [audit_blocks(weighted, core, shell, args.bits, choice)
               for choice in (False, True)]
    results.append(audit_blocks(weighted, core, shell, args.bits, True, coupled=True))
    report = dict(matrix=str(args.matrix), base=str(args.base), center=str(args.center),
                  component_weights=[str(w) for w in component_weights],
                  core_cutoff=[args.core_m, args.core_n], results=results,
                  infinite_tail_bounds=None, physical_certificate=False,
                  limitation="blocos da mesma matriz finita; exterior infinito nao estimado")
    output = json.dumps(report, indent=2)+"\n"
    print(output, end="")
    if args.out:
        args.out.write_text(output, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
