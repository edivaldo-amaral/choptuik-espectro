#!/usr/bin/env python3
"""Analisa o truncamento de L_A(s)=L_A(0)+sI; não certifica o espectro."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
from scipy.linalg import eig

K_RT = 1.7227262011139106749857559727881918622210
S_TARGET = K_RT * 2.674 / (2.0 * math.pi)


def read_matrix(path: Path, expected_magic: str = "CHOPTUIK_RT_SPECTRAL_MATRIX_V1") -> tuple[np.ndarray, list[tuple[int, int, int, int]]]:
    with path.open(encoding="ascii") as stream:
        if stream.readline().strip() != expected_magic:
            raise ValueError("formato de matriz desconhecido")
        key, raw_dimension = stream.readline().split()
        if key != "dimension":
            raise ValueError("cabeçalho sem dimensão")
        dimension = int(raw_dimension)
        stream.readline()  # sector
        stream.readline()  # background
        if stream.readline().strip() != "BEGIN_ADDRESSES":
            raise ValueError("endereços ausentes")
        addresses: list[tuple[int, int, int, int]] = []
        for expected_index in range(dimension):
            index, component, part, fourier, chebyshev = map(int, stream.readline().split())
            if index != expected_index:
                raise ValueError("índice de endereço fora de ordem")
            addresses.append((component, part, fourier, chebyshev))
        if stream.readline().strip() != "BEGIN_ENTRIES":
            raise ValueError("entradas ausentes")
        matrix = np.zeros((dimension, dimension), dtype=float)
        seen: set[tuple[int, int]] = set()
        for line_number, line in enumerate(stream, dimension + 8):
            fields = list(map(int, line.split()))
            if len(fields) == 3:  # o C do RT omite 2^0
                row, column, coefficient = fields
                exponent = 0
            elif len(fields) == 4:
                row, column, coefficient, exponent = fields
            else:
                raise ValueError(f"entrada malformada na linha {line_number}")
            if not (0 <= row < dimension and 0 <= column < dimension):
                raise ValueError(f"entrada fora da matriz na linha {line_number}")
            if (row, column) in seen:
                raise ValueError(f"entrada duplicada na linha {line_number}")
            seen.add((row, column))
            matrix[row, column] = math.ldexp(float(coefficient), exponent)
    return matrix, addresses


def read_constraints(path: Path, columns: int) -> tuple[np.ndarray, np.ndarray]:
    with path.open(encoding="ascii") as stream:
        if stream.readline().strip() != "CHOPTUIK_RT_SPECTRAL_CONSTRAINTS_V1":
            raise ValueError("formato de constraints desconhecido")
        rows = int(stream.readline().split()[1])
        actual_columns = int(stream.readline().split()[1])
        if actual_columns != columns:
            raise ValueError("número de colunas incompatível nas constraints")
        if stream.readline().strip() != "BEGIN_ADDRESSES":
            raise ValueError("endereços das constraints ausentes")
        for expected_index in range(rows):
            index, *_ = map(int, stream.readline().split())
            if index != expected_index:
                raise ValueError("índice de constraint fora de ordem")
        if stream.readline().strip() != "BEGIN_ENTRIES":
            raise ValueError("entradas das constraints ausentes")
        constant = np.zeros((rows, columns), dtype=float)
        linear = np.zeros((rows, columns), dtype=float)
        seen: set[tuple[int, int, int]] = set()
        for line_number, line in enumerate(stream, rows + 7):
            fields = list(map(int, line.split()))
            if len(fields) == 4:  # o C do RT omite 2^0
                degree, row, column, coefficient = fields
                exponent = 0
            elif len(fields) == 5:
                degree, row, column, coefficient, exponent = fields
            else:
                raise ValueError(f"entrada de constraint malformada na linha {line_number}")
            if degree not in (0, 1):
                raise ValueError(f"grau inválido na linha {line_number}")
            key = (degree, row, column)
            if key in seen:
                raise ValueError(f"entrada de constraint duplicada na linha {line_number}")
            seen.add(key)
            target = constant if degree == 0 else linear
            target[row, column] = math.ldexp(float(coefficient), exponent)
    return constant, linear


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("matrix", type=Path)
    parser.add_argument("--count", type=int, default=12)
    parser.add_argument("--show-b", action="store_true", help="mostra o setor B via deslocamento exato -i/2")
    args = parser.parse_args()

    matrix, _ = read_matrix(args.matrix)
    # Se q é autovalor de L_A(0), s=-q resolve det(L_A(s))=0.
    eigenvalues, eigenvectors = eig(matrix)
    exponents_s = -eigenvalues
    scale = 2.0 * math.pi / K_RT
    candidate = min(exponents_s, key=lambda value: abs(value - S_TARGET))
    candidate_index = int(np.argmin(np.abs(exponents_s - S_TARGET)))
    phase_candidate = min(exponents_s, key=abs)

    print(f"dimensão real: {matrix.shape[0]}")
    print("AVISO: resultado de Galerkin em ponto flutuante; NÃO É PROVA.")
    print(
        "candidato físico: "
        f"s={candidate.real:+.12f}{candidate.imag:+.12f}i, "
        f"lambda={scale*candidate.real:+.12f}{scale*candidate.imag:+.12f}i"
    )
    constraint_path = Path(f"{args.matrix}.constraints")
    if constraint_path.exists():
        constant, linear = read_constraints(constraint_path, matrix.shape[0])
        vector = eigenvectors[:, candidate_index]
        residual = np.linalg.norm((constant + candidate * linear) @ vector) / np.linalg.norm(vector)
        print(f"resíduo l2 das constraints do candidato: {residual:.12e}")
    else:
        print("constraints: arquivo auxiliar ausente")
    print(
        "candidato mais próximo do modo de fase: "
        f"s={phase_candidate.real:+.12f}{phase_candidate.imag:+.12f}i"
    )
    print("raízes mais próximas do alvo histórico:")
    for value in sorted(exponents_s, key=lambda item: abs(item - S_TARGET))[: args.count]:
        physical = scale * value
        print(
            f"  s={value.real:+.9f}{value.imag:+.9f}i  "
            f"lambda={physical.real:+.9f}{physical.imag:+.9f}i"
        )
    if args.show_b:
        print("setor B (conjugação exata M=exp(i*tau_R/2)): s_B=s_A-i/2")
        for value in sorted(exponents_s, key=lambda item: abs(item - S_TARGET))[: args.count]:
            shifted = value - 0.5j
            physical = scale * shifted
            print(
                f"  s_B={shifted.real:+.9f}{shifted.imag:+.9f}i  "
                f"lambda={physical.real:+.9f}{physical.imag:+.9f}i"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
