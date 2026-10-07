#!/usr/bin/env python3
"""Caudas EXATAS de RefA e bound condicional de acoplamento baixo->distante.

Usa a desigualdade bilinear ||S Xi Gamma2(v,.)|| <= 3||v|| de RT.
A incerteza do fundo e adicionada separadamente, nunca tomada como zero.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import re

from check_free_preconditioner import input_tail

PAIR = re.compile(r"\((-?\d+),(-?\d+)\)")
RT_BACKGROUND_ERROR = Q(1, 2**25) + Q(1, 2**277)


def weighted_coefficients(path: Path) -> list[tuple[int, int, Q]]:
    text = path.read_text(encoding="ascii")
    blocks = re.split(r"(?m)^Field\s*$", text)[1:]
    if len(blocks) != 4:
        raise ValueError("esperados quatro campos RT")
    result = []
    for block in blocks:
        header = {}
        for name in ("TwoExp", "off_m", "off_n", "num_m", "num_n"):
            match = re.search(rf"(?m)^{name}\s+(-?\d+)\s*$", block)
            if match is None:
                raise ValueError(f"cabecalho ausente: {name}")
            header[name] = int(match.group(1))
        if min(header["off_m"], header["off_n"]) < 0 or min(header["num_m"], header["num_n"]) <= 0:
            raise ValueError("setor invalido")
        pairs = list(PAIR.finditer(block))
        if len(pairs) != header["num_m"] * header["num_n"]:
            raise ValueError("quantidade de coeficientes incorreta")
        wm = [(2-(m == 0))*Q(65, 64)**m
              for m in range(header["off_m"], header["off_m"]+header["num_m"])]
        wn = [(2-(n == 0))*Q(5, 4)**n
              for n in range(header["off_n"], header["off_n"]+header["num_n"])]
        scale = Q(2)**header["TwoExp"]
        for index, pair in enumerate(pairs):
            amplitude = abs(int(pair[1])) + abs(int(pair[2]))
            if amplitude:
                i, j = divmod(index, header["num_n"])
                result.append((header["off_m"]+i, header["off_n"]+j,
                               amplitude*scale*wm[i]*wn[j]))
    return result


def tail(coefficients: list[tuple[int, int, Q]], m: int, n: int) -> Q:
    if m < 1 or n < 1:
        raise ValueError("cortes devem ser positivos")
    return sum((weight for fm, cn, weight in coefficients if fm >= m or cn >= n), Q(0))


def guard_tail(coefficients: list[tuple[int, int, Q]], core: tuple[int, int],
               outer: tuple[int, int]) -> Q:
    if min(core) < 1 or any(b < a for a, b in zip(core, outer)):
        raise ValueError("a caixa externa deve conter o nucleo")
    # Indices mantidos sao estritamente menores que os cortes.
    return tail(coefficients, outer[0]-core[0]+1, outer[1]-core[1]+1)


def reverse_coupling(coefficients: list[tuple[int, int, Q]], core: tuple[int, int],
                     middle: tuple[int, int], outer: tuple[int, int],
                     base: Q = Q(0), distance: Q = Q(5,4)) -> dict:
    """Bound de P D Q0 (I-G), antes da multiplicacao pela parametriz V."""
    if (base < 0 or distance < 0 or min(core) < 1
            or any(not p <= u <= g for p, u, g in zip(core, middle, outer))):
        raise ValueError("exige P contido em U contido em G, base e distancia nao negativas")
    r = Q(4,5)
    q_middle_from_far = 4/(Q(1,6)+base) * r**(outer[1]-middle[1]+1)/(1-r)
    gamma1 = 4*Q(17,100)*r**(middle[1]-core[1]+1)/(1-r*r)
    nonlinear = 6*(guard_tail(coefficients, core, middle)+RT_BACKGROUND_ERROR)
    result = ((23+distance)*q_middle_from_far
              + (gamma1+nonlinear)*input_tail(*outer, base))
    return dict(core=core, middle=middle, outer=outer,
                base=str(base), max_distance=str(distance),
                q_middle_from_far=str(q_middle_from_far),
                reverse_H_bound=str(result), reverse_H_bound_decimal=float(result))


def crown_exterior(core: tuple[int, int], outer: tuple[int, int],
                   base: Q = Q(0), distance: Q = Q(5, 4)) -> dict:
    """Majorantes uniformes condicionais a beta<=23, nao certificados fisicos."""
    if distance < 0 or any(g <= p for p, g in zip(core, outer)):
        raise ValueError("exige coroa estrita e distancia nao negativa")
    tp, tg = input_tail(*core, base), input_tail(*outer, base)
    # T Q0 S = 0, pois Q0 preserva G.
    bounds = dict(T_from_S=23*tp, S_from_T=(23+distance)*tg,
                  TT_defect=(23+distance)*tg)
    return dict(core=core, outer=outer, base=str(base), max_distance=str(distance),
                bounds={k: str(v) for k, v in bounds.items()},
                decimals={k: float(v) for k, v in bounds.items()},
                exterior_contraction=bounds["TT_defect"] < 1,
                physical_certificate=False)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("background", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    coefficients = weighted_coefficients(args.background)
    norm = sum((w for _, _, w in coefficients), Q(0))
    rows = []
    for core, outer in (((4,12),(6,18)), ((4,12),(12,36)), ((4,12),(24,72)),
                        ((4,12),(45,113)), ((12,36),(24,72)), ((12,36),(53,137))):
        exact = guard_tail(coefficients, core, outer)
        bound = 6*(exact+RT_BACKGROUND_ERROR)
        rows.append(dict(core=core, outer=outer,
                         reference_tail=str(exact), reference_tail_decimal=float(exact),
                         conditional_B_far_from_core=str(bound),
                         conditional_B_far_from_core_decimal=float(bound)))
    report = dict(background=str(args.background),
                  sha256=hashlib.sha256(args.background.read_bytes()).hexdigest(),
                  reference_norm=str(norm), reference_norm_decimal=float(norm),
                  assumed_background_error=str(RT_BACKGROUND_ERROR),
                  assumption="bound publicado RT relativo a RefA, nao validado pelo hash",
                  guards=rows, physical_certificate=False)
    report["reverse_couplings"] = [reverse_coupling(coefficients, *boxes)
        for boxes in (((4,12),(12,36),(24,72)), ((4,12),(24,72),(45,113)),
                      ((4,12),(24,72),(64,192)), ((12,36),(32,96),(64,192)))]
    report["crown_exterior"] = [crown_exterior((12,36), box)
        for box in ((24,72), (64,192), (874,2621))]
    output = json.dumps(report, indent=2)+"\n"
    print(output, end="")
    if args.out:
        args.out.write_text(output, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
