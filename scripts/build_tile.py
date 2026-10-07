#!/usr/bin/env python3
"""Constroi um `Tile` de Neumann a partir da matriz espectral real.

Le `build/spectrum/rt-A-*.dat` (formato `CHOPTUIK_RT_SPECTRAL_MATRIX_V1`,
entradas diadicas exatas) e emite o bloco Lean de um
`Choptuik.RationalCertificate.Tile` cujo `centerMatrix` e
`F_N(z_j) = I + delta * R_N`, com `R_N = (A_N + z_0 I)^-1`,
`A_N` a matriz do arquivo e `delta = z_j - z_0`.
Essa resolvente finita NAO e automaticamente a compressao da resolvente
infinita: a diferenca exige um bound adicional (docs/ANALYTIC_AUDIT.md).
`--raw-coefficient` preserva o antigo benchmark I + delta*A_N, sem
interpretacao como certificado do pencil espectral.

O parser exato e reaproveitado de `scripts/build_contour.py` (nao editado).

PONTO FLUTUANTE, E POR QUE ELE NAO ESTA NO CAMINHO DA PROVA
-----------------------------------------------------------

O campo `approxInverse` (`V`) e calculado aqui em ponto flutuante
(`numpy.linalg.inv`) e depois **arredondado para racionais de denominador
`2^bits`**. Isso nao fere a regra epistemologica do projeto, e vale explicar
por que, porque e o unico lugar do repositorio em que um `float` toca um dado
que entra num certificado:

* `V` e um **chute**. O certificado nao afirma nada sobre ele e nao confia
  nele. Ele nao precisa ser o inverso, nem estar perto do inverso.
* O que o certificado afirma e `inverseDefect < 1`, onde
  `inverseDefect = ||I - V F|| + ||V|| * matrixVariation`. Essa quantidade
  **mede** o quao ruim `V` e, e e recomputada do zero pelo kernel em `Rat`
  exato a partir dos racionais emitidos.
* Se o `float` mentir, o unico efeito possivel e o certificado ser **recusado**
  (defeito grande demais). Nao ha como um `V` ruim fazer um tile falso passar.

Em outras palavras: o ponto flutuante entra como heuristica de busca, nunca
como testemunha. E exatamente o mesmo papel que ele tem no programa de
Reiterer-Trubowitz.

ARREDONDAR O CENTRO (`--center-bits`)
-------------------------------------

`--center-bits b` emite como `centerMatrix` o arredondado de `F_N(z_j)` para
multiplos de `2^-b`, e emite **a parte descartada** como `centerResidual`. A
soma das duas e o centro exato, entrada por entrada, que e exatamente a
obrigacao declarada no campo `centerResidual` do `Tile`.

Isto e solido, e nao por convencao: `Tile.matrixVariation` passou a ser
`raio * ||C_N|| + ||centerResidual||_1`, e o kernel **recomputa** a segunda
parcela a partir da matriz de residuo escrita no arquivo. O arredondamento nao
e estimado nem prometido -- e pago, em `inverseDefect`. Um arredondamento
grosso demais nao produz certificado otimista: produz certificado **recusado**
(ver `rejectedRoundingTile` em `RationalCertificate.lean`).

O inverso aproximado `V` e calculado contra o centro **arredondado**, que e a
matriz que aparece no produto `I - V F`.

USO

    scripts/build_tile.py build/spectrum/rt-A-4x12.dat --dim 16 --out /tmp/T.lean

`--dim k` usa a submatriz principal `k x k` (para medir a curva de custo);
omitido, usa a dimensao inteira.
"""
from __future__ import annotations

import argparse
import math
import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_contour import read_spectral_matrix  # noqa: E402  (parser exato reaproveitado)


def exact_inverse(matrix: list[list[Fraction]]) -> list[list[Fraction]]:
    """Gauss-Jordan racional com pivotamento; sem arredondamento."""
    size = len(matrix)
    work = [list(row) + [Fraction(i == j) for j in range(size)]
            for i, row in enumerate(matrix)]
    for column in range(size):
        pivot = next((row for row in range(column, size)
                      if work[row][column]), None)
        if pivot is None:
            raise ValueError("ponto-base singular: A_N + z_0 I nao e invertivel")
        work[column], work[pivot] = work[pivot], work[column]
        divisor = work[column][column]
        work[column] = [value / divisor for value in work[column]]
        for row in range(size):
            if row != column:
                factor = work[row][column]
                if factor:
                    work[row] = [a - factor * b
                                 for a, b in zip(work[row], work[column])]
    return [row[size:] for row in work]


def coefficient_matrix(matrix: list[list[int]], exponent: int,
                       base_point: Fraction, raw: bool = False,
                       weights: list[Fraction] | None = None
                       ) -> tuple[list[list[int]], int]:
    """Numerador e denominador comum de R_N (ou A_N no benchmark)."""
    scale = Fraction(2) ** exponent
    exact = [[value * scale for value in row] for row in matrix]
    if not raw:
        for index in range(len(exact)):
            exact[index][index] += base_point
        exact = exact_inverse(exact)
    if weights is not None:
        if len(weights) != len(exact) or any(weight <= 0 for weight in weights):
            raise ValueError("pesos devem ser positivos, um por coordenada")
        exact = [[value * weights[i] / weights[j] for j, value in enumerate(row)]
                 for i, row in enumerate(exact)]
    denominator = math.lcm(*(value.denominator for row in exact for value in row))
    return ([[value.numerator * (denominator // value.denominator)
              for value in row] for row in exact], denominator)


def read_rt_weights(path: Path, dimension: int) -> list[Fraction]:
    """Diagonal dos pesos RT na ordem dos DOFs reais exportados."""
    addresses = {}
    active = False
    with path.open(encoding="ascii") as source:
        for line in source:
            if line.strip() == "BEGIN_ADDRESSES":
                active = True
                continue
            if line.strip() == "BEGIN_ENTRIES":
                break
            if active:
                index, component, part, m, n = map(int, line.split())
                if (index in addresses or not 0 <= index < dimension
                        or component not in range(4) or part not in (0, 1)
                        or m < 0 or n < 0):
                    raise ValueError("enderecos RT invalidos")
                addresses[index] = ((2 - (m == 0)) * (2 - (n == 0))
                                    * Fraction(65, 64) ** m * Fraction(5, 4) ** n)
    if len(addresses) != dimension:
        raise ValueError("faltam enderecos para calcular os pesos RT")
    return [addresses[index] for index in range(dimension)]


def approximate_inverse(center: list[list[int]], scale: int, bits: int) -> list[list[int]]:
    """`V ~ F^-1` em ponto flutuante, arredondado para multiplos de `2^-bits`.

    `center` e o numerador inteiro de `F` na escala `scale` (isto e,
    `F = center / scale`). Devolve o numerador inteiro de `V` na escala
    `2^bits`.
    """
    import numpy as np

    size = len(center)
    dense = np.empty((size, size), dtype=float)
    for row in range(size):
        for col in range(size):
            dense[row][col] = center[row][col] / scale
    inverse = np.linalg.inv(dense)
    grid = 1 << bits
    return [[int(round(float(inverse[row][col]) * grid)) for col in range(size)]
            for row in range(size)]


def round_to_grid(numerator: int, step: int) -> int:
    """Arredonda `numerator/step` para o inteiro mais proximo (metade para longe
    de zero). Exato: so aritmetica inteira."""
    if numerator >= 0:
        return (2 * numerator + step) // (2 * step)
    return -((-2 * numerator + step) // (2 * step))


def one_norm(numerators: list[list[int]], denominator: int) -> Fraction:
    """Norma induzida por l^1: maxima soma de valores absolutos por COLUNA."""
    size = len(numerators)
    best = 0
    for col in range(size):
        total = 0
        for row in range(size):
            total += abs(numerators[row][col])
        best = max(best, total)
    return Fraction(best, denominator)


def block_one_norm(numerators: list[list[int]], denominator: int,
                   start: int, count: int) -> Fraction:
    """Norma `l^1` restrita as colunas `[start, start + count)`.

    E exatamente o que `Choptuik.TileBlocks.residualBlockNorm` calcula: maximo,
    sobre as colunas do bloco, da soma de valores absolutos. O maximo comeca em
    zero, como o `maxList` do Lean."""
    size = len(numerators)
    best = 0
    for col in range(start, start + count):
        total = 0
        for row in range(size):
            total += abs(numerators[row][col])
        best = max(best, total)
    return Fraction(best, denominator)


def column_blocks(size: int, width: int) -> list[tuple[int, int]]:
    """Blocos consecutivos de `width` colunas; o ultimo leva o resto."""
    blocks = []
    start = 0
    while start < size:
        blocks.append((start, min(width, size - start)))
        start += width
    return blocks


def matrix_product(left: list[list[int]], right: list[list[int]]) -> list[list[int]]:
    size = len(left)
    result = [[0] * size for _ in range(size)]
    for row in range(size):
        left_row = left[row]
        out_row = result[row]
        for inner in range(size):
            coefficient = left_row[inner]
            if coefficient:
                right_row = right[inner]
                for col in range(size):
                    out_row[col] += coefficient * right_row[col]
    return result


def lean_rational(numerator: int, denominator: int) -> str:
    if denominator == 1:
        return f"({numerator})" if numerator < 0 else str(numerator)
    head = f"({numerator})" if numerator < 0 else str(numerator)
    return f"mkRat {head} {denominator}"


def emit_matrix(name: str, numerators: list[list[int]], denominator: int) -> str:
    size = len(numerators)
    cells = []
    for row in range(size):
        for col in range(size):
            cells.append(lean_rational(numerators[row][col], denominator))
    body = ",\n    ".join(", ".join(cells[i:i + 4]) for i in range(0, len(cells), 4))
    return (f"def {name} : RatMatrix where\n"
            f"  rows := {size}\n  cols := {size}\n"
            f"  data :=\n   [{body}]\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("matrix", type=Path)
    parser.add_argument("--dim", type=int, default=None,
                        help="submatriz principal k x k (padrao: dimensao inteira)")
    parser.add_argument("--delta", type=str, default="1",
                        help="z_j - z_0, racional exato (padrao 1)")
    parser.add_argument("--base-point", default="1",
                        help="z_0 real racional da resolvente finita (padrao 1)")
    parser.add_argument("--raw-coefficient", action="store_true",
                        help="benchmark legado: usa A_N como coeficiente, nao resolvente")
    parser.add_argument("--rt-weights", action="store_true",
                        help="conjuga pela diagonal dos pesos RT antes da norma de colunas")
    parser.add_argument("--bits", type=int, default=40,
                        help="V arredondado para multiplos de 2^-bits")
    parser.add_argument("--center-bits", type=int, default=0,
                        help="centerMatrix arredondado para multiplos de 2^-b, com a "
                             "parte descartada emitida em centerResidual (0 = exato)")
    parser.add_argument("--max-heartbeats", type=int, default=None,
                        help="emite set_option maxHeartbeats N no arquivo gerado. "
                             "Omitido, o tile e verificado no orcamento PADRAO -- que e "
                             "a metrica que interessa, porque e a que um auditor externo "
                             "reproduz sem afrouxar limite nenhum")
    parser.add_argument("--block-cols", type=int, default=None,
                        help="parte a verificacao em blocos de tantas COLUNAS, uma "
                             "declaracao por bloco, e recompoe por qmax "
                             "(Choptuik.TileBlocks). Omitido, emite o tile monolitico")
    parser.add_argument("--tile-radius", type=str, default="1/1048576")
    parser.add_argument("--max-distance", type=str, default="2")
    parser.add_argument("--tail-error", type=str, default="1/1000000")
    parser.add_argument("--operator-bound", type=str, default="1")
    parser.add_argument("--name", type=str, default="realTile")
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()
    if args.bits < 0 or args.center_bits < 0:
        parser.error("--bits e --center-bits devem ser nao negativos")
    if args.block_cols is not None and args.block_cols <= 0:
        parser.error("--block-cols deve ser positivo")

    dimension, integer_matrix, exponent = read_spectral_matrix(args.matrix)
    size = dimension if args.dim is None else args.dim
    if not 1 <= size <= dimension:
        raise SystemExit(f"--dim fora de [1, {dimension}]")

    delta = Fraction(args.delta)
    base_point = Fraction(args.base_point)
    principal = [row[:size] for row in integer_matrix[:size]]
    try:
        weights = read_rt_weights(args.matrix, dimension)[:size] if args.rt_weights else None
        integer_matrix, shift = coefficient_matrix(
            principal, exponent, base_point, args.raw_coefficient, weights)
    except ValueError as error:
        raise SystemExit(str(error)) from error
    print("BENCHMARK I+delta*A_N" if args.raw_coefficient else
          f"RESOLVENTE FINITA (A_N+({base_point})I)^-1; cauda PDE nao certificada",
          file=sys.stderr)
    print("norma: l1 com pesos RT" if args.rt_weights else "norma: l1 sem pesos",
          file=sys.stderr)
    # F = I + delta * R_N, denominador comum exato (nao necessariamente diadico).
    scale = shift * delta.denominator
    center = [[delta.numerator * integer_matrix[row][col] for col in range(size)]
              for row in range(size)]
    for index in range(size):
        center[index][index] += scale

    # Arredondamento do centro. `center_residual` guarda a parte descartada, de
    # modo que `centro_arredondado + residuo = centro_exato` entrada por entrada
    # -- a obrigacao declarada no campo `centerResidual` do `Tile`.
    center_residual = None
    center_scale = scale
    if args.center_bits:
        target = 1 << args.center_bits
        residual_denominator = math.lcm(scale, target)
        rounded = [[round_to_grid(center[row][col] * target, scale) for col in range(size)]
                   for row in range(size)]
        center_residual = [[center[row][col] * (residual_denominator // scale)
                            - rounded[row][col] * (residual_denominator // target)
                            for col in range(size)] for row in range(size)]
        scale = residual_denominator
        center = rounded
        center_scale = target

    inverse = approximate_inverse(center, center_scale, args.bits)
    inverse_scale = 1 << args.bits

    # residuo = I - V F, numeradores na escala center_scale * inverse_scale
    product = matrix_product(inverse, center)
    residual_scale = center_scale * inverse_scale
    residual = [[(residual_scale if row == col else 0) - product[row][col]
                 for col in range(size)] for row in range(size)]

    residual_norm = one_norm(residual, residual_scale)
    inverse_norm = one_norm(inverse, inverse_scale)
    matrix_norm = one_norm([[integer_matrix[r][c] for c in range(size)]
                            for r in range(size)], shift)
    center_error = (one_norm(center_residual, scale) if center_residual is not None
                    else Fraction(0))

    tile_radius = Fraction(args.tile_radius)
    max_distance = Fraction(args.max_distance)
    tail_error = Fraction(args.tail_error)
    operator_bound = Fraction(args.operator_bound)

    variation = tile_radius * matrix_norm + center_error
    defect = residual_norm + inverse_norm * variation
    bound = inverse_norm / (1 - defect) if defect != 1 else None
    majorant = (max_distance * tail_error * (1 + max_distance * bound * operator_bound)
                if bound is not None else None)

    print(f"dim={size}  ||I-VF||_1={float(residual_norm):.6e}  ||V||_1={float(inverse_norm):.6e}"
          f"  ||C_N||_1={float(matrix_norm):.6e}", file=sys.stderr)
    print(f"  Vbits={args.bits}  centerBits={args.center_bits or 'exato'}"
          f"  centerError={float(center_error):.6e}", file=sys.stderr)
    print(f"  matrixVariation={float(variation):.6e}  inverseDefect={float(defect):.6e}"
          f"  (<1: {defect < 1})", file=sys.stderr)
    if bound is not None:
        print(f"  inverseBound={float(bound):.6e}  transferMajorant={float(majorant):.6e}"
              f"  (<1: {majorant < 1})", file=sys.stderr)
    valid = (all(value >= 0 for value in
                 (tile_radius, max_distance, tail_error, operator_bound))
             and defect < 1 and majorant is not None and majorant < 1)
    print(f"  tile valido (aritmetica em Fraction): {valid}", file=sys.stderr)

    name = args.name
    blocks = (column_blocks(size, args.block_cols) if args.block_cols else None)
    if blocks is not None:
        widths = ", ".join(str(count) for _, count in blocks)
        print(f"  blocos de colunas ({len(blocks)}): {widths}", file=sys.stderr)
    text = [
        ("import ChoptuikFormal.TileBlocks" if blocks is not None
         else "import ChoptuikFormal.RationalCertificate"),
        "",
        "/-! Gerado por scripts/build_tile.py a partir da matriz espectral real.",
        ("Benchmark legado I+delta*A_N; nao e a resolvente espectral." if args.raw_coefficient
         else f"Resolvente FINITA no ponto-base {base_point}; transferencia PDE aberta."),
        ("Matrizes conjugadas pelos pesos RT." if args.rt_weights else "Norma l1 sem pesos."),
        "NAO adicionar ao repositorio: e um artefato de medicao. -/",
        "",
        "namespace Choptuik.RealTile",
        "open Choptuik.RationalCertificate",
    ]
    if blocks is not None:
        text += ["open Choptuik.TileBlocks"]
    text += [
        "",
        "unseal Rat.mul Rat.add Rat.sub Rat.inv",
        "",
        "set_option maxRecDepth 1000000",
    ]
    if args.max_heartbeats is not None:
        text += [f"set_option maxHeartbeats {args.max_heartbeats}"]
    text += [
        "",
        emit_matrix(f"{name}Center", center, center_scale),
        "",
        emit_matrix(f"{name}Inverse", inverse, inverse_scale),
        "",
    ]
    if center_residual is not None:
        text += [emit_matrix(f"{name}CenterResidual", center_residual, scale), ""]
    text += [
        f"def {name} : Tile where",
        f"  centerMatrix := {name}Center",
        f"  approxInverse := {name}Inverse",
        f"  tileRadius := {lean_rational(tile_radius.numerator, tile_radius.denominator)}",
        f"  matrixNorm := {lean_rational(matrix_norm.numerator, matrix_norm.denominator)}",
        f"  maxDistance := {lean_rational(max_distance.numerator, max_distance.denominator)}",
        f"  tailError := {lean_rational(tail_error.numerator, tail_error.denominator)}",
        f"  operatorBound := {lean_rational(operator_bound.numerator, operator_bound.denominator)}",
    ]
    if center_residual is not None:
        text += [f"  centerResidual := {name}CenterResidual"]
    if blocks is None:
        text += [
            "",
            f"example : {name}.valid = true := by decide",
        ]
    else:
        # Uma declaracao por bloco de colunas: e aqui que o custo se distribui.
        # Os valores sao exatos (Fraction), e o kernel os recomputa.
        block_norms = [block_one_norm(residual, residual_scale, start, count)
                       for start, count in blocks]
        assert max(block_norms) == residual_norm, "composicao dos blocos != norma total"
        text += ["", "/-! Um bloco de colunas por declaracao. Cada uma cabe no orcamento",
                 "padrao; a soma delas e o mesmo trabalho que o monolitico faria de uma",
                 "vez so. -/", ""]
        text += [f"theorem {name}Shape : {name}.shapeCorrect := by decide", ""]
        for index, ((start, count), value) in enumerate(zip(blocks, block_norms)):
            text += [f"theorem {name}Block{index} : residualBlockNorm {name} {start} {count}"
                     f" = {lean_rational(value.numerator, value.denominator)} := by decide",
                     ""]
        # Recomposicao por qmax: puro remanejo de teoremas, sem recomputar nada.
        text += ["/-! A recomposicao: `qmax` dos blocos. Nenhuma entrada e reexaminada. -/",
                 ""]
        lines = []
        remaining = size
        for index, (start, count) in enumerate(blocks[:-1]):
            rest = remaining - count
            nxt = start + count
            lines.append(f"      residualBlockNorm_split {name} {start} {remaining} "
                         f"{count} {rest} {nxt} rfl rfl,")
            remaining = rest
        for index in range(len(blocks)):
            lines.append(f"      {name}Block{index},")
        lines[-1] = lines[-1].rstrip(",")
        text += [f"theorem {name}Residual : residualBlockNorm {name} 0 {size}"
                 f" = {lean_rational(residual_norm.numerator, residual_norm.denominator)} := by",
                 "  rw [" + "\n".join(lines).lstrip() + "]",
                 "  decide",
                 ""]
        text += [f"theorem {name}InverseNorm : {name}.approxInverse.oneNorm"
                 f" = {lean_rational(inverse_norm.numerator, inverse_norm.denominator)}"
                 " := by decide", ""]
        text += [f"/-- Mesma proposicao que `by decide` no tile inteiro produziria. -/",
                 f"theorem {name}Valid : {name}.valid = true :=",
                 f"  tileValidOfBlockNorms {name} "
                 f"({lean_rational(residual_norm.numerator, residual_norm.denominator)}) "
                 f"({lean_rational(inverse_norm.numerator, inverse_norm.denominator)})",
                 f"    {name}Shape {name}Residual {name}InverseNorm",
                 "    (by decide) (by decide) (by decide) (by decide) (by decide)",
                 "    (by decide) (by decide)",
                 "",
                 f"#print axioms {name}Valid"]
    text += [
        "",
        "end Choptuik.RealTile",
        "",
    ]
    output = "\n".join(text)
    if args.out is None:
        sys.stdout.write(output)
    else:
        args.out.write_text(output, encoding="ascii")
        print(f"  escrito: {args.out} ({len(output)} bytes)", file=sys.stderr)
    return 0 if valid else 1


if __name__ == "__main__":
    raise SystemExit(main())
