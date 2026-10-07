#!/usr/bin/env python3
"""Constrói o arquivo de nós do contorno `Gamma_A` a partir da matriz espectral.

Este é o primeiro programa do repositório que vai da matriz de Galerkin real
até dados consumíveis pelo verificador de winding. Ele é **diagnóstico**, não
prova: a matriz `4x12` é um truncamento cuja origem é numérica, e o winding
deste polinômio é o winding do **truncamento**, não o do operador infinito.
Ver `docs/WINDING_DATA.md` e `docs/PROOF_OBLIGATIONS.md` (C1, C3).

O que é exato aqui, e o que não é
---------------------------------

* Exato: a leitura das entradas diádicas (`coeficiente · 2^expoente`), o
  polinômio `p(z) = det(A + zI)`, os valores nos nós, os coeficientes de Taylor
  em cada nó, os raios de inflação das arestas e todo arredondamento (sempre
  para fora). Nenhum número emitido passa por ponto flutuante.
* Não exato, e usado **apenas** para escolher parâmetros (quantos nós, onde
  cortar o contorno) e para imprimir diagnóstico: uma diagonalização em
  `float`. Nada disso entra no arquivo de saída.

Arquitetura, e por que ela é a prescrita em `docs/WINDING_DATA.md` §3.4
----------------------------------------------------------------------

A prescrição é a cadeia relativa

    h(z_0) = 1,   h(z_{k+1}) = h(z_k) · det[I + (z_{k+1}-z_k) C_N F(z_k)^{-1}] .

No pencil afim `P(s) = L_A(0) + s·I` vale `C_N = (A + z_0 I)^{-1}` e
`F(z) = (A + zI)(A + z_0 I)^{-1}`, de modo que o fator relativo é
`det(A + z_{k+1}I)/det(A + z_k I)` e a cadeia **telescopa exatamente** para

    h(z_k) = p(z_k) / p(z_0) ,   p(z) = det(A + zI) .

Ou seja: a cadeia relativa e o polinômio característico produzem os **mesmos
números**; o que muda é o custo de obtê-los. Medido nesta máquina, em `n = 138`:

* um determinante exato (Bareiss fraction-free) custa **42 s**, logo a cadeia
  avaliada nó a nó custa 42 s por nó — cerca de 6 h para 500 nós;
* o polinômio característico inteiro, por Hessenberg modular + CRT, custa
  **~300 s uma única vez**, e depois cada nó é uma avaliação de Horner
  (milissegundos).

Além do custo, o polinômio dá de graça o que a forma de Jacobi mais precisa:
`tr(F(z)^{-1}) = p'(z)/p(z)` **exatamente**, sem inverso aproximado `V`, sem
resíduo de Frobenius e sem o fator `n` que `WINDING_DATA.md` §3.2 mostrou ser
proibitivo. A prescrição não foi trocada; foi implementada pelo caminho barato,
e a igualdade dos dois caminhos é **verificada** contra Bareiss (`--verify`).

Bound de aresta
---------------

Em cada nó calcula-se o deslocamento de Taylor exato
`p(z_k + u) = Σ_m c_m u^m`. Para `|u| ≤ L` saem, com `|c| ≤ |Re c| + |Im c|`:

    delta_taylor = Σ_{m≥1} |c_m| L^m                    (enclosure direto)
    sup|p'| ≤ Σ_{m≥1} m |c_m| L^{m-1}                   (Jacobi, numerador)
    inf|p|  ≥ |c_0| - Σ_{m≥1} |c_m| L^m                 (Jacobi, denominador)
    K = sup|p'| / inf|p| ,  S ≤ m0/(1-K·L) ,  delta_jacobi = K·L·m0/(1-K·L)

Os dois são rigorosos; o programa relata os dois e emite o menor. `K·L < 1/2` é
conferido e reportado.
"""
from __future__ import annotations

import argparse
import math
import sys
import time
from fractions import Fraction
from pathlib import Path

import numpy as np

HEADER_WARNING = """# GERADO POR scripts/build_contour.py -- NAO E PROVA.
#
# AVISO DE DIAGNOSTICO. Os valores abaixo sao exatos como racionais, mas o
# objeto que eles descrevem e o determinante do TRUNCAMENTO de Galerkin
# {label}, cuja origem e numerica. O winding deste contorno e o winding da
# matriz truncada; ele NAO estabelece o winding do operador infinito, que
# depende das obrigacoes C1 (bound de cauda), S3 (constraints) e S4 (gauge),
# todas abertas. Ver docs/PROOF_OBLIGATIONS.md e docs/WINDING_DATA.md.
"""


# --------------------------------------------------------------------------
# 1. Leitura exata da matriz
# --------------------------------------------------------------------------

def read_spectral_matrix(path: Path, expected_magic: str = "CHOPTUIK_RT_SPECTRAL_MATRIX_V1") -> tuple[int, list[list[int]], int]:
    """Devolve `(n, M, E)` com `A = M · 2^E` e `M` inteira.

    As entradas do formato `CHOPTUIK_RT_SPECTRAL_MATRIX_V1` sao diadicas
    (`coeficiente · 2^expoente`), portanto racionais exatos. Aqui elas sao
    lidas como inteiros e um expoente comum, sem tocar em `float`.
    """
    lines = path.read_text(encoding="ascii").splitlines()
    if lines[0].strip() != expected_magic:
        raise SystemExit(f"{path}: formato de matriz desconhecido")
    key, raw_dimension = lines[1].split()
    if key != "dimension":
        raise SystemExit(f"{path}: cabecalho sem dimensao")
    dimension = int(raw_dimension)
    start = lines.index("BEGIN_ENTRIES") + 1

    entries: list[tuple[int, int, int, int]] = []
    exponents: set[int] = set()
    seen: set[tuple[int, int]] = set()
    for offset, line in enumerate(lines[start:], start=start + 1):
        fields = line.split()
        if len(fields) == 3:
            row, column, coefficient = fields
            exponent = "0"
        elif len(fields) == 4:
            row, column, coefficient, exponent = fields
        else:
            raise SystemExit(f"{path}:{offset}: entrada malformada")
        key_rc = (int(row), int(column))
        if key_rc in seen:
            raise SystemExit(f"{path}:{offset}: entrada duplicada")
        seen.add(key_rc)
        entries.append((int(row), int(column), int(coefficient), int(exponent)))
        exponents.add(int(exponent))

    smallest = min(exponents)
    matrix = [[0] * dimension for _ in range(dimension)]
    for row, column, coefficient, exponent in entries:
        matrix[row][column] = coefficient << (exponent - smallest)
    return dimension, matrix, smallest


# --------------------------------------------------------------------------
# 2. Polinomio caracteristico inteiro: Hessenberg modular + CRT
# --------------------------------------------------------------------------

def _charpoly_mod(matrix: list[list[int]], dimension: int, prime: int) -> list[int]:
    """Coeficientes de `det(yI + M)` modulo `prime`, indice = grau.

    Reducao a Hessenberg por transformacoes de SEMELHANCA (operacao de linha
    seguida da operacao de coluna inversa), que preservam o polinomio
    caracteristico, e depois a recorrencia classica de Hessenberg.
    """
    hessenberg = np.array([[(-value) % prime for value in row] for row in matrix],
                          dtype=np.int64)
    for column in range(dimension - 2):
        pivot = -1
        for row in range(column + 1, dimension):
            if hessenberg[row, column]:
                pivot = row
                break
        if pivot < 0:
            continue
        if pivot != column + 1:
            hessenberg[[column + 1, pivot], :] = hessenberg[[pivot, column + 1], :]
            hessenberg[:, [column + 1, pivot]] = hessenberg[:, [pivot, column + 1]]
        inverse = pow(int(hessenberg[column + 1, column]), prime - 2, prime)
        for row in range(column + 2, dimension):
            if not hessenberg[row, column]:
                continue
            factor = (int(hessenberg[row, column]) * inverse) % prime
            hessenberg[row, :] = (hessenberg[row, :]
                                  - factor * hessenberg[column + 1, :]) % prime
            hessenberg[:, column + 1] = (hessenberg[:, column + 1]
                                         + factor * hessenberg[:, row]) % prime

    polynomials = [np.zeros(dimension + 1, dtype=np.int64)
                   for _ in range(dimension + 1)]
    polynomials[0][0] = 1
    for size in range(1, dimension + 1):
        current = np.zeros(dimension + 1, dtype=np.int64)
        current[1:size + 1] = polynomials[size - 1][0:size]
        diagonal = int(hessenberg[size - 1, size - 1])
        current[0:size] = (current[0:size]
                           - diagonal * polynomials[size - 1][0:size]) % prime
        product = 1
        for step in range(1, size):
            product = (product * int(hessenberg[size - step, size - step - 1])) % prime
            coefficient = (product * int(hessenberg[size - step - 1, size - 1])) % prime
            if coefficient:
                current[0:size - step] = (
                    current[0:size - step]
                    - coefficient * polynomials[size - step - 1][0:size - step]) % prime
        polynomials[size] = current % prime
    return [int(value) for value in polynomials[dimension] % prime]


def _coefficient_bound(matrix: list[list[int]], dimension: int) -> int:
    """Majorante de Hadamard para os coeficientes do polinomio caracteristico.

    O coeficiente de grau `n-k` e uma soma de menores principais `k x k`, logo
    `|c| <= C(n,k) · (sqrt(k) · B)^k` com `B` o maior modulo de entrada.
    """
    largest = max((abs(value) for row in matrix for value in row), default=0)
    if largest == 0:
        return 1
    bound = 1
    for size in range(1, dimension + 1):
        # Arredonde o majorante INTEIRO completo para cima. Arredondar
        # sqrt(k**k) para baixo antes de multiplicar por B**k perde
        # uma quantidade arbitrariamente grande quando B cresce.
        squared = (math.comb(dimension, size) ** 2
                   * size ** size * largest ** (2 * size))
        root = math.isqrt(squared)
        candidate = root + (root * root < squared)
        bound = max(bound, candidate)
    return bound


def _primes_above(start: int, count: int) -> list[int]:
    primes: list[int] = []
    candidate = start | 1
    while len(primes) < count:
        if _is_prime(candidate):
            primes.append(candidate)
        candidate += 2
    return primes


def _is_prime(number: int) -> bool:
    if number < 2:
        return False
    for small in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if number % small == 0:
            return number == small
    remainder, exponent = number - 1, 0
    while remainder % 2 == 0:
        remainder //= 2
        exponent += 1
    for witness in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        value = pow(witness, remainder, number)
        if value in (1, number - 1):
            continue
        for _ in range(exponent - 1):
            value = value * value % number
            if value == number - 1:
                break
        else:
            return False
    return True


def charpoly_exact(matrix: list[list[int]], dimension: int,
                   verbose: bool = True) -> list[int]:
    """Coeficientes inteiros de `q(y) = det(yI + M)`, indice = grau."""
    bound = _coefficient_bound(matrix, dimension)
    needed = 2 * bound + 1
    primes: list[int] = []
    modulus = 1
    pool = _primes_above(1 << 30, 8)
    index = 0
    while modulus <= needed:
        if index >= len(pool):
            pool.extend(_primes_above(pool[-1] + 2, 32))
        primes.append(pool[index])
        modulus *= pool[index]
        index += 1

    if verbose:
        print(f"  polinomio caracteristico exato: {len(primes)} primos de 31 bits "
              f"(cota de Hadamard: {bound.bit_length()} bits)", flush=True)

    residues: list[list[int]] = []
    started = time.time()
    for position, prime in enumerate(primes, start=1):
        residues.append(_charpoly_mod(matrix, dimension, prime))
        if verbose and (position % 50 == 0 or position == len(primes)):
            print(f"    primo {position}/{len(primes)} "
                  f"({time.time() - started:.0f}s)", flush=True)

    coefficients: list[int] = []
    for degree in range(dimension + 1):
        value = 0
        running = 1
        for prime, residue in zip(primes, residues):
            target = residue[degree]
            difference = (target - value) % prime
            value += running * (difference * pow(running % prime, prime - 2, prime) % prime)
            running *= prime
        if value > running // 2:
            value -= running
        coefficients.append(value)
    if coefficients[dimension] != 1:
        raise SystemExit("polinomio caracteristico nao saiu monico: erro de CRT")
    return coefficients


def bareiss_determinant(matrix: list[list[int]], dimension: int) -> int:
    """Determinante inteiro exato, sem fracoes (Bareiss). Verificacao lenta."""
    work = [row[:] for row in matrix]
    sign = 1
    previous = 1
    for step in range(dimension - 1):
        if work[step][step] == 0:
            for row in range(step + 1, dimension):
                if work[row][step]:
                    work[step], work[row] = work[row], work[step]
                    sign = -sign
                    break
            else:
                return 0
        pivot = work[step][step]
        for row in range(step + 1, dimension):
            factor = work[row][step]
            target, source = work[row], work[step]
            for column in range(step + 1, dimension):
                target[column] = (pivot * target[column]
                                  - factor * source[column]) // previous
            target[step] = 0
        previous = pivot
    return sign * work[dimension - 1][dimension - 1]


# --------------------------------------------------------------------------
# 3. Aritmetica de inteiros de Gauss
# --------------------------------------------------------------------------

def gauss_multiply(left: tuple[int, int], right: tuple[int, int]) -> tuple[int, int]:
    return (left[0] * right[0] - left[1] * right[1],
            left[0] * right[1] + left[1] * right[0])


def gauss_add(left: tuple[int, int], right: tuple[int, int]) -> tuple[int, int]:
    return (left[0] + right[0], left[1] + right[1])


def gauss_modulus_upper(value: tuple[int, int]) -> int:
    """`|w| <= |Re w| + |Im w|`, majorante inteiro sem raiz quadrada."""
    return abs(value[0]) + abs(value[1])


def gauss_modulus_lower(value: tuple[int, int]) -> int:
    """`|w| >= max(|Re w|, |Im w|)`, minorante inteiro sem raiz quadrada."""
    return max(abs(value[0]), abs(value[1]))


def horner(coefficients: list[int], point: tuple[int, int]) -> tuple[int, int]:
    value = (0, 0)
    for coefficient in reversed(coefficients):
        value = gauss_add(gauss_multiply(value, point), (coefficient, 0))
    return value


def taylor_shift(coefficients: list[int],
                 point: tuple[int, int]) -> list[tuple[int, int]]:
    """Coeficientes de `R(u + point)`, por divisoes sinteticas repetidas."""
    work: list[tuple[int, int]] = [(value, 0) for value in coefficients]
    degree = len(work) - 1
    shifted: list[tuple[int, int]] = []
    for order in range(degree + 1):
        carry = (0, 0)
        for index in range(degree - order, -1, -1):
            carry = gauss_add(gauss_multiply(carry, point), work[index])
            work[index] = carry
        shifted.append(work[0])
        work = work[1:degree - order + 1] if degree - order > 0 else []
        if not work:
            break
    while len(shifted) < degree + 1:
        shifted.append((0, 0))
    return shifted


# --------------------------------------------------------------------------
# 4. Contorno
# --------------------------------------------------------------------------

def rectangle_nodes(left: Fraction, right: Fraction, height: Fraction,
                    horizontal: int, vertical: int,
                    imag_center: Fraction = Fraction(0)
                    ) -> list[tuple[Fraction, Fraction]]:
    """Nos do retangulo `[left,right] x [-height,height]`, sentido anti-horario.

    Cada canto aparece uma unica vez. `vertical` deve ser IMPAR, para que
    `Im s = 0` caia no interior de uma aresta e nao sobre um no: e a regra 1 de
    `docs/WINDING_DATA.md` §2.6. A regra 2 (nos em `±h`) sai de graca, porque a
    subdivisao uniforme de um intervalo simetrico e simetrica.
    """
    if left >= right or height <= 0 or horizontal <= 0 or vertical <= 0:
        raise SystemExit("retangulo exige left < right, height > 0 e subdivisoes positivas")
    if vertical % 2 == 0:
        raise SystemExit("subdivisao vertical tem de ser IMPAR "
                         "(regra 1 de docs/WINDING_DATA.md §2.6)")
    nodes: list[tuple[Fraction, Fraction]] = []
    width = right - left
    span = 2 * height
    for index in range(horizontal):
        nodes.append((left + width * index / horizontal, -height))
    for index in range(vertical):
        nodes.append((right, -height + span * index / vertical))
    for index in range(horizontal):
        nodes.append((right - width * index / horizontal, height))
    for index in range(vertical):
        nodes.append((left, height - span * index / vertical))
    seen = set(nodes)
    if len(seen) != len(nodes):
        raise SystemExit("contorno degenerado: pontos z repetidos")
    return [(real, imag + imag_center) for real, imag in nodes]


def integer_lattice_nodes(nodes: list[tuple[Fraction, Fraction]], exponent: int
                          ) -> tuple[int, list[tuple[int, int]]]:
    """y=2^-E*z; u=scale*y inteiro, inclusive para E positivo."""
    factor = Fraction(2) ** (-exponent)
    scaled = [(real * factor, imag * factor) for real, imag in nodes]
    scale = math.lcm(*(part.denominator for point in scaled for part in point))
    return scale, [(int(real * scale), int(imag * scale)) for real, imag in scaled]


# --------------------------------------------------------------------------
# 5. Arredondamento para fora
# --------------------------------------------------------------------------

def round_up(value: Fraction, bits: int) -> Fraction:
    """Menor racional curto >= `value`. Usado nos raios, onde arredondar para
    cima e sempre o lado seguro: a caixa so cresce."""
    if value <= 0:
        return Fraction(0)
    exponent = value.numerator.bit_length() - value.denominator.bit_length()
    grid = Fraction(2) ** (exponent - bits)
    return Fraction(math.ceil(value / grid)) * grid


def outward_box(real: Fraction, imag: Fraction, bits: int
                ) -> tuple[Fraction, Fraction, Fraction]:
    """Centro e raio racionais curtos cuja caixa contem `real + i·imag`."""
    magnitude = max(abs(real), abs(imag))
    if magnitude == 0:
        return Fraction(0), Fraction(0), Fraction(0)
    exponent = magnitude.numerator.bit_length() - magnitude.denominator.bit_length()
    grid = Fraction(2) ** (exponent - bits)
    centre_real = Fraction(math.floor(real / grid)) * grid
    centre_imag = Fraction(math.floor(imag / grid)) * grid
    return centre_real + grid / 2, centre_imag + grid / 2, grid


# --------------------------------------------------------------------------
# 6. Classificacao (replica do Lean, so para o relatorio)
# --------------------------------------------------------------------------

def classify(start: tuple[Fraction, ...], end: tuple[Fraction, ...],
             edge: tuple[Fraction, ...]) -> int | None:
    if edge[0] > 0:
        return 0
    if edge[2] > 0 or edge[3] < 0:
        return 0
    if edge[1] < 0:
        if start[2] > 0 and end[2] > 0:
            return 0
        if start[3] < 0 and end[3] < 0:
            return 0
        if start[2] > 0 and end[3] < 0:
            return 1
        if start[3] < 0 and end[2] > 0:
            return -1
        return None
    return None


# --------------------------------------------------------------------------
# 7. Programa
# --------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("matrix", type=Path)
    parser.add_argument("--sigma0", default="1/8",
                        help="borda esquerda do retangulo em Re(s)")
    parser.add_argument("--radius", default="1",
                        help="borda direita do retangulo em Re(s)")
    parser.add_argument("--height", default="1/4",
                        help="meia-altura em Im(s); a faixa e [-h, h]")
    parser.add_argument("--imag-center", default="0",
                        help="centro imaginario; use 1/2 para o contorno reduzido de B")
    parser.add_argument("--horizontal", type=int, default=64)
    parser.add_argument("--vertical", type=int, default=33,
                        help="numero IMPAR de subdivisoes nos lados verticais")
    parser.add_argument("--bits", type=int, default=64,
                        help="bits significativos do arredondamento para fora")
    parser.add_argument("--name", default="contourA")
    parser.add_argument("--label", default="4x12")
    parser.add_argument("--out", type=Path)
    parser.add_argument("--cache", type=Path,
                        help="arquivo onde guardar/reler os coeficientes "
                             "inteiros do polinomio caracteristico")
    parser.add_argument("--verify", action="store_true",
                        help="confere o polinomio contra um determinante de "
                             "Bareiss independente (lento: ~42 s por ponto)")
    arguments = parser.parse_args()

    dimension, matrix, exponent = read_spectral_matrix(arguments.matrix)
    print(f"matriz: n={dimension}, A = M·2^({exponent}), "
          f"{max(abs(v).bit_length() for row in matrix for v in row)} bits no maior coeficiente")

    started = time.time()
    if arguments.cache is not None and arguments.cache.exists():
        charpoly = [int(token) for token in
                    arguments.cache.read_text(encoding="ascii").split()]
        if len(charpoly) != dimension + 1 or charpoly[dimension] != 1:
            raise SystemExit(f"{arguments.cache}: cache incompativel")
        print(f"  polinomio caracteristico relido de {arguments.cache}")
    else:
        charpoly = charpoly_exact(matrix, dimension)
        if arguments.cache is not None:
            arguments.cache.write_text(
                "\n".join(str(value) for value in charpoly) + "\n",
                encoding="ascii")
    print(f"  pronto em {time.time() - started:.0f}s; "
          f"maior coeficiente com {max(abs(c).bit_length() for c in charpoly)} bits")

    if arguments.verify:
        for probe in (1, 1 << 71):
            shifted = [[matrix[i][j] + (probe if i == j else 0)
                        for j in range(dimension)] for i in range(dimension)]
            started = time.time()
            direct = bareiss_determinant(shifted, dimension)
            through = horner(charpoly, (probe, 0))[0]
            status = "OK" if direct == through else "DIVERGENCIA"
            print(f"  verificacao Bareiss em y={probe}: {status} "
                  f"({time.time() - started:.0f}s)")
            if direct != through:
                return 1

    # Diagnostico em ponto flutuante: SO para escolher parametros e relatar.
    float_matrix = np.array([[float(value) for value in row] for row in matrix])
    float_matrix *= 2.0 ** exponent
    roots = -np.linalg.eigvals(float_matrix)
    left = Fraction(arguments.sigma0)
    right = Fraction(arguments.radius)
    height = Fraction(arguments.height)
    imag_center = Fraction(arguments.imag_center)
    inside = [root for root in roots
              if float(left) < root.real < float(right)
              and float(imag_center - height) < root.imag < float(imag_center + height)]
    print(f"diagnostico (ponto flutuante, NAO entra no certificado): "
          f"{len(inside)} raizes dentro do retangulo")
    for root in sorted(inside, key=lambda value: -value.real):
        print(f"    s = {root.real:+.9f}{root.imag:+.9f}i")

    nodes = rectangle_nodes(left, right, height,
                            arguments.horizontal, arguments.vertical, imag_center)
    print(f"contorno: {len(nodes)} nos "
          f"({arguments.horizontal} por lado horizontal, "
          f"{arguments.vertical} por lado vertical, impar)")

    # Denominador comum: leva todos os `y = 2^-E · z` a inteiros de Gauss.
    scale, integer_nodes = integer_lattice_nodes(nodes, exponent)

    # Polinomio homogeneizado em `u = scale · 2^-E · z`: coeficientes inteiros.
    homogeneous = [charpoly[degree] * scale ** (dimension - degree)
                   for degree in range(dimension + 1)]

    started = time.time()
    values = [horner(homogeneous, point) for point in integer_nodes]
    print(f"valores exatos nos {len(values)} nos em {time.time() - started:.1f}s; "
          f"maior valor com {max(gauss_modulus_upper(v).bit_length() for v in values)} bits")
    for index, value in enumerate(values):
        if value == (0, 0):
            raise SystemExit(f"no {index}: p(z) = 0 exatamente; o contorno "
                             "passa por uma raiz do truncamento")

    reference = values[0]
    reference_lower = gauss_modulus_lower(reference)

    started = time.time()
    deltas: list[Fraction] = []
    local_ratio: list[Fraction] = []
    jacobi_report: list[tuple[Fraction, Fraction, bool]] = []
    for index, point in enumerate(integer_nodes):
        following = integer_nodes[(index + 1) % len(integer_nodes)]
        length = max(abs(following[0] - point[0]), abs(following[1] - point[1]))
        coefficients = taylor_shift(homogeneous, point)
        variation = 0
        derivative = 0
        power = 1
        for order in range(1, dimension + 1):
            magnitude = gauss_modulus_upper(coefficients[order])
            if magnitude == 0:
                power *= length
                continue
            power_here = power * length
            variation += magnitude * power_here
            derivative += order * magnitude * power
            power = power_here
        constant = gauss_modulus_lower(coefficients[0])
        deltas.append(round_up(Fraction(variation, reference_lower),
                               arguments.bits))
        local_ratio.append(Fraction(variation, constant) if constant else Fraction(-1))
        lower = constant - variation
        if lower > 0:
            gain = Fraction(derivative, lower)
            product = gain * length
            jacobi_report.append((product,
                                  product * Fraction(gauss_modulus_upper(coefficients[0]),
                                                     reference_lower) / (1 - product)
                                  if product < 1 else Fraction(-1),
                                  product < Fraction(1, 2)))
        else:
            jacobi_report.append((Fraction(-1), Fraction(-1), False))
    print(f"bounds de aresta (Taylor exato + Jacobi) em {time.time() - started:.1f}s")

    boxes: list[tuple[Fraction, Fraction, Fraction]] = []
    for value in values:
        denominator = reference[0] ** 2 + reference[1] ** 2
        real = Fraction(value[0] * reference[0] + value[1] * reference[1], denominator)
        imag = Fraction(value[1] * reference[0] - value[0] * reference[1], denominator)
        boxes.append(outward_box(real, imag, arguments.bits))

    increments: list[int | None] = []
    for index in range(len(values)):
        following = (index + 1) % len(values)
        start_real, start_imag, start_radius = boxes[index]
        end_real, end_imag, end_radius = boxes[following]
        start_box = (start_real - start_radius, start_real + start_radius,
                     start_imag - start_radius, start_imag + start_radius)
        end_box = (end_real - end_radius, end_real + end_radius,
                   end_imag - end_radius, end_imag + end_radius)
        inflation = deltas[index]
        edge_box = (min(start_box[0], end_box[0]) - inflation,
                    max(start_box[1], end_box[1]) + inflation,
                    min(start_box[2], end_box[2]) - inflation,
                    max(start_box[3], end_box[3]) + inflation)
        increments.append(classify(start_box, end_box, edge_box))

    classified = [value for value in increments if value is not None]
    print(f"classificacao: {len(classified)}/{len(increments)} arestas classificam")
    print(f"delta/|p(z_k)| (Taylor, quantidade que decide a classificacao): "
          f"maior {float(max(local_ratio)):.3e}, menor {float(min(local_ratio)):.3e}")
    valid_products = [product for product, _, _ in jacobi_report if product >= 0]
    if valid_products:
        print(f"Jacobi K·L: maior {float(max(valid_products)):.3e}; "
              f"arestas com K·L < 1/2: "
              f"{sum(1 for _, _, ok in jacobi_report if ok)}/{len(jacobi_report)}")
    if len(classified) != len(increments):
        failed = [index for index, value in enumerate(increments) if value is None]
        print(f"ARESTAS QUE NAO CLASSIFICAM: {failed[:20]}"
              f"{' ...' if len(failed) > 20 else ''}")
        print("  refine o contorno (mais subdivisoes) ou aumente --bits.")
        return 2
    total = sum(classified)
    print(f"winding do truncamento sobre este contorno: {total}")

    if arguments.out is None:
        return 0
    lines = [HEADER_WARNING.format(label=arguments.label),
             f"name {arguments.name}",
             f"target {total}",
             ""]
    for index, ((real, imag, radius), (point_real, point_imag)) in enumerate(
            zip(boxes, nodes)):
        lines.append(f"node {real} {imag} {radius}")
        lines.append(f"point {point_real} {point_imag}")
        start_box = (real - radius, real + radius, imag - radius, imag + radius)
        following = (index + 1) % len(boxes)
        next_real, next_imag, next_radius = boxes[following]
        end_box = (next_real - next_radius, next_real + next_radius,
                   next_imag - next_radius, next_imag + next_radius)
        inflation = deltas[index]
        lines.append(f"edge {min(start_box[0], end_box[0]) - inflation} "
                     f"{max(start_box[1], end_box[1]) + inflation} "
                     f"{min(start_box[2], end_box[2]) - inflation} "
                     f"{max(start_box[3], end_box[3]) + inflation}")
    arguments.out.write_text("\n".join(lines) + "\n", encoding="ascii")
    print(f"escrito: {arguments.out} ({arguments.out.stat().st_size} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
