#!/usr/bin/env python3
"""As constraints discriminam as raizes instaveis do pencil truncado?

DIAGNOSTICO EM PONTO FLUTUANTE. Nada aqui e prova. O objetivo e responder,
com numeros, a uma pergunta que orienta as obrigacoes S3 (propagacao das
constraints) e S4 (gauge): das raizes que o winding conta em `Gamma_A`, quais
carregam um autovetor que satisfaz as constraints?

Tres medidas, por modo:

1. `residuo` -- `||C(s) v|| / ||v||`, a mesma quantidade que
   `analyze_spectrum.py` imprime para o candidato. Depende da escala das linhas
   de `C`, logo so e comparavel entre modos do MESMO truncamento.
2. `q` -- o residuo dividido pelo que um vetor generico produziria,
   `||C(s)||_F / sqrt(n)`. Adimensional: `q ~ 1` significa "o autovetor viola as
   constraints tanto quanto um vetor aleatorio"; `q << 1` significa que ele e
   especial.
3. `sigma_min` do sistema empilhado `[L_A(s); C(s)]`, normalizado por
   `||C(s)||_2`. Esta e a medida robusta: ela pergunta se existe ALGUM vetor
   quase-nulo para o pencil E para as constraints ao mesmo tempo, sem depender
   de o autovetor estar bem condicionado.

Tambem imprime o numero de condicionamento do autovalor, porque um residuo
grande num autovalor mal condicionado nao significa nada.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import eig, svd, svdvals

sys.path.insert(0, str(Path(__file__).resolve().parent))
from analyze_spectrum import K_RT, read_constraints, read_matrix  # noqa: E402


def constrained_gap(matrix: np.ndarray, constant: np.ndarray,
                    linear: np.ndarray, s: complex,
                    tolerance: float = 1e-10) -> tuple[float, int]:
    """`min ||(A+sI)v|| / ||A+sI||` sobre `v` unitario no NUCLEO de `C(s)`.

    Esta e a pergunta de S3 posta diretamente: existe vetor que anula o pencil
    *e* satisfaz as constraints? Ao contrario do residuo do autovetor, ela nao
    depende de o autovetor calculado estar bem condicionado, e ao contrario do
    `sigma_min` empilhado ela nao mistura as duas condicoes numa soma de
    quadrados: as constraints entram como vinculo exato (ate a tolerancia do
    posto numerico) e o pencil como quantidade minimizada.
    """
    dimension = matrix.shape[0]
    operator = constant + s * linear
    _, values, right = svd(operator)
    rank = int(np.sum(values > tolerance * values[0]))
    kernel = right[rank:].conj().T
    pencil = matrix + s * np.eye(dimension)
    projected = pencil @ kernel
    return float(svdvals(projected)[-1] / np.linalg.norm(pencil, 2)), kernel.shape[1]


def mode_table(matrix: np.ndarray, constant: np.ndarray, linear: np.ndarray,
               selected: list[complex], label: str) -> list[dict]:
    dimension = matrix.shape[0]
    values, right = eig(matrix)
    left_values, left = eig(matrix.T)
    exponents = -values
    left_exponents = -left_values

    rows: list[dict] = []
    for target in selected:
        index = int(np.argmin(np.abs(exponents - target)))
        s = exponents[index]
        vector = right[:, index]
        vector = vector / np.linalg.norm(vector)

        # Condicionamento do autovalor: 1/|w . x|, com `w` autovetor a ESQUERDA.
        # `eig(A.T)` devolve `w` com `w^T A = lambda w^T`, de modo que o produto
        # correto e o bilinear `w . x`, SEM conjugar (usar `vdot` conjugaria e
        # daria um numero sem sentido para modos complexos).
        partner = int(np.argmin(np.abs(left_exponents - s)))
        left_vector = left[:, partner]
        left_vector = left_vector / np.linalg.norm(left_vector)
        overlap = abs(np.dot(left_vector, vector))
        condition = 1.0 / overlap if overlap > 0 else np.inf

        operator = constant + s * linear
        residual = np.linalg.norm(operator @ vector)
        generic = np.linalg.norm(operator, "fro") / np.sqrt(dimension)
        stacked = np.vstack([matrix + s * np.eye(dimension), operator])
        smallest = svdvals(stacked)[-1] / np.linalg.norm(operator, 2)

        gap, kernel_dimension = constrained_gap(matrix, constant, linear, s)
        rows.append({
            "s": s,
            "gap": gap,
            "kernel": kernel_dimension,
            "residual": residual,
            "q": residual / generic if generic else np.inf,
            "stacked": smallest,
            "condition": condition,
            "label": label,
        })
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("matrix", type=Path)
    parser.add_argument("--scan", action="store_true",
                        help="varre s real e imprime sigma_min empilhado")
    arguments = parser.parse_args()

    matrix, _ = read_matrix(arguments.matrix)
    constraint_path = Path(f"{arguments.matrix}.constraints")
    if not constraint_path.exists():
        raise SystemExit(f"{constraint_path}: ausente")
    constant, linear = read_constraints(constraint_path, matrix.shape[0])

    print("AVISO: diagnostico em ponto flutuante; NAO E PROVA.")
    print(f"matriz {arguments.matrix.name}: dimensao {matrix.shape[0]}, "
          f"constraints {constant.shape[0]}x{constant.shape[1]}")

    values, _ = eig(matrix)
    exponents = -values
    unstable = sorted([value for value in exponents if value.real > 1e-9],
                      key=lambda value: -value.real)
    reals = [value for value in unstable if abs(value.imag) < 1e-6]
    complexes = [value for value in unstable if abs(value.imag) >= 1e-6]
    print(f"raizes com Re s > 0: {len(unstable)} "
          f"({len(reals)} reais, {len(complexes)} complexas)")

    # Controle interno: o modo de fase (o mais proximo de s=0) e um modo de
    # gauge GENUINO -- a fase da solucao DSS. Ele calibra o piso: qualquer
    # discriminacao tem de coloca-lo do lado "satisfaz as constraints".
    phase = min(exponents, key=abs)
    selected = reals + sorted(complexes, key=lambda value: -value.real)[:4]
    if not any(abs(value - phase) < 1e-12 for value in selected):
        selected = selected + [phase]
    rows = mode_table(matrix, constant, linear, selected, "")

    scale = 2.0 * np.pi / K_RT
    print()
    print(f"{'s':>26} {'lambda':>12} {'residuo':>12} {'q':>10} "
          f"{'sigma_min':>11} {'gap_S3':>11} {'cond':>9}")
    for row in rows:
        s = row["s"]
        tag = ("FASE (gauge conhecido)" if abs(s) < 0.05
               else "real" if abs(s.imag) < 1e-6 else "complexa")
        print(f"{s.real:+12.9f}{s.imag:+12.9f}i {scale*s.real:+12.6f} "
              f"{row['residual']:12.4e} {row['q']:10.4f} "
              f"{row['stacked']:11.4e} {row['gap']:11.4e} "
              f"{row['condition']:9.2e}  {tag}")

    # contexto: distribuicao de q sobre TODOS os modos
    everything = mode_table(matrix, constant, linear, list(exponents), "")
    qs = np.array([row["q"] for row in everything])
    print()
    print(f"distribuicao de q sobre os {len(qs)} modos: "
          f"min {qs.min():.4f}, mediana {np.median(qs):.4f}, max {qs.max():.4f}")
    for row in rows:
        rank = float((qs < row["q"]).mean())
        print(f"   s={row['s'].real:+.6f}{row['s'].imag:+.6f}i "
              f"esta no percentil {100*rank:5.1f} de q")

    # Linha de base em `s` GENERICO: pontos reais que ficam longe de qualquer
    # raiz. Sem esse afastamento a comparacao seria viciada, porque perto de uma
    # raiz o bloco do pencil ja e quase singular por si so.
    dimension = matrix.shape[0]
    probes = [value for value in np.linspace(0.05, 1.0, 40)
              if np.min(np.abs(exponents - value)) > 0.05]
    baseline = []
    gaps = []
    for probe in probes:
        operator = constant + probe * linear
        stacked = np.vstack([matrix + probe * np.eye(dimension), operator])
        baseline.append(svdvals(stacked)[-1] / np.linalg.norm(operator, 2))
        gaps.append(constrained_gap(matrix, constant, linear, probe)[0])
    print(f"linha de base sobre {len(probes)} pontos genericos "
          f"(a mais de 0,05 de qualquer raiz):")
    print(f"   sigma_min: {min(baseline):.4e} a {max(baseline):.4e} "
          f"(mediana {np.median(baseline):.4e})")
    print(f"   gap_S3:    {min(gaps):.4e} a {max(gaps):.4e} "
          f"(mediana {np.median(gaps):.4e}); dim ker C(s) = "
          f"{constrained_gap(matrix, constant, linear, float(probes[0]))[1]}"
          f" de {matrix.shape[0]}")

    if arguments.scan:
        print()
        print("varredura em s real: sigma_min empilhado normalizado")
        dimension = matrix.shape[0]
        for step in np.linspace(0.05, 1.0, 39):
            operator = constant + step * linear
            stacked = np.vstack([matrix + step * np.eye(dimension), operator])
            print(f"   s={step:5.3f}  {svdvals(stacked)[-1]/np.linalg.norm(operator,2):.6e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
