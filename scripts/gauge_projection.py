#!/usr/bin/env python3
"""O modo `lambda~0,61` e uma direcao de gauge?

DIAGNOSTICO EM PONTO FLUTUANTE. Nao e prova.

Satisfazer as constraints nao atesta ser fisico: o modo de fase tambem as
satisfaz, e e gauge conhecido. O teste certo, portanto, nao e o residuo -- e
projetar os autovetores sobre as direcoes de gauge construidas de forma
INDEPENDENTE do espectro, a partir do fundo `omega*` lido de `RefA.dat`.

Construcao da direcao de fase (translacao em tau), fixada pelo fonte de RT:

* `Field_U_0_m_times_im_unit` multiplica o coeficiente de indice de Fourier `m`
  por `i*m`: e o `d/dtau` do codigo. Logo a direcao de fase e
  `(i*m) * omega*_{d,m,n}`, um fator global irrelevante para sobreposicao.
* `IndexDOF.h`: `_Re_ = 0`, `_Im_ = 1`; o endereco de um DOF e `(d, k, m, n)`.
* Se `omega*_{d,m,n} = a + i b`, entao `(i m)(a + i b) = -m b + i (m a)`, de modo
  que o DOF `_Re_` recebe `-m b` e o DOF `_Im_` recebe `+m a`.

A FAMILIA DE REPARAMETRIZACAO (docs/MODE_DISCRIMINATION.md §9.6.1). Sob
`tau -> tau + eps f(tau)` o fundo varia por `f d/dtau omega* + f' W omega*`,
com `W` a matriz (desconhecida aqui) dos pesos tensoriais. Com `f = e^{s tau}`
isso da `h = d/dtau omega* + s W omega*`: o termo dominante NAO depende de `s`,
de modo que a fase (`s = 0`) e o tempo de acumulacao (`lambda = 1`, isto e
`s = 0,274`) tem a MESMA direcao dominante. Duas consequencias operacionais:

* projetar sobre `d/dtau omega*` e cego a `s` por construcao -- alta
  sobreposicao diz "vive no setor de reparametrizacao", nunca "este `s` e
  admissivel" (isso e S4);
* o teste da leitura "gauge com correcoes de ordem `s`" e projetar sobre o
  ESPACO que contem toda escolha possivel de `W`:
  `G = span{d/dtau omega*} + span{E_{d<-e} omega*_e}`, onde `E_{d<-e}` poe a
  componente `e` do fundo no slot `d`. As combinacoes proibidas por paridade
  caem sozinhas: os enderecos correspondentes nao existem na matriz.

Limitacao, registrada em docs/MODE_DISCRIMINATION.md §9 antes de rodar: em
`s = 0` a direcao de gauge da fase E exatamente `d/dtau omega*`, o que torna a
calibracao legitima; em `s != 0` a direcao verdadeira e `d/dtau omega*` mais
termos de ordem `s`. Para `s != 0` a sobreposicao com o vetor unico e um PROXY;
a projecao sobre `G` e uma cota SUPERIOR do carater gauge, porque `G` contem
mais do que as reparametrizacoes admissiveis.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import eig

sys.path.insert(0, str(Path(__file__).resolve().parent))
from analyze_spectrum import K_RT, read_constraints, read_matrix  # noqa: E402

PAIR = re.compile(r"\((-?\d+),(-?\d+)\)")
GAUGED_DOF = (3, 0, 1, 0)  # o unico DOF que `Field_set_VGauged` de RT zera


def read_mu(path: Path) -> float:
    """`mu` de RT, o diadico EXATO no cabecalho `DyadicQ` de `RefA.dat`.

    O arquivo comeca com `mu_MultiField / DyadicQ / TwoExp e / coeff c`, de modo
    que `mu = c * 2^e`. `SPECTRAL_PROBLEM.md` §1 diz que esse parametro controla
    a reescala dos nulos `u_±`; ele nao e o periodo de eco.
    """
    head = path.read_text(encoding="ascii")[:400]
    exponent = re.search(r"TwoExp\s+(-?\d+)", head)
    coefficient = re.search(r"coeff\s+(-?\d+)", head)
    if exponent is None or coefficient is None:
        raise SystemExit(f"{path}: cabecalho sem `mu`")
    return float(int(coefficient.group(1))) * 2.0 ** int(exponent.group(1))


def read_background(path: Path) -> tuple[list[np.ndarray], list[tuple[int, int]]]:
    """Le `RefA.dat`. Devolve, por componente, a matriz complexa de
    coeficientes `omega*_{m,n} = 2^TwoExp * (Re + i Im)` e o par `(off_m, off_n)`.
    """
    text = path.read_text(encoding="ascii")
    # Cuidado: "mu_MultiField" e "MultiField" contem "Field" como substring, de
    # modo que um split ingenuo quebra no lugar errado. Aqui o corte e feito
    # apenas em linhas que sejam exatamente "Field".
    lines = text.splitlines(keepends=True)
    starts = [index for index, line in enumerate(lines)
              if line.strip() == "Field"]
    blocks = ["".join(lines[start + 1:end])
              for start, end in zip(starts, starts[1:] + [len(lines)])]
    fields: list[np.ndarray] = []
    offsets: list[tuple[int, int]] = []
    for block in blocks:
        header = {}
        for key in ("TwoExp", "off_m", "off_n", "num_m", "num_n"):
            match = re.search(rf"{key}\s+(-?\d+)", block)
            if match is None:
                raise SystemExit(f"{path}: cabecalho sem '{key}'")
            header[key] = int(match.group(1))
        pairs = PAIR.findall(block)
        expected = header["num_m"] * header["num_n"]
        if len(pairs) != expected:
            raise SystemExit(f"{path}: esperados {expected} coeficientes, "
                             f"lidos {len(pairs)}")
        scale = 2.0 ** header["TwoExp"]
        data = np.empty((header["num_m"], header["num_n"]), dtype=complex)
        position = 0
        for u in range(header["num_m"]):
            for v in range(header["num_n"]):
                real, imag = pairs[position]
                position += 1
                data[u, v] = complex(float(int(real)) * scale,
                                     float(int(imag)) * scale)
        fields.append(data)
        offsets.append((header["off_m"], header["off_n"]))
    if len(fields) != 4:
        raise SystemExit(f"{path}: esperados 4 campos, lidos {len(fields)}")
    return fields, offsets


def phase_direction(background: list[np.ndarray], offsets: list[tuple[int, int]],
                    addresses: list[tuple[int, int, int, int]]) -> np.ndarray:
    """`d/dtau omega*` no espaco de DOFs, na ordem dos enderecos da matriz."""
    vector = np.zeros(len(addresses), dtype=float)
    for index, (component, part, fourier, chebyshev) in enumerate(addresses):
        field = background[component]
        off_m, off_n = offsets[component]
        u, v = fourier - off_m, chebyshev - off_n
        if not (0 <= u < field.shape[0] and 0 <= v < field.shape[1]):
            continue
        value = field[u, v]
        vector[index] = (-fourier * value.imag if part == 0
                         else fourier * value.real)
    return vector


def embedding(background: list[np.ndarray], offsets: list[tuple[int, int]],
              addresses: list[tuple[int, int, int, int]],
              slot: int, source: int) -> np.ndarray:
    """`E_{slot<-source} omega*`: a componente `source` do fundo colocada no
    slot `slot` do espaco de DOFs. Enderecos que a paridade do setor nao admite
    simplesmente nao existem na lista, de modo que a restricao e automatica."""
    vector = np.zeros(len(addresses), dtype=float)
    field = background[source]
    off_m, off_n = offsets[source]
    for index, (component, part, fourier, chebyshev) in enumerate(addresses):
        if component != slot:
            continue
        u, v = fourier - off_m, chebyshev - off_n
        if not (0 <= u < field.shape[0] and 0 <= v < field.shape[1]):
            continue
        value = field[u, v]
        vector[index] = value.real if part == 0 else value.imag
    return vector


def orthonormal(columns: list[np.ndarray], tolerance: float = 1e-10) -> np.ndarray:
    """Base ortonormal do espaco gerado, com deteccao de posto."""
    usable = [column for column in columns if np.linalg.norm(column) > 0]
    if not usable:
        return np.zeros((0, 0))
    stacked = np.column_stack([column / np.linalg.norm(column)
                               for column in usable])
    left, values, _ = np.linalg.svd(stacked, full_matrices=False)
    rank = int(np.sum(values > tolerance * values[0]))
    return left[:, :rank]


def subspace_overlap(basis: np.ndarray, vector: np.ndarray) -> float:
    """`||Q^T v|| / ||v||`: cosseno do angulo entre `v` e o subespaco."""
    norm = np.linalg.norm(vector)
    if norm == 0 or basis.size == 0:
        return float("nan")
    return float(np.linalg.norm(basis.T @ vector) / norm)


def overlap(left: np.ndarray, right: np.ndarray) -> float:
    left_norm = np.linalg.norm(left)
    right_norm = np.linalg.norm(right)
    if left_norm == 0 or right_norm == 0:
        return float("nan")
    return float(abs(np.vdot(left, right)) / (left_norm * right_norm))


def gauge_fixed(matrix: np.ndarray,
                addresses: list[tuple[int, int, int, int]],
                reference: list[float], mu: float) -> None:
    """Apaga linha E coluna do DOF que `Field_set_VGauged` zera, e reexamina.

    Diagnostico com controle embutido: o modo de FASE tem de sair do espectro.
    Se ele nao sair, esta reducao nao realiza a fixacao de gauge e nada se
    conclui dela (docs/MODE_DISCRIMINATION.md §9.6.4, F2).
    """
    if GAUGED_DOF not in addresses:
        print("  o DOF de `VGauged` nao existe neste truncamento; teste pulado")
        return
    position = addresses.index(GAUGED_DOF)
    keep = [index for index in range(matrix.shape[0]) if index != position]
    reduced = matrix[np.ix_(keep, keep)]
    exponents = -eig(reduced, right=False)
    reals = sorted([value.real for value in exponents
                    if abs(value.imag) < 1e-6 and value.real > -0.05],
                   reverse=True)
    print(f"  fixacao de gauge (DOF {GAUGED_DOF} removido, dimensao "
          f"{reduced.shape[0]}): raizes reais em Re s > -0,05")
    for value in reals:
        near = min(reference, key=lambda other: abs(other - value))
        print(f"    s = {value:+.9f}   (mais proximo do original "
              f"{near:+.9f}, distancia {abs(near - value):.2e})")
    for original in reference:
        distance = min((abs(value - original) for value in reals),
                       default=float("inf"))
        verdict = "SOBREVIVE" if distance < 5e-3 else "SAIU (ou moveu muito)"
        print(f"    modo original {original:+.9f}: {verdict} "
              f"(distancia {distance:.2e})")
    close = int(np.sum(np.abs(exponents - mu) < 1e-3))
    print(f"    autovalores a menos de 1e-3 de +mu apos a fixacao: {close}")


def analyse(path: Path, background: list[np.ndarray],
            offsets: list[tuple[int, int]], want_constraints: bool,
            mu: float) -> dict:
    matrix, addresses = read_matrix(path)
    dimension = matrix.shape[0]
    direction = phase_direction(background, offsets, addresses)
    generic = 1.0 / np.sqrt(dimension)

    # Checagem estrutural independente do espectro: a condicao `VGauged` de RT
    # zera exatamente UM DOF, `Re` da componente 3 em `(m=1, n=0)`. Se a direcao
    # de fase construida for mesmo a direcao de gauge que essa condicao remove,
    # ela tem de ter componente nao nula ali.
    gauge_dof = [index for index, address in enumerate(addresses)
                 if address == (3, 0, 1, 0)]
    marker = (abs(direction[gauge_dof[0]]) / np.linalg.norm(direction)
              if gauge_dof else float("nan"))

    # A familia de reparametrizacao: `G_diag` com `W` diagonal (o caso tensorial
    # usual, um peso por componente) e `G_full` com mistura arbitraria, que e
    # uma cota SUPERIOR e paga o preco de um nivel generico maior.
    diagonal = [embedding(background, offsets, addresses, d, d) for d in range(4)]
    mixing = [embedding(background, offsets, addresses, d, e)
              for d in range(4) for e in range(4)]
    basis_weights = orthonormal(diagonal)
    basis_diag = orthonormal([direction] + diagonal)
    basis_full = orthonormal([direction] + mixing)

    constant = linear = None
    if want_constraints:
        constraint_path = Path(f"{path}.constraints")
        if constraint_path.exists():
            constant, linear = read_constraints(constraint_path, dimension)

    values, vectors = eig(matrix)
    left_values, left_vectors = eig(matrix.T)
    exponents = -values
    left_exponents = -left_values
    reals = sorted([value for value in exponents
                    if abs(value.imag) < 1e-6 and value.real > -0.05],
                   key=lambda value: -value.real)

    unit = direction / np.linalg.norm(direction)
    modes = []
    for target in reals:
        index = int(np.argmin(np.abs(exponents - target)))
        raw = vectors[:, index]
        if np.linalg.norm(raw.imag) > 1e-8 * np.linalg.norm(raw.real):
            print(f"  AVISO: autovetor de s={target.real:+.6f} tem parte "
                  f"imaginaria nao desprezivel; usando a parte real")
        vector = np.real(raw)
        partner = int(np.argmin(np.abs(left_exponents - target)))
        left_vector = left_vectors[:, partner]
        overlap_lr = abs(np.dot(left_vector / np.linalg.norm(left_vector),
                                raw / np.linalg.norm(raw)))
        unit_vector = vector / np.linalg.norm(vector)
        residual = unit_vector - np.dot(unit_vector, unit) * unit
        record = {
            "s": target.real,
            "mu_distance": abs(target.real - mu),
            "lambda": 2.0 * np.pi / K_RT * target.real,
            "vector": vector,
            "residual": residual,
            "cos_phase": overlap(direction, vector),
            "cos_diag": subspace_overlap(basis_diag, vector),
            "cos_full": subspace_overlap(basis_full, vector),
            "cos_weights_of_residual": subspace_overlap(basis_weights, residual),
            "condition": 1.0 / overlap_lr if overlap_lr > 0 else np.inf,
            "constraint": float("nan"),
        }
        if constant is not None:
            operator = constant + target.real * linear
            record["constraint"] = float(
                np.linalg.norm(operator @ unit_vector))
        modes.append(record)
    near_mu = int(np.sum(np.abs(exponents - mu) < 1e-3))
    return {
        "name": path.name,
        "dimension": dimension,
        "matrix": matrix,
        "addresses": addresses,
        "near_mu": near_mu,
        "generic": generic,
        "marker": marker,
        "rank_diag": basis_diag.shape[1] if basis_diag.size else 0,
        "rank_full": basis_full.shape[1] if basis_full.size else 0,
        "modes": modes,
    }


def report(result: dict) -> None:
    generic = result["generic"]
    k_diag, k_full = result["rank_diag"], result["rank_full"]
    print()
    print(f"=== {result['name']}: dimensao {result['dimension']}, "
          f"sobreposicao generica ~{generic:.4f} ===")
    print(f"peso da direcao de fase no DOF que `VGauged` remove "
          f"(comp 3, Re, m=1, n=0): {result['marker']:.4f}")
    print(f"subespaco de reparametrizacao: dim {k_diag} (diagonal) e {k_full} "
          f"(mistura); niveis genericos {np.sqrt(k_diag/result['dimension']):.4f}"
          f" e {np.sqrt(k_full/result['dimension']):.4f}")
    print(f"autovalores de L_A(0) a menos de 1e-3 de -mu: {result['near_mu']}")
    print(f"{'s':>14} {'lambda':>10} {'|cos| fase':>11} {'x gen':>7} "
          f"{'|cos| G_diag':>13} {'|cos| G_full':>13} {'resid->W':>9} "
          f"{'constr':>10} {'|s-mu|':>10} {'cond':>9}")
    for mode in result["modes"]:
        print(f"{mode['s']:+14.9f} {mode['lambda']:+10.6f} "
              f"{mode['cos_phase']:11.4f} {mode['cos_phase']/generic:7.1f} "
              f"{mode['cos_diag']:13.4f} {mode['cos_full']:13.4f} "
              f"{mode['cos_weights_of_residual']:9.4f} "
              f"{mode['constraint']:10.3e} {mode['mu_distance']:10.3e} "
              f"{mode['condition']:9.2e}")
    print("  sobreposicao mutua dos autovetores (|cos|):")
    print("        " + " ".join(f"{mode['s']:+9.4f}" for mode in result["modes"]))
    for mode in result["modes"]:
        row = " ".join(f"{overlap(mode['vector'], other['vector']):9.4f}"
                       for other in result["modes"])
        print(f"{mode['s']:+7.4f} {row}")


def match(first: dict, second: dict) -> None:
    """Pareamento entre truncamentos por SOBREPOSICAO, nao por ordenacao.

    A §9.1 mostrou que no 4x12 o rotulo "fase" estava trocado; parear por ordem
    de `s` propaga esse tipo de erro em silencio. Aqui cada modo do truncamento
    menor e casado com o modo do maior que maximiza `|cos|` nos DOFs comuns, e o
    pareamento escolhido e impresso junto com o segundo melhor, para que se veja
    se a escolha foi folgada ou apertada.
    """
    shared = sorted(set(first["addresses"]) & set(second["addresses"]))
    index_a = {address: position
               for position, address in enumerate(first["addresses"])}
    index_b = {address: position
               for position, address in enumerate(second["addresses"])}
    take_a = [index_a[address] for address in shared]
    take_b = [index_b[address] for address in shared]

    print()
    print(f"=== {first['name']} -> {second['name']}: pareamento por |cos| "
          f"em {len(shared)} DOFs comuns ===")
    print(f"{'s (menor)':>12} {'s (maior)':>12} {'|cos|':>8} {'2o melhor':>10} "
          f"{'ortogonal ao gauge':>19}")
    used: set[int] = set()
    for mode in first["modes"]:
        left = mode["vector"][take_a]
        scores = [(overlap(left, other["vector"][take_b]), position)
                  for position, other in enumerate(second["modes"])]
        scores.sort(reverse=True)
        if not scores:
            continue
        best, position = scores[0]
        runner = scores[1][0] if len(scores) > 1 else float("nan")
        partner = second["modes"][position]
        tag = "  (JA USADO)" if position in used else ""
        used.add(position)
        residual_cos = overlap(mode["residual"][take_a],
                               partner["residual"][take_b])
        print(f"{mode['s']:+12.6f} {partner['s']:+12.6f} {best:8.4f} "
              f"{runner:10.4f} {residual_cos:19.4f}{tag}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("matrices", type=Path, nargs="+")
    parser.add_argument("--background", type=Path,
                        default=Path(".cache/rt-1203.3766v1/sourcecode/RefA.dat"))
    parser.add_argument("--constraints", action="store_true",
                        help="tambem reporta ||C(s)v||/||v|| por modo")
    parser.add_argument("--gauge-fixed", action="store_true",
                        help="apaga o DOF que `VGauged` zera e reexamina (F2)")
    arguments = parser.parse_args()

    background, offsets = read_background(arguments.background)
    mu = read_mu(arguments.background)
    print("AVISO: diagnostico em ponto flutuante; NAO E PROVA.")
    print(f"fundo: {arguments.background}, 4 componentes "
          f"{'x'.join(str(dimension) for dimension in background[0].shape)}")
    print(f"mu de RT (diadico exato do cabecalho): {mu!r}")

    results = []
    for path in arguments.matrices:
        result = analyse(path, background, offsets, arguments.constraints, mu)
        report(result)
        if arguments.gauge_fixed:
            gauge_fixed(result["matrix"], result["addresses"],
                        [mode["s"] for mode in result["modes"]], mu)
        results.append(result)

    for position in range(len(results) - 1):
        match(results[position], results[position + 1])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
