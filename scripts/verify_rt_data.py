#!/usr/bin/env python3
"""Verifica exatamente os dois certificados INV distribuídos por RT.

O programa não confia nos marcadores ``TRUE`` gravados nos arquivos: ele
recalcula as duas desigualdades com inteiros e verifica a cobertura de todos
os graus de liberdade permitidos pelas paridades Fourier--Chebyshev.
"""

from __future__ import annotations

import argparse
import hashlib
import re
from fractions import Fraction
from pathlib import Path

LINE = re.compile(
    rb"([0-9]+)_([01])_([0-9]+)_([0-9]+) "
    rb"(-?[0-9]+)\*2\^\((-?[0-9]+)\) (TRUE|FALSE) "
    rb"(-?[0-9]+)\*2\^\((-?[0-9]+)\) (TRUE|FALSE) *\n"
)

EXPECTED = {
    "INV_0_0_250_750.dat": (
        654_375,
        "df3c2b1aa022ed6e1a5998822e7f014819e669d5f4ff49f9345740cee84f1efc",
        4,
    ),
    "SHARP_INV_0_0_250_750.dat": (
        280_125,
        "d198799b93d214a2b3e7c51b170e613667bae6fb025231c1c56a566383cb91b1",
        2,
    ),
}


def dyadic(coefficient: int, exponent: int) -> Fraction:
    if exponent >= 0:
        return Fraction(coefficient * (1 << exponent), 1)
    return Fraction(coefficient, 1 << (-exponent))


def gothic_a(fourier: int, chebyshev: int) -> Fraction:
    if fourier < 50 and chebyshev < 150:
        return Fraction(256)
    if fourier < 100 and chebyshev < 300:
        return Fraction(16)
    if fourier < 150 and chebyshev < 450:
        return Fraction(4)
    if fourier < 200 and chebyshev < 600:
        return Fraction(2)
    return Fraction(25, 16)


def valid_address(component: int, part: int, fourier: int, chebyshev: int, components: int) -> bool:
    if not (0 <= component < components and part in (0, 1)):
        return False
    if not (0 <= fourier < 250 and 0 <= chebyshev < 750):
        return False
    # Realidade da série: o modo Fourier zero não tem parte imaginária.
    if fourier == 0 and part == 1:
        return False
    # componentes 0,...,components-1 são periódicas; no sistema principal a
    # componente 3 é antiperiódica.
    periodic = not (components == 4 and component == 3)
    if periodic != (fourier % 2 == 0):
        return False
    # A primeira componente é ímpar em xi.
    if component == 0 and chebyshev % 2 == 0:
        return False
    return True


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify(path: Path) -> None:
    try:
        expected_lines, expected_hash, components = EXPECTED[path.name]
    except KeyError as error:
        raise ValueError(f"nome de certificado não reconhecido: {path.name}") from error
    actual_hash = sha256(path)
    if actual_hash != expected_hash:
        raise ValueError(f"SHA-256 incorreto para {path.name}: {actual_hash}")

    seen: set[tuple[int, int, int, int]] = set()
    with path.open("rb") as stream:
        for line_number, raw_line in enumerate(stream, 1):
            match = LINE.fullmatch(raw_line)
            if match is None:
                raise ValueError(f"linha {line_number}: formato inválido")
            d, k, m, n, c1, e1, flag1, c2, e2, flag2 = match.groups()
            address = tuple(map(int, (d, k, m, n)))
            if address in seen:
                raise ValueError(f"linha {line_number}: endereço duplicado {address}")
            if not valid_address(*address, components):
                raise ValueError(f"linha {line_number}: grau de liberdade inválido {address}")
            seen.add(address)

            first = dyadic(int(c1), int(e1))
            second = dyadic(int(c2), int(e2))
            first_true = first <= Fraction(1, 1 << 16)
            second_true = second <= gothic_a(address[2], address[3])
            if (flag1 == b"TRUE") != first_true:
                raise ValueError(f"linha {line_number}: primeiro marcador falso")
            if (flag2 == b"TRUE") != second_true:
                raise ValueError(f"linha {line_number}: segundo marcador falso")
            if not (first_true and second_true):
                raise ValueError(f"linha {line_number}: desigualdade não satisfeita")

    if len(seen) != expected_lines:
        raise ValueError(f"cobertura incompleta: {len(seen)} de {expected_lines}")
    print(f"OK {path.name}: {len(seen)} linhas, SHA-256 e desigualdades exatas")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("files", type=Path, nargs="+")
    args = parser.parse_args()
    for path in args.files:
        verify(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

