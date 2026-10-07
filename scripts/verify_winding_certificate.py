#!/usr/bin/env python3
"""Verifica o formato inteiro usado por WindingCertificate no Lean.

Este programa não calcula winding a partir de dados flutuantes: ele somente
confere um certificado já produzido por um estimador intervalar independente.
"""
from __future__ import annotations

import argparse
from pathlib import Path


def verify(path: Path, target: int) -> int:
    values = [int(line.strip()) for line in path.read_text().splitlines()
              if line.strip() and not line.lstrip().startswith("#")]
    if not values:
        raise ValueError("contorno vazio")
    invalid = [value for value in values if value not in (-1, 0, 1)]
    if invalid:
        raise ValueError(f"incrementos inválidos: {invalid}")
    total = sum(values)
    if total != target:
        raise ValueError(f"soma {total} diferente do alvo {target}")
    print(f"OK: {len(values)} arestas, winding={total}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("certificate", type=Path)
    parser.add_argument("target", type=int)
    args = parser.parse_args()
    return verify(args.certificate, args.target)


if __name__ == "__main__":
    raise SystemExit(main())
