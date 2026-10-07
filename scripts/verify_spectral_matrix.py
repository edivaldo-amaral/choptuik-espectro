#!/usr/bin/env python3
"""Auditoria estrutural exata de uma matriz espectral RT em formato dyádico."""
from __future__ import annotations
import argparse
from pathlib import Path

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("matrix", type=Path)
    args = ap.parse_args()
    with args.matrix.open(encoding="ascii") as f:
        if f.readline().strip() != "CHOPTUIK_RT_SPECTRAL_MATRIX_V1":
            raise ValueError("magic inválido")
        key, raw = f.readline().split(); assert key == "dimension"
        n = int(raw); assert n > 0
        f.readline(); f.readline()
        assert f.readline().strip() == "BEGIN_ADDRESSES"
        addresses = []
        for i in range(n):
            fields = f.readline().split()
            if len(fields) != 5 or int(fields[0]) != i:
                raise ValueError(f"endereço inválido na posição {i}")
            addresses.append(tuple(map(int, fields[1:])))
        assert f.readline().strip() == "BEGIN_ENTRIES"
        seen: set[tuple[int, int]] = set()
        count = 0
        for lineno, line in enumerate(f, n + 8):
            fields = list(map(int, line.split()))
            if len(fields) == 3:
                row, col, coeff = fields; exponent = 0
            elif len(fields) == 4:
                row, col, coeff, exponent = fields
            else:
                raise ValueError(f"entrada inválida na linha {lineno}")
            if not (0 <= row < n and 0 <= col < n):
                raise ValueError(f"índice fora da matriz na linha {lineno}")
            if (row, col) in seen:
                raise ValueError(f"entrada duplicada na linha {lineno}")
            if coeff == 0:
                raise ValueError(f"coeficiente zero explícito na linha {lineno}")
            seen.add((row, col)); count += 1
    print(f"OK: dimensão={n}, endereços={len(addresses)}, entradas={count}; aritmética dyádica preservada")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
