#!/usr/bin/env python3
"""Conta usos REAIS da tática `native_decide` em arquivos Lean.

Menções em comentários não contam: o que importa para o ledger de confiança é
a tática efetivamente elaborada, não o texto que fala sobre ela. Por isso os
comentários de linha (`--`) e de bloco aninhado (`/- ... -/`, incluindo `/-! -/`)
são removidos antes da contagem. Ver docs/KERNEL_TRUST.md.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

TOKEN = re.compile(r"\bnative_decide\b")


def strip_comments(source: str) -> str:
    """Remove comentários Lean preservando as quebras de linha (para numerar)."""
    out: list[str] = []
    index = 0
    depth = 0
    length = len(source)
    while index < length:
        if depth == 0:
            if source.startswith("--", index):
                end = source.find("\n", index)
                if end == -1:
                    break
                out.append("\n")
                index = end + 1
            elif source.startswith("/-", index):
                depth = 1
                index += 2
            else:
                out.append(source[index])
                index += 1
        else:
            if source.startswith("/-", index):
                depth += 1
                index += 2
            elif source.startswith("-/", index):
                depth -= 1
                index += 2
            else:
                if source[index] == "\n":
                    out.append("\n")
                index += 1
    return "".join(out)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    args = parser.parse_args()

    total = 0
    for path in sorted(args.root.rglob("*.lean")):
        if ".lake" in path.parts:
            continue
        stripped = strip_comments(path.read_text(encoding="utf-8"))
        for number, line in enumerate(stripped.splitlines(), 1):
            if TOKEN.search(line):
                total += 1
                print(f"  {path}:{number}", file=sys.stderr)
    print(total)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
