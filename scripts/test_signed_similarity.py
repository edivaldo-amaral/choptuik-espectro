#!/usr/bin/env python3
"""Auditoria: o operador com sinal (signed_operator, base do certificado de S3a) e o
montador real validado (sharp_operator_builder) sao o MESMO operador, entrada a entrada:
K_sinal = T K_real T^-1, com T a mudanca de base explicita

    c_{+m} = x_(m, parte 0) + i x_(m, parte 1),   c_{-m} = x_(m, parte 0) - i x_(m, parte 1),
    c_0    = x_(0, parte 0),

(parte 0/1 = parte real/imaginaria do coeficiente em +m, convencao de mult() no montador
real). Nao usa o espectro: compara as matrizes."""
from __future__ import annotations

from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import sharp_operator_builder as sb
import signed_operator as so


def T_matriz(M, N):
    addr = sb.enderecos(M, N)
    pos_real = {a: i for i, a in enumerate(addr)}
    ms = so.modos(M)
    n1, n2 = so.radial(N)
    k = len(n1) + len(n2)
    loc = [(0, n) for n in n1] + [(1, n) for n in n2]     # (componente, n) na ordem local
    T = np.zeros((len(ms)*k, len(addr)), complex)
    for i, m in enumerate(ms):
        for l, (d, n) in enumerate(loc):
            linha = i*k + l
            if m == 0:
                T[linha, pos_real[(d, 0, 0, n)]] = 1
            else:
                T[linha, pos_real[(d, 0, abs(m), n)]] = 1
                T[linha, pos_real[(d, 1, abs(m), n)]] = 1j if m > 0 else -1j
    return T


def main():
    for M, N in ((8, 24), (12, 36), (16, 50)):
        Kr, _ = sb.operador(M, N)
        Ks, _ = so.operador(M, N)
        T = T_matriz(M, N)
        dif = Ks - T @ Kr @ np.linalg.inv(T)
        print(f'caixa {M}x{N} (dim {len(Kr)}): max |K_sinal - T K_real T^-1| = {np.abs(dif).max():.2e}'
              f'   (max |K| = {np.abs(Ks).max():.2e})')


if __name__ == '__main__':
    main()
