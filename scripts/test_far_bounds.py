#!/usr/bin/env python3
"""Conferencia numerica das cotas externas usadas pelo nivel far de tile_perron:

  * cauda_Rinv(sigma, N)  >= ||R_m^-1 P_{n >= N}||      (pesada, radial)
  * Rinv_uniforme(M, s0)  >= sup_{|m| >= M} ||R_m^-1||
  * descida(J, Nc)        >= ||P_{n < J} R_m^-1 P_{n >= Nc}||
  * NORMA_B = 1,64967     >= ||B|| numa caixa (compressao do operador infinito)

A "verdade" e a norma numa truncagem radial grande (Nbig), que para R_m^-1 e exata no
bloco (triangular superior) e cresce para a do operador infinito com Nbig.
"""
from __future__ import annotations

from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fourier_tail_bound as ft
import tile_perron as tp
from signed_operator import Rinv, operador, modos, MU, K1


def pesado(m, s, N):
    md = tp.Modo(N)
    return (md.w[:, None]*Rinv(m, s, N))/md.w[None, :], md


def main():
    Nbig = 1500
    pior = 0.0
    print('cauda_Rinv (radial):')
    for s in (0.0, 0.2, 1.0 + 0.25j):
        for m in (0, 10, 40, -100):
            for N in (160, 400):
                R, md = pesado(m, s, Nbig)
                sig = s + 1j*m/2
                for passo, sel in ((2, md.n % 2 == 1), (1, None)):
                    # componente c1 (n impar) = primeira metade; c2 = segunda
                    n1 = len(md.n1)
                    blk = slice(0, n1) if passo == 2 else slice(n1, len(md.n))
                    Rc = R[blk, blk]; nc = md.n[blk]
                    verdade = np.linalg.norm(Rc[:, nc >= N], 2)
                    cota = ft.cauda_Rinv(sig, N, passo, MU)
                    pior = max(pior, verdade/cota)
                    print(f'  s {s!s:>10} m {m:+4d} N {N} passo {passo}: {verdade:.4f} <= {cota:.4f}')
    print('Rinv_uniforme (Fourier distante):')
    for s in (0.0, 0.25j, 1.0):
        for M in (60, 200):
            verdade = max(np.linalg.norm(pesado(m, s, 800)[0], 2) for m in (M, M + 2, -M, M + 40))
            cota = ft.Rinv_uniforme(M, MU, s)
            pior = max(pior, verdade/cota)
            print(f'  s {s!s:>6} |m| >= {M}: {verdade:.4f} <= {cota:.4f}')
    print('descida:')
    for (J, Nc) in ((181, 400), (201, 400), (121, 300)):
        verdade = 0.0
        for s in (0.0, 1.0):
            for m in (0, 20, -60):
                R, md = pesado(m, s, Nbig)
                verdade = max(verdade, np.linalg.norm(R[np.ix_(md.n < J, md.n >= Nc)], 2))
        cota = tp.descida(J, Nc)
        pior = max(pior, verdade/cota)
        print(f'  J {J} Nc {Nc}: {verdade:.3e} <= {cota:.3e}')
    print('||B|| numa caixa (fundo completo) contra o simbolo:')
    K, _ = operador(24, 60, s=0.0)
    ms = modos(24); md = tp.Modo(60); k = len(md.n)
    Kf = np.zeros_like(K)
    from signed_operator import Jrad
    for i, m in enumerate(ms):
        Kf[i*k:(i+1)*k, i*k:(i+1)*k] = Jrad(m, 0.0, 60)
    peso = np.concatenate([K1**m*md.w for m in ms])
    B = (peso[:, None]*(K - Kf))/peso[None, :]
    verdade = np.linalg.norm(B, 2)
    pior = max(pior, verdade/tp.NORMA_B)
    print(f'  caixa 24x60: ||B|| = {verdade:.4f} <= {tp.NORMA_B}')
    print(f'pior razao verdade/cota = {pior:.4f}  ->  ' + ('OK' if pior <= 1 else 'COTA VIOLADA'))


if __name__ == '__main__':
    main()
