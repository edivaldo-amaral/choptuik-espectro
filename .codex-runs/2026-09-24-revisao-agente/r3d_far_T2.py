#!/usr/bin/env python3
"""R3d -- bloco (far, T2) EXATO na caixa pequena: P_far B Q0(s2) P_T2, com o fundo COMPLETO.
Q0 preserva F, e B completo leva F no maximo a |m| + 40, n + 100: as linhas far alcancaveis
estao todas na caixa (Mc + 42, Nc + 102). Gram acumulado por modo de saida; verdade =
sqrt(lam_max). Compara com N0[3][2] = alfa_f2 de tile_perron (modo rigoroso)."""
from __future__ import annotations

from pathlib import Path
import sys
import time

import numpy as np

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI.parents[1]/'scripts')); sys.path.insert(0, str(AQUI))
import fourier_tail_bound as ft
import tile_perron as tp
import signed_operator as so
from r3_ataque import Caixa, ZB, GP, FB, DD
from r1_forma_fechada import forma_fechada, lw

MU = ft.MU


def main():
    t0 = time.time()
    tp.RIG = True
    ctx = tp.Ctx(ZB, GP, FB, DD)
    g0 = ft.modos_campos()
    (MP, NP), (Mc, Nc) = GP, FB
    grande = Caixa(Mc + 42, Nc + 102, g0)
    peq = Caixa(Mc, Nc, g0)
    for s2 in (0.1 + 0.1j, 0.0 + 0.25j, 1.125 + 0.125j):
        q2, _, _ = tp.parte_T2(ctx, s2, lambda *a, **k: None)
        # colunas T2 por modo de entrada (na base local da caixa pequena)
        colsT2 = {}; Qs = {}
        for m in peq.ms:
            n_loc = peq.nloc
            sel = np.where(n_loc >= NP)[0] if abs(m) < MP else np.arange(peq.k)
            colsT2[m] = sel
            sg = s2 + 1j*m/2
            R1 = forma_fechada(peq.n1, sg); R2 = forma_fechada(peq.n2, sg)
            Q = np.block([[R1, np.zeros((len(peq.n1), len(peq.n2)))], [np.zeros((len(peq.n2), len(peq.n1))), R2]])
            Qs[m] = (Q*np.exp(lw(n_loc)[:, None] - lw(n_loc)[None, :]))[:, sel]   # pesado (radial)
        ncols = sum(len(v) for v in colsT2.values())
        off = {}; o = 0
        for m in peq.ms:
            off[m] = o; o += len(colsT2[m])
        G = np.zeros((ncols, ncols), complex)
        cacheB = {}
        for mo in grande.ms:
            n_out = grande.nloc
            linsel = np.where(n_out >= Nc)[0] if abs(mo) < Mc else np.arange(grande.k)
            X = np.zeros((len(linsel), ncols), complex)
            for mi in peq.ms:
                d = mo - mi
                if abs(d) > 40:
                    continue
                if d not in cacheB:
                    Bd = ft.B_rad(g0, d, peq.n1, peq.n2, grande.n1, grande.n2, MU)
                    # peso: kappa1^d sqrt(omega(n_out))/sqrt(omega(n_in))
                    cacheB[d] = (ft.K1**d)*Bd*np.exp(lw(grande.nloc)[:, None] - lw(peq.nloc)[None, :])
                X[:, off[mi]:off[mi] + len(colsT2[mi])] = cacheB[d][linsel] @ Qs[mi]
            G += X.conj().T @ X
        verdade = float(np.sqrt(max(np.linalg.eigvalsh(G)[-1], 0.0)))
        cota = q2['alfa_f2']*tp.montar.__defaults__[0] if False else q2['alfa_f2']
        print(f's2 = {s2}: ||P_far B Q0 P_T2|| = {verdade:.6e}   cota alfa_f2 = {cota:.6e}   razao {verdade/cota:.4f}'
              f'   [{time.time()-t0:.0f}s]', flush=True)
        del G


if __name__ == '__main__':
    main()
