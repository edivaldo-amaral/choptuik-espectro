#!/usr/bin/env python3
"""R3c -- ||P_F B P_F|| (compressao do fundo completo, pesado) numa caixa: cota INFERIOR de ||B||,
que tem de ficar abaixo de NORMA_B = 1,64967 (entrada externa, simbolo)."""
from pathlib import Path
import sys
import numpy as np
AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI.parents[1]/'scripts')); sys.path.insert(0, str(AQUI))
import fourier_tail_bound as ft
import tile_perron as tp
from r3_ataque import Caixa, pesar
from scipy.sparse.linalg import svds
g0 = ft.modos_campos()
for (M, N) in ((16, 200), (40, 60), (24, 120)):
    cx = Caixa(M, N, g0)
    B = pesar(cx.Bdenso(), cx.lw, cx.lw)
    v = float(svds(B, k=1, return_singular_vectors=False, tol=1e-10)[0])
    print(f'caixa {M}x{N} (dim {cx.dim}): ||P_F B P_F|| = {v:.5f}  contra NORMA_B = {tp.NORMA_B}')
    del B
