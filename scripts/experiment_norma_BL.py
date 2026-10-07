#!/usr/bin/env python3
"""DIAGNOSTICO (ponto flutuante): cota INFERIOR de ||B_L|| (fundo de L, norma do toro (65/64, 5/4),
m com sinal) pela compressao P B P numa caixa |m| < M, n < N (||P B P|| <= ||B||), por potencia
implicita. Confere a consistencia da cota do simbolo em Arb, NORMA_BL = 3,964043."""
import argparse, sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import nk_L as nk

ap = argparse.ArgumentParser(); ap.add_argument('--M', type=int, default=60); ap.add_argument('--N', type=int, default=200)
ap.add_argument('--iters', type=int, default=80)
a = ap.parse_args()
ctx = nk.Contexto(0.4010247330, a.N)
ms = list(range(-(a.M - 1), a.M)); dims = [ctx.rad(m, a.N).dim for m in ms]
off = np.cumsum([0] + dims); n = int(off[-1]); t0 = time.time()
def mv(x, adj=False):
    y = np.zeros(n, complex)
    for p, mo in enumerate(ms):
        for q, mi in enumerate(ms):
            if abs(mo - mi) <= nk.BAND_M:
                if not adj:
                    y[off[p]:off[p+1]] += ctx.B(mo, mi, a.N, a.N) @ x[off[q]:off[q+1]]
                else:
                    y[off[q]:off[q+1]] += ctx.B(mo, mi, a.N, a.N).conj().T @ x[off[p]:off[p+1]]
    return y
v = np.random.default_rng(0).standard_normal(n) + 0j; v /= np.linalg.norm(v)
for k in range(a.iters):
    u = mv(mv(v), True); v = u/np.linalg.norm(u)
    if k % 20 == 19:
        print(f'  iter {k+1}: ||P B P|| >= {np.linalg.norm(mv(v)):.5f}  [{time.time()-t0:.0f}s]', flush=True)
print(f'caixa |m| < {a.M}, n < {a.N} (dim {n}): ||P B_L P|| ~ {np.linalg.norm(mv(v)):.5f}  (cota do simbolo 3,964043)')
