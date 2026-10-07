#!/usr/bin/env python3
"""Diagnóstico (ponto flutuante): posto numérico do acoplamento Z <- T do Rouché de L. Com G = Y Y^*
(Y = P_Z (B + E) Q0(s) P_T, Gram de nk_L.acoplamentos_Z, com deflação), imprime o decaimento dos
autovalores de G e a(k) = ||V_Z Y_k|| com Y_k a projeção de posto k, contra a(completo)."""
from __future__ import annotations
import argparse, math, sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import nk_L as nk
from experiment_tile_variacao import gram_Y

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--s', type=complex, default=0j)
    ap.add_argument('--Z', default='32x128'); ap.add_argument('--F', default='200x400')
    ap.add_argument('--deflacao', default='0,0.16831')
    a = ap.parse_args(); t0 = time.time()
    MZ, NZ = map(int, a.Z.split('x')); Mc, Nc = map(int, a.F.split('x'))
    ctx = nk.Contexto(a.s, Nc)
    z = nk.nivel_Z(ctx, MZ, NZ, sem_borda=True, deflacao=[complex(x) for x in a.deflacao.split(',') if x.strip()])
    print(f'nivel Z [{time.time()-t0:.0f}s]', flush=True)
    G = gram_Y(ctx, z, MZ, NZ, Mc, Nc); print(f'Gram [{time.time()-t0:.0f}s]', flush=True)
    lam, U = np.linalg.eigh(G); lam = lam[::-1]; U = U[:, ::-1]; del G
    print('autovalores de Y Y^* (raiz): ' + ', '.join(f'{k}:{math.sqrt(max(lam[k], 0)):.2e}' for k in (0, 10, 50, 100, 200, 400, 800, 1600, 3200) if k < len(lam)), flush=True)
    V = z['Vb']
    def a_k(k):
        Wk = V @ (U[:, :k]*np.sqrt(np.maximum(lam[:k], 0)))
        return float(np.linalg.norm(Wk, 2))
    for k in (50, 100, 200, 400, 800, 1600):
        if k < len(lam):
            print(f'  a(posto {k}) = {a_k(k):.5f};  resto ||Y - Y_k|| = {math.sqrt(max(lam[k], 0)):.2e}   [{time.time()-t0:.0f}s]', flush=True)
    print(f'  a(completo) = {a_k(len(lam)):.5f}   [{time.time()-t0:.0f}s]', flush=True)

if __name__ == '__main__':
    main()
