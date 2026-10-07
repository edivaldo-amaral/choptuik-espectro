#!/usr/bin/env python3
"""Diagnóstico: custo da decomposição de Schur da truncagem deflacionada L~_Z (pesada) e das
resoluções triangulares com k colunas, base do esquema rápido de tiles do Rouché de L."""
from __future__ import annotations
import argparse, sys, time
from pathlib import Path
import numpy as np
from scipy.linalg import schur, solve_triangular
sys.path.insert(0, str(Path(__file__).resolve().parent))
import nk_L as nk
import signed_operator_L as sl

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--Z', default='32x128'); ap.add_argument('--k', type=int, default=400)
    a = ap.parse_args(); t0 = time.time()
    M, N = map(int, a.Z.split('x'))
    L0, rs, offs = sl.operador(M, N, 0.0)
    w = np.concatenate([nk.pesos(r) for r in rs]); A = (w[:, None]*L0)/w[None, :]; del L0
    n = len(w); print(f'dim {n}, montagem [{time.time()-t0:.0f}s]', flush=True)
    T, U = schur(A, output='complex', overwrite_a=False, check_finite=False)
    print(f'Schur [{time.time()-t0:.0f}s]', flush=True)
    x = np.random.default_rng(1).standard_normal((n, 4)) + 0j          # conferencias por vetores (memoria)
    r1 = np.linalg.norm(A @ x - U @ (T @ (U.conj().T @ x)))/np.linalg.norm(A @ x)
    r2 = np.linalg.norm(U.conj().T @ (U @ x) - x)/np.linalg.norm(x)
    print(f'||(A - U T U^*) x||/||A x|| = {r1:.2e}; ||(U^*U - I) x||/||x|| = {r2:.2e}  [{time.time()-t0:.0f}s]', flush=True)
    ev = -np.diag(T)
    sel = sorted([e for e in ev if -0.3 < e.real < 3.1 and abs(e.imag) <= 0.3], key=lambda z: z.real)
    print('raizes s na faixa: ' + ', '.join(f'{e.real:.6f}{e.imag:+.1e}i' for e in sel), flush=True)
    Bk = np.random.default_rng(0).standard_normal((n, a.k)) + 0j
    for s in (0.0, 0.25j):
        t1 = time.time(); X = solve_triangular(T + s*np.eye(n), Bk, check_finite=False)
        print(f'solve_triangular k={a.k} em s={s}: {time.time()-t1:.1f}s', flush=True)

if __name__ == '__main__':
    main()
