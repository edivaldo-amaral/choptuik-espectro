#!/usr/bin/env python3
"""Diagnóstico (ponto flutuante, NÃO é prova) para a reavaliação de C1/T1 pela rota de Rouché de operadores:
ao longo de um contorno candidato na faixa A, mede na caixa Z = M x N, nos pesos de L:

  res(s) = ||L_Z(s)^-1||           (resolvente da truncagem; o raio de tile é ~ margem/res)
  tV(s)  = ||J_Z(s) L_Z(s)^-1||    (= ||Hh_ZZ(s)^-1|| com Q0(s): o nível Z do Perron)

L_Z(s) = L_Z(0) + s I; normas por iteração de potência sobre LU (5 a 8 passos; subestimam pouco).
Pontos: metade superior (o espectro é simétrico por conjugação)."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import time

import numpy as np
from scipy.linalg import lu_factor, lu_solve

sys.path.insert(0, str(Path(__file__).resolve().parent))
import nk_L as nk
import signed_operator_L as sl


def norma_inv(lu, n, esq=None, esq_h=None, it=8, rng=np.random.default_rng(0)):
    """||E A^-1||_2 por potência em (E A^-1)^*(E A^-1); E aplicado por esq (função) e E^* por esq_h."""
    x = rng.standard_normal(n) + 1j*rng.standard_normal(n); x /= np.linalg.norm(x)
    v = 0.0
    for _ in range(it):
        y = lu_solve(lu, x)
        if esq is not None:
            y = esq(y)
        v = np.linalg.norm(y)
        z = esq_h(y) if esq is not None else y
        x = lu_solve(lu, z, trans=2); x /= np.linalg.norm(x)
    return float(v)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--Z', default='16x64')
    ap.add_argument('--pontos', required=True, help='lista "re+imj,..." ou arquivo .json')
    ap.add_argument('--saida', default=None)
    a = ap.parse_args()
    M, N = map(int, a.Z.split('x'))
    t0 = time.time()
    L0, rs, offs = sl.operador(M, N, 0.0)
    w = np.concatenate([nk.pesos(r) for r in rs])
    Lw0 = (w[:, None]*L0)/w[None, :]; del L0
    Jb = [(offs[j], offs[j+1], (w[offs[j]:offs[j+1], None]*sl.J_livre(r, 0.0))/w[None, offs[j]:offs[j+1]])
          for j, r in enumerate(rs)]                   # J em blocos (diagonal em m)
    n = len(w)
    pts = json.loads(Path(a.pontos).read_text()) if a.pontos.endswith('.json') else [complex(p) for p in a.pontos.split(',')]
    pts = [complex(p) if not isinstance(p, list) else complex(*p) for p in pts]
    print(f'Z = {a.Z}, dim {n}, montagem {time.time() - t0:.0f}s', flush=True)
    out = []
    for s in pts:
        A = Lw0.copy(); A[np.diag_indices(n)] += s
        lu = lu_factor(A, overwrite_a=True, check_finite=False); del A
        res = norma_inv(lu, n)
        def J(y, s=s):
            z = s*y
            for o0, o1, B in Jb:
                z[o0:o1] += B @ y[o0:o1]
            return z
        def Jh(y, s=s):
            z = np.conj(s)*y
            for o0, o1, B in Jb:
                z[o0:o1] += B.conj().T @ y[o0:o1]
            return z
        tV = norma_inv(lu, n, esq=J, esq_h=Jh)
        out.append(dict(s=[s.real, s.imag], res=res, tV=tV))
        print(f'  s = {s.real:.4f}{s.imag:+.4f}i: ||L_Z^-1|| = {res:9.2f}   ||Hh^-1|| = {tV:9.2f}   [{time.time() - t0:.0f}s]', flush=True)
    if a.saida:
        Path(a.saida).write_text(json.dumps(dict(Z=a.Z, pontos=out), indent=1) + '\n')


if __name__ == '__main__':
    main()
