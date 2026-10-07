#!/usr/bin/env python3
"""R4 (revisor): ||V_J|| = ||(I + P_J B_16 Q0(c) P_J)^-1|| para a janela 192..199 (e, se der tempo, 80..95) da
cauda de 0,5i, por SVD densa. Codigo proprio para os pesos (kappa1^m sqrt(omega(n))), para Q0(c) (J triangular
superior: diag mu(n+1) + c + i m/2, 2 mu n acima) e para a montagem/inversao; os blocos NAO pesados de B vem de
signed_operator_L.bloco_B (validado pelo autor contra o exportador de RT; nao reimplementado aqui).
Compara com nV do T3.json (cota superior rigorosa)."""
import json, sys, time
from pathlib import Path
import numpy as np
from scipy.linalg import solve_triangular
sys.path.insert(0, '<repo>/scripts')
import signed_operator_L as sl

MU = 722873400/2**32; K1 = 65/64; K2 = 5/4; Nc = 400; BANDA = 16
c = 0.5j
t0 = time.time()
def omega(n):
    return np.where(n == 0, 1.0, K2**(2.0*n) + K2**(-2.0*n))
def ns_modo(m):
    return [np.arange(1, Nc, 2), np.arange(Nc), np.arange(Nc)] if m % 2 == 0 else [np.arange(Nc)]
def w_modo(m):
    return K1**m*np.sqrt(omega(np.concatenate(ns_modo(m))))
def Q0w(m):
    bl = []
    for ns in ns_modo(m):
        J = np.diag(MU*(ns + 1.0) + c + 1j*m/2).astype(complex) + np.triu(np.tile(2*MU*ns, (len(ns), 1)), 1)
        Q = solve_triangular(J, np.eye(len(ns)))
        w = np.sqrt(omega(ns))
        bl.append((w[:, None]*Q)/w[None, :])
    n = sum(len(b) for b in bl); out = np.zeros((n, n), complex); o = 0
    for b in bl:
        out[o:o+len(b), o:o+len(b)] = b; o += len(b)
    return out
g = sl.campos()
res = {}
for nome, modos in (('192..199', list(range(192, 200))), ('-199..-192', list(range(-199, -191)))):
    dims = [sum(len(x) for x in ns_modo(m)) for m in modos]; off = np.cumsum([0] + dims)
    X = np.zeros((off[-1], off[-1]), complex)
    Qs = {m: Q0w(m) for m in modos}
    for j, mi in enumerate(modos):
        ri = sl.Radial(mi, Nc)
        for i, mo in enumerate(modos):
            if abs(mo - mi) > BANDA:
                continue
            ro = sl.Radial(mo, Nc)
            Bb = sl.bloco_B(g, ro, ri)
            Bw = (w_modo(mo)[:, None]*Bb)/w_modo(mi)[None, :]
            X[off[i]:off[i+1], off[j]:off[j+1]] = Bw @ Qs[mi]
    H = np.eye(off[-1]) + X
    s = np.linalg.svd(H, compute_uv=False)
    nV = 1/s.min()
    res[nome] = dict(dim=int(off[-1]), nV_svd=float(nV), sigma_min=float(s.min()), sigma_max=float(s.max()))
    print(f'janela {nome} (dim {off[-1]}): ||V_J|| (SVD) = {nV:.10f}  [{time.time()-t0:.0f}s]', flush=True)
    del X, H
t3 = json.load(open('<repo>/build/faixaB/cauda/T3_0_5.json'))
for nome in res:
    i = t3['nome'].index(nome)
    res[nome]['nV_T3'] = t3['nV'][i]
    print(f'{nome}: T3 nV = {t3["nV"][i]:.10f}; razao T3/SVD = {t3["nV"][i]/res[nome]["nV_svd"]:.8f}')
json.dump(res, open('<repo>/.codex-runs/2026-10-02-revisao-B/r4_nV_resultado.json', 'w'), indent=1)
