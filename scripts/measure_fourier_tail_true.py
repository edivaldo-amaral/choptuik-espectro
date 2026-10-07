#!/usr/bin/env python3
# DIAGNOSTICO EM PONTO FLUTUANTE -- NAO E PROVA. Ver S3_CONDICIONAMENTO_SHARP.md 10.4.
"""Norma VERDADEIRA (ponto flutuante, truncada) da cauda de Fourier ||P_F B Q0 P_F||,
montando os blocos (m+d, m) = kappa1^d-pesados B_d R_m^-1 para m em [M, Mtop]."""
import sys; sys.path.insert(0, 'scripts')
import numpy as np, time
import fourier_tail_bound as ft
g = ft.modos_campos()
NR = 240
n1 = np.arange(1, NR, 2); n2 = np.arange(NR)
w = np.sqrt(np.concatenate([ft.omega(n1), ft.omega(n2)]))
def Rm(m, s0):
    sig = s0 + 1j*m/2
    R1 = ft.livre_inv(n1, 2, sig, ft.MU); R2 = ft.livre_inv(n2, 1, sig, ft.MU)
    return np.block([[R1, np.zeros((len(n1), len(n2)))], [np.zeros((len(n2), len(n1))), R2]])
Bcache = {}
def Bd(d):
    if d not in Bcache:
        Bcache[d] = ft.B_rad(g, d, n1, n2, n1, n2, ft.MU)       # saidas truncadas no mesmo intervalo
    return Bcache[d]
for M, span in ((20, 60), (30, 60)):
    t = time.time()
    ms = list(range(M, M + span + 1, 2))
    k = len(n1) + len(n2)
    X = np.zeros((len(ms)*k, len(ms)*k), complex)
    for j, m in enumerate(ms):
        Rw = (w[:, None]*Rm(m, 0.0))/w[None, :]
        for i, mo in enumerate(ms):
            d = mo - m
            if abs(d) > 40: continue
            Bw = (w[:, None]*Bd(d))/w[None, :]
            X[i*k:(i+1)*k, j*k:(j+1)*k] = ft.K1**d*(Bw @ Rw)
    v = np.linalg.norm(X, 2) if X.shape[0] < 9000 else ft.norma2_sup(X)
    print(f'M={M}: ||P_F B Q0 P_F|| (m em [{M},{M+span}], n<{NR}) = {v:.4f}   ({time.time()-t:.0f}s, dim {X.shape[0]})', flush=True)
print('comparar: soma modo a modo 0,879 (M=20), 0,626 (M=30)')
