#!/usr/bin/env python3
"""A resolvente tem decaimento exponencial fora da diagonal? (premissa Combes-Thomas)

A cauda em uso, t(G)~81/(N+1), e ALGEBRICA e vem so do crescimento diagonal de
J_mu. Combes-Thomas daria decaimento EXPONENCIAL, explorando que B e bandada.
Aqui se mede o decaimento real de V=H^-1 em funcao da distancia de indice
radial, na norma conjugada pelos pesos.
"""
import sys, time
from fractions import Fraction as Q
from pathlib import Path
import numpy as np
sys.path.insert(0,'scripts'); sys.path.insert(0,'<scratch>')
from kappa import enderecos, pesos, MU
from check_free_preconditioner import (read_spectral_matrix, read_addresses,
                                       free_matrix, exact_inverse, multiply)

MALHA = sys.argv[1] if len(sys.argv) > 1 else '6x18'
CENTRO = Q(sys.argv[2]) if len(sys.argv) > 2 else Q(1)

p = Path(f'build/spectrum/rt-A-{MALHA}.dat')
dim, raw, expo = read_spectral_matrix(p)
addr_raw, addr = read_addresses(p, dim), enderecos(p, dim)
w = np.array([float(x) for x in pesos(addr, Q(65,64), Q(5,4))])

t0 = time.time()
principal = free_matrix(addr_raw, MU)
inv_free = exact_inverse(principal)
pencil = [[Q(x)*Q(2)**expo for x in row] for row in raw]
for i in range(dim): pencil[i][i] += CENTRO
H = np.array([[float(x) for x in row] for row in multiply(pencil, inv_free)])
print(f"H montada ({dim}x{dim}) em {time.time()-t0:.0f}s; centro s={float(CENTRO)}")

V = np.linalg.inv(H)
Vc = V*w[:, None]/w[None, :]          # conjugada pelos pesos
print(f"||V||_1 conjugada = {np.abs(Vc).sum(0).max():.1f}")

n_idx = np.array([a[3] for a in addr])
m_idx = np.array([a[2] for a in addr])
print(f"\n{'dist n':>7} {'pares':>8} {'max |V_ij|':>13} {'razao':>8}")
print('-'*42)
D = np.abs(n_idx[:, None]-n_idx[None, :])
ant = None
for d in range(0, int(D.max())+1):
    mask = D == d
    if not mask.any(): continue
    v = np.abs(Vc[mask]).max()
    r = v/ant if ant else float('nan')
    print(f"{d:>7} {int(mask.sum()):>8} {v:>13.3e} {r:>8.3f}")
    ant = v

print(f"\n{'dist m':>7} {'pares':>8} {'max |V_ij|':>13} {'razao':>8}")
print('-'*42)
Dm = np.abs(m_idx[:, None]-m_idx[None, :])
ant = None
for d in range(0, int(Dm.max())+1):
    mask = Dm == d
    if not mask.any(): continue
    v = np.abs(Vc[mask]).max()
    r = v/ant if ant else float('nan')
    print(f"{d:>7} {int(mask.sum()):>8} {v:>13.3e} {r:>8.3f}")
    ant = v
