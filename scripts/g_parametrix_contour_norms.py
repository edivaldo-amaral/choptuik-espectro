#!/usr/bin/env python3
"""DIAGNOSTICO em float (nao e prova) para docs/G_PARAMETRIX_FEASIBILITY.md, itens 3-4.

(1) ||Q H(s)^-1|| e ||H(s)^-1|| (norma eta x pesos RT, soma por coluna) ao longo de Gamma_A e
    do retangulo B, na matriz 12x36 -> numero de centros de uma parametriz por centro.
(2) Bloco de casca H_OO(s), O = 12x36 - nucleo (6x18 ou 8x24): ||H_OO^-1||, defeito da
    Neumann amortecida de grau 32 e raio projetado (1-a)/(||H_OO^-1|| ||Q||).
(3) Fatia de beta_eta vinda dos coeficientes de RefA com m!=0 (majorante EXATO por componentes).
(4) Micro-benchmark: uma aplicacao matriz-livre completa de B (18 convolucoes FFT) por grade.
"""
from __future__ import annotations

import json
import time
from fractions import Fraction as Q

import numpy as np
from scipy.signal import fftconvolve

import g_parametrix_prototype as G
import independent_column_diagnostic as D
from independent_component_beta import parse_refa, selected_majorant, gamma1_majorant, eta_norm


def cn(M):
    return float(np.max(np.sum(np.abs(M.real)+np.abs(M.imag), axis=0)))


def main():
    A, addr = G.read(G.ROOT/'build/spectrum/rt-A-12x36.dat')
    n = len(addr)
    Qm = np.linalg.inv(G.free(addr))
    w = np.array([G.ETA[d]*(2-(m == 0))*(2-(nn == 0))*(65/64)**m*(5/4)**nn for d, p, m, nn in addr])
    Qw = w[:, None]*Qm/w[None, :]
    out = dict(note='DIAGNOSTICO float; nao e prova', resolvent={}, shell={})
    pts = [(0.125, 0), (0.125, 0.125), (0.125, 0.25), (0.5, 0.25), (1, 0.25), (1, 0.125), (1, 0),
           (0.125, 0.5), (0.125, 0.75), (0.5, 0.75), (1, 0.75), (1, 0.5), (0.5, 0.5)]
    for x, y in pts:
        s = x+1j*y
        R = np.linalg.inv(A+s*np.eye(n))
        Hinv = np.linalg.inv((A+s*np.eye(n))@Qm)
        out['resolvent'][f'{x}+{y}i'] = dict(norm_QHinv=cn(w[:, None]*R/w[None, :]),
                                             norm_Hinv=cn(w[:, None]*Hinv/w[None, :]))
    for cm, cnn in ((6, 18), (8, 24)):
        Oi = np.where(np.array([not (m < cm and nn < cnn) for d, p, m, nn in addr]))[0]
        for x, y in pts:
            s = x+1j*y
            H = (w[:, None]*((A+s*np.eye(n))@Qm)/w[None, :]).astype(complex)
            HOO = H[np.ix_(Oi, Oi)]
            inv = np.linalg.inv(HOO)
            best = min(cn(np.linalg.matrix_power(np.eye(len(Oi))-om*HOO, 32)) for om in (1.0, 0.5))
            out['shell'][f'core{cm}x{cnn} s={x}+{y}i'] = dict(
                norm_HOO_inv=cn(inv), neumann32_best_defect=best,
                radius_proj=(1-min(best, 0.99))/(cn(inv)*cn(Qw)))
    fields, scale = parse_refa()
    lin = gamma1_majorant(Q(17, 100), Q(5, 4))
    tot = lambda M: eta_norm([[a+b for a, b in zip(x, y)] for x, y in zip(M, lin)])
    Af, _ = selected_majorant(fields, scale)
    An, _ = selected_majorant([{k: v for k, v in f.items() if k[0] != 0} for f in fields], scale)
    out['beta_eta_total'] = float(tot(Af))
    out['beta_eta_only_m_nonzero_plus_gamma1'] = float(tot(An))
    coefs = D.coef_arrays()
    kern = [a for arrs in coefs.values() for a in arrs if np.any(a)]
    bench = {}
    for M, N in ((24, 72), (64, 200), (107, 384)):
        x = np.random.rand(2*M+1, 2*N+1)+1j*np.random.rand(2*M+1, 2*N+1)
        t = time.time()
        for k in kern:
            fftconvolve(k, x, mode='full')
        bench[f'{M}x{N}'] = time.time()-t
    out['benchmark_seconds_per_full_application_1core'] = bench
    out['convolutions_per_application'] = len(kern)
    (G.ROOT/'build/spectrum/g-parametrix-contour-norms.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps({k: v for k, v in out.items() if k not in ('resolvent', 'shell')}, indent=2))


if __name__ == '__main__':
    main()
