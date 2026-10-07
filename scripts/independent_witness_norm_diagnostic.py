#!/usr/bin/env python3
"""DIAGNOSTICO EM PONTO FLUTUANTE (revisao cruzada, rodada 2) -- NAO E PROVA.

Pergunta (i) do orquestrador: na raiz ~0.401 de L_A, a razao ||C(s)h||/||h|| cai com
o truncamento na norma ponderada (129/128,9/8). Isso e C h -> 0 ou ||h|| ponderado
crescendo? E em que norma o testemunho deve valer?

Normas consideradas (todas l1 ponderadas, soma por coluna/componente, |Re|+|Im| ja
separado nos DOFs reais):
  weak(v)    pesos (129/128,9/8)            -- saida de C e dominio sharp (SHARP_DOMAIN)
  strongE(v) pesos (65/64,5/4) x eta        -- Y+ da realizacao H=L Q (A1/A2/A3)
  x = J_ref h                               -- coordenada em que tiles/caudas enclausuram
Parsers proprios; nao importa scripts do Fisico nem do ChatGPT.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SPEC = ROOT/'build/spectrum'
MU = 722873400/2**32
ETA = np.array([916., 4096., 243., 1233.])


def read_selected(path):
    with open(path) as fh:
        assert fh.readline().strip() == 'CHOPTUIK_RT_SPECTRAL_MATRIX_V1'
        dim = int(fh.readline().split()[1])
        line = fh.readline()
        while line.strip() != 'BEGIN_ADDRESSES':
            line = fh.readline()
        addr = [tuple(map(int, fh.readline().split()))[1:] for _ in range(dim)]
        assert fh.readline().strip() == 'BEGIN_ENTRIES'
        A = np.zeros((dim, dim))
        for line in fh:
            f = line.split()
            if f:
                A[int(f[0]), int(f[1])] += math.ldexp(float(int(f[2])), int(f[3]) if len(f) > 3 else 0)
    return A, addr


def read_constraints(path, columns):
    with open(path) as fh:
        assert fh.readline().strip() == 'CHOPTUIK_RT_SPECTRAL_CONSTRAINTS_V1'
        rows = int(fh.readline().split()[1])
        assert int(fh.readline().split()[1]) == columns
        assert fh.readline().strip() == 'BEGIN_ADDRESSES'
        addr = [tuple(map(int, fh.readline().split()))[1:] for _ in range(rows)]
        assert fh.readline().strip() == 'BEGIN_ENTRIES'
        C0, C1 = np.zeros((rows, columns)), np.zeros((rows, columns))
        for line in fh:
            f = line.split()
            if f:
                (C0 if int(f[0]) == 0 else C1)[int(f[1]), int(f[2])] = math.ldexp(
                    float(int(f[3])), int(f[4]) if len(f) > 4 else 0)
    return C0, C1, addr


def weights(addr, k1, k2, eta=None):
    w = np.array([(2-(m == 0))*(2-(n == 0))*k1**m*k2**n for d, p, m, n in addr])
    if eta is not None:
        w = w*np.array([eta[d] for d, p, m, n in addr])
    return w


def free(addr, mu=MU):
    look = {a: i for i, a in enumerate(addr)}
    J = np.zeros((len(addr), len(addr)))
    for c, (d, p, m, n) in enumerate(addr):
        J[c, c] += mu*(n+1)
        if m:
            J[look[(d, 1-p, m, n)], c] += (m/2 if p == 0 else -m/2)
        step = 2 if d == 0 else 1
        for j in range(n-step, -1, -step):
            J[look[(d, p, m, j)], c] += 2*mu*n
    return J


def analyse(M, N, target):
    A, addr = read_selected(SPEC/f'rt-A-{M}x{N}.dat')
    C0, C1, caddr = read_constraints(SPEC/f'rt-A-{M}x{N}.dat.constraints', len(addr))
    lam, vec = np.linalg.eig(A)
    s_all = -lam
    k = int(np.argmin(np.abs(s_all-target)))
    s = s_all[k]
    h = vec[:, k]
    h = h/h[np.argmax(np.abs(h))]
    Ch = (C0+s*C1)@h
    J = free(addr)
    x = J@h
    wweak = weights(addr, 129/128, 9/8)
    wstrE = weights(addr, 65/64, 5/4, ETA)
    wstr = weights(addr, 65/64, 5/4)
    cweak = weights(caddr, 129/128, 9/8)
    nrm = lambda v, w: float(np.sum(w*np.abs(v)))
    small = np.array([m < 4 and n < 12 for d, p, m, n in addr])
    outer = np.array([m >= M-2 or n >= N-6 for d, p, m, n in addr])
    csmall = np.array([m < 4 and n < 12 for d, p, m, n in caddr])
    return dict(
        box=[M, N], s=complex(s).real, s_imag=complex(s).imag,
        ratio_l2=float(np.linalg.norm(Ch)/np.linalg.norm(h)),
        ratio_weak_weak=nrm(Ch, cweak)/nrm(h, wweak),
        ratio_weak_over_strongEta_h=nrm(Ch, cweak)/nrm(h, wstrE),
        ratio_weak_over_strongEta_x=nrm(Ch, cweak)/nrm(x, wstrE),
        ratio_weak_over_graph_strong=nrm(Ch, cweak)/(nrm(h, wstr)+nrm(x, wstr)),
        h_weak_norm_over_l2=nrm(h, wweak)/np.linalg.norm(h),
        h_weak_fraction_outer_band=nrm(h*outer, wweak)/nrm(h, wweak),
        x_strongEta_fraction_outer_band=nrm(x*outer, wstrE)/nrm(x, wstrE),
        Ch_weak_over_l2=nrm(Ch, cweak)/np.linalg.norm(Ch),
        box4x12_ratio_weak=nrm(Ch*csmall, cweak[:]*1)/max(nrm(h*small, wweak), 1e-300),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--target', type=complex, default=0.401+0j)
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    rows = [analyse(M, N, args.target) for M, N in ((6, 18), (8, 24), (10, 30), (12, 36))]
    text = json.dumps(dict(note='DIAGNOSTICO float; nao e prova', target=str(args.target), rows=rows), indent=2)
    if args.out:
        args.out.write_text(text+'\n')
    print(text)


if __name__ == '__main__':
    main()
