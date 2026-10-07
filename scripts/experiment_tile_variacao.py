#!/usr/bin/env python3
"""Diagnóstico (ponto flutuante) para o raio de tile do Rouché de L (C1_REAVALIACAO.md): como o
acoplamento a(s) = ||V_Z(s) Y(s)|| varia com s, e quais quantidades no centro c controlam a variação.
Com R_c = L~_Z(c)^-1 (pesado) e V_Z = J_Z R_c, vale exatamente
    V_Z(s) Y = V_Z(c) Y - eps (V_Z(c) - I) R_c (I + eps R_c)^-1 Y,   eps = s - c,
e (V_Z(c) - I) R_c = -(B + E)_Z R_c^2. Mede: a(c), ||R_c Y||, ||R_c^2 Y||, ||(V_Z(c) - I) R_c Y||, ||R_c||,
e o a(c + eps) verdadeiro (Y com Q0(c + eps)). Gram Y Y^* (como nk_L.acoplamentos_Z)."""
from __future__ import annotations

import argparse
import math
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import nk_L as nk
import signed_operator_L as sl


def gram_Y(ctx, z, MZ, NZ, Mc, Nc):
    nZ = z['n']; ms, offs = z['ms'], z['offs']; idx = {m: i for i, m in enumerate(ms)}
    Ga = np.zeros((nZ, nZ), complex)
    for m in range(-(min(MZ + nk.BAND_M, Mc) - 1), min(MZ + nk.BAND_M, Mc)):
        c = nk.colsT(ctx, m, MZ, NZ, Nc)
        RT = ctx.Q(m, Nc)[:, c]
        Y = np.zeros((nZ, len(c)), complex)
        for mo in ms:
            if abs(mo - m) <= nk.BAND_M:
                i = idx[mo]
                Y[offs[i]:offs[i+1]] = ctx.B(mo, m, NZ, Nc) @ RT
        if m in idx:
            i = idx[m]
            for dk, vw, uw in z['defl']:
                Y += dk*np.outer(vw, uw[offs[i]:offs[i+1]].conj() @ RT[nk.nZ_local(ctx, m, NZ, Nc)])
        Ga += Y @ Y.conj().T
    return Ga


def nrm(M, Ga):
    return math.sqrt(max(np.linalg.eigvalsh(M @ Ga @ M.conj().T)[-1], 0.0))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--Z', default='12x48'); ap.add_argument('--F', default='40x160')
    ap.add_argument('--c', type=complex, default=0j)
    ap.add_argument('--eps', default='1e-4,1e-3,1e-2')
    ap.add_argument('--deflacao', default='0,0.16831')
    a = ap.parse_args()
    MZ, NZ = map(int, a.Z.split('x')); Mc, Nc = map(int, a.F.split('x'))
    dfl = [complex(x) for x in a.deflacao.split(',') if x.strip()]
    def nivel(s):
        ctx = nk.Contexto(s, Nc)
        z = nk.nivel_Z(ctx, MZ, NZ, log=lambda *x, **k: None, sem_borda=True, deflacao=dfl)
        return ctx, z
    ctx, z = nivel(a.c)
    V = z['Vb']; n = z['n']
    QZ = np.zeros((n, n), complex)
    for j, r in enumerate(z['rs']):
        QZ[z['offs'][j]:z['offs'][j+1], z['offs'][j]:z['offs'][j+1]] = ctx.Q(r.m, NZ)
    R = QZ @ V                                                   # R_c = J_Z(c)^-1 V_Z(c) = L~_Z(c)^-1
    Ga = gram_Y(ctx, z, MZ, NZ, Mc, Nc)
    I = np.eye(n)
    print(f'c = {a.c}: ||V_Z|| {np.linalg.norm(V, 2):.2f}, ||R_c|| {np.linalg.norm(R, 2):.2f}, ||Q_ZZ|| {np.linalg.norm(QZ, 2):.3f}')
    print(f'  a(c) = ||V Y|| {nrm(V, Ga):.4f};  ||Y|| {math.sqrt(np.linalg.eigvalsh(Ga)[-1]):.3f};  ||R Y|| {nrm(R, Ga):.3f};'
          f'  ||R^2 Y|| {nrm(R @ R, Ga):.3f};  ||(V - I) R Y|| {nrm((V - I) @ R, Ga):.3f}')
    for e in map(float, a.eps.split(',')):
        for d in (e, 1j*e):
            ctx2, z2 = nivel(a.c + d)
            a2 = nrm(z2['Vb'], gram_Y(ctx2, z2, MZ, NZ, Mc, Nc))
            print(f'  a(c + {d}) = {a2:.4f}')


if __name__ == '__main__':
    main()
