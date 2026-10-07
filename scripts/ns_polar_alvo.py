#!/usr/bin/env python3
"""Polar l >= 2: autovalores perto de um alvo sigma por shift-invert (LU densa de A - sigma B + Arnoldi), com o residuo
relativo dos vinculos (cpp, cmm) de cada autovetor. Exploracao numerica (NAO rigorosa); ver docs/NAO_ESFERICO.md."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys
import time

import numpy as np
import scipy.linalg as sla
from scipy.sparse.linalg import LinearOperator, eigs
from numpy.polynomial import chebyshev as C

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ns_axial as na
import ns_polar as npl


def operadores(l, Nt, Nx, xs, eqs, campos, mu):
    tau, xi, Dt, Dx = na.matrizes(Nt, Nx, xs)
    bg, _ = npl.fundo(campos, mu, tau, xi)
    Xg = np.broadcast_to(xi[None, :], (Nt, Nx)); env = dict(bg); env.update(xi=Xg, mu=mu, l=l)
    COEF = npl.coeficientes(); It, Ix = np.eye(Nt), np.eye(Nx)
    DT = np.kron(Dt, Ix); DX = np.kron(It, Dx); R = np.kron(It, np.eye(Nx)[::-1]); n = Nt*Nx; I = np.eye(n)
    unk = ['Phi', 'k', 'A']; M = [np.zeros((len(eqs)*n, 3*n)) for _ in range(3)]
    for ie, e in enumerate(eqs):
        for (u, i, j), expr in COEF[e].items():
            c = np.broadcast_to(np.asarray(npl.avalia_coef(expr, env), float), (Nt, Nx)).ravel()
            Dj = np.linalg.matrix_power(DX, j) if j else I
            bl = [(0, Dj)] if i == 0 else ([(0, DT @ Dj), (1, Dj)] if i == 1 else [(0, DT @ DT), (1, 2*DT), (2, I)])
            col = unk.index(u) if u != 'B' else 2
            for pw, Mat in bl:
                op = c[:, None]*Mat
                if u == 'B':
                    op = (op @ R)*(-1)**l
                M[pw][ie*n:(ie + 1)*n, col*n:(col + 1)*n] += op
    Peven = np.array([xi**l*C.chebval(xi/xs, np.eye(Nx)[2*k]) for k in range(Nx//2)]).T
    Pfull = np.array([xi**l*C.chebval(xi/xs, np.eye(Nx)[k]) for k in range(Nx)]).T
    P = sla.block_diag(np.kron(It, Peven), np.kron(It, Peven), np.kron(It, Pfull))
    return M, P


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--l', type=int, default=2)
    ap.add_argument('--res', default='25x40')
    ap.add_argument('--xs', type=float, default=1.02)
    ap.add_argument('--alvo', type=complex, default=complex(-0.005, 0.147))
    ap.add_argument('--k', type=int, default=10)
    a = ap.parse_args()
    campos, mu = na.campos_refA()
    for r in a.res.split(','):
        Nt, Nx = map(int, r.split('x'))
        t0 = time.time()
        M, P = operadores(a.l, Nt, Nx, a.xs, ['wave', 'trace', 'EAm'], campos, mu)
        Mc, _ = operadores(a.l, Nt, Nx, a.xs, ['cpp', 'cmm'], campos, mu)
        Pi = np.linalg.pinv(P)
        Mr = [Pi @ Mk @ P for Mk in M]; nr = Mr[0].shape[0]
        Z = np.zeros((nr, nr)); Id = np.eye(nr)
        A = np.block([[Z, Id], [Mr[0], Mr[1]]]); B = np.block([[Id, Z], [Z, -Mr[2]]])
        sig = a.alvo
        lu = sla.lu_factor(A - sig*B)
        op = LinearOperator(A.shape, matvec=lambda v: sla.lu_solve(lu, B @ v), dtype=complex)
        th, V = eigs(op, k=a.k, which='LM', tol=1e-12, maxiter=5000)
        s_ = sig + 1/th
        print(f'{Nt}x{Nx} (n = {A.shape[0]}, {time.time() - t0:.0f}s), alvo {sig}: perto do alvo, residuo de vinculo:', flush=True)
        ordem = np.argsort(np.abs(s_ - sig))
        for i in ordem:
            s = s_[i]; u = P @ V[:nr, i]
            cr = (Mc[0] + s*Mc[1] + s*s*Mc[2]) @ u
            esc = np.linalg.norm(Mc[0] @ u) + abs(s)*np.linalg.norm(Mc[1] @ u) + abs(s)**2*np.linalg.norm(Mc[2] @ u)
            ev_res = (Mr[0] + s*Mr[1] + s*s*Mr[2]) @ V[:nr, i]
            print(f'   s = {s.real:+.7f} {s.imag:+.7f}i   kappa*Delta = {2*np.pi*s.real/na.K_RT*2*na.K_RT:+.4f}   '
                  f'vinculo {np.linalg.norm(cr)/esc:.1e}', flush=True)


if __name__ == '__main__':
    main()
