#!/usr/bin/env python3
"""O operador de constraints K (setor A) na base de Fourier COM SINAL: modos
m par em (-M, M), coeficientes complexos, componente c1 (n impar) e c2 (todo n).
E a base da norma do toro de raios (129/128, 9/8),

    ||x||^2 = sum_{m, n} kappa1^(2m) omega2(n) |x_mn|^2,   m com sinal,

que e a norma em que o simbolo certificou ||B#|| <= 1,65. O montador real
(sharp_operator_builder) usa partes (cos, sen) com peso simetrico omega1(|m|);
complexificado, isso e OUTRA norma (equivalente), e as cotas de blocos de um
certificado precisam de uma so.

Bloco (m + d, m) de K: delta_d0 (J_rad + s + i m/2) + B_d, com B_d da
fourier_tail_bound (inclui mu R_x em d = 0). Validado contra o montador real pelo
espectro (--validar): a mudanca (x_cos, x_sen) -> (x_cos + i x_sen, x_cos - i x_sen)
e uma semelhanca.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fourier_tail_bound as ft

MU = ft.MU
K1, K2 = ft.K1, ft.K2


def radial(N, n0=0):
    """Indices radiais locais de um modo: (c1 impar em [n0, N), c2 em [n0, N))."""
    n1 = np.arange(n0 + (1 - n0 % 2), N, 2)
    n2 = np.arange(n0, N)
    return n1, n2


def peso_radial(n1, n2):
    return np.sqrt(np.concatenate([ft.omega(n1), ft.omega(n2)]))


def Rinv(m, s, N, mu=MU):
    """(J_rad + s + i m/2)^-1 truncado em n < N (exato: J e triangular superior)."""
    sig = s + 1j*m/2
    n1, n2 = radial(N)
    R1 = ft.livre_inv(n1, 2, sig, mu); R2 = ft.livre_inv(n2, 1, sig, mu)
    Z = np.zeros((len(n1), len(n2)))
    return np.block([[R1, Z], [Z.T, R2]])


def Jrad(m, s, N, mu=MU):
    sig = s + 1j*m/2
    blocos = []
    for ns in radial(N):
        J = np.diag(mu*(ns + 1.0)).astype(complex) + sig*np.eye(len(ns))
        for b in range(len(ns)):
            J[:b, b] += 2*mu*ns[b]
        blocos.append(J)
    Z = np.zeros((blocos[0].shape[0], blocos[1].shape[1]))
    return np.block([[blocos[0], Z], [Z.T, blocos[1]]])


def modos(M):
    return list(range(-(M - 2), M - 1, 2))


def operador(M, N, s=0.0, mu=MU, g=None):
    """K(s) na caixa |m| < M, n < N (sem pesos). Devolve (K, lista de (m, k0, k1))."""
    g = ft.modos_campos() if g is None else g
    ms = modos(M)
    n1, n2 = radial(N)
    k = len(n1) + len(n2)
    K = np.zeros((len(ms)*k, len(ms)*k), complex)
    for j, mi in enumerate(ms):
        K[j*k:(j+1)*k, j*k:(j+1)*k] += Jrad(mi, s, N, mu)
        for i, mo in enumerate(ms):
            d = mo - mi
            if abs(d) <= 40:
                K[i*k:(i+1)*k, j*k:(j+1)*k] += ft.B_rad(g, d, n1, n2, n1, n2, mu)
    return K, [(m, i*k, (i+1)*k) for i, m in enumerate(ms)]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--validar', default='8x24')
    a = ap.parse_args()
    M, N = map(int, a.validar.split('x'))
    import sharp_operator_builder as sb
    Kr, _ = sb.operador(M, N)
    Ks, _ = operador(M, N)
    from scipy.optimize import linear_sum_assignment
    er = np.linalg.eigvals(Kr); es = np.linalg.eigvals(Ks)
    C = np.abs(er[:, None] - es[None, :])
    i, j = linear_sum_assignment(C)
    print(f'caixa {M}x{N}: dim real {Kr.shape[0]}, dim assinada {Ks.shape[0]}')
    print(f'  espectros emparelhados: max {C[i, j].max():.2e}, mediana {np.median(C[i, j]):.2e}'
          f'  (max |autovalor| = {np.abs(er).max():.2f})')
    print(f'  menor |autovalor|: {es[np.argmin(np.abs(es))]:.8f}')


if __name__ == '__main__':
    main()
