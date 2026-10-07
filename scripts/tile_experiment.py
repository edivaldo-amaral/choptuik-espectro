#!/usr/bin/env python3
"""EXPERIMENTO EM PONTO FLUTUANTE -- NAO E PROVA. Mede, para um tile em s0, as
quantidades que decidem se um certificado com bloco finito verificado Z e cauda
T fecha, na norma do toro (129/128, 9/8).

H(s0) = K(s0) Q0 = I + B Q0,  Q0 = (K_livre + s0)^-1,  B = mu R_x + Gamma2sharp.
Aproximacao da cauda: uma faixa finita T' alem de Z (caixa grande G').

Inverso aproximado A = diag(V, I), V = H_ZZ^-1. Entao I - A H tem blocos
    [[I - V H_ZZ,  -V H_ZT],
     [ -H_TZ,       I - H_TT]]
e o certificado fecha se a norma espectral da matriz de normas de blocos for < 1.
Mede tambem ||V|| e ||V P_faixa||: o inverso so e grande nos modos baixos?
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path
import sys
import time

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from sharp_operator_builder import operador, MU

K1, K2 = 129/128, 9/8


def norma2(X, iters=60):
    """Norma espectral por iteracao de potencia (estimativa)."""
    if X.size == 0:
        return 0.0
    rng = np.random.default_rng(0)
    v = rng.standard_normal(X.shape[1]) + 0j; v /= np.linalg.norm(v)
    s = 0.0
    for _ in range(iters):
        w = X @ v; u = X.conj().T @ w
        s = math.sqrt(np.linalg.norm(u)); nu = np.linalg.norm(u)
        if nu == 0:
            return 0.0
        v = u/nu
    return float(np.linalg.norm(X @ v))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--s0', type=complex, default=1.0)
    ap.add_argument('--Z', default='12x36')
    ap.add_argument('--G', default='20x60', help='caixa grande (Z + faixa que aproxima a cauda)')
    ap.add_argument('--faixa', type=int, default=41, help='largura radial da faixa de Z junto a T')
    a = ap.parse_args()
    MZ, NZ = map(int, a.Z.split('x')); MG, NG = map(int, a.G.split('x'))
    t = time.time()
    K, addr = operador(MG, NG)
    Kl, _ = operador(MG, NG, mu=MU)          # reusa enderecos; livre montado abaixo
    n = len(addr)
    print(f'caixa grande {MG}x{NG}: {n} DOFs, montada em {time.time()-t:.0f}s', flush=True)
    # parte livre separada (mesma montagem, sem fundo): K_livre = K - B, com B = mu R_x + Gamma
    from sharp_operator_builder import enderecos
    idx = {x: i for i, x in enumerate(addr)}
    Kf = np.zeros((n, n))
    for col, (d, part, m, nn) in enumerate(addr):
        Kf[col, col] += MU*(nn + 1)
        if m:
            Kf[idx[(d, 1-part, m, nn)], col] += (m/2)*(1 if part == 0 else -1)
        passo = 2 if d == 0 else 1
        for low in range(nn - passo, -1, -passo):
            Kf[idx[(d, part, m, low)], col] += 2*MU*nn
    del Kl
    s0 = a.s0
    om = lambda k, j: 1.0 if j == 0 else k**(2*j) + k**(-2*j)
    dvec = np.sqrt(np.array([om(K1, m)*om(K2, nn) for _, _, m, nn in addr]))
    B = K - Kf
    Kfs = Kf.astype(complex) + s0*np.eye(n)
    Q0 = np.linalg.inv(Kfs)
    H = np.eye(n) + B @ Q0
    Hw = (dvec[:, None]*H)/dvec[None, :]
    del H, Q0, K, Kf, Kfs
    Zi = np.array([i for i, (_, _, m, nn) in enumerate(addr) if m < MZ and nn < NZ])
    Ti = np.array([i for i, (_, _, m, nn) in enumerate(addr) if not (m < MZ and nn < NZ)])
    fx = np.array([k for k, i in enumerate(Zi) if addr[i][3] >= NZ - a.faixa or addr[i][2] >= MZ - 21])
    HZZ = Hw[np.ix_(Zi, Zi)]
    V = np.linalg.inv(HZZ)
    rho = norma2(np.eye(len(Zi)) - V @ HZZ)
    VHzt = V @ Hw[np.ix_(Zi, Ti)]
    q = dict(dimZ=len(Zi), dimT=len(Ti),
             normV=norma2(V), normV_faixa=norma2(V[:, fx]), normV_linhas_faixa=norma2(V[fx, :]),
             rho=rho, VHzt=norma2(VHzt), Htz=norma2(Hw[np.ix_(Ti, Zi)]),
             Hzt=norma2(Hw[np.ix_(Zi, Ti)]), alfa=norma2(Hw[np.ix_(Ti, Ti)] - np.eye(len(Ti))))
    Mj = np.array([[q['rho'], q['VHzt']], [q['Htz'], q['alfa']]])
    q['majorante'] = float(np.linalg.norm(Mj, 2))
    print(f"s0={s0}  Z={a.Z} ({q['dimZ']})  T'={a.G} menos Z ({q['dimT']})")
    for k in ('normV', 'normV_faixa', 'normV_linhas_faixa', 'rho', 'Hzt', 'VHzt', 'Htz', 'alfa', 'majorante'):
        print(f'   {k:>20} = {q[k]:.4f}')
    print(f'   tempo total {time.time()-t:.0f}s')


if __name__ == '__main__':
    main()
