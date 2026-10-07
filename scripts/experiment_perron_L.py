#!/usr/bin/env python3
"""DIAGNOSTICO EM PONTO FLUTUANTE -- NAO E PROVA. Viabilidade do tile de Perron para L (S3b):
numa caixa G', H = I + B Q0 (pesos do toro (65/64, 5/4), m com sinal), bloco Z invertido,
T = G' \\ Z. Mede ||V||, ||V Q_ZZ||, a = ||V H_ZT||, b = ||H_TZ||, alfa = ||H_TT - I|| e o
criterio de Schur rho + a b/(1 - alfa)."""
import argparse, math, sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import signed_operator_L as sl
import fourier_tail_bound as ft
from tile_certificate import norma2


def omega2(n):
    return np.where(n == 0, 1.0, sl.K2**(2.0*n) + sl.K2**(-2.0*n))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--s0', type=complex, default=0.276)
    ap.add_argument('--Z', default='16x64'); ap.add_argument('--G', default='28x112')
    ap.add_argument('--gamma1-livre', action='store_true', help='mu S Gamma1 na parte livre (Q0 = (O_mu + mu S Gamma1 + s)^-1)')
    a = ap.parse_args()
    MZ, NZ = map(int, a.Z.split('x')); MG, NG = map(int, a.G.split('x'))
    t0 = time.time()
    g = sl.campos()
    L0, rs, offs = sl.operador(MG, NG, 0.0, g=g)
    Kf = np.zeros_like(L0)
    for i, r in enumerate(rs):
        Kf[offs[i]:offs[i+1], offs[i]:offs[i+1]] = sl.J_livre(r, 0.0)
        if a.gamma1_livre:     # bloco d = 0 de mu S Gamma1 (fundo com campos nulos)
            nulo = [{} for _ in range(4)]
            Kf[offs[i]:offs[i+1], offs[i]:offs[i+1]] += sl.bloco_B(nulo, r, r)
    B = L0 - Kf
    del L0
    peso = np.concatenate([sl.K1**r.m*np.sqrt(np.concatenate([omega2(ns) for _, ns in r.comps])) for r in rs])
    Q = np.zeros_like(Kf)
    for i, r in enumerate(rs):
        blk = slice(offs[i], offs[i+1])
        Q[blk, blk] = np.linalg.inv(Kf[blk, blk] + a.s0*np.eye(r.dim))
    del Kf
    W = lambda X: (peso[:, None]*X)/peso[None, :]
    H = W(B @ Q); H[np.diag_indices(len(H))] += 1
    Qw = W(Q); del B, Q
    em_Z = np.zeros(len(H), bool)
    for i, r in enumerate(rs):
        if abs(r.m) < MZ:
            for c, ns in r.comps:
                o0, _ = r.off[c]
                em_Z[offs[i] + o0 + np.where(ns < NZ)[0]] = True
    iZ = np.where(em_Z)[0]; iT = np.where(~em_Z)[0]
    V = np.linalg.inv(H[np.ix_(iZ, iZ)])
    print(f'[gamma1 livre: {a.gamma1_livre}] L em s0 = {a.s0}: G\' = {MG}x{NG} ({len(H)}), Z = {MZ}x{NZ} ({len(iZ)}), T ({len(iT)})  [{time.time()-t0:.0f}s]', flush=True)
    nV = norma2(V); VQ = norma2(V @ Qw[np.ix_(iZ, iZ)])
    A_ = norma2(V @ H[np.ix_(iZ, iT)]); B_ = norma2(H[np.ix_(iT, iZ)])
    X = H[np.ix_(iT, iT)]; X[np.diag_indices(len(iT))] -= 1
    alfa = norma2(X)
    # separacao da cauda: Fourier (|m| >= MZ) e radial (|m| < MZ, n >= NZ)
    mT = np.concatenate([np.full(r.dim, abs(r.m)) for r in rs])[iT]
    iF = np.where(mT >= MZ)[0]; iR = np.where(mT < MZ)[0]
    blocos = {nome: norma2(X[np.ix_(I, J)]) for nome, (I, J) in
              dict(FF=(iF, iF), RR=(iR, iR), FR=(iF, iR), RF=(iR, iF)).items()}
    print('  alfa por bloco (linhas<-colunas; F Fourier, R radial): ' +
          '  '.join(f'{k} {v:.4f}' for k, v in blocos.items()))
    print(f'  ||V|| {nV:.2f}  ||V Q_ZZ|| {VQ:.2f}  a {A_:.4f}  b {B_:.4f}  alfa {alfa:.4f}'
          f'  ->  Schur ab/(1-alfa) = {A_*B_/(1-alfa) if alfa < 1 else float("inf"):.4f}   [{time.time()-t0:.0f}s]')


if __name__ == '__main__':
    main()
