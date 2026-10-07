#!/usr/bin/env python3
"""DIAGNOSTICO EM PONTO FLUTUANTE -- NAO E PROVA. Viabilidade do certificado de
Newton-Kantorovich (NK) BORDEJADO de L no autovalor ~0,401 (S3b, rota local).

Incognitas (lambda, y), h = Q0 y, Q0 = (O_mu + lambda0)^-1; F = (L(lambda) h, l(h) - 1),
L(s) = L0 + s. Jacobiano pre-condicionado (pesos do toro (65/64, 5/4), m com sinal):

    M_b = [[H, h0], [l Q0, 0]],   H = L(lambda0) Q0 = I + B Q0.

Niveis: Z' = Z + {lambda} (bloco bordejado invertido: V_b) e T = G' \\ Z (banda).
Mede ||V_b||, a = ||V_b [H_ZT; l Q0_T]||, b = ||H_TZ||, alfa = ||H_TT - I||, o criterio
ab/(1 - alfa) e o residuo r = L(lambda0) h0 fora de Z (h0 = autovetor da truncagem Z).
Raio NK previsto: eps_y ~ ||A r||/(1 - theta).
"""
import argparse, math, sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import signed_operator_L as sl
from tile_certificate import norma2


def omega2(n):
    return np.where(n == 0, 1.0, sl.K2**(2.0*n) + sl.K2**(-2.0*n))


def chaves(rs, offs):
    """(m, componente, n) -> indice global."""
    out = {}
    for i, r in enumerate(rs):
        for c, ns in r.comps:
            o0, _ = r.off[c]
            for k, n in enumerate(ns):
                out[(r.m, c, int(n))] = offs[i] + o0 + k
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--alvo', type=float, default=0.401)
    ap.add_argument('--Z', default='24x96'); ap.add_argument('--G', default='32x128')
    a = ap.parse_args()
    MZ, NZ = map(int, a.Z.split('x')); MG, NG = map(int, a.G.split('x'))
    t0 = time.time()
    g = sl.campos()
    # autopar da truncagem Z: (L0 + s) h = 0, s perto do alvo
    LZ, rZ, oZ = sl.operador(MZ, NZ, 0.0, g=g)
    ev, evec = np.linalg.eig(LZ)
    k = np.argmin(np.abs(-ev - a.alvo))
    lam0 = complex(-ev[k]); lam0 = lam0.real if abs(lam0.imag) < 1e-10 else lam0
    hZ = evec[:, k]
    del LZ, evec
    print(f'autovalor da truncagem Z = {MZ}x{NZ}: lambda0 = {lam0:.10f}  [{time.time()-t0:.0f}s]', flush=True)
    # operador na caixa G', no lugar: L0 -> B -> B Q0 -> H pesado
    H, rs, offs = sl.operador(MG, NG, 0.0, g=g)
    ch = chaves(rs, offs); chZ = chaves(rZ, oZ)
    h0 = np.zeros(len(H), complex)
    for key, j in chZ.items():
        h0[ch[key]] = hZ[j]
    Lh = H @ h0 + lam0*h0                      # residuo L(lambda0) h0, antes de sobrescrever H
    Qb = []
    for i, r in enumerate(rs):
        blk = slice(offs[i], offs[i+1])
        J = sl.J_livre(r, 0.0)
        H[blk, blk] -= J
        Qb.append(np.linalg.inv(J + lam0*np.eye(r.dim)))
    for i in range(len(rs)):
        blk = slice(offs[i], offs[i+1])
        H[:, blk] = H[:, blk] @ Qb[i]
    peso = np.concatenate([sl.K1**r.m*np.sqrt(np.concatenate([omega2(ns) for _, ns in r.comps])) for r in rs])
    H *= peso[:, None]; H /= peso[None, :]
    H[np.diag_indices(len(H))] += 1
    Qw = [(peso[offs[i]:offs[i+1], None]*Qb[i])/peso[None, offs[i]:offs[i+1]] for i in range(len(rs))]
    h0w = peso*h0; h0w /= np.linalg.norm(h0w)
    em_Z = np.zeros(len(H), bool); em_Z[[ch[key] for key in chZ]] = True
    iZ = np.where(em_Z)[0]; iT = np.where(~em_Z)[0]
    # residuo r = L(lambda0) h0 = (B + O_mu + lambda0) h0: no pesado, r_w = W r
    rw = peso*Lh/np.linalg.norm(peso*h0)
    rZn, rTn = np.linalg.norm(rw[iZ]), np.linalg.norm(rw[iT])
    print(f'  G\' = {MG}x{NG} ({len(H)}), Z = {len(iZ)}, T = {len(iT)};  residuo relativo: em Z {rZn:.2e}, fora de Z {rTn:.2e}'
          f'  [{time.time()-t0:.0f}s]', flush=True)
    # funcional l(h) = <h0w, W h> (pesado); linha l Q0 em coordenadas y pesadas
    lQ = np.zeros(len(H), complex)
    for i, r in enumerate(rs):
        blk = slice(offs[i], offs[i+1])
        lQ[blk] = h0w[blk].conj() @ Qw[i]
    # bloco Z bordejado e V_b
    nZ = len(iZ)
    Mb = np.zeros((nZ + 1, nZ + 1), complex)
    Mb[:nZ, :nZ] = H[np.ix_(iZ, iZ)]
    Mb[:nZ, nZ] = h0w[iZ]
    Mb[nZ, :nZ] = lQ[iZ]
    Vb = np.linalg.inv(Mb)
    nVb = norma2(Vb)
    C = np.vstack([H[np.ix_(iZ, iT)], lQ[iT][None, :]])
    a_b = norma2(Vb @ C)
    b = norma2(H[np.ix_(iT, iZ)])
    X = H[np.ix_(iT, iT)]; X[np.diag_indices(len(iT))] -= 1
    alfa = norma2(X)
    crit = a_b*b/(1 - alfa) if alfa < 1 else float('inf')
    print(f'  ||V_b|| {nVb:.2f}   a_b {a_b:.4f}   b {b:.4f}   alfa {alfa:.4f}   ->  ab/(1-alfa) = {crit:.4f}'
          f'  [{time.time()-t0:.0f}s]', flush=True)
    if crit < 1:
        # Perron 2x2 [[0, a],[b, alfa]] com o residuo (Z: V_b r_Z; T: r_T)
        N = np.array([[0.0, a_b], [b, alfa]])
        lam, vec = np.linalg.eig(N); th = max(abs(lam))
        rZb = norma2(Vb[:, :nZ] @ rw[iZ][:, None]) if rZn > 0 else 0.0
        eps = max(rZb, rTn)/(1 - th)
        print(f'  theta {th:.4f}; raio NK previsto (norma de Perron, y) ~ {eps:.2e}')


if __name__ == '__main__':
    main()
