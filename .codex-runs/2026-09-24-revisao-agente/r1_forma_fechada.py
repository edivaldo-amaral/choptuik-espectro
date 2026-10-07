#!/usr/bin/env python3
"""R1 -- forma fechada de (J + sigma)^-1 e teste das cotas cauda_Rinv, Rinv_uniforme, descida.

J (signed_operator.Jrad, um componente): lista ns (passo 1: n = 0, 1, 2, ...; passo 2: n impar),
    diagonal  d_k = mu (n_k + 1) + sigma,
    acima da diagonal, coluna c: v_c = 2 mu n_c em todas as linhas j < c.

Substituicao regressiva de (J + sigma) x = e_c:
    x_c = 1/d_c;  para j < c:  d_j x_j + sum_{k=j+1}^{c} v_k x_k = 0.
Com S_j = sum_{k=j}^{c} v_k x_k:  x_j = -S_{j+1}/d_j  e  S_j = S_{j+1} (1 - v_j/d_j), S_c = v_c/d_c. Logo

    (J + sigma)^-1_{jc} = 1/d_c                                   (j = c)
                        = -(v_c/(d_j d_c)) prod_{j<k<c} (1 - v_k/d_k)   (j < c)
                        = 0                                        (j > c).

E  |1 - v_k/d_k| = |mu(1 - n_k) + sigma| / |mu(1 + n_k) + sigma| <= 1  se Re sigma >= -mu.
"""
from __future__ import annotations

import math
from pathlib import Path
import sys
import time

import numpy as np

RAIZ = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(RAIZ/'scripts'))
import fourier_tail_bound as ft          # cotas sob teste
import tile_perron as tp                 # descida
from signed_operator import Jrad, radial

MU = ft.MU
K2 = ft.K2


def lista(N, passo):
    return np.arange(N) if passo == 1 else np.arange(1, N, 2)


def forma_fechada(ns, sigma, mu=MU):
    """Matriz (J + sigma)^-1 pela formula fechada (independente de np.linalg.inv)."""
    n = len(ns)
    d = mu*(ns + 1.0) + sigma
    v = 2*mu*ns.astype(float)
    f = 1 - v/d
    R = np.zeros((n, n), complex)
    for c in range(n):
        R[c, c] = 1/d[c]
        if c == 0:
            continue
        # p_j = prod_{k=j+1}^{c-1} f_k, j = c-1, ..., 0
        p = np.ones(c, complex)
        if c >= 2:
            p[:c-1] = np.cumprod(f[1:c][::-1])[::-1]
        R[:c, c] = -(v[c]/(d[:c]*d[c]))*p
    return R


def J_de(ns, sigma, mu=MU):
    J = np.diag(mu*(ns + 1.0)).astype(complex) + sigma*np.eye(len(ns))
    for b in range(len(ns)):
        J[:b, b] += 2*mu*ns[b]
    return J


def lw(ns):
    """log sqrt(omega(n)) sem estouro."""
    ns = np.asarray(ns, float)
    return np.where(ns == 0, 0.0, ns*math.log(K2) + 0.5*np.log1p(K2**(-4.0*ns)))


def pesada(R, ns_lin, ns_col):
    return R*np.exp(lw(ns_lin)[:, None] - lw(ns_col)[None, :])


def parte1():
    print('== R1a: formula fechada contra np.linalg.inv (e contra signed_operator.Jrad)')
    rng = np.random.default_rng(7)
    pior = 0.0
    sigmas = [0.0, 1e-3, 0.25j, 0.2 + 0.1j, 1.765 + 0.25j, 3.0 - 7.0j, 0.0 + 20.0j, MU, 5*MU]
    sigmas += list(rng.uniform(0, 2, 6) + 1j*rng.uniform(-60, 60, 6))
    for passo in (1, 2):
        for N in (1, 2, 3, 5, 10, 25, 40):
            ns = lista(N if passo == 1 else 2*N + 1, passo)
            for sg in sigmas:
                J = J_de(ns, sg)
                Rf = forma_fechada(ns, sg)
                Ri = np.linalg.inv(J)
                err = np.abs(Rf - Ri).max()/max(np.abs(Ri).max(), 1e-300)
                pior = max(pior, err)
                assert err < 1e-11, (passo, N, sg, err)
    # consistencia com o J do certificado: Jrad(m, s, N) = blocos (c1 impar, c2 todo n)
    for (m, s, N) in ((0, 0.0, 24), (-6, 0.25j, 36), (10, 0.3 + 0.1j, 40)):
        Jr = Jrad(m, s, N)
        n1, n2 = radial(N)
        sg = s + 1j*m/2
        Jmeu = np.block([[J_de(n1, sg), np.zeros((len(n1), len(n2)))],
                         [np.zeros((len(n2), len(n1))), J_de(n2, sg)]])
        assert np.abs(Jr - Jmeu).max() == 0.0
    print(f'   pior erro relativo (max entrada) formula/inv: {pior:.2e}  (passos 1 e 2, n ate 40/81, '
          f'{len(sigmas)} sigmas)  -> OK')
    # |1 - v/d| <= 1 para Re sigma >= 0
    worst = 0.0
    for sg in list(rng.uniform(0, 3, 200) + 1j*rng.uniform(-200, 200, 200)) + [0.0, 1e-15]:
        ns = np.arange(0, 3000)
        f = np.abs(1 - 2*MU*ns/(MU*(ns + 1.0) + sg))
        worst = max(worst, f.max())
    print(f'   max |1 - v_k/d_k| com Re sigma >= 0: {worst:.15f} (<= 1)')
    return pior


def sigma1(X):
    """maior valor singular: Lanczos (svds) -- valor de Ritz, logo cota INFERIOR da verdade."""
    from scipy.sparse.linalg import svds
    if min(X.shape) <= 400:
        return float(np.linalg.norm(X, 2))
    try:
        return float(svds(X, k=1, return_singular_vectors=False, tol=1e-10)[0])
    except Exception:
        return float(np.linalg.norm(X, 2))


def pesada_cheia(sigma, passo, Ntr):
    ns = lista(Ntr, passo)
    return ns, pesada(forma_fechada(ns, sigma), ns, ns)


def parte2(Ntr=1400):
    print(f'\n== R1b: cotas da forma fechada contra a norma verdadeira pesada (truncagem {Ntr})')
    print('   (colunas N <= n < Ntr com linhas completas: J e triangular superior, entao a norma truncada')
    print('    e cota INFERIOR da verdadeira; uma violacao aqui seria real)')
    linhas = []
    pior = 0.0
    for s in (0.0, 0.25j, 0.2 + 0.1j, 0.401 + 0.125j, 1.0 + 0.25j, 1.765):
        for m in (0, -2, 2, -20, 38, -100, 198):
            sg = s + 1j*m/2
            for passo in (1, 2):
                ns, W = pesada_cheia(sg, passo, Ntr)
                for N in (24, 80, 160, 400):
                    cota = ft.cauda_Rinv(complex(sg), N, passo, MU)
                    v = sigma1(W[:, ns >= N])
                    pior = max(pior, v/cota)
                    linhas.append(('cauda_Rinv', s, m, N, passo, v, cota))
                    assert v <= cota, ('VIOLACAO cauda_Rinv', s, m, N, passo, v, cota)
                if m in (0, 198) and s in (0.0, 1.765):
                    pass
    rc = max(l[5]/l[6] for l in linhas)
    print(f'   cauda_Rinv: {len(linhas)} casos, max verdade/cota = {rc:.4f}')
    for l in sorted(linhas, key=lambda l: -l[5]/l[6])[:4]:
        print(f'      s={l[1]}, m={l[2]}, N={l[3]}, passo={l[4]}: verdade {l[5]:.5f}  cota {l[6]:.5f}  razao {l[5]/l[6]:.4f}')
    # Rinv_uniforme: ||R_m^-1|| inteiro
    lu = []
    for s in (0.0, 0.25j, 0.1 + 0.25j, 1.0):
        for m in (20, -40, 60, -120, 200, -202):
            sg = s + 1j*m/2
            cota = ft.Rinv_uniforme(m, MU, s)
            v = max(sigma1(pesada_cheia(sg, 1, Ntr)[1]), sigma1(pesada_cheia(sg, 2, Ntr)[1]))
            lu.append((s, m, v, cota)); pior = max(pior, v/cota)
            assert v <= cota, ('VIOLACAO Rinv_uniforme', s, m, v, cota)
    ru = max(l[2]/l[3] for l in lu)
    print(f'   Rinv_uniforme: {len(lu)} casos, max verdade/cota = {ru:.4f}')
    for l in sorted(lu, key=lambda l: -l[2]/l[3])[:3]:
        print(f'      s={l[0]}, m={l[1]}: verdade {l[2]:.5f}  cota {l[3]:.5f}')
    # descida
    ld = []
    cache = {}
    for s in (0.0, 0.25j, 1.0):
        for m in (0, -2, 18, -38):
            for passo in (1, 2):
                cache[(s, m, passo)] = pesada_cheia(s + 1j*m/2, passo, Ntr)
    for (J, Nc) in ((24, 77), (36, 77), (65, 77), (61, 400), (121, 400), (160, 400), (201, 400), (380, 400)):
        cota = tp.descida(J, Nc)
        v = 0.0
        for (ns, W) in cache.values():
            v = max(v, sigma1(W[np.ix_(ns < J, ns >= Nc)]))
        ld.append((J, Nc, v, cota)); pior = max(pior, v/cota)
        assert v <= cota, ('VIOLACAO descida', J, Nc, v, cota)
    print('   descida(J, Nc): ' + '; '.join(f'({J},{Nc}) verdade {v:.2e} cota {c:.2e}' for J, Nc, v, c in ld))
    return pior, rc, ru, ld


def parte3():
    """Rederivacao independente (em Python puro) das tres cotas, conferida contra as do codigo."""
    print('\n== R1c: rederivacao independente das cotas (Schur com majorantes da forma fechada)')
    cw = math.sqrt(1 + K2**-4)
    # Rinv_uniforme: |d_k| >= h = |m|/2 - |Im s|; |v_c/d_c| < 2; prod <= 1; pesos <= cw K2^(n_j - n_c)
    # coluna: 1/h + (1/h) 2 cw sum_{k>=1} K2^-k = (1 + 2 cw/(K2 - 1))/h ; linha idem.
    for m, s in ((200, 0.125j), (60, 0.25j), (-400, 0.0)):
        h = abs(m)/2 - abs(complex(s).imag)
        meu = (1 + 2*cw/(K2 - 1))/h
        assert abs(meu - ft.Rinv_uniforme(m, MU, s)) < 1e-15*meu
    # descida: |ent| <= (2 cw/(mu (n_j+1))) K2^(n_j - n_c); colunas <= (2cw/mu) K2^(J-Nc) H_J;
    # linhas <= (2cw/mu) K2^(J-Nc) K2/(K2-1)
    for J, Nc in ((121, 400), (201, 400)):
        HJ = sum(1/(j + 1) for j in range(J))
        base = (2*cw/MU)*K2**(J - Nc)
        meu = math.sqrt(base*HJ*base*K2/(K2 - 1))
        assert abs(meu*1.0001 - tp.descida(J, Nc)) < 1e-12*meu
    print('   Rinv_uniforme e descida: formulas conferem com a minha derivacao (ver RELATORIO)')


if __name__ == '__main__':
    t0 = time.time()
    parte1()
    pior, rc, ru, ld = parte2()
    parte3()
    print(f'\nRESUMO R1: nenhuma violacao; pior verdade/cota global {pior:.4f}  [{time.time()-t0:.0f}s]')
