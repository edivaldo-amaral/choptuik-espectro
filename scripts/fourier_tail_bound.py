#!/usr/bin/env python3
"""Cota da cauda infinita beta_f = ||(B + ds) Q0 P_far|| do operador de constraints K,
na norma do toro (129/128, 9/8), modo a modo em Fourier.

Q0 = (K_livre + s0)^-1 e diagonal em m (m par, com sinal: o setor A complexificado).
O fundo se decompoe em modos B = sum_d S_d B_d; o bloco (m+d, m) de B Q0 e a
composicao RADIAL B_d R_m^-1, R_m = J_rad + s0 + i m/2. Teste de Schur nos blocos
de Fourier (o peso do toro da o fator kappa1^|d|):

    ||B Q0 P_far|| <= sum_d kappa1^|d| sup_m ||B_d R_m^-1 P_far,m||.

far = {|m| >= M+} (todo n) U {|m| < M+, n >= N+}.

Cada composicao: truncamento radial em Nr, mais a cauda ||B_d|| ||R_m^-1 P_{>Nr}||.
  * ||B_d||: Young por modo (norma l1 radial com pesos simetricos), estrutura 2x3.
  * ||R_m^-1 P_{>N}||: teste de Schur com majorantes da forma fechada
    |(R^-1)_jn| <= (1/|d_j|)(v_n/|d_n|) prod |1 - v_k/d_k|, fatores <= 1 se Re s0 >= 0.
  * |m| >= Mmax: ||R_m^-1|| <= (2 + 4 c_w/(1-1/k))/|m| uniforme, vezes ||B_d||.
  * |d| > D: produto ||B_d|| sup ||R^-1 P_far||.

ESTADO: prototipo; normas finitas por cota superior verificada (Cholesky) em ponto
flutuante com margem.
"""
from __future__ import annotations

import math
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from component_exterior import read_fields
from tile_certificate import norma2_sup

REF = Path('.cache/rt-1203.3766v1/sourcecode/RefA.dat')
MU = 722873400/2**32
MU_MIN, MU_MAX = 1/6, 0.17
K1, K2 = 129/128, 9/8


def modos_campos():
    """g[j][d] = coeficientes simetricos (n = 0..100) do modo de Fourier d >= 0 do campo j."""
    F = read_fields(REF)
    g = []
    for f in F[:3]:
        dd = {}
        for (m, n, part), v in f.items():
            dd.setdefault(m, np.zeros(101, complex))[n] += float(v)*(1j if part else 1)
        g.append(dd)
    return g


def modo(g, j, d):
    a = g[j].get(abs(d))
    if a is None:
        return np.zeros(101, complex)
    return a.conj() if d < 0 else a


def refl(a):
    return a*(-1.0)**np.arange(len(a))


def combos(g, d):
    v1, v2, v3 = (modo(g, j, d) for j in range(3))
    P1, P2, P3 = refl(v1), refl(v2), refl(v3)
    return ((v2+P2-P3-v3, v1-v2+P2), (-v1+P3-v3, -v1+v2+3*P2), (v2+P2+P3+v3, -v1+v2+3*P2))


def mult(a, nin, nout):
    """multiplicacao radial por a (coeficientes simetricos), entradas nin -> saidas nout."""
    L = len(a)
    D = nout[:, None] - nin[None, :]
    S = nout[:, None] + nin[None, :]
    A = np.where(np.abs(D) < L, a[np.minimum(np.abs(D), L-1)], 0)
    Hk = np.where((S < L) & (nin[None, :] > 0), a[np.minimum(S, L-1)], 0)
    return A + Hk


def omega(j):
    return np.where(j == 0, 1.0, K2**(2.0*j) + K2**(-2.0*j))


def B_rad(g, d, n1i, n2i, n1o, n2o, mu):
    """Bloco radial B_d: entradas (c1 impar n1i, c2 n2i) -> saidas (c1 n1o, c2 n2o)."""
    F1, F2, F3 = combos(g, d)
    ev = (n2i % 2 == 0)
    b11 = 0.5*mult(F1[0], n1i, n1o)
    b12 = 0.5*mult(F2[0], n2i, n1o)*ev[None, :] + 0.5*mult(F3[0], n2i, n1o)*(~ev)[None, :]
    b21 = 0.5*mult(F1[1], n1i, n2o)
    b22 = 0.5*mult(F2[1], n2i, n2o)*ev[None, :] + 0.5*mult(F3[1], n2i, n2o)*(~ev)[None, :]
    if d == 0:
        pos = {n: i for i, n in enumerate(n2o)}
        for c, n in enumerate(n1i):
            k = (n - 1)//2
            for jj in range(k + 1):
                if 2*jj in pos:
                    b21[pos[2*jj], c] += mu*2*(-1)**(k - jj)
    return np.block([[b11, b12], [b21, b22]])


def norma_Bd_young(g, d, mu):
    """||B_d|| <= norma espectral da matriz 2x3 de normas l1 radiais (pesos simetricos)."""
    def l1(a):
        n = np.arange(len(a))
        return float(np.sum(np.abs(a)*np.where(n == 0, 1.0, 2*K2**n)))
    F1, F2, F3 = combos(g, d)
    M = np.array([[0.5*l1(F1[0]), 0.5*l1(F2[0]), 0.5*l1(F3[0])],
                  [0.5*l1(F1[1]), 0.5*l1(F2[1]), 0.5*l1(F3[1])]])
    if d == 0:
        M[1, 0] += mu*2*(1/K2)/(1 - 1/K2**2)
    return float(np.linalg.norm(M, 2))*1.0001


_RCACHE = {}


def livre_inv(ns, passo, sigma, mu):
    chave = (len(ns), passo, sigma, mu)
    if chave in _RCACHE:
        return _RCACHE[chave]
    J = np.diag(mu*(ns + 1.0)).astype(complex) + sigma*np.eye(len(ns))
    for b in range(len(ns)):
        J[:b, b] += 2*mu*ns[b]
    R = np.linalg.inv(J)
    res = np.linalg.norm(J @ R - np.eye(len(ns)), 2)
    assert res < 1e-8, res
    _RCACHE[chave] = R
    return R


from functools import lru_cache


@lru_cache(maxsize=None)
def cauda_Rinv(sigma, N, passo, mu, Nbig=3000):
    """Cota rigorosa (Schur com majorantes da forma fechada) de ||R^-1 P_{n >= N}||,
    radial, para um componente (passo 1: todo n; passo 2: n impar)."""
    ns = np.arange(Nbig + 1) if passo == 1 else np.arange(1, 2*Nbig + 2, 2)
    dk = mu*(ns + 1.0) + sigma
    vk = 2*mu*ns
    rho = np.abs(1 - vk/dk)*(1 + 1e-12)
    rho[0] = 1.0
    # log sqrt(omega(n)) sem estouro: omega(6001) = (9/8)^12002 passa de 1e308, e a
    # razao inf/inf virava nan -- que max() descartava em silencio
    lw = np.where(ns == 0, 0.0, ns*math.log(K2) + 0.5*np.log1p(K2**(-4.0*ns)))
    logr = np.log(np.maximum(rho, 1e-300))
    S = np.cumsum(logr)                     # S[i] = sum_{k<=i} log rho_k
    cols = np.where(ns >= N)[0]
    Cmax = 0.0
    rows_acc = np.zeros(len(ns))
    for c in cols:
        j = np.arange(c + 1)
        prod = np.exp(np.where(j < c, S[c-1] - S[j], 0.0)) if c > 0 else np.ones(1)
        ent = (1/np.abs(dk[j]))*(vk[c]/np.abs(dk[c]))*prod
        ent[c] = 1/np.abs(dk[c])
        wt = ent*np.exp(lw[j] - lw[c])
        cs = float(wt.sum())
        assert math.isfinite(cs), (sigma, N, passo, c)
        Cmax = max(Cmax, cs)
        rows_acc[:c+1] += wt
    # cauda n > Nbig: majorante uniforme (fatores <= 1), decrescente em n
    cw = math.sqrt(1 + K2**-4)
    nb = ns[-1]
    Ctail = (1/(mu*(nb + 1)) + 2*cw*(1/((nb/2 + 1)*(K2 - 1)) + K2**(-nb/2)*(1 + math.log(nb + 1)))/mu)
    Cmax = max(Cmax, Ctail)
    jj = np.arange(len(ns))
    rows_acc += (2/(mu*(ns + 1)))*cw*K2**(ns - nb)/(K2 - 1)
    assert np.all(np.isfinite(rows_acc))
    Rmax = max(float(rows_acc.max()), 2*cw*K2/(mu*(nb + 1)*(K2 - 1)))
    return math.sqrt(Cmax*Rmax)*1.0001


def Rinv_uniforme(m, mu, s0=0.0):
    """||R_m^-1|| uniforme em n: fatores do produto <= 1 (Re s0 >= 0), |d_k| >= h,
    v_n/|d_n| < 2; Schur: (1 + 2 c_w/(k-1))/h, h = |m|/2 - |Im s0|."""
    cw = math.sqrt(1 + K2**-4)
    h = abs(m)/2 - abs(complex(s0).imag)
    assert h > 0
    return (1 + 2*cw/(K2 - 1))/h


def composicao(g, d, m, s0, far_n, Nr, mu):
    """||B_d R_m^-1 P_far|| (sem o fator de Fourier), truncando em Nr + cauda."""
    sigma = s0 + 1j*m/2
    n1 = np.arange(1, Nr, 2); n2 = np.arange(Nr)
    R1 = livre_inv(n1, 2, sigma, mu); R2 = livre_inv(n2, 1, sigma, mu)
    c1 = n1 >= far_n; c2 = n2 >= far_n
    n1o = np.arange(1, Nr + 101, 2); n2o = np.arange(Nr + 101)
    Bd = B_rad(g, d, n1, n2, n1o, n2o, mu)
    Rp = np.block([[R1[:, c1], np.zeros((len(n1), c2.sum()))],
                   [np.zeros((len(n2), c1.sum())), R2[:, c2]]])
    X = Bd @ Rp
    wo = np.sqrt(np.concatenate([omega(n1o), omega(n2o)]))
    wi = np.sqrt(np.concatenate([omega(n1[c1]), omega(n2[c2])]))
    Xw = (wo[:, None]*X)/wi[None, :]
    finito = norma2_sup(Xw) if Xw.size else 0.0
    cauda = norma_Bd_young(g, d, mu)*max(cauda_Rinv(sigma, Nr, 1, mu), cauda_Rinv(sigma, Nr, 2, mu))
    return finito + cauda


def norma_Rinv_far(m, s0, far_n, Nr, mu):
    sigma = s0 + 1j*m/2
    n1 = np.arange(1, Nr, 2); n2 = np.arange(Nr)
    a = []
    for ns, passo in ((n1, 2), (n2, 1)):
        R = livre_inv(ns, passo, sigma, mu)
        c = ns >= far_n
        w = np.sqrt(omega(ns))
        X = (w[:, None]*R[:, c])/w[c][None, :]
        a.append((norma2_sup(X) if X.size else 0.0) + cauda_Rinv(sigma, Nr, passo, mu))
    return max(a)


def beta_far(s0, r, Mp, Np, Nr=500, Mmax=260, D=6, passo_m=2, mu=MU, log=print):
    """Cota de ||(B + ds) Q0 P_far||, far = complemento de G+ = (Mp, Np)."""
    g = modos_campos()
    ms = [m for m in range(-Mmax, Mmax + 1, passo_m)]
    total = 0.0
    supQ = 0.0
    for d in range(-D, D + 1, 2):
        sup_d = 0.0
        for m in ms:
            far_n = 0 if abs(m) >= Mp else Np
            if abs(m) < Mp and far_n >= Nr:
                v = norma_Bd_young(g, d, mu)*max(cauda_Rinv(s0 + 1j*m/2, far_n, 1, mu),
                                                  cauda_Rinv(s0 + 1j*m/2, far_n, 2, mu))
            else:
                v = composicao(g, d, m, s0, far_n, Nr, mu)
            sup_d = max(sup_d, v)
        sup_d = max(sup_d, norma_Bd_young(g, d, mu)*Rinv_uniforme(Mmax + passo_m, mu, s0))
        total += K1**abs(d)*sup_d
        log(f'   d={d:+d}: sup_m = {sup_d:.4f}')
    resto = sum(K1**abs(d)*norma_Bd_young(g, d, mu) for d in range(-40, 41, 2) if abs(d) > D)
    for m in ms:
        far_n = 0 if abs(m) >= Mp else Np
        supQ = max(supQ, norma_Rinv_far(m, s0, far_n, Nr, mu) if not (abs(m) < Mp and far_n >= Nr)
                   else max(cauda_Rinv(s0 + 1j*m/2, far_n, 1, mu), cauda_Rinv(s0 + 1j*m/2, far_n, 2, mu)))
    supQ = max(supQ, Rinv_uniforme(Mmax + passo_m, mu, s0))
    beta = total + resto*supQ + r*supQ
    log(f'   soma |d|<={D}: {total:.4f}; |d|>{D}: {resto*supQ:.2e}; r*||Q0 P_far||: {r*supQ:.4f} (||Q0 P_far|| <= {supQ:.4f})')
    return beta, supQ
