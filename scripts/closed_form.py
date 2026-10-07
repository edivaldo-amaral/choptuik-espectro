#!/usr/bin/env python3
"""Cotas da forma fechada de (J + sigma)^-1, PARAMETRIZADAS em kappa2 (o raio radial do toro),
para o certificado de L (S3b, kappa2 = 5/4). As versoes de K (kappa2 = 9/8) ficam em
fourier_tail_bound.py e tile_perron.py, auditadas; este modulo as reproduz exatamente com
kappa2 = 9/8 (teste de regressao em main).

J: triangular superior na lista de n de uma componente (passo 1: todo n; passo 2: n impar),
diagonal d_k = mu(n_k + 1) + sigma, acima da diagonal v_c = 2 mu n_c na coluna c. Forma fechada:
    (J + sigma)^-1_{jc} = -(1/d_j)(v_c/d_c) prod_{j<k<c} (1 - v_k/d_k),   diagonal 1/d_c,
com |1 - v_k/d_k| <= 1 se Re sigma >= 0. Pesos sqrt(omega(n)), omega(n) = k2^2n + k2^-2n (n > 0),
e sqrt(omega(j)/omega(c)) <= c_w k2^(j - c), c_w = sqrt(1 + k2^-4).
"""
from __future__ import annotations

from functools import lru_cache
import math

import numpy as np


def omega(j, k2):
    return np.where(j == 0, 1.0, k2**(2.0*j) + k2**(-2.0*j))


@lru_cache(maxsize=None)
def cauda_Rinv(sigma, N, passo, mu, k2, Nbig=3000):
    """Cota (Schur com majorantes da forma fechada) de ||R^-1 P_{n >= N}|| pesada, radial, numa
    componente. Identica a fourier_tail_bound.cauda_Rinv com k2 = 9/8."""
    ns = np.arange(Nbig + 1) if passo == 1 else np.arange(1, 2*Nbig + 2, 2)
    dk = mu*(ns + 1.0) + sigma
    vk = 2*mu*ns
    rho = np.abs(1 - vk/dk)*(1 + 1e-12)
    rho[0] = 1.0
    lw = np.where(ns == 0, 0.0, ns*math.log(k2) + 0.5*np.log1p(k2**(-4.0*ns)))
    logr = np.log(np.maximum(rho, 1e-300))
    S = np.cumsum(logr)
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
        assert math.isfinite(cs)
        Cmax = max(Cmax, cs)
        rows_acc[:c+1] += wt
    cw = math.sqrt(1 + k2**-4)
    nb = ns[-1]
    Ctail = (1/(mu*(nb + 1)) + 2*cw*(1/((nb/2 + 1)*(k2 - 1)) + k2**(-nb/2)*(1 + math.log(nb + 1)))/mu)
    Cmax = max(Cmax, Ctail)
    rows_acc += (2/(mu*(ns + 1)))*cw*k2**(ns - nb)/(k2 - 1)
    assert np.all(np.isfinite(rows_acc))
    Rmax = max(float(rows_acc.max()), 2*cw*k2/(mu*(nb + 1)*(k2 - 1)))
    return math.sqrt(Cmax*Rmax)*1.0001


def cauda_Rinv_disco(sigma_p, r, N, passo, mu, k2, Nbig=3000):
    """Como cauda_Rinv, mas cota sup ||R^-1 P_{n >= N}|| sobre sigma = sigma_p + e, |e| <= r, Re sigma >= 0 (rouche_L).
    Entradas da forma fechada: |(R^-1)_{jc}| = (1/|d_j|)(v_c/|d_c|) prod_{j<k<c} |1 - v_k/d_k|, d_k = mu(n_k + 1) + sigma,
    v_k = 2 mu n_k, todas crescentes em 1/|d_k| e em |1 - v_k/d_k|. No disco:
      |d_k| >= max(|d_k(p)| - r, mu(n_k + 1))                       (Re sigma >= 0),
      |1 - v_k/d_k| = |d_k - v_k|/|d_k| <= min(1, (|d_k(p) - v_k| + r)/|d_k|_inf),
    e |d_k|^2 - |d_k - v_k|^2 = 2 mu n_k (2 mu + 2 Re sigma) >= 0 da o 1. A cauda n > Nbig e os majorantes de linha/
    coluna sao os de cauda_Rinv (so dependem de mu, k2 e de Re sigma >= 0)."""
    ns = np.arange(Nbig + 1) if passo == 1 else np.arange(1, 2*Nbig + 2, 2)
    dkp = mu*(ns + 1.0) + sigma_p
    vk = 2*mu*ns
    dinf = np.maximum(np.abs(dkp) - r, mu*(ns + 1.0))*(1 - 1e-12)
    rho = np.minimum(1.0, (np.abs(dkp - vk) + r)/dinf)*(1 + 1e-12)
    rho[0] = 1.0
    lw = np.where(ns == 0, 0.0, ns*math.log(k2) + 0.5*np.log1p(k2**(-4.0*ns)))
    S = np.cumsum(np.log(np.maximum(rho, 1e-300)))
    Cmax = 0.0
    rows_acc = np.zeros(len(ns))
    for c in np.where(ns >= N)[0]:
        j = np.arange(c + 1)
        prod = np.exp(np.where(j < c, S[c-1] - S[j], 0.0)) if c > 0 else np.ones(1)
        ent = (1/dinf[j])*(vk[c]/dinf[c])*prod
        ent[c] = 1/dinf[c]
        wt = ent*np.exp(lw[j] - lw[c])
        cs = float(wt.sum())
        assert math.isfinite(cs)
        Cmax = max(Cmax, cs)
        rows_acc[:c+1] += wt
    cw = math.sqrt(1 + k2**-4)
    nb = ns[-1]
    Ctail = (1/(mu*(nb + 1)) + 2*cw*(1/((nb/2 + 1)*(k2 - 1)) + k2**(-nb/2)*(1 + math.log(nb + 1)))/mu)
    Cmax = max(Cmax, Ctail)
    rows_acc += (2/(mu*(ns + 1)))*cw*k2**(ns - nb)/(k2 - 1)
    assert np.all(np.isfinite(rows_acc))
    Rmax = max(float(rows_acc.max()), 2*cw*k2/(mu*(nb + 1)*(k2 - 1)))
    return math.sqrt(Cmax*Rmax)*1.0001


def Rinv_uniforme(m, mu, s0, k2):
    """||R_m^-1|| <= (1 + 2 c_w/(k2 - 1))/h, h = |m|/2 - |Im s0| > 0 (Re s0 >= 0); x 1,0001."""
    cw = math.sqrt(1 + k2**-4)
    h = abs(m)/2 - abs(complex(s0).imag)
    assert h > 0
    return (1 + 2*cw/(k2 - 1))/h*1.0001


def descida(J, Nc, mu, k2):
    """||P_{n < J} R^-1 P_{n >= Nc}|| <= (2 c_w/mu) k2^(J - Nc) sqrt(H_J k2/(k2 - 1)); x 1,0001."""
    cw = math.sqrt(1 + k2**-4)
    HJ = sum(1/(j + 1) for j in range(J))
    base = (2*cw/mu)*k2**(J - Nc)
    return math.sqrt(base*HJ*base*k2/(k2 - 1))*1.0001


def descida_divxi(J, Nc, k2, fator=1.0):
    """||P_{n < J} (fator mu R) Q0 P_{n >= Nc}|| para uma divisao regularizada por xi R (T_n/xi ->
    2(-1)^L em n - 1 - 2L), com |entrada pesada de mu R| <= 2 mu c_w k2^(j - n). Composicao
    (j <- c) <= 4 c_w^2 fator (1 + ln(c+1)) k2^(j - c); Schur como em tile_perron.descida_rx."""
    cw = math.sqrt(1 + k2**-4)
    f = (1 + math.log(Nc + 1))*k2**(-Nc)
    q = (1 + 1/((Nc + 1)*(1 + math.log(Nc + 1))))/k2
    assert q < 1
    base = 4*cw*cw*abs(fator)*f*k2**J
    return math.sqrt(base/(1 - q)*base/(1 - k2**-2))*1.0001


def norma_divxi(k2):
    """||R|| (sem mu) pelo teste de Schur: |entrada| <= 2 c_w k2^(j - n), passo 2 em n:
    linhas e colunas <= 2 c_w k2^-1/(1 - k2^-2)."""
    cw = math.sqrt(1 + k2**-4)
    return 2*cw/(k2*(1 - k2**-2))*1.0001


def _regressao():
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import fourier_tail_bound as ft
    import tile_perron as tp
    k2 = 9/8
    pior = 0.0
    for sig in (0.0, 0.3 + 0.2j, 1.0 + 20j):
        for N in (80, 160, 400):
            for passo in (1, 2):
                a = cauda_Rinv(sig, N, passo, ft.MU, k2); b = ft.cauda_Rinv(sig, N, passo, ft.MU)
                pior = max(pior, abs(a - b)/b)
    pior = max(pior, abs(Rinv_uniforme(200, ft.MU, 0.25j, k2)/1.0001 - ft.Rinv_uniforme(200, ft.MU, 0.25j))/ft.Rinv_uniforme(200, ft.MU, 0.25j))
    pior = max(pior, abs(descida(181, 400, ft.MU, k2) - tp.descida(181, 400))/tp.descida(181, 400))
    pior = max(pior, abs(descida_divxi(80, 400, k2) - tp.descida_rx(80, 400))/tp.descida_rx(80, 400))
    print(f'regressao contra as versoes de K (kappa2 = 9/8): maior diferenca relativa {pior:.2e}')
    assert pior < 1e-12


if __name__ == '__main__':
    _regressao()
