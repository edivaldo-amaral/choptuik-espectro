#!/usr/bin/env python3
"""Montagem (estimativa em ponto flutuante) do Perron de 3 níveis do Rouché de L em pontos p do
contorno (C1_REAVALIACAO.md §7). Para cada p, usa o centro de cauda c do seu fator:
  Z' <- T   a_esc(p) (rouche_L_Z.py pontos; escala D_Z = J_Z(c) Q_ZZ(p)) / min_{SZ} w;
  T  <- Z'  TZ(c) (exato na escala);  T <- T  ||W (N(c) + diag rho_J(p)) W^-1|| / (1 - |p - c| q_TT),
            rho_J(p) = kappa |p - c| (sensibilidade medida: kappa <= 0,2, experiment_cauda_s);
  T <- far  max(w nVhi) beta(p) + sqrt(sum (w nV)^2) desc;  far <- T  ||Nf W^-1|| / (1 - |p - c| q_TT);
  far <- far beta(p) = ||B|| ||Q0(p) P_far|| (formas fechadas em p).
Imprime theta(p) com os pesos 1 e sqrt(u/z) do centro."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import closed_form as cf
import nk_L as nk

MU, K2 = nk.MU, nk.K2


def beta_em(p, Mc, Nc):
    rad = max(max(cf.cauda_Rinv(p + 1j*m/2, Nc, q, MU, K2, 1500) for q in ((2, 1, 1) if m % 2 == 0 else (1,)))
              for m in range(-(Mc - 1), Mc))
    return nk.NORMA_BL*max(rad, cf.Rinv_uniforme(Mc, MU, p, K2))


def theta(pz, cauda, p, c, qTT, kappa, Mc, Nc):
    N = np.array(cauda['N']); TZ = np.array(cauda['TZ']); nV = np.array(cauda['nV']); nVhi = np.array(cauda['nVhi'])
    Nf = np.array(cauda['Nf']); desc = cauda['desc']; SZ = cauda['SZ']
    k = len(nV); d = abs(p - c); fq = 1/(1 - d*qTT)
    beta = beta_em(p, Mc, Nc)
    out = {}
    for nomew, r in cauda['res'].items():
        w = np.array(r['w'])
        Nd = N*fq + np.diag(np.full(k, kappa*d))
        nT = float(np.linalg.norm((w[:, None]*Nd)/w[None, :], 2))
        TZw = math.sqrt(float(np.sum((w*TZ)**2)))
        ZT = pz['a_esc']/min(w[j] for j in SZ)
        Tf = max(w*nVhi)*beta + math.sqrt(float(np.sum((w*nV)**2)))*desc
        fT = float(np.linalg.norm(Nf/w[None, :], 2))*fq
        M = np.array([[1e-6, ZT, 1e-9], [TZw, nT, Tf], [0.0, fT, beta]])
        out[nomew] = float(max(abs(np.linalg.eigvals(M))))
    return out, beta


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--pontos', required=True, help='JSON de rouche_L_Z.py pontos')
    ap.add_argument('--fator', required=True, help='JSON do fator (qTT)')
    ap.add_argument('--cauda', required=True, help='JSON de nk_L3.py no centro c')
    ap.add_argument('--kappa', type=float, default=0.2)
    ap.add_argument('--F', default='200x400')
    a = ap.parse_args()
    Mc, Nc = map(int, a.F.split('x'))
    P = json.load(open(a.pontos)); fat = json.load(open(a.fator)); cauda = json.load(open(a.cauda))
    c = complex(*P['c']); qTT = fat['qTT']
    pior = 0.0
    for pz in P['pontos']:
        p = complex(*pz['p'])
        th, beta = theta(pz, cauda, p, c, qTT, a.kappa, Mc, Nc)
        m = min(th.values()); pior = max(pior, m)
        print(f'p = {p.real:.4f}{p.imag:+.4f}i (|p - c| {abs(p - c):.3f}): Z<-T {pz["a_esc"]:.4f}, beta {beta:.4f}; '
              + ', '.join(f'theta[{k}] {v:.4f}' for k, v in th.items()))
    print(f'pior theta (melhor dos pesos por ponto): {pior:.4f}')


if __name__ == '__main__':
    main()
