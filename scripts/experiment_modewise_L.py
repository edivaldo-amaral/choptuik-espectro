#!/usr/bin/env python3
"""DIAGNOSTICO EM PONTO FLUTUANTE -- NAO E PROVA. A cota modo a modo da cauda de L
(S3b): alfa_modo = sum_d kappa1^d sup_m ||B_d R_m^-1 P_T|| (linhas restritas a T), para
T = F \\ G, G = (M, N), F = (M + 40, N + 150), no autovalor lambda0 ~ 0,401. Compare com o
alfa verdadeiro da banda (experiment_bordered_L.py) para medir a perda da soma sobre d."""
import argparse, sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import signed_operator_L as sl
from tile_certificate import norma2


def om(n):
    return np.where(n == 0, 1.0, sl.K2**(2.0*n) + sl.K2**(-2.0*n))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--G', default='20x80'); ap.add_argument('--D', type=int, default=10)
    ap.add_argument('--lam', type=float, default=0.4010247330)
    ap.add_argument('--janela', type=int, default=24, help='so |m| < M + janela (o sup mora perto da borda)')
    a = ap.parse_args()
    M, N = map(int, a.G.split('x')); Mc, Nc = M + 40, N + 150
    g = sl.campos(); t0 = time.time()
    sup = {}
    for m in range(-(M + a.janela - 1), M + a.janela):
        rin = sl.Radial(m, Nc)
        J = sl.J_livre(rin, 0.0) + a.lam*np.eye(rin.dim)
        win = np.sqrt(np.concatenate([om(ns) for _, ns in rin.comps]))
        R = (win[:, None]*np.linalg.inv(J))/win[None, :]
        nin = np.concatenate([ns for _, ns in rin.comps])
        cols = np.where(nin >= N)[0] if abs(m) < M else np.arange(rin.dim)
        RT = R[:, cols]
        for d in range(-a.D, a.D + 1):
            mo = m + d
            if abs(mo) >= Mc:
                continue
            rout = sl.Radial(mo, Nc)
            Bd = sl.bloco_B(g, rout, rin)
            if not np.any(Bd):
                continue
            wout = np.sqrt(np.concatenate([om(ns) for _, ns in rout.comps]))
            X = sl.K1**d*((wout[:, None]*Bd)/win[None, :]) @ RT
            nout = np.concatenate([ns for _, ns in rout.comps])
            linhas = np.where(nout >= N)[0] if abs(mo) < M else np.arange(rout.dim)
            v = norma2(X[linhas])
            sup[d] = max(sup.get(d, 0.0), v)
    tot = sum(sup.values())
    print(f'G = {M}x{N}: alfa modo a modo (|d| <= {a.D}) = {tot:.4f}   [{time.time()-t0:.0f}s]')
    print('   por d: ' + ' '.join(f'{d:+d}:{v:.3f}' for d, v in sorted(sup.items())))


if __name__ == '__main__':
    main()
