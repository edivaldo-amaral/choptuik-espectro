#!/usr/bin/env python3
"""DIAGNOSTICO EM PONTO FLUTUANTE -- NAO E PROVA. A cauda de Fourier de L (S3b) por JANELAS
de W modos: normas EXATAS dos blocos janela <- janela de P_T B Q0 P_T (T = |m| >= M, n < Nc),
e o raio de Perron da matriz de blocos (e a soma de Schur por deslocamento de janela).
Compare com a banda verdadeira e com a cota modo a modo."""
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
    ap.add_argument('--M', type=int, default=20); ap.add_argument('--Nc', type=int, default=120)
    ap.add_argument('--W', type=int, default=8); ap.add_argument('--nwin', type=int, default=4)
    ap.add_argument('--lam', type=float, default=0.4010247330)
    a = ap.parse_args()
    g = sl.campos(); t0 = time.time()
    # janelas: [M + k W, M + (k+1) W) e as simetricas negativas
    jan = []
    for k in range(a.nwin):
        lo = a.M + k*a.W
        jan.append(list(range(lo, lo + a.W))); jan.append(list(range(-(lo + a.W - 1), -lo + 1)))
    Rs = {}; ws = {}; rads = {}
    for J in jan:
        for m in J:
            r = sl.Radial(m, a.Nc); rads[m] = r
            w = np.sqrt(np.concatenate([om(ns) for _, ns in r.comps])); ws[m] = w
            Jm = sl.J_livre(r, 0.0) + a.lam*np.eye(r.dim)
            Rs[m] = (w[:, None]*np.linalg.inv(Jm))/w[None, :]
    def bloco(I, J):
        di = [rads[m].dim for m in I]; dj = [rads[m].dim for m in J]
        X = np.zeros((sum(di), sum(dj)), complex)
        oi = np.cumsum([0] + di); oj = np.cumsum([0] + dj)
        for q, mj in enumerate(J):
            for p, mi in enumerate(I):
                d = mi - mj
                if abs(d) > sl.BAND_M:
                    continue
                Bd = sl.bloco_B(g, rads[mi], rads[mj])
                if np.any(Bd):
                    Bw = sl.K1**d*((ws[mi][:, None]*Bd)/ws[mj][None, :])
                    X[oi[p]:oi[p+1], oj[q]:oj[q+1]] = Bw @ Rs[mj]
        return X
    n = len(jan)
    Nmat = np.zeros((n, n))
    for i, I in enumerate(jan):
        for j, J in enumerate(jan):
            if min(abs(x - y) for x in I for y in J) > sl.BAND_M:
                continue
            Nmat[i, j] = norma2(bloco(I, J))
    rho = max(abs(np.linalg.eigvals(Nmat)))
    # tudo junto (verdade na banda)
    todos = [m for J in jan for m in J]
    verdade = norma2(bloco(todos, todos))
    print(f'M = {a.M}, W = {a.W}, {a.nwin} janelas de cada lado, Nc = {a.Nc}: verdade (banda) {verdade:.4f};'
          f'  Perron das janelas {rho:.4f}   [{time.time()-t0:.0f}s]')
    np.set_printoptions(precision=3, suppress=True, linewidth=160)
    print(Nmat[:6, :6])


if __name__ == '__main__':
    main()
