#!/usr/bin/env python3
"""DIAGNOSTICO EM PONTO FLUTUANTE -- NAO E PROVA. Bloco-Jacobi nas janelas de T para L.

A = P_T B Q0 P_T tem norma 0,637 mas raio espectral pequeno (0,17 no bloco da borda de Z):
e muito nao normal. Com V_J = (I + A_JJ)^-1 nas janelas fortes (as de norma > --limiar) e I nas
demais, o nivel T da certificacao tem E_JJ = 0 (densas) ou -A_JJ, e E_IJ = -V_I A_IJ. Mede a
matriz de Perron por janelas e o seu raio. Blocos densos montados inteiros (F pequeno)."""
import argparse, sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import nk_L as nk
from tile_certificate import norma2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--Z', default='32x128'); ap.add_argument('--F', default='80x200')
    ap.add_argument('--W', type=int, default=16)
    ap.add_argument('--limiar', type=float, default=0.3)
    ap.add_argument('--lam0', type=float, default=0.4010247330)
    a = ap.parse_args()
    MZ, NZ = map(int, a.Z.split('x')); Mc, Nc = map(int, a.F.split('x'))
    ctx = nk.Contexto(a.lam0, Nc)
    t0 = time.time()
    ms = list(range(-(Mc - 1), Mc))
    cols = {m: nk.colsT(ctx, m, MZ, NZ, Nc) for m in ms}
    QT = {m: ctx.Q(m, Nc)[:, cols[m]] for m in ms}
    ctx.Q.cache_clear()
    Bs = {}
    def B(mo, mi):
        k = (mo - mi, mi % 2)
        if k not in Bs:
            Bs[k] = ctx.B(mo, mi, Nc, Nc); ctx._B.cache_clear()
        return Bs[k]
    jan = [ms[k:k+a.W] for k in range(0, len(ms), a.W)]
    nj = len(jan)
    def bloco(I, J):
        mI, mJ = jan[I], jan[J]
        oi = np.cumsum([0] + [len(cols[m]) for m in mI]); oj = np.cumsum([0] + [len(cols[m]) for m in mJ])
        X = np.zeros((int(oi[-1]), int(oj[-1])), complex)
        for p, mo in enumerate(mI):
            for q, mi in enumerate(mJ):
                if abs(mo - mi) <= nk.BAND_M:
                    X[oi[p]:oi[p+1], oj[q]:oj[q+1]] = B(mo, mi)[cols[mo]] @ QT[mi]
        return X
    perto = lambda I, J: min(abs(x - y) for x in jan[I] for y in jan[J]) <= nk.BAND_M
    N0 = np.zeros((nj, nj)); N = np.zeros((nj, nj))
    for I in range(nj):                       # um V_I por vez (memoria)
        A = bloco(I, I)
        if not A.size:
            continue
        N0[I, I] = norma2(A)
        VI = None
        if N0[I, I] > a.limiar:
            VI = np.linalg.inv(np.eye(len(A)) + A)
            ra = float(max(abs(np.linalg.eigvals(A))))
            print(f'  janela {I} (modos {jan[I][0]}..{jan[I][-1]}, {len(A)}): ||A_II|| {N0[I, I]:.4f}, raio espectral'
                  f' {ra:.4f}, ||V_I|| {norma2(VI):.4f}  [{time.time()-t0:.0f}s]', flush=True)
        else:
            N[I, I] = N0[I, I]
        del A
        for J in range(nj):
            if J == I or not perto(I, J):
                continue
            A = bloco(I, J)
            if not A.size:
                continue
            N0[I, J] = norma2(A)
            N[I, J] = norma2(VI @ A) if VI is not None else N0[I, J]
            del A
        del VI
    np.set_printoptions(precision=3, suppress=True, linewidth=160)
    raio = lambda M: float(max(abs(np.linalg.eigvals(M))))
    print(f'  sem precondicionar: ||N0||_2 {np.linalg.norm(N0, 2):.4f}, raio {raio(N0):.4f}')
    print(f'  bloco-Jacobi nas janelas com ||A_JJ|| > {a.limiar}: raio de Perron {raio(N):.4f}, ||N||_2 {np.linalg.norm(N, 2):.4f}'
          f'  [{time.time()-t0:.0f}s]')
    print(N)


if __name__ == '__main__':
    main()
