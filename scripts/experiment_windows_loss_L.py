#!/usr/bin/env python3
"""DIAGNOSTICO EM PONTO FLUTUANTE -- NAO E PROVA. Quanto a cota por janelas de
||P_T B Q0 P_T|| (L, S3b) perde contra a norma verdadeira, no Z do certificado.

Em Z = 32x128, F = 200x400, as janelas de W = 16 deram alfa_T = 0,952 (theta0 1,143). Aqui:
alfa_T verdadeiro por iteracao de potencia com produtos IMPLICITOS sobre todo T, e o das
janelas, num F menor (a borda de Z, onde mora o maior bloco, fica igual). Tambem o raio de
Perron da matriz de normas das janelas (se as janelas virarem niveis de Perron)."""
import argparse, sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import nk_L as nk


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--Z', default='32x128'); ap.add_argument('--F', default='80x200')
    ap.add_argument('--W', default='16', help='larguras de janela, separadas por virgula')
    ap.add_argument('--lam0', type=float, default=0.4010247330)
    ap.add_argument('--iters', type=int, default=50)
    ap.add_argument('--so-vetor', action='store_true')
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
            Bs[k] = ctx.B(mo, mi, Nc, Nc)
            ctx._B.cache_clear()
        return Bs[k]
    oj = np.cumsum([0] + [len(cols[m]) for m in ms]); n = int(oj[-1])
    viz = {mo: [mi for mi in ms if abs(mo - mi) <= nk.BAND_M] for mo in ms}
    def mv(v):
        u = {m: QT[m] @ v[oj[q]:oj[q+1]] for q, m in enumerate(ms)}
        out = np.zeros(n, complex)
        for p, mo in enumerate(ms):
            acc = np.zeros(len(cols[mo]), complex)
            for mi in viz[mo]:
                acc += B(mo, mi)[cols[mo]] @ u[mi]
            out[oj[p]:oj[p+1]] = acc
        return out
    def rmv(w):
        t = {m: np.zeros(QT[m].shape[0], complex) for m in ms}
        for p, mo in enumerate(ms):
            wm = w[oj[p]:oj[p+1]]
            for mi in viz[mo]:
                t[mi] += B(mo, mi)[cols[mo]].conj().T @ wm
        return np.concatenate([QT[m].conj().T @ t[m] for m in ms])
    rng = np.random.default_rng(0)
    v = rng.standard_normal(n) + 1j*rng.standard_normal(n); v /= np.linalg.norm(v)
    for _ in range(a.iters):
        u = rmv(mv(v)); v = u/np.linalg.norm(u)
    Av = mv(v); alfa = float(np.linalg.norm(Av))
    print(f'Z = {MZ}x{NZ}, F = {Mc}x{Nc} (|T| = {n}): alfa_T verdadeiro {alfa:.4f}  [{time.time()-t0:.0f}s]', flush=True)
    # onde vivem os vetores singulares dominantes: massa por faixa de |m| e de n
    def massa(x, nome):
        fm = [(0, 16), (16, 32), (32, 40), (40, 48), (48, 56), (56, 64), (64, 999)]
        fn = [(0, 64), (64, 128), (128, 160), (160, 200), (200, 9999)]
        tab = np.zeros((len(fm), len(fn)))
        for q, m in enumerate(ms):
            r = ctx.rad(m, Nc); nn = np.concatenate([ns for _, ns in r.comps])[cols[m]]
            xm = np.abs(x[oj[q]:oj[q+1]])**2
            i = next(k for k, (lo, hi) in enumerate(fm) if lo <= abs(m) < hi)
            for k, (lo, hi) in enumerate(fn):
                tab[i, k] += xm[(nn >= lo) & (nn < hi)].sum()
        tab /= tab.sum()
        print(f'  massa de {nome} (linhas |m| em {fm}; colunas n em {fn}):')
        print(np.array2string(tab, precision=3, suppress_small=True))
    massa(v, 'v (entrada)'); massa(Av/alfa, 'u (saida)')
    if a.so_vetor:
        return
    for W in map(int, a.W.split(',')):
        aT, afT, QTn, Nw = nk.janelas(ctx, MZ, NZ, Mc, Nc, W, log=lambda *x, **k: None)
        rP = float(max(abs(np.linalg.eigvals(Nw))))
        print(f'  W = {W}: janelas ||N_win||_2 {aT:.4f} (perda {aT/alfa:.2f}x), raio de Perron {rP:.4f}'
              f' ({rP/alfa:.2f}x), maior bloco {Nw.max():.4f}  [{time.time()-t0:.0f}s]', flush=True)
        np.set_printoptions(precision=3, suppress=True, linewidth=160)
        print(Nw)


if __name__ == '__main__':
    main()
