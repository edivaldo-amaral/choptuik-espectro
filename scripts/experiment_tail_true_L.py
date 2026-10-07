#!/usr/bin/env python3
"""DIAGNOSTICO EM PONTO FLUTUANTE -- NAO E PROVA. Norma VERDADEIRA do nivel de cauda de L em
F = 200x400, sem e com bloco-Jacobi (os niveis D de nk_L2), por potencia com produtos implicitos
sobre toda a cauda T = F \\ Z. Separa perda de agregacao (normas de blocos) de acoplamento
intrinseco. Q0 por modo em O(n): J e diagonal mais triangular superior de colunas constantes,
(J x)_j = d_j x_j + sum_{b > j} 2 mu n_b x_b; a inversa sai por substituicao com soma acumulada."""
import argparse, math, sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import closed_form as cf
import nk_L as nk
import nk_L2
import signed_operator_L as sl


class Q0op:
    """Q0 pesado (sqrt omega) do modo m truncado em n < N, aplicado em O(n) por componente."""
    def __init__(self, m, N, lam0):
        r = sl.Radial(m, N)
        self.partes = []
        for c, ns in r.comps:
            o0, o1 = r.off[c]
            d = sl.MU*(ns + 1.0) + lam0 + 1j*m/2
            v = 2*sl.MU*ns
            w = np.sqrt(cf.omega(ns, sl.K2))
            self.partes.append((o0, o1, d, v, w))
        self.dim = r.dim

    def __call__(self, y):
        out = np.empty_like(y)
        for o0, o1, d, v, w in self.partes:
            b = y[o0:o1]/w
            x = np.empty_like(b); s = 0.0
            for k in range(len(b) - 1, -1, -1):
                x[k] = (b[k] - s)/d[k]; s += v[k]*x[k]
            out[o0:o1] = w*x
        return out

    def adj(self, y):
        out = np.empty_like(y)
        for o0, o1, d, v, w in self.partes:
            # (W J^-1 W^-1)^* = W^-1 J^-* W: J^* = diag(conj d) + inferior de linhas constantes
            b = y[o0:o1]*w
            x = np.empty_like(b); s = 0.0
            dc = np.conj(d)
            for k in range(len(b)):
                x[k] = (b[k] - v[k]*s)/dc[k]; s += x[k]
            out[o0:o1] = x/w
        return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--Z', default='32x128'); ap.add_argument('--F', default='200x400')
    ap.add_argument('--W', type=int, default=16); ap.add_argument('--densas', default='16,80')
    ap.add_argument('--N1', type=int, default=200); ap.add_argument('--lam0', type=float, default=0.4010247330)
    ap.add_argument('--iters', type=int, default=40); ap.add_argument('--dados-Z', default='0,0,0,0,0')
    a = ap.parse_args()
    cert = nk_L2.Certificado(a)
    ctx, Nc = cert.ctx, cert.Nc
    t0 = time.time()
    # Q0 por modo, O(n)
    ms = list(range(-(cert.Mc - 1), cert.Mc))
    Qs = {m: Q0op(m, Nc, a.lam0) for m in ms}
    # teste de Q0op contra o inverso denso num modo
    m0 = 5; Qd = ctx.Q(m0, Nc); y = np.random.default_rng(0).standard_normal(Qd.shape[0]) + 0j
    print(f'Q0op x inverso denso (m = {m0}): {np.abs(Qs[m0](y) - Qd @ y).max():.1e}; adjunto {np.abs(Qs[m0].adj(y) - Qd.conj().T @ y).max():.1e}', flush=True)
    ctx.Q.cache_clear()
    # T por modo: uniao dos niveis
    selT = {m: nk.colsT(ctx, m, cert.MZ, cert.NZ, Nc) for m in ms}
    oT = np.cumsum([0] + [len(selT[m]) for m in ms]); nT = int(oT[-1]); pos = {m: q for q, m in enumerate(ms)}
    Bc = {}
    def B(mo, mi):
        k = (mo - mi, mi % 2)
        if k not in Bc:
            Bc[k] = ctx.B(mo, mi, Nc, Nc); ctx._B.cache_clear()
        return Bc[k]
    def A(x):          # P_T B Q0 P_T (banda 16)
        u = {}
        for m in ms:
            y = np.zeros(Qs[m].dim, complex); y[selT[m]] = x[oT[pos[m]]:oT[pos[m]+1]]
            u[m] = Qs[m](y)
        out = np.zeros(nT, complex)
        for mo in ms:
            acc = np.zeros(Qs[mo].dim, complex)
            for mi in range(max(-(cert.Mc - 1), mo - nk_L2.BANDA), min(cert.Mc, mo + nk_L2.BANDA + 1)):
                acc += B(mo, mi) @ u[mi]
            out[oT[pos[mo]]:oT[pos[mo]+1]] = acc[selT[mo]]
        return out
    def At(w):
        t = {m: np.zeros(Qs[m].dim, complex) for m in ms}
        for mo in ms:
            wm = np.zeros(Qs[mo].dim, complex); wm[selT[mo]] = w[oT[pos[mo]]:oT[pos[mo]+1]]
            for mi in range(max(-(cert.Mc - 1), mo - nk_L2.BANDA), min(cert.Mc, mo + nk_L2.BANDA + 1)):
                t[mi] += B(mo, mi).conj().T @ wm
        out = np.zeros(nT, complex)
        for m in ms:
            out[oT[pos[m]]:oT[pos[m]+1]] = Qs[m].adj(t[m])[selT[m]]
        return out
    alfa = nk._potencia(A, At, nT, iters=a.iters)
    print(f'|T| = {nT}: ||P_T B Q0 P_T|| verdadeiro {alfa:.4f}  [{time.time()-t0:.0f}s]', flush=True)
    # bloco-Jacobi: V_D nas linhas densas
    cert.densos()
    idxD = []
    for lv in cert.niveis:
        if lv['tipo'] == 'D':
            ii = np.concatenate([oT[pos[m]] + np.searchsorted(selT[m], lv['sel'][m]) for m in lv['modos'] if len(lv['sel'][m])])
            idxD.append((lv['V'], ii))
    # E = I - A_T (I + A) = -(A_T A + A_T - I); com A_T = diag(V_D, I): nas linhas D, E = I - V (I + A) = -V A_off
    # (V inverte I + A_DD). Implementacao: y = A x; nas linhas D: V (A x)_D - V A_DD x_D ... mais simples:
    def Ecomp(x):
        y = A(x)
        for V, ii in idxD:
            yD = y[ii] + x[ii]              # ((I + A) x)_D
            y[ii] = V @ yD - x[ii]          # V (I + A) x - x = -(E x)_D
        return y
    def Ecomp_t(w):
        w2 = w.copy()
        extra = np.zeros(nT, complex)
        for V, ii in idxD:
            u = V.conj().T @ w[ii]
            w2[ii] = u
            extra[ii] += u - w[ii]
        return At(w2) + extra
    th = nk._potencia(Ecomp, Ecomp_t, nT, iters=a.iters)
    print(f'||E_TT|| verdadeiro com bloco-Jacobi (8 janelas densas em n < {a.N1}): {th:.4f}  [{time.time()-t0:.0f}s]', flush=True)


if __name__ == '__main__':
    main()
