#!/usr/bin/env python3
"""DIAGNOSTICO EM PONTO FLUTUANTE -- NAO E PROVA. Nivel T1 com inverso aproximado DENSO para L.

Em Z = 32x128 o singular dominante de A = P_T B Q0 P_T (0,637) mora em 32 <= |m| < 40,
64 <= n < 160; as janelas perdem 1,44x por somar acoplamentos vizinhos. Aqui T1 = {M1a <= |m|
< M1b, n < N1} (dois lados, desacoplados se 2 M1a > banda) recebe V_1 = (I + A_11)^-1, e o
nivel T da certificacao vira (T1, T2) com

    E_11 = 0,  E_12 = -V_1 A_12,  E_21 = -A_21,  E_22 = -A_22.

Mede ||V_1||, ||V_1 A_12||, ||A_21||, ||A_22|| (verdadeiro, potencia implicita) e por janelas,
e o raio de Perron da parte T."""
import argparse, sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import nk_L as nk
from tile_certificate import norma2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--Z', default='32x128'); ap.add_argument('--F', default='80x200')
    ap.add_argument('--T1', default='32,48,160', help='M1a,M1b,N1')
    ap.add_argument('--W', type=int, default=16)
    ap.add_argument('--lam0', type=float, default=0.4010247330)
    ap.add_argument('--iters', type=int, default=40)
    a = ap.parse_args()
    MZ, NZ = map(int, a.Z.split('x')); Mc, Nc = map(int, a.F.split('x'))
    M1a, M1b, N1 = map(int, a.T1.split(','))
    ctx = nk.Contexto(a.lam0, Nc)
    t0 = time.time()
    ms = list(range(-(Mc - 1), Mc))
    nn = {m: np.concatenate([ns for _, ns in ctx.rad(m, Nc).comps]) for m in ms}
    colsT = {m: nk.colsT(ctx, m, MZ, NZ, Nc) for m in ms}
    em1 = lambda m: M1a <= abs(m) < M1b
    c1 = {m: (colsT[m][nn[m][colsT[m]] < N1] if em1(m) else np.array([], int)) for m in ms}
    c2 = {m: np.setdiff1d(colsT[m], c1[m]) for m in ms}
    Q = {m: ctx.Q(m, Nc) for m in ms}
    ctx.Q.cache_clear()
    Bs = {}
    def B(mo, mi):
        k = (mo - mi, mi % 2)
        if k not in Bs:
            Bs[k] = ctx.B(mo, mi, Nc, Nc); ctx._B.cache_clear()
        return Bs[k]
    viz = {mo: [mi for mi in ms if abs(mo - mi) <= nk.BAND_M] for mo in ms}

    def operador(lin, col, modos_lin=None, modos_col=None):
        """mv, rmv, n_col, n_lin de P_lin B Q0 P_col (produtos implicitos modo a modo)."""
        ml = [m for m in (modos_lin or ms) if len(lin[m])]
        mc = [m for m in (modos_col or ms) if len(col[m])]
        oc = np.cumsum([0] + [len(col[m]) for m in mc]); ol = np.cumsum([0] + [len(lin[m]) for m in ml])
        QC = {m: Q[m][:, col[m]] for m in mc}
        setc = set(mc)
        def mv(v):
            u = {m: QC[m] @ v[oc[q]:oc[q+1]] for q, m in enumerate(mc)}
            out = np.zeros(int(ol[-1]), complex)
            for p, mo in enumerate(ml):
                acc = np.zeros(len(lin[mo]), complex)
                for mi in viz[mo]:
                    if mi in setc:
                        acc += B(mo, mi)[lin[mo]] @ u[mi]
                out[ol[p]:ol[p+1]] = acc
            return out
        def rmv(w):
            t = {m: np.zeros(Q[m].shape[0], complex) for m in mc}
            for p, mo in enumerate(ml):
                wm = w[ol[p]:ol[p+1]]
                for mi in viz[mo]:
                    if mi in setc:
                        t[mi] += B(mo, mi)[lin[mo]].conj().T @ wm
            return np.concatenate([QC[m].conj().T @ t[m] for m in mc])
        return mv, rmv, int(oc[-1]), int(ol[-1]), ml, mc, ol, oc

    # A_11 denso, por lado
    res = {}
    V1 = {}
    for lado, sgn in (('+', 1), ('-', -1)):
        m1 = [m for m in ms if em1(m) and np.sign(m) == sgn]
        sel = {m: c1[m] if m in m1 else np.array([], int) for m in ms}
        mv, rmv, nc, nl, ml, mc, ol, oc = operador(sel, sel, m1, m1)
        A11 = np.column_stack([mv(e) for e in np.eye(nc)])
        H = np.eye(nc) + A11
        V = np.linalg.inv(H)
        V1[lado] = (V, sel, m1)
        res[lado] = dict(n=nc, A11=norma2(A11), rA11=float(max(abs(np.linalg.eigvals(A11)))), V1=norma2(V))
        print(f'  T1{lado} ({nc}): ||A_11|| {res[lado]["A11"]:.4f}, raio espectral {res[lado]["rA11"]:.4f},'
              f' ||V_1|| {res[lado]["V1"]:.4f}  [{time.time()-t0:.0f}s]', flush=True)
    # V_1 A_12 (linhas T1, colunas T2) e A_21, conjuntos nos dois lados
    mv12, rmv12, n2, n1, ml12, _, ol12, _ = operador(c1, c2)
    # V_1 bloco-diagonal nos dois lados, na ordem das linhas de T1 (modos crescentes)
    pos = {m: (ol12[p], ol12[p+1]) for p, m in enumerate(ml12)}
    blocos = []
    for lado in ('+', '-'):
        V, sel, m1 = V1[lado]
        blocos.append((V, np.concatenate([np.arange(*pos[m]) for m in m1])))
    def aplica(x, adj=False):
        y = np.zeros_like(x)
        for V, idx in blocos:
            y[idx] = (V.conj().T if adj else V) @ x[idx]
        return y
    a12 = nk._potencia(lambda v: aplica(mv12(v)), lambda w: rmv12(aplica(w, True)), n2, iters=a.iters)
    b12 = nk._potencia(mv12, rmv12, n2, iters=a.iters)
    print(f'  T1 <- T2: ||A_12|| {b12:.4f}, ||V_1 A_12|| {a12:.4f}  [{time.time()-t0:.0f}s]', flush=True)
    mv21, rmv21, n1c, _, *_ = operador(c2, c1)
    a21 = nk._potencia(mv21, rmv21, n1c, iters=a.iters)
    print(f'  T2 <- T1: ||A_21|| {a21:.4f}  [{time.time()-t0:.0f}s]', flush=True)
    mv22, rmv22, n22, *_ = operador(c2, c2)
    a22 = nk._potencia(mv22, rmv22, n22, iters=a.iters)
    print(f'  T2 <- T2: ||A_22|| verdadeiro {a22:.4f}  [{time.time()-t0:.0f}s]', flush=True)
    # janelas em T2
    jan = [ms[k:k+a.W] for k in range(0, len(ms), a.W)]
    Nw = np.zeros((len(jan), len(jan)))
    for J, mJ in enumerate(jan):
        for I, mI in enumerate(jan):
            if min(abs(x - y) for x in mI for y in mJ) > nk.BAND_M:
                continue
            mv, rmv, ncJ, nlI, *_ = operador(c2, c2, mI, mJ)
            if ncJ and nlI:
                Nw[I, J] = nk._potencia(mv, rmv, ncJ, iters=a.iters)
    a22w = float(np.linalg.norm(Nw, 2))
    print(f'  T2 <- T2 por janelas (W = {a.W}): {a22w:.4f} (perda {a22w/a22:.2f}x)  [{time.time()-t0:.0f}s]', flush=True)
    np.set_printoptions(precision=3, suppress=True, linewidth=160)
    print(Nw)
    for nome, x22 in (('verdadeiro', a22), ('janelas', a22w)):
        N = np.array([[0.0, a12], [a21, x22]])
        print(f'  Perron da parte T (T1, T2) com A_22 {nome}: {max(abs(np.linalg.eigvals(N))):.4f}')


if __name__ == '__main__':
    main()
