#!/usr/bin/env python3
"""Protótipo em ponto flutuante do nível Z do Rouché de L no contorno (C1_REAVALIACAO.md §7).

Confere, contra o cálculo direto (experiment_tile_variacao), as identidades do desenho:
  (1) V_Z(s) Y(s) = [V_Z(s) Y~ - J_ZT] Q0_TT(s),   Y~ = P_Z L~(0) P_T,  J_ZT = P_Z J P_T;
  (2) a escala D_Z = J_Z(c) Q_ZZ(s): entrada Z <- T = ||J_Z(c)[R(s) Y~ - Q_ZZ(s) J_ZT] Q0_TT(s)||,
      que em s = c coincide com a(c);
  (3) o posto baixo: Gram de [[Y~ Q0_TT(c)], [J_ZT Q0_TT(c)]] = F F^* + resto, e a entrada via F com
      R(s) aplicado pela decomposição de Schur de L~_Z (A = U T U^*).
Tudo pesado nos pesos de L; L~ = L + deflação de Brauer (nk_L.nivel_Z)."""
from __future__ import annotations

import argparse
import math
from pathlib import Path
import sys
import time

import numpy as np
from scipy.linalg import schur, solve_triangular

sys.path.insert(0, str(Path(__file__).resolve().parent))
import nk_L as nk
import signed_operator_L as sl
from experiment_tile_variacao import gram_Y


def montar(MZ, NZ, Mc, Nc, c, dfl):
    """Contexto em c, nível Z (V_Z(c), deflação) e as peças constantes: A = W L~_Z(0) W^-1,
    blocos de Y~ e J_ZT por modo de T, Q0_TT(c) por modo."""
    ctx = nk.Contexto(c, Nc)
    z = nk.nivel_Z(ctx, MZ, NZ, log=lambda *x, **k: None, sem_borda=True, deflacao=dfl)
    ms, offs, rs = z['ms'], z['offs'], z['rs']; n = z['n']; idx = {m: i for i, m in enumerate(ms)}
    # A = W L~_Z(0) W^-1 (pesado): J(0) + B + deflação
    A = np.zeros((n, n), complex)
    for j, ri in enumerate(rs):
        w = nk.pesos(ri)
        A[offs[j]:offs[j+1], offs[j]:offs[j+1]] += (w[:, None]*sl.J_livre(ri, 0.0))/w[None, :]
        for i, ro in enumerate(rs):
            if abs(ro.m - ri.m) <= nk.BAND_M:
                A[offs[i]:offs[i+1], offs[j]:offs[j+1]] += ctx.B(ro.m, ri.m, NZ, NZ)
    for dk, vw, uw in z['defl']:
        A += dk*np.outer(vw, uw.conj())
    # colunas de T por modo: Y~_m = P_Z (J + B) P_{T,m} (E P_T = 0: phi em Z), J_ZT,m (só o próprio modo), Q0_TT,m(c)
    blocos = []
    for m in range(-(min(MZ + nk.BAND_M, Mc) - 1), min(MZ + nk.BAND_M, Mc)):
        cT = nk.colsT(ctx, m, MZ, NZ, Nc)
        Yt = np.zeros((n, len(cT)), complex)
        for mo in ms:
            if abs(mo - m) <= nk.BAND_M:
                i = idx[mo]
                Yt[offs[i]:offs[i+1]] = ctx.B(mo, m, NZ, Nc)[:, cT]
        Jzt = np.zeros((n, len(cT)), complex)
        if m in idx:
            i = idx[m]; r = ctx.rad(m, Nc); w = nk.pesos(r)
            Jw = (w[:, None]*sl.J_livre(r, 0.0))/w[None, :]
            Jzt[offs[i]:offs[i+1]] = Jw[np.ix_(nk.nZ_local(ctx, m, NZ, Nc), cT)]
        Yt += Jzt                                   # Y~ = P_Z (J + B + E) P_T inclui J_ZT
        QTT = ctx.Q(m, Nc)[np.ix_(cT, cT)]
        blocos.append((m, cT, Yt, Jzt, QTT))
    return ctx, z, A, blocos


def QZZ_de(s, z, NZ):
    ctx = nk.Contexto(s, NZ)
    n = z['n']; out = np.zeros((n, n), complex)
    for j, r in enumerate(z['rs']):
        out[z['offs'][j]:z['offs'][j+1], z['offs'][j]:z['offs'][j+1]] = ctx.Q(r.m, NZ)
    return out


def JZ_de(s, z):
    n = z['n']; out = np.zeros((n, n), complex)
    for j, r in enumerate(z['rs']):
        w = nk.pesos(r)
        out[z['offs'][j]:z['offs'][j+1], z['offs'][j]:z['offs'][j+1]] = (w[:, None]*sl.J_livre(r, s))/w[None, :]
    return out


def norma_bloco(Mleft, blocos, Qs=None):
    """||sum_m Mleft(m) Q0_TT,m|| por Gram: Mleft(m) devolve a matriz n x |T_m|."""
    G = None
    for bl in blocos:
        X = Mleft(bl) @ (bl[4] if Qs is None else Qs[bl[0]])
        G = X @ X.conj().T if G is None else G + X @ X.conj().T
    return math.sqrt(max(np.linalg.eigvalsh(G)[-1], 0.0))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--Z', default='12x48'); ap.add_argument('--F', default='40x160')
    ap.add_argument('--c', type=complex, default=0j)
    ap.add_argument('--eps', default='0.002,0.01j')
    ap.add_argument('--k', type=int, default=400)
    ap.add_argument('--deflacao', default='0,0.16831')
    a = ap.parse_args(); t0 = time.time()
    MZ, NZ = map(int, a.Z.split('x')); Mc, Nc = map(int, a.F.split('x'))
    dfl = [complex(x) for x in a.deflacao.split(',') if x.strip()]
    ctx, z, A, blocos = montar(MZ, NZ, Mc, Nc, a.c, dfl)
    n = z['n']; I = np.eye(n)
    Jc = JZ_de(a.c, z); Qc = QZZ_de(a.c, z, NZ)
    Vc = z['Vb']
    print(f'n = {n}, |T acoplado| = {sum(len(b[1]) for b in blocos)}   [{time.time() - t0:.0f}s]', flush=True)
    # (1) identidade em s = c: direto (gram_Y) contra [V Y~ - J_ZT] Q0_TT
    a_dir = math.sqrt(max(np.linalg.eigvalsh(Vc @ gram_Y(ctx, z, MZ, NZ, Mc, Nc) @ Vc.conj().T)[-1], 0.0))
    a_id = norma_bloco(lambda b: Vc @ b[2] - b[3], blocos)
    print(f'(1) a(c) direto {a_dir:.6f}  via [V Y~ - J_ZT] Q0_TT {a_id:.6f}', flush=True)
    # Schur de A (solver)
    T, U = schur(A, output='complex', check_finite=False)
    def Rsolve(s, X):                     # (A + s)^-1 X via Schur
        return U @ solve_triangular(T + s*I, U.conj().T @ X, check_finite=False)
    # (3) posto baixo do Gram empilhado em c
    G2 = np.zeros((2*n, 2*n), complex)
    for (m, cT, Yt, Jzt, QTT) in blocos:
        S = np.vstack([Yt @ QTT, Jzt @ QTT]); G2 += S @ S.conj().T
    lam, W = np.linalg.eigh(G2); lam = lam[::-1]; W = W[:, ::-1]
    k = min(a.k, 2*n); F = W[:, :k]*np.sqrt(np.maximum(lam[:k], 0)); resto = max(lam[k] if k < 2*n else 0.0, 0.0)
    print(f'(3) posto {k}: resto do Gram {resto:.2e} (sqrt {math.sqrt(resto):.2e})   [{time.time() - t0:.0f}s]', flush=True)
    Ft, Fb = F[:n], F[n:]
    for e in [0j] + [complex(x) for x in a.eps.split(',')]:
        s = a.c + e
        Qs_ZZ = QZZ_de(s, z, NZ)
        ctx_s = nk.Contexto(s, Nc)
        Qs = {b[0]: ctx_s.Q(b[0], Nc)[np.ix_(b[1], b[1])] for b in blocos}
        # (2) escala: ||J_Z(c)[R(s) Y~ - Q_ZZ(s) J_ZT] Q0_TT(s)|| (direto, sem posto baixo)
        esc = norma_bloco(lambda b: Jc @ (Rsolve(s, b[2]) - Qs_ZZ @ b[3]), blocos, Qs)
        # via posto baixo (com Q0_TT(c): o fator 1/(1 - |e| q_TT) cobre Q0_TT(s))
        X = Jc @ (Rsolve(s, Ft) - Qs_ZZ @ Fb)
        lr = float(np.linalg.norm(X, 2))
        # a(s) no nível Z não escalado (referência)
        ctx2, z2 = ctx_s, nk.nivel_Z(ctx_s, MZ, NZ, log=lambda *x, **k: None, sem_borda=True, deflacao=dfl)
        a_s = math.sqrt(max(np.linalg.eigvalsh(z2['Vb'] @ gram_Y(ctx2, z2, MZ, NZ, Mc, Nc) @ z2['Vb'].conj().T)[-1], 0.0))
        print(f'  s = c + {e}: Z<-T escalado {esc:.5f}  (posto baixo com Q0_TT(c): {lr:.5f});  a(s) sem escala {a_s:.5f}   [{time.time() - t0:.0f}s]', flush=True)


if __name__ == '__main__':
    main()
