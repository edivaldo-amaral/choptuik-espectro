#!/usr/bin/env python3
"""S4: o vetor de gauge g = (Z + 1) omega na base de L (Fourier com sinal, Chebyshev simetrico),
Z = mu^-1 d_tau + xi d_xi (GAUGE_RT_ACTION.md). Para qualquer fundo w,

    L_w(mu) (Z + 1) w = (Z + 2) F_w,

com F_w o defeito de w nas equacoes selecionadas: para o fundo exato, g e autovetor de L em s = mu.

Convencoes: modo de Fourier m com d_tau -> i m/2; xi d_xi nos coeficientes simetricos
(signed_constraints.xidxi: n na diagonal e 2n nas linhas n - 2, n - 4, ...). Os quatro campos de
RefA (signed_operator_L.campos) na ordem das componentes de L (0: n impar, m par; 1, 2: m par;
3: m impar).

Teste (--testar): monta L(mu) numa caixa grande (sem pesos) e mede ||L(mu) g|| relativo a ||g||
na parte que a caixa enxerga inteira; deve ser pequeno (o defeito de RefA e o truncamento), e o
cosseno com o autovetor numerico da truncagem perto de mu deve ser ~1."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import signed_operator_L as sl
from signed_constraints import xidxi

MU = sl.MU


def g_modo(g, m, N, mu=MU):
    """Coeficientes (sem pesos) de (Z + 1) omega no modo m, no layout Radial(m, N)."""
    r = sl.Radial(m, N)
    out = np.zeros(r.dim, complex)
    for c, ns in r.comps:
        a = sl.modo(g, c, m)                    # coeficientes simetricos n = 0..100 do campo c no modo m
        full = np.zeros(max(N, len(a)), complex); full[:len(a)] = a
        v = full[ns] if ns.max() < len(full) else np.array([full[n] if n < len(full) else 0 for n in ns])
        # (Z + 1) = 1 + i m/(2 mu) + xi d_xi
        X = xidxi(ns, ns)
        o0, o1 = r.off[c]
        out[o0:o1] = (1 + 1j*m/(2*mu))*v + X @ v
    return out, r


def vetor(M, N, mu=MU, g=None):
    g = sl.campos() if g is None else g
    partes = []; layout = []
    for m in sl.modos(M):
        v, r = g_modo(g, m, N, mu)
        partes.append(v); layout.append(r)
    return np.concatenate(partes), layout


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--caixa', default='12x48', help='caixa do teste (M x N)')
    a = ap.parse_args()
    M, N = map(int, a.caixa.split('x'))
    g = sl.campos()
    h, lay = vetor(M, N, MU, g)
    L, _, _ = sl.operador(M, N, s=MU, g=g)
    r = L @ h
    # a parte que a caixa enxerga inteira: linhas longe das bordas (|m| < M - 40 nao existe em caixa
    # pequena; use o total e o interior radial n < N - 101)
    off = np.cumsum([0] + [x.dim for x in lay])
    inter = []
    for i, x in enumerate(lay):
        nn = np.concatenate([ns for _, ns in x.comps])
        inter.append(off[i] + np.where(nn < N - 110)[0] if N > 110 else off[i] + np.arange(x.dim))
    inter = np.concatenate(inter)
    print(f'caixa {M}x{N}: ||g|| = {np.linalg.norm(h):.4e}; ||L(mu) g|| / ||g|| = {np.linalg.norm(r)/np.linalg.norm(h):.3e}'
          f' (interior radial: {np.linalg.norm(r[inter])/np.linalg.norm(h):.3e})')
    w, V = np.linalg.eig(L + MU*0)
    i0 = np.argmin(np.abs(w))                    # L(mu) singular: autovalor de L(mu) perto de 0
    v = V[:, i0]
    cos = abs(np.vdot(v, h))/(np.linalg.norm(v)*np.linalg.norm(h))
    print(f'autovalor de L(mu) mais perto de 0: {w[i0]:.3e} (s = mu - isso: {MU - w[i0]:.10f}); |cos(g, autovetor)| = {cos:.8f}')


if __name__ == '__main__':
    main()


def residuo_completo(Mg=50, Ng=110, mu=MU, g=None):
    """||L(mu) g|| sem truncar a saida: g nos modos |m| < Mg, n < Ng (o suporte de RefA cabe), e as
    linhas em todos os modos que a banda alcanca (|m| < Mg + 40) e n < Ng + 101. Sem pesos e com os
    pesos de L."""
    import nk_L as nk
    g = sl.campos() if g is None else g
    gm = {m: g_modo(g, m, Ng, mu)[0] for m in sl.modos(Mg)}
    No = Ng + sl.BAND_N + 1
    tot2 = 0.0; totw2 = 0.0; gn2 = 0.0; gnw2 = 0.0
    for m, v in gm.items():
        r = sl.Radial(m, Ng); w = nk.pesos(r)
        gn2 += float(np.sum(np.abs(v)**2)); gnw2 += float(np.sum(np.abs(w*v)**2))
    for mo in range(-(Mg + sl.BAND_M - 1), Mg + sl.BAND_M):
        ro = sl.Radial(mo, No)
        acc = np.zeros(ro.dim, complex)
        for mi, v in gm.items():
            if abs(mo - mi) <= sl.BAND_M:
                acc += sl.bloco_B(g, ro, sl.Radial(mi, Ng), mu) @ v
        if mo in gm:
            ri = sl.Radial(mo, Ng)
            J = sl.J_livre(ri, mu, mu)                     # J(mu) + mu: s = mu
            nn_o = np.concatenate([ns for _, ns in ro.comps]); pos = {}
            o = 0
            for c, ns in ro.comps:
                for k, n in enumerate(ns):
                    pos[(c, int(n))] = o + k
                o += len(ns)
            idx = np.array([pos[(c, int(n))] for c, ns in ri.comps for n in ns])
            acc[idx] += J @ gm[mo]
        w = nk.pesos(ro)
        tot2 += float(np.sum(np.abs(acc)**2)); totw2 += float(np.sum(np.abs(w*acc)**2))
    return np.sqrt(tot2/gn2), np.sqrt(totw2/gnw2)


def salvar_Z(MZ, NZ, caminho, mu=MU, Mg=50, Ng=110):
    """Grava o vetor de gauge pesado (pesos de L) restrito a Z = MZ x NZ, normalizado, e devolve
    ||P_T g_w||/||P_Z g_w|| (cauda fora de Z, com g no suporte de RefA: |m| < Mg, n < Ng)."""
    import nk_L as nk
    g = sl.campos()
    partes = []; cauda2 = 0.0
    for m in sl.modos(max(MZ, Mg)):
        N = max(NZ, Ng)
        v, r = g_modo(g, m, N, mu)
        vw = nk.pesos(r)*v
        nn = np.concatenate([ns for _, ns in r.comps])
        if abs(m) < MZ:
            partes.append(vw[nn < NZ]); cauda2 += float(np.sum(np.abs(vw[nn >= NZ])**2))
        else:
            cauda2 += float(np.sum(np.abs(vw)**2))
    hz = np.concatenate(partes)
    nz = float(np.linalg.norm(hz))
    np.save(caminho, hz/nz)
    return np.sqrt(cauda2)/nz, nz
