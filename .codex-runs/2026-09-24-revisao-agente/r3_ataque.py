#!/usr/bin/env python3
"""R3 -- ataque as cotas N[i][j] de tile_perron.montar com montagem PROPRIA de E.

Caixa pequena: Z = 8x24, G+ = 12x36, F = 32x77, D = 4 (banda 20,40 como no certificado).
Modo rigoroso do tile_perron (RIG = True) para as cotas; a verdade e montada aqui:

  * B completo (campos sem truncar, todos os d da caixa) via ft.B_rad -- o mesmo operador que
    signed_operator.operador (conferido contra ele: K(0) - blocos de J);
  * Q0(sigma) por modo pela forma fechada do R1 (nao usa livre_inv);
  * pesos kappa1^m sqrt(omega(n)); V e o do certificado (parte_Z), reordenado para a caixa.

Para linhas e colunas dentro de F, P_i E P_j do operador infinito e o da caixa (Q0 preserva F).
E = I - A H, H - I = sum_l (B + delta_l) Q_l P_l, delta_l = s - s_l. Cada bloco (i, j) e afim num
unico delta_j:  E_ij = E0_ij + delta_j C_ij. Testes:
  (a) ||E0_ij|| contra N0[i][j] e ||C_ij|| contra N_j[i][j] (as duas pecas da cota);
  (b) max sobre s na borda do tile (32 direcoes) de ||E_ij(s)|| contra N0 + r_j N_j;
  (c) max sobre |delta_j| = r_j (32 direcoes; mais forte que o tile) contra N0 + r_j N_j.
Niveis far: blocos (far, Z) e (far, T1) EXATOS numa caixa maior (a banda completa do fundo
leva Z/G+ no maximo a |m| + 40, n + 100); (Z, far) e (T1, far) PARCIAIS (colunas far dentro de
uma caixa maior: cota inferior da verdade).
"""
from __future__ import annotations

import json
import math
from pathlib import Path
import sys
import time

import numpy as np

AQUI = Path(__file__).resolve().parent
RAIZ = AQUI.parents[1]
sys.path.insert(0, str(RAIZ/'scripts')); sys.path.insert(0, str(AQUI))
import fourier_tail_bound as ft
import tile_perron as tp
import signed_operator as so
from r1_forma_fechada import forma_fechada

MU, K1, K2 = ft.MU, ft.K1, ft.K2
ZB, GP, FB, DD = (8, 24), (12, 36), (32, 77), 4
NIV = ['Z', 'T1', 'T2', 'far']


def norma(X):
    if X.size == 0:
        return 0.0
    if min(X.shape) <= 700:
        return float(np.linalg.norm(X, 2))
    from scipy.sparse.linalg import svds
    return float(svds(X, k=1, return_singular_vectors=False, tol=1e-10)[0])


class Caixa:
    """Indices, pesos e B completo numa caixa (M, N) com modos m em modos(M)."""
    def __init__(self, M, N, g):
        self.M, self.N = M, N
        self.ms = so.modos(M)
        self.n1, self.n2 = so.radial(N)
        self.k = len(self.n1) + len(self.n2)
        self.nloc = np.concatenate([self.n1, self.n2])
        self.comp = np.concatenate([np.ones(len(self.n1), int), 2*np.ones(len(self.n2), int)])
        self.dim = len(self.ms)*self.k
        self.m_de = np.repeat(self.ms, self.k)
        self.n_de = np.tile(self.nloc, len(self.ms))
        self.c_de = np.tile(self.comp, len(self.ms))
        lw = np.where(self.n_de == 0, 0.0, self.n_de*math.log(K2) + 0.5*np.log1p(K2**(-4.0*self.n_de)))
        self.lw = self.m_de*math.log(K1) + lw          # log do peso kappa1^m sqrt(omega(n))
        self.g = g

    def idx(self, m, comp, n):
        i = self.ms.index(m)
        loc = (n - 1)//2 if comp == 1 else len(self.n1) + n
        return i*self.k + loc

    def Bdenso(self):
        B = np.zeros((self.dim, self.dim), complex)
        cache = {}
        for j, mi in enumerate(self.ms):
            for i, mo in enumerate(self.ms):
                d = mo - mi
                if d not in cache:
                    cache[d] = ft.B_rad(self.g, d, self.n1, self.n2, self.n1, self.n2, MU)
                B[i*self.k:(i+1)*self.k, j*self.k:(j+1)*self.k] = cache[d]
        return B

    def Qmodo(self, m, s):
        sg = s + 1j*m/2
        R1 = forma_fechada(self.n1, sg); R2 = forma_fechada(self.n2, sg)
        Z = np.zeros((len(self.n1), len(self.n2)))
        return np.block([[R1, Z], [Z.T, R2]])


def pesar(X, lw_lin, lw_col):
    return X*np.exp(lw_lin[:, None] - lw_col[None, :])


def niveis(cx, Z, Gp):
    (MZ, NZ), (MP, NP) = Z, Gp
    am = np.abs(cx.m_de)
    emZ = (am < MZ) & (cx.n_de < NZ)
    emG = (am < MP) & (cx.n_de < NP)
    return {'Z': np.where(emZ)[0], 'T1': np.where(emG & ~emZ)[0], 'T2': np.where(~emG)[0]}


def ordem_V(cx, Z):
    """Posicao na caixa de cada indice de V (ordem de tile_perron.parte_Z)."""
    MZ, NZ = Z
    mZ = tp.Modo(NZ)
    pos = []
    for m in so.modos(MZ):
        for n in mZ.n1:
            pos.append(cx.idx(m, 1, int(n)))
        for n in mZ.n2:
            pos.append(cx.idx(m, 2, int(n)))
    return np.array(pos)


def monta_E(cx, B, V, posV, lv, centros):
    """Blocos E0_ij e C_ij (pesados) para i, j em Z, T1, T2."""
    lw = cx.lw
    # Q_l por modo nas colunas do nivel l; H - I = B Q + delta Q
    BQ = np.zeros((cx.dim, cx.dim), complex)
    Qm = np.zeros((cx.dim, cx.dim), complex)
    for j, m in enumerate(cx.ms):
        sl = slice(j*cx.k, (j+1)*cx.k)
        cols = np.arange(j*cx.k, (j+1)*cx.k)
        for l in ('Z', 'T1', 'T2'):
            cl = np.intersect1d(cols, lv[l])
            if not len(cl):
                continue
            Q = cx.Qmodo(m, centros[l])[:, cl - j*cx.k]
            Qm[sl, cl] = Q
            BQ[:, cl] = B[:, sl] @ Q
    BQ = pesar(BQ, lw, lw); Qm = pesar(Qm, lw, lw)
    # linhas de Z na ordem de V: iZ da caixa deve ser a MESMA lista (conjunto) que posV
    assert set(posV.tolist()) == set(lv['Z'].tolist())
    lvZ = posV                              # linhas/colunas de Z na ordem de V
    ordem = {'Z': lvZ, 'T1': lv['T1'], 'T2': lv['T2']}
    E0, C = {}, {}
    for i in ('Z', 'T1', 'T2'):
        for j in ('Z', 'T1', 'T2'):
            ri, cj = ordem[i], ordem[j]
            HB = BQ[np.ix_(ri, cj)]; HQ = Qm[np.ix_(ri, cj)]
            if i == 'Z':
                Id = np.eye(len(ri)) if j == 'Z' else 0.0
                E0[i, j] = Id - V @ (Id + HB) if j == 'Z' else -(V @ HB)
                C[i, j] = -(V @ HQ)
            else:
                E0[i, j] = -HB
                C[i, j] = -HQ
    del BQ, Qm
    return E0, C


def ataque(cfg, ctx, cx, B, cacheT, log):
    sZ, rZ, s1, r1, s2, r2 = cfg
    assert abs(sZ - s1) + rZ <= r1 and abs(sZ - s2) + rZ <= r2
    silencio = lambda *a, **k: None
    qZ, V = tp.parte_Z(ctx, sZ, silencio)
    if ('T1', s1) not in cacheT:
        cacheT['T1', s1] = tp.parte_T1(ctx, s1, silencio)
    if ('T2', s2) not in cacheT:
        cacheT['T2', s2] = tp.parte_T2(ctx, s2, silencio)
    q1, HZ1, QZ1 = cacheT['T1', s1]
    q2, SZ, SQZ = cacheT['T2', s2]
    N0, NZm, N1m, N2m, qa = tp.montar(ctx, qZ, V, q1, HZ1, QZ1, q2, SZ, SQZ, rZ, r1, r2, silencio)
    Nj = [NZm, N1m, N2m, N2m]
    rj = [rZ, r1, r2, r2]
    lv = niveis(cx, ZB, GP)
    posV = ordem_V(cx, ZB)
    E0, C = monta_E(cx, B, V, posV, lv, {'Z': sZ, 'T1': s1, 'T2': s2})
    cen = {'Z': sZ, 'T1': s1, 'T2': s2}
    res = {}
    fi = np.linspace(0, 2*np.pi, 32, endpoint=False)
    for a, i in enumerate(('Z', 'T1', 'T2')):
        for b, j in enumerate(('Z', 'T1', 'T2')):
            e0 = norma(E0[i, j]); c = norma(C[i, j])
            grande = min(E0[i, j].shape) > 700
            dirs = fi[::4] if grande else fi
            # (b) s na borda do tile
            vt = max(norma(E0[i, j] + (sZ + rZ*np.exp(1j*f) - cen[j])*C[i, j]) for f in dirs) if c > 0 else e0
            # (c) |delta_j| = r_j
            vc = max(norma(E0[i, j] + rj[b]*np.exp(1j*f)*C[i, j]) for f in dirs) if c > 0 else e0
            cota = N0[a][b] + rj[b]*Nj[b][a][b]
            res[i, j] = dict(E0=e0, N0=N0[a][b], C=c, Nj=Nj[b][a][b], tile=vt, circ=vc, cota=cota)
    return res, dict(N0=N0, NZ=NZm, N1=N1m, N2=N2m, V=qZ['V'], VQ_ZZ=qZ['VQ_ZZ']), V, posV


def far_exatos(cfg, ctx, V, cacheT, g0):
    """(far, Z) e (far, T1) exatos: colunas de Z (em sZ) e T1 (em s1); linhas fora de F numa caixa
    que contem todo o alcance do fundo COMPLETO (|d| <= 40, subida radial <= 100)."""
    sZ, rZ, s1, r1, s2, r2 = cfg
    (MZ, NZ), (MP, NP), (Mc, Nc) = ZB, GP, FB
    out = {}
    for nome, (Mcol, Ncol, s, fora) in (('far<-Z', (MZ, NZ, sZ, None)), ('far<-T1', (MP, NP, s1, (MZ, NZ)))):
        cg = Caixa(Mcol + 42, Ncol + 102, g0)       # alcance completo das colunas
        am = np.abs(cg.m_de)
        colsel = (am < Mcol) & (cg.n_de < Ncol)
        if fora is not None:
            colsel &= ~((am < fora[0]) & (cg.n_de < fora[1]))
        linsel = ~((am < Mc) & (cg.n_de < Nc))
        cols = np.where(colsel)[0]; lins = np.where(linsel)[0]
        X = np.zeros((len(lins), len(cols)), complex)
        cacheB = {}
        # B Q nas colunas: por modo de coluna
        for j, mi in enumerate(cg.ms):
            cl = cols[(cols >= j*cg.k) & (cols < (j+1)*cg.k)]
            if not len(cl):
                continue
            Q = cg.Qmodo(mi, s)[:, cl - j*cg.k]
            # Q preserva a caixa da coluna: so linhas n < Ncol importam, mas usamos o modo inteiro
            for i, mo in enumerate(cg.ms):
                d = mo - mi
                li = np.where((lins >= i*cg.k) & (lins < (i+1)*cg.k))[0]
                if not len(li) or abs(d) > 40:
                    continue
                if d not in cacheB:
                    cacheB[d] = ft.B_rad(g0, d, cg.n1, cg.n2, cg.n1, cg.n2, MU)
                Bd = cacheB[d]
                X[li[:, None], np.searchsorted(cols, cl)[None, :]] = (Bd @ Q)[lins[li] - i*cg.k]
        Xw = pesar(X, cg.lw[lins], cg.lw[cols])
        out[nome] = norma(Xw)
    return out


def far_parcial(cfg, V, posV_small, g0, Nfar=260, Mfar=48):
    """(Z, far) e (T1, far) PARCIAIS: colunas far dentro da caixa (Mfar, Nfar); linhas Z (com V) e T1.
    Cota INFERIOR da verdade (restricao a um subespaco de colunas)."""
    sZ, rZ, s1, r1, s2, r2 = cfg
    (MZ, NZ), (MP, NP), (Mc, Nc) = ZB, GP, FB
    cg = Caixa(Mfar, Nfar, g0)
    am = np.abs(cg.m_de)
    colsel = ~((am < Mc) & (cg.n_de < Nc))
    emZ = (am < MZ) & (cg.n_de < NZ); emG = (am < MP) & (cg.n_de < NP)
    lZ = np.array([cg.idx(m, c, n) for (m, c, n) in posV_small])
    l1 = np.where(emG & ~emZ)[0]
    lins = np.concatenate([lZ, l1])
    cols = np.where(colsel)[0]
    HB = np.zeros((len(lins), len(cols)), complex); HQ = np.zeros_like(HB)
    cacheB = {}
    for j, mi in enumerate(cg.ms):
        cl = cols[(cols >= j*cg.k) & (cols < (j+1)*cg.k)]
        if not len(cl):
            continue
        Q = cg.Qmodo(mi, s2)[:, cl - j*cg.k]
        cpos = np.searchsorted(cols, cl)
        for i, mo in enumerate(cg.ms):
            d = mo - mi
            li = np.where((lins >= i*cg.k) & (lins < (i+1)*cg.k))[0]
            if not len(li) or abs(d) > 40:
                continue
            if d not in cacheB:
                cacheB[d] = ft.B_rad(g0, d, cg.n1, cg.n2, cg.n1, cg.n2, MU)
            Bd = cacheB[d]
            HB[li[:, None], cpos[None, :]] = (Bd @ Q)[lins[li] - i*cg.k]
            if d == 0:
                HQ[li[:, None], cpos[None, :]] = Q[lins[li] - i*cg.k]
    HB = pesar(HB, cg.lw[lins], cg.lw[cols]); HQ = pesar(HQ, cg.lw[lins], cg.lw[cols])
    nz = len(lZ)
    EZ0 = -(V @ HB[:nz]); CZ = -(V @ HQ[:nz])
    E10 = -HB[nz:]; C1 = -HQ[nz:]
    fi = np.linspace(0, 2*np.pi, 16, endpoint=False)
    return dict(Zfar_E0=norma(EZ0), Zfar_circ=max(norma(EZ0 + r2*np.exp(1j*f)*CZ) for f in fi),
                T1far_E0=norma(E10), T1far_circ=max(norma(E10 + r2*np.exp(1j*f)*C1) for f in fi))


def main():
    t0 = time.time()
    tp.RIG = True
    ctx = tp.Ctx(ZB, GP, FB, DD)
    g0 = ft.modos_campos()
    cx = Caixa(*FB, g0)
    B = cx.Bdenso()
    # conferencia: B = operador(F, s=0) - blocos de J
    K, _ = so.operador(*FB, s=0.0)
    for j, m in enumerate(cx.ms):
        K[j*cx.k:(j+1)*cx.k, j*cx.k:(j+1)*cx.k] -= so.Jrad(m, 0.0, FB[1])
    dif = float(np.abs(K - B).max()); del K
    print(f'B proprio contra signed_operator.operador - J: max |dif| = {dif:.1e}  (dim {cx.dim})')
    # posicoes de Z na ordem de V, como (m, comp, n), para a caixa maior do teste far
    mZ = tp.Modo(ZB[1])
    posV_small = [(m, 1, int(n)) for m in so.modos(ZB[0]) for n in mZ.n1]
    posV_small = [(m, c, int(n)) for m in so.modos(ZB[0]) for (c, ns) in ((1, mZ.n1), (2, mZ.n2)) for n in ns]
    cfgs = [
        # (sZ, rZ, s1, r1, s2, r2)
        (0.0 + 0.1j, 0.02, 0.02 + 0.1j, 0.05, 0.1 + 0.1j, 0.15),        # Re sZ = 0
        (0.0 + 0.25j, 0.02, 0.02 + 0.23j, 0.05, 0.125 + 0.125j, 0.2),   # Re sZ = 0, Im sZ = 1/4
        (0.0, 0.015, 0.03 + 0.03j, 0.06, 0.125 + 0.125j, 0.2),          # canto s = 0
        (0.25 + 0.25j, 0.03, 0.23 + 0.22j, 0.07, 0.2 + 0.2j, 0.12),     # Im sZ = 1/4
        (0.401 + 0.13j, 0.02, 0.42 + 0.15j, 0.05, 0.375 + 0.125j, 0.1768),  # junto a D
        (0.276, 0.02, 0.276, 0.02, 0.276, 0.02),                        # eixo real, borda de D
        (1.0 + 0.2j, 0.1, 1.05 + 0.15j, 0.18, 1.125 + 0.125j, 0.25),    # Re s grande
    ]
    cacheT = {}
    pior = {}
    saida = []
    for cfg in cfgs:
        res, N, V, posV = ataque(cfg, ctx, cx, B, cacheT, print)
        fx = far_exatos(cfg, ctx, V, cacheT, g0)
        fp = far_parcial(cfg, V, posV_small, g0)
        sZ, rZ, s1, r1, s2, r2 = cfg
        print(f'\ncfg sZ={sZ} rZ={rZ} | s1={s1} r1={r1} | s2={s2} r2={r2}  ||V|| {N["V"]:.3f}'
              f'  [{time.time()-t0:.0f}s]')
        print('   bloco      ||E0||/N0     ||C||/Nj   tile/cota  circ/cota')
        for (i, j), r in res.items():
            rE = r['E0']/r['N0'] if r['N0'] > 0 else (0.0 if r['E0'] == 0 else math.inf)
            rC = r['C']/r['Nj'] if r['Nj'] > 0 else (0.0 if r['C'] < 1e-300 else math.inf)
            rt, rc = r['tile']/r['cota'], r['circ']/r['cota']
            print(f'   {i:>2} <- {j:<3} {rE:10.4f}  {rC:10.4f}  {rt:9.4f}  {rc:9.4f}'
                  f'    (E0 {r["E0"]:.3e} N0 {r["N0"]:.3e}; C {r["C"]:.3e} Nj {r["Nj"]:.3e})')
            for chave, val in (('E0', rE), ('C', rC), ('tile', rt), ('circ', rc)):
                pior[(i, j, chave)] = max(pior.get((i, j, chave), 0.0), val)
        N0, N2 = np.array(N['N0']), np.array(N['N2'])
        cz = N0[3, 0]; c1 = N0[3, 1]
        print(f'   far <- Z  : verdade {fx["far<-Z"]:.3e}  cota {cz:.3e}  razao {fx["far<-Z"]/cz:.4f}')
        print(f'   far <- T1 : verdade {fx["far<-T1"]:.3e}  cota {c1:.3e}  razao {fx["far<-T1"]/c1:.4f}')
        czf = N0[0, 3] + r2*N2[0, 3]; c1f = N0[1, 3] + r2*N2[1, 3]
        print(f'   Z <- far  (parcial): E0 {fp["Zfar_E0"]:.3e} (N0 {N0[0,3]:.3e}); |delta|=r2: {fp["Zfar_circ"]:.3e}'
              f'  cota {czf:.3e}  razao {fp["Zfar_circ"]/czf:.4f}')
        print(f'   T1 <- far (parcial): E0 {fp["T1far_E0"]:.3e} (N0 {N0[1,3]:.3e}); |delta|=r2: {fp["T1far_circ"]:.3e}'
              f'  cota {c1f:.3e}  razao {fp["T1far_circ"]/c1f:.4f}')
        for chave, val in (('far<-Z', fx['far<-Z']/cz), ('far<-T1', fx['far<-T1']/c1),
                           ('Z<-far parcial', fp['Zfar_circ']/czf), ('T1<-far parcial', fp['T1far_circ']/c1f)):
            pior[('far', chave, '')] = max(pior.get(('far', chave, ''), 0.0), val)
        saida.append(dict(cfg=[str(x) for x in cfg], res={f'{i}<-{j}': r for (i, j), r in res.items()},
                          far_exato=fx, far_parcial=fp, N0=N['N0'], NZ=N['NZ'], N1=N['N1'], N2=N['N2']))
        (AQUI/'r3_resultados.json').write_text(json.dumps(saida, indent=1, default=str))
    print('\n== maior razao verdade/cota por par (todas as configuracoes)')
    for k, v in sorted(pior.items()):
        print(f'   {k}: {v:.4f}')
    viol = {k: v for k, v in pior.items() if v > 1.0}
    print('VIOLACOES:', viol if viol else 'nenhuma', f'  [{time.time()-t0:.0f}s]')


if __name__ == '__main__':
    main()
