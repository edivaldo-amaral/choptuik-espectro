#!/usr/bin/env python3
"""S3b, rota local: certificado de Newton-Kantorovich (NK) BORDEJADO do autovalor de L perto
de 0,401. Ver docs/S3B_ROTA.md.

Incognitas x = (y, lambda), h = Q0 y, Q0 = (O_mu + lambda0)^-1 (bloco-diagonal em m, triangular
superior em n); F(x) = (L(lambda) h, l(h) - 1), L(s) = L0 + s. Na norma do toro de L (raios
(65/64, 5/4), m com sinal), pesos kappa1^m sqrt(omega2(n)).

Jacobiano pre-condicionado: DF(x0) = [[H, h0], [l Q0, 0]], H = I + B Q0 (B = mu S Gamma1 +
2 S Xi Gamma2). Niveis:
  Z'   = Z + {lambda}: inverso bordejado denso V_b;
  T    = F \\ Z, em JANELAS de W modos: blocos exatos janela <- janela; alfa_T = ||(normas)||_2;
  far  = complemento de F: ||B_L|| (simbolo, em Arb) vezes ||Q0 P_far|| (forma fechada, kappa2 = 5/4).
A = diag(V_b, I, I); Perron nos 3 niveis; NK: ||T(x) - x0||_v <= Y + (theta0 + z1 r) r <= r.

ESTADO: estimativa em ponto flutuante (parte 1: infraestrutura e nivel Z').
"""
from __future__ import annotations

import argparse
from functools import lru_cache
import math
from pathlib import Path
import sys
import time

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import signed_operator_L as sl
import closed_form as cf
from tile_certificate import norma2

K1, K2, MU = sl.K1, sl.K2, sl.MU
NORMA_BL = 3.964043             # ||B_L(RefA)|| no toro (65/64, 5/4): t_norma = 1937/500 mais a folga do simbolo
                                # em Arb (0,0900420665) da 3,9640420665 <= 3,964043. E a cota do fundo DOS DADOS: o
                                # erro de RT entra a parte (rt, CHOPTUIK_EPS_FUNDO_L), em beta = (NORMA_BL + rt) Qfar
                                # etc. Nao confundir com 3,9641333 de symbol-L-2048x4096.json, que recombina as
                                # faixas com a folga gravada (usada no L-real). Revisao do NK em 0,7332, 09/10, I-1.
BAND_M, BAND_N = 40, 100


def pesos(r: sl.Radial):
    return K1**r.m*np.sqrt(np.concatenate([cf.omega(ns, K2) for _, ns in r.comps]))


class Contexto:
    def __init__(self, lam0, N, g=None):
        self.g = sl.campos() if g is None else g
        self.lam0 = lam0
        self.N = N                  # truncagem radial comum (o F e Mc x N)
        self._rad = {}

    def rad(self, m, N=None):
        N = self.N if N is None else N
        k = (m, N)
        if k not in self._rad:
            self._rad[k] = sl.Radial(m, N)
        return self._rad[k]

    @lru_cache(maxsize=40)
    def Q(self, m, N):
        """(O_mu + lambda0)^-1 do modo m truncado em n < N, pesado (sem o fator kappa1^m, que cancela)."""
        r = self.rad(m, N)
        w = np.sqrt(np.concatenate([cf.omega(ns, K2) for _, ns in r.comps]))
        J = sl.J_livre(r, self.lam0)
        return (w[:, None]*np.linalg.inv(J))/w[None, :]

    def B(self, mo, mi, No, Ni):
        """Bloco pesado (modo mo, n < No) <- (modo mi, n < Ni) do fundo. So depende de
        d = mo - mi, da paridade de mi e das truncagens (os pesos de Fourier dao kappa1^d)."""
        return self._B(mo - mi, mi % 2, No, Ni)

    @lru_cache(maxsize=90)
    def _B(self, d, par, No, Ni):
        ri = sl.Radial(par, Ni); ro = sl.Radial(par + d, No)
        X = sl.bloco_B(self.g, ro, ri)
        return (pesos(ro)[:, None]*X)/pesos(ri)[None, :]


def nivel_Z(ctx: Contexto, MZ, NZ, log=print, sem_borda=False, deflacao=(), phi='biortogonal'):
    """Autopar da truncagem Z e o bloco bordejado. Devolve V_b, h0 (pesado, normalizado),
    layout de Z e a linha l Q0 restrita a Z.
    sem_borda (Rouche de operadores no contorno, C1_REAVALIACAO.md): V = Hh_ZZ^-1, sem autopar nem borda.
    deflacao = raizes s_k da truncagem a deflacionar (Brauer/Wielandt): L0 + delta_k v_k u_k^*, com
    v_k, u_k autovetores direito e esquerdo de L0 (u^* v = 1) e delta_k = 1 - lambda_k (a raiz vai
    para s = -1). So mexe nas linhas de Z: Hh_ZZ ganha delta (W v)(W^-1 u)^* Q0, e o acoplamento
    Z <- T ganha delta (W v) (W^-1 u)^* Q0 P_T (em acoplamentos_Z).
    phi = 'biortogonal' (padrao, 26/09): os phi_k em span{v_j} com phi_i^* v_j = delta_ij (nos pesos). Os
    autovetores de fase e de gauge sao quase paralelos (cos 0,87); a deflacao conjunta biortogonal leva as
    duas raizes exatamente a -1 com ||E|| ~ 1,3 e resolvente ~40 no contorno (contra ~190 com phi = v, em que
    o acoplamento das duas deflacoes deixa uma raiz em -0,056).
    phi = 'direito': u_w = v_w/||v_w||^2 nos pesos (Brauer vale para qualquer phi
    com phi^* v = 1; ||delta v phi^*|| = delta, e o resolvente no contorno fica menor que com o autovetor
    esquerdo, que da ||E|| ~ 36); phi = 'esquerdo': autovetor esquerdo (o piloto de 26/09 usou este)."""
    ms = sl.modos(MZ)
    rs = [ctx.rad(m, NZ) for m in ms]
    offs = np.cumsum([0] + [r.dim for r in rs])
    n = offs[-1]
    # L0 na caixa Z (sem pesos) para o autopar
    L0 = np.zeros((n, n), complex)
    for j, ri in enumerate(rs):
        L0[offs[j]:offs[j+1], offs[j]:offs[j+1]] += sl.J_livre(ri, 0.0)
        for i, ro in enumerate(rs):
            if abs(ro.m - ri.m) <= BAND_M:
                L0[offs[i]:offs[i+1], offs[j]:offs[j+1]] += sl.bloco_B(ctx.g, ro, ri)
    from scipy.linalg import lu_factor, lu_solve
    w = np.concatenate([pesos(r) for r in rs])
    defl = []
    for sk in deflacao:                            # autovetores direito e esquerdo de L0 em -s_k
        lu = lu_factor(L0 + (sk + 1e-7)*np.eye(n))
        v = np.random.default_rng(2).standard_normal(n) + 0j; u = v.copy()
        for _ in range(8):
            v = lu_solve(lu, v); v /= np.linalg.norm(v)
            u = lu_solve(lu, u, trans=2); u /= np.linalg.norm(u)
        del lu
        u = u/np.conj(np.vdot(u, v))               # u^* v = 1
        lamk = np.vdot(u, L0 @ v)                  # autovalor de L0 (raiz s = -lamk)
        vw = w*v
        uw = u/w if phi == 'esquerdo' else vw/np.vdot(vw, vw).real          # 'biortogonal': refeito abaixo
        defl.append((1.0 - lamk, vw, uw))
        log(f'  deflacao: raiz {-lamk.real:.8f}{-lamk.imag:+.2e}i -> -1; ||v_w|| {np.linalg.norm(w*v):.3g}, ||u_w|| {np.linalg.norm(u/w):.3g}')
    if phi == 'biortogonal' and len(defl) > 1:      # phi_i^* v_j = delta_ij em span{v_j} (pesos)
        X = np.stack([vw for _, vw, _ in defl], 1)
        Phi = X @ np.linalg.inv(X.conj().T @ X)
        defl = [(dk, vw, Phi[:, k]) for k, (dk, vw, _) in enumerate(defl)]
    if sem_borda:
        del L0; x = None; lam = complex(np.nan)
    else:
        # iteracao inversa em L0 + lambda0 (autovalor simples e isolado perto de lambda0)
        lu = lu_factor(L0 + ctx.lam0*np.eye(n))
        x = np.random.default_rng(1).standard_normal(n) + 0j
        for _ in range(6):
            x = lu_solve(lu, x); x /= np.linalg.norm(x)
        lam = -(x.conj() @ (L0 @ x))/(x.conj() @ x)
        del lu, L0
    h0w = None
    if x is not None:
        h0w = w*x
        h0w /= np.linalg.norm(h0w)
    # H_ZZ = I + B Q0 (pesado) e l Q0
    H = np.eye(n, dtype=complex)
    for j, ri in enumerate(rs):
        Qj = ctx.Q(ri.m, NZ)
        for i, ro in enumerate(rs):
            if abs(ro.m - ri.m) <= BAND_M:
                H[offs[i]:offs[i+1], offs[j]:offs[j+1]] += ctx.B(ro.m, ri.m, NZ, NZ) @ Qj
    for dk, vw, uw in defl:                        # Brauer: + delta (W v)(W^-1 u)^* Q0 nas linhas de Z
        uQ = np.concatenate([uw[offs[j]:offs[j+1]].conj() @ ctx.Q(r.m, NZ) for j, r in enumerate(rs)])
        H += dk*np.outer(vw, uQ)
    if sem_borda:
        Mb = H
    else:
        lQ = np.concatenate([h0w[offs[j]:offs[j+1]].conj() @ ctx.Q(r.m, NZ) for j, r in enumerate(rs)])
        Mb = np.zeros((n + 1, n + 1), complex)
        Mb[:n, :n] = H; del H
        Mb[:n, n] = h0w; Mb[n, :n] = lQ
    Vb = np.linalg.inv(Mb)
    X = Vb @ Mb; del Mb
    X[np.diag_indices(X.shape[0])] -= 1
    rho = norma2(X); del X
    log(f'  Z = {MZ}x{NZ} ({n}): lambda da truncagem {lam.real:.12f}; ||V_b|| {norma2(Vb):.3f}; rho {rho:.1e}'
        + ('  [sem borda]' if sem_borda else ''))
    return dict(Vb=Vb, h0w=h0w, ms=ms, rs=rs, offs=offs, lam=lam, rho=rho, n=n, defl=defl, borda=not sem_borda)


def colsT(ctx, m, MZ, NZ, N):
    """Indices locais (no modo m truncado em N) das colunas/linhas de T = F \\ Z."""
    r = ctx.rad(m, N)
    n = np.concatenate([ns for _, ns in r.comps])
    return np.where(n >= NZ)[0] if abs(m) < MZ else np.arange(r.dim)


def nZ_local(ctx, m, NZ, N):
    r = ctx.rad(m, N)
    n = np.concatenate([ns for _, ns in r.comps])
    return np.where(n < NZ)[0]


def acoplamentos_Z(ctx, z, MZ, NZ, Mc, Nc, log=print):
    """a = ||V_b [P_Z B Q0 P_T; l Q0 P_T]|| (Gram conjunto) e b = ||P_T B Q0 P_Z||."""
    nZ = z['n']; ms, offs = z['ms'], z['offs']; idx = {m: i for i, m in enumerate(ms)}
    nb = nZ + 1 if z.get('borda', True) else nZ
    Ga = np.zeros((nb, nb), complex)
    for m in range(-(min(MZ + BAND_M, Mc) - 1), min(MZ + BAND_M, Mc)):
        c = colsT(ctx, m, MZ, NZ, Nc)
        RT = ctx.Q(m, Nc)[:, c]
        Y = np.zeros((nb, len(c)), complex)
        for mo in ms:
            if abs(mo - m) <= BAND_M:
                i = idx[mo]
                Y[offs[i]:offs[i+1]] = ctx.B(mo, m, NZ, Nc) @ RT
        if m in idx:                                  # l Q0: descida radial para dentro de Z
            i = idx[m]
            if z.get('borda', True):
                Y[nZ] = z['h0w'][offs[i]:offs[i+1]].conj() @ RT[nZ_local(ctx, m, NZ, Nc)]
            for dk, vw, uw in z.get('defl', []):      # Brauer: delta (W v) (W^-1 u)^* Q0 P_T
                Y[:nZ] += dk*np.outer(vw, uw[offs[i]:offs[i+1]].conj() @ RT[nZ_local(ctx, m, NZ, Nc)])
        Ga += Y @ Y.conj().T
    Vb = z['Vb']
    rng = np.random.default_rng(0); x = rng.standard_normal(nb) + 0j; x /= np.linalg.norm(x); lam = 0.0
    for _ in range(80):                        # maior autovalor de V_b Ga V_b^* sem forma-lo
        y = Vb @ (Ga @ (Vb.conj().T @ x)); lam = float(np.linalg.norm(y)); x = y/lam
    a = math.sqrt(lam)
    del Ga
    Gb = np.zeros((nZ, nZ), complex)
    NZb = NZ + BAND_N + 1
    for r_ in range(-(min(MZ + BAND_M, Mc) - 1), min(MZ + BAND_M, Mc)):
        linhas = colsT(ctx, r_, MZ, NZ, NZb)
        if not len(linhas):
            continue
        X = np.zeros((len(linhas), nZ), complex)
        for j, mi in enumerate(ms):
            if abs(r_ - mi) <= BAND_M:
                X[:, offs[j]:offs[j+1]] = ctx.B(r_, mi, NZb, NZ)[linhas] @ ctx.Q(mi, NZ)
        Gb += X.conj().T @ X
    b = math.sqrt(max(np.linalg.eigvalsh(Gb)[-1], 0.0))
    log(f'  Z\' <- T: a = {a:.4f};  T <- Z: b = {b:.4f}')
    return a, b


def _potencia(mv, rmv, n, iters=60, seed=0):
    """Maior valor singular por iteracao de potencia com produtos implicitos (estimativa)."""
    rng = np.random.default_rng(seed)
    v = rng.standard_normal(n) + 1j*rng.standard_normal(n); v /= np.linalg.norm(v)
    s = 0.0
    for _ in range(iters):
        u = rmv(mv(v)); nu = np.linalg.norm(u)
        if nu == 0:
            return 0.0
        v = u/nu
    return float(np.linalg.norm(mv(v)))


def janelas(ctx, MZ, NZ, Mc, Nc, W, log=print):
    """Blocos janela <- janela de P_T B Q0 P_T (linhas em F) e far <- janela, por PRODUTOS
    IMPLICITOS modo a modo (sem montar o bloco denso: cada um teria ~2 GB).
    Devolve (alfa_T = ||N_win||_2, alfa_fT = ||F_win||_2, ||Q0 P_T||)."""
    ms = list(range(-(Mc - 1), Mc))
    jan = [ms[k:k+W] for k in range(0, len(ms), W)]
    nj = len(jan)
    Nw = np.zeros((nj, nj))
    Nf = np.zeros((nj + 2, nj))
    No = Nc + BAND_N + 1
    QT = 0.0
    def linhas_T(mo):
        return colsT(ctx, mo, MZ, NZ, Nc)
    def linhas_far(mo):
        rr = ctx.rad(mo, No); n = np.concatenate([ns for _, ns in rr.comps])
        return np.where(n >= Nc)[0] if abs(mo) < Mc else np.arange(rr.dim)
    for J, modosJ in enumerate(jan):
        cols = {m: colsT(ctx, m, MZ, NZ, Nc) for m in modosJ}
        QTm = {m: ctx.Q(m, Nc)[:, cols[m]] for m in modosJ}
        for m in modosJ:
            QT = max(QT, norma2(QTm[m]))
        oj = np.cumsum([0] + [len(cols[m]) for m in modosJ]); nJ = int(oj[-1])
        def blocos_linhas(modosI, fsel, NoI):
            sel = {mo: fsel(mo) for mo in modosI}
            oi = np.cumsum([0] + [len(sel[mo]) for mo in modosI])
            def mv(v):
                u = {mi: QTm[mi] @ v[oj[q]:oj[q+1]] for q, mi in enumerate(modosJ)}
                out = np.zeros(int(oi[-1]), complex)
                for p_, mo in enumerate(modosI):
                    acc = np.zeros(len(sel[mo]), complex)
                    for mi in modosJ:
                        if abs(mo - mi) <= BAND_M:
                            acc += ctx.B(mo, mi, NoI, Nc)[sel[mo]] @ u[mi]
                    out[oi[p_]:oi[p_+1]] = acc
                return out
            def rmv(w):
                t = {mi: np.zeros(QTm[mi].shape[0], complex) for mi in modosJ}
                for p_, mo in enumerate(modosI):
                    wm = w[oi[p_]:oi[p_+1]]
                    for mi in modosJ:
                        if abs(mo - mi) <= BAND_M:
                            t[mi] += ctx.B(mo, mi, NoI, Nc)[sel[mo]].conj().T @ wm
                return np.concatenate([QTm[mi].conj().T @ t[mi] for mi in modosJ])
            return mv, rmv, int(oi[-1])
        for I, modosI in enumerate(jan):
            if min(abs(x - y) for x in modosI for y in modosJ) > BAND_M:
                continue
            mv, rmv, nI = blocos_linhas(modosI, linhas_T, Nc)
            if nI:
                Nw[I, J] = _potencia(mv, rmv, nJ)
            mv, rmv, nI = blocos_linhas(modosI, linhas_far, No)
            if nI:
                Nf[I, J] = _potencia(mv, rmv, nJ)
        for lado, alvo in ((0, list(range(Mc, Mc + BAND_M))), (1, list(range(-(Mc + BAND_M - 1), -Mc + 1)))):
            alvo = [mo for mo in alvo if min(abs(mo - mi) for mi in modosJ) <= BAND_M]
            if alvo:
                mv, rmv, nI = blocos_linhas(alvo, linhas_far, No)
                Nf[nj + lado, J] = _potencia(mv, rmv, nJ)
        log(f'    janela {J+1}/{nj} (modos {modosJ[0]}..{modosJ[-1]}): diagonal {Nw[J, J]:.4f}', flush=True)
    aT = float(np.linalg.norm(Nw, 2)); afT = float(np.linalg.norm(Nf, 2))
    log(f'  T <- T por janelas (W = {W}): {aT:.4f}  (maior bloco {Nw.max():.4f});  far <- T: {afT:.4f};  ||Q0 P_T|| {QT:.3f}')
    return aT, afT, QT, Nw


def far(ctx, MZ, NZ, Mc, Nc, log=print):
    rad = 0.0
    for m in range(-(Mc - 1), Mc):
        sig = ctx.lam0 + 1j*m/2
        passos = (2, 1, 1) if m % 2 == 0 else (1,)
        rad = max(rad, max(cf.cauda_Rinv(sig, Nc, p, MU, K2, 1500) for p in passos))
    fou = cf.Rinv_uniforme(Mc, MU, ctx.lam0, K2)
    Qf = max(rad, fou)
    log(f'  ||Q0 P_far|| <= max(radial {rad:.4f}, Fourier {fou:.4f});  beta = {NORMA_BL*Qf:.4f}')
    return Qf


def residuo(ctx, z, MZ, NZ, Mc, Nc):
    """F(x0) = (L(lambda0) h0, 0): partes em Z (pesada) e fora de Z (pesada), h0 com ||W h0|| = 1."""
    ms, offs, rs = z['ms'], z['offs'], z['rs']
    h0w = z['h0w']
    rZ = np.zeros(z['n'], complex); rT2 = 0.0
    NZb = NZ + BAND_N + 1
    for mo in range(-(min(MZ + BAND_M, Mc) - 1), min(MZ + BAND_M, Mc)):
        acc = np.zeros(ctx.rad(mo, NZb).dim, complex)
        for j, mi in enumerate(ms):
            if abs(mo - mi) <= BAND_M:
                acc += ctx.B(mo, mi, NZb, NZ) @ h0w[offs[j]:offs[j+1]]
        # parte livre + lambda0 so no proprio modo (dentro de Z)
        if mo in ms:
            j = ms.index(mo)
            r = ctx.rad(mo, NZ); w = np.sqrt(np.concatenate([cf.omega(ns, K2) for _, ns in r.comps]))
            Jw = (w[:, None]*sl.J_livre(r, ctx.lam0))/w[None, :]
            loc = nZ_local(ctx, mo, NZ, NZb)
            acc[loc] += Jw @ h0w[offs[j]:offs[j+1]]
            rZ[offs[j]:offs[j+1]] = acc[loc]
        t = colsT(ctx, mo, MZ, NZ, NZb)
        rT2 += float(np.sum(np.abs(acc[t])**2))
    return rZ, math.sqrt(rT2)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--lam0', type=complex, default=0.4010247330)
    ap.add_argument('--sem-borda', action='store_true', help='nivel Z sem borda (Rouche no contorno)')
    ap.add_argument('--deflacao', default='', help='raizes s_k a deflacionar (Brauer), ex.: "0,0.16831"')
    ap.add_argument('--Z', default='12x48')
    ap.add_argument('--F', default='40x200')
    ap.add_argument('--W', type=int, default=16)
    ap.add_argument('--dados-Z', default=None,
                    help='a,b,||V_b||,rho,Y_Z,r_T de uma rodada anterior (pula o nivel Z; so estimativa)')
    ap.add_argument('--so-Z', action='store_true', help='so o nivel Z (imprime DADOS-Z e sai)')
    a = ap.parse_args()
    MZ, NZ = map(int, a.Z.split('x')); Mc, Nc = map(int, a.F.split('x'))
    ctx = Contexto(a.lam0, Nc)
    t0 = time.time()
    tempo = lambda: f'[{time.time()-t0:.0f}s]'
    if a.dados_Z:
        aZ, bZ, nVb, rho, YZ, rT = map(float, a.dados_Z.split(','))
        z = dict(rho=rho)
        print(f'  nivel Z reaproveitado: a {aZ}, b {bZ}, ||V_b|| {nVb}, rho {rho}, Y_Z {YZ}, r_T {rT}', flush=True)
    else:
        dfl = [complex(x) for x in a.deflacao.split(',') if x.strip()]
        z = nivel_Z(ctx, MZ, NZ, sem_borda=a.sem_borda, deflacao=dfl); print('  ', tempo(), flush=True)
        aZ, bZ = acoplamentos_Z(ctx, z, MZ, NZ, Mc, Nc); print('  ', tempo(), flush=True)
        # tudo o que usa V_b antes das janelas, e depois libera V_b (memoria)
        nVb = norma2(z['Vb'])
        if a.sem_borda:                              # Rouche: nao ha residuo
            rZ, rT, YZ = None, 0.0, 0.0
        else:
            rZ, rT = residuo(ctx, z, MZ, NZ, Mc, Nc)
            YZ = norma2(z['Vb'][:, :z['n']] @ rZ[:, None]) if np.any(rZ) else 0.0
        print(f'  ||V_b|| {nVb:.3f};  residuo: Z {YZ:.2e}, T {rT:.2e}   {tempo()}', flush=True)
        del z['Vb']
        ctx.Q.cache_clear(); ctx._B.cache_clear()
        print(f"  DADOS-Z {aZ:.6g},{bZ:.6g},{nVb:.6g},{z['rho']:.3g},{YZ:.3g},{rT:.3g}  (lambda da truncagem {z['lam'].real:.12f})", flush=True)
        if a.so_Z:
            return
    aT, afT, QT, Nw = janelas(ctx, MZ, NZ, Mc, Nc, a.W); print('  ', tempo(), flush=True)
    Qf = far(ctx, MZ, NZ, Mc, Nc); beta = NORMA_BL*Qf
    epsZ = nVb*NORMA_BL*(cf.descida(NZ + BAND_N + 1, Nc, MU, K2) + 6*cf.descida_divxi(NZ, Nc, K2))
    N0 = np.array([[z['rho'], aZ, epsZ], [bZ, aT, beta], [0.0, afT, beta]])
    lam, vec = np.linalg.eig(N0); i = int(np.argmax(np.abs(lam))); v = np.abs(vec[:, i]); v /= v.max()
    theta0 = float(max((N0 @ v)/v))
    Y = max(YZ/v[0], rT/v[1])
    np.set_printoptions(precision=4, suppress=True)
    print(f'  N0 (Zb, T, far) =\n{N0}\n  v = {v};  theta0 = {theta0:.4f}')
    print(f'  residuo: Z\' {YZ:.2e}, T {rT:.2e};  Y = {Y:.2e};  raio NK ~ Y/(1 - theta0) = '
          f'{Y/(1 - theta0) if theta0 < 1 else float("inf"):.2e}   {tempo()}')


if __name__ == '__main__':
    main()
