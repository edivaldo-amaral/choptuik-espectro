#!/usr/bin/env python3
"""S3b, L, modo rigoroso, ETAPA A: o nivel Z' (bordejado) do certificado NK de L.

x = (y, lambda), h = Q0 y, F(x) = (L(lambda) h, l(h) - 1), Q0 = (J_mu + lambda0)^-1 (mu = mu_RefA,
fixo). y0 e escolhido com Q0 y0 = h0w EXATAMENTE (h0w: o vetor de ponto flutuante gravado), e
l(h) = <h0w, h> (pesado). O bloco Z' de DF(x0) e M = [[H_ZZ, h0w], [l Q0 P_Z, 0]], e o
certificado usa o inverso EXATO da matriz calculada Mh (um operador fixo). Cotas:

  t_V   >= ||Mh^-1||                         (Loewner: t^2 Mh Mh^* - I >= 0);
  a     >= ||Mh^-1 [P_Z B Q0 P_T; l Q0 P_T]|| (Loewner com Ga = Y Y^*, banda |d| <= BANDA) mais
           t_V (erro de Y + (dB_L + eps_mu_L) ||Q0 P_T||);
  rhoZ  >= ||Mh^-1 (M - Mh)|| = t_V ||M - Mh|| (arredondamento de H_ZZ e de l Q0, erro de Q0,
           erro de RT, mu; em Z a banda e a completa, |d| <= 40);
  Y_Z   >= ||Mh^-1 P_Z' F(x0)|| (solucao com correcao);
  epsZ  >= ||Mh^-1 P_Z B Q0 P_far|| (descidas; a perturbacao global vezes ||Q0 P_far||).
Grava h0w (npy) para as etapas seguintes. Pico de memoria ~10 GB (Z = 32x128: 14 016).
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys
import time

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import closed_form as cf
import fundo_L
import nk_L as nk
import signed_operator_L as sl
from tile_certificate import norma2
from verified_inverse import pd_hermitiana, gram_linhas
from verified_norms import cota_norma, gamma, U

FATOR_PESO = 1 + 1e-12
MU, K2 = sl.MU, sl.K2


def err_mm(A, B):
    return gamma(A.shape[1] + 8)*float(np.linalg.norm(A))*float(np.linalg.norm(B))


def Q_cert(m, N, lam0):
    """(J_mu + lambda0)^-1 do modo m truncado em n < N, pesado (sqrt omega), e cota de
    ||Q_exato - Q|| pelo residuo: ||Q_ex - Q|| <= ||Q|| e/(1 - e), e >= ||I - Jw Q||."""
    r = sl.Radial(m, N)
    w = np.sqrt(np.concatenate([cf.omega(ns, K2) for _, ns in r.comps]))
    Jw = (w[:, None]*sl.J_livre(r, lam0))/w[None, :]
    Q = np.linalg.inv(Jw)
    E = np.eye(len(Q)) - Jw @ Q
    e = cota_norma(E, err_mm(Jw, Q) + U*float(np.linalg.norm(E)))
    assert e < 0.5
    nQ = cota_norma(Q)
    return Q, nQ*e/(1 - e), nQ


class Log:
    def __init__(self):
        self.t0 = time.time()
    def __call__(self, s):
        print(f'{s}   [{time.time()-self.t0:.0f}s]', flush=True)



def cota_tQ(a, M, QZ, dQZ, ms, offs, n, log):
    """tQ >= ||Mh^-1 [Q_ZZ; 0]|| (a parte de Z do termo de 2a ordem de nk_L3_C, que antes era t_V ||Q_ZZ||):
    Loewner t^2 Mh Mh^* - X X^* >= 0, X = [Q_ZZ; 0] (Q_ZZ calculada, bloco-diagonal; X X^* subtraida por blocos).
    O Q_ZZ exato difere do calculado por <= dQZ: ||Mh^-1 [Q_e; 0]|| <= t + t_V dQZ (t_V de Z.json). Grava em Z.json."""
    from scipy.linalg import lu_factor, lu_solve
    z = json.loads((a.dir/'Z.json').read_text())
    lu = lu_factor(M, check_finite=False)
    v = np.random.default_rng(3).standard_normal(n) + 0j; v /= np.linalg.norm(v)
    def aplX(u):
        out = np.zeros(n + 1, complex)
        for j, m in enumerate(ms):
            out[offs[j]:offs[j+1]] = QZ[m] @ u[offs[j]:offs[j+1]]
        return out
    def aplXh(u):
        out = np.zeros(n, complex)
        for j, m in enumerate(ms):
            out[offs[j]:offs[j+1]] = QZ[m].conj().T @ u[offs[j]:offs[j+1]]
        return out
    for _ in range(40):
        u = lu_solve(lu, aplX(v)); v = aplXh(lu_solve(lu, u, trans=2)); v /= np.linalg.norm(v)
    est = float(np.linalg.norm(lu_solve(lu, aplX(v))))
    del lu
    log(f'tQ: estimativa ||Mh^-1 [Q_ZZ; 0]|| ~ {est:.4f} (t_V ||Q_ZZ|| = {z["tV"]*z["nQZ"]:.2f})')
    P, eP = gram_linhas(M); del M
    fP = float(np.linalg.norm(P))
    Sb = []; eS = 0.0; fS2 = 0.0
    for m in ms:
        Sm = QZ[m] @ QZ[m].conj().T
        Sb.append(Sm); fS2 += float(np.linalg.norm(Sm))**2
        eS = max(eS, gamma(len(Sm))*float(np.linalg.norm(QZ[m]))**2)
    fS = math.sqrt(fS2)
    t = est*(1 + 1e-3)
    for _ in range(60):
        G = (t*t)*P
        for j in range(len(ms)):
            G[offs[j]:offs[j+1], offs[j]:offs[j+1]] -= Sb[j]
        if pd_hermitiana([G], t*t*eP + eS + 2*U*(t*t*fP + fS)):
            break
        t *= 1.003
    else:
        raise RuntimeError('tQ: Loewner nao verificou')
    tQ = (t + z['tV']*dQZ)*FATOR_PESO
    z['tQ'] = tQ
    (a.dir/'Z.json').write_text(json.dumps(z, indent=1) + '\n')
    log(f'tQ RIGOROSO: ||Mh^-1 [Q_ZZ; 0]|| <= {tQ:.5f} (Loewner {t:.5f}); gravado em Z.json')

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--lam0', type=complex, default=0.4010247330)
    ap.add_argument('--Z', default='32x128'); ap.add_argument('--F', default='200x400')
    ap.add_argument('--banda', type=int, default=16)
    ap.add_argument('--chutes', default='20.7,0.375', help='estimativas de ||Mh^-1|| e de a (ponto flutuante)')
    ap.add_argument('--dir', type=Path, default=Path('build/s3b/rig'))
    ap.add_argument('--h0', type=Path, default=None, help='centro h0 (pesado, layout de Z) em .npy')
    ap.add_argument('--so-tQ', dest='so_tQ', action='store_true',
                    help='so tQ >= ||Mh^-1 [Q_ZZ; 0]|| com o h0w.npy de --dir (grava em Z.json)')
    a = ap.parse_args()
    if a.lam0.imag == 0:
        a.lam0 = a.lam0.real                          # S3b e S4: lambda0 real, como antes
    log = Log()
    MZ, NZ = map(int, a.Z.split('x')); Mc, Nc = map(int, a.F.split('x'))
    a.dir.mkdir(parents=True, exist_ok=True)
    ctx = nk.Contexto(a.lam0, Nc)
    g = ctx.g
    ms = sl.modos(MZ); rs = [ctx.rad(m, NZ) for m in ms]
    offs = np.cumsum([0] + [r.dim for r in rs]); n = int(offs[-1])
    fb = fundo_L.dB_L(g, a.banda)
    fb40 = fundo_L.dB_L(g, sl.BAND_M)                   # em Z a banda e completa: sem resto
    emu = fundo_L.eps_mu_L(complex(a.lam0))
    nBL = nk.NORMA_BL + fb['total']                   # ||B|| (simbolo) + perturbacao
    log(f'Z = {MZ}x{NZ} ({n}); dB_L banda {a.banda}: {fb["total"]:.2e}; sem resto: {fb40["total"]:.2e}; eps_mu {emu:.1e}')
    # 1. Q0 nos modos de Z (n < NZ) e nas colunas de T perto de Z (n < Nc)
    QZ = {}; dQZ = 0.0; nQZ = 0.0
    for m in ms:
        Q, dq, nq = Q_cert(m, NZ, a.lam0); QZ[m] = Q; dQZ = max(dQZ, dq); nQZ = max(nQZ, nq + dq)
    log(f'Q0 em Z: ||Q_ZZ|| <= {nQZ:.4f}, erro <= {dQZ:.1e}')
    # 2. autovetor aproximado (ponto flutuante: e so um vetor escolhido) e h0w
    L0 = np.zeros((n, n), complex)
    for j, ri in enumerate(rs):
        L0[offs[j]:offs[j+1], offs[j]:offs[j+1]] += sl.J_livre(ri, 0.0)
        for i, ro in enumerate(rs):
            if abs(ro.m - ri.m) <= sl.BAND_M:
                L0[offs[i]:offs[i+1], offs[j]:offs[j+1]] += sl.bloco_B(g, ro, ri)
    from scipy.linalg import lu_factor, lu_solve
    if a.so_tQ:
        # so tQ: o MESMO h0w.npy das etapas seguintes, sem renormalizar (Mh identica a da etapa A)
        del L0
        h0w = np.load(a.dir/'h0w.npy')
        assert len(h0w) == n
    elif a.h0 is not None:
        # centro dado (S4: o vetor de gauge (Z + 1) omega_RefA restrito a Z, pesado)
        del L0
        h0w = np.load(a.h0).astype(complex); h0w /= np.linalg.norm(h0w)
        assert len(h0w) == n
    else:
        # iteracao inversa NA MATRIZ PESADA (bem escalada): a sem pesos tem erro absoluto ~u ||x|| em todas as
        # componentes, e os pesos sqrt(omega(n)) kappa1^m (ate ~1e12 em n = 128) o amplificam nas de n alto
        # (faixa B: Y_Z 2e-4 em vez de ~1e-10)
        del L0
        Lw = np.zeros((n, n), complex)
        for j, ri in enumerate(rs):
            wi = np.sqrt(np.concatenate([cf.omega(ns, K2) for _, ns in ri.comps]))
            Lw[offs[j]:offs[j+1], offs[j]:offs[j+1]] += (wi[:, None]*sl.J_livre(ri, a.lam0))/wi[None, :]
            for i, ro in enumerate(rs):
                if abs(ro.m - ri.m) <= sl.BAND_M:
                    Lw[offs[i]:offs[i+1], offs[j]:offs[j+1]] += ctx.B(ro.m, ri.m, NZ, NZ)
        ctx._B.cache_clear()
        lu = lu_factor(Lw, overwrite_a=True, check_finite=False); del Lw
        x = np.random.default_rng(1).standard_normal(n) + 0j
        for _ in range(6):
            x = lu_solve(lu, x); x /= np.linalg.norm(x)
        del lu
        h0w = x/np.linalg.norm(x)
    if not a.so_tQ:
        np.save(a.dir/'h0w.npy', h0w)
        log('autovetor aproximado gravado (h0w.npy)')
    # 3. Mh = [[H_ZZ, h0w], [l Q0, 0]] e o erro ||M - Mh||
    M = np.zeros((n + 1, n + 1), complex)
    M[np.diag_indices(n)] = 1.0
    e2 = 0.0
    for j, ri in enumerate(rs):
        Qj = QZ[ri.m]
        for i, ro in enumerate(rs):
            if abs(ro.m - ri.m) <= sl.BAND_M:
                Bb = ctx.B(ro.m, ri.m, NZ, NZ)
                M[offs[i]:offs[i+1], offs[j]:offs[j+1]] += Bb @ Qj
                e2 += err_mm(Bb, Qj)**2
    ctx._B.cache_clear()
    M[:n, n] = h0w
    for j, r in enumerate(rs):
        M[n, offs[j]:offs[j+1]] = h0w[offs[j]:offs[j+1]].conj() @ QZ[r.m]
    fM = float(np.linalg.norm(M))
    dM = (math.sqrt(e2) + U*fM + (fb40['arred'] + fb40['rt'])*nQZ + nBL*dQZ + emu
          + dQZ*(1 + n*U) + gamma(3*NZ)*nQZ)         # linha l Q0: erro de Q e do produto
    log(f'Mh montada: ||M - Mh|| <= {dM:.2e}')
    if a.so_tQ:
        cota_tQ(a, M, QZ, dQZ, ms, offs, n, log)
        return
    # 4. residuo F(x0) em Z': (J + lambda0) h0 + B h0 nas linhas de Z, e l(h0) - 1
    rZ = np.zeros(n + 1, complex); er2 = 0.0
    for i, ro in enumerate(rs):
        wo = np.sqrt(np.concatenate([cf.omega(ns, K2) for _, ns in ro.comps]))
        Jw = (wo[:, None]*sl.J_livre(ro, a.lam0))/wo[None, :]
        acc = Jw @ h0w[offs[i]:offs[i+1]]; er2 += err_mm(Jw, h0w[offs[i]:offs[i+1], None])**2
        for j, ri in enumerate(rs):
            if abs(ro.m - ri.m) <= sl.BAND_M:
                Bb = ctx.B(ro.m, ri.m, NZ, NZ)
                acc += Bb @ h0w[offs[j]:offs[j+1]]; er2 += err_mm(Bb, h0w[offs[j]:offs[j+1], None])**2
        rZ[offs[i]:offs[i+1]] = acc
    ctx._B.cache_clear()
    rZ[n] = np.vdot(h0w, h0w) - 1
    # mu: ||delta_mu (J1 + S Gamma1) h0|| <= DELTA_MU (||J1 h0|| + 2 norma_divxi ||h0||)
    J1h = 0.0
    for i, r in enumerate(rs):
        wo = np.sqrt(np.concatenate([cf.omega(ns, K2) for _, ns in r.comps]))
        J1 = (sl.J_livre(r, 0.0, mu=1.0) - 1j*r.m/2*np.eye(r.dim))
        J1h += float(np.linalg.norm(((wo[:, None]*J1)/wo[None, :]) @ h0w[offs[i]:offs[i+1]]))**2
    from tile_perron import DELTA_MU
    dr = (math.sqrt(er2) + U*float(np.linalg.norm(rZ)) + (fb40['arred'] + fb40['rt'])*1.0001
          + DELTA_MU*(math.sqrt(J1h)*1.0001 + 2*cf.norma_divxi(K2)) + 4*n*U)
    lu = lu_factor(M, check_finite=False)
    xs = lu_solve(lu, rZ)
    res = M @ xs - rZ
    eres = gamma(n + 1)*fM*float(np.linalg.norm(xs)) + U*float(np.linalg.norm(res))
    nres = float(np.linalg.norm(res))*(1 + 1e-6) + eres
    nxs = float(np.linalg.norm(xs))*(1 + 1e-6)
    # chute de ||Mh^-1|| por potencia com o LU
    v = np.random.default_rng(2).standard_normal(n + 1) + 0j; v /= np.linalg.norm(v)
    for _ in range(40):
        u = lu_solve(lu, v); v = lu_solve(lu, u, trans=2); v /= np.linalg.norm(v)
    chV = float(np.linalg.norm(lu_solve(lu, v)))
    del lu
    log(f'residuo: ||Mh^-1 r|| ~ {nxs:.2e} (residuo da solucao {nres:.1e}, erro de r {dr:.1e}); ||Mh^-1|| ~ {chV:.4f}')
    # 5. P = Mh Mh^* (gravado) e t_V
    P, eP = gram_linhas(M); del M
    np.save(a.dir/'P.npy', P)
    fP = float(np.linalg.norm(P))
    t = chV*(1 + 1e-4)
    while True:
        P *= t*t
        P[np.diag_indices(n + 1)] -= 1.0
        ok = pd_hermitiana(P, t*t*eP + 2*U*(t*t*fP + math.sqrt(n + 1)))
        if ok:
            break
        P = np.load(a.dir/'P.npy'); t *= 1.002
    tV = t*FATOR_PESO
    del P
    log(f't_V = ||Mh^-1|| <= {tV:.5f}')
    YZ = (nxs + tV*(nres + dr))*FATOR_PESO
    rhoZ = tV*dM*FATOR_PESO
    # 6. a: t_a^2 P - S >= 0, S = Y Y^*, Y = [P_Z B Q0 P_T; l Q0 P_T] (banda |d| <= banda).
    #    S e acumulado uma vez (lotes de modos) e gravado; cada tentativa monta t^2 P - S no lugar.
    chA = float(a.chutes.split(',')[1])
    idx = {m: i for i, m in enumerate(ms)}
    eY2 = 0.0; fY2 = 0.0; eS = 0.0; nQT = 0.0; dQT = 0.0
    S = np.zeros((n + 1, n + 1), complex)
    lote = []
    def descarrega():
        nonlocal S
        if lote:
            Yb = np.concatenate(lote, axis=1); lote.clear()
            S += Yb @ Yb.conj().T
            del Yb
    for mi in range(-(MZ + a.banda - 1), MZ + a.banda):
        c = nk.colsT(ctx, mi, MZ, NZ, Nc)
        Q, dq, nq = Q_cert(mi, Nc, a.lam0)
        RT = Q[:, c]
        nQT = max(nQT, nq + dq); dQT = max(dQT, dq)
        Y = np.zeros((n + 1, len(c)), complex)
        for mo in ms:
            if abs(mo - mi) <= sl.BAND_M:          # banda COMPLETA nas colunas perto de Z: sem resto
                i = idx[mo]; Bb = ctx.B(mo, mi, NZ, Nc)
                Y[offs[i]:offs[i+1]] = Bb @ RT; eY2 += err_mm(Bb, RT)**2
        if mi in idx:
            i = idx[mi]
            Y[n] = h0w[offs[i]:offs[i+1]].conj() @ RT[nk.nZ_local(ctx, mi, NZ, Nc)]
            eY2 += err_mm(h0w[None, offs[i]:offs[i+1]], RT)**2
        fy = float(np.linalg.norm(Y))**2; fY2 += fy
        lote.append(Y)
        if sum(x.shape[1] for x in lote) >= 6000:
            descarrega()
        ctx.Q.cache_clear(); ctx._B.cache_clear()
    descarrega()
    eS = gamma(8000)*fY2                              # produtos em lotes (<= 8000 colunas) e somas
    eS += U*float(np.linalg.norm(S))*(2*MZ + 2*a.banda)
    np.save(a.dir/'S.npy', S); del S
    log(f'S = Y Y^* acumulado (||Y||_F^2 = {fY2:.3e})')
    ta = chA*1.01
    while True:
        G = np.load(a.dir/'P.npy'); G *= ta*ta
        Sm = np.load(a.dir/'S.npy', mmap_mode='r')
        for i0 in range(0, n + 1, 2000):
            G[i0:i0+2000] -= Sm[i0:i0+2000]
        del Sm
        ok = pd_hermitiana(G, ta*ta*eP + eS + 2*U*(ta*ta*fP + fY2))
        del G
        log(f'  a: t = {ta:.5f} -> ' + ('verificado' if ok else 'nao'))
        if ok:
            break
        ta *= 1.01
    # colunas com MZ + banda <= |m| < MZ + 40: so |d| > banda as liga a Z; entram em quadratura,
    # ||Mh^-1 [Y1 Y2]|| <= sqrt(||Mh^-1 Y1||^2 + ||Mh^-1 Y2||^2)
    nQx = max(Q_cert(m, Nc, a.lam0)[2] for m in (MZ + a.banda, -(MZ + a.banda)))
    a_extra = tV*fb['resto']*nQx*1.01                 # ||Q0|| decresce com |m|: o maior e o da borda
    dY = math.sqrt(eY2) + nBL*dQT + fb40['arred']*nQT + dQT    # produto, erro de Q, montagem de B, linha l
    aZ = (math.sqrt(ta*ta + a_extra*a_extra) + tV*(dY + (fb40['rt'] + emu)*nQT))*FATOR_PESO
    # 7. Z' <- far
    Qf = max(max(cf.cauda_Rinv(a.lam0 + 1j*m/2, Nc, p, MU, K2, 1500) for p in ((2, 1, 1) if m % 2 == 0 else (1,)))
             for m in range(-(Mc - 1), Mc))
    Qf = max(Qf, cf.Rinv_uniforme(Mc, MU, a.lam0, K2))
    epsZ = tV*(nBL*(cf.descida(NZ + nk.BAND_N + 1, Nc, MU, K2) + 6*cf.descida_divxi(NZ, Nc, K2))
               + (fb['total'] + emu)*Qf)*FATOR_PESO
    out = dict(Z=[MZ, NZ], F=[Mc, Nc], lam0=([a.lam0.real, a.lam0.imag] if isinstance(a.lam0, complex) else a.lam0), banda=a.banda, n=n, tV=tV, aZ=aZ, a_loewner=ta, a_extra=a_extra,
               J1h=math.sqrt(J1h)*1.0001,
               erro_Y=dY, rhoZ=rhoZ, YZ=YZ, epsZ=epsZ, dM=dM, dB_L=fb, eps_mu=emu, Qfar=Qf,
               nQZ=nQZ, dQZ=dQZ, nQT=nQT, dQT=dQT)
    (a.dir/'Z.json').write_text(json.dumps(out, indent=1) + '\n')
    log(f'Z\' RIGOROSO: ||Mh^-1|| <= {tV:.5f}, a <= {aZ:.5f}, rho <= {rhoZ:.2e}, Y_Z <= {YZ:.2e}, eps_far <= {epsZ:.2e}')
    (a.dir/'P.npy').unlink(); (a.dir/'S.npy').unlink()


if __name__ == '__main__':
    main()
