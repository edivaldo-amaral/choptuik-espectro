#!/usr/bin/env python3
"""Nível Z do Rouché de L no contorno (C1_REAVALIACAO.md §7), em ponto flutuante (estimativa).

Subcomandos:
  schur    monta A = W L~_Z(0) W^-1 (deflação de Brauer da fase e de mu) e grava U, T (A = U T U^*),
           os autovalores (raízes s = -diag T) e os vetores de deflação;
  fator    para um centro de cauda c: range finder aleatório (em fluxo por modo de T) da matriz
           empilhada S = [[Y~ Q0_TT(c)], [J_ZT Q0_TT(c)]] (2n x |T|): S ~ Qb Qb^* S, com
           F = Qb (C C^*)^(1/2), C = Qb^* S, e resto_F = ||S - Qb Qb^* S||_F;
  pontos   para pontos p (e o centro c do fator): a_esc(p) = ||J_Z(c)[R(p) F_t - Q_ZZ(p) F_b]||
           (+ resto), b1(p) = ||J_Z(c)[R(p)^2 F_t - Q_ZZ(p)^2 F_b]|| (primeira ordem) e
           ||R(p)|| (potência). R(p) = (A + p)^-1 via Schur."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys
import time

import numpy as np
from scipy.linalg import schur, solve_triangular

sys.path.insert(0, str(Path(__file__).resolve().parent))
import nk_L as nk
import signed_operator_L as sl


def log(t0, *x):
    print(*x, f'[{time.time() - t0:.0f}s]', flush=True)


def layout(MZ, NZ, Nc, c=0.0):
    ctx = nk.Contexto(c, Nc)
    ms = sl.modos(MZ); rs = [ctx.rad(m, NZ) for m in ms]
    offs = np.cumsum([0] + [r.dim for r in rs])
    return ctx, ms, rs, offs


def JZ_blocos(rs, s):
    return [(nk.pesos(r)[:, None]*sl.J_livre(r, s))/nk.pesos(r)[None, :] for r in rs]


def aplica_blocos(blocos, offs, X):
    out = np.empty_like(X)
    for j, M in enumerate(blocos):
        out[offs[j]:offs[j+1]] = M @ X[offs[j]:offs[j+1]]
    return out


def vetores_refA(g, rs, offs):
    """Vetores de deflacao EXPLICITOS de RefA na caixa Z (pesados): fase v = d_tau omega_A (autovalor
    lambda = 0 de L0) e gauge g = (Z_A + 1) omega_A (lambda = -mu_A). O certificado usa estes vetores; a
    diferenca para os autovetores exatos do operador verdadeiro (d_tau omega*, (Z* + 1) omega*) nas linhas de
    Z e ||P_Z d_tau (RefB + Corr)|| e ||P_Z (g* - g_A)|| (etapa E de S4), e fora de Z entra nas linhas de T."""
    from gauge_vector_L import g_modo
    n = int(offs[-1]); vf = np.zeros(n, complex); vg = np.zeros(n, complex)
    for j, r in enumerate(rs):
        w = nk.pesos(r)
        for c, ns in r.comps:
            a_ = sl.modo(g, c, r.m); full = np.zeros(max(int(ns.max()) + 1, len(a_)), complex); full[:len(a_)] = a_
            o0, o1 = r.off[c]
            vf[offs[j] + o0:offs[j] + o1] = (1j*r.m/2)*full[ns]
        gm, rr = g_modo(g, r.m, max(int(np.concatenate([ns for _, ns in r.comps]).max()) + 1, 1))
        assert rr.dim == r.dim
        vf[offs[j]:offs[j+1]] *= w; vg[offs[j]:offs[j+1]] = w*gm
    return [(0.0, vf), (-nk.MU, vg)]


def caudas_refA(g, MZ, NZ, Mg=41):
    """Normas (pesos de L) das partes de d_tau omega_A e (Z_A + 1) omega_A FORA da caixa Z (modos
    MZ <= |m| < Mg; RefA tem n <= 100 < NZ): entram nas linhas de T da deflacao verdadeira."""
    from gauge_vector_L import g_modo
    cv = 0.0; cg = 0.0
    for m in [m for m in sl.modos(Mg) if abs(m) >= MZ]:
        r = sl.Radial(m, max(NZ, 110)); w = nk.pesos(r)
        v = np.zeros(r.dim, complex)
        for c, ns in r.comps:
            a_ = sl.modo(g, c, m); full = np.zeros(max(int(ns.max()) + 1, len(a_)), complex); full[:len(a_)] = a_
            o0, o1 = r.off[c]; v[o0:o1] = (1j*m/2)*full[ns]
        gm, _ = g_modo(g, m, max(NZ, 110))
        cv += float(np.sum(np.abs(w*v)**2)); cg += float(np.sum(np.abs(w*gm)**2))
    return math.sqrt(cv)*(1 + 1e-9), math.sqrt(cg)*(1 + 1e-9)


def monta_A(Zs, Fs, deflacao, phi, vetores):
    """Matriz de referencia A = W L~_Z(0) W^-1 (deflacao de Brauer da fase e de mu), a partir dos dados de RT.
    Usada por cmd_schur e por scripts/reconstroi_A.py (nivel 2 da reproducao)."""
    MZ, NZ = map(int, Zs.split('x')); Mc, Nc = map(int, Fs.split('x'))
    ctx = nk.Contexto(0.0, Nc)
    dfl = [complex(x) for x in deflacao.split(',') if x.strip()]
    z = nk.nivel_Z(ctx, MZ, NZ, log=lambda *x, **k: None, sem_borda=True,
                   deflacao=(dfl if vetores == 'truncagem' else []), phi=phi)
    del z['Vb']
    if vetores == 'refA':                      # deflacao por vetores explicitos de RefA, phi = v/||v||^2
        vs = vetores_refA(ctx.g, z['rs'], z['offs'])
        X = np.stack([v for _, v in vs], 1)
        Phi = X @ np.linalg.inv(X.conj().T @ X) if phi == 'biortogonal' else X/np.sum(np.abs(X)**2, 0)
        z['defl'] = [(1.0 - lam, v, Phi[:, k]) for k, (lam, v) in enumerate(vs)]
        cauda_v, cauda_g = caudas_refA(ctx.g, MZ, NZ)
        extra_defl = dict(delta=[float((1.0 - lam).real) for lam, _ in vs],
                          norma_x=[float(np.linalg.norm(v)) for _, v in vs],
                          norma_phi=[float(np.linalg.norm(Phi[:, k])) for k in range(len(vs))],
                          cauda_x=[cauda_v, cauda_g])
    else:
        extra_defl = {}
    ms, offs, rs = z['ms'], z['offs'], z['rs']; n = z['n']
    A = np.zeros((n, n), complex)
    for j, (ri, Jb) in enumerate(zip(rs, JZ_blocos(rs, 0.0))):
        A[offs[j]:offs[j+1], offs[j]:offs[j+1]] += Jb
        for i, ro in enumerate(rs):
            if abs(ro.m - ri.m) <= nk.BAND_M:
                A[offs[i]:offs[i+1], offs[j]:offs[j+1]] += ctx.B(ro.m, ri.m, NZ, NZ)
    ctx._B.cache_clear()
    for dk, vw, uw in z['defl']:
        A += dk*np.outer(vw, uw.conj())
    return A, z, extra_defl


def cmd_schur(a):
    t0 = time.time()
    MZ, NZ = map(int, a.Z.split('x')); Mc, Nc = map(int, a.F.split('x'))
    A, z, extra_defl = monta_A(a.Z, a.F, a.deflacao, a.phi, a.vetores)
    n = z['n']
    log(t0, f'A montada (n = {n})')
    np.save(a.dir/'A.npy', A)
    T, U = schur(A, output='complex', overwrite_a=True, check_finite=False); del A
    log(t0, 'Schur')
    np.save(a.dir/'T.npy', T); np.save(a.dir/'U.npy', U)
    ev = -np.diag(T)
    dentro = [e for e in ev if -1e-9 < e.real < a.R and abs(e.imag) < 0.25]
    json.dump(dict(Z=[MZ, NZ], F=[Mc, Nc], deflacao=a.deflacao, phi=a.phi, vetores=a.vetores, n=int(n),
                   raizes_em_Omega=[[e.real, e.imag] for e in sorted(dentro, key=lambda z: z.real)],
                   defl=[dict(delta=[complex(d).real, complex(d).imag]) for d, _, _ in z['defl']], defl_refA=extra_defl),
              open(a.dir/'schur.json', 'w'), indent=1)
    np.save(a.dir/'defl.npy', np.array([np.concatenate([[d], v, u]) for d, v, u in z['defl']]))
    log(t0, 'raizes de L~_Z em Omega: ' + ', '.join(f'{e.real:.6f}{e.imag:+.1e}i' for e in sorted(dentro, key=lambda z: z.real)))


def blocos_T(ctx, ms, offs, MZ, NZ, Mc, Nc):
    """Por modo m de T acoplado a Z: (m, cT, Y~_m = P_Z (J + B) P_{T,m}, J_ZT,m)."""
    idx = {m: i for i, m in enumerate(ms)}; n = int(offs[-1])
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
        Yt += Jzt
        yield m, cT, Yt, Jzt


def cmd_fator(a):
    t0 = time.time()
    MZ, NZ = map(int, a.Z.split('x')); Mc, Nc = map(int, a.F.split('x'))
    c = a.c
    ctx, ms, rs, offs = layout(MZ, NZ, Nc, c); n = int(offs[-1])
    rng = np.random.default_rng(7)
    k = a.k
    Ysk = np.zeros((2*n, k + 20), complex)
    for m, cT, Yt, Jzt in blocos_T(ctx, ms, offs, MZ, NZ, Mc, Nc):   # passada 1: esboço S Omega
        Q = ctx.Q(m, Nc)[np.ix_(cT, cT)]
        Om = rng.standard_normal((len(cT), k + 20)) + 1j*rng.standard_normal((len(cT), k + 20))
        QO = Q @ Om
        Ysk[:n] += Yt @ QO; Ysk[n:] += Jzt @ QO
    ctx.Q.cache_clear(); ctx._B.cache_clear()
    Qb, _ = np.linalg.qr(Ysk); Qb = Qb[:, :k]; del Ysk
    log(t0, f'passada 1 (esboço, posto {k})')
    CC = np.zeros((k, k), complex); fro2 = 0.0; proj2 = 0.0
    for m, cT, Yt, Jzt in blocos_T(ctx, ms, offs, MZ, NZ, Mc, Nc):   # passada 2: C C^* e resto
        Q = ctx.Q(m, Nc)[np.ix_(cT, cT)]
        S = np.vstack([Yt @ Q, Jzt @ Q])
        Cm = Qb.conj().T @ S
        CC += Cm @ Cm.conj().T
        fro2 += float(np.sum(np.abs(S)**2)); proj2 += float(np.sum(np.abs(S - Qb @ Cm)**2))
    ctx.Q.cache_clear(); ctx._B.cache_clear()
    lam, W = np.linalg.eigh(CC)
    F = Qb @ (W*np.sqrt(np.maximum(lam, 0)))
    resto = math.sqrt(proj2)
    np.save(a.dir/f'F_{rotulo(c)}.npy', F)
    json.dump(dict(c=[c.real, c.imag], k=k, resto_F=resto, fro=math.sqrt(fro2), qTT=qTT(ctx, MZ, NZ, Mc, Nc)),
              open(a.dir/f'F_{rotulo(c)}.json', 'w'), indent=1)
    log(t0, f'fator em c = {c}: posto {k}, ||S||_F {math.sqrt(fro2):.3e}, resto_F {resto:.2e}')


def blocos_T2(ctx, ms, offs, MZ, NZ, Mc, Nc):
    """Por modo m de T acoplado a Z: (m, cT, P_Z B P_{T,m}, J_ZT,m restrito as linhas de Z do modo m, i, Q_TT,m(c))."""
    idx = {m: i for i, m in enumerate(ms)}; n = int(offs[-1])
    for m in range(-(min(MZ + nk.BAND_M, Mc) - 1), min(MZ + nk.BAND_M, Mc)):
        cT = nk.colsT(ctx, m, MZ, NZ, Nc)
        Yb = np.zeros((n, len(cT)), complex)
        for mo in ms:
            if abs(mo - m) <= nk.BAND_M:
                i = idx[mo]
                Yb[offs[i]:offs[i+1]] = ctx.B(mo, m, NZ, Nc)[:, cT]
        Jm = None; i = idx.get(m)
        if i is not None:
            r = ctx.rad(m, Nc); w = nk.pesos(r)
            Jw = (w[:, None]*sl.J_livre(r, 0.0))/w[None, :]
            Jm = Jw[np.ix_(nk.nZ_local(ctx, m, NZ, Nc), cT)]
        yield m, cT, Yb, Jm, i, ctx.Q(m, Nc)[np.ix_(cT, cT)]


def cmd_fator2(a):
    """Fator com a parte de J EXATA: S = [[(B + J)_ZT Q_TT], [J_ZT Q_TT]] = [[S_B + S_J], [S_J]], S_J = J_ZT Q_TT
    bloco-diagonal por modo, de posto <= ~8 por modo (as colunas de J_ZT nas linhas de Z sao constantes por
    componente e paridade): SVD por modo, S_J = F_J W^* com W de colunas ortonormais. S_B por range finder com
    iteracao de potencia: S_B = Q_B C_B + Resto_B, ||Resto_B||_F <= resto. Com K = [C_B; W^*][C_B; W^*]^* = L L^*,
    grava F = [[Q_B, F_J] L; [0, F_J] L]: para qualquer M1, M2, ||M1 (S_B + S_J) - M2 S_J|| <= ||[M1, M2](F)|| +
    ||M1|| resto (o mesmo formato de `fator`, e a avaliacao por ponto nao muda)."""
    t0 = time.time()
    MZ, NZ = map(int, a.Z.split('x')); Mc, Nc = map(int, a.F.split('x'))
    c = a.c
    ctx, ms, rs, offs = layout(MZ, NZ, Nc, c); n = int(offs[-1])
    rng = np.random.default_rng(7); k = a.k; kk = k + 20
    # parte J: SVD por modo (exata a menos de sigma relativos < 1e-14)
    FJ_cols = []; W_blocos = {}; eJ2 = 0.0
    Ysk = np.zeros((n, kk), complex)
    for m, cT, Yb, Jm, i, Q in blocos_T2(ctx, ms, offs, MZ, NZ, Mc, Nc):
        Om = rng.standard_normal((len(cT), kk)) + 1j*rng.standard_normal((len(cT), kk))
        Ysk += Yb @ (Q @ Om)
        if Jm is not None:
            X = Jm @ Q
            Uu, sv, Wh = np.linalg.svd(X, full_matrices=False)
            keep = sv > sv[0]*1e-14 if sv.size and sv[0] > 0 else np.zeros(0, bool)
            eJ2 += float(np.sum(sv[~keep]**2))
            cols = np.zeros((n, int(keep.sum())), complex)
            lz = nk.nZ_local(ctx, m, NZ, Nc)
            cols[offs[i]:offs[i+1]] = (Uu[:, keep]*sv[keep])
            FJ_cols.append(cols); W_blocos[m] = Wh[keep].conj().T      # |cT| x r_m, colunas ortonormais
    FJ = np.hstack(FJ_cols) if FJ_cols else np.zeros((n, 0), complex); rJ = FJ.shape[1]
    Qb, _ = np.linalg.qr(Ysk); del Ysk
    log(t0, f'passada 1: esboco de S_B (posto {k}), parte J exata de posto {rJ}')
    for it in range(a.potencia):                      # iteracao de potencia: Qb <- qr(S_B S_B^* Qb)
        Z = np.zeros((n, Qb.shape[1]), complex)
        for m, cT, Yb, Jm, i, Q in blocos_T2(ctx, ms, offs, MZ, NZ, Mc, Nc):
            SB = Yb @ Q
            Z += SB @ (SB.conj().T @ Qb)
        Qb, _ = np.linalg.qr(Z); del Z
        log(t0, f'iteracao de potencia {it + 1}')
    Qb = Qb[:, :k]
    CC = np.zeros((k, k), complex); CW = np.zeros((k, rJ), complex); fro2 = 0.0; proj2 = 0.0
    oj = np.cumsum([0] + [W_blocos[m].shape[1] for m in W_blocos]); posJ = {m: (oj[t], oj[t+1]) for t, m in enumerate(W_blocos)}
    for m, cT, Yb, Jm, i, Q in blocos_T2(ctx, ms, offs, MZ, NZ, Mc, Nc):
        SB = Yb @ Q
        Cm = Qb.conj().T @ SB
        CC += Cm @ Cm.conj().T
        if m in W_blocos:
            j0, j1 = posJ[m]; CW[:, j0:j1] += Cm @ W_blocos[m]
        R_ = SB - Qb @ Cm
        fro2 += float(np.sum(np.abs(SB)**2)); proj2 += float(np.sum(np.abs(R_)**2))
    ctx.Q.cache_clear(); ctx._B.cache_clear()
    K = np.block([[CC, CW], [CW.conj().T, np.eye(rJ)]])
    lam, V = np.linalg.eigh(K)
    L = V*np.sqrt(np.maximum(lam, 0))                # K = L L^*
    Ft = np.hstack([Qb, FJ]) @ L; Fb = np.hstack([np.zeros((n, k), complex), FJ]) @ L
    resto = math.sqrt(proj2)*(1 + 1e-6) + math.sqrt(eJ2)
    np.save(a.dir/f'F_{rotulo(c)}.npy', np.vstack([Ft, Fb]))
    json.dump(dict(c=[c.real, c.imag], k=k, rJ=rJ, potencia=a.potencia, resto_F=resto, fro=math.sqrt(fro2), qTT=qTT(ctx, MZ, NZ, Mc, Nc),
                   tipo='fator2 (parte J exata)'), open(a.dir/f'F_{rotulo(c)}.json', 'w'), indent=1)
    log(t0, f'fator2 em c = {c}: posto {k} + J {rJ}, ||S_B||_F {math.sqrt(fro2):.3e}, resto_F(S_B) {resto:.2e}')


def qTT(ctx, MZ, NZ, Mc, Nc):
    q = 0.0
    for m in range(-(Mc - 1), Mc):
        cT = nk.colsT(ctx, m, MZ, NZ, Nc)
        q = max(q, float(np.linalg.norm(ctx.Q(m, Nc)[np.ix_(cT, cT)], 2)))
    ctx.Q.cache_clear()
    return q


def rotulo(c):
    return f'{c.real:+.4f}{c.imag:+.4f}i'.replace('.', '_')


def cmd_pontos(a):
    t0 = time.time()
    MZ, NZ = map(int, a.Z.split('x')); Mc, Nc = map(int, a.F.split('x'))
    c = a.c
    _, ms, rs, offs = layout(MZ, NZ, Nc); n = int(offs[-1])
    T = np.load(a.dir/'T.npy'); U = np.load(a.dir/'U.npy')
    F = np.load(a.dir/f'F_{rotulo(c)}.npy'); info = json.load(open(a.dir/f'F_{rotulo(c)}.json'))
    Ft, Fb = F[:n], F[n:]
    UFt = U.conj().T @ Ft
    Jc = JZ_blocos(rs, c)
    I_n = np.eye(n)
    pts = [complex(p) for p in a.pontos.split(',')]
    out = []
    for p in pts:
        Tp = T + p*I_n
        X1 = solve_triangular(Tp, UFt, check_finite=False)          # (T + p)^-1 U^* F_t
        R1 = U @ X1                                                  # R(p) F_t
        R2 = U @ solve_triangular(Tp, X1, check_finite=False)        # R(p)^2 F_t
        del Tp
        ctxp = nk.Contexto(p, NZ)
        Qp = [ctxp.Q(r.m, NZ) for r in rs]
        QF = aplica_blocos(Qp, offs, Fb); Q2F = aplica_blocos(Qp, offs, QF)
        E0 = aplica_blocos(Jc, offs, R1 - QF); E1 = aplica_blocos(Jc, offs, R2 - Q2F)
        a0 = float(np.linalg.norm(E0, 2)); b1 = float(np.linalg.norm(E1, 2))
        # ||R(p)|| e ||J_Z(c) R(p)|| por potência (estimativa) para o resto e a segunda ordem
        x = np.random.default_rng(1).standard_normal(n) + 0j
        for _ in range(6):
            y = U @ solve_triangular(T + p*I_n, U.conj().T @ x, check_finite=False); res = np.linalg.norm(y)
            x = U @ solve_triangular((T + p*I_n).conj().T, U.conj().T @ (y/res), lower=True, check_finite=False)
            x /= np.linalg.norm(x)
        fator_q = 1/(1 - abs(p - c)*info['qTT'])
        a_esc = (a0 + 0.0)*fator_q                                  # resto_F entra na versão rigorosa
        out.append(dict(p=[p.real, p.imag], a_esc=a_esc, b1=b1, res=float(res)))
        log(t0, f'  p = {p.real:.4f}{p.imag:+.4f}i: Z<-T escalado {a_esc:.4f} (b1 {b1:.3f}, ||R|| {res:.1f};'
                f' raio 1a ordem p/ +0,02: {0.02/max(b1, 1e-9):.2e}, raio resolvente 0,3/||R||: {0.3/res:.2e})')
    json.dump(dict(c=[c.real, c.imag], pontos=out), open(a.dir/f'pontos_{rotulo(c)}.json', 'w'), indent=1)


def cmd_certF(a):
    """Certificado do fator F de um centro c, sem identidades exatas sobre quantidades de ponto flutuante.
    S = [[(B + J)_ZT Q0_TT(c)], [J_ZT Q0_TT(c)]] EXATA (B = os blocos calculados, tomados como dados: a
    diferenca para o B verdadeiro esta em dB na montagem; Q0_TT(c) exata, com ||Q0 - Q|| <= dq pelo residuo
    do modo). Para QUALQUER matriz Phi: ||X S|| <= ||X F|| ||Phi|| + ||X|| ||S - F Phi||. Phi_m = V f(Sigma) U^* S_m
    (F ~ U Sigma V^* em ponto flutuante; f = sigma/(sigma^2 + lambda^2), Tikhonov; Phi e so um dado), por modo m
    de T acoplado a Z (colunas disjuntas), para varios lambda na mesma passada:
      ||Phi||^2 <= ||sum_m Phi_m Phi_m^*|| + arredondamento,  rho = ||S - F Phi||_F (com o erro de S_calc),
    e q_TT = max_m ||Q0_TT,m(c)|| certificado nos mesmos modos (as colunas dos demais modos de S sao nulas).
    Cada par (||Phi||, rho) e um certificado; a montagem usa o melhor por disco.
    Arredondamentos por Frobenius: |fl(X Y) - X Y| <= gamma_k |X||Y|, || |X||Y| ||_F <= ||X||_F ||Y||_F."""
    from verified_norms import U as UNIT, cota_norma, gamma
    t0 = time.time()
    MZ, NZ = map(int, a.Z.split('x')); Mc, Nc = map(int, a.F.split('x'))
    c = a.c
    ctx, ms, rs, offs = layout(MZ, NZ, Nc, c); n = int(offs[-1])
    idx = {m: i for i, m in enumerate(ms)}
    import hashlib
    with open(a.dir/f'F_{rotulo(c)}.npy', 'rb') as fh:
        shaF = hashlib.sha256(fh.read()).hexdigest()
    F = np.load(a.dir/f'F_{rotulo(c)}.npy'); r = F.shape[1]
    fFt = float(np.linalg.norm(F[:n])); fFb = float(np.linalg.norm(F[n:]))
    Uf, sv, Vh = np.linalg.svd(F, full_matrices=False)
    Uh = Uf.conj().T; del Uf                          # r x 2n
    V = Vh.conj().T; del Vh
    lams = [float(x)*sv[0] for x in a.tikhonov_rel.split(',')]
    filtros = [sv/(sv**2 + l**2) if l > 0 else np.where(sv > sv[0]*1e-12, 1/sv, 0.0) for l in lams]
    log(t0, f'F: {F.shape}, sigma {sv[0]:.3e} .. {sv[-1]:.3e}; lambdas {", ".join(f"{l:.2e}" for l in lams)}')
    nl = len(lams)
    G = [np.zeros((r, r), complex) for _ in lams]; eG = [0.0]*nl; d2 = [0.0]*nl; e2 = [0.0]*nl; eS2 = 0.0; qTT = 0.0
    for m in range(-(min(MZ + nk.BAND_M, Mc) - 1), min(MZ + nk.BAND_M, Mc)):
        cT = nk.colsT(ctx, m, MZ, NZ, Nc)
        rr = ctx.rad(m, Nc)
        w = np.sqrt(np.concatenate([cf_omega(ns) for _, ns in rr.comps]))
        Jw = (w[:, None]*sl.J_livre(rr, c))/w[None, :]
        Qw = ctx.Q(m, Nc); ctx.Q.cache_clear()       # cada modo e usado uma vez (memoria)
        E = np.eye(len(Qw)) - Jw @ Qw
        fE = float(np.linalg.norm(E)); fQ = float(np.linalg.norm(Qw))
        fJQ = float(np.linalg.norm(np.abs(Jw) @ np.abs(Qw)))*(1 + 1e-6)   # elemento a elemento: produto e pesos
        e = fE*(1 + 1e-6) + (gamma(len(Qw) + 8) + 4*UNIT)*fJQ + UNIT*fE
        assert e < 0.5
        dq = fQ*(1 + 1e-6)*e/(1 - e)                 # ||Q0_m - Qw|| (espectral <= Frobenius)
        del E, Jw
        QT = Qw[np.ix_(cT, cT)]
        qTT = max(qTT, cota_norma(QT, dq))
        Yb = np.zeros((n, len(cT)), complex)
        for mo in ms:
            if abs(mo - m) <= nk.BAND_M:
                i = idx[mo]
                Yb[offs[i]:offs[i+1]] = ctx.B(mo, m, NZ, Nc)[:, cT]
        Jm = None
        if m in idx:
            i = idx[m]; wz = nk.pesos(rr)
            Jm = ((wz[:, None]*sl.J_livre(rr, 0.0))/wz[None, :])[np.ix_(nk.nZ_local(ctx, m, NZ, Nc), cT)]
            Yb[offs[i]:offs[i+1]] += Jm                  # (B + J)_ZT
        fY = float(np.linalg.norm(Yb)); aQT = np.abs(QT)
        top = Yb @ QT
        # Q exata (dq), produto e soma B + J / pesos (elemento a elemento: gamma |Y||Q| + 3u |Y||Q|)
        em = fY*dq + (gamma(len(cT)) + 3*UNIT)*float(np.linalg.norm(np.abs(Yb) @ aQT))*(1 + 1e-6)
        del Yb
        Wm = Uh[:, :n] @ top
        if Jm is not None:
            bot = Jm @ QT; fJm = float(np.linalg.norm(Jm))
            em = math.hypot(em, fJm*dq + (gamma(len(cT)) + 2*UNIT)*float(np.linalg.norm(np.abs(Jm) @ aQT))*(1 + 1e-6))
            rows = np.arange(offs[i], offs[i+1]) + n
            Wm += Uh[:, rows] @ bot
        eS2 += em**2
        nTop = float(np.linalg.norm(top))
        for t, f in enumerate(filtros):
            Phi = V @ (f[:, None]*Wm)
            fPhi = float(np.linalg.norm(Phi))
            FP = F[:n] @ Phi; FP -= top                  # = -(S - F Phi), linhas de cima
            dm2 = float(np.linalg.norm(FP))**2
            er = gamma(r)*fFt*fPhi + UNIT*(nTop + float(np.linalg.norm(FP)) + nTop)
            FP = F[n:] @ Phi
            if Jm is not None:
                FP[offs[i]:offs[i+1]] -= bot
            dm2 += float(np.linalg.norm(FP))**2
            er = math.hypot(er, gamma(r)*fFb*fPhi + UNIT*2*float(np.linalg.norm(FP)))
            del FP
            d2[t] += dm2; e2[t] += er**2
            G[t] += Phi @ Phi.conj().T; eG[t] += gamma(len(cT))*fPhi**2
        del top, Wm
        if m % 8 == 0:
            log(t0, f'  modo {m}: rho parcial ' + ', '.join(f'{math.sqrt(x):.3e}' for x in d2))
    ctx.Q.cache_clear(); ctx._B.cache_clear()
    pares = []
    for t, l in enumerate(lams):
        phi = math.sqrt(cota_norma(G[t], eG[t]))*(1 + 1e-9)
        rho = math.sqrt(d2[t])*(1 + 1e-6) + math.sqrt(e2[t]) + math.sqrt(eS2)
        pares.append(dict(lam=l, phi=phi, rho=rho))
        log(t0, f'certF em c = {c}, lambda {l:.2e}: ||Phi|| <= {phi:.9f}, rho <= {rho:.4e}')
    fat = json.load(open(a.dir/f'F_{rotulo(c)}.json'))
    json.dump(dict(c=[c.real, c.imag], sha256_F=shaF, pares=pares, qTT=qTT, posto=r, erro_S=math.sqrt(eS2), sigma=[float(sv[0]), float(sv[-1])],
                   resto_F_antigo=fat['resto_F'], qTT_antigo=fat['qTT']),
              open(a.dir/f'certF_{rotulo(c)}.json', 'w'), indent=1)
    log(t0, f'q_TT <= {qTT:.6f} (fator: {fat["qTT"]:.6f}); resto_F do fator {fat["resto_F"]:.3e}; erro de S_calc {math.sqrt(eS2):.1e}')


def cf_omega(ns):
    import closed_form as cf
    return cf.omega(ns, nk.K2)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('cmd', choices=['schur', 'fator', 'fator2', 'pontos', 'certF'])
    ap.add_argument('--Z', default='32x128'); ap.add_argument('--F', default='200x400')
    ap.add_argument('--dir', type=Path, default=Path('build/rouche_L'))
    ap.add_argument('--deflacao', default='0,0.16831')
    ap.add_argument('--phi', default='biortogonal', choices=['biortogonal', 'direito', 'esquerdo'])
    ap.add_argument('--vetores', default='truncagem', choices=['truncagem', 'refA'])
    ap.add_argument('--R', type=float, default=3.0)
    ap.add_argument('--c', type=complex, default=0j)
    ap.add_argument('--k', type=int, default=800)
    ap.add_argument('--potencia', type=int, default=1)
    ap.add_argument('--pontos', default='0')
    ap.add_argument('--tikhonov-rel', dest='tikhonov_rel', default='1e-7,3e-7,1e-6',
                    help='certF: lambdas de Tikhonov relativos a sigma_max (0 = pseudo-inversa)')
    a = ap.parse_args()
    a.dir.mkdir(parents=True, exist_ok=True)
    dict(schur=cmd_schur, fator=cmd_fator, fator2=cmd_fator2, pontos=cmd_pontos, certF=cmd_certF)[a.cmd](a)


if __name__ == '__main__':
    main()
