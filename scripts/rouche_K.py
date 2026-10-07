#!/usr/bin/env python3
"""S3b, K: exatamente uma raiz de K (simples) em |s - c| < rho (c = 0,401, rho = 0,05, o
disco que a cobertura de S3a deixou de fora). Duas partes.

(1) ROUCHE DE OPERADORES (Gohberg-Sigal) na circunferencia |s - c| = rho, subcomando circ.
    P FIXO: P = Q0(c) nos tres niveis, H(s) = K(s) P = I + (B + s - c) Q0(c). Referencia
        A(s) = diag(Hc_ZZ(s), I, I, I),   Hc_ZZ(s) = I + P_Z (B_t + s - c) Q0(c) P_Z,
    com o fundo TRUNCADO B_t (dados de RT em aritmetica exata) e mu = mu_RefA (diadico).
    Tile j da circunferencia (centro p_j = c + dV_j, raio rZ): V_j ~ Hc_ZZ(p_j)^-1 e a matriz
    de Perron N de tile_perron (parte_Z com dV = dV_j; T1 e T2 uma vez, em c, raio R). Com
    e = N_ZZ >= ||I - V Hc_ZZ(s)|| < 1 (rho de parte_Z cobre Hc: dH >= dB ||Q|| + arredondamento):
        ||Hc_ZZ^-1 H_Zk|| <= N_Zk/(1 - e)                                    (k != Z),
        ||Hc_ZZ^-1 (H_ZZ - Hc_ZZ)|| <= ||V|| (dB ||Q_ZZ|| + eps_mu(c))/(1 - e),
    e as linhas de T1, T2 e far sao as de N (a referencia e I nelas). Se N' v <= theta v com
    theta < 1 (racionais), ||A^-1 (H - A)|| < 1 na norma max_i ||x_i||/v_i. H = I + compacto,
    A = I + posto finito: H tem no disco tantos zeros (com multiplicidade) quanto A, isto e,
    quanto det Hc_ZZ(s) = det Kc_Z(s)/det J_Z(c) -- os s com Kc_Z(s) = J_Z(s) + P_Z B_t P_Z
    singular (K_livre triangular superior em n: Q0 P_Z = P_Z Q0 P_Z e Q_ZZ = J_Z(c)^-1).
    A norma com sinal NAO e simetrica pela conjugacao: a circunferencia inteira e coberta.

(2) CONTAGEM FINITA, subcomando contagem: Kc_Z tem exatamente um autovalor (simples) no disco.
    Bordejado M(s) = [[Kc_Z(s), r], [l^H, 0]] (r, l autovetores aproximados, fixos);
    t(s) = e^T M(s)^-1 e = det Kc_Z(s)/det M(s). Com V_b ~ M(c)^-1, E(s) = I - V_b M(s) =
    E_c - (s - c) V_b D (D = diag(I, 0)), ||E_c|| <= eps_c, ||V_b D|| <= beta, q = eps_c +
    rho beta < 1: M invertivel no disco fechado e
        t(s) = sum_{j <= K} (-(s - c))^j tau_j + R,  tau_j = e^T (V_b D)^j V_b e,
        |R| <= ||V_b e|| (sum_{j <= K} (q^j - (rho beta)^j) + q^(K+1)/(1 - q)).
    Rouche escalar contra g(s) = tau_0 - (s - c) tau_1 (um zero): basta
        sum_{j >= 2} rho^j |tau_j| + sum_j rho^j dtau_j + |R| < rho |tau_1| - |tau_0|.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Fr
import json
import math
from pathlib import Path
import sys
import time

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import tile_perron as tp
from signed_operator import Jrad
from verified_norms import cota_norma, gamma, U

ARRED = 1 + Fr(1, 10**9)


def perron_q(N):
    """v > 0 (autovetor de Perron em ponto flutuante) e max_i (N v)_i/v_i em racionais exatos."""
    M = np.array([[float(x) for x in linha] for linha in N])
    lam, vec = np.linalg.eig(M + 1e-9)
    v = np.maximum(np.abs(vec[:, np.argmax(np.abs(lam))]), 1e-9)
    v = v/v.max()
    vq = [Fr(x) for x in v.tolist()]
    th = max(sum(N[i][k]*vq[k] for k in range(len(vq)))/vq[i] for i in range(len(vq)))
    return v, th


def centros(c, rho, n):
    """Centros dV_j = rho e^(2 pi i j/n) (ponto flutuante: o centro EXATO do tile e c + dV_j),
    raio rZ = 2 rho sin(pi/(2n)) (1 + 1e-6) e raio R dos niveis T, com folga."""
    dV = [complex(rho*math.cos(2*math.pi*j/n), rho*math.sin(2*math.pi*j/n)) for j in range(n)]
    rZ = 2*rho*math.sin(math.pi/(2*n))*(1 + 1e-6)
    R = (rho + rZ)*(1 + 1e-6)
    return dV, rZ, R


def confere_cobertura(rho, n, dV, rZ, R):
    """Rigoroso (Arb + racionais): cada arco phi em [pi(2j-1)/n, pi(2j+1)/n] da circunferencia
    |s - c| = rho esta no disco |s - c - dV_j| <= rZ, e |dV_j| + rZ <= R. A distancia
    |rho e^(i phi) - dV_j|^2 = rho^2 + |dV_j|^2 - 2 rho |dV_j| cos(phi - arg dV_j) e
    quase-convexa em phi num arco de comprimento < 2 pi: o maximo esta nas pontas."""
    import flint
    flint.ctx.prec = 200
    rho_a, rZ_a = flint.arb(rho), flint.arb(rZ)
    for j, d in enumerate(dV):
        for k in (2*j - 1, 2*j + 1):
            f = flint.arb(flint.fmpq(k, n))
            x = rho_a*f.cos_pi() - flint.arb(d.real)
            y = rho_a*f.sin_pi() - flint.arb(d.imag)
            if not (x*x + y*y < rZ_a*rZ_a):
                return False, f'arco {j}, ponta {k}/{n}'
        if Fr(R) < Fr(rZ) or (Fr(R) - Fr(rZ))**2 < Fr(d.real)**2 + Fr(d.imag)**2:
            return False, f'tile {j} fora do disco R dos niveis T'
    return True, 'ok'


def circ(a):
    tp.RIG = a.rigoroso
    tp.BAND_M, tp.BAND_N = map(int, a.banda.split(','))
    cx = lambda t: tuple(map(int, t.split('x')))
    Z, Gp, F = cx(a.Z), cx(a.Gp), cx(a.F)
    c = complex(a.c)
    dV, rZ, R = centros(a.c, a.rho, a.n)
    ok_cob, msg = confere_cobertura(a.rho, a.n, dV, rZ, R)
    print(f'circunferencia |s - {a.c}| = {a.rho}: {a.n} tiles de raio rZ = {rZ:.6f}; niveis T em c, '
          f'raio R = {R:.6f}; Z = {Z}, G+ = {Gp}, F = {F}, rigoroso = {tp.RIG}', flush=True)
    print(f'  cobertura da circunferencia (Arb): {msg}', flush=True)
    ctx = tp.Ctx(Z, Gp, F, a.D)
    print(f'  fundo: resto da banda {ctx.resto_banda:.1e}, dB {ctx.dB:.1e}', flush=True)
    q1, HZ1, QZ1 = tp.com_cache(a.cache1, c, ctx, lambda: tp.parte_T1(ctx, c), ('HZ1', 'QZ1'))
    q2, SZ, SQZ = tp.com_cache(a.cache2, c, ctx, lambda: tp.parte_T2(ctx, c), ('SZ', 'SQZ'))
    e_mu = Fr(tp.eps_mu(c))*ARRED
    js = range(a.n) if a.so is None else [int(x) for x in a.so.split(',')]
    saida = dict(c=a.c, rho=a.rho, n=a.n, rZ=rZ, R=R, Z=Z, Gp=Gp, F=F, D=a.D, banda=a.banda,
                 rigoroso=tp.RIG, cobertura=msg, dB=ctx.dB, eps_mu=float(e_mu), tiles=[])
    pior = Fr(0)
    for j in js:
        t0 = time.time()
        qZ, V = tp.parte_Z(ctx, c, dV=dV[j])
        N0, NZm, N1m, N2m, qa = tp.montar(ctx, qZ, V, q1, HZ1, QZ1, q2, SZ, SQZ, rZ, R, R)
        del V
        fr, fR = Fr(rZ), Fr(R)
        N = [[ARRED*(Fr(N0[i, k]) + fr*Fr(NZm[i, k]) + fR*(Fr(N1m[i, k]) + Fr(N2m[i, k])))
              for k in range(4)] for i in range(4)]
        nV = Fr(qZ['V'])
        for k in range(4):
            N[0][k] += nV*e_mu
            for i in (1, 2, 3):
                N[i][k] += e_mu
        v, th = perron_q(N)
        den = 1 - N[0][0]
        assert den > 0, f'tile {j}: N_ZZ = {float(N[0][0]):.3f} >= 1'
        Nr = [linha[:] for linha in N]
        Nr[0] = [nV*(Fr(ctx.dB)*Fr(qZ['nQZ']) + e_mu)*ARRED/den] + [N[0][k]/den for k in (1, 2, 3)]
        vr, thr = perron_q(Nr)
        pior = max(pior, thr)
        print(f'  tile {j:2d} (p = {c + dV[j]:.4f}): ||V|| {qZ["V"]:.2f}, N_ZZ {float(N[0][0]):.4f},'
              f' Perron do tile {float(th):.4f}, ROUCHE {float(thr):.4f}'
              + ('' if thr < 1 else '  <-- FALHA') + f'  [{time.time()-t0:.0f}s]', flush=True)
        saida['tiles'].append(dict(j=j, dV=[dV[j].real, dV[j].imag], V=qZ['V'], nQZ=qZ['nQZ'],
                                   rho=qZ['rho'], VQ_ZZ=qZ['VQ_ZZ'],
                                   N=[[float(x) for x in l] for l in N],
                                   N_rouche=[[str(x) for x in l] for l in Nr],
                                   v=vr.tolist(), theta_tile=float(th), theta_rouche=float(thr),
                                   ok=bool(thr < 1)))
    tudo = ok_cob and all(t['ok'] for t in saida['tiles']) and len(saida['tiles']) == a.n
    saida['pior_theta_rouche'] = float(pior); saida['certificado'] = bool(tudo and tp.RIG)
    print(f'  pior theta de Rouche {float(pior):.6f};  ' + (
        'CIRCUNFERENCIA CERTIFICADA: zeros de H no disco = autovalores de Kc_Z no disco'
        if saida['certificado'] else ('todos os tiles fecham (estimativa)' if tudo else 'NAO fecha')))
    saida['T1'] = {k: (str(x) if isinstance(x, (complex, tuple)) else x) for k, x in q1.items()}
    saida['T2'] = {k: (str(x) if isinstance(x, (complex, tuple)) else x) for k, x in q2.items()}
    if a.saida:
        a.saida.write_text(json.dumps(saida, indent=1, default=str) + '\n')


def matriz_Z(ctx, s):
    """Kc_Z(s) pesada (pesos radiais e kappa1^d de tile_perron) e cota do erro de montagem
    ||fl(Kc_Z) - Kc_Z|| (arredondamento de B_t: dB o majora; pesos e kappa1^d: <= 25u por
    entrada; J: <= 5u por entrada)."""
    (MZ, NZ) = ctx.Z
    mZ, msZ, kZ = ctx.mZ, ctx.msZ, ctx.kZ
    n = len(msZ)*kZ
    K = np.zeros((n, n), complex)
    Bc = {}
    fB2 = 0.0; fJ2 = 0.0
    for j, mi in enumerate(msZ):
        for i, mo in enumerate(msZ):
            d = mo - mi
            if abs(d) <= tp.BAND_M:
                if d not in Bc:
                    Bc[d] = tp.Bw(ctx.g, d, mZ, mZ)
                K[i*kZ:(i+1)*kZ, j*kZ:(j+1)*kZ] += Bc[d]
                fB2 += tp.fro(Bc[d])**2
        Jw = (mZ.w[:, None]*Jrad(mi, s, NZ))/mZ.w[None, :]
        K[j*kZ:(j+1)*kZ, j*kZ:(j+1)*kZ] += Jw
        fJ2 += tp.fro(Jw)**2
    dK = ctx.dB + 25*U*math.sqrt(fB2) + 5*U*math.sqrt(fJ2) + U*tp.fro(K)
    return K, dK


def contagem(a):
    import scipy.linalg as sla
    tp.RIG = True
    tp.BAND_M, tp.BAND_N = map(int, a.banda.split(','))
    cx = lambda t: tuple(map(int, t.split('x')))
    Z = cx(a.Z); Mz, Nz = Z
    ctx = tp.Ctx(Z, (Mz, Nz), (Mz + tp.BAND_M, Nz + tp.BAND_N + 1), 12)
    c, rho = complex(a.c), a.rho
    t0 = time.time()
    K, dK = matriz_Z(ctx, c)
    n = len(K)
    w, vl, vr = sla.eig(K, left=True, right=True)
    lam = c - w                                   # Kc_Z(s) = Kc_Z(c) + (s - c) I
    dentro = np.where(np.abs(lam - c) < rho)[0]
    print(f'Kc_Z em Z = {Z} (dim {n}), s em |s - {a.c}| < {rho}: autovalores em ponto flutuante '
          + ', '.join(f'{lam[i]:.10f}' for i in dentro) + f'  [{time.time()-t0:.0f}s]', flush=True)
    i0 = int(np.argmin(np.abs(w)))
    r = vr[:, i0]/np.linalg.norm(vr[:, i0]); l = vl[:, i0]/np.linalg.norm(vl[:, i0])
    del vl, vr
    M = np.zeros((n + 1, n + 1), complex)
    M[:n, :n] = K; M[:n, n] = r; M[n, :n] = l.conj()
    Vb = np.linalg.inv(M)
    X = np.eye(n + 1) - Vb @ M
    nVb = cota_norma(Vb)
    eps_c = cota_norma(X, gamma(n + 1)*tp.fro(Vb)*tp.fro(M) + U*tp.fro(X)) + nVb*dK
    del X
    VD = Vb[:, :n]
    beta = cota_norma(VD)
    nVbe = cota_norma(Vb[:, n:])
    qq = eps_c + rho*beta
    print(f'  ||V_b|| {nVb:.4f}, eps_c {eps_c:.2e}, beta = ||V_b D|| {beta:.4f}, q = eps_c + rho beta'
          f' {qq:.4f}, ||V_b e|| {nVbe:.4f}  (erro de montagem de Kc_Z {dK:.1e})', flush=True)
    assert qq < 1
    gVD = gamma(n)*tp.fro(VD)
    wk = Vb[:, n].copy(); err = 0.0
    tau = [complex(wk[n])]; dtau = [0.0]
    for k in range(1, a.K + 1):
        nw = float(np.linalg.norm(wk))*(1 + 1e-6)
        wk = VD @ wk[:n]
        err = beta*err + gVD*nw
        tau.append(complex(wk[n])); dtau.append(err)
    Q = lambda x: Fr(float(x))
    # ||E^j - (pura)^j|| <= (eps_c + x)^j - x^j, crescente em x: x = rho beta (cota superior), exato
    frho = Q(rho); x = frho*Q(beta); fq = Q(eps_c) + x
    assert fq < 1
    R1 = sum(fq**j - x**j for j in range(a.K + 1))
    R2 = fq**(a.K + 1)/(1 - fq)
    resto = Q(nVbe)*(R1 + R2)
    # |tau| em ponto flutuante: abs() com erro relativo <= 2u; inflar e deflacionar
    up = lambda z: Q(abs(z))*(1 + Fr(1, 10**12))
    lo = lambda z: Q(abs(z))*(1 - Fr(1, 10**12))
    esq = sum(frho**j*up(tau[j]) for j in range(2, a.K + 1)) + sum(frho**j*Q(dtau[j]) for j in range(a.K + 1)) + resto
    dir_ = frho*lo(tau[1]) - up(tau[0])
    raiz = c + tau[0]/tau[1]
    print(f'  tau_0 {tau[0]:.3e}, tau_1 {tau[1]:.6f}, |tau_2..| ' + ' '.join(f'{abs(x):.1e}' for x in tau[2:6])
          + f'; resto {float(resto):.1e}', flush=True)
    ok = esq < dir_
    print(f'  Rouche escalar em |s - c| = {rho}: {float(esq):.3e} < {float(dir_):.6f}?  -> '
          + (f'EXATAMENTE UM autovalor (simples) de Kc_Z no disco, em ~{raiz:.10f}' if ok else 'NAO certificado'))
    if a.saida:
        a.saida.write_text(json.dumps(dict(
            Z=Z, banda=a.banda, c=a.c, rho=rho, dim=n, autovalores=[str(lam[i]) for i in dentro],
            dK=dK, nVb=nVb, eps_c=eps_c, beta=beta, nVbe=nVbe, q=qq, K=a.K,
            tau=[str(x) for x in tau], dtau=dtau, resto=float(resto), esquerda=float(esq),
            direita=float(dir_), raiz_aprox=str(raiz), certificado=bool(ok)), indent=1) + '\n')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('circ', help='tiles da circunferencia (Rouche de operadores)')
    p.add_argument('--n', type=int, default=24)
    p.add_argument('--Gp', default='40x160'); p.add_argument('--F', default='200x400')
    p.add_argument('--D', type=int, default=12)
    p.add_argument('--rigoroso', action='store_true')
    p.add_argument('--cache1', type=Path); p.add_argument('--cache2', type=Path)
    p.add_argument('--so', default=None, help='so estes tiles (j1,j2,...), para teste')
    p2 = sub.add_parser('contagem', help='autovalores de Kc_Z no disco (bordejamento, rigoroso)')
    p2.add_argument('--K', type=int, default=12, help='termos de Taylor de t')
    for q in (p, p2):
        q.add_argument('--c', type=float, default=0.401); q.add_argument('--rho', type=float, default=0.05)
        q.add_argument('--Z', default='20x80'); q.add_argument('--banda', default='20,40')
        q.add_argument('--saida', type=Path)
    a = ap.parse_args()
    (circ if a.cmd == 'circ' else contagem)(a)


if __name__ == '__main__':
    main()
