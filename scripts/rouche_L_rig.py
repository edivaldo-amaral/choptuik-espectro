#!/usr/bin/env python3
"""Nível Z do Rouché de L no contorno, RIGOROSO, ponto a ponto (C1_REAVALIACAO.md §7).

Dados (de rouche_L_Z.py, em --dir): A = W L~_Z(0) W^-1 (float, exata como dado; a diferença para o
operador verdadeiro entra em dA), a decomposição de Schur A ~ U T U^* (float; só serve de solver) e o
fator F = [F_t; F_b] do centro c, tomado como dado: a relação com S = [[Y~ Q0_TT(c)], [J_ZT Q0_TT(c)]]
(||X S|| <= ||X F|| ||Phi|| + ||X|| rho) é certificada à parte (rouche_L_Z certF).

Por ponto p (subcomando `ponto`):
  t_p >= ||(A + p)^-1||, por Loewner (verified_inverse.cota_inv);
  X_k ~ R(p)^k F_t (k = 1, 2, 3) pelo Schur; Res_k = (A + p) X_k - X_{k-1} com cota de arredondamento
    elemento a elemento (|fl(A X) - A X| <= gamma_n |A||X|); então Y_k = R(p)^k F_t exato satisfaz
    ||Y_k - X_k|| <= t_p ||Res_k|| + t_p ||Y_{k-1} - X_{k-1}||;
  a0, b1, b2 = cotas de ||J_Z(c)[Y_{j+1} - Q_ZZ(p)^{j+1} F_b]|| (j = 0, 1, 2) e ||Y_3||, ||Q_ZZ(p)^3 F_b||.
Versao 2 (revisao independente de 01/10, L4 e L5): grava sha256 de A e de F; t_p corrigido pelo arredondamento
de fl(A + p); J_Z(0) e J_Z(c) calculadas em ponto flutuante (4u elemento a elemento) em ||B^|| e em a0, b1, b2;
arredondamento de Q_ZZ^k F_b e de J_Z(c) D elemento a elemento (gamma_dim || |M| |X| ||_F).
Para |s - p| = |e|: pela identidade do resolvente (a mesma para R e para Q_ZZ),
  ||J_Z(c)[R(s) F_t - Q_ZZ(s) F_b]|| <= a0 + |e| b1 + |e|^2 b2
                                       + |e|^3 (tau(s) ||Y_3|| + ||J_Z(c) Q_ZZ(s)|| ||Q_ZZ(p)^3 F_b||),
  tau(s) = ||J_Z(c) R(s)|| <= 1 + (|c - s| + ||B^||) t(s),  t(s) <= t_p/(1 - |e| t_p),  B^ = A - J_Z(0)."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys
import time

import numpy as np
from scipy.linalg import solve_triangular

sys.path.insert(0, str(Path(__file__).resolve().parent))
import nk_L as nk
from nk_L2_Z import Q_cert
from rouche_L_Z import JZ_blocos, aplica_blocos, layout, rotulo
from verified_inverse import cota_inv
from verified_norms import U as UNIT, cota_norma, gamma


def log(t0, *x):
    print(*x, f'[{time.time() - t0:.0f}s]', flush=True)


def prod_cert_desloc(A, p, X, bloco=2048):
    """fl((A + p) X) = fl(fl(A X) + fl(p X)) sem copiar A, e cota de ||fl - (A + p) X||_F:
    elemento a elemento, gamma_n |A||X| + 3u |p||X| + u |P| (|A||X| por blocos de linhas de A)."""
    P = A @ X
    P += p*X
    aX = np.abs(X)
    s2 = 0.0
    for i0 in range(0, A.shape[0], bloco):
        s2 += float(np.linalg.norm(np.abs(A[i0:i0 + bloco]) @ aX))**2
    del aX
    E = gamma(A.shape[1])*math.sqrt(s2)*(1 + 1e-6) + 3*UNIT*abs(p)*float(np.linalg.norm(X)) + UNIT*float(np.linalg.norm(P))
    return P, E


def norma_blocos_cert(blocos):
    return max(cota_norma(M) for M in blocos)


def sha256_arquivo(caminho):
    import hashlib
    h = hashlib.sha256()
    with open(caminho, 'rb') as fh:
        for bloco in iter(lambda: fh.read(1 << 24), b''):
            h.update(bloco)
    return h.hexdigest()


def blocos_abs_fro(blocos, offs, X):
    """|| |M| |X| ||_F com M bloco-diagonal (blocos) e X empilhado por blocos de linhas (offs): a cota
    elemento a elemento do arredondamento de fl(M X) e gamma_dim || |M| |X| ||_F."""
    s2 = 0.0
    for j, M in enumerate(blocos):
        s2 += float(np.linalg.norm(np.abs(M) @ np.abs(X[offs[j]:offs[j+1]])))**2
    return math.sqrt(s2)*(1 + 1e-6)


def tp_precond(A, p, rs, offs, NZ, nBhat, fJ0, est, log, t0):
    """||(A + p)^-1|| por Loewner PRE-CONDICIONADA (faixa B: resolvente ate ~200, onde a forma direta nao fecha:
    a margem de Rump ~ (2n) u tr(t^2 A A^*) cresce com ||A||_F^2 ~ 50^2 por linha).
    Com Q = Q_ZZ(p) = J_Z(p)^-1 EXATA e B^ = A - J_Z(0): Q (A + p) = I + Q B^ =: K, logo (A + p)^-1 = K^-1 Q e
    ||K^-1 Q|| <= t <=> t^2 K K^* - Q Q^* >= 0 (Loewner, como cota_inv_X). Calculados: K_c = I + Q_c B^_c (por blocos
    de linhas, um modo por vez) com ||K_e - K_c|| <= dK (Q_cert, J_Z(0) calculada, arredondamento elemento a
    elemento) e Q_c com ||Q_e - Q_c|| <= dq. Se ||K_c^-1 Q_c|| <= t:
      ||K_c^-1|| <= t ||Q_c^-1|| <= t nJ/(1 - nJ dq) =: kc   (Q_c^-1 = (I + J_Z(p)(Q_c - Q_e))^-1 J_Z(p)),
      ||(A + p)^-1|| = ||K_e^-1 Q_e|| <= t/(1 - kc dK) + dq kc/(1 - kc dK).
    A e alterada (vira B^) e liberada pelo chamador."""
    from verified_inverse import pd_hermitiana
    n = A.shape[0]
    J0 = JZ_blocos(rs, 0.0)
    Jp = JZ_blocos(rs, p)
    nJ = norma_blocos_cert(Jp) + 4*UNIT*math.sqrt(sum(float(np.linalg.norm(M))**2 for M in Jp))*1.0001
    Qs = []; dq = 0.0
    for r in rs:
        Qm, dqm, _ = Q_cert(r.m, NZ, p); Qs.append(Qm); dq = max(dq, dqm)
    dimz = max(int(offs[j+1] - offs[j]) for j in range(len(rs)))
    K = np.empty_like(A)
    e2 = 0.0
    for j in range(len(rs)):
        i0, i1 = int(offs[j]), int(offs[j+1])
        X = A[i0:i1].copy()
        X[:, i0:i1] -= J0[j]                       # B^ = A - J_Z(0): |fl - | <= u |B^_c| + 4u |J0| (J0 calculada)
        eX = UNIT*float(np.linalg.norm(X)) + 4*UNIT*float(np.linalg.norm(J0[j]))
        Y = Qs[j] @ X
        eY = gamma(dimz)*float(np.linalg.norm(np.abs(Qs[j]) @ np.abs(X)))*(1 + 1e-6)
        Y[:, i0:i1] += np.eye(i1 - i0)
        eY += UNIT*float(np.linalg.norm(Y))
        K[i0:i1] = Y
        e2 += (eY + float(np.linalg.norm(Qs[j], 2))*1.0001*eX)**2
        del X, Y
    dK = math.sqrt(e2)*(1 + 1e-6) + dq*(nBhat + fJ0*4*UNIT)    # + (Q_e - Q_c) B^_e
    del A                                            # memoria: o chamador nao guarda outra referencia
    aK = np.abs(K)
    eP = gamma(n)*float(aK.sum(axis=0).max())*float(aK.sum(axis=1).max())*1.01
    del aK
    P = K @ K.conj().T
    fP = float(np.linalg.norm(P))
    del K
    Sb = []; eS = 0.0; fS2 = 0.0
    for Qm in Qs:
        Sb.append(Qm @ Qm.conj().T)
        eS = max(eS, gamma(dimz)*float(np.linalg.norm(np.abs(Qm) @ np.abs(Qm).T))*(1 + 1e-6))
        fS2 += float(np.linalg.norm(Sb[-1]))**2
    fS = math.sqrt(fS2)
    log(t0, f'  p = {p}: K = I + Q B^ montada (dK {dK:.1e}, dq {dq:.1e}, ||J_Z(p)|| {nJ:.2f}); erro de P {eP:.1e}')
    t = est*1.003
    for _ in range(60):
        G = (t*t)*P
        for j, Sj in enumerate(Sb):
            i0, i1 = int(offs[j]), int(offs[j+1])
            G[i0:i1, i0:i1] -= Sj
        extra = t*t*eP + eS + 2*UNIT*(t*t*fP + fS)
        if pd_hermitiana([G], extra):
            break
        t *= 1.01
    else:
        raise RuntimeError('Loewner pre-condicionada nao verificou')
    del P
    kc = t*nJ/(1 - nJ*dq)
    assert kc*dK < 1e-3
    return t/(1 - kc*dK) + dq*kc/(1 - kc*dK), t


def cmd_globais(a):
    """Uma vez: ||B^|| = ||A - J_Z(0)|| (certificada) e ||A||_F."""
    t0 = time.time()
    info = json.load(open(a.dir/'schur.json')); MZ, NZ = info['Z']; Mc, Nc = info['F']
    _, ms, rs, offs = layout(MZ, NZ, Nc)
    A = np.load(a.dir/'A.npy')
    J0 = JZ_blocos(rs, 0.0)
    for j, M in enumerate(J0):
        A[offs[j]:offs[j+1], offs[j]:offs[j+1]] -= M
    nB = cota_norma(A, UNIT*float(np.linalg.norm(A))*4)
    log(t0, f'||B^|| = ||A - J_Z(0)|| <= {nB:.4f}')
    json.dump(dict(nBhat=nB), open(a.dir/'globais.json', 'w'), indent=1)


def cmd_ponto(a):
    t0 = time.time()
    info = json.load(open(a.dir/'schur.json')); MZ, NZ = info['Z']; Mc, Nc = info['F']
    g = json.load(open(a.dir/'globais.json'))
    c = a.c
    _, ms, rs, offs = layout(MZ, NZ, Nc); n = int(offs[-1])
    fat = json.load(open(a.dir/f'F_{rotulo(c)}.json'))
    shaA = sha256_arquivo(a.dir/'A.npy'); shaF = sha256_arquivo(a.dir/f'F_{rotulo(c)}.npy')
    log(t0, f'sha256 A {shaA[:16]}..., F {shaF[:16]}...')
    F = np.load(a.dir/f'F_{rotulo(c)}.npy'); Ft, Fb = F[:n], F[n:]
    A = np.load(a.dir/'A.npy')
    # A + p e formada em ponto flutuante na fase 1: fl(A_ii + p) = A_ii + p + e_i, |e_i| <= u |A_ii + p|, entao a
    # Loewner cota ||fl(A + p)^-1|| e ||(A + p)^-1|| <= t/(1 - t delta), delta = u max_i |A_ii + p|
    dA = np.abs(A.diagonal())
    # J_Z(0) e calculada em ponto flutuante (pesos): |J_calc - J| <= 4u |J_calc| elemento a elemento, entao a
    # B^ verdadeira (A - J_Z(0) exata) difere da de globais em <= 4u ||J_Z(0)||_F
    eJ0 = 4*UNIT*math.sqrt(sum(float(np.linalg.norm(M))**2 for M in JZ_blocos(rs, 0.0)))*1.0001
    nBhat = g['nBhat'] + eJ0
    pts = [complex(x) for x in a.pontos.split(',')]
    # t_p ja certificados (de um log anterior, arredondados PARA CIMA): '' = calcular
    dados = (a.tps.split(',') if a.tps else []) + [''] * len(pts)
    conhecidos = [float(x) if x else None for x in dados[:len(pts)]]
    # fase 1 (memoria: so A): t_p por Loewner, chute por LU + potencia
    from scipy.linalg import lu_factor, lu_solve
    from verified_inverse import gram_linhas, cota_loewner
    tps = []; precs = []
    for p, tconh in zip(pts, conhecidos):
        if tconh is not None:
            tps.append(tconh); precs.append(False)
            log(t0, f'p = {p}: t_p = ||(A + p)^-1|| <= {tconh:.6f}  (do log anterior)')
            continue
        A[np.diag_indices(n)] += p
        lu = lu_factor(A, check_finite=False)
        # chute preciso: Lanczos (ARPACK) no maior autovalor de (A + p)^-* (A + p)^-1, pela LU
        from scipy.sparse.linalg import LinearOperator, eigsh
        op = LinearOperator((n, n), matvec=lambda v: lu_solve(lu, lu_solve(lu, v), trans=2), dtype=complex)
        lam_ = eigsh(op, k=1, which='LM', tol=1e-8, return_eigenvectors=False)
        est = float(np.sqrt(np.max(np.abs(lam_))))
        del lu, op
        if a.precond or est > a.limiar_precond:
            del A                                    # memoria: tp_precond recebe a unica referencia e a libera
            tp, tl = tp_precond(np.load(a.dir/'A.npy'), p, rs, offs, NZ, nBhat,
                                math.sqrt(sum(float(np.linalg.norm(M))**2 for M in JZ_blocos(rs, 0.0))), est, log, t0)
            A = np.load(a.dir/'A.npy')
            tps.append(tp)
            log(t0, f'p = {p}: t_p = ||(A + p)^-1|| <= {tp:.6f}  (estimativa {est:.3f}; Loewner pre-condicionada {tl:.6f})')
            precs.append(True)
            continue
        # P = A A^* com cota ELEMENTO A ELEMENTO: |fl(A A^*) - A A^*| <= gamma_n |A| |A|^*, e
        # || |A| |A|^* ||_2 <= || |A| ||_1 || |A| ||_inf (a de Frobenius, gamma_n ||A||_F^2, e grande demais
        # para A sem pre-condicionamento: J tem entradas ate ~50 e n = 14 mil)
        P = A @ A.conj().T
        aA = np.abs(A)
        eP = gamma(n)*float(aA.sum(axis=0).max())*float(aA.sum(axis=1).max())*1.01
        del aA
        A[np.diag_indices(n)] -= p                   # A de volta (a soma e desfeita exatamente? nao: recarrega)
        # A sem pre-condicionamento: a margem de Rump (~ n u ||t^2 A A^*||) cresce com t^2 e obriga t acima do
        # valor verdadeiro (medido: 1,0025 em t ~ 5,5; 1,033 em t ~ 23; 1,12 em t ~ 44); comecar perto do
        # necessario poupa dezenas de Cholesky de n
        f0 = 1.005 + 0.13*(est/45.0)**2
        tp = cota_loewner(P, eP, None, 0.0, est*f0, passo=1.02, tentativas=40)
        del P
        A = np.load(a.dir/'A.npy')                   # sem erro de arredondamento acumulado
        tps.append(tp); precs.append(False)
        log(t0, f'p = {p}: t_p = ||(A + p)^-1|| <= {tp:.6f}  (estimativa {est:.3f})')
    # so grava se algum t_p foi calculado aqui: a fase 2 (todos de --tps) nao sobrescreve o tp da fase 1
    # (revisao da faixa B, L2: a origem de cada t_p, direta ou pre-condicionada, fica registrada)
    if any(x is None for x in conhecidos):
        json.dump(dict(c=[c.real, c.imag], pontos=[[p.real, p.imag] for p in pts], t_p=tps, precond=precs),
                  open(a.dir/f'tp_{rotulo(c)}_{a.rotulo}.json', 'w'), indent=1)
    if a.so_fase1:                                   # a fase 2 roda noutro processo (memoria: sem swap)
        return
    # correcao do deslocamento arredondado (vale tambem para os t_p reaproveitados: vieram da mesma fase 1)
    tps_ef = []
    for p, tp, pc in zip(pts, tps, precs):
        if pc:                                       # pre-condicionada: A + p nao e formada
            tps_ef.append(tp); continue
        delta = UNIT*float(np.max(np.abs(dA + p)))*1.0001
        assert tp*delta < 1e-6
        tps_ef.append(tp/(1 - tp*delta)*(1 + 1e-15))
    del dA
    # fase 2: Schur, SEM copias n x n (A, T e U ja ocupam ~9,4 GB em n = 14016): a diagonal de T e trocada
    # no lugar e restaurada exatamente; (A + p) X = A X + p X
    T = np.load(a.dir/'T.npy'); Us = np.load(a.dir/'U.npy')
    dT = T.diagonal().copy()
    UFt = Us.conj().T @ Ft
    Jc = JZ_blocos(rs, c)
    fJc = math.sqrt(sum(float(np.linalg.norm(M))**2 for M in Jc))
    nJc = norma_blocos_cert(Jc) + 4*UNIT*fJc*1.0001        # >= ||J_Z(c)|| exata
    saida = []
    for p, tp in zip(pts, tps_ef):
        # (b) Schur: X_k e resíduos certificados
        T[np.diag_indices(n)] = dT + p
        Xs = []; Z = UFt; erros = []; eprev = 0.0; Xprev = Ft
        for k in range(3):
            Z = solve_triangular(T, Z, check_finite=False)
            Xk = Us @ Z
            P, eP = prod_cert_desloc(A, p, Xk)
            Res = P - Xprev
            nRes = float(np.linalg.norm(Res))*(1 + 1e-6) + eP + UNIT*(float(np.linalg.norm(P)) + float(np.linalg.norm(Xprev)))
            ek = tp*nRes + tp*eprev                  # ||Y_k - X_k||_F
            Xs.append(Xk); erros.append(ek); eprev = ek; Xprev = Xk
            del P, Res
        T[np.diag_indices(n)] = dT
        del Z
        # (c) Q_ZZ(p) em blocos, certificado (Q_cert: erro dq por modo)
        Qp = []; dQ = 0.0; nQ = 0.0
        for r in rs:
            Q, dq, nq = Q_cert(r.m, NZ, p); Qp.append(Q); dQ = max(dQ, dq); nQ = max(nQ, nq + dq)
        # ||Q_ZZ^k F_b calc - exato|| (Frobenius) por recorrencia: Q exata com ||Q|| <= nQ, ||Q - Q_calc|| <= dQ por
        # bloco, e o arredondamento de cada produto elemento a elemento, gamma_dim || |Q| |X| ||_F
        dimz = max(int(offs[j+1] - offs[j]) for j in range(len(rs)))
        QF = [Fb]; dQF = [0.0]
        for k in range(3):
            X = QF[-1]
            eR = gamma(dimz)*blocos_abs_fro(Qp, offs, X)
            QF.append(aplica_blocos(Qp, offs, X))
            dQF.append(nQ*dQF[-1] + dQ*float(np.linalg.norm(X))*(1 + 1e-6) + eR)
        cotas = []
        for j in range(3):
            D = Xs[j] - QF[j + 1]                    # |fl(D) - D| <= u |D|
            fD = float(np.linalg.norm(D))*(1 + 1e-6)
            M = aplica_blocos(Jc, offs, D)
            # J_Z(c) exata contra a calculada (4u |Jc|) e o produto (gamma_dim |Jc| |D|)
            eM = nJc*(erros[j] + dQF[j + 1] + UNIT*fD) + (gamma(dimz) + 4*UNIT)*blocos_abs_fro(Jc, offs, D)
            del D
            cotas.append(cota_norma(M, eM))
        a0, b1, b2 = cotas
        nY3 = cota_norma(Xs[2], erros[2])
        nQ3F = cota_norma(QF[3], dQF[3])
        resto = fat['resto_F']
        r_ = dict(p=[p.real, p.imag], t_p=tp, a0=a0, b1=b1, b2=b2, nY3=nY3, nQ3F=nQ3F, nQZ=nQ, nJc=nJc,
                  resto_F=resto, qTT=fat['qTT'], nBhat=nBhat, erros=erros, dQF=dQF)
        saida.append(r_)
        log(t0, f'  a0 {a0:.4f}  b1 {b1:.3f}  b2 {b2:.2f}  ||Y3|| {nY3:.2f}  ||Q^3 F_b|| {nQ3F:.2f}  (erros Schur {erros[0]:.1e}, {erros[2]:.1e})')
    json.dump(dict(c=[c.real, c.imag], sha256_A=shaA, sha256_F=shaF, versao=2, pontos=saida),
              open(a.dir/f'rig_{rotulo(c)}_{a.rotulo}.json', 'w'), indent=1)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('cmd', choices=['globais', 'ponto'])
    ap.add_argument('--dir', type=Path, default=Path('build/rouche_L'))
    ap.add_argument('--c', type=complex, default=0j)
    ap.add_argument('--pontos', default='0')
    ap.add_argument('--chute', type=float, default=None, help='chute para ||(A + p)^-1|| (Loewner)')
    ap.add_argument('--rotulo', default='a')
    ap.add_argument('--tps', default='', help="t_p ja certificados por ponto (cotas superiores), '' = calcular")
    ap.add_argument('--precond', action='store_true', help='fase 1 sempre pela Loewner pre-condicionada')
    ap.add_argument('--so-fase1', dest='so_fase1', action='store_true', help='so a fase 1 (grava tp_*.json e sai)')
    ap.add_argument('--limiar-precond', dest='limiar_precond', type=float, default=40.0,
                    help='estimativa de ||(A + p)^-1|| acima da qual a fase 1 usa a Loewner pre-condicionada')
    a = ap.parse_args()
    dict(globais=cmd_globais, ponto=cmd_ponto)[a.cmd](a)


if __name__ == '__main__':
    main()
