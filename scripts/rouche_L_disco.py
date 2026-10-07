#!/usr/bin/env python3
"""Montagem RIGOROSA do Rouché de L num pedaço do contorno: para cada ponto p (nível Z de rouche_L_rig)
e o centro de cauda c (T3.json de nk_L3_T --rouche), o maior raio r tal que o Perron de 3 níveis
(Z', T, far) fica <= alvo para TODO s do contorno com |s - p| <= r (C1_REAVALIACAO.md §7).

Norma em s (escala s-dependente no nível Z, D_Z = J_Z(c) Q_ZZ(s); pesos w nas janelas, os de nk_L3_C):
  Z'<-T  [E ||Phi|| + (tau + JQ) rho  (E = a0 + r b1 + r^2 b2 + r^3 (tau ||Y3|| + JQ ||Q^3 F_b||); (||Phi||, rho) de certF)
          + tau ((dB + dEZ) nQx + emu_s)] / (1 - d q_TT) / min_SZ w,
         tau = 1 + (d + ||B^||) t_p/(1 - r t_p) >= ||J_Z(c) R(s)||,  JQ = 1 + d nQZ_s >= ||J_Z(c) Q_ZZ(s)||;
  Z'<-Z' tau ((dB + dEZ) ||Q_ZZ(c)|| + emu_c);  Z'<-far  tau (||B|| descidas + (dB + dEZ) Qf_s + emu_s + ||E|| q_Zf);
  T<-Z'  TZ_J(c) exato na escala (+ variação de mu e linhas de T da deflação verdadeira);
  T<-T   N(c) com sensibilidade a s (d = |p - c| + r), pela identidade do resolvente e pelas identidades
         triangulares Q0(s) P_J = P_Z Q_ZJ(s) + P_J Q_JJ(s), Q0(c) P_Z = P_Z Q_ZZ(c),
         (I + (s - c) Q_JJ(c))^-1 = I - (s - c) Q_JJ(s), fJ = 1/(1 - d nQ_JJ) >= 1 + d ||Q_JJ(s)||:
           V_J P_J B (Q0(s) - Q0(c)) P_J = -(s - c) [V_J P_J B P_Z Q_ZZ(c) Q_ZJ(s)
                                  + (V_J P_J B Q0(c)^2 P_J - V_J P_J B P_Z Q_ZZ(c) Q_ZJ(c)) (I - (s - c) Q_JJ(s))];
         janela J central:     rho_J^s = d [TZ_J q_ZJ(s) + (K_J + TZ_J q_ZJ(c)) fJ];
         janela J não central: rho_J^s = d K_J fJ (sem linhas de Z);
         fora da diagonal:     N_ij(s) <= N_ij(c) + d [TZ_i q_Zj(s) + N_ij(c) nQ_jj fJ_j] (q_Zj = 0 se j não central);
         q_ZJ(s) <= qz_J(p)/(1 - r nQZ(p)) (1 + d nQ_JJ fJ), qz_J(p) = ||Q_ZZ(p) J_ZJ Q_JJ(c)||, q_ZJ(c) = ||Q_ZJ(c)||;
  far    beta_s = (||B|| + rt) Qf_s + emu_s, Qf_s = sup no disco de ||Q0(s) P_far|| (todas as linhas; formas fechadas
         com majorantes no disco, Qfar_disco); Nf(s) = Nf(c)/(1 - d nQ_JJ) (P_far B P_Z = 0 na banda: o termo de Z
         da inversa em blocos nao aparece; o resto da banda esta na perturbacao global);
  perturbação global eg_s = rt + emu_s em toda entrada, como em nk_L3_C (L1, L2 da revisão de S3b incluídas), com
         ||Q0(s) P_J|| <= (||Q0(c) P_J|| + d ||Q_ZZ(s)|| ||Q_ZJ(c)||) fJ e fator >= 1 (eps_mu ja inclui Q0(s));
  deflação verdadeira fora de Z (revisão de 01/10, L1): dET (nQx nas colunas de T, Qf_s nas do far, ||Q_ZZ(c)|| nas
         de Z') nas linhas de T (max w nV) e do far.
Só pontos do CONTORNO (Re s >= 0): as formas fechadas de Q0 valem lá."""
from __future__ import annotations

import argparse
from fractions import Fraction as Fr
import json
import math
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import closed_form as cf
import fundo_L
import nk_L as nk
import signed_operator_L as sl
from nk_L2_Z import Q_cert
from nk_L3_C import Q, raiz_sup, norma2_cota
from verified_norms import U as UNIT, cota_norma

MU, K2 = nk.MU, nk.K2


def dados_centrais(t3, c, MZ, NZ, Nc):
    """Por janela central (modos |m| < MZ): matrizes X_m = J_ZJ,m Q_JJ,m(c) por modo, jq = max ||X_m||,
    nQJJ = max ||Q_JJ,m(c)|| (certificados)."""
    ctx = nk.Contexto(c, Nc)
    out = {}
    for i, modos in enumerate(t3['modos']):
        if min(abs(m) for m in modos) >= MZ:
            continue
        Xs = {}; jq = 0.0; nq = 0.0; qzc = 0.0
        for m in modos:
            if abs(m) >= MZ:
                continue
            Qm, dq, nqm = Q_cert(m, Nc, c)
            r = ctx.rad(m, Nc); w = np.sqrt(np.concatenate([cf.omega(ns, K2) for _, ns in r.comps]))
            Jw = (w[:, None]*sl.J_livre(r, 0.0))/w[None, :]
            nn = np.concatenate([ns for _, ns in r.comps])
            iz = np.where(nn < NZ)[0]; iT = np.where(nn >= NZ)[0]
            QJJ = Qm[np.ix_(iT, iT)]
            qzc = max(qzc, cota_norma(Qm[np.ix_(iz, iT)], dq))
            X = Jw[np.ix_(iz, iT)] @ QJJ
            eX = cota_norma(Jw[np.ix_(iz, iT)])*dq + 4*UNIT*float(np.linalg.norm(X))*len(iT)
            Xs[m] = (X, eX)
            jq = max(jq, cota_norma(X, eX)); nq = max(nq, cota_norma(QJJ, dq))
        out[i] = dict(X=Xs, jq=jq, nQJJ=nq, qZc=qzc)
    return out


def qz_ponto(cent, p, NZ):
    """qz_J(p) = max_m ||Q_ZZ,m(p) X_m|| por janela central, certificado."""
    res = {}
    cache = {}
    for i, d in cent.items():
        v = 0.0
        for m, (X, eX) in d['X'].items():
            if m not in cache:
                Qz, dqz, nqz = Q_cert(m, NZ, p); cache[m] = (Qz, dqz, nqz)
            Qz, dqz, nqz = cache[m]
            Y = Qz @ X
            eY = (nqz + dqz)*eX + dqz*cota_norma(X) + 4*UNIT*float(np.linalg.norm(Y))*X.shape[0]
            v = max(v, cota_norma(Y, eY))
        res[i] = v
    return res


def Qfar_disco(p, r, Mc, Nc):
    """sup ||Q0(s) P_far|| (TODAS as linhas) sobre s do contorno com |s - p| <= r (Re s >= 0): formas fechadas
    com majorantes no disco (closed_form.cauda_Rinv_disco) nos modos |m| < Mc, e Rinv_uniforme com
    h = Mc/2 - |Im p| - r nos modos |m| >= Mc. Substitui Qf(c)/(1 - d Qf(c)), que so vale no bloco far<-far
    (revisao independente de 01/10, L2: as linhas de Z/T de Q0(s) P_far trazem (I + (s - c) Q_NN)^-1)."""
    q = 0.0
    for m in range(-(Mc - 1), Mc):
        for passo in ((2, 1) if m % 2 == 0 else (1,)):
            q = max(q, cf.cauda_Rinv_disco(p + 1j*m/2, r, Nc, passo, MU, K2, 1500))
    cw = math.sqrt(1 + K2**-4)
    h = Mc/2 - abs(p.imag) - r
    assert h > 0
    return max(q, (1 + 2*cw/(K2 - 1))/h*1.0001)


def theta_disco(p, r, zp, t3, c, cent, qz, g, exato=False, Qf_disco=None):
    """Perron majorante (entrada a entrada) em todos os s do contorno com |s - p| <= r. Devolve (theta, detalhes)."""
    Fl = (lambda x: Q(x)) if exato else float
    d = abs(p - c) + r
    tp, nQZp = zp['t_p'], zp['nQZ']
    if r*tp >= 0.9 or r*nQZp >= 0.9:
        return float('inf'), None
    if Qf_disco is None:                                  # so na bissecao (ponto flutuante, heuristica)
        assert not exato, 'o Perron rigoroso exige Qf_disco'
        if d*t3['Qfar'] >= 0.9:
            return float('inf'), None
        Qf_s = t3['Qfar']/(1 - d*t3['Qfar'])
    else:
        Qf_s = Qf_disco
    t_s = tp/(1 - r*tp)
    tau = 1 + (d + zp['nBhat'])*t_s
    nQZs = nQZp/(1 - r*nQZp)
    JQ = 1 + d*nQZs
    cF = g.get('certF')
    qTT = cF['qTT'] if cF else zp['qTT']
    if d*qTT >= 0.9:
        return float('inf'), None
    # mu: eps_mu_L crescente em |Re s| e |Im s|: pior canto do disco
    emu_s = fundo_L.eps_mu_L(complex(abs(p.real) + r, abs(p.imag) + r))
    emu_c = t3['eps_mu']
    rt, arred, fb_tot = g['rt'], g['arred40'], t3['dB_L']['total']
    dB = arred + rt
    dEZ = g['dEZ']; dET = g['dET']
    nQx = g['nQx_s']
    E = zp['a0'] + r*zp['b1'] + r*r*zp['b2'] + r**3*(tau*zp['nY3'] + JQ*zp['nQ3F'])
    # fator: ||X S|| <= ||X F|| ||Phi|| + ||X|| rho (certF, o melhor par); sem certF, a forma antiga com resto_F
    if cF:
        EF = min(E*pr['phi'] + (tau + JQ)*pr['rho'] for pr in cF['pares'])
    else:
        EF = E + (tau + JQ)*zp['resto_F']
    ZT0 = (EF + tau*((dB + dEZ)*nQx + emu_s))/(1 - d*qTT)
    ZZ = tau*((dB + dEZ)*t3['nQZ'] + emu_c)
    Zf = tau*(nk.NORMA_BL*(cf.descida(g['NZ'] + nk.BAND_N + 1, g['Nc'], MU, K2) + 6*cf.descida_divxi(g['NZ'], g['Nc'], K2))
              + (dB + dEZ)*Qf_s + emu_s + g['nE']*cf.descida(g['NZ'], g['Nc'], MU, K2))
    # cauda com sensibilidade
    N = np.array(t3['N']); TZ = np.array(t3['TZ']); nV = np.array(t3['nV']); nVhi = np.array(t3['nVhi'])
    rho = np.diag(N).copy(); K = np.array(t3['K']); nQm = np.array(t3['nQm']); Nf = np.array(t3['Nf'])   # rho na diagonal de N
    k = len(nV)
    nQJJ = np.array([cent[i]['nQJJ'] if i in cent else nQm[i] for i in range(k)])
    if np.any(d*nQJJ >= 0.9):
        return float('inf'), None
    fJ = 1/(1 - d*nQJJ)
    qZ = {i: qz[i]/(1 - r*nQZp)*(1 + d*nQJJ[i]*fJ[i]) for i in cent}
    # ||Q0(s) P_J|| <= (||Q0(c) P_J|| + d ||Q_ZZ(s)|| ||Q_ZJ(c)||) fJ (centrais; pela inversa em blocos de
    # I + (s - c) Q0(c)); nas nao centrais ||Q0(c) P_J|| fJ
    nQ_c = np.array(t3['nQ'])
    nQs = np.array([(nQ_c[i] + (d*nQZs*cent[i]['qZc'] if i in cent else 0.0))*fJ[i] for i in range(k)])
    Ns = np.zeros((k, k))
    dpert = max(0.0, emu_s - emu_c)
    # sensibilidade a s (d = |s - c|; formulas e identidades na docstring; so janelas centrais tem a parte Z):
    #   diagonal: rho_J^s <= d [TZ_J q_ZJ(s) + (K_J + TZ_J q_ZJ(c)) fJ]  (nao centrais: d K_J fJ)
    #   fora:     N_ij(s) <= N_ij(c) + d [TZ_i q_Zj(s) + N_ij(c) nQ_jj fJ_j]  (q_Zj = 0 nas nao centrais)
    for i in range(k):
        if i in cent:
            rs_ = d*(TZ[i]*qZ[i] + (K[i] + TZ[i]*cent[i]['qZc'])*fJ[i])
        else:
            rs_ = d*K[i]*fJ[i]
        Ns[i, i] = rho[i] + rs_ + nV[i]*dpert*max(1.0, nQs[i])   # eps_mu ja inclui Q0(s): fator >= 1
        for j in range(k):
            if j == i or N[i, j] == 0:
                continue
            Ns[i, j] = N[i, j] + d*((TZ[i]*qZ[j] if j in cent else 0.0) + N[i, j]*nQJJ[j]*fJ[j])
            Ns[i, j] += nV[i]*dpert*max(1.0, nQs[j])
    TZs = TZ + nV*(dpert*t3['nQZ'] + dET*t3['nQZ'])
    Nfs = Nf*fJ[None, :]
    beta_s = (nk.NORMA_BL + rt)*Qf_s + emu_s
    lim = g['Nc'] - nk.BAND_N - t3.get('folga_far', 60)
    nB = nk.NORMA_BL + fb_tot
    desc_s = nB*(cf.descida(lim + nk.BAND_N + 1, g['Nc'], MU, K2) + 6*cf.descida_divxi(lim, g['Nc'], K2)) + (fb_tot + emu_s)*Qf_s
    eg = rt + emu_s
    SZ = t3['SZ']
    Mc = g['Mc']                                        # L1, L2 da revisao de S3b (160, 176 com Mc = 200)
    borda = [max(abs(m) for m in t3['modos'][i]) >= Mc - nk.BAND_M for i in range(k)]
    nQb = max(nQs[i] for i in range(k) if max(abs(m) for m in t3['modos'][i]) >= Mc - 24)
    nQmax = max(1.0, float(nQs.max()))                  # eg = rt + emu_s: rt pede ||Q0(s) P_J||, emu ja inclui Q0
    # deflacao verdadeira fora de Z (revisao de 01/10, L1): delta_k (x_k - P_Z x_k^A) phi_k^* Q0(s), phi_k em Z,
    # ||(I - P_Z) x_k|| <= cauda_x + erro, soma dET; colunas: Z' (escala: Q_ZZ(c)), T (||P_Z Q0(s) P_T|| <= nQx),
    # far (Qf_s); linhas T (V_J, max w nV) e far
    dfT_T = dET*nQx; dfT_f = dET*Qf_s
    # pesos: 1 e sqrt(u/z) de N(c)
    lam, Vr = np.linalg.eig(N); zv = np.abs(Vr[:, np.argmax(np.abs(lam))])
    lam, Vl = np.linalg.eig(N.T); uv = np.abs(Vl[:, np.argmax(np.abs(lam))])
    melhor = (float('inf'), None)
    for nomew, w in (('1', np.ones(k)), ('sqrt(u/z)', np.sqrt(np.maximum(uv, 1e-9)/np.maximum(zv, 1e-9)))):
        w = w/w.max()
        if exato:
            wq = [Fr(float(x)) for x in w]
            nT = norma2_cota((w[:, None]*Ns)/w[None, :]); fT = norma2_cota(Nfs/w[None, :])
            TZw = raiz_sup(sum((wq[i]*Q(TZs[i]))**2 for i in range(k)))
            wmin = min(wq[j] for j in SZ); wmin0 = min(wq)
            maxwV0 = max(wq[i]*Q(nV[i]) for i in range(k))
            Tf = max(wq[i]*Q(nV[i] if borda[i] else nVhi[i]) for i in range(k))*Q(beta_s) \
                + raiz_sup(sum((wq[i]*Q(nV[i]))**2 for i in range(k)))*Q(desc_s)
            fT += Q(t3['dB_L']['resto'])*Q(nQb)/wmin0
            egq = Q(eg)
            nT += maxwV0*egq*Q(nQmax)/wmin0; TZw += maxwV0*egq*Q(t3['nQZ']); Tf += maxwV0*egq*Q(Qf_s)
            ZT = Q(ZT0)/wmin + Q(tau)*egq*Q(nQmax)/wmin0; fT += egq*Q(nQmax)/wmin0
            nT += maxwV0*Q(dfT_T)/wmin0; Tf += maxwV0*Q(dfT_f); fT += Q(dfT_T)/wmin0           # L1
            fZ = Q(fb_tot + emu_s)*Q(t3['nQZ']) + Q(dET)*Q(t3['nQZ']); ff = Q(beta_s) + Q(dfT_f)
            M = [[Q(ZZ), ZT, Q(Zf)], [TZw, nT, Tf], [fZ, fT, ff]]
            Mf = np.array([[float(x) for x in l] for l in M])
            lamM, vec = np.linalg.eig(Mf); v = np.abs(vec[:, np.argmax(np.abs(lamM))]); v = np.maximum(v/v.max(), 1e-9)
            vq = [Fr(float(x)) for x in v]
            th = max(sum(M[i][j]*vq[j] for j in range(3))/vq[i] for i in range(3))
        else:
            nT = float(np.linalg.norm((w[:, None]*Ns)/w[None, :], 2)); fT = float(np.linalg.norm(Nfs/w[None, :], 2))
            TZw = math.sqrt(float(np.sum((w*TZs)**2))); wmin = min(w[j] for j in SZ); wmin0 = float(w.min())
            maxwV0 = float(np.max(w*nV))
            Tf = max(float(w[i]*(nV[i] if borda[i] else nVhi[i])) for i in range(k))*beta_s + math.sqrt(float(np.sum((w*nV)**2)))*desc_s
            fT += t3['dB_L']['resto']*nQb/wmin0
            nT += maxwV0*eg*nQmax/wmin0; TZw += maxwV0*eg*t3['nQZ']; Tf += maxwV0*eg*Qf_s
            ZT = ZT0/wmin + tau*eg*nQmax/wmin0; fT += eg*nQmax/wmin0
            nT += maxwV0*dfT_T/wmin0; Tf += maxwV0*dfT_f; fT += dfT_T/wmin0                        # L1
            M = np.array([[ZZ, ZT, Zf], [TZw, nT, Tf], [(fb_tot + emu_s + dET)*t3['nQZ'], fT, beta_s + dfT_f]])
            th = float(max(abs(np.linalg.eigvals(M))))
        if float(th) < float(melhor[0]):
            melhor = (th, dict(pesos=nomew, M=[[float(x) for x in l] for l in M], ZT=float(ZT), nT=float(nT), TZw=float(TZw),
                                tau=tau, E=E, d=d, Qf_s=Qf_s, nQmax_s=nQmax))
    return melhor


def contexto(dir_, zrig, cauda, E, dtauB, certF, shaA):
    """Preparacao do disco (usada por main e por scripts/verifica_T2.py): le schur.json, o nivel Z, a cauda e o
    certF, confere os sha256 e monta as constantes globais g e os dados centrais da cauda."""
    dir_ = Path(dir_)
    info = json.load(open(dir_/'schur.json')); MZ, NZ = info['Z']; Mc, Nc = info['F']
    zr = json.load(open(zrig)); t3 = json.load(open(cauda))
    c = complex(*zr['c']); lc = t3['lam0']; assert complex(*lc) == c, 'centro do nivel Z e da cauda diferentes'
    dfr = info['defl_refA']
    eE = json.load(open(E))['e'] if Path(E).exists() else 2.51e-8
    fb40 = fundo_L.dB_L(sl.campos(), nk.BAND_M)
    errZ = [dtauB + 1e-60 + UNIT*dfr['norma_x'][0]*10, eE + UNIT*dfr['norma_x'][1]*10]
    g = dict(rt=fb40['rt'], arred40=fb40['arred'], NZ=NZ, Nc=Nc, Mc=Mc,
             dEZ=sum(dk*nphi*e for dk, nphi, e in zip(dfr['delta'], dfr['norma_phi'], errZ)),
             dET=sum(dk*nphi*(cx + e) for dk, nphi, cx, e in zip(dfr['delta'], dfr['norma_phi'], dfr['cauda_x'], errZ)),
             nE=sum(dk*nphi*nx for dk, nphi, nx in zip(dfr['delta'], dfr['norma_phi'], dfr['norma_x'])),
             nQx_s=(1 + 2*math.sqrt(1 + K2**-4)/(K2 - 1))/MU*1.0001)
    if certF:
        cF = json.load(open(certF)); assert complex(*cF['c']) == c, 'centro do certF diferente'
        g['certF'] = cF
        if 'sha256_F' in zr:                        # rig versao 2: o F do nivel Z e o do certF; A e a de referencia
            assert zr['sha256_F'] == cF['sha256_F'], 'F do rig diferente do F do certF'
            assert zr['sha256_A'] == shaA, 'A do rig diferente da A de referencia'
        else:
            print('AVISO: rig sem sha256 (versao 1)')
    else:
        print('AVISO: sem --certF, o fator usa resto_F e q_TT de ponto flutuante (nao rigoroso)')
    cent = dados_centrais(t3, c, MZ, NZ, Nc)
    return zr, t3, c, cent, g, MZ, NZ, Mc, Nc


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--dir', type=Path, default=Path('build/rouche_L'))
    ap.add_argument('--zrig', required=True, help='JSON de rouche_L_rig ponto (centro c)')
    ap.add_argument('--cauda', required=True, help='T3.json de nk_L3_T --rouche no centro c')
    ap.add_argument('--E', default='build/s4/rig/E.json', help='etapa E de S4: ||g* - g_A||')
    ap.add_argument('--dtauB', type=float, default=1.4563e-09, help='||d_tau RefB||_L (gauge_refB)')
    ap.add_argument('--certF', type=Path, default=None, help='certF_<c>.json de rouche_L_Z certF (obrigatorio para o rigoroso)')
    ap.add_argument('--alvo', type=float, default=0.99)
    ap.add_argument('--shaA', default='70791ce0c30f8065f93770c36ff89c1d8290978c5abf84ab0570313bc09fd002',
                    help='sha256 da A.npy de referencia (a mesma nas tres maquinas)')
    ap.add_argument('--saida', type=Path, default=None)
    a = ap.parse_args()
    zr, t3, c, cent, g, MZ, NZ, Mc, Nc = contexto(a.dir, a.zrig, a.cauda, a.E, a.dtauB, a.certF, a.shaA)
    out = []
    for zp in zr['pontos']:
        p = complex(*zp['p'])
        qz = qz_ponto(cent, p, NZ)
        lo, hi = 0.0, 0.5
        for _ in range(40):                         # maior r com theta <= alvo (ponto flutuante)
            mid = (lo + hi)/2
            th, _ = theta_disco(p, mid, zp, t3, c, cent, qz, g)
            lo, hi = (mid, hi) if th <= a.alvo else (lo, mid)
        r = lo*0.999
        for _ in range(20):                         # rigoroso com sup ||Q0(s) P_far|| no disco; recua se falhar
            qf = Qfar_disco(p, r, Mc, Nc)
            th, det = theta_disco(p, r, zp, t3, c, cent, qz, g, exato=True, Qf_disco=qf)
            if det is not None and float(th) < 1:
                break
            r *= 0.95
        ok = det is not None and float(th) < 1
        out.append(dict(p=[p.real, p.imag], r=r, t_p=zp['t_p'], theta=str(th), theta_f=float(th), ok=bool(ok), det=det))
        print(f'p = {p.real:.4f}{p.imag:+.4f}i: raio {r:.3e}, theta (racional) {float(th):.6f}  '
              f'[Z<-T {det["ZT"]:.4f}, T<-T {det["nT"]:.4f}, T<-Z {det["TZw"]:.4f}, tau {det["tau"]:.1f}, pesos {det["pesos"]}]'
              + ('' if ok else '  FALHA'), flush=True)
    if a.saida:
        a.saida.write_text(json.dumps(dict(c=[c.real, c.imag], certF=str(a.certF) if a.certF else None, versao=2,
                                           argumentos=dict(alvo=a.alvo, cauda=str(a.cauda), dtauB=a.dtauB, E=str(a.E),
                                                           dir=str(a.dir), shaA_ref=a.shaA),
                                           zrig=str(a.zrig), sha256_A=zr.get('sha256_A'), sha256_F=zr.get('sha256_F'),
                                           pontos=out), indent=1) + '\n')


if __name__ == '__main__':
    main()
