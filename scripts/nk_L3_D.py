#!/usr/bin/env python3
"""S3b, L, modo rigoroso, ETAPA D: o TESTEMUNHO de constraints no autovetor enclausurado.

Com o certificado NK (C3.json): h* = Q0 y*, ||x* - x0|| <= r na norma do certificado, h0 = h0w
(Z.json), |lambda* - lambda0| <= v_0 r. Mostra Pi_F C(lambda*) h* != 0, com Pi_F a projecao
numa caixa F de saida (espaco de K: m par, |m| < MF, n < NF; pesos (129/128, 9/8)); entrada no
espaco de L (pesos (65/64, 5/4)). C(s) = C(0) + s C1, C1 h = (0, -xi h2).

  ||Pi_F C(lambda*) h*|| >= c_ref - v_0 r c_Z - v_T r sqrt(sum_J (c_J/w_J)^2) - v_f r c_far
                           - v_0 r ||Pi_F C1|| (1 + q r) - eps_C,
  c_ref = ||Pi_F C(lambda0) h0|| (vetor finito; cota inferior), c_Z = ||Pi_F C Q0 P_Z||,
  c_J = ||Pi_F C Q0 P_J|| (janelas de nk_L3), c_far: Pi_F C so ve entradas de n >= 250 com
  pesos (9/8)^NF/(5/4)^n, e Q0 P_far so desce ate la por (5/4)^-150; eps_C: dados de RT e mu.
Se o lado direito e > 0, o autovetor de L em lambda* VIOLA as constraints: lambda* e raiz de K
(K C = M L), logo lambda* = s_K (S3B_ROTA.md, 1 e 8)."""
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
import fourier_tail_bound as ft
import nk_L as nk
import nk_L3
import signed_constraints as sc
import signed_operator_L as sl
from nk_L2_Z import Q_cert, err_mm
from verified_norms import cota_norma, gamma, U

K1K, K2K = ft.K1, ft.K2          # saida: espaco de K


def pesos_K(so: sc.Saida):
    n0, n1 = so.comps[0][1], so.comps[1][1]
    return K1K**so.m*np.sqrt(np.concatenate([ft.omega(n0), ft.omega(n1)]))


def bloco_Cw(g, mo, NF, mi, Nin, s):
    so = sc.Saida(mo, NF); ri = sl.Radial(mi, Nin)
    X = sc.bloco_C(g, so, ri, s)
    return (pesos_K(so)[:, None]*X)/nk.pesos(ri)[None, :]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--dir', type=Path, default=Path('build/s3b/rig'))
    ap.add_argument('--F', default='12x48', help='caixa de saida Pi_F (espaco de K)')
    ap.add_argument('--W', type=int, default=16)
    a = ap.parse_args()
    z = json.loads((a.dir/'Z.json').read_text()); c3 = json.loads((a.dir/'C3.json').read_text())
    t3 = json.loads((a.dir/'T3.json').read_text())
    MF, NF = map(int, a.F.split('x'))
    MZ, NZ = z['Z']; Mc, Nc = z['F']; lam0 = z['lam0']
    lam0 = complex(*lam0) if isinstance(lam0, list) else lam0      # faixa B: lambda0 complexo
    g = sl.campos()
    h0w = np.load(a.dir/'h0w.npy')
    msZ = sl.modos(MZ); rs = [sl.Radial(m, NZ) for m in msZ]
    offs = np.cumsum([0] + [r.dim for r in rs])
    saidas = [m for m in sl.modos(MF) if m % 2 == 0]
    BM = sl.BAND_M
    # c_ref: Pi_F C(lambda0) h0
    partes = []; e2 = 0.0
    for mo in saidas:
        acc = 0
        for j, mi in enumerate(msZ):
            if abs(mo - mi) <= BM:
                Cb = bloco_Cw(g, mo, NF, mi, NZ, lam0); h = h0w[offs[j]:offs[j+1]]
                acc = acc + Cb @ h; e2 += err_mm(Cb, h[:, None])**2
        partes.append(acc)
    c = np.concatenate(partes)
    # dados de RT e mu em C (nota das normas do testemunho, em l1: (81/8) eps_omega + 45,19 eps_mu);
    # aqui, conservador: 100 eps_omega + 100 eps_mu, vezes ||h0|| = 1, e o fator de norma l2 <= l1
    from tile_perron import DELTA_MU
    # mu: C e afim em mu; ||Pi_F (C(mu=1) - C(mu=0)) h0|| calculado, vezes |mu_true - mu_RefA|
    cmu = []
    for mo in saidas:
        acc = 0
        for j, mi in enumerate(msZ):
            if abs(mo - mi) <= BM:
                so = sc.Saida(mo, NF); ri = sl.Radial(mi, NZ)
                X = sc.bloco_C(g, so, ri, lam0, mu=1.0) - sc.bloco_C(g, so, ri, lam0, mu=0.0)
                acc = acc + ((pesos_K(so)[:, None]*X)/nk.pesos(ri)[None, :]) @ h0w[offs[j]:offs[j+1]]
        cmu.append(acc)
    epsC = 100*(2.0**-25 + 2.0**-277) + DELTA_MU*float(np.linalg.norm(np.concatenate(cmu)))*1.01
    c_ref = float(np.linalg.norm(c))*(1 - 1e-6) - math.sqrt(e2) - U*float(np.linalg.norm(c)) - epsC
    print(f'c_ref = ||Pi_F C(lambda0) h0|| >= {c_ref:.5f}   (F = {MF}x{NF}, |h0| = {np.linalg.norm(h0w):.6f})')
    # c_Z = ||Pi_F C Q0 P_Z||
    def bloco_linhas(cols_por_modo, Nin, Qs):
        """[Pi_F C Q0] nas colunas dadas: blocos (saidas) x (modos de entrada), com erro."""
        mi_l = [m for m in cols_por_modo if len(cols_por_modo[m])]
        oc = np.cumsum([0] + [len(cols_por_modo[m]) for m in mi_l])
        rows = []; e2 = 0.0; dQ = 0.0
        for mo in saidas:
            linha = np.zeros((sc.Saida(mo, NF).dim, int(oc[-1])), complex)
            for q, mi in enumerate(mi_l):
                if abs(mo - mi) <= BM:
                    Cb = bloco_Cw(g, mo, NF, mi, Nin, lam0)
                    Q, dq, _ = Qs(mi)
                    R = Q[:, cols_por_modo[mi]]
                    linha[:, oc[q]:oc[q+1]] = Cb @ R
                    e2 += err_mm(Cb, R)**2 + (float(np.linalg.norm(Cb))*dq)**2
            rows.append(linha)
        return np.concatenate(rows), math.sqrt(e2)
    QZc = {}
    def QsZ(m):
        if m not in QZc:
            QZc[m] = Q_cert(m, NZ, lam0)
        return QZc[m]
    XZ, eZ = bloco_linhas({m: np.arange(sl.Radial(m, NZ).dim) for m in msZ}, NZ, QsZ)
    cZ = cota_norma(XZ, eZ); del XZ
    print(f'c_Z = ||Pi_F C Q0 P_Z|| <= {cZ:.4f}')
    # c_J nas janelas de nk_L3 alcancadas (|m| < MF + BAND_M)
    import nk_L as nkL
    ctx = nkL.Contexto(lam0, Nc)
    js = nk_L3.janelas_L3(MZ, Mc, a.W)
    assert [f'{w[0]}..{w[-1]}' for w in js] == t3['nome']
    QTc = {}
    def QsT(m):
        if m not in QTc:
            if len(QTc) > 40:
                QTc.clear()
            QTc[m] = Q_cert(m, Nc, lam0)
        return QTc[m]
    cJ = []
    for w in js:
        if min(abs(m) for m in w) >= MF + BM:
            cJ.append(0.0); continue
        cols = {m: nkL.colsT(ctx, m, MZ, NZ, Nc) for m in w}
        X, e = bloco_linhas(cols, Nc, QsT)
        cJ.append(cota_norma(X, e)); del X
        print(f'  c_J[{w[0]}..{w[-1]}] <= {cJ[-1]:.4f}', flush=True)
    # far: Pi_F C P_(n >= 250) e P_(n < 250) Q0 P_far (desprezaveis; cotas explicitas)
    #  |entrada pesada de Pi_F C na coluna n| <= (sup de C em Chebyshev ~ 4 n + 40) (9/8)^NF/(5/4)^n,
    #  soma de Schur em n >= 250 e linhas < 2 NF: majorante folgado
    #  majorante generoso de |entrada pesada| na coluna n: (1e3 n + 1e5) (9/8)^NF (65/64)^52 / (5/4)^n;
    #  Frobenius sobre <= 2 NF MF linhas e 4 componentes x 2 (MF + BM) modos de entrada
    assert Nc >= 350, 'a descida P_(n<250) Q0 P_(n>=Nc) precisa de Nc >= 350'
    col = lambda n: (1e3*n + 1e5)*(9/8)**NF*(65/64)**52/(5/4)**n
    # termos decrescem com razao <= ((n + 101)/n)^2 (4/5)^2 < 0,7 para n >= 250: cauda geometrica
    fro2 = (sum(col(n)**2 for n in range(250, 3000)) + col(3000)**2/(1 - 0.7))*2*NF*MF*4*2*(MF + BM)
    c_far = math.sqrt(fro2) + 1e3*cf.descida(250, Nc, sl.MU, sl.K2)
    print(f'c_far <= {c_far:.1e}')
    # C1: (0, -xi h2); ||xi||_{L -> K} <= (9/8 + 1)/2 folgado
    nC1 = 1.5                      # ||xi||_{L->K} <= (5/4) max razao de pesos em m (<= 1,1), folgado
    w = [float(Fr(x)) for x in c3['w']]; v = [float(Fr(x)) for x in c3['v']]; r = float(Fr(c3['r']))
    SJ = math.sqrt(sum((cJ[i]/w[i])**2 for i in range(len(cJ))))
    # perturbacao global de C (dados de RT e mu) sobre h* - h0 = Q0 (y* - y0): nao respeita a banda
    wmin = min(w); nQmax = max(t3['nQ'] + [z['nQZ'], z['Qfar']])
    epsCg = 100*(2.0**-25 + 2.0**-277) + DELTA_MU*1e3
    perda = v[0]*r*cZ + v[1]*r*SJ + v[2]*r*c_far + v[0]*r*nC1*(1 + 10*r) \
            + epsCg*nQmax*(v[0] + v[1]/wmin + v[2])*r
    print(f'perdas: Z {v[0]*r*cZ:.2e}, cauda {v[1]*r*SJ:.2e}, far {v[2]*r*c_far:.1e}, lambda {v[0]*r*nC1:.2e};'
          f' total {perda:.3e}')
    folga = c_ref - perda
    print(f'||Pi_F C(lambda*) h*|| >= {c_ref:.5f} - {perda:.3e} = {folga:.5f}  ->  '
          + ('TESTEMUNHO: o autovetor de L em lambda* VIOLA as constraints' if folga > 0 else 'NAO certificado'))
    (a.dir/'D.json').write_text(json.dumps(dict(F=[MF, NF], c_ref=c_ref, cZ=cZ, cJ=cJ, c_far=c_far, nC1=nC1,
                                                perda=perda, folga=folga, certificado=folga > 0), indent=1) + '\n')


if __name__ == '__main__':
    main()
