#!/usr/bin/env python3
"""S3b, L, modo rigoroso, ETAPA C da arquitetura nk_L3: montagem do certificado NK (racionais).

Le Z.json (etapa A) e T3.json (nk_L3_T.py). Norma: max(||x_Z'||/v_0, ||x_T||_W/v_T, ||x_f||/v_f),
||x_T||_W^2 = sum_J w_J^2 ||x_J||^2 (w escolhido em ponto flutuante; tudo o mais e cota). Entradas
de ||E|| por nivel (E = I - A DF(x0), A = diag(Mh^-1, V_J, I)):
  Z'<-Z' rhoZ;  Z'<-T a/min_{SZ} w;  Z'<-far epsZ;
  T<-Z' sqrt(sum (w_J TZ_J)^2);  T<-T ||W N W^-1||_2 (Rump na matriz pequena);
  T<-far max_J (w_J ||V_J P_alto||) beta + sqrt(sum (w_J ||V_J||)^2) desc;
  far<-Z' pert ||Q_ZZ||;  far<-T ||Nf W^-1||_2;  far<-far beta.
Perron 3x3 exato; 2a ordem: DF(x) - DF(x0) = [(lambda - lambda0) Q0 dy + dlambda Q0 (y - y0); 0],
Q0 diagonal em m e triangular superior em n (P_T Q0 P_Z = 0, P_far Q0 P_(Z,T) = 0):
  Z2 = 2 v_0 max( (tQ v_0 + t_V (qZT v_T + Qf v_f))/v_0,   [tQ >= ||Mh^-1 [Q_ZZ; 0]||, ou t_V nQZ]  maxwV (qTT v_T + Qf v_f)/v_T,  Qf ),
com qTT = max_J ||Q0 P_J||/min w, qZT = (o mesmo nas janelas dos modos de Z), maxwV = max w_J ||V_J||.
NK: Y + theta r + Z2 r^2 <= r e theta + 2 Z2 r < 1 => zero unico, DF(x*) invertivel, lambda* simples."""
from __future__ import annotations

import argparse
from fractions import Fraction as Fr
import json
import math
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from verified_norms import cota_norma, U

ARRED = 1 + Fr(1, 10**9)
Q = lambda x: Fr(float(x))*ARRED


def raiz_sup(x: Fr) -> Fr:
    s = Fr(math.sqrt(float(x)))*(1 + Fr(1, 10**12))
    while s*s < x:
        s *= 1 + Fr(1, 10**9)
    return s


def norma2_cota(M):
    return Q(cota_norma(M, 4*U*float(np.linalg.norm(M))))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--dir', type=Path, default=Path('build/s3b/rig'))
    a = ap.parse_args()
    z = json.loads((a.dir/'Z.json').read_text()); t = json.loads((a.dir/'T3.json').read_text())
    qj = json.loads((a.dir/'Q.json').read_text())               # nk_L3_q.py: normas restritas de Q0
    # rastreabilidade (revisao de S4, L5): um so lambda0 e um so erro do fundo em todas as etapas
    for nome, d in (('T3.json', t), ('Q.json', qj)):
        l = d.get('lam0', z['lam0']); l = complex(*l) if isinstance(l, list) else complex(l)
        lz = z['lam0']; lz = complex(*lz) if isinstance(lz, list) else complex(lz)
        assert l == lz, f'{nome}: lam0 difere do de Z.json'
    assert t['dB_L']['rt'] == z['dB_L']['rt'] and t['eps_mu'] == z['eps_mu'], 'T3.json e Z.json: fundo ou mu diferentes'
    # L5 (revisao): a_extra com o MAXIMO calculado de ||Q0|| nos modos MZ + banda <= |m| < MZ + 40
    tV_ = z['tV']; a_extra = tV_*z['dB_L']['resto']*qj['nQx_extra']
    aZ_rig = (math.sqrt(z['a_loewner']**2 + a_extra**2)
              + tV_*(z['erro_Y'] + (z['dB_L']['rt'] + z['eps_mu'])*z['nQT']))*(1 + 1e-12)
    # L4: linha l Q0 P_far (descida de n >= Nc ate Z) em Z'<-far
    epsZ = z['epsZ'] + tV_*qj['qZf']*1.0001
    N = np.array(t['N']); TZ = np.array(t['TZ']); nV = np.array(t['nV']); nVhi = np.array(t['nVhi'])
    Nf = np.array(t['Nf']); Y = np.array(t['Y']); nQ = np.array(t['nQ']); SZ = t['SZ']; k = len(N)
    MZ = z['Z'][0]
    beta, desc, pert = Q(t['beta']), Q(t['desc']), Q(t['pert'])
    # perturbacao GLOBAL (erro dos dados de RT e mu): nao respeita a banda, acopla tudo; as etapas
    # A e B a somaram so dentro da banda: aqui ela entra de novo, conservadora, em toda entrada
    eg = Q(t['dB_L']['rt']) + Q(t['eps_mu'])
    nQmax = Q(max(nQ))
    # pesos: sqrt(u/z) dos vetores de Perron de N (so propoem)
    lam, Vr = np.linalg.eig(N); zv = np.abs(Vr[:, np.argmax(np.abs(lam))])
    lam, Vl = np.linalg.eig(N.T); uv = np.abs(Vl[:, np.argmax(np.abs(lam))])
    melhor = None
    for nomew, w in (('1', np.ones(k)), ('sqrt(u/z)', np.sqrt(np.maximum(uv, 1e-9)/np.maximum(zv, 1e-9)))):
        w = w/w.max()
        wq = [Fr(float(x)) for x in w]
        nT = norma2_cota((w[:, None]*N)/w[None, :])
        TZw = raiz_sup(sum((wq[i]*Q(TZ[i]))**2 for i in range(k)))
        wmin = min(wq[j] for j in SZ)
        ZT = Q(aZ_rig)/wmin
        # L1 (revisao): o far de Fourier (|m| >= 200, todo n) alcanca pela banda |d| <= 40 TODAS as
        # linhas das janelas com |m| >= 160: nelas vale ||V_J||, nao ||V_J P_alto||
        borda = [max(abs(m) for m in t['modos'][i]) >= 160 for i in range(k)]
        Tf = max(wq[i]*Q(nV[i] if borda[i] else nVhi[i]) for i in range(k))*beta \
            + raiz_sup(sum((wq[i]*Q(nV[i]))**2 for i in range(k)))*desc
        fT = norma2_cota(Nf/w[None, :])
        # L2 (revisao): linhas do far em 216 <= |m| < 240 recebem das janelas |m| >= 176 pelo resto
        nQb = max(Q(nQ[i]) for i in range(k) if max(abs(m) for m in t['modos'][i]) >= 176)
        fT += Q(t['dB_L']['resto'])*nQb/min(wq)
        maxwV0 = max(wq[i]*Q(nV[i]) for i in range(k)); wmin0 = min(wq)
        nT += maxwV0*eg*nQmax/wmin0; TZw += maxwV0*eg*Q(z['nQZ']); Tf += maxwV0*eg*Q(z['Qfar'])
        ZT += Q(z['tV'])*eg*nQmax/wmin0; fT += eg*nQmax/wmin0
        M = [[Q(z['rhoZ']), ZT, Q(epsZ)], [TZw, nT, Tf], [pert*Q(z['nQZ']), fT, beta]]
        Mf = np.array([[float(x) for x in l] for l in M])
        lamM, vec = np.linalg.eig(Mf); v = np.abs(vec[:, np.argmax(np.abs(lamM))]); v = np.maximum(v/v.max(), 1e-9)
        vq = [Fr(float(x)) for x in v]
        th = max(sum(M[i][j]*vq[j] for j in range(3))/vq[i] for i in range(3))
        print(f'pesos {nomew}: T<-T {float(nT):.4f}, T<-Z\' {float(TZw):.4f}, Z\'<-T {float(ZT):.4f}, T<-far {float(Tf):.4f},'
              f' far<-T {float(fT):.4f}, beta {float(beta):.4f}  ->  theta = {float(th):.6f}')
        if melhor is None or th < melhor[0]:
            melhor = (th, M, vq, wq, nomew)
    th, M, vq, wq, nomew = melhor
    assert th < 1, f'theta = {float(th):.6f} >= 1'
    theta = th*(1 + Fr(1, 10**9))
    for i in range(3):
        assert sum(M[i][j]*vq[j] for j in range(3)) <= theta*vq[i]
    print(f'PERRON (3 niveis, pesos {nomew}) verificado em racionais: theta = {float(theta):.6f}')
    # 2a ordem
    wmin = min(wq); tV = Q(z['tV']); Qf = Q(z['Qfar']); nQZ = Q(z['nQZ'])
    # normas RESTRITAS de Q0 (nk_L3_q.py): ||P_Z Q0 P_T|| (descida da cauda para Z), ||P_T Q0 P_T||,
    # ||P_Z Q0 P_far|| (descida de n >= Nc); Q0 e diagonal em m e triangular superior em n
    qTT = Q(qj['qTT'])/wmin; qZT = Q(qj['qZT'])/wmin; qZf = Q(qj['qZf'])
    maxwV = max(wq[i]*Q(nV[i]) for i in range(k))
    v0, vT, vf = vq
    # parte de Z': ||Mh^-1 [P_Z Q0 v; 0]|| <= tQ ||v_Z|| + t_V (qZT ||v_T|| + qZf ||v_f||), tQ >= ||Mh^-1 [Q_ZZ; 0]||
    # (nk_L2_Z --so-tQ; sem ele, tQ = t_V ||Q_ZZ||)
    tQ = Q(z['tQ']) if 'tQ' in z else tV*nQZ
    Z2 = 2*v0*max((tQ*v0 + tV*(qZT*vT + qZf*vf))/v0, maxwV*(qTT*vT + Qf*vf)/vT, Qf)*ARRED
    # residuo
    YT = raiz_sup(sum((wq[i]*Q(Y[i]))**2 for i in range(k)))
    # L3 (revisao): o erro de RT em B h0 nao respeita a banda: soma em T (todas as janelas) e no far
    rtY = Q(t['dB_L']['rt'])*ARRED
    YT += maxwV*rtY
    Ymax = max(Q(z['YZ'])/v0, YT/vT, rtY/vf)
    d = (1 - theta)**2 - 4*Z2*Ymax
    assert d > 0, 'NK: discriminante <= 0'
    r = Fr((float(1 - theta) - math.sqrt(float(d)))/(2*float(Z2)))*(1 + Fr(1, 10**6))
    assert Ymax + theta*r + Z2*r*r <= r, 'NK: a bola nao e invariante'
    assert theta + 2*Z2*r < 1, 'NK: nao contrai'
    print(f'Z2 = {float(Z2):.3e}; Y = {float(Ymax):.3e} (Z\' {float(Q(z["YZ"])/v0):.2e}, T {float(YT/vT):.2e});'
          f' raio r = {float(r):.3e}')
    print(f'NK VERIFICADO: zero unico (y*, lambda*) com ||x* - x0|| <= r; |lambda* - {z["lam0"]}| <= {float(v0*r):.3e};'
          f' DF(x*) invertivel: lambda* algebricamente simples.')
    out = dict(theta=str(theta), v=[str(x) for x in vq], pesos=nomew, w=[str(x) for x in wq], M=[[str(x) for x in l] for l in M],
               Z2=str(Z2), Y=str(Ymax), r=str(r), lambda0=z['lam0'], erro_lambda=str(v0*r),
               erro_yZ=str(v0*r), erro_yT_W=str(vT*r), erro_yT_l2=str(vT*r/wmin), erro_yfar=str(vf*r), certificado=True)
    (a.dir/'C3.json').write_text(json.dumps(out, indent=1) + '\n')


if __name__ == '__main__':
    main()
