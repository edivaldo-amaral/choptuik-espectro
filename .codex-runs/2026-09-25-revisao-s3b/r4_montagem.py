#!/usr/bin/env python3
"""R4 (revisao independente de S3b): remontagem PROPRIA do certificado NK de L (3 niveis) a partir
de Z.json e T3.json, com as desigualdades derivadas pelo revisor, e o efeito das lacunas achadas:

  L1  T <- far: a cota do codigo usa ||V_J P_alto|| x beta nas linhas altas e so a DESCIDA radial nas
      baixas; mas as colunas do far de Fourier (|m| >= Mc = 200, TODO n) alcancam pela banda |d| <= 40
      as linhas BAIXAS das janelas com |m| >= 160. Conserto: nessas janelas usar ||V_J|| no lugar de
      ||V_J P_alto||.
  L2  far <- T: os grupos de linhas do far cobrem |m| < 216; as linhas do far em 216 <= |m| < 240 sao
      alcancadas pelo resto da banda (17 <= |d| <= 40) das janelas com |m| >= 176 e nao entram.
      Conserto: grupos extras com entradas resto x ||Q0 P_J||.
  L3  Y: falta a componente far do residuo (erro de RT e mu em h0: ||P_far dB h0|| <= rt ||h0||).
  L4  Z' <- far: falta a linha de bordejamento l Q0 P_far (descida de n >= 400 para n < 128).
  L5  a_extra: ||Q0 P_(48 <= |m| < 72)|| tomado no modo 48 ('decresce com |m|'): conferido aqui.
Tambem: a sensibilidade (quanto theta, Y, Z2 podem crescer)."""
import json
import math
import sys
from fractions import Fraction as Fr
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path('scripts').resolve()))
import closed_form as cf
import signed_operator_L as sl
from verified_norms import cota_norma, U

MU, K2 = sl.MU, sl.K2


def monta(z, t, correcoes=False, log=print):
    N = np.array(t['N']); TZ = np.array(t['TZ']); nV = np.array(t['nV']); nVhi = np.array(t['nVhi'])
    Nf = np.array(t['Nf']); Y = np.array(t['Y']); nQ = np.array(t['nQ']); SZ = t['SZ']; k = len(N)
    w = np.ones(k)                                  # pesos '1' (os do certificado)
    beta, desc, pert = t['beta'], t['desc'], t['pert']
    eg = t['dB_L']['rt'] + t['eps_mu']
    nQmax = max(nQ)
    nT = cota_norma(N, 4*U*np.linalg.norm(N))
    TZw = math.sqrt(sum(TZ**2))
    ZT = z['aZ']
    # L1: janelas a menos de 40 modos do far de Fourier: ||V_J|| (linhas baixas alcancadas)
    borda = [i for i, ms in enumerate(t['modos']) if max(abs(m) for m in ms) >= 200 - 40]
    fator_T = [nV[i] if (correcoes and i in borda) else nVhi[i] for i in range(k)]
    Tf = max(fator_T)*beta + math.sqrt(sum(nV**2))*desc
    fT_M = Nf.copy()
    if correcoes:                                   # L2: grupos 216..239 e -239..-216
        extra = np.zeros((2, k))
        for i, ms in enumerate(t['modos']):
            if max(ms) >= 216 - 40:
                extra[0, i] = t['dB_L']['resto']*nQ[i]
            if min(ms) <= -(216 - 40):
                extra[1, i] = t['dB_L']['resto']*nQ[i]
        fT_M = np.vstack([Nf, extra])
    fT = cota_norma(fT_M, 4*U*np.linalg.norm(fT_M))
    maxwV = max(nV)
    nT += maxwV*eg*nQmax; TZw += maxwV*eg*z['nQZ']; Tf += maxwV*eg*z['Qfar']
    ZT += z['tV']*eg*nQmax; fT += eg*nQmax
    epsZ = z['epsZ']
    if correcoes:                                   # L4: l Q0 P_far, ||h0|| = 1
        epsZ += z['tV']*cf.descida(z['Z'][1], z['F'][1], MU, K2)
    M = np.array([[z['rhoZ'], ZT, epsZ], [TZw, nT, Tf], [pert*z['nQZ'], fT, beta]])*(1 + 1e-9)
    lam, vec = np.linalg.eig(M); v = np.abs(vec[:, np.argmax(np.abs(lam))]); v = v/v.max()
    Mq = [[Fr(float(x)) for x in l] for l in M]; vq = [Fr(float(x)) for x in v]
    th = max(sum(Mq[i][j]*vq[j] for j in range(3))/vq[i] for i in range(3))
    v0, vT, vf = (float(x) for x in vq)
    th = float(th)
    zw = [i for i, ms in enumerate(t['modos']) if min(abs(m) for m in ms) < z['Z'][0]]
    Z2 = 2*v0*max(z['tV']*(z['nQZ']*v0 + max(nQ[i] for i in zw)*vT + z['Qfar']*vf)/v0,
                  maxwV*(nQmax*vT + z['Qfar']*vf)/vT, z['Qfar'])
    YT = math.sqrt(sum(Y**2))
    Ys = [z['YZ']/v0, YT/vT]
    if correcoes:                                   # L3: residuo no far (RT e mu em h0) e janelas nao alcancadas
        from tile_perron import DELTA_MU
        rt = t['dB_L']['rt'] + DELTA_MU*(z['J1h'] + 2*cf.norma_divxi(K2))
        Ys.append(rt/vf)
        Ys[1] = (YT + maxwV*rt)/vT
    Ym = max(Ys)
    d = (1 - th)**2 - 4*Z2*Ym
    r = ((1 - th) - math.sqrt(d))/(2*Z2) if d > 0 else float('nan')
    log(f'  M = {np.array2string(M, precision=6)}')
    log(f'  v = {np.round([v0, vT, vf], 5)};  theta = {th:.6f};  Z2 = {Z2:.3f};  Y = {Ym:.4e} (partes: '
        + ', '.join(f'{x:.2e}' for x in Ys) + f');  discriminante = {d:.3e};  r = {r:.4e};  |lambda - lambda0| <= {v0*r:.3e}')
    return th, Z2, Ym, d


def main():
    rig = Path(sys.argv[1]) if len(sys.argv) > 1 else Path('build/s3b/rig')
    z = json.loads((rig/'Z.json').read_text()); t = json.loads((rig/'T3.json').read_text())
    print('Remontagem propria (pesos 1), sem correcoes (deve reproduzir C3.json):')
    th0, Z20, Y0, d0 = monta(z, t, False)
    print('Com as correcoes L1-L4:')
    th1, Z21, Y1, d1 = monta(z, t, True)
    print(f'  => o NK {"CONTINUA FECHANDO" if d1 > 0 else "NAO FECHA"} com as correcoes (theta {th0:.6f} -> {th1:.6f})')
    # sensibilidade
    print('Sensibilidade (sem correcoes): theta maximo tolerado = '
          f'{1 - math.sqrt(4*Z20*Y0):.6f} (folga {1 - math.sqrt(4*Z20*Y0) - th0:.2e}); Y e Z2 podem crescer '
          f'{((1 - th0)**2/(4*Z20*Y0) - 1)*100:.2f}% (produto)')
    print(f'  Y e dominado pelo erro de RT: t_V x 40 eps_omega = {z["tV"]*t["dB_L"]["rt"]:.4e} de Y_Z = {z["YZ"]:.4e}'
          f' ({z["tV"]*t["dB_L"]["rt"]/z["YZ"]*100:.2f}%); a constante 40 pode ir no maximo a '
          f'{40*(1 + ((1 - th0)**2/(4*Z20*Y0) - 1)*z["YZ"]/(z["tV"]*t["dB_L"]["rt"])):.2f}')
    # L5: ||Q_m|| para 48 <= |m| < 72
    nQx = z['a_extra']/(z['tV']*z['dB_L']['resto']*1.01)
    lam0 = z['lam0']; Nc = z['F'][1]
    pior = (0, 0.0)
    for m in list(range(48, 72)) + list(range(-71, -47)):
        r = sl.Radial(m, Nc)
        w = np.sqrt(np.concatenate([cf.omega(ns, K2) for _, ns in r.comps]))
        Q = np.linalg.solve(sl.J_livre(r, lam0), np.eye(r.dim))
        v = np.linalg.norm((w[:, None]*Q)/w[None, :], 2)
        if v > pior[1]:
            pior = (m, v)
    print(f'L5: max_(48 <= |m| < 72) ||Q_m|| = {pior[1]:.5f} (modo {pior[0]}); usado em a_extra: {nQx:.5f} -> '
          + ('OK' if pior[1] <= nQx else 'SUBESTIMADO'))
    # Qfar: secao finita de ||R_m^-1 P_(n >= 400)|| (cota inferior da verdadeira) contra Qfar
    pior = 0.0
    for m in (0, 1, -1, 2, 5, 20, 100, 199):
        r = sl.Radial(m, 1400)
        J = sl.J_livre(r, lam0)
        # por componente (J e bloco-diagonal nas componentes)
        for c, ns in r.comps:
            o0, o1 = r.off[c]
            Jc = J[o0:o1, o0:o1]
            lw = np.where(ns == 0, 0.0, ns*math.log(K2) + 0.5*np.log1p(K2**(-4.0*ns)))   # log sqrt(omega)
            Qc = np.linalg.solve(Jc, np.eye(len(ns)))
            cols = ns >= Nc
            # pesos relativos exp(lw_j - lw_c) sem estouro
            Wt = np.exp(lw[:, None] - lw[None, cols])
            v = np.linalg.norm(Qc[:, cols]*Wt, 2)
            pior = max(pior, v)
    print(f'Qfar: secao finita (n < 1400) max_m ||R_m^-1 P_(n>=400)|| ~ {pior:.5f}  <=  Qfar usado {z["Qfar"]:.5f}: '
          + ('OK' if pior <= z['Qfar'] else 'PROBLEMA'))


if __name__ == '__main__':
    main()
