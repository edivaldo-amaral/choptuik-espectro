#!/usr/bin/env python3
"""R5 (revisor): remontagem da matriz 3x3 (Z', T, far) e de theta em RACIONAIS para 3 discos da faixa B.

Codigo proprio para a montagem (a partir das formulas, derivadas de novo: escala D_Z = J_Z(c) Q_ZZ(s) no nivel Z,
logo Z'<-Z' = J_Z(c) R(s) dB Q_ZZ(c), Z'<-T = J_Z(c) R(s)[...] etc.) e para o Perron: theta <= max_i (M v)_i / v_i
com v racional (Collatz-Wielandt, exato para matriz >= 0), e ||X||_2^2 <= max_i (X^T X u)_i/u_i com u racional
(X >= 0 entrada a entrada) em vez de cota_norma. Reaproveita do autor so as constantes fechadas (eps_mu_L, dB_L,
descida, descida_divxi, Qfar_disco) e as normas por modo das janelas centrais (dados_centrais, qz_ponto), que sao
normas de matrizes pequenas ja revisadas na faixa A. Compara entrada a entrada com det['M'] do JSON do disco."""
import json, math, sys
from fractions import Fraction as Fr
from pathlib import Path
import numpy as np
RAIZ = Path('<repo>'); sys.path.insert(0, str(RAIZ/'scripts'))
import closed_form as cf, fundo_L, signed_operator_L as sl
from rouche_L_disco import dados_centrais, qz_ponto, Qfar_disco

MU = 722873400/2**32; K2 = 5/4; NORMA_BL = 3.964043; BAND_N, BAND_M = 100, 40
UP = 1 + Fr(1, 10**12)
def q(x):                      # racional >= x (folga 1e-12)
    return Fr(float(x))*UP

def cw_norma2(X):
    """cota racional de ||X||_2 para X >= 0: sqrt(max_i (X^T X u)_i/u_i), u = vetor de Perron (float -> racional)."""
    Xq = [[q(x) for x in row] for row in X]
    G = np.array(X, float).T @ np.array(X, float)
    lam, V = np.linalg.eigh(G); u = np.maximum(np.abs(V[:, -1]), 1e-6); u = u/u.max()
    uq = [Fr(float(x)) for x in u]
    m, n = len(Xq), len(Xq[0])
    Xu = [sum(Xq[i][j]*uq[j] for j in range(n)) for i in range(m)]
    GTu = [sum(Xq[i][j]*Xu[i] for i in range(m)) for j in range(n)]
    r2 = max(GTu[j]/uq[j] for j in range(n))
    s = Fr(math.sqrt(float(r2)))*(1 + Fr(1, 10**9))
    while s*s < r2:
        s *= 1 + Fr(1, 10**9)
    return s

def sqrt_sup(x):
    s = Fr(math.sqrt(float(x)))*(1 + Fr(1, 10**9))
    while s*s < x:
        s *= 1 + Fr(1, 10**9)
    return s

def theta_proprio(p, r, zp, t3, c, cF, schur, eE, dtauB=1.4563e-09):
    MZ, NZ = schur['Z']; Mc, Nc = schur['F']
    dfr = schur['defl_refA']
    U = 2.0**-53
    errZ = [dtauB + 1e-60 + U*dfr['norma_x'][0]*10, eE + U*dfr['norma_x'][1]*10]
    dEZ = sum(dk*nphi*e for dk, nphi, e in zip(dfr['delta'], dfr['norma_phi'], errZ))
    dET = sum(dk*nphi*(cx + e) for dk, nphi, cx, e in zip(dfr['delta'], dfr['norma_phi'], dfr['cauda_x'], errZ))
    nE = sum(dk*nphi*nx for dk, nphi, nx in zip(dfr['delta'], dfr['norma_phi'], dfr['norma_x']))
    cwt = math.sqrt(1 + K2**-4)
    nQx = (1 + 2*cwt/(K2 - 1))/MU*1.0001            # ||P_Z Q0(s) P_T|| uniforme (Re s >= 0)
    fb40 = fundo_L.dB_L(sl.campos(), BAND_M)
    rt, arred = fb40['rt'], fb40['arred']
    fb_tot = t3['dB_L']['total']
    # escalares do ponto
    d = abs(p - c) + r
    tp, nQZp = zp['t_p'], zp['nQZ']
    assert r*tp < 0.9 and r*nQZp < 0.9
    ts = tp/(1 - r*tp)                               # ||R(s)||, Neumann a partir de p
    tau = 1 + (d + zp['nBhat'])*ts                    # ||J_Z(c) R(s)|| <= 1 + ||(B^ + s - c) R(s)||
    nQZs = nQZp/(1 - r*nQZp)
    JQ = 1 + d*nQZs                                   # ||J_Z(c) Q_ZZ(s)||
    qTT = cF['qTT']; assert d*qTT < 0.9
    emu_s = fundo_L.eps_mu_L(complex(abs(p.real) + r, abs(p.imag) + r))
    emu_c = t3['eps_mu']
    dB = arred + rt
    Qf_s = Qfar_disco(p, r, Mc, Nc)
    E = zp['a0'] + r*zp['b1'] + r*r*zp['b2'] + r**3*(tau*zp['nY3'] + JQ*zp['nQ3F'])
    EF = min(E*pr['phi'] + (tau + JQ)*pr['rho'] for pr in cF['pares'])
    ZT0 = (EF + tau*((dB + dEZ)*nQx + emu_s))/(1 - d*qTT)
    ZZ = tau*((dB + dEZ)*t3['nQZ'] + emu_c)
    Zf = tau*(NORMA_BL*(cf.descida(NZ + BAND_N + 1, Nc, MU, K2) + 6*cf.descida_divxi(NZ, Nc, K2))
              + (dB + dEZ)*Qf_s + emu_s + nE*cf.descida(NZ, Nc, MU, K2))
    # janelas
    N = np.array(t3['N']); TZ = np.array(t3['TZ']); nV = np.array(t3['nV']); nVhi = np.array(t3['nVhi'])
    K = np.array(t3['K']); nQm = np.array(t3['nQm']); Nf = np.array(t3['Nf']); k = len(nV)
    cent = dados_centrais(t3, c, MZ, NZ, Nc); qz = qz_ponto(cent, p, NZ)
    nQJJ = np.array([cent[i]['nQJJ'] if i in cent else nQm[i] for i in range(k)])
    assert np.all(d*nQJJ < 0.9)
    fJ = 1/(1 - d*nQJJ)
    qZs = {i: qz[i]/(1 - r*nQZp)*(1 + d*nQJJ[i]*fJ[i]) for i in cent}
    nQs = np.array([(t3['nQ'][i] + (d*nQZs*cent[i]['qZc'] if i in cent else 0.0))*fJ[i] for i in range(k)])
    dpert = max(0.0, emu_s - emu_c)
    Ns = np.zeros((k, k))
    for i in range(k):
        sens = d*(TZ[i]*qZs[i] + (K[i] + TZ[i]*cent[i]['qZc'])*fJ[i]) if i in cent else d*K[i]*fJ[i]
        Ns[i, i] = N[i, i] + sens + nV[i]*dpert*max(1.0, nQs[i])
        for j in range(k):
            if j != i and N[i, j] != 0:
                Ns[i, j] = N[i, j] + d*((TZ[i]*qZs[j] if j in cent else 0.0) + N[i, j]*nQJJ[j]*fJ[j]) + nV[i]*dpert*max(1.0, nQs[j])
    TZs = TZ + nV*(dpert + dET)*t3['nQZ']
    Nfs = Nf*fJ[None, :]
    beta = (NORMA_BL + rt)*Qf_s + emu_s
    lim = Nc - BAND_N - t3.get('folga_far', 60)
    desc = (NORMA_BL + fb_tot)*(cf.descida(lim + BAND_N + 1, Nc, MU, K2) + 6*cf.descida_divxi(lim, Nc, K2)) + (fb_tot + emu_s)*Qf_s
    eg = rt + emu_s
    borda = [max(abs(m) for m in t3['modos'][i]) >= Mc - BAND_M for i in range(k)]
    nQb = max(nQs[i] for i in range(k) if max(abs(m) for m in t3['modos'][i]) >= Mc - 24)
    nQmax = max(1.0, float(nQs.max()))
    SZ = t3['SZ']
    # racionais (pesos 1)
    print(f'   ||Ns||_2 (SVD float) = {np.linalg.norm(Ns, 2):.7f}, CW = {float(cw_norma2(Ns)):.7f}; ||Nfs||_2 (SVD) = {np.linalg.norm(Nfs, 2):.7f}, CW = {float(cw_norma2(Nfs)):.7f}', flush=True)
    nT = cw_norma2(Ns) + q(max(nV))*q(eg)*q(nQmax) + q(max(nV))*q(dET*nQx)
    fT = cw_norma2(Nfs) + q(t3['dB_L']['resto'])*q(nQb) + q(eg)*q(nQmax) + q(dET*nQx)
    TZw = sqrt_sup(sum(q(x)**2 for x in TZs)) + q(max(nV))*q(eg)*q(t3['nQZ'])
    Tf = max(q(nV[i] if borda[i] else nVhi[i]) for i in range(k))*q(beta) + sqrt_sup(sum(q(x)**2 for x in nV))*q(desc) \
        + q(max(nV))*q(eg)*q(Qf_s) + q(max(nV))*q(dET*Qf_s)
    ZT = q(ZT0) + q(tau)*q(eg)*q(nQmax)
    fZ = q(fb_tot + emu_s)*q(t3['nQZ']) + q(dET)*q(t3['nQZ'])
    ff = q(beta) + q(dET*Qf_s)
    M = [[q(ZZ), ZT, q(Zf)], [TZw, nT, Tf], [fZ, fT, ff]]
    Mf = np.array([[float(x) for x in l] for l in M])
    lam, V = np.linalg.eig(Mf); v = np.abs(V[:, np.argmax(np.abs(lam))]); v = np.maximum(v/v.max(), 1e-9)
    vq = [Fr(float(x)) for x in v]
    th = max(sum(M[i][j]*vq[j] for j in range(3))/vq[i] for i in range(3))
    return th, Mf, dict(d=d, tau=tau, E=E, Qf_s=Qf_s, emu_s=emu_s, nQmax=nQmax)

casos = [('build/faixaB/discos/disco_0_5.json', (0.0, 0.502449)),
         ('build/faixaB/discos/disco_1_0.json', (1.0, 0.75)),
         ('build/faixaB/discos/disco_0_75.json', (0.0, 0.75))]
schur = json.load(open(RAIZ/'build/rouche_L_rig/z/lab2/schur.json'))
eE = json.load(open(RAIZ/'build/s4/rig/E.json'))['e']
saida = []
import os
if os.environ.get('R5_SO'):
    casos = [casos[int(os.environ['R5_SO'])]]
for arq, pp in casos:
    dd = json.load(open(RAIZ/arq))
    zr = json.load(open(RAIZ/dd['zrig'])); cF = json.load(open(RAIZ/dd['certF']))
    c = complex(*zr['c'])
    rot = {0.5j: 'build/faixaB/cauda/T3_0_5.json', 0.75j: 'build/faixaB/cauda/T3_0_75.json', 1.0: 'build/rouche_L_rig/cauda/T3_1_0.json'}[c]
    t3 = json.load(open(RAIZ/rot)); assert complex(*t3['lam0']) == c
    assert zr['sha256_F'] == cF['sha256_F']
    x = [x for x in dd['pontos'] if tuple(x['p']) == pp][0]
    zp = [z for z in zr['pontos'] if tuple(z['p']) == pp][0]
    p = complex(*pp); r = x['r']
    th, Mf, aux = theta_proprio(p, r, zp, t3, c, cF, schur, eE)
    Ma = np.array(x['det']['M'])
    rel = np.abs(Mf - Ma)/np.maximum(np.abs(Ma), 1e-300)
    print(f'{arq} p = {p}: r = {r:.6e}, d = {aux["d"]:.4f}, tau = {aux["tau"]:.2f}')
    print('   M (revisor)  = ' + np.array2string(Mf, precision=6))
    print('   M (autor)    = ' + np.array2string(Ma, precision=6))
    print(f'   maior diferenca relativa por entrada {rel.max():.2e}; theta revisor (racional) = {float(th):.7f} < 1: {th < 1};'
          f' theta autor = {x["theta_f"]:.7f}')
    saida.append(dict(disco=arq, p=pp, r=r, theta_revisor=str(th), theta_revisor_f=float(th), theta_autor=x['theta_f'],
                      M_revisor=Mf.tolist(), M_autor=Ma.tolist(), dif_rel_max=float(rel.max()), ok=bool(th < 1)))
if not os.environ.get('R5_SO'):
    json.dump(saida, open(RAIZ/'.codex-runs/2026-10-02-revisao-B/r5_resultado.json', 'w'), indent=1)
