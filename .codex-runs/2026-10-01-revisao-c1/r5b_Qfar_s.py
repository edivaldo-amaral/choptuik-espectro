#!/usr/bin/env python3
"""R5b (revisao independente C1): a sensibilidade do far, Qf_s = Qf(c)/(1 - d Qf(c)) (rouche_L_disco), e uma
cota valida de ||Q0(s) P_far||?  Q0(s) P_far tem linhas de Z/T (Q0 e triangular SUPERIOR em n):
  P_far Q0(s) P_far = Q_ff(c) (I + (s - c) Q_ff(c))^-1                    (norma <= Qf/(1 - d Qf): ok)
  P_N   Q0(s) P_far = (I + (s - c) Q_NN(c))^-1 Q_Nf(c) (I + (s - c) Q_ff(c))^-1   (fator extra (I + (s-c)Q_NN(c))^-1!)
Teste numerico: ||Q0_m(s) P_{n >= Nc}|| (pesos sqrt omega, por componente; truncagem n < 1200 = cota INFERIOR
do valor verdadeiro) contra Qf_s, para o disco de pior d (centro 0,25i, p = 0,1258i) e outros.
Para comparacao, a cota do autor do far em c (T3 Qfar) e a cota fechada cauda_Rinv avaliada DIRETO em s."""
import json, math, sys
import numpy as np
from scipy.linalg import solve_triangular
sys.path.insert(0, 'scripts')
import closed_form as cf

MU = 722873400/2**32; K2 = 5/4; Nc = 400; NB = 1200
def far_norm(sig, passo):
    ns = np.arange(NB) if passo == 1 else np.arange(1, 2*NB, 2)
    k = len(ns)
    J = np.diag(MU*(ns + 1.0) + sig).astype(complex)
    for b in range(k):
        J[:b, b] += 2*MU*ns[b]
    cols = np.where(ns >= Nc)[0]
    E = np.zeros((k, len(cols)), complex); E[cols, np.arange(len(cols))] = 1
    X = solve_triangular(J, E)
    lw = np.where(ns == 0, 0.0, ns*math.log(K2) + 0.5*np.log1p(K2**(-4.0*ns)))
    Xw = X*np.exp(lw[:, None] - lw[cols][None, :])
    return float(np.linalg.norm(Xw, 2))
def modo(s, m):
    sig = s + 1j*m/2
    return max(far_norm(sig, 2), far_norm(sig, 1)) if m % 2 == 0 else far_norm(sig, 1)

casos = [('0_25', 0.25j, 'build/rouche_L_rig/discos/disco_0_25.json', 0), ('0', 0j, 'build/rouche_L_rig/discos/disco_0.json', 10),
         ('1_0', 1+0j, 'build/rouche_L_rig/discos/disco_1_0.json', 9)]
pior = 0
for rot, c, arq, ip in casos:
    t3 = json.load(open(f'build/rouche_L_rig/cauda/T3_{rot}.json'))
    dp = json.load(open(arq))['pontos'][ip]; p = complex(*dp['p']); r = dp['r']; d = abs(p - c) + r
    Qf = t3['Qfar']; Qfs = Qf/(1 - d*Qf)
    # pontos s do contorno no disco (extremos), e o modo com maior norma
    ss = [p]
    if p.real == 0: ss += [complex(0, p.imag + r*0.999), complex(0, p.imag - r*0.999)]
    if abs(abs(p.imag) - 0.25) < 1e-12: ss += [complex(p.real + 0.999*r, p.imag), complex(max(p.real - 0.999*r, 0), p.imag)]
    for s in ss:
        vals = {m: modo(s, m) for m in range(-6, 7)}
        mm = max(vals, key=vals.get)
        vc = modo(c, mm)
        fechada = max(cf.cauda_Rinv(s + 1j*mm/2, Nc, pp, MU, K2, 1500) for pp in ((2, 1) if mm % 2 == 0 else (1,)))
        print(f'c={c} p={p:.4f} r={r:.3e} d={d:.3f} s={s:.4f}: max_m ||Q0_m(s)P_far|| (trunc.) = {vals[mm]:.5f} (m={mm}; em c: {vc:.5f});'
              f' Qf_s do autor = {Qfs:.5f} (Qf(c) = {Qf:.5f}); cauda_Rinv fechada em s = {fechada:.5f}'
              + ('   *** VIOLA ***' if vals[mm] > Qfs else ''))
        pior = max(pior, vals[mm]/Qfs)
print(f'pior razao verdadeiro/Qf_s = {pior:.4f}')
