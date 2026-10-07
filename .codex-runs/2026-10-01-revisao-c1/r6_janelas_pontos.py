#!/usr/bin/env python3
"""R6 (revisao independente C1): recalculo PROPRIO (ponto flutuante + Neumann; as constantes da cauda
rho, TZ, K, nV, qZc vem do T3, que sao dados certificados) da cota de ||I - V_J(c) H^_J(s)|| no bordo de
um retangulo Omega_k, em amostras q com raio r (disco que cobre o subsegmento), para a janela que deu o
maximo no relatorio do autor:
  v = rho_J + d [TZ_J q_ZJ + (K_J + TZ_J q_ZJ(c)) (1 + d q_JJ)] + nV_J pert (q_ZJ + q_JJ),  d = |q - c| + r,
  q_JJ = max_m ||Q_JJ,m(q)||/(1 - r ||Q_JJ,m(q)||),
  q_ZJ = max_m ||Q_ZJ,m(q)||/((1 - r ||Q_ZZ,m(q)||)(1 - r ||Q_JJ,m(q)||)),
  pert = rt + eps_mu(|Re q| + r, |Im q| + r) + max(0, eps_mu(s) - eps_mu(c)).
Os ||Q||: inversa de J_mu + q pesado (sqrt omega), modo a modo, norma espectral em ponto flutuante."""
import json, sys
import numpy as np
sys.path.insert(0, 'scripts')
import closed_form as cf
import fundo_L
import signed_operator_L as sl

MZ, NZ, Nc, K2 = 32, 128, 400, 5/4
rt = 4.7e-08
def blocos(m, s):
    r = sl.Radial(m, Nc); w = np.sqrt(np.concatenate([cf.omega(ns, K2) for _, ns in r.comps]))
    Q = np.linalg.inv((w[:, None]*sl.J_livre(r, s))/w[None, :])
    nn = np.concatenate([ns for _, ns in r.comps]); iz = np.where(nn < NZ)[0]; iT = np.where(nn >= NZ)[0]
    n2 = lambda X: float(np.linalg.norm(X, 2))
    return n2(Q[np.ix_(iz, iz)]), n2(Q[np.ix_(iz, iT)]), n2(Q[np.ix_(iT, iT)])
casos = [(0j, '0', 0.3 + 0.1j, 0.05, 0.1488), (0.25j, '0_25', 0.3 + 0.2j, 0.05, 0.1051),
         (0j, '0', 0.0 + 0.0j, 0.05, None)]
for c, rot, q, r, autor in casos:
    t3 = json.load(open(f'build/rouche_L_rig/cauda/T3_{rot}.json'))
    i = 0; modos = [m for m in t3['modos'][i] if abs(m) < MZ]          # janela -15..0 (central)
    rho = t3['N'][i][i]; TZ = t3['TZ'][i]; K = t3['K'][i]; nV = t3['nV'][i]
    qZc = max(blocos(m, c)[1] for m in modos)
    qjj = qzj = 0.0
    for m in modos:
        a, b, e = blocos(m, q)
        assert r*a < 0.9 and r*e < 0.9
        qjj = max(qjj, e/(1 - r*e)); qzj = max(qzj, b/((1 - r*a)*(1 - r*e)))
    d = abs(q - c) + r
    emu_c = t3['eps_mu']; emu_s = fundo_L.eps_mu_L(complex(abs(q.real) + r, abs(q.imag) + r))
    pert = rt + emu_s + max(0.0, emu_s - emu_c)
    v = rho + d*(TZ*qzj + (K + TZ*qZc)*(1 + d*qjj)) + nV*pert*(qzj + qjj)
    print(f'centro {c}, janela {t3["modos"][i][0]}..{t3["modos"][i][-1]}, q = {q}, r = {r}: rho {rho:.3e}, TZ {TZ:.4f}, K {K:.4f}, '
          f'qZJ(c) {qZc:.4f}, q_ZJ {qzj:.4f}, q_JJ {qjj:.4f}, d {d:.4f}  ->  cota {v:.4f}' + (f'  (autor: {autor})' if autor else ''))
