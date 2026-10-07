#!/usr/bin/env python3
"""Rouché de L: os blocos de janela da referência não têm zeros em Omega (C1_REAVALIACAO.md §9).

A referência do Rouché é A(s) = diag(H^_ZZ(s), H^_J(s) nas janelas, I), com H^_J(s) = P_J H(s) P_J e
H(s) = L~(s) Q0(s). O Perron no contorno dá H_t = A + t(H - A) invertível no contorno para t em [0, 1],
logo N(H, Omega) = N(A, Omega) = N(H^_ZZ, Omega) + soma_J N(H^_J, Omega). Este script certifica
N(H^_J, Omega) = 0 para toda janela J: H^_J(s) invertível em TODO s de Omega, e não só no contorno.

Partição de Omega = [0, 3] x [-1/4, 1/4] em retângulos Omega_k, um por centro de cauda c_k. Para s em
Omega_k, F_k(s) = I - V_J(c_k) H^_J(s) é analítica (Q0(s) é analítica em Re s > -mu: J_livre é
triangular com diagonal mu (n + 1) + s + i m/2). Pelo princípio do máximo (||F|| é subharmônica),
sup_{Omega_k} ||F_k|| <= max_{bordo de Omega_k} ||F_k||, e ||F_k|| < 1 dá H^_J invertível.

No bordo, com d = |s - c| e as identidades triangulares (Q0(s) P_J = P_Z Q_ZJ(s) + P_J Q_JJ(s),
Q0(c) P_Z = P_Z Q_ZZ(c), (I + (s - c) Q_JJ(c))^-1 = I - (s - c) Q_JJ(s)):
  V_J P_J B (Q0(s) - Q0(c)) P_J = -(s - c) [V_J P_J B P_Z Q_ZZ(c) Q_ZJ(s)
                                          + (V_J P_J B Q0(c)^2 P_J - V_J P_J B P_Z Q_ZZ(c) Q_ZJ(c)) (I - (s - c) Q_JJ(s))],
  janela central:     rho_J(s) <= rho_J(c) + d [TZ_J q_ZJ(s) + (K_J + TZ_J q_ZJ(c)) (1 + d q_JJ(s))],
  janela não central: rho_J(s) <= rho_J(c) + d K_J/(1 - d nQ_J(c)),
mais a perturbação nV_J [rt ||Q0(s) P_J|| + (eps_mu(s) + max(0, eps_mu(s) - eps_mu(c))) max(1, ||Q0(s) P_J||)]
(B exato contra o calculado; eps_mu_L já inclui Q0(s)) e, nas janelas centrais, a deflação verdadeira fora de Z,
nV_J dET q_ZJ(s) (revisão independente de 01/10, L1), como em rouche_L_disco. q_XY(s) = max por modo de ||Q_XY,m(s)||, cotados num disco |s - q| <= r
em torno de pontos de amostra q por Neumann estruturado (bloco triangular):
  Q_JJ(s) = Q_JJ(q) (I + (s - q) Q_JJ(q))^-1,  Q_ZJ(s) = (I + (s - q) Q_ZZ(q))^-1 Q_ZJ(q) (I + (s - q) Q_JJ(q))^-1.
A soma final de cada janela é feita em racionais."""
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
import fundo_L
import nk_L as nk
import signed_operator_L as sl
import closed_form as cf
from nk_L2_Z import err_mm
from nk_L3_C import Q
from verified_norms import U as UNIT, cota_norma

K2 = nk.K2

# (centro, rótulo do T3 em --caudas, Re0, Re1, Im0, Im1): cobre [0, 3] x [-1/4, 1/4]
REGIOES = [
    (0.25j, '0_25', 0.0, 0.3, 0.15, 0.25),
    (0j, '0', 0.0, 0.3, -0.15, 0.15),
    (-0.25j, '-0_25', 0.0, 0.3, -0.25, -0.15),
    (0.55 + 0.25j, '0_55p0_25', 0.3, 0.8, 0.0, 0.25),
    (0.55 - 0.25j, '0_55-0_25', 0.3, 0.8, -0.25, 0.0),
    (1.0, '1_0', 0.8, 1.2, -0.25, 0.25),
    (1.5, '1_5', 1.2, 1.8, -0.25, 0.25),
    (2.0, '2_0', 1.8, 2.2, -0.25, 0.25),
    (2.5, '2_5', 2.2, 2.8, -0.25, 0.25),
    (3.0, '3_0', 2.8, 3.0, -0.25, 0.25),
]

# faixa B: [0, 3] x [1/4, 3/4] (centros novos em --caudas-B, os da faixa A em --caudas); preenchida em REGIOES_B
REGIOES_B = [
    (0.25j, '0_25', 0.0, 0.3, 0.25, 0.375),
    (0.5j, '0_5', 0.0, 0.3, 0.375, 0.625),
    (0.75j, '0_75', 0.0, 0.3, 0.625, 0.75),
    (0.55 + 0.25j, '0_55p0_25', 0.3, 0.8, 0.25, 0.5),
    (0.55 + 0.75j, '0_55p0_75', 0.3, 0.8, 0.5, 0.75),
    (1.0, '1_0', 0.8, 1.2, 0.25, 0.75),
    (1.5, '1_5', 1.2, 1.8, 0.25, 0.75),
    (2.0, '2_0', 1.8, 2.2, 0.25, 0.75),
    (2.5, '2_5', 2.2, 2.8, 0.25, 0.75),
    (3.0, '3_0', 2.8, 3.0, 0.25, 0.75),
]


def amostras(re0, re1, im0, im1):
    """Pontos médios de subsegmentos do bordo do retângulo e o raio que cobre cada subsegmento.
    Passo 0,1 em Re <= 0,8 e 0,2 à direita."""
    out = []
    for a, b in (((re0, im0), (re1, im0)), ((re1, im0), (re1, im1)), ((re1, im1), (re0, im1)), ((re0, im1), (re0, im0))):
        za, zb = complex(*a), complex(*b)
        L = abs(zb - za)
        h = 0.1 if max(a[0], b[0]) <= 0.8 + 1e-12 else 0.2
        n = max(1, round(L/h))
        for k in range(n):
            q = za + (zb - za)*(k + 0.5)/n
            out.append((complex(round(q.real, 9), round(q.imag, 9)), L/n/2*(1 + 1e-9)))
    return out


def Q_rapido(m, N, s):
    """Como nk_L2_Z.Q_cert (mesma matriz pesada), mas com cotas de Frobenius (||X||_2 <= ||X||_F) no
    resíduo e em ||Q||: o resíduo é ~1e-13, e a cota espectral certificada custava 10x mais."""
    r = sl.Radial(m, N)
    w = np.sqrt(np.concatenate([cf.omega(ns, K2) for _, ns in r.comps]))
    Jw = (w[:, None]*sl.J_livre(r, s))/w[None, :]
    Q = np.linalg.inv(Jw)
    E = np.eye(len(Q)) - Jw @ Q
    fE = float(np.linalg.norm(E))
    e = fE*(1 + 1e-6) + err_mm(Jw, Q) + UNIT*fE
    assert e < 0.5
    nQ = float(np.linalg.norm(Q))*(1 + 1e-6)
    return Q, nQ*e/(1 - e)


class Blocos:
    """Por ponto q: por modo central m, cotas certificadas de ||Q_ZZ,m(q)||, ||Q_ZJ,m(q)||, ||Q_JJ,m(q)||."""
    def __init__(self, MZ, NZ, Nc):
        self.MZ, self.NZ, self.Nc = MZ, NZ, Nc
        self.cache = {}
        self.idx = {}

    def __call__(self, q):
        if q in self.cache:
            return self.cache[q]
        res = {}
        for m in range(-(self.MZ - 1), self.MZ):
            Qm, dq = Q_rapido(m, self.Nc, q)
            if m not in self.idx:
                r = sl.Radial(m, self.Nc); nn = np.concatenate([ns for _, ns in r.comps])
                self.idx[m] = (np.where(nn < self.NZ)[0], np.where(nn >= self.NZ)[0])
            iz, iT = self.idx[m]
            res[m] = (cota_norma(Qm[np.ix_(iz, iz)], dq), cota_norma(Qm[np.ix_(iz, iT)], dq),
                      cota_norma(Qm[np.ix_(iT, iT)], dq))
        self.cache[q] = res
        return res


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--caudas', type=Path, default=Path('build/rouche_L_rig/cauda'))
    ap.add_argument('--Z', default='32x128'); ap.add_argument('--F', default='200x400')
    ap.add_argument('--regioes', default=None, help='índices das regiões (padrão: todas)')
    ap.add_argument('--faixa', default='A', choices=['A', 'B'])
    ap.add_argument('--caudas-B', dest='caudas_B', type=Path, default=Path('build/faixaB/cauda'))
    ap.add_argument('--saida', type=Path, default=Path('build/rouche_L_rig/janelas.json'))
    ap.add_argument('--schur', type=Path, default=Path('build/rouche_L_rig/z/lab2/schur.json'), help='defl_refA')
    ap.add_argument('--E', default='build/s4/rig/E.json', help='etapa E de S4: ||g* - g_A||')
    ap.add_argument('--dtauB', type=float, default=1.4563e-09, help='||d_tau RefB||_L (gauge_refB)')
    a = ap.parse_args()
    # deflacao verdadeira fora de Z (revisao de 01/10, L1), como em rouche_L_disco: na janela central J entra
    # delta_k P_J x_k phi_k^* P_Z Q_ZJ(s), norma <= dET q_ZJ(s); nas nao centrais Q_ZJ = 0
    dfr = json.load(open(a.schur))['defl_refA']
    eE = json.load(open(a.E))['e']
    errZ = [a.dtauB + 1e-60 + UNIT*dfr['norma_x'][0]*10, eE + UNIT*dfr['norma_x'][1]*10]
    dET = sum(dk*nphi*(cx + e) for dk, nphi, cx, e in zip(dfr['delta'], dfr['norma_phi'], dfr['cauda_x'], errZ))
    t0 = time.time()
    MZ, NZ = map(int, a.Z.split('x')); Mc, Nc = map(int, a.F.split('x'))
    rt = fundo_L.dB_L(sl.campos(), nk.BAND_M)['rt']
    blocos = Blocos(MZ, NZ, Nc)
    regs = REGIOES if a.faixa == 'A' else REGIOES_B
    idx = range(len(regs)) if a.regioes is None else [int(x) for x in a.regioes.split(',')]
    saida = json.loads(a.saida.read_text()) if a.saida.exists() else {}
    for k in idx:
        c, rot, re0, re1, im0, im1 = regs[k]
        c = complex(c)
        arq = a.caudas/f'T3_{rot}.json'
        if not arq.exists():
            arq = a.caudas_B/f'T3_{rot}.json'
        t3 = json.load(open(arq))
        assert complex(*t3['lam0']) == c, f'centro da cauda {t3["lam0"]} != {c}'
        N = np.array(t3['N']); rho = np.diag(N); TZ = t3['TZ']; K = t3['K']; nQm = t3['nQm']; nV = t3['nV']
        nj = len(K)
        cent = {i: [m for m in t3['modos'][i] if abs(m) < MZ] for i in range(nj) if min(abs(m) for m in t3['modos'][i]) < MZ}
        bc = blocos(c)
        qZc = {i: max(bc[m][1] for m in ms) for i, ms in cent.items()}
        emu_c = t3['eps_mu']
        pior = (Fr(0), None, None)
        for q, r in amostras(re0, re1, im0, im1):
            b = blocos(q)
            d = abs(q - c) + r
            emu_s = fundo_L.eps_mu_L(complex(abs(q.real) + r, abs(q.imag) + r))
            emu_t = emu_s + max(0.0, emu_s - emu_c)        # eps_mu ja inclui Q0(s): fator max(1, ||Q0(s) P_J||)
            for i in range(nj):
                if i in cent:
                    qzj = qjj = 0.0
                    for m in cent[i]:
                        qzz_m, qzj_m, qjj_m = b[m]
                        assert r*qzz_m < 0.9 and r*qjj_m < 0.9, (q, r, m, qzz_m, qjj_m)
                        qjj = max(qjj, qjj_m/(1 - r*qjj_m))
                        qzj = max(qzj, qzj_m/((1 - r*qzz_m)*(1 - r*qjj_m)))
                    nQs = Q(qzj) + Q(qjj)
                    v = Q(rho[i]) + Q(d)*(Q(TZ[i])*Q(qzj) + (Q(K[i]) + Q(TZ[i])*Q(qZc[i]))*(1 + Q(d)*Q(qjj))) \
                        + Q(nV[i])*(Q(rt)*nQs + Q(emu_t)*max(Fr(1), nQs) + Q(dET)*Q(qzj))
                else:
                    assert d*nQm[i] < 0.9
                    fq = 1/(1 - Q(d)*Q(nQm[i]))
                    nQs = Q(nQm[i])*fq
                    v = Q(rho[i]) + Q(d)*Q(K[i])*fq + Q(nV[i])*(Q(rt)*nQs + Q(emu_t)*max(Fr(1), nQs))
                if v > pior[0]:
                    pior = (v, q, i)
        v, q, i = pior
        ok = v < 1
        m0, m1 = t3['modos'][i][0], t3['modos'][i][-1]
        print(f'regiao {k} (centro {c}, [{re0}, {re1}] x [{im0}, {im1}]): max rho_J no bordo <= {float(v):.4f} '
              f'(janela {m0}..{m1}, s = {q})  {"OK" if ok else "FALHA"}  [{time.time() - t0:.0f}s, {len(blocos.cache)} pontos]',
              flush=True)
        saida[str(k)] = dict(centro=[c.real, c.imag], retangulo=[re0, re1, im0, im1], max_rho=str(v), max_rho_f=float(v),
                             janela=[m0, m1], s=[q.real, q.imag], ok=bool(ok))
        a.saida.write_text(json.dumps(saida, indent=1) + '\n')
    print('todas as janelas sem zeros em Omega' if all(x['ok'] for x in saida.values()) and len(saida) == len(regs)
          else 'incompleto ou com falha')


if __name__ == '__main__':
    main()
