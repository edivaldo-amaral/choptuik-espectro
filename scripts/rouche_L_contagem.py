#!/usr/bin/env python3
"""Rouché de L, contagem finita: quantos zeros a referência de nível Z tem em Omega (C1_REAVALIACAO.md §9).

Os zeros de H^_ZZ(s) = L~_Z(s) J_Z(s)^-1 em Omega são os s com A + s singular (A = W L~_Z(0) W^-1, a
matriz de ponto flutuante tomada exata; J_Z(s) é invertível em Re s > -mu), isto é, s = -lambda com
lambda autovalor de A.

Subcomando `residuo` (na máquina com A.npy, T.npy e U.npy de rouche_L_Z schur): com a decomposição de
Schur em ponto flutuante A ~ U T U^* (T triangular superior, exata como dado),
  eps_R >= ||A U - U T||_F   (por blocos de colunas; arredondamento pela cota de Frobenius
                              |fl(X Y) - X Y| <= gamma_n |X||Y|, ||(|X||Y|)||_F <= ||X||_F ||Y||_F),
  eps_G >= ||I - U^* U||_F,  ||U^-1|| <= 1/sqrt(1 - eps_G),
  ||A - U T U^-1|| <= eps_R/sqrt(1 - eps_G) =: eps_B.
Os autovalores de B = U T U^-1 são EXATAMENTE a diagonal de T. Grava contagem.json com eps_B e os
-T_ii perto de Omega.

Subcomando `rouche` (local, com os discos de rouche_L_disco): no contorno, ||(A + s)^-1|| <= t(s) com
t(s) <= t_p/(1 - r t_p) em cada disco (p, r) da cobertura. Se max t(s) eps_B < 1, A + s + t (B - A) é
invertível no contorno para t em [0, 1] e A tem em Omega tantos autovalores (com multiplicidade)
quanto B: #{i : -T_ii em Omega}."""
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
from verified_norms import U as UNIT, gamma

RE0, RE1, IM0, IM1 = 0.0, 3.0, -0.25, 0.25         # faixa A (padrao); --regiao troca


def dentro(z):
    return RE0 < z.real < RE1 and IM0 < z.imag < IM1


def dist_bordo(z):
    """Distância de z ao bordo do retângulo (para z fora ou dentro)."""
    dx = max(RE0 - z.real, 0.0, z.real - RE1); dy = max(IM0 - z.imag, 0.0, z.imag - IM1)
    if dx > 0 or dy > 0:
        return math.hypot(dx, dy)
    return min(z.real - RE0, RE1 - z.real, z.imag - IM0, IM1 - z.imag)


def cmd_residuo(a):
    t0 = time.time()
    A = np.load(a.dir/'A.npy'); n = len(A)
    fA = float(np.linalg.norm(A))
    U = np.load(a.dir/'U.npy'); fU = float(np.linalg.norm(U))
    T = np.load(a.dir/'T.npy')
    assert np.all(np.tril(T, -1) == 0), 'T nao e triangular superior'
    lam = np.diag(T).copy()
    eR2 = 0.0; eG2 = 0.0; arR = 0.0; arG = 0.0
    b = a.bloco
    for j0 in range(0, n, b):
        j1 = min(n, j0 + b)
        Ub = U[:, j0:j1]; Tb = T[:j1, j0:j1]
        P1 = A @ Ub
        P2 = U[:, :j1] @ Tb
        R = P1 - P2
        eR2 += float(np.linalg.norm(R))**2
        fUb = float(np.linalg.norm(Ub)); fTb = float(np.linalg.norm(Tb))
        arR += 4*gamma(n)*fA*fUb + 4*gamma(j1)*float(np.linalg.norm(U[:, :j1]))*fTb \
            + 2*UNIT*(float(np.linalg.norm(P1)) + float(np.linalg.norm(P2)))
        del P1, P2, R
        G = U.conj().T @ Ub
        G[j0:j1, :] -= np.eye(j1 - j0)
        eG2 += float(np.linalg.norm(G))**2
        arG += 4*gamma(n)*fU*fUb + 2*UNIT*float(np.linalg.norm(G))
        del G
        print(f'  colunas {j0}..{j1}  [{time.time() - t0:.0f}s]', flush=True)
    # soma de cotas de erro por bloco: ||E||_F <= sum_b ||E_b||_F (grosseiro, mas folgado)
    eR = math.sqrt(eR2)*(1 + 1e-6) + arR
    eG = math.sqrt(eG2)*(1 + 1e-6) + arG
    assert eG < 0.5
    eB = eR/math.sqrt(1 - eG)*(1 + 1e-9)
    zs = -lam
    em = sorted([complex(z) for z in zs if dentro(z)], key=lambda z: z.real)
    perto = sorted([complex(z) for z in zs if not dentro(z)], key=dist_bordo)[:12]
    print(f'||A U - U T||_F <= {eR:.3e}, ||I - U^*U||_F <= {eG:.3e}, ||A - U T U^-1|| <= {eB:.3e}')
    print(f'-T_ii em Omega ({len(em)}): ' + ', '.join(f'{z.real:.6f}{z.imag:+.2e}i' for z in em))
    print('fora, mais perto do bordo: ' + ', '.join(f'{z.real:.4f}{z.imag:+.4f}i (dist {dist_bordo(z):.3e})' for z in perto))
    json.dump(dict(eps_R=eR, eps_G=eG, eps_B=eB, n=n, em_Omega=[[z.real, z.imag] for z in em],
                   fora_perto=[[z.real, z.imag, dist_bordo(z)] for z in perto]),
              open(a.dir/'contagem.json', 'w'), indent=1)


def cmd_rouche(a):
    cont = json.load(open(a.contagem))
    eB = Fr(cont['eps_B'])*(1 + Fr(1, 10**9))
    if a.diag:                                       # contagem em outra regiao, pela diagonal de T (diag.npy)
        zs = -np.load(a.diag)
        em = sorted([complex(z) for z in zs if dentro(z)], key=lambda z: z.real)
        cont['em_Omega'] = [[z.real, z.imag] for z in em]
        perto = sorted([complex(z) for z in zs if not dentro(z)], key=dist_bordo)[:4]
        print('fora, mais perto do bordo: ' + ', '.join(f'{z.real:.4f}{z.imag:+.4f}i (dist {dist_bordo(z):.3e})' for z in perto))
    tmax = Fr(0); pior = None
    for f in a.discos:
        dd = json.load(open(f))
        assert dd.get('certF'), f'{f}: disco sem certF'
        for x in dd['pontos']:
            if not x['ok']:
                continue
            tp = Fr(x['t_p']); r = Fr(x['r'])
            assert r*tp < 1
            t = tp/(1 - r*tp)
            if t > tmax:
                tmax, pior = t, (f, x['p'])
    marg = tmax*eB
    ok = marg < 1
    print(f'max t(s) no contorno <= {float(tmax):.3f} (disco {pior}); t eps_B <= {float(marg):.3e}  {"OK" if ok else "FALHA"}')
    print(f'zeros de A + s em Omega = {len(cont["em_Omega"])}: '
          + ', '.join(f'{z[0]:.6f}{z[1]:+.2e}i' for z in cont['em_Omega']))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('cmd', choices=['residuo', 'rouche'])
    ap.add_argument('--dir', type=Path, default=Path('build/rouche_L'))
    ap.add_argument('--bloco', type=int, default=1024)
    ap.add_argument('--contagem', type=Path, default=Path('build/rouche_L_rig/contagem.json'))
    ap.add_argument('--discos', type=Path, nargs='*', default=[])
    ap.add_argument('--regiao', default='0,3,-0.25,0.25', help='re0,re1,im0,im1 (faixa B: 0,3,0.25,0.75)')
    ap.add_argument('--diag', type=Path, default=None, help='rouche: diag(T) em .npy (contagem em outra regiao)')
    a = ap.parse_args()
    global RE0, RE1, IM0, IM1
    RE0, RE1, IM0, IM1 = (float(x) for x in a.regiao.split(','))
    dict(residuo=cmd_residuo, rouche=cmd_rouche)[a.cmd](a)


if __name__ == '__main__':
    main()
