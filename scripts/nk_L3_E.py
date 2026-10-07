#!/usr/bin/env python3
"""S4, ETAPA E: IDENTIFICACAO do autopar do NK de L em lambda0 = mu_RefA com o autopar de gauge
exato. Depois de A (Z.json, centro h0 = P_Z g_A normalizado), B (T3.json) e C (C3.json: r, v, w).

O fundo verdadeiro w* = RefA + RefB + Corr (bola de RT) tem o autovetor exato g* = (Z* + 1) w* de
L_*(mu*), Z* = mu*^-1 d_tau + xi d_xi (GAUGE_RT_ACTION.md). Normalizado, g_n = g*/<h0, g*>, e
x_g = (y_g, mu*) com y_g = Q0^-1 g_n e zero de F (g* no dominio forte: RADIUS_RECOVERY.md). O NK
da o zero x* em B(x0, r) e, como F e quadratica, ||DT(x)|| <= theta + 2 Z2 ||x - x0||_v para todo x:
se theta + 2 Z2 max(||x_g - x0||_v, r) < 1, T e contracao nessa bola e x_g = x*. Entao mu* e o
lambda* certificado, raiz ALGEBRICAMENTE SIMPLES de L com autoespaco span(g*).

Pela equacao de autovalor (sem derivadas do erro de RT):
    y_g - y0 = -B_*(g_n - h0) - L_*(lambda0) h0 + (lambda0 - mu*) g_n - dmu J1 (g_n - h0),
    ||y_g - y0|| <= (||B_*|| ||g_n - h0|| + ||L_*(lambda0) h0|| + |lambda0 - mu*| ||g_n||)
                    / (1 - |dmu| ||J1 Q0||).
Com g_A = (Z_A + 1) RefA = nz h0 + t (t = parte fora de Z) e e = g* - g_A:
    ||g_n - h0|| <= (||t|| + 2||e||)/(nz - ||e||),
    ||e|| <= ||(Z_A + 1) RefB|| + |mu*^-1 - mu_A^-1| (||d_tau RefA|| + ||d_tau RefB||) + ||(Z* + 1) Corr||.
||(Z* + 1) Corr||_L <= 1,1e-75 (corr_gauge_bound.py: equacao da correcao de RT, O_mu* Corr = G com
||G|| <= 5,5e-80, e cotas de banda de RT); aqui entra como 1e-60. Norma do certificado: ||x||_v = max(||x_Z'||/v0,
||x_T||_W/vT, ||x_f||/vf), com w <= 1: ||x_g - x0||_v <= sqrt(delta^2 + dlambda^2)/min(v)."""
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
import gauge_vector_L as G
import nk_L as nk
import signed_operator_L as sl
from tile_perron import DELTA_MU, MU_80
from verified_norms import gamma


def residuo_h0(h0w, MZ, NZ, lam0):
    """||L_A(lam0) h0|| (pesos de L), saida completa: h0 em Z (pesado), linhas em |m| < MZ + 40,
    n < NZ + 101. Devolve o valor e a cota do arredondamento: por componente, |fl(acc) - acc| <=
    gamma(K + 1) (|B| |hu| + |J| |hu|), K = comprimento dos produtos mais as somas de blocos (o erro
    relativo u de hu = h0w/pesos e das entradas de J entra no +1; o de B esta em dB_L['arred'])."""
    g = sl.campos()
    ms = sl.modos(MZ); rs = [sl.Radial(m, NZ) for m in ms]
    off = np.cumsum([0] + [r.dim for r in rs])
    hu = {m: h0w[off[i]:off[i+1]]/nk.pesos(rs[i]) for i, m in enumerate(ms)}     # sem pesos
    No = NZ + sl.BAND_N + 1
    K = int(off[-1]) + len(ms) + 2
    tot = 0.0; tabs = 0.0
    for mo in range(-(MZ + sl.BAND_M - 1), MZ + sl.BAND_M):
        ro = sl.Radial(mo, No); acc = np.zeros(ro.dim, complex); ab = np.zeros(ro.dim)
        for i, mi in enumerate(ms):
            if abs(mo - mi) <= sl.BAND_M:
                Bb = sl.bloco_B(g, ro, rs[i])
                acc += Bb @ hu[mi]; ab += np.abs(Bb) @ np.abs(hu[mi])
        if mo in hu:
            i = ms.index(mo)
            pos = G_pos(ro)
            idx = np.array([pos[(c, int(n))] for c, ns in rs[i].comps for n in ns])
            Jb = sl.J_livre(rs[i], lam0)
            acc[idx] += Jb @ hu[mo]; ab[idx] += np.abs(Jb) @ np.abs(hu[mo])
        w = nk.pesos(ro)
        tot += float(np.sum(np.abs(w*acc)**2)); tabs += float(np.sum((w*ab)**2))
    gK = gamma(K + 1)
    r = math.sqrt(tot)
    return r, gK*math.sqrt(tabs)*(1 + gK)*1.01 + gK*r


def G_pos(r):
    pos = {}; o = 0
    for c, ns in r.comps:
        for k, n in enumerate(ns):
            pos[(c, int(n))] = o + k
        o += len(ns)
    return pos


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--dir', type=Path, default=Path('build/s4/rig'))
    ap.add_argument('--h0', type=Path, default=None, help='centro (padrao: dir/h0w.npy, o da etapa A)')
    ap.add_argument('--ZB', type=float, default=2.4542e-08, help='||(Z_A + 1) RefB||_L (gauge_refB.py)')
    ap.add_argument('--B', type=float, default=2.3506e-10, help='||RefB||_L (gauge_refB.py)')
    a = ap.parse_args()
    z = json.loads((a.dir/'Z.json').read_text())
    c3 = json.loads((a.dir/'C3.json').read_text()) if (a.dir/'C3.json').exists() else None
    MZ, NZ = z['Z']; lam0 = z['lam0']
    muA = sl.MU; assert lam0 == muA
    # g_A = (Z_A + 1) RefA nos pesos de L: parte em Z (nz) e cauda (t)
    cauda_rel, nz = G.salvar_Z(MZ, NZ, a.dir/'_g_tmp.npy'); (a.dir/'_g_tmp.npy').unlink()
    # + 1e-12 nz: o centro h0 (h0w.npy) e P_Z g_A/nz so a menos de arredondamento (relativo ~1e-15)
    t = cauda_rel*nz*(1 + 1e-9) + 1e-12*nz
    # ||d_tau RefA|| (pesos de L): i m/2 nos coeficientes
    g = sl.campos(); dt2 = 0.0
    for m in sl.modos(50):
        r = sl.Radial(m, 110); v = np.zeros(r.dim, complex)
        for c, ns in r.comps:
            aa = sl.modo(g, c, m); o0, o1 = r.off[c]
            full = np.zeros(110, complex); full[:len(aa)] = aa
            v[o0:o1] = full[ns]*(m/2)
        dt2 += float(np.sum(np.abs(nk.pesos(r)*v)**2))
    dtA = math.sqrt(dt2)*1.0001
    # |mu*^-1 - mu_A^-1|, com |mu* - MU_80| <= 1e-80 (RT)
    dmu_inv = float(abs(1/Fr(722873400, 2**32) - 1/MU_80) + Fr(1, 10**70))*1.0001
    e = a.ZB*1.0001 + dmu_inv*(dtA + a.B*1e3) + 1e-60
    gn_h0 = (t + 2*e)/(nz - e)
    gn = (nz + t + e)/(nz - e)
    # ||L_*(lam0) h0||: parte RefA calculada; RT (RefB + Corr em B) e mu
    h0w = np.load(a.h0 or a.dir/'h0w.npy').astype(complex)             # o centro EXATO da etapa A
    nh = float(np.linalg.norm(h0w)); assert abs(nh - 1) < 1e-12
    rA, er = residuo_h0(h0w, MZ, NZ, lam0)
    fb = fundo_L.dB_L(g, sl.BAND_M)
    rL = rA*(1 + 1e-9) + er + fb['rt'] + fb['arred'] + DELTA_MU*(z['J1h'] + 2*cf.norma_divxi(sl.K2))
    nB = nk.NORMA_BL + fb['rt'] + DELTA_MU*2*cf.norma_divxi(sl.K2)       # mu* S Gamma1 - mu_A S Gamma1
    dlam = DELTA_MU                                                  # |lambda0 - mu*| (mu_A diadico)
    J1Q0 = fundo_L.eps_mu_L(complex(lam0))/DELTA_MU                  # ||(J1 + S Gamma1) Q0|| majora ||J1 Q0||
    delta = (nB*gn_h0 + rL + dlam*gn)/(1 - DELTA_MU*J1Q0)*(1 + 1e-9)
    print(f'g_A: ||P_Z g|| = {nz:.5f}, cauda {t:.2e}; ||g* - g_A|| <= {e:.2e}; ||g_n - h0|| <= {gn_h0:.2e}')
    print(f'||L_A(lambda0) h0|| = {rA:.3e} (arredondamento {er:.1e}); com RT e mu: {rL:.3e}; ||B_*|| <= {nB:.4f}')
    print(f'||y_g - y0|| <= {delta:.3e}')
    if c3 is None:
        print('sem C3.json: so a previa (etapa C ainda nao rodou)')
        return
    # ||x||_v = max(||x_Z'||/v0, ||W x_T||/vT, ||x_f||/vf) <= max(1, max w) ||x||/min v
    v = [float(Fr(x)) for x in c3['v']]; r = float(Fr(c3['r'])); wmax = max(1.0, max(float(Fr(x)) for x in c3['w']))
    dist = math.sqrt(delta**2 + dlam**2)*wmax/min(v)*(1 + 1e-12)
    # unicidade: o zero x* (em B(x0, r)) e o unico em B(x0, rho) sempre que theta + 2 Z2 rho < 1 (dois
    # zeros la dariam ||T x1 - T x2|| <= (theta + 2 Z2 rho) ||x1 - x2||); rho = max(dist, r)
    th, Z2 = Fr(c3['theta']), Fr(c3['Z2'])
    rho = max(Fr(dist), Fr(c3['r']))
    lip = th + 2*Z2*rho
    rho_max = (1 - th)/(2*Z2)
    print(f'||x_g - x0||_v <= {dist:.3e};  raio NK r = {r:.3e}, unicidade ate (1 - theta)/(2 Z2) = {float(rho_max):.3e}'
          f'  (theta + 2 Z2 rho = {float(lip):.6f})')
    ok = lip < 1
    print('IDENTIFICADO: o autopar do NK e o de gauge; mu* e raiz ALGEBRICAMENTE SIMPLES de L, autoespaco span(g*)'
          if ok else 'NAO identificado')
    (a.dir/'E.json').write_text(json.dumps(dict(nz=nz, cauda=t, e=e, gn_h0=gn_h0, rA=rA, arred=er, rL=rL,
                                                rt=fb['rt'], delta=delta, dist=dist, r=r, rho_max=float(rho_max),
                                                lip=str(lip), identificado=bool(ok)), indent=1) + '\n')


if __name__ == '__main__':
    main()
