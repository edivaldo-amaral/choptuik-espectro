#!/usr/bin/env python3
"""R5 (revisao S4): recalculo independente das grandezas da etapa E e da cota ||x_g - x0||_v.

Proprio: g_A = (Z_A + 1)RefA (xi d_xi e d_tau proprios, r2_identidade), pesos kappa1^m sqrt(omega2(n))
(m com sinal), nz = ||P_Z g_A||, cauda ||P_T g_A||, h0 contra h0w.npy, ||d_tau RefA||, e = ||g* - g_A||
(com ||(Z_A+1)RefB|| e ||RefB|| de r1), ||L_A(lambda0) h0|| na saida completa (aplica_L proprio de r2),
e a formula de ||y_g - y0|| e dist com v de C3.json. Compara com E.json."""
from __future__ import annotations
from fractions import Fraction as Fr
import json, math, sys, time
from pathlib import Path
import numpy as np

AQUI = Path(__file__).resolve().parent
RAIZ = AQUI.parents[1]
sys.path.insert(0, str(RAIZ/'scripts')); sys.path.insert(0, str(AQUI))
import signed_operator_L as sl, closed_form as cf, fundo_L
from tile_perron import MU_80, DELTA_MU
import r2_identidade as R2

out = []
def P(s=''):
    print(s, flush=True); out.append(s)
T0 = time.time()
K1, K2, MU = sl.K1, sl.K2, sl.MU
rig = RAIZ/'build/s4/rig'
z = json.loads((rig/'Z.json').read_text()); c3 = json.loads((rig/'C3.json').read_text()); E = json.loads((rig/'E.json').read_text())
MZ, NZ = z['Z']; lam0 = z['lam0']

def pesos(m, N):
    r = sl.Radial(m, N)
    return K1**m*np.sqrt(np.concatenate([np.where(ns == 0, 1.0, K2**(2.0*ns) + K2**(-2.0*ns)) for _, ns in r.comps]))

g = sl.campos()
Ng = 101
wA = R2.campos_para_modos(g, 41, Ng)
gA = R2.aplica_Z(wA, Ng, MU, 1)                                  # (Z_A + 1) RefA, sem pesos
# nz e cauda (pesos de L), com a caixa Z = MZ x NZ
nz2 = 0.0; t2 = 0.0; partes = []
for m in sl.modos(41):
    r = sl.Radial(m, Ng); nn = np.concatenate([ns for _, ns in r.comps]); vw = pesos(m, Ng)*gA[m]
    if abs(m) < MZ:
        nz2 += float(np.sum(np.abs(vw[nn < NZ])**2)); t2 += float(np.sum(np.abs(vw[nn >= NZ])**2))
    else:
        t2 += float(np.sum(np.abs(vw)**2))
nz = math.sqrt(nz2); t = math.sqrt(t2)
P(f'nz = ||P_Z g_A|| = {nz:.10f} (E.json {E["nz"]:.10f});  cauda ||g_A - P_Z g_A|| = {t:.6e} (E.json {E["cauda"]:.6e})')
# h0 contra h0w.npy e contra gauge_h0w_32x128.npy
h0w = np.load(rig/'h0w.npy').astype(complex)
hz = []
for m in sl.modos(MZ):
    rN = sl.Radial(m, NZ); w_ = pesos(m, NZ)
    full = np.zeros(rN.dim, complex)
    rG = sl.Radial(m, Ng)
    posG = {(c, int(n)): k for k, (c, n) in enumerate((c, n) for c, ns in rG.comps for n in ns)}
    for k, (c, n) in enumerate((c, n) for c, ns in rN.comps for n in ns):
        if (c, int(n)) in posG:
            full[k] = gA[m][posG[(c, int(n))]]
    hz.append(w_*full)
hz = np.concatenate(hz)
P(f'||h0w.npy|| = {np.linalg.norm(h0w):.15f};  ||h0w - P_Z g_A/nz|| = {np.linalg.norm(h0w - hz/nz):.2e}')
hg = np.load(RAIZ/'build/s4/gauge_h0w_32x128.npy').astype(complex)
P(f'gauge_h0w_32x128.npy (entrada da etapa A) normalizado == rig/h0w.npy: {np.array_equal(hg/np.linalg.norm(hg), h0w)} (dif {np.linalg.norm(hg/np.linalg.norm(hg) - h0w):.1e})')
# ||d_tau RefA|| (pesos), e
dt2 = sum(float(np.sum(np.abs(pesos(m, Ng)*(m/2)*wA[m])**2)) for m in wA)
dtA = math.sqrt(dt2)
dmu_inv = float(abs(1/Fr(722873400, 2**32) - 1/MU_80) + Fr(1, 10**70))
ZB, nB_ = 2.454216e-08*(1 + 1e-6), 2.350632e-10*(1 + 1e-6)      # r1 (proprios)
e = ZB + dmu_inv*(dtA + nB_*449/2) + 1e-60
P(f'||d_tau RefA||_L = {dtA:.6f}; |1/mu* - 1/mu_A| <= {dmu_inv:.4e}; e = ||g* - g_A|| <= {e:.6e} (E.json {E["e"]:.6e}; o autor usa 1000 ||RefB|| para ||d_tau RefB||, eu 449/2)')
gn_h0 = (t + 2*e)/(nz - e); gn = (nz + t + e)/(nz - e)
P(f'||g_n - h0|| <= {gn_h0:.6e} (E.json {E["gn_h0"]:.6e});  ||g_n|| <= {gn:.8f}')
# ||L_A(lambda0) h0|| (pesos), saida completa
hu = {}; off = 0
for m in sl.modos(MZ):
    rN = sl.Radial(m, NZ); hu[m] = h0w[off:off + rN.dim]/pesos(m, NZ); off += rN.dim
Mo, No = MZ + 40, NZ + 102
res = R2.aplica_L(g, hu, NZ, No, Mo, lam0, MU)
rA = math.sqrt(sum(float(np.sum(np.abs(pesos(m, No)*x)**2)) for m, x in res.items()))
P(f'||L_A(lambda0) h0|| = {rA:.6e} (E.json rA {E["rA"]:.6e})   [{time.time()-T0:.0f}s]')
# formula de delta e dist (margens do autor para arredondamento; as constantes de erro sao as de E.json)
import nk_L as nk
fb = fundo_L.dB_L(g, sl.BAND_M)
rL = rA + E['arred'] + fb['rt'] + fb['arred'] + DELTA_MU*(z['J1h'] + 2*cf.norma_divxi(K2))
nBst = nk.NORMA_BL + fb['rt'] + DELTA_MU*2*cf.norma_divxi(K2)     # inclui delta_mu S Gamma1 (o autor nao)
J1Q0 = fundo_L.eps_mu_L(complex(lam0))/DELTA_MU
delta = (nBst*gn_h0 + rL + DELTA_MU*gn)/(1 - DELTA_MU*J1Q0)
v = [float(Fr(x)) for x in c3['v']]
dist = math.sqrt(delta**2 + DELTA_MU**2)/min(v)
th, Z2, r = Fr(c3['theta']), Fr(c3['Z2']), Fr(c3['r'])
lip = float(th + 2*Z2*max(Fr(dist), r))
P(f'||y_g - y0|| <= {delta:.6e} (E.json {E["delta"]:.6e}); dist = {dist:.6e} (E.json {E["dist"]:.6e}); r = {float(r):.6e}; theta + 2 Z2 max(dist, r) = {lip:.6f}')
# sensibilidade: quanto dist pode crescer ate perder a unicidade
P(f'folga: unicidade vale ate rho = (1 - theta)/(2 Z2) = {float((1 - th)/(2*Z2)):.4e}, isto e {float((1 - th)/(2*Z2))/dist:.1f} vezes dist')
(AQUI/'r5_saida.txt').write_text('\n'.join(out) + '\n')
