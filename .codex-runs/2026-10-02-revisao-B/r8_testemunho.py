#!/usr/bin/env python3
"""R8 (revisor): o testemunho de constraints da raiz da faixa B.

(a) Caminho INDEPENDENTE do signed_constraints do autor: as matrizes exportadas pelo proprio codigo de RT
    (build/spectrum/rt-A-<malha>.dat e .dat.constraints, base real de RT, C(s) = C0 + s C1, L(s) = L0 + s):
    autopares de L0 (s = -autovalor) perto de 0,0414 + 0,5i (raiz da faixa B), 0,733 (modo fisico), 0,401 (s_K),
    0 (fase) e 0,168 (gauge), e a razao ||C(s) h|| / ||h|| (euclidiana na base de RT). Se o metodo esta certo, a raiz
    da faixa B deve ter razao O(1), como 0,401, e 0,733/0/0,168 devem ter razao ~ 0.
(b) c_ref = ||Pi_F C(lambda0) h0|| recalculado com montagem e pesos proprios (pesos de K (129/128, 9/8) na saida,
    de L (65/64, 5/4) na entrada), blocos nao pesados de signed_constraints.bloco_C, h0w.npy do NK, F = 12x48.
(c) lambda0 do NK contra o autovalor de A (diag de T do laboratorio2) e contra a regiao Omega_B."""
import json, math, sys
from pathlib import Path
import numpy as np
RAIZ = Path('<repo>'); sys.path.insert(0, str(RAIZ/'scripts'))
import independent_numeric_audit as NA
import signed_operator_L as sl, signed_constraints as sc, closed_form as cf
res = {}
alvos = [0.0414194105 + 0.5j, 0.7331796597, 0.4010247330, 0.0, 0.1683070790]
for malha in ('10x30', '12x36'):
    L0, addr = NA.read_rt_matrix(NA.SPEC/f'rt-A-{malha}.dat', 'CHOPTUIK_RT_SPECTRAL_MATRIX_V1')
    C0, C1, caddr = NA.read_constraints(NA.SPEC/f'rt-A-{malha}.dat.constraints', len(L0))
    L0 = np.array(L0, float) if not isinstance(L0, np.ndarray) else L0
    lam, V = np.linalg.eig(L0)
    ss = -lam
    linha = {}
    for a in alvos:
        i = int(np.argmin(np.abs(ss - a)))
        s = ss[i]; h = V[:, i]
        r = float(np.linalg.norm((C0 + s*C1) @ h)/np.linalg.norm(h))
        rl = float(np.linalg.norm(L0 @ h + s*h)/np.linalg.norm(h))
        linha[str(a)] = dict(s=[s.real, s.imag], razao_C=r, residuo_L=rl)
        print(f'{malha}: s = {s.real:+.6f}{s.imag:+.6f}i (alvo {a}): ||C(s)h||/||h|| = {r:.4e}  (residuo de L {rl:.1e})')
    res[malha] = linha
# (b)
MZ, NZ = 32, 128; MF, NF = 12, 48
z = json.load(open(RAIZ/'build/faixaB/nk/Z.json'))
lam0 = complex(*z['lam0'])
h0 = np.load(RAIZ/'build/faixaB/nk/h0w.npy')
print(f'lambda0 = {lam0}, |h0w| = {np.linalg.norm(h0):.12f}, dim {len(h0)}')
g = sl.campos()
K1L, K2L = 65/64, 5/4; K1K, K2K = 129/128, 9/8
def om(n, k2):
    return np.where(n == 0, 1.0, k2**(2.0*n) + k2**(-2.0*n))
ms = sl.modos(MZ); rs = [sl.Radial(m, NZ) for m in ms]; offs = np.cumsum([0] + [r.dim for r in rs])
partes = []
for mo in [m for m in sl.modos(MF) if m % 2 == 0]:
    so = sc.Saida(mo, NF)
    wout = K1K**mo*np.sqrt(om(np.concatenate([ns for _, ns in so.comps]), K2K))
    acc = np.zeros(so.dim, complex)
    for j, mi in enumerate(ms):
        if abs(mo - mi) <= sl.BAND_M:
            win = K1L**mi*np.sqrt(om(np.concatenate([ns for _, ns in rs[j].comps]), K2L))
            X = sc.bloco_C(g, so, rs[j], lam0)
            acc += wout*(X @ (h0[offs[j]:offs[j+1]]/win))
    partes.append(acc)
c = np.concatenate(partes)
cref = float(np.linalg.norm(c))
D = json.load(open(RAIZ/'build/faixaB/nk/D.json'))
print(f'(b) c_ref (revisor, ponto flutuante) = {cref:.6f}; D.json c_ref (cota inferior) = {D["c_ref"]:.6f}; folga D.json = {D["folga"]:.5f}')
res['c_ref_revisor'] = cref; res['c_ref_autor'] = D['c_ref']
# sensibilidade a caixa F (Pi_F e projecao: ||Pi_F x|| cresce com F)
for MF2, NF2 in ((8, 32), (16, 64)):
    pp = []
    for mo in [m for m in sl.modos(MF2) if m % 2 == 0]:
        so = sc.Saida(mo, NF2)
        wout = K1K**mo*np.sqrt(om(np.concatenate([ns for _, ns in so.comps]), K2K))
        acc = np.zeros(so.dim, complex)
        for j, mi in enumerate(ms):
            if abs(mo - mi) <= sl.BAND_M:
                win = K1L**mi*np.sqrt(om(np.concatenate([ns for _, ns in rs[j].comps]), K2L))
                acc += wout*(sc.bloco_C(g, so, rs[j], lam0) @ (h0[offs[j]:offs[j+1]]/win))
        pp.append(acc)
    v = float(np.linalg.norm(np.concatenate(pp)))
    print(f'    F = {MF2}x{NF2}: ||Pi_F C(lambda0) h0|| = {v:.6f}')
    res[f'c_ref_F{MF2}x{NF2}'] = v
# (c)
dT = np.load(RAIZ/'build/rouche_L_rig/diagT_lab2.npy')
sA = -dT[np.argmin(np.abs(-dT - lam0))]
C3 = json.load(open(RAIZ/'build/faixaB/nk/C3.json'))
from fractions import Fraction as Fr
el = float(Fr(C3['erro_lambda']))
print(f'(c) autovalor de A (Schur do lab2) mais proximo: s = {sA}; |s_A - lambda0| = {abs(sA - lam0):.2e}; '
      f'|lambda* - lambda0| <= {el:.3e}; Re lambda* >= {lam0.real - el:.6f} > 0; Im lambda* em [{lam0.imag - el:.6f}, {lam0.imag + el:.6f}]')
res['lambda'] = dict(lam0=[lam0.real, lam0.imag], sA=[sA.real, sA.imag], erro=el)
json.dump(res, open(RAIZ/'.codex-runs/2026-10-02-revisao-B/r8_resultado.json', 'w'), indent=1)
