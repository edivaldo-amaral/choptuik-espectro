#!/usr/bin/env python3
"""N4 (nível Z, a): estimativa INDEPENDENTE (ARPACK, ponto flutuante) de
  a = ||Mh^-1 [P_Z B Q0 P_T; l Q0 P_T]||, banda completa |d| <= 40, com as colunas de T em |m| < 48 (o que o Loewner
  de Z.json cobre, a_loewner) e depois em |m| < 72 (todas as que alcançam Z; comparar com sqrt(a_loewner^2 + a_extra^2)).
Texto original de n4_Z.py: refaz, com o h0w.npy GRAVADO, a matriz bordejada Mh de nk_L2_Z, o resíduo Y_Z e a cota de
Loewner de ||Mh^-1||, e acrescenta controles do revisor:
  (a) estimativa independente de ||Mh^-1|| por ARPACK (svds no operador Mh^-1 via LU), em vez da potência;
  (b) Loewner na cota GRAVADA t = tV de Z.json (tem de passar);
  (c) controle negativo: Loewner em t = 0,999 * estimativa (tem de FALHAR: mostra que o teste não é vazio);
  (d) Y_Z decomposto (||Mh^-1 r||, resíduo da solução, erro de r) e comparado com o gravado.
Roda no laboratorio2 a partir de uma cópia do projeto (~/revisao_nk733/proj), com o dir de dados copiado."""
import json, math, sys, time
from pathlib import Path
import numpy as np

import os
OUT = Path(os.environ.get('N4_OUT', str(Path.home()/'revisao_nk733')))
PROJ = Path(os.environ.get('N4_PROJ', str(OUT/'proj')))
sys.path.insert(0, str(PROJ/'scripts'))
import closed_form as cf
import fundo_L
import nk_L as nk
import signed_operator_L as sl
from nk_L2_Z import Q_cert, err_mm, FATOR_PESO
from verified_inverse import pd_hermitiana, gram_linhas
from verified_norms import gamma, U
from tile_perron import DELTA_MU
from scipy.linalg import lu_factor, lu_solve
from scipy.sparse.linalg import LinearOperator, svds

t0 = time.time()
def log(s):
    print(f'{s}   [{time.time()-t0:.0f}s]', flush=True)

D = Path(os.environ.get('N4_DIR', str(PROJ/'build'/'nk733copia')))
z = json.loads((D/'Z.json').read_text())
lam0 = z['lam0']; MZ, NZ = z['Z']; Mc, Nc = z['F']; banda = z['banda']
log(f'EPS_FUNDO_L = {fundo_L.EPS_FUNDO_L}; lam0 = {lam0!r}; Z = {MZ}x{NZ}')
ctx = nk.Contexto(lam0, Nc)
g = ctx.g
ms = sl.modos(MZ); rs = [ctx.rad(m, NZ) for m in ms]
offs = np.cumsum([0] + [r.dim for r in rs]); n = int(offs[-1])
fb = fundo_L.dB_L(g, banda); fb40 = fundo_L.dB_L(g, sl.BAND_M); emu = fundo_L.eps_mu_L(complex(lam0))
log(f'n = {n}; dB_L banda {banda} {fb}; eps_mu {emu!r}  (Z.json: {z["dB_L"]}, {z["eps_mu"]!r})')
nBL = nk.NORMA_BL + fb['total']
QZ = {}; dQZ = 0.0; nQZ = 0.0
for m in ms:
    Q, dq, nq = Q_cert(m, NZ, lam0); QZ[m] = Q; dQZ = max(dQZ, dq); nQZ = max(nQZ, nq + dq)
log(f'nQZ {nQZ!r} (Z.json {z["nQZ"]!r}); dQZ {dQZ!r} (Z.json {z["dQZ"]!r})')
h0w = np.load(D/'h0w.npy'); assert len(h0w) == n
M = np.zeros((n + 1, n + 1), complex); M[np.diag_indices(n)] = 1.0; e2 = 0.0
for j, ri in enumerate(rs):
    Qj = QZ[ri.m]
    for i, ro in enumerate(rs):
        if abs(ro.m - ri.m) <= sl.BAND_M:
            Bb = ctx.B(ro.m, ri.m, NZ, NZ)
            M[offs[i]:offs[i+1], offs[j]:offs[j+1]] += Bb @ Qj
            e2 += err_mm(Bb, Qj)**2
ctx._B.cache_clear()
M[:n, n] = h0w
for j, r in enumerate(rs):
    M[n, offs[j]:offs[j+1]] = h0w[offs[j]:offs[j+1]].conj() @ QZ[r.m]
fM = float(np.linalg.norm(M))
dM = (math.sqrt(e2) + U*fM + (fb40['arred'] + fb40['rt'])*nQZ + nBL*dQZ + emu + dQZ*(1 + n*U) + gamma(3*NZ)*nQZ)
log(f'dM {dM!r} (Z.json {z["dM"]!r})')

from scipy.sparse.linalg import eigsh
lu = lu_factor(M, overwrite_a=True, check_finite=False); del M
idx = {m: i for i, m in enumerate(ms)}
S = np.zeros((n + 1, n + 1), complex)
def acumula(mis):
    lote = []
    def desc():
        if lote:
            Yb = np.concatenate(lote, axis=1); lote.clear()
            YbH = Yb.conj().T
            for i0 in range(0, n + 1, 2000):
                S[i0:i0+2000] += Yb[i0:i0+2000] @ YbH
            del Yb, YbH
    for mi in mis:
        c = nk.colsT(ctx, mi, MZ, NZ, Nc)
        Q, dq, nq = Q_cert(mi, Nc, lam0)
        RT = Q[:, c]
        Y = np.zeros((n + 1, len(c)), complex)
        for mo in ms:
            if abs(mo - mi) <= sl.BAND_M:
                i = idx[mo]; Y[offs[i]:offs[i+1]] = ctx.B(mo, mi, NZ, Nc) @ RT
        if mi in idx:
            i = idx[mi]
            Y[n] = h0w[offs[i]:offs[i+1]].conj() @ RT[nk.nZ_local(ctx, mi, NZ, Nc)]
        lote.append(Y)
        if sum(x.shape[1] for x in lote) >= 6000:
            desc()
        ctx.Q.cache_clear(); ctx._B.cache_clear()
    desc()
def amax():
    op = LinearOperator((n + 1, n + 1), matvec=lambda v: lu_solve(lu, S @ lu_solve(lu, v, trans=2)), dtype=complex)
    return math.sqrt(float(eigsh(op, k=1, which='LM', return_eigenvectors=False, tol=1e-10)[0].real))
acumula(range(-(MZ + banda - 1), MZ + banda))
a48 = amax(); log(f'a (colunas |m| < 48, banda 40) ~ {a48!r};  a_loewner de Z.json {z["a_loewner"]!r}')
acumula([m for m in range(-(MZ + sl.BAND_M - 1), MZ + sl.BAND_M) if abs(m) >= MZ + banda])
a72 = amax()
a_extra = z['tV']*z['dB_L']['resto']*json.loads((D/'Q.json').read_text())['nQx_extra']
log(f'a (colunas |m| < 72, banda 40) ~ {a72!r};  sqrt(a_loewner^2 + a_extra^2) = {math.hypot(z["a_loewner"], a_extra)!r}; aZ de Z.json {z["aZ"]!r}')
(OUT/'n4_aZ.json').write_text(json.dumps(dict(a48=a48, a72=a72, a_loewner=z['a_loewner'], a_extra_C3=a_extra, aZ=z['aZ']), indent=1) + '\n')
log('FIM n4_aZ')
