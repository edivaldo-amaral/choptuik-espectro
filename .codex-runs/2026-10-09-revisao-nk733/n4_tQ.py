#!/usr/bin/env python3
"""N4 (nível Z, tQ): refaz Mh (h0w GRAVADO) e testa a cota de Loewner de tQ >= ||Mh^-1 [Q_ZZ; 0]||
(t^2 Mh Mh^* - X X^* >= 0, X = [Q_ZZ; 0]) no valor GRAVADO, com controle negativo em 0,999 * estimativa.
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

from verified_norms import gamma as gam
lu = lu_factor(M, check_finite=False)
def aplX(u):
    out = np.zeros(n + 1, complex)
    for j, m in enumerate(ms):
        out[offs[j]:offs[j+1]] = QZ[m] @ u[offs[j]:offs[j+1]]
    return out
def aplXh(u):
    out = np.zeros(n, complex)
    for j, m in enumerate(ms):
        out[offs[j]:offs[j+1]] = QZ[m].conj().T @ u[offs[j]:offs[j+1]]
    return out
op = LinearOperator((n + 1, n), matvec=lambda v: lu_solve(lu, aplX(v)), rmatvec=lambda v: aplXh(lu_solve(lu, v, trans=2)), dtype=complex)
est = float(max(svds(op, k=2, return_singular_vectors=False, tol=1e-12)))
log(f'ARPACK: ||Mh^-1 [Q_ZZ; 0]|| ~ {est!r}')
del lu
P, eP = gram_linhas(M); del M
fP = float(np.linalg.norm(P))
np.save(OUT/'P_tQ.npy', P); del P
Sb = []; eS = 0.0; fS2 = 0.0
for m in ms:
    Sm = QZ[m] @ QZ[m].conj().T
    Sb.append(Sm); fS2 += float(np.linalg.norm(Sm))**2
    eS = max(eS, gam(len(Sm))*float(np.linalg.norm(QZ[m]))**2)
fS = math.sqrt(fS2)
t_grav = z['tQ']/FATOR_PESO - z['tV']*dQZ
out = dict(arpack=est, tQ_gravado=z['tQ'], t_loewner_gravado=t_grav)
for nome, t in (('gravado', t_grav), ('controle_negativo', est*0.999)):
    G = np.load(OUT/'P_tQ.npy'); G *= t*t
    for j in range(len(ms)):
        G[offs[j]:offs[j+1], offs[j]:offs[j+1]] -= Sb[j]
    ok = pd_hermitiana([G], t*t*eP + eS + 2*U*(t*t*fP + fS))
    out[f'loewner_{nome}'] = dict(t=t, verificado=bool(ok))
    log(f'Loewner tQ em t = {t:.6f} ({nome}): {"VERIFICADO" if ok else "nao verificado"}')
(OUT/'P_tQ.npy').unlink()
(OUT/'n4_tQ.json').write_text(json.dumps(out, indent=1) + '\n')
log('FIM n4_tQ')
