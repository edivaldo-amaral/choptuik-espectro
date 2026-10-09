#!/usr/bin/env python3
"""N4 (nível Z): refaz, com o h0w.npy GRAVADO, a matriz bordejada Mh de nk_L2_Z, o resíduo Y_Z e a cota de
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
# residuo
rZ = np.zeros(n + 1, complex); er2 = 0.0
for i, ro in enumerate(rs):
    wo = np.sqrt(np.concatenate([cf.omega(ns, sl.K2) for _, ns in ro.comps]))
    Jw = (wo[:, None]*sl.J_livre(ro, lam0))/wo[None, :]
    acc = Jw @ h0w[offs[i]:offs[i+1]]; er2 += err_mm(Jw, h0w[offs[i]:offs[i+1], None])**2
    for j, ri in enumerate(rs):
        if abs(ro.m - ri.m) <= sl.BAND_M:
            Bb = ctx.B(ro.m, ri.m, NZ, NZ)
            acc += Bb @ h0w[offs[j]:offs[j+1]]; er2 += err_mm(Bb, h0w[offs[j]:offs[j+1], None])**2
    rZ[offs[i]:offs[i+1]] = acc
ctx._B.cache_clear()
rZ[n] = np.vdot(h0w, h0w) - 1
J1h = 0.0
for i, r in enumerate(rs):
    wo = np.sqrt(np.concatenate([cf.omega(ns, sl.K2) for _, ns in r.comps]))
    J1 = (sl.J_livre(r, 0.0, mu=1.0) - 1j*r.m/2*np.eye(r.dim))
    J1h += float(np.linalg.norm(((wo[:, None]*J1)/wo[None, :]) @ h0w[offs[i]:offs[i+1]]))**2
termo_mu = DELTA_MU*(math.sqrt(J1h)*1.0001 + 2*cf.norma_divxi(sl.K2))
dr = (math.sqrt(er2) + U*float(np.linalg.norm(rZ)) + (fb40['arred'] + fb40['rt'])*1.0001 + termo_mu + 4*n*U)
log(f'||r_Z|| (calculado) {np.linalg.norm(rZ):.3e}; J1h {math.sqrt(J1h)*1.0001!r} (Z.json {z["J1h"]!r}); '
    f'dr {dr:.4e} = produto {math.sqrt(er2):.2e} + RT/arred {(fb40["arred"] + fb40["rt"])*1.0001:.3e} + mu {termo_mu:.2e}')
lu = lu_factor(M, check_finite=False)
xs = lu_solve(lu, rZ); res = M @ xs - rZ
eres = gamma(n + 1)*fM*float(np.linalg.norm(xs)) + U*float(np.linalg.norm(res))
nres = float(np.linalg.norm(res))*(1 + 1e-6) + eres; nxs = float(np.linalg.norm(xs))*(1 + 1e-6)
log(f'||Mh^-1 r|| ~ {nxs:.3e}; residuo da solucao {nres:.2e}')
# (a) ARPACK: maior valor singular de Mh^-1
op = LinearOperator((n + 1, n + 1), matvec=lambda v: lu_solve(lu, v), rmatvec=lambda v: lu_solve(lu, v, trans=2),
                    dtype=complex)
s = svds(op, k=3, return_singular_vectors=False, tol=1e-12)
est = float(max(s))
log(f'ARPACK: 3 maiores valores singulares de Mh^-1: {sorted(s, reverse=True)}')
del lu
P, eP = gram_linhas(M); del M
fP = float(np.linalg.norm(P))
np.save(OUT/'P_Z.npy', P)
tV = z['tV']
YZ = (nxs + tV*(nres + dr))*FATOR_PESO; rhoZ = tV*dM*FATOR_PESO
log(f'Y_Z refeito {YZ!r} (Z.json {z["YZ"]!r}); rho_Z refeito {rhoZ!r} (Z.json {z["rhoZ"]!r})')
out = dict(dM=dM, dM_gravado=z['dM'], nxs=nxs, nres=nres, dr=dr, YZ=YZ, YZ_gravado=z['YZ'], rhoZ=rhoZ,
           rhoZ_gravado=z['rhoZ'], nQZ=nQZ, dQZ=dQZ, J1h=math.sqrt(J1h)*1.0001, arpack=est, tV_gravado=tV)
def teste(t):
    G = np.load(OUT/'P_Z.npy')
    G *= t*t; G[np.diag_indices(n + 1)] -= 1.0
    return pd_hermitiana([G], t*t*eP + 2*U*(t*t*fP + math.sqrt(n + 1)))
del P
for nome, t in (('gravado', tV/FATOR_PESO), ('controle_negativo', est*0.999)):
    ok = teste(t)
    out[f'loewner_{nome}'] = dict(t=t, verificado=bool(ok))
    log(f'Loewner t^2 Mh Mh^* - I >= 0 em t = {t:.6f} ({nome}): {"VERIFICADO" if ok else "nao verificado"}')
(OUT/'P_Z.npy').unlink()
(OUT/'n4_Z.json').write_text(json.dumps(out, indent=1) + '\n')
log('FIM n4_Z')
