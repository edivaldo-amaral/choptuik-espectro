#!/usr/bin/env python3
"""N4 (janelas): controles do revisor para janelas escolhidas de T3.json.

Para cada janela J: monta Hh_J = I + fl(P_J B Q0 P_J) com a MESMA montagem do certificado (Cauda3.produto_cert),
e estima em ponto flutuante, por ARPACK (svds com a LU de Hh_J; método diferente da potência e do Loewner):
  ||V_J||, ||V_J P_alto||, ||V_J X_JI|| (vizinhas a <= 16 modos), ||V_J X_JZ|| e ||V_J x_res|| (resolução direta).
As cotas rigorosas (T3.json gravado e a reexecução de nk_L3_T.py) têm de ficar ACIMA destas estimativas e perto
delas. Controle negativo: Loewner de ||V_J|| em t = 0,999 * estimativa tem de FALHAR."""
import argparse, json, math, sys, time
from pathlib import Path
import numpy as np

import os
OUT = Path(os.environ.get('N4_OUT', str(Path.home()/'revisao_nk733')))
PROJ = Path(os.environ.get('N4_PROJ', str(OUT/'proj')))
sys.path.insert(0, str(PROJ/'scripts'))
import nk_L as nk
from nk_L3_T import Cauda3
from verified_inverse import gram_linhas, pd_hermitiana
from verified_norms import U
from scipy.linalg import lu_factor, lu_solve
from scipy.sparse.linalg import LinearOperator, svds

ap = argparse.ArgumentParser()
ap.add_argument('janelas', nargs='+')
a0 = ap.parse_args()
D = Path(os.environ.get('N4_DIR', str(PROJ/'build'/'nk733copia')))
z = json.loads((D/'Z.json').read_text()); T = json.loads((D/'T3.json').read_text())
a = argparse.Namespace(lam0=z['lam0'], Z='32x128', F='200x400', W=16, folga_far=60, dir=D, rouche=False)
C = Cauda3(a)
ctx, Nc = C.ctx, C.Nc
niv = C.niveis
h0w = np.load(D/'h0w.npy')
Yr = C.residuo_janelas(z, h0w)
colZ = {m: nk.nZ_local(ctx, m, C.NZ, Nc) for m in nk.sl.modos(C.MZ)}
lim = Nc - nk.BAND_N - 60
saida = {}
t0 = time.time()

def sv(lu, X, k=1):
    """maior valor singular de H^-1 X (X denso ou None = identidade), ARPACK."""
    if X is None:
        mv = lambda v: lu_solve(lu, v); rmv = lambda v: lu_solve(lu, v, trans=2); shp = (lu[0].shape[0],)*2
    else:
        mv = lambda v: lu_solve(lu, X @ v); rmv = lambda v: X.conj().T @ lu_solve(lu, v, trans=2)
        shp = (lu[0].shape[0], X.shape[1])
    op = LinearOperator(shp, matvec=mv, rmatvec=rmv, dtype=complex)
    return float(max(svds(op, k=k, return_singular_vectors=False, tol=1e-10)))

for nome in a0.janelas:
    i = [lv['nome'] for lv in niv].index(nome); r = niv[i]
    X, erro, nQ = C.produto_cert(C.lin_nivel(r, Nc), r['sel'], Nc)
    H = X; H[np.diag_indices(len(H))] += 1.0; del X
    lu = lu_factor(H, check_finite=False)
    res = dict(dim=r['dim'])
    res['nV'] = sv(lu, None)
    nn = np.concatenate([np.concatenate([ns for _, ns in ctx.rad(m, Nc).comps])[r['sel'][m]]
                         for m in r['modos'] if len(r['sel'][m])])
    alto = np.where(nn >= lim)[0]
    E = np.zeros((r['dim'], len(alto))); E[alto, np.arange(len(alto))] = 1.0
    res['nVhi'] = sv(lu, E) if len(alto) else 0.0; del E
    res['N'] = {}
    for j, c in enumerate(niv):
        if j != i and C.perto(r['modos'], c['modos']):
            Xj, _, _ = C.produto_cert(C.lin_nivel(r, Nc), c['sel'], Nc)
            res['N'][c['nome']] = sv(lu, Xj); del Xj
    if C.perto(r['modos'], list(colZ)):
        Xz, _, _ = C.produto_cert(C.lin_nivel(r, Nc), colZ, Nc)
        res['TZ'] = sv(lu, Xz); del Xz
    if Yr[i] is not None:
        x, e = Yr[i]
        res['Y_sem_erro'] = float(np.linalg.norm(lu_solve(lu, x)))
        res['Y_erro_x'] = e
    del lu
    # controle negativo do Loewner (||V_J||): t abaixo da estimativa tem de falhar
    P, eP = gram_linhas(H); del H
    fP = float(np.linalg.norm(P)); n = P.shape[0]
    tneg = res['nV']*0.999
    P *= tneg*tneg; P[np.diag_indices(n)] -= 1.0
    res['loewner_nV_controle_negativo'] = dict(t=tneg, verificado=bool(pd_hermitiana([P], tneg*tneg*eP + 2*U*(tneg*tneg*fP + math.sqrt(n)))))
    del P
    grav = dict(nV=T['nV'][i], nVhi=T['nVhi'][i], TZ=T['TZ'][i], Y=T['Y'][i],
                N={T['nome'][j]: T['N'][i][j] for j in range(26) if j != i and T['nome'][j] in res['N']})
    res['gravado_T3'] = grav
    saida[nome] = res
    print(f'{nome} [{time.time()-t0:.0f}s]: ' + json.dumps(res), flush=True)
    ctx.Q.cache_clear(); ctx._B.cache_clear()
    (OUT/'n4_janelas_float.json').write_text(json.dumps(saida, indent=1) + '\n')
print('FIM n4_janelas_float', flush=True)
