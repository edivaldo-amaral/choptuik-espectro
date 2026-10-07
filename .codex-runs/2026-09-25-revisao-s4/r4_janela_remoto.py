#!/usr/bin/env python3
"""R4 (revisao S4): normas VERDADEIRAS (ponto flutuante, ARPACK) de blocos de E = I - A DF(x0) em mu,
com montagem propria, para comparar com as cotas de build/s4/rig/T3.json.

Janela J = modos 32..47 (linhas e colunas de T: todo n < Nc = 400, pois |m| >= MZ = 32):
  H_J = I + P_J B16 Q0 P_J (a matriz do certificado: banda |d| <= 16, pesos de L), V_J = H_J^-1;
  ||V_J||                           (T3: nV)
  ||V_J P_J B Q0 P_Z||              (T3: TZ; aqui com banda 16 E com a banda completa 40)
  ||V_J P_J B16 Q0 P_{48..63}||     (T3: N[32..47 <- 48..63])
  ||V_J P_J B16 Q0 P_{16..31}||     (T3: N[32..47 <- 16..31]; colunas de T de 16..31 sao n >= 128)
Reusa so J_livre/bloco_B de signed_operator_L e os pesos."""
import json, math, sys, time
from pathlib import Path
import numpy as np
from scipy.linalg import lu_factor, lu_solve
from scipy.sparse.linalg import LinearOperator, svds

sys.path.insert(0, 'scripts')
import signed_operator_L as sl
import closed_form as cf

T0 = time.time()
def log(s):
    print(f'{s}   [{time.time()-T0:.0f}s]', flush=True)

K1, K2, MU = sl.K1, sl.K2, sl.MU
lam0 = MU
MZ, NZ, Nc = 32, 128, 400
g = sl.campos()

def pesos(r):
    return K1**r.m*np.sqrt(np.concatenate([cf.omega(ns, K2) for _, ns in r.comps]))

Qc = {}
def Qw(m):
    if m not in Qc:
        r = sl.Radial(m, Nc)
        w = np.sqrt(np.concatenate([cf.omega(ns, K2) for _, ns in r.comps]))
        Qc[m] = (w[:, None]*np.linalg.inv(sl.J_livre(r, lam0)))/w[None, :]
    return Qc[m]

Bc = {}
def Bw(mo, mi):
    k = (mo - mi, mi % 2)
    if k not in Bc:
        ri = sl.Radial(mi % 2, Nc); ro = sl.Radial(mi % 2 + mo - mi, Nc)
        Bc[k] = (pesos(ro)[:, None]*sl.bloco_B(g, ro, ri))/pesos(ri)[None, :]
    return Bc[k]

def cols(m, parte):
    r = sl.Radial(m, Nc); nn = np.concatenate([ns for _, ns in r.comps])
    if parte == 'Z':
        return np.where(nn < NZ)[0]
    return np.where(nn >= NZ)[0] if abs(m) < MZ else np.arange(r.dim)

def produto(lin, col, banda):
    """P_lin B Q0 P_col (linhas: todo o layout de T do modo)"""
    ml = list(lin); mc = [m for m in col if len(col[m])]
    ol = np.cumsum([0] + [len(lin[m]) for m in ml]); oc = np.cumsum([0] + [len(col[m]) for m in mc])
    X = np.zeros((int(ol[-1]), int(oc[-1])), complex)
    for q, mi in enumerate(mc):
        Qi = Qw(mi)[:, col[mi]]
        for p, mo in enumerate(ml):
            if abs(mo - mi) <= banda:
                X[ol[p]:ol[p+1], oc[q]:oc[q+1]] = Bw(mo, mi)[lin[mo]] @ Qi
    return X

J = list(range(32, 48))
linJ = {m: cols(m, 'T') for m in J}
H = produto(linJ, linJ, 16)
H[np.diag_indices(len(H))] += 1.0
nJ = len(H)
log(f'H_J montada ({nJ})')
lu = lu_factor(H, overwrite_a=True, check_finite=False); del H

def norma_VX(X):
    if X is None:
        op = LinearOperator((nJ, nJ), matvec=lambda v: lu_solve(lu, v), rmatvec=lambda v: lu_solve(lu, v, trans=2), dtype=complex)
    else:
        op = LinearOperator((nJ, X.shape[1]), matvec=lambda v: lu_solve(lu, X @ v),
                            rmatvec=lambda v: X.conj().T @ lu_solve(lu, v, trans=2), dtype=complex)
    return float(max(svds(op, k=2, return_singular_vectors=False, tol=1e-9, maxiter=3000)))

res = {}
res['nV'] = norma_VX(None); log(f"||V_J|| = {res['nV']:.5f}")
colZ = {m: cols(m, 'Z') for m in sl.modos(MZ)}
X = produto(linJ, colZ, 16); res['TZ_b16'] = norma_VX(X); del X; log(f"||V_J P_J B16 Q0 P_Z|| = {res['TZ_b16']:.5f}")
X = produto(linJ, colZ, 40); res['TZ_b40'] = norma_VX(X); del X; log(f"||V_J P_J B40 Q0 P_Z|| = {res['TZ_b40']:.5f}")
Bc.clear()
X = produto(linJ, {m: cols(m, 'T') for m in range(48, 64)}, 16); res['N_48_63'] = norma_VX(X); del X
log(f"||V_J A(32..47 <- 48..63)|| = {res['N_48_63']:.5f}")
X = produto(linJ, {m: cols(m, 'T') for m in range(16, 32)}, 16); res['N_16_31'] = norma_VX(X); del X
log(f"||V_J A(32..47 <- 16..31)|| = {res['N_16_31']:.5f}")
t = json.loads(Path('T3.json').read_text()); i = t['nome'].index('32..47')
aut = dict(nV=t['nV'][i], TZ=t['TZ'][i], N_48_63=t['N'][i][t['nome'].index('48..63')], N_16_31=t['N'][i][t['nome'].index('16..31')])
out = dict(estimativa=res, cota_T3=aut)
Path('r4_janela_resultado.json').write_text(json.dumps(out, indent=1) + '\n')
print(json.dumps(out, indent=1))
for k_est, k_aut in (('nV', 'nV'), ('TZ_b40', 'TZ'), ('N_48_63', 'N_48_63'), ('N_16_31', 'N_16_31')):
    print(f'{k_est}: verdadeira {res[k_est]:.5f} <= cota {aut[k_aut]:.5f}: {res[k_est] <= aut[k_aut]}')
