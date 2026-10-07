#!/usr/bin/env python3
"""R3 (revisao S4): etapa A em mu, em PONTO FLUTUANTE, com montagem propria do bloco bordejado.

Roda numa maquina remota, em ~/revisao_s4 (copia de scripts/ e dos dados). Reusa so os blocos de
signed_operator_L (J_livre, bloco_B) e os pesos; o resto (Q0, M^, residuo, S = Y Y^*, normas por
ARPACK) e montado aqui, sem nk_L2_Z.

  M^ = [[I + P_Z B Q0 P_Z, h0w], [h0w^* Q0 P_Z, 0]]    (pesos de L, banda completa |d| <= 40)
  ||M^-1||           : svds (ARPACK) no operador M^-1 via LU
  Y_Z ~ ||M^-1 F_Z'(x0)||, F_Z' = ((J + lam0) h0 + B h0 nas linhas de Z, ||h0w||^2 - 1)
  a   ~ ||M^-1 [P_Z B Q0 P_T; l Q0 P_T]||, colunas de T nos modos |m| < MZ + 16 (banda completa)
Compara com Z.json do autor (tV 52,14972; a_loewner 0,404; YZ 3,0e-6 com as margens)."""
import json, math, sys, time
from pathlib import Path
import numpy as np
from scipy.linalg import lu_factor, lu_solve
from scipy.sparse.linalg import LinearOperator, svds, eigsh

sys.path.insert(0, 'scripts')
import signed_operator_L as sl
import closed_form as cf

T0 = time.time()
def log(s):
    print(f'{s}   [{time.time()-T0:.0f}s]', flush=True)

K1, K2, MU = sl.K1, sl.K2, sl.MU
lam0 = MU                          # mu_RefA (o lambda0 alegado)
MZ, NZ, Nc, BANDA_A = 32, 128, 400, 16
g = sl.campos()

def pesos(r):
    return K1**r.m*np.sqrt(np.concatenate([cf.omega(ns, K2) for _, ns in r.comps]))

def Qw(m, N):                      # (J + lam0)^-1 pesado (o fator K1^m cancela no bloco diagonal)
    r = sl.Radial(m, N)
    w = np.sqrt(np.concatenate([cf.omega(ns, K2) for _, ns in r.comps]))
    return (w[:, None]*np.linalg.inv(sl.J_livre(r, lam0)))/w[None, :]

Bc = {}
def Bw(mo, mi, No, Ni):
    k = (mo - mi, mi % 2, No, Ni)
    if k not in Bc:
        ri = sl.Radial(mi % 2, Ni); ro = sl.Radial(mi % 2 + mo - mi, No)
        Bc[k] = (pesos(ro)[:, None]*sl.bloco_B(g, ro, ri))/pesos(ri)[None, :]
    return Bc[k]

ms = sl.modos(MZ); rs = [sl.Radial(m, NZ) for m in ms]
offs = np.cumsum([0] + [r.dim for r in rs]); n = int(offs[-1])
h0w = np.load('h0w.npy').astype(complex); h0w /= np.linalg.norm(h0w)
assert len(h0w) == n
QZ = {m: Qw(m, NZ) for m in ms}
M = np.zeros((n + 1, n + 1), complex)
M[np.diag_indices(n)] = 1.0
for j, ri in enumerate(rs):
    for i, ro in enumerate(rs):
        if abs(ro.m - ri.m) <= sl.BAND_M:
            M[offs[i]:offs[i+1], offs[j]:offs[j+1]] += Bw(ro.m, ri.m, NZ, NZ) @ QZ[ri.m]
M[:n, n] = h0w
for j, r in enumerate(rs):
    M[n, offs[j]:offs[j+1]] = h0w[offs[j]:offs[j+1]].conj() @ QZ[r.m]
log(f'M^ montada (n+1 = {n+1}), ||M||_F = {np.linalg.norm(M):.4e}')
# residuo F_Z'(x0): (J + lam0) h0 + B h0 nas linhas de Z; ultima linha ||h0w||^2 - 1
rZ = np.zeros(n + 1, complex)
for i, ro in enumerate(rs):
    wo = np.sqrt(np.concatenate([cf.omega(ns, K2) for _, ns in ro.comps]))
    acc = ((wo[:, None]*sl.J_livre(ro, lam0))/wo[None, :]) @ h0w[offs[i]:offs[i+1]]
    for j, ri in enumerate(rs):
        if abs(ro.m - ri.m) <= sl.BAND_M:
            acc += Bw(ro.m, ri.m, NZ, NZ) @ h0w[offs[j]:offs[j+1]]
    rZ[offs[i]:offs[i+1]] = acc
rZ[n] = np.vdot(h0w, h0w) - 1
Bc.clear()
log(f'||F_Z\'(x0)|| = {np.linalg.norm(rZ):.4e}')
lu = lu_factor(M, overwrite_a=True, check_finite=False); del M
log('LU pronta')
xs = lu_solve(lu, rZ)
log(f'Y_Z (sem margens) ~ ||M^-1 F_Z\'|| = {np.linalg.norm(xs):.4e}')
op = LinearOperator((n + 1, n + 1), matvec=lambda v: lu_solve(lu, v),
                    rmatvec=lambda v: lu_solve(lu, v, trans=2), dtype=complex)
s = svds(op, k=3, return_singular_vectors=False, tol=1e-10, maxiter=3000)
tV_est = float(max(s))
log(f'||M^-1|| (ARPACK, 3 maiores valores singulares) = {sorted(s, reverse=True)}')
# S = Y Y^*, Y = [P_Z B Q0 P_T; l Q0 P_T], colunas de T nos modos |mi| < MZ + 16, banda completa nas linhas
S = np.zeros((n + 1, n + 1), complex)
idx = {m: i for i, m in enumerate(ms)}
lote = []; ncol = 0
def descarrega():
    global lote
    if lote:
        Yb = np.concatenate(lote, axis=1); lote = []
        S[...] += Yb @ Yb.conj().T
for mi in range(-(MZ + BANDA_A - 1), MZ + BANDA_A):
    r = sl.Radial(mi, Nc)
    nn = np.concatenate([ns for _, ns in r.comps])
    c = np.where(nn >= NZ)[0] if abs(mi) < MZ else np.arange(r.dim)
    RT = Qw(mi, Nc)[:, c]
    Y = np.zeros((n + 1, len(c)), complex)
    for mo in ms:
        if abs(mo - mi) <= sl.BAND_M:
            i = idx[mo]
            Y[offs[i]:offs[i+1]] = Bw(mo, mi, NZ, Nc) @ RT
    if mi in idx:
        i = idx[mi]
        Y[n] = h0w[offs[i]:offs[i+1]].conj() @ RT[np.where(nn < NZ)[0]]
    lote.append(Y); ncol += len(c)
    if sum(x.shape[1] for x in lote) >= 5000:
        descarrega()
    Bc.clear()
descarrega()
log(f'S = Y Y^* acumulado ({ncol} colunas de T)')
opS = LinearOperator((n + 1, n + 1), matvec=lambda v: lu_solve(lu, S @ lu_solve(lu, v, trans=2)), dtype=complex)
ev = eigsh(opS, k=2, which='LA', return_eigenvectors=False, tol=1e-10, maxiter=3000)
a_est = math.sqrt(float(max(ev.real)))
log(f'a ~ ||M^-1 Y|| (ARPACK) = {a_est:.6f}')
z = json.loads(Path('Z.json').read_text())
out = dict(tV_est=tV_est, sv3=[float(x) for x in s], a_est=a_est, YZ_est=float(np.linalg.norm(xs)),
           rZ=float(np.linalg.norm(rZ)), ncolT=ncol, autor=dict(tV=z['tV'], a_loewner=z['a_loewner'], aZ=z['aZ'], YZ=z['YZ']))
Path('r3_resultado.json').write_text(json.dumps(out, indent=1) + '\n')
print(json.dumps(out, indent=1))
print(f"VEREDITO: tV_est {tV_est:.5f} <= {z['tV']:.5f}: {tV_est <= z['tV']};  a_est {a_est:.5f} <= a_loewner {z['a_loewner']}: {a_est <= z['a_loewner']};"
      f"  YZ_est {np.linalg.norm(xs):.3e} <= YZ {z['YZ']:.3e}: {np.linalg.norm(xs) <= z['YZ']}")
