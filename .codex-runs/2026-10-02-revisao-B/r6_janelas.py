#!/usr/bin/env python3
"""R6 (revisor): janelas sem zeros na faixa B.

(1) a particao REGIOES_B cobre [0,3] x [1/4,3/4] (racionais; cada faixa vertical coberta em Im);
(2) a cota do autor em 2 amostras, refeita com codigo proprio para as normas por modo (SVD, ponto flutuante) e a formula:
    (a) regiao 5 (centro 1,0), q = 1,1 + 0,75i, r = 0,1, janela -47..-32 (nao central, d ~ 0,86);
    (b) regiao 1 (centro 0,5i), q = 0,3 + 0,5625i, r = 0,0625, janela -15..0 (central);
(3) estimativa DIRETA de ||F(s)|| = ||I - V_J(c) H^_J(s)|| = ||(I + X(c))^-1 (X(c) - X(s))||, X(s) = P_J B_16 Q0(s) P_J,
    na janela -47..-32, c = 1, s = 1,1 + 0,75i e s = 0,8 + 0,75i (LU densa 11200 + Lanczos). Roda no PC 3.
argv: pasta com T3_1_0.json e T3_0_5.json; caminho dos scripts do projeto."""
import json, math, sys, time
from fractions import Fraction as Fr
from pathlib import Path
import numpy as np
P = Path(sys.argv[1]); sys.path.insert(0, sys.argv[2])
import signed_operator_L as sl, fundo_L
from scipy.linalg import solve_triangular, lu_factor, lu_solve
from scipy.sparse.linalg import LinearOperator, eigsh
t0 = time.time()
def log(*x):
    print(*x, f'[{time.time()-t0:.0f}s]', flush=True)
res = {}
# (1)
R = [(0.0, 0.3, 0.25, 0.375), (0.0, 0.3, 0.375, 0.625), (0.0, 0.3, 0.625, 0.75), (0.3, 0.8, 0.25, 0.5), (0.3, 0.8, 0.5, 0.75),
     (0.8, 1.2, 0.25, 0.75), (1.2, 1.8, 0.25, 0.75), (1.8, 2.2, 0.25, 0.75), (2.2, 2.8, 0.25, 0.75), (2.8, 3.0, 0.25, 0.75)]
Rq = [tuple(Fr(str(x)) for x in r) for r in R]
xs = sorted({r[0] for r in Rq} | {r[1] for r in Rq})
ok = xs[0] == 0 and xs[-1] == 3
for a, b in zip(xs, xs[1:]):
    ys = sorted((r[2], r[3]) for r in Rq if r[0] <= a and r[1] >= b)
    y = Fr(1, 4)
    for lo, hi in ys:
        if lo > y:
            ok = False
        y = max(y, hi)
    ok = ok and y >= Fr(3, 4)
log(f'(1) particao REGIOES_B cobre [0,3] x [1/4,3/4]: {ok}')
res['particao'] = ok
# normas por modo de Q0 (pesos sqrt omega), proprias
MU, K1, K2 = sl.MU, sl.K1, sl.K2; MZ, NZ, Nc = 32, 128, 400
def om(n):
    return np.where(n == 0, 1.0, K2**(2.0*n) + K2**(-2.0*n))
def Qmodo(m, s, N=Nc):
    comps = [np.arange(1, N, 2), np.arange(N), np.arange(N)] if m % 2 == 0 else [np.arange(N)]
    bl = []; nn = []
    for ns in comps:
        J = np.diag(MU*(ns + 1.0) + s + 1j*m/2).astype(complex) + np.triu(np.tile(2*MU*ns, (len(ns), 1)), 1)
        Q = solve_triangular(J, np.eye(len(ns))); w = np.sqrt(om(ns))
        bl.append((w[:, None]*Q)/w[None, :]); nn.append(ns)
    n = sum(len(b) for b in bl); X = np.zeros((n, n), complex); o = 0
    for b in bl:
        X[o:o+len(b), o:o+len(b)] = b; o += len(b)
    return X, np.concatenate(nn)
n2 = lambda X: float(np.linalg.norm(X, 2)) if X.size else 0.0
rt = fundo_L.EPS_FUNDO_L
# (2a) nao central
t3 = json.load(open(P/'T3_1_0.json')); i = t3['nome'].index('-47..-32')
c = 1.0; q = 1.1 + 0.75j; r = 0.1
d = abs(q - c) + r
nQm = max(n2(Qmodo(m, c)[0]) for m in range(-47, -31))
emu_s = fundo_L.eps_mu_L(complex(abs(q.real) + r, abs(q.imag) + r)); emu_c = t3['eps_mu']
emu_t = emu_s + max(0.0, emu_s - emu_c)
fq = 1/(1 - d*t3['nQm'][i]); nQs = t3['nQm'][i]*fq
v = t3['N'][i][i] + d*t3['K'][i]*fq + t3['nV'][i]*(rt*nQs + emu_t*max(1.0, nQs))
log(f'(2a) janela -47..-32, c = 1, q = {q}, r = {r}: d = {d:.4f}; nQm (T3) = {t3["nQm"][i]:.6f}, max ||Q0(c) P_m|| (SVD, revisor) = {nQm:.6f}; '
    f'd nQm = {d*t3["nQm"][i]:.4f}; cota rho_J = {v:.4f} (autor: 0,1517)')
res['2a'] = dict(d=d, nQm_T3=t3['nQm'][i], nQm_revisor=nQm, cota=v)
# (2b) central
t3b = json.load(open(P/'T3_0_5.json')); ib = t3b['nome'].index('-15..0')
cb = 0.5j; qb = 0.3 + 0.5625j; rb = 0.0625
db = abs(qb - cb) + rb
qzj = qjj = 0.0; qzc = 0.0
for m in range(-15, 1):
    Qq, nn = Qmodo(m, qb); Qc_, _ = Qmodo(m, cb)
    iz, iT = np.where(nn < NZ)[0], np.where(nn >= NZ)[0]
    qzz_m, qzj_m, qjj_m = n2(Qq[np.ix_(iz, iz)]), n2(Qq[np.ix_(iz, iT)]), n2(Qq[np.ix_(iT, iT)])
    qjj = max(qjj, qjj_m/(1 - rb*qjj_m)); qzj = max(qzj, qzj_m/((1 - rb*qzz_m)*(1 - rb*qjj_m)))
    qzc = max(qzc, n2(Qc_[np.ix_(iz, iT)]))
emu_s = fundo_L.eps_mu_L(complex(abs(qb.real) + rb, abs(qb.imag) + rb)); emu_t = emu_s + max(0.0, emu_s - t3b['eps_mu'])
sch = json.load(open(P/'schur.json'))['defl_refA']; eE = json.load(open(P/'E.json'))['e']; U = 2.0**-53
errZ = [1.4563e-09 + 1e-60 + U*sch['norma_x'][0]*10, eE + U*sch['norma_x'][1]*10]
dET = sum(dk*nphi*(cx + e) for dk, nphi, cx, e in zip(sch['delta'], sch['norma_phi'], sch['cauda_x'], errZ))
nQs = qzj + qjj
vb = t3b['N'][ib][ib] + db*(t3b['TZ'][ib]*qzj + (t3b['K'][ib] + t3b['TZ'][ib]*qzc)*(1 + db*qjj)) \
    + t3b['nV'][ib]*(rt*nQs + emu_t*max(1.0, nQs) + dET*qzj)
log(f'(2b) janela -15..0, c = 0,5i, q = {qb}, r = {rb}: d = {db:.4f}, q_ZJ(s) <= {qzj:.4f}, q_JJ(s) <= {qjj:.4f}, '
    f'q_ZJ(c) = {qzc:.4f}; cota rho_J = {vb:.4f} (autor: 0,1486 no pior ponto da regiao 1)')
res['2b'] = dict(d=db, qzj=qzj, qjj=qjj, qzc=qzc, cota=vb)
json.dump(res, open(P/'r6_resultado.json', 'w'), indent=1)
# (3) direta
g = sl.campos()
modos = list(range(-47, -31))
dims = [sl.Radial(m, Nc).dim for m in modos]; off = np.cumsum([0] + dims); n = int(off[-1])
def w_modo(m):
    rr = sl.Radial(m, Nc)
    return K1**m*np.sqrt(om(np.concatenate([ns for _, ns in rr.comps])))
def Bw(mo, mi):
    return (w_modo(mo)[:, None]*sl.bloco_B(g, sl.Radial(mo, Nc), sl.Radial(mi, Nc)))/w_modo(mi)[None, :]
def X_de(s):
    X = np.zeros((n, n), complex)
    for b_, mi in enumerate(modos):
        Q, _ = Qmodo(mi, s)
        for a_, mo in enumerate(modos):
            if abs(mo - mi) <= 16:
                X[off[a_]:off[a_+1], off[b_]:off[b_+1]] = Bw(mo, mi) @ Q
    return X
Xc = X_de(1.0)
H = Xc.copy(); H[np.diag_indices(n)] += 1.0
lu = lu_factor(H, overwrite_a=True, check_finite=False); del H
log(f'(3) janela -47..-32 (dim {n}): LU de I + X(1) pronta')
res['3'] = {}
for s in (1.1 + 0.75j, 0.8 + 0.75j, 1.2 + 0.75j):
    D = Xc - X_de(s)
    op = LinearOperator((n, n), matvec=lambda v: D.conj().T @ lu_solve(lu, lu_solve(lu, D @ v), trans=2), dtype=complex)
    lam = eigsh(op, k=1, which='LM', tol=1e-8, return_eigenvectors=False)
    e = float(math.sqrt(abs(lam[0])))
    log(f'    s = {s}: ||(I + X(c))^-1 (X(c) - X(s))|| ~ {e:.4f}  (d = |s - c| = {abs(s - 1.0):.4f})')
    res['3'][str(s)] = e
    del D, op
    json.dump(res, open(P/'r6_resultado.json', 'w'), indent=1)
log('fim')
