#!/usr/bin/env python3
"""R2/R3/R7 (revisor, codigo proprio): contas densas na A.npy do laboratorio2.

R2: para pontos p perto de 0,5i, ||(A + p)^-1|| por LU densa (scipy) + Lanczos (eigsh) no operador
    (A + p)^-* (A + p)^-1, e uma cota INFERIOR independente do solver: ||x|| / ||(A + p) x|| com x = LU^-1 v
    (o produto (A + p) x e feito com a A original, por matvec).
R3: no primeiro ponto, a0 ~ ||J_Z(c) [(A + p)^-1 F_t - Q_ZZ(p) F_b]|| com LU (e nao Schur), J_Z e Q_ZZ montadas
    por codigo proprio (pesos kappa1^m sqrt(omega(n)), J = diag mu(n+1) + s + i m/2 + 2 mu n acima).
R7: autovalores de A em -Omega_B = (-3, 0) x (-3/4, -1/4) por shift-invert (ARPACK) em dois deslocamentos.
Tudo em ponto flutuante: e estimativa, nao certificado. Grava JSON parcial a cada passo."""
import hashlib, json, sys, time
from pathlib import Path
import numpy as np
from scipy.linalg import lu_factor, lu_solve, solve_triangular
from scipy.sparse.linalg import LinearOperator, eigsh, eigs

D = Path.home()/'Choptuik/build/rouche_L_v2'
OUT = Path.home()/'revisao_B/r2_resultado.json'
MU = 722873400/2**32; K1 = 65/64; K2 = 5/4; MZ, NZ = 32, 128
t0 = time.time()
res = {}
def log(*x):
    print(*x, f'[{time.time()-t0:.0f}s]', flush=True)
def salva():
    OUT.write_text(json.dumps(res, indent=1) + '\n')

def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 24), b''):
            h.update(b)
    return h.hexdigest()

def comps(m, N):
    if m % 2 == 0:
        return [np.arange(1, N, 2), np.arange(N), np.arange(N)]
    return [np.arange(N)]

def omega(n):
    return np.where(n == 0, 1.0, K2**(2.0*n) + K2**(-2.0*n))

def blocos_J(s):
    """lista por modo (m = -31..31) de listas por componente de J pesada (sqrt omega)."""
    out = []
    for m in range(-(MZ - 1), MZ):
        bl = []
        for ns in comps(m, NZ):
            J = np.diag(MU*(ns + 1.0) + s + 1j*m/2).astype(complex)
            J += np.triu(np.tile(2*MU*ns, (len(ns), 1)), 1)
            w = np.sqrt(omega(ns))
            bl.append((w[:, None]*J)/w[None, :])
        out.append(bl)
    return out

def aplica(bl, X, inv=False):
    Y = np.empty_like(X); o = 0
    for bm in bl:
        for Jc in bm:
            k = Jc.shape[0]
            Y[o:o+k] = solve_triangular(Jc, X[o:o+k]) if inv else Jc @ X[o:o+k]
            o += k
    assert o == X.shape[0]
    return Y

shaA = sha(D/'A.npy'); res['sha256_A'] = shaA; log('sha A', shaA)
assert shaA.startswith('70791ce0')
A0 = np.load(D/'A.npy'); n = A0.shape[0]; assert n == 14016
pontos = [complex(x) for x in sys.argv[1].split(',')] if len(sys.argv) > 1 else [0.502449j, 0.47849j, 0.464083j]
res['R2'] = []
for ip, p in enumerate(pontos):
    M = A0.copy(); M[np.diag_indices(n)] += p
    lu = lu_factor(M, overwrite_a=True, check_finite=False); del M
    log('LU', p)
    op = LinearOperator((n, n), matvec=lambda v: lu_solve(lu, lu_solve(lu, v), trans=2), dtype=complex)
    lam, V = eigsh(op, k=3, which='LM', tol=1e-12)
    o = np.argsort(-lam); lam = lam[o]; V = V[:, o]
    sv = np.sqrt(np.abs(lam))
    v = V[:, 0]; x = lu_solve(lu, v)
    Ax = A0 @ x + p*x
    inf = float(np.linalg.norm(x)/np.linalg.norm(Ax))
    resid = float(np.linalg.norm(Ax - v)/np.linalg.norm(v))
    r = dict(p=[p.real, p.imag], sigma_max_inv=sv.tolist(), cota_inferior=inf, residuo_solve=resid)
    log(f'p = {p}: ||(A+p)^-1|| ~ {sv[0]:.6f} (2o, 3o: {sv[1]:.3f}, {sv[2]:.3f}); cota inferior ||x||/||(A+p)x|| = {inf:.6f}; residuo {resid:.1e}')
    if ip == 0:   # R3: a0 com LU
        c = 0.5j
        F = np.load(D/'F_+0_0000+0_5000i.npy'); res['sha256_F'] = sha(D/'F_+0_0000+0_5000i.npy')
        Ft, Fb = F[:n], F[n:]
        X1 = lu_solve(lu, Ft)
        QF = aplica(blocos_J(p), Fb, inv=True)
        Mt = aplica(blocos_J(c), X1 - QF)
        G = Mt.conj().T @ Mt
        a0 = float(np.sqrt(np.linalg.eigvalsh(G).max()))
        # residuo da solucao LU
        rr = float(np.linalg.norm(A0 @ X1 + p*X1 - Ft)/np.linalg.norm(Ft))
        r['a0_LU'] = a0; r['residuo_X1'] = rr; r['k_F'] = int(F.shape[1])
        log(f'R3: a0 (LU, codigo proprio) = {a0:.6f}; residuo relativo de X1 {rr:.1e}; F com {F.shape[1]} colunas')
        del F, Ft, Fb, X1, QF, Mt
    res['R2'].append(r); salva()
    del lu, op
# R7: autovalores de A em -Omega_B
res['R7'] = []
for sig in (-0.75 - 0.5j, -2.25 - 0.5j):
    M = A0.copy(); M[np.diag_indices(n)] -= sig
    lu = lu_factor(M, overwrite_a=True, check_finite=False); del M
    op = LinearOperator((n, n), matvec=lambda v: lu_solve(lu, v), dtype=complex)
    nu = eigs(op, k=30, which='LM', tol=1e-10, return_eigenvectors=False)
    lam = sig + 1/nu
    dist = np.abs(lam - sig)
    o = np.argsort(dist); lam = lam[o]; dist = dist[o]
    dentro = [complex(z) for z in lam if -3 < z.real < 0 and -0.75 < z.imag < -0.25]
    r = dict(sigma=[sig.real, sig.imag], lam=[[z.real, z.imag] for z in lam], dist_max=float(dist.max()),
             dentro=[[z.real, z.imag] for z in dentro])
    log(f'R7 sigma = {sig}: {len(dentro)} autovalores em -Omega_B: {dentro}; raio coberto {dist.max():.3f}')
    log('   mais proximos: ' + ', '.join(f'{z.real:.5f}{z.imag:+.5f}i' for z in lam[:12]))
    res['R7'].append(r); salva()
    del lu, op
res['fim'] = True; salva(); log('fim')
