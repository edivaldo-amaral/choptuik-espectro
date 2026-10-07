#!/usr/bin/env python3
"""R2 (revisao independente C1): valores VERDADEIROS (ponto flutuante, LU densa) de
  ||(A + s)^-1||_2  e  ||J_Z(c)[R(s) F_t - Q_ZZ(s) F_b]||_2
em s do contorno dentro de discos |s - p| <= r do centro c = 0,25i, comparados com as cotas
  t(s) = t_p/(1 - |s - p| t_p)  e  E(s) = a0 + e b1 + e^2 b2 + e^3 (tau ||Y3|| + JQ ||Q^3 F_b||)
(e = |s - p|; tau = 1 + (|s - c| + ||B^||) t(s); JQ = 1 + |s - c| nQZ(p)/(1 - e nQZ(p))).
Tambem ||J_Z(c) R(s)||_2 (contra tau) e, em s = p, o valor contra a0 (cota certificada).
Codigo proprio: Q_ZZ(s) por inversao do J pesado de cada modo; a LU e a do scipy (estimativa,
nao certificado: o objetivo e flagrar uma cota MENOR que o valor verdadeiro).
Uso: python r2_lu_densa.py --dir <com A.npy, F_+0_0000+0_2500i.npy> --rig <rig json> --disco <disco json>"""
import argparse, json, math, sys, time, hashlib
from pathlib import Path
import numpy as np
from scipy.linalg import lu_factor, lu_solve
from scipy.sparse.linalg import LinearOperator, eigsh

sys.path.insert(0, str(Path.home()/'Choptuik/scripts'))
import signed_operator_L as sl   # so a DEFINICAO de J_livre e Radial
import closed_form as cf          # so omega

MZ, NZ = 32, 128
K2 = 5/4
t0 = time.time()
def log(*x):
    print(*x, f'[{time.time()-t0:.0f}s]', flush=True)

def blocos(s):
    """J pesado (sqrt omega; kappa1^m cancela no bloco diagonal) de cada modo |m| < MZ, n < NZ."""
    out = []
    for m in range(-(MZ - 1), MZ):
        r = sl.Radial(m, NZ)
        w = np.sqrt(np.concatenate([cf.omega(ns, K2) for _, ns in r.comps]))
        out.append((w[:, None]*sl.J_livre(r, s))/w[None, :])
    return out

def aplica(Bs, X):
    out = np.empty_like(X); o = 0
    for M in Bs:
        k = M.shape[0]; out[o:o+k] = M @ X[o:o+k]; o += k
    return out

def norma2(M):
    G = M.conj().T @ M if M.shape[0] >= M.shape[1] else M @ M.conj().T
    return float(np.sqrt(np.max(np.linalg.eigvalsh(G))))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dir', type=Path, required=True)
    ap.add_argument('--rig', type=Path, required=True)
    ap.add_argument('--disco', type=Path, required=True)
    ap.add_argument('--ip', default='0,8', help='indices dos pontos p do rig')
    ap.add_argument('--saida', type=Path, required=True)
    ap.add_argument('--c', type=complex, default=0.25j)
    ap.add_argument('--F', default='F_+0_0000+0_2500i.npy')
    ap.add_argument('--extra', default='', help='ip:s;ip:s pontos s extras')
    a = ap.parse_args()
    c = a.c
    rig = json.load(open(a.rig)); dis = json.load(open(a.disco)); nB = 3.850313347957818
    with open(a.dir/'A.npy', 'rb') as fh:
        log('sha256 A', hashlib.sha256(fh.read()).hexdigest())
    with open(a.dir/a.F, 'rb') as fh:
        log('sha256 F', hashlib.sha256(fh.read()).hexdigest())
    A = np.load(a.dir/'A.npy'); n = A.shape[0]
    F = np.load(a.dir/a.F); Ft, Fb = F[:n], F[n:]; del F
    assert abs(complex(*rig['c']) - c) < 1e-12 and abs(complex(*dis['c']) - c) < 1e-12
    Jc = blocos(c)
    res = json.load(open(a.saida)) if a.saida.exists() else []
    for ip in map(int, a.ip.split(',')):
        zp = rig['pontos'][ip]; dp = dis['pontos'][ip]
        p = complex(*zp['p']); assert abs(complex(*dp['p']) - p) < 1e-12
        r = dp['r']; tp = zp['t_p']
        # s no contorno dentro do disco: bordo esquerdo (Re s = 0) e superior (Im s = 1/4)
        ss = [p]
        if p.real == 0:
            for f in (0.999, 0.5, -0.5, -0.999):
                y = p.imag + f*r
                if -0.25 <= y <= 0.25: ss.append(complex(0, y))
        if abs(p.imag - 0.25) < 1e-12:
            for f in (0.999, 0.5):
                ss.append(complex(p.real + f*r, 0.25))
        for item in [x for x in a.extra.split(';') if x]:
            i_, s_ = item.split(':')
            if int(i_) == ip: ss.append(complex(s_))
        ss = list(dict.fromkeys(ss))
        for s in ss:
            e = abs(s - p); d = abs(s - c)
            As = A.copy(); As[np.diag_indices(n)] += s
            lu = lu_factor(As, overwrite_a=True, check_finite=False); del As
            op = LinearOperator((n, n), matvec=lambda v: lu_solve(lu, lu_solve(lu, v), trans=2), dtype=complex)
            nR = float(np.sqrt(np.max(np.abs(eigsh(op, k=1, which='LM', tol=1e-10, return_eigenvectors=False)))))
            # ||J_Z(c) R(s)|| : maior autovalor de (J R)^* (J R)
            JcH = [M.conj().T for M in Jc]
            op2 = LinearOperator((n, n), matvec=lambda v: lu_solve(lu, aplica(JcH, aplica(Jc, lu_solve(lu, v))), trans=2), dtype=complex)
            nJR = float(np.sqrt(np.max(np.abs(eigsh(op2, k=1, which='LM', tol=1e-8, return_eigenvectors=False)))))
            G = lu_solve(lu, Ft, check_finite=False); del lu
            Qs = [np.linalg.inv(M) for M in blocos(s)]
            M = aplica(Jc, G - aplica(Qs, Fb)); del G
            val = norma2(M); del M
            ts = tp/(1 - e*tp) if e*tp < 1 else float('inf')
            nQZp = zp['nQZ']; nQZs = nQZp/(1 - e*nQZp)
            tau = 1 + (d + nB)*ts; JQ = 1 + d*nQZs
            E = zp['a0'] + e*zp['b1'] + e*e*zp['b2'] + e**3*(tau*zp['nY3'] + JQ*zp['nQ3F'])
            # tau e d como na montagem (d = |p - c| + r, raio inteiro): a cota do disco inteiro
            dd = abs(p - c) + r; tsr = tp/(1 - r*tp); taur = 1 + (dd + nB)*tsr
            Er = zp['a0'] + r*zp['b1'] + r*r*zp['b2'] + r**3*(taur*zp['nY3'] + (1 + dd*nQZp/(1 - r*nQZp))*zp['nQ3F'])
            ok = nR <= ts and val <= E and nJR <= tau
            reg = dict(p=[p.real, p.imag], s=[s.real, s.imag], e=e, r=r, R_verd=nR, t_s=ts, JR_verd=nJR, tau_s=tau,
                       tau_disco=taur, val_verd=val, E_s=E, E_disco=Er, a0=zp['a0'], ok=ok)
            res.append(reg)
            log(f'p={p:.4f} s={s:.5f} |s-p|={e:.2e}: ||R(s)|| {nR:.3f} <= t(s) {ts:.2f}; ||J_Z(c)R(s)|| {nJR:.2f} <= tau {tau:.1f} (disco {taur:.1f});'
                f' ||J(c)[RF_t - QF_b]|| {val:.6f} <= E(s) {E:.6f} (disco {Er:.6f}; a0 {zp["a0"]:.6f})' + ('' if ok else '  *** VIOLACAO ***'))
            json.dump(res, open(a.saida, 'w'), indent=1)
    log('fim')

if __name__ == '__main__':
    main()
