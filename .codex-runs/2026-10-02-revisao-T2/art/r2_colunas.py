#!/usr/bin/env python3
"""R2: conferencia direta (ponto flutuante) das somas de colunas pesadas de Q = J_mu^-1 na norma l1 de RT
(pesos (2 - delta_n0) kappa2^n, kappa2 = 5/4; o peso de Fourier cancela, Q e diagonal em m), contra as
cotas de STRUCTURED_EXTERIOR: colunas n >= N: 9/(mu (n+1)) (parte J, k = 1) e (1 + 32n/(9(n-1)))/(mu(n+1))
(parte K, k = 2); modos |m| >= M: 1/(b0 (1 - r)) com b0 = M/2.  (Sem o fator 3/2 dos blocos reais.)"""
import numpy as np, math, sys, json
MU = 722873400/2**32
K2 = 5/4

def colsum(m, passo, cols, Ntot):
    ns = np.arange(Ntot) if passo == 1 else np.arange(1, Ntot, 2)
    sig = 1j*m/2
    d = MU*(ns + 1.0) + sig; v = 2*MU*ns
    lw = np.where(ns == 0, 0.0, math.log(2) + ns*math.log(K2))
    logrho = np.log(np.abs(1 - v/d)); logrho[0] = 0.0
    S = np.cumsum(logrho)
    out = []
    for c in cols:
        i = int(np.searchsorted(ns, c))
        j = np.arange(i + 1)
        prod = np.exp(np.where(j < i, S[i-1] - S[j], 0.0)) if i > 0 else np.ones(1)
        ent = (1/np.abs(d[j]))*(v[i]/np.abs(d[i]))*prod
        ent[i] = 1/np.abs(d[i])
        out.append(float(np.sum(ent*np.exp(lw[j] - lw[i]))))
    return np.array(out)

res = {}
N, M = 4096, 1024
cols = np.arange(N, N + 400, 1)
for m in (0, 1, 2, 7, 1023):
    for passo in (1, 2):
        cs = colsum(m, passo, cols if passo == 1 else cols[cols % 2 == 1], N + 401)
        nn = cols if passo == 1 else cols[cols % 2 == 1]
        cota = 9/(MU*(nn + 1)) if passo == 1 else (1 + 32*nn/(9*(nn - 1)))/(MU*(nn + 1))
        res[f'radial m={m} passo={passo}'] = dict(max_soma=float(cs.max()), max_razao=float((cs/cota).max()))
        print(f'colunas n >= {N}, m = {m}, passo {passo}: max soma {cs.max():.6f}; max soma/cota {(cs/cota).max():.4f}')
r = 4/5; b0 = M/2
for m in (1024, 1025, 3000):
    for passo in (1, 2):
        cols2 = np.arange(1 if passo == 2 else 0, 3000, 7)
        cols2 = cols2[cols2 % 2 == 1] if passo == 2 else cols2
        cs = colsum(m, passo, cols2, 3001)
        t = 1/(b0*(1 - r)) if passo == 1 else 1/(b0*(1 - r*r)) + 2*0.17*r*r/(b0*b0*(1 - r*r))
        res[f'fourier m={m} passo={passo}'] = dict(max_soma=float(cs.max()), cota=t)
        print(f'modo m = {m}, passo {passo}: max soma de coluna {cs.max():.6f}  cota {t:.6f}')
json.dump(res, open(sys.argv[1], 'w'), indent=1)
