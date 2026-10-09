#!/usr/bin/env python3
"""N3 (estrutura usada em Z2 e nos níveis): J_mu + lambda0 é diagonal em m e 'triangular superior em n'
(J[i, j] = 0 se n_i > n_j), logo Q0 também; e o bloco de Q0 em n < NZ é o inverso do bloco truncado.
Checagem numérica em alguns modos (código do autor só para montar J; a checagem é do revisor)."""
import sys
sys.path.insert(0, 'scripts')
import numpy as np
import signed_operator_L as sl
import closed_form as cf
lam0 = 0.7331796596606924
for m in (0, 1, 5, 31, 32, 48, 100, 199, -7, -150):
    for N in (128, 400):
        r = sl.Radial(m, N)
        n = np.concatenate([ns for _, ns in r.comps])
        J = sl.J_livre(r, lam0)
        viol = np.abs(J[n[:, None] > n[None, :]]).max() if np.any(n[:, None] > n[None, :]) else 0.0
        Q = np.linalg.inv(J)
        violQ = np.abs(Q[n[:, None] > n[None, :]]).max()
        print(f'm {m:5d} N {N}: dim {r.dim:5d}; max |J| abaixo (n_i > n_j) = {viol:.1e}; max |Q| abaixo = {violQ:.1e}', end='')
        if N == 400:
            r1 = sl.Radial(m, 128); n1 = np.concatenate([ns for _, ns in r1.comps])
            # bloco n < 128 de Q(400) contra inv(J(128)), mesmos índices (componente, n)
            pos = {(c, int(k)): i for i, (c, k) in enumerate((c, k) for c, ns in r.comps for k in ns)}
            idx = [pos[(c, int(k))] for c, ns in r1.comps for k in ns]
            Q1 = np.linalg.inv(sl.J_livre(r1, lam0))
            print(f';  ||Q(400)[n<128] - inv(J(128))|| / ||.|| = {np.linalg.norm(Q[np.ix_(idx, idx)] - Q1)/np.linalg.norm(Q1):.1e}', end='')
        print()
