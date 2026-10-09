#!/usr/bin/env python3
"""N4 (Q.json e constantes de Q0): recálculo em ponto flutuante, com SVD, das normas restritas de Q0 pesado
(q_ZT, q_TT, max ||Q0|| em 48 <= |m| < 72, ||Q_ZZ||) para TODOS os modos |m| < 200, e comparação com o gravado.
Só a montagem de J (sl.J_livre) e os pesos vêm do código do autor; inversão e normas são do revisor."""
import json, sys
sys.path.insert(0, 'scripts')
import numpy as np
import signed_operator_L as sl
import closed_form as cf
lam0 = 0.7331796596606924; MZ, NZ, Mc, Nc = 32, 128, 200, 400
qj = json.load(open('build/nk733/Q.json')); z = json.load(open('build/nk733/Z.json')); t = json.load(open('build/nk733/T3.json'))
def Qw(m, N):
    r = sl.Radial(m, N); n = np.concatenate([ns for _, ns in r.comps])
    w = np.sqrt(cf.omega(n, sl.K2))
    J = sl.J_livre(r, lam0)
    return (w[:, None]*np.linalg.inv(J))/w[None, :], n
qZT = qTT = nQx = nQZ = 0.0; nQmax = 0.0; argZT = argTT = None
for m in range(-(Mc - 1), Mc):
    Q, n = Qw(m, Nc)
    colT = np.where(n >= NZ)[0] if abs(m) < MZ else np.arange(len(n))
    if abs(m) < MZ:
        v = np.linalg.norm(Q[np.ix_(np.where(n < NZ)[0], colT)], 2)
        if v > qZT: qZT, argZT = v, m
        QZ, _ = Qw(m, NZ); nQZ = max(nQZ, np.linalg.norm(QZ, 2))
    v = np.linalg.norm(Q[np.ix_(colT, colT)], 2)
    if v > qTT: qTT, argTT = v, m
    nm = np.linalg.norm(Q, 2); nQmax = max(nQmax, nm)
    if MZ + 16 <= abs(m) < MZ + 40:
        nQx = max(nQx, nm)
res = dict(qZT=(qZT, qj['qZT'], argZT), qTT=(qTT, qj['qTT'], argTT), nQx_extra=(nQx, qj['nQx_extra']),
           nQZ=(nQZ, z['nQZ']), max_nQ_F=(nQmax, max(t['nQ'])))
for k, v in res.items():
    print(f'{k:10s} float {v[0]:.10f}   gravado {v[1]:.10f}   gravado >= float: {v[1] >= v[0]}   razao {v[1]/v[0]:.8f}' + (f'   (modo {v[2]})' if len(v) > 2 else ''))
print('qZf gravado', qj['qZf'], '; Qfar gravado', z['Qfar'])
