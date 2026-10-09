#!/usr/bin/env python3
"""N4: ||Q0 P_far|| refeito das formas fechadas para este lambda0 (rodar da raiz do projeto)."""
import sys, json; sys.path.insert(0, 'scripts')
import closed_form as cf, signed_operator_L as sl
lam0 = 0.7331796596606924; Mc, Nc = 200, 400
Qf = max(max(cf.cauda_Rinv(lam0 + 1j*m/2, Nc, p, sl.MU, sl.K2, 1500) for p in ((2, 1, 1) if m % 2 == 0 else (1,))) for m in range(-(Mc - 1), Mc))
U = cf.Rinv_uniforme(Mc, sl.MU, lam0, sl.K2)
z = json.load(open('build/nk733/Z.json'))
print('Qfar radial', repr(Qf), ' uniforme |m|>=200', repr(U), ' max', repr(max(Qf, U)), ' Z.json', repr(z['Qfar']))
