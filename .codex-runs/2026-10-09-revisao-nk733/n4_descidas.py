#!/usr/bin/env python3
"""N4: epsZ, desc e qZf refeitos das formas fechadas com os dados de Z.json (rodar da raiz do projeto)."""
import sys, json; sys.path.insert(0, 'scripts')
import closed_form as cf, signed_operator_L as sl
z = json.load(open('build/nk733/Z.json')); t = json.load(open('build/nk733/T3.json'))
NB = 3.964043; MU, K2 = sl.MU, sl.K2; NZ, Nc = 128, 400; lim = Nc - 100 - 60
fb = z['dB_L']; emu = z['eps_mu']; Qf = z['Qfar']; tV = z['tV']
nBL = NB + fb['total']
epsZ = tV*(nBL*(cf.descida(NZ + 101, Nc, MU, K2) + 6*cf.descida_divxi(NZ, Nc, K2)) + (fb['total'] + emu)*Qf)*(1 + 1e-12)
desc = nBL*(cf.descida(lim + 101, Nc, MU, K2) + 6*cf.descida_divxi(lim, Nc, K2)) + (fb['total'] + emu)*Qf
qZf = cf.descida(NZ, Nc, MU, K2)
print(f'epsZ refeito {epsZ!r}  gravado {z["epsZ"]!r}  igual {epsZ == z["epsZ"]}')
print(f'desc refeito {desc!r}  gravado {t["desc"]!r}  igual {desc == t["desc"]}')
print(f'qZf  refeito {qZf!r}  gravado {json.load(open("build/nk733/Q.json"))["qZf"]!r}')
