#!/usr/bin/env python3
"""N4: tabela float (ARPACK) x rigoroso (T3.json) das janelas refeitas, com a perturbação global esperada."""
import json
from pathlib import Path
R = Path(__file__).parent
f = json.load(open(R/'lab2/n4_janelas_float.json')); T = json.load(open('<repo>/build/nk733/T3.json'))
pert = T['pert']; nQ = dict(zip(T['nome'], T['nQ'])); nQZ_F = max(nQ[n] for n in ('-15..0', '1..15'))
print('janela        qtd          float(ARPACK)    rigoroso(T3)     rig - float      nV*pert*nQ(col)  rig>=float')
for nm, r in f.items():
    g = r['gravado_T3']; nV = g['nV']
    linhas = [('nV', r['nV'], g['nV'], None), ('nVhi', r['nVhi'], g['nVhi'], None)]
    for c, v in r['N'].items():
        linhas.append(('N<-' + c, v, g['N'][c], nV*pert*nQ[c]))
    if 'TZ' in r:
        linhas.append(('TZ', r['TZ'], g['TZ'], nV*pert*nQZ_F))
    if 'Y_sem_erro' in r:
        linhas.append(('Y', r['Y_sem_erro'], g['Y'], nV*r['Y_erro_x']))
    for q, a, b, p in linhas:
        print(f'{nm:12s}  {q:12s} {a:.10e} {b:.10e} {b - a:+.4e}  ' + (f'{p:.4e}' if p is not None else ' ' * 10) + f'   {b >= a}')
    print(f'{nm:12s}  controle negativo Loewner ||V|| em t = 0,999*ARPACK = {r["loewner_nV_controle_negativo"]["t"]:.6f}: '
          + ('verificado (MAU)' if r['loewner_nV_controle_negativo']['verificado'] else 'rejeitado (esperado)'))
