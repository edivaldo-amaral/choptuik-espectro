#!/usr/bin/env python3
"""R4 (revisor): a cauda de 0,5i feita em duas maquinas e juntada.

(a) janelas em comum entre T3.parcial.json (Joao, ordem direta; = T3.joao.json no lab2) e T3.lab2.json (lab2,
    ordem reversa): todos os campos (linha N, TZ, nV, nVhi, rho, Y, K, nQm) devem coincidir;
(b) o T3.json juntado (lab2; = build/faixaB/cauda/T3_0_5.json) tem as 26 janelas e cada entrada vem de um dos
    checkpoints (linha N fora da diagonal, diagonal = rho, TZ, nV, nVhi, K, nQm, Y);
(c) a copia antiga de T3.lab2.json que esta no PC 3 (sha diferente) e um subconjunto da do lab2.
argv[1] = pasta com lab2/ e joao/ (copias)."""
import json, sys
from pathlib import Path
import numpy as np

S = Path(sys.argv[1])
RAIZ = Path('<repo>')
jo = json.loads((S/'joao/T3.parcial.json').read_text())
jo2 = json.loads((S/'lab2/T3.joao.json').read_text())
lb = json.loads((S/'lab2/T3.lab2.json').read_text())
lb_old = json.loads((S/'joao/T3.lab2.json').read_text())
fin = json.loads((S/'lab2/T3.json').read_text())
loc = json.loads((RAIZ/'build/faixaB/cauda/T3_0_5.json').read_text())
out = {}
print('T3.parcial (Joao) == T3.joao.json (lab2):', jo == jo2)
print('T3.json (lab2) == build/faixaB/cauda/T3_0_5.json:', fin == loc)
print(f'janelas: Joao {len(jo)}, lab2 {len(lb)}, copia antiga do lab2 no Joao {len(lb_old)}')
comuns = sorted(set(jo) & set(lb), key=lambda x: fin['nome'].index(x))
print(f'em comum ({len(comuns)}):', comuns)
campos = ['N', 'TZ', 'nV', 'nVhi', 'rho', 'Y', 'K', 'nQm']
difmax = 0.0; iguais = True
for w in comuns:
    for c in campos:
        a, b = np.array(jo[w][c], float), np.array(lb[w][c], float)
        if not np.array_equal(a, b):
            iguais = False
        d = float(np.max(np.abs(a - b)/np.maximum(np.abs(a), 1e-300))) if np.any(a != b) else 0.0
        difmax = max(difmax, d)
    print(f'  {w:12s}: nV {jo[w]["nV"]:.10f} / {lb[w]["nV"]:.10f}  TZ {jo[w]["TZ"]:.6g} / {lb[w]["TZ"]:.6g}  '
          f'K {jo[w]["K"]:.6g} / {lb[w]["K"]:.6g}  identicos: {all(jo[w][c] == lb[w][c] for c in campos)}')
print(f'(a) janelas em comum bit a bit iguais: {iguais}; maior diferenca relativa {difmax:.2e}')
out['a'] = dict(comuns=comuns, iguais=iguais, difmax=difmax)
# (b)
nomes = fin['nome']; k = len(nomes)
assert k == 26, k
N = np.array(fin['N'])
cob = set(jo) | set(lb)
falta = [w for w in nomes if w not in cob]
origem = {}; ruins = []
for i, w in enumerate(nomes):
    fontes = []
    for nome_ck, ck in (('joao', jo), ('lab2', lb)):
        if w not in ck:
            continue
        e = ck[w]
        Nr = np.array(e['N'], float); Nr[i] = e['rho']
        ok = (np.array_equal(Nr, N[i]) and e['TZ'] == fin['TZ'][i] and e['nV'] == fin['nV'][i] and
              e['nVhi'] == fin['nVhi'][i] and e['K'] == fin['K'][i] and e['nQm'] == fin['nQm'][i] and e['Y'] == fin['Y'][i])
        if ok:
            fontes.append(nome_ck)
    origem[w] = fontes
    if not fontes:
        ruins.append(w)
print(f'(b) T3.json: {k} janelas; sem checkpoint: {falta}; entradas que nao batem com nenhum checkpoint: {ruins}')
print('    origem: ' + ', '.join(f'{w}:{"+".join(v)}' for w, v in origem.items()))
out['b'] = dict(n=k, sem_checkpoint=falta, nao_batem=ruins, origem=origem)
# (c)
sub = all(w in lb and lb_old[w] == lb[w] for w in lb_old)
print(f'(c) a copia antiga de T3.lab2.json no Joao ({len(lb_old)} janelas) e subconjunto identico da do lab2: {sub}')
out['c'] = dict(subconjunto=sub, n_antiga=len(lb_old))
# outros campos globais do T3 juntado
print(f"lam0 {fin['lam0']}, rouche {fin['rouche']}, eps_fundo_L {fin['eps_fundo_L']}, nQZ {fin['nQZ']:.4f}, Qfar {fin['Qfar']:.4f}, eps_mu {fin['eps_mu']:.3e}")
json.dump(out, open(RAIZ/'.codex-runs/2026-10-02-revisao-B/r4_resultado.json', 'w'), indent=1)
