#!/usr/bin/env python3
"""N3: a constante NORMA_BL = 3,964043 (nk_L.py) usada em beta, desc, epsZ e nas descidas do NK de L,
contra o certificado do símbolo de L (build/spectrum/symbol-L-faixa*.json, grade 2048x4096) com a folga refeita
em Arb (symbol_slack_arb.folga_arb). Imprime t_norma + folga (dados RefA) e + 40 eps (correção genérica)."""
import json, sys
from fractions import Fraction as Q
from pathlib import Path
sys.path.insert(0, 'scripts')
import certify_symbol_numerical_range as cs
import symbol_slack_arb as ss
from flint import arb
cfg = cs.OPER['L']
_, lips = cs.carrega(cfg['k1'], cfg['k2'])
fs = [json.loads(Path(f'build/spectrum/symbol-L-faixa{i}.json').read_text()) for i in range(4)]
assert all(f['falhas'] == 0 and f['t_norma'] == fs[0]['t_norma'] for f in fs)
cob = sorted(tuple(f['linhas']) for f in fs); print('faixas', cob, 'nphi', fs[0]['nphi'], 'ntheta', fs[0]['ntheta'])
fl, fn = ss.folga_arb(cfg['forma'], cfg['dim'], lips, fs[0]['ntheta'], fs[0]['nphi'], cfg['k2'], cs.MUS)
tn = ss.A(fs[0]['t_norma'])
dados = (tn + arb(fn)).upper(); corr = fs[0]['correcao']
print(f't_norma = {fs[0]["t_norma"]}; folga Arb (norma) = {fn!r}')
print(f'||B_L(RefA)|| <= t + folga = {float(dados):.12f}')
print(f'  + 40 eps_omega ({corr:.4e}) = {float(dados) + corr:.12f}')
print(f'  + 4,7e-8 (CHOPTUIK_EPS_FUNDO_L) = {float(dados) + 4.7e-8:.12f}')
N = 3.964043
print(f'NORMA_BL = {N}: NORMA_BL >= t + folga? {Q(N) >= Q(float(dados))} (folga {N - float(dados):.3e});'
      f'  NORMA_BL + 4,7e-8 >= t + folga + 4,7e-8 ? sim se o anterior')
print(f'  NORMA_BL >= t + folga + 40 eps? {N >= float(dados) + corr} (diferenca {N - float(dados) - corr:.3e})')
