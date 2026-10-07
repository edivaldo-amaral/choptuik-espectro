#!/usr/bin/env python3
"""Controle negativo de r4_curvatura_rt.py: altera UM coeficiente de cada formula de RT e confirma que a
diferenca deixa de ser zero (o teste nao e vacuo). Requer sympy (PYTHONPATH do scratchpad)."""
import sys
import sympy as sp
import r4_curvatura_rt as R

x = R.x
alter = {
    'rr com 3*Omega3 em vez de 2*Omega3': (R.LHS['rr'], (R.A + R.B + R.Om_m['O3'] + R.Om_p['O3'])/(4*x)),
    '00 sem o termo -Omega3': (R.LHS['00'], R.RHS['00'] + (R.Om_m['O3'] + R.Om_p['O3'])/(2*x)),
    'box com fator 1/(4 xi)': (R.LHS['box'], (R.Om_m['O4'] + R.Om_p['O4'])/(4*x)),
    'theta trocando o sinal de Omega2sharp': (R.LHS['theta'], (R.A + 2*R.Om_m['O2s'] + 2*R.Om_p['O2s'])/(4*x)),
}
subs, pt = R.jet_subs(1)
nz = 0
for k, (l, r) in alter.items():
    d = R.evaluate(l - r, subs, pt)
    print(f'{k}: LHS - RHS alterado = {sp.N(d, 10)}')
    nz += d != 0
print('CONTROLE OK (todas as alteracoes detectadas)' if nz == len(alter) else 'CONTROLE FALHOU')
sys.exit(0 if nz == len(alter) else 1)
