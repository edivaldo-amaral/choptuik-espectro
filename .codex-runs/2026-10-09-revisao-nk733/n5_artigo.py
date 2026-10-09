#!/usr/bin/env python3
"""N5: a linha de s* da Proposição 6.2 (paper/sec/roots.tex) e o intervalo do Teorema Principal
(paper/sec/main_theorem.tex) contra os racionais de build/nk733/C3.json, com lambda0 = o double EXATO."""
import json, re
from fractions import Fraction as Fr
from decimal import Decimal, getcontext
from pathlib import Path
getcontext().prec = 40
RAIZ = Path('<repo>')
c = json.loads((RAIZ/'build/nk733/C3.json').read_text())
roots = (RAIZ/'paper/sec/roots.tex').read_text(); mt = (RAIZ/'paper/sec/main_theorem.tex').read_text()
linha = [l for l in roots.splitlines() if l.startswith('$\\sphys$ & $0.7331796596606924$')][0]
print('artigo, roots.tex:', linha)
cel = [x.strip().strip('$\\').replace('\\\\', '').strip() for x in linha.split('&')]
lam0 = Fr(c['lambda0'])                              # double exato
print('lambda0 (double exato) =', Decimal(lam0.numerator)/Decimal(lam0.denominator))
def num(s):
    s = s.replace('\\cdot10^{', 'e').replace('}', '').replace('$', '').strip()
    return Fr(Decimal(s))
pares = [('theta', cel[2]), ('Z2', cel[3]), ('Y', cel[4]), ('r', cel[5]), ('erro_lambda', cel[6])]
res = {}
for chave, txt in pares:
    v = Fr(c[chave]); a = num(txt)
    ok = a >= v
    res[chave] = dict(certificado=float(v), artigo=txt, artigo_ge_certificado=ok, folga_rel=float((a - v)/v))
    print(f'  {chave:12s} certificado {float(v):.10e}  artigo {txt:>18s}  artigo >= certificado: {ok}  (folga relativa {float((a - v)/v):.2e})')
print('  lambda0 do artigo == C3.lambda0:', Fr(Decimal(cel[1])) == lam0, '(decimal do artigo vs double:',
      float(Fr(Decimal(cel[1])) - lam0), ')')
m = re.search(r'\\sphys\\in\[([\d.]+),\\ ([\d.]+)\]', mt)
lo, hi = Fr(Decimal(m.group(1))), Fr(Decimal(m.group(2)))
e_art = Fr(Decimal('1.0163e-5')); e_cert = Fr(c['erro_lambda'])
print(f'Teorema Principal: [{m.group(1)}, {m.group(2)}]')
for nome, e in (('erro do certificado', e_cert), ('1,0163e-5 do artigo', e_art)):
    a, b = lam0 - e, lam0 + e
    print(f'  com {nome}: [{float(a):.13f}, {float(b):.13f}];  lo <= a: {lo <= a} (folga {float(a - lo):.2e});'
          f'  b <= hi: {b <= hi} (folga {float(hi - b):.2e})')
    res[nome] = dict(lo_ok=lo <= a, hi_ok=b <= hi, folga_lo=float(a - lo), folga_hi=float(hi - b))
# arredondamento "para fora com oito dígitos": o menor intervalo de 8 casas que contém [a, b]
import math
a, b = lam0 - e_art, lam0 + e_art
lo8 = Fr(math.floor(a*10**8), 10**8); hi8 = Fr(math.ceil(b*10**8), 10**8)
print(f'  arredondamento para fora de lambda0 -+ 1,0163e-5 em 8 casas: [{float(lo8):.8f}, {float(hi8):.8f}] == artigo: {(lo8, hi8) == (lo, hi)}')
res['intervalo_8_casas_igual'] = (lo8, hi8) == (lo, hi)
Path(__file__).with_suffix('.json').write_text(json.dumps(res, indent=1, default=str) + '\n')
