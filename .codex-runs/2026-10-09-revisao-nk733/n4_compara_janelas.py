#!/usr/bin/env python3
"""N4: compara as janelas refeitas no laboratorio2 (nk_L3_T.py, código do snapshot, h0w gravado) com T3.json."""
import json
from pathlib import Path
R = Path(__file__).parent
ck = json.loads((R/'lab2/n4_ckpt.json').read_text()); T = json.loads(Path('<repo>/build/nk733/T3.json').read_text())
for nm in ('-15..0', '32..47', '-191..-176'):
    i = T['nome'].index(nm); f = ck[nm]
    lin = list(T['N'][i]); rho = lin[i]; lin[i] = 0.0
    pares = dict(nV=(f['nV'], T['nV'][i]), nVhi=(f['nVhi'], T['nVhi'][i]), rho=(f['rho'], rho), TZ=(f['TZ'], T['TZ'][i]),
                 Y=(f['Y'], T['Y'][i]))
    for j in range(26):
        if lin[j] or f['N'][j]:
            pares['N<-' + T['nome'][j]] = (f['N'][j], lin[j])
    ig = sum(a == b for a, b in pares.values())
    print(f'{nm}: {ig}/{len(pares)} quantidades IDENTICAS bit a bit; maior |dif| relativa '
          f'{max(abs(a - b)/max(abs(b), 1e-300) for a, b in pares.values()):.2e}')
    for k, (a, b) in pares.items():
        print(f'    {k:18s} refeito {a!r:24s} gravado {b!r:24s} {"=" if a == b else "DIFERE"}')
