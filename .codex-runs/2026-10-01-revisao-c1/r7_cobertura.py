#!/usr/bin/env python3
"""R7 (revisao independente C1): verificador PROPRIO de cobertura do contorno de
Omega = [0,3] x [-1/4,1/4] pelos discos certificados, em aritmetica racional exata.

Metodo diferente do autor: em vez de cotar sqrt por baixo, cada aresta e varrida por um ponto "x atual"
(racional) e, a cada passo, procura-se um disco que contenha o ponto (teste exato |z - p|^2 <= r^2) e
avanca-se ate um ponto racional DENTRO do disco obtido por bissecao exata (40 bissecoes) sobre a aresta.
Se nenhum disco contem o ponto atual, ha lacuna. Os 4 cantos sao testados explicitamente.
Tambem confere, por disco: ok == True, theta (string racional) < 1 recalculado como Fraction,
r > 0, r*t_p < 1, que o JSON registra certF, e que (p) e um ponto do contorno."""
import json, sys
from fractions import Fraction as Fr
from pathlib import Path

base = Path('build/rouche_L_rig/discos')
arqs = sorted(base.glob('disco_*.json'))
discos = []; problemas = []
for f in arqs:
    d = json.load(open(f))
    if not d.get('certF'):
        problemas.append(f'{f.name}: sem certF')
    for x in d['pontos']:
        th = Fr(x['theta'])
        px, py = Fr(x['p'][0]), Fr(x['p'][1]); r = Fr(x['r']); tp = Fr(x['t_p'])
        no_contorno = (px == 0 or px == 3) and -Fr(1, 4) <= py <= Fr(1, 4) or (py == Fr(1, 4) or py == -Fr(1, 4)) and 0 <= px <= 3
        if not (x['ok'] and th < 1 and r > 0 and r*tp < 1 and no_contorno):
            problemas.append(f'{f.name}: p={x["p"]} ok={x["ok"]} theta={float(th)} r={x["r"]} contorno={no_contorno}')
            continue
        discos.append((px, py, r, th, f.name))
print(f'{len(arqs)} arquivos, {len(discos)} discos validos; problemas: {problemas if problemas else "nenhum"}')
print(f'max theta (racional) = {float(max(d[3] for d in discos)):.6f} ({max(discos, key=lambda d: d[3])[4]}, p = {[float(v) for v in max(discos, key=lambda d: d[3])[:2]]})')

def dentro(z, dsc):
    px, py, r = dsc[:3]
    return (z[0] - px)**2 + (z[1] - py)**2 <= r*r

arestas = [((Fr(0), -Fr(1, 4)), (Fr(0), Fr(1, 4))), ((Fr(0), Fr(1, 4)), (Fr(3), Fr(1, 4))),
           ((Fr(3), Fr(1, 4)), (Fr(3), -Fr(1, 4))), ((Fr(3), -Fr(1, 4)), (Fr(0), -Fr(1, 4)))]
tudo = True
for a, b in arestas:
    ponto = lambda t: (a[0] + t*(b[0] - a[0]), a[1] + t*(b[1] - a[1]))
    t = Fr(0); passos = 0
    while t < 1:
        z = ponto(t)
        cands = [d for d in discos if dentro(z, d)]
        if not cands:
            print(f'  LACUNA em {float(z[0]):.6f}{float(z[1]):+.6f}i'); tudo = False; break
        melhor = t
        for dsc in cands:                 # maior t' tal que ponto(t') esta no disco (bissecao exata)
            if dentro(ponto(Fr(1)), dsc):
                melhor = Fr(1); break
            lo, hi = t, Fr(1)
            for _ in range(40):
                mid = (lo + hi)/2
                if dentro(ponto(mid), dsc):
                    lo = mid
                else:
                    hi = mid
            melhor = max(melhor, lo)
        if melhor == t:
            print(f'  sem avanco em t = {float(t)}'); tudo = False; break
        # o segmento [t, melhor] esta no disco (convexo, extremos dentro); o proximo ponto e 'melhor'
        t = melhor; passos += 1
    print(f'  aresta {[float(v) for v in a]} -> {[float(v) for v in b]}: {"coberta" if t >= 1 else "NAO coberta"} em {passos} passos')
cantos = [(Fr(0), -Fr(1, 4)), (Fr(0), Fr(1, 4)), (Fr(3), Fr(1, 4)), (Fr(3), -Fr(1, 4))]
for z in cantos:
    ds = [d[4] + f' p=({float(d[0])},{float(d[1])})' for d in discos if dentro(z, d)]
    print(f'  canto ({z[0]}, {z[1]}): {len(ds)} discos: {ds[:3]}')
    tudo &= bool(ds)
print('CONTORNO COBERTO (verificador proprio)' if tudo else 'CONTORNO NAO COBERTO')
