#!/usr/bin/env python3
"""R7 (revisor): verificador PROPRIO de cobertura de d([0,3] x [1/4,3/4]) pelos discos certificados.

Sem raiz quadrada: cada aresta [0, L] e percorrida por uma cadeia de pontos RACIONAIS 0 = t_0 < t_1 < ... < t_m = L
tal que t_i e t_{i+1} estao (exatamente, |z - p|^2 <= r^2 em racionais) no MESMO disco; pela convexidade do disco
o segmento [t_i, t_{i+1}] fica coberto. Os cantos sao t_0 e t_m de cada aresta. Discos: ok e r > 0 dos JSONs dados
(faixa A em build/rouche_L_rig/discos_v2 e faixa B em build/faixaB/discos). Tambem: max t(s) = t_p/(1 - r t_p) so nos
discos usados na cadeia, e t eps_B com o eps_B de contagem.json."""
import glob, json, math, sys
from fractions import Fraction as Fr
from pathlib import Path
RAIZ = Path('<repo>')
arqs = sorted(glob.glob(str(RAIZ/'build/rouche_L_rig/discos_v2/*.json'))) + sorted(glob.glob(str(RAIZ/'build/faixaB/discos/*.json')))
discos = []
for f in arqs:
    d = json.load(open(f))
    assert d.get('certF'), f
    for x in d['pontos']:
        if x['ok'] and x['r'] > 0:
            assert float(x['theta_f']) < 1 and Fr(x['theta']) < 1
            discos.append((Fr(x['p'][0]), Fr(x['p'][1]), Fr(x['r']), Fr(x['t_p']), Path(f).parent.name + '/' + Path(f).name))
print(f'{len(discos)} discos (ok, r > 0) em {len(arqs)} arquivos')
re0, re1, im0, im1 = Fr(0), Fr(3), Fr(1, 4), Fr(3, 4)
arestas = [('inferior Im=1/4', (re0, im0), (1, 0), re1 - re0), ('direita Re=3', (re1, im0), (0, 1), im1 - im0),
           ('superior Im=3/4', (re0, im1), (1, 0), re1 - re0), ('esquerda Re=0', (re0, im0), (0, 1), im1 - im0)]
def dentro(D, z):
    return (z[0] - D[0])**2 + (z[1] - D[1])**2 <= D[2]**2
usados = set(); ok_tudo = True; relato = []
for nome, a, u, L in arestas:
    pt = lambda t: (a[0] + t*u[0], a[1] + t*u[1])
    t = Fr(0); cadeia = 0
    while True:
        cands = [D for D in discos if dentro(D, pt(t))]
        if not cands:
            print(f'  LACUNA na aresta {nome} em t = {float(t):.6f}'); ok_tudo = False; break
        # para cada candidato: maior t' racional (bissecao em racionais) com pt(t') ainda no disco
        melhor = None
        for D in cands:
            if dentro(D, pt(L)):
                melhor = (L, D); break
            lo, hi = t, L
            for _ in range(60):
                mid = (lo + hi)/2
                mid = Fr(math.floor(mid*10**15), 10**15) if mid.denominator > 10**15 else mid
                if mid <= lo:
                    break
                if dentro(D, pt(mid)):
                    lo = mid
                else:
                    hi = mid
            if melhor is None or lo > melhor[0]:
                melhor = (lo, D)
        t2, D = melhor
        usados.add(D); cadeia += 1
        if t2 <= t:
            print(f'  SEM AVANCO na aresta {nome} em t = {float(t):.6f}'); ok_tudo = False; break
        t = t2
        if t == L:
            break
    relato.append((nome, cadeia))
    print(f'  aresta {nome}: coberta por uma cadeia de {cadeia} discos (cantos incluidos)' if t == L else f'  aresta {nome}: FALHA')
print('CONTORNO COBERTO (verificador do revisor)' if ok_tudo else 'contorno NAO coberto')
cont = json.load(open(RAIZ/'build/rouche_L_rig/contagem.json'))
eB = Fr(cont['eps_B'])*(1 + Fr(1, 10**9))
tmax = max(D[3]/(1 - D[2]*D[3]) for D in usados)
pior = max(usados, key=lambda D: D[3]/(1 - D[2]*D[3]))
fontes = sorted({D[4] for D in usados})
print(f'max t(s) nos {len(usados)} discos da cadeia <= {float(tmax):.3f} ({pior[4]}, p = {float(pior[0])}+{float(pior[1])}i); '
      f't eps_B <= {float(tmax*eB):.4e} < 1: {tmax*eB < 1}')
print('arquivos usados:', fontes)
json.dump(dict(coberto=ok_tudo, n_discos=len(discos), usados=len(usados), arestas=relato, tmax=float(tmax),
               t_epsB=float(tmax*eB), arquivos=fontes), open(RAIZ/'.codex-runs/2026-10-02-revisao-B/r7_cobertura_resultado.json', 'w'), indent=1)
