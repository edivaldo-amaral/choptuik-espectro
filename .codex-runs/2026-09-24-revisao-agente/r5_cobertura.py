#!/usr/bin/env python3
"""R5 -- cobertura da regiao por pontos, em racionais exatos.

Regiao: 0 <= Re s <= 1,765, 0 <= Im s <= 1/4, fora do disco ABERTO |s - 0,401| < 1/8.
Discos: build/cobertura_lab/grande-*.json cuja reverificacao passa. Reverifico com
check_cover.reverifica E com uma reimplementacao propria (Perron e contencoes em Fraction).

Pontos (todos racionais exatos):
  * 10^5 uniformes na regiao (denominador 2^40);
  * 2000 em cada lado do retangulo e 2000 pontos racionais EXATAMENTE sobre a circunferencia de D
    (parametrizacao racional do circulo);
  * 10^4 perto dos vertices das malhas (tiles, quadrados medios e grandes): o vertice (float ->
    racional exato) e deslocamentos de ate 1e-7 hZ.
"""
from __future__ import annotations

from fractions import Fraction as Fr
import glob
import json
import math
from pathlib import Path
import random
import sys
import time

import numpy as np

AQUI = Path(__file__).resolve().parent
RAIZ = AQUI.parents[1]
sys.path.insert(0, str(RAIZ/'scripts'))
import check_cover as cc

RE, IM = Fr(1765, 1000), Fr(1, 4)
CD, RD = Fr(401, 1000), Fr(1, 8)


def meu_reverifica(t):
    ARRED = 1 + Fr(1, 10**9)
    rZ, r1, r2 = Fr(t['rZ']), Fr(t['r1']), Fr(t['r2'])
    N = [[ARRED*(Fr(t['N0'][i][j]) + rZ*Fr(t['NZ'][i][j]) + r1*Fr(t['N1'][i][j]) + r2*Fr(t['N2'][i][j]))
          for j in range(4)] for i in range(4)]
    v = [Fr(x) for x in t['v']]
    if min(v) <= 0 or any(N[i][j] < 0 for i in range(4) for j in range(4)):
        return False, None
    th = max(sum(N[i][j]*v[j] for j in range(4))/v[i] for i in range(4))
    sZ = (Fr(t['sZ'][0]), Fr(t['sZ'][1]))
    for c, r in ((t['s1'], r1), (t['s2'], r2)):
        d2 = (sZ[0] - Fr(c[0]))**2 + (sZ[1] - Fr(c[1]))**2
        if r < rZ or (r - rZ)**2 < d2:
            return False, th
    # centros com Re >= 0 (a forma fechada pede Re s_l >= 0)
    if min(sZ[0], Fr(t['s1'][0]), Fr(t['s2'][0])) < 0:
        return False, th
    return th < 1, th


def carrega():
    discos = []; todos = []; thetas = []; divergem = 0; recusados = 0
    for f in sorted(glob.glob(str(RAIZ/'build/cobertura_lab/grande-*.json'))):
        r = json.loads(Path(f).read_text())
        assert r['rigoroso'] is True and tuple(r['Z']) == (20, 80) and tuple(r['Gp']) == (40, 160) \
            and tuple(r['F']) == (200, 400)
        for t in r['tiles']:
            todos.append(t)
            ok1, _ = cc.reverifica(t)
            ok2, th = meu_reverifica(t)
            if ok1 != ok2:
                divergem += 1
            if t.get('ok') and ok1 and ok2:
                discos.append((Fr(t['sZ'][0]), Fr(t['sZ'][1]), Fr(t['rZ']), t))
                thetas.append(float(th))
            elif t.get('ok'):
                recusados += 1
    return discos, todos, thetas, divergem, recusados


def fora_de_D(x, y):
    return (x - CD)**2 + y*y >= RD*RD


def na_regiao(x, y):
    return 0 <= x <= RE and 0 <= y <= IM and fora_de_D(x, y)


class Indice:
    def __init__(self, discos):
        self.d = discos
        self.cx = np.array([float(a) for a, b, r, t in discos]); self.cy = np.array([float(b) for a, b, r, t in discos])
        self.r = np.array([float(r) for a, b, r, t in discos])

    def coberto(self, x, y):
        xf, yf = float(x), float(y)
        cand = np.where((self.cx - xf)**2 + (self.cy - yf)**2 <= (self.r*(1 + 1e-9) + 1e-12)**2)[0]
        for k in cand:
            a, b, r, _ = self.d[k]
            if (x - a)**2 + (y - b)**2 <= r*r:
                return True
        return False


def main():
    t0 = time.time()
    discos, todos, thetas, divergem, recusados = carrega()
    print(f'tiles gravados: {len(todos)}; marcados ok: {sum(bool(t.get("ok")) for t in todos)}; '
          f'discos aceitos (check_cover.reverifica E a minha): {len(discos)}; recusados: {recusados}; '
          f'divergencias entre as duas reverificacoes: {divergem}')
    print(f'theta: max {max(thetas):.6f}, min {min(thetas):.6f}; rZ min {min(float(d[2]) for d in discos):.5f}')
    idx = Indice(discos)
    rng = random.Random(24092026)
    DEN = 2**40
    faltam = []
    # 1) uniformes
    n = 0; testados = 0
    while n < 100_000:
        x = Fr(rng.randrange(0, int(RE*DEN) + 1), DEN); y = Fr(rng.randrange(0, int(IM*DEN) + 1), DEN)
        if not na_regiao(x, y):
            continue
        n += 1
        if not idx.coberto(x, y):
            faltam.append(('uniforme', x, y))
    print(f'1) 10^5 uniformes: descobertos {len(faltam)}  [{time.time()-t0:.0f}s]')
    # 2) bordas do retangulo e circunferencia de D (pontos racionais exatos sobre o circulo)
    nb = 0
    for k in range(2000):
        u = Fr(rng.randrange(0, DEN + 1), DEN)
        for (x, y) in ((Fr(0), u*IM), (RE, u*IM), (u*RE, Fr(0)), (u*RE, IM)):
            if na_regiao(x, y):
                nb += 1
                if not idx.coberto(x, y):
                    faltam.append(('borda', x, y))
        # circulo: ((1 - t^2)/(1 + t^2), 2t/(1 + t^2)), t racional; y >= 0 <=> t >= 0
        tt = Fr(rng.randrange(0, 4*DEN), DEN)
        x = CD + RD*(1 - tt*tt)/(1 + tt*tt); y = RD*2*tt/(1 + tt*tt)
        assert (x - CD)**2 + y*y == RD*RD
        if na_regiao(x, y):
            nb += 1
            if not idx.coberto(x, y):
                faltam.append(('circulo D', x, y))
    # os dois pontos do circulo no eixo real e o topo
    for (x, y) in ((CD - RD, Fr(0)), (CD + RD, Fr(0)), (CD, RD)):
        nb += 1
        if not idx.coberto(x, y):
            faltam.append(('circulo D (especial)', x, y))
    print(f'2) {nb} pontos nas bordas e sobre a circunferencia de D: descobertos (acumulado) {len(faltam)}'
          f'  [{time.time()-t0:.0f}s]')
    # 3) perto dos vertices das malhas
    vert = []
    for (a, b, r, t) in discos:
        h = float(r)/(1 + 1e-6)*math.sqrt(2)
        for dx in (-0.5, 0.5):
            for dy in (-0.5, 0.5):
                vert.append((float(a) + dx*h, float(b) + dy*h, h))
        for c, rr in ((t['s1'], t['r1']), (t['s2'], t['r2'])):
            h = rr/(1 + 1e-6)*math.sqrt(2)
            for dx in (-0.5, 0.5):
                for dy in (-0.5, 0.5):
                    vert.append((c[0] + dx*h, c[1] + dy*h, h))
    vert = list({(round(x, 15), round(y, 15)): (x, y, h) for x, y, h in vert}.values())
    nv = 0
    alvo = 10_000
    i = 0
    while nv < alvo:
        x0, y0, h = vert[i % len(vert)]; i += 1
        if i <= len(vert):
            x, y = Fr(x0), Fr(y0)                  # o proprio vertice, exato
        else:
            x = Fr(x0) + Fr(rng.randrange(-10**6, 10**6 + 1), 10**13)*Fr(h)
            y = Fr(y0) + Fr(rng.randrange(-10**6, 10**6 + 1), 10**13)*Fr(h)
        if not na_regiao(x, y):
            continue
        nv += 1
        if not idx.coberto(x, y):
            faltam.append(('vertice', x, y))
    print(f'3) {nv} pontos em/perto de {len(vert)} vertices: descobertos (acumulado) {len(faltam)}'
          f'  [{time.time()-t0:.0f}s]')
    for f_ in faltam[:10]:
        print('   DESCOBERTO:', f_[0], float(f_[1]), float(f_[2]))
    print('RESUMO R5:', 'todos os pontos cobertos' if not faltam else f'{len(faltam)} pontos DESCOBERTOS')


if __name__ == '__main__':
    main()
