#!/usr/bin/env python3
"""Cobertura da regiao finita de S3a por tiles de Perron (tile_perron), em tres escalas.

Regiao: 0 <= Re s <= RE_MAX, 0 <= Im s <= IM_MAX (K e real: injetividade no espaco com
sinal para Im s >= 0 da a do espaco simetrico, invariante pela conjugacao), menos o
disco D de S3b em torno do autovalor 0,401.

  * quadrado GRANDE de lado h2 -> centro s2 das colunas de T2 e far, raio r2 = h2/sqrt2
    (uma parte T2 por quadrado grande: o item caro);
  * quadrados MEDIOS -> centro s1 das colunas de T1, raio r1 = h1/sqrt2. Comeca com
    k1 x k1 e subdivide enquanto a base (rZ = 0) nao fechar;
  * quadrados PEQUENOS -> o tile |s - sZ| <= rZ = hZ/sqrt2. O lado inicial vem do rZ
    admissivel no centro do medio; cada tile que falha e subdividido (ate --prof).

O disco circunscrito de um subquadrado cabe no do quadrado-pai (|s1 - s2| + r1 <= r2 com
igualdade no canto), entao as condicoes de centro do tile valem por construcao. Um
quadrado inteiro dentro de D e pulado; um que corta D e subdividido.

Uma invocacao processa UM quadrado grande (--grande i,j), para rodar em paralelo.
Saida: JSON por quadrado grande com cada tile (centro, raio, raio de Perron, resultado).
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys
import time

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import tile_perron as tp

RE_MAX, IM_MAX = 1.765, 0.25
DISCO = (0.401, 0.125)           # (centro, raio) do disco de S3b: docs/S3B_DISCO.md
# raios inflados: o disco circunscrito exato poe os vertices da malha SOBRE os circulos, e
# h/sqrt2 arredondado para baixo os deixaria descobertos. Com o mesmo fator em r1 e r2 as
# contencoes |s1 - s2| + r1 <= r2 e |sZ - s1| + rZ <= r1 continuam valendo.
FOLGA = 1 + 1e-6
raio = lambda h: h/math.sqrt(2)*FOLGA


def dentro_do_disco(c, h, disco):
    """Quadrado de centro c e lado h inteiro dentro do disco?"""
    cD, rD = disco
    return max(abs(c + (dx + 1j*dy)*h/2 - cD) for dx in (-1, 1) for dy in (-1, 1)) <= rD


def corta_regiao(c, h):
    return (c.real + h/2 > 0 and c.real - h/2 < RE_MAX and c.imag + h/2 > 0 and c.imag - h/2 < IM_MAX)


def filhos(c, h):
    return [(c + complex(dx, dy)*h/4, h/2) for dx in (-1, 1) for dy in (-1, 1)]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--grande', required=True, help='indices i,j do quadrado grande')
    ap.add_argument('--h2', type=float, default=0.25)
    ap.add_argument('--k1', type=int, default=2, help='divisoes iniciais do quadrado grande')
    ap.add_argument('--prof1', type=int, default=3, help='subdivisoes extras do quadrado medio')
    ap.add_argument('--prof', type=int, default=3, help='subdivisoes extras do tile pequeno')
    ap.add_argument('--alvo', type=float, default=0.97, help='raio de Perron visado ao escolher rZ')
    ap.add_argument('--Z', default='20x80'); ap.add_argument('--Gp', default='40x160')
    ap.add_argument('--F', default='200x400'); ap.add_argument('--D', type=int, default=12)
    ap.add_argument('--disco', default=f'{DISCO[0]},{DISCO[1]}')
    ap.add_argument('--rigoroso', action='store_true')
    ap.add_argument('--dir', type=Path, default=Path('build/cobertura'))
    a = ap.parse_args()
    tp.RIG = a.rigoroso
    cx = lambda t: tuple(map(int, t.split('x')))
    Z, Gp, F = cx(a.Z), cx(a.Gp), cx(a.F)
    disco = tuple(map(float, a.disco.split(',')))
    i, j = map(int, a.grande.split(','))
    h2 = a.h2
    s2 = complex((i + 0.5)*h2, (j + 0.5)*h2); r2 = raio(h2)
    a.dir.mkdir(parents=True, exist_ok=True)
    saida = a.dir/f'grande-{i}-{j}.json'
    registro = dict(s2=[s2.real, s2.imag], r2=r2, h2=h2, Z=Z, Gp=Gp, F=F, D=a.D,
                    rigoroso=a.rigoroso, disco=disco, tiles=[], medios=[])
    salvar = lambda: saida.write_text(json.dumps(registro, indent=1) + '\n')
    if not corta_regiao(s2, h2) or dentro_do_disco(s2, h2, disco):
        registro['pulado'] = True; salvar(); return
    ctx = tp.Ctx(Z, Gp, F, a.D)
    silencio = lambda *x, **y: None
    cache2 = a.dir/f'T2-{i}-{j}.npz'
    t0 = time.time()
    agora = lambda: f'[{time.time()-t0:.0f}s]'
    h1 = h2/a.k1
    medios = [(s2 + complex((u + 0.5)*h1 - h2/2, (v + 0.5)*h1 - h2/2), h1, 0)
              for u in range(a.k1) for v in range(a.k1)]
    n_cache = 0
    while medios:
        s1, h1, p1 = medios.pop()
        if not corta_regiao(s1, h1) or dentro_do_disco(s1, h1, disco):
            continue
        r1 = raio(h1)
        cache1 = a.dir/f'T1-{i}-{j}-{n_cache}.npz'; n_cache += 1
        roda = lambda sZ, rZ: tp.tile(sZ, rZ, Z, Gp, F, D=a.D, s1=s1, r1=r1, s2=s2, r2=r2,
                                      cache1=cache1, cache2=cache2, log=silencio, ctx=ctx, alvo=a.alvo)
        q = roda(s1, 0.0)
        med = dict(s1=[s1.real, s1.imag], r1=r1, prof=p1, rZ_adm=q['r_adm'], falhas=0)
        registro['medios'].append(med)
        if q['r_adm'] <= 0.0:
            base = q['perron']
            print(f'medio {s1:.4f} h1 {h1:.4f}: base {base:.4f} nao fecha  {agora()}', flush=True)
            if p1 < a.prof1:
                medios += [(c, h, p1 + 1) for c, h in filhos(s1, h1)]
                med['subdividido'] = True
            else:
                med['falhas'] = 1
            salvar(); continue
        kZ = max(1, math.ceil(h1/(math.sqrt(2)*0.9*q['r_adm'])))
        hZ0 = h1/kZ
        print(f'medio {s1:.4f} h1 {h1:.4f}: rZ admissivel {q["r_adm"]:.4f} -> {kZ}x{kZ} tiles  {agora()}',
              flush=True)
        fila = [(s1 + complex((u + 0.5)*hZ0 - h1/2, (v + 0.5)*hZ0 - h1/2), hZ0, 0)
                for u in range(kZ) for v in range(kZ)]
        while fila:
            sZ, hZ, prof = fila.pop()
            if not corta_regiao(sZ, hZ) or dentro_do_disco(sZ, hZ, disco):
                continue
            rZ = raio(hZ)
            q = roda(sZ, rZ)
            ok = q['certificado'] if a.rigoroso else q['perron'] < 1.0
            registro['tiles'].append(dict(sZ=[sZ.real, sZ.imag], rZ=rZ, prof=prof, perron=q['perron'],
                                          theta=q.get('perron_theta'), ok=bool(ok), s1=[s1.real, s1.imag],
                                          VQ_ZZ=q['Z']['VQ_ZZ'], V=q['Z']['V'], r1=r1, s2=[s2.real, s2.imag],
                                          r2=r2, v=q.get('perron_v'),
                                          N0=q['N0'], NZ=q['NZ'], N1=q['N1'], N2=q['N2']))
            if not ok:
                print(f'  tile {sZ:.4f} r {rZ:.4f}: Perron {q["perron"]:.4f} falhou  {agora()}', flush=True)
                if prof < a.prof:
                    fila += [(c, h, prof + 1) for c, h in filhos(sZ, hZ)]
                else:
                    med['falhas'] += 1
            salvar()
        salvar()
    n_ok = sum(t['ok'] for t in registro['tiles'])
    registro['resumo'] = dict(tiles=len(registro['tiles']), certificados=n_ok,
                              falhas_finais=sum(m['falhas'] for m in registro['medios']),
                              medios=len(registro['medios']), tempo=time.time() - t0)
    salvar()
    print(f'quadrado ({i},{j}) em {s2}: {registro["resumo"]}')


if __name__ == '__main__':
    main()
