#!/usr/bin/env python3
"""Verificacao independente da cobertura de S3a: a uniao dos discos CERTIFICADOS
(|s - sZ| <= rZ com Perron verificado) cobre a regiao

    R = [0, RE_MAX] x [0, IM_MAX]  menos o disco aberto D de S3b.

Le so os JSON de saida (grande-*.json), nao a logica do planejador. Quadtree propria em
aritmetica RACIONAL exata. Uma celula esta resolvida se:
  * esta inteira no disco ABERTO D (os 4 cantos estritamente dentro; D e convexo), ou
  * esta inteira num disco certificado (os 4 cantos nele; o disco e convexo).
Senao, subdivide ate --prof; o que sobrar e reportado. A circunferencia de D faz parte
da regiao a cobrir. Celulas que passam da borda da regiao sao exigidas inteiras (mais
forte que o necessario).
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Fr
import json
from pathlib import Path


def reverifica(t, caixas=None):
    """Refaz, em racionais exatos, o que certifica o tile: N v <= theta v com theta < 1 para
    N = N0 + rZ NZ + r1 N1 + r2 N2 e o v gravado; e as contencoes dos discos dos centros.
    Nao confia na marca 'ok' da rodada."""
    if not all(c in t for c in ('N0', 'NZ', 'N1', 'N2', 'v', 's1', 'r1', 's2', 'r2')) or t['v'] is None:
        return False, 'sem dados'
    rZ, r1, r2 = Fr(t['rZ']), Fr(t['r1']), Fr(t['r2'])
    # 1 + 1e-9: as entradas gravadas sao somas de cotas em ponto flutuante (arredondamento)
    ARRED = 1 + Fr(1, 10**9)
    N = [[ARRED*(Fr(t['N0'][i][j]) + rZ*Fr(t['NZ'][i][j]) + r1*Fr(t['N1'][i][j]) + r2*Fr(t['N2'][i][j]))
          for j in range(4)] for i in range(4)]
    if caixas is not None:
        # far -> Z e far -> T1 pelo caminho de mu R_x (desce de qualquer n), que os tiles
        # calculados antes de 24/09 (tarde) nao incluiam; somar de novo nos novos e conservador
        from tile_perron import descida_rx
        (MZ, NZ), (MP, NP), (Mc, Nc) = caixas
        N[0][3] += Fr(t['V'])*Fr(descida_rx(NZ, Nc))*ARRED
        N[1][3] += Fr(descida_rx(NP, Nc))*ARRED
        # mu verdadeiro de RT (|mu - mu_RefA| <= 3,9e-11): coluna do nivel l ganha eps_mu(s_l),
        # vezes ||V|| na linha de Z (revisao independente, 24/09)
        from tile_perron import eps_mu
        cent = [complex(*t['sZ']), complex(*t['s1']), complex(*t['s2']), complex(*t['s2'])]
        for j in range(4):
            e = Fr(eps_mu(cent[j]))*ARRED
            N[0][j] += Fr(t['V'])*e
            for i in (1, 2, 3):
                N[i][j] += e
    v = [Fr(x) for x in t['v']]
    if min(v) <= 0:
        return False, 'v nao positivo'
    theta = max(sum(N[i][j]*v[j] for j in range(4))/v[i] for i in range(4))
    if theta >= 1:
        return False, f'theta {float(theta):.4f}'
    # contencoes: |sZ - s1| + rZ <= r1 e |sZ - s2| + rZ <= r2, sem raiz: (r - rZ)^2 >= |.|^2
    sZ = (Fr(t['sZ'][0]), Fr(t['sZ'][1]))
    for c, r in ((t['s1'], r1), (t['s2'], r2)):
        d2 = (sZ[0] - Fr(c[0]))**2 + (sZ[1] - Fr(c[1]))**2
        if r < rZ or (r - rZ)**2 < d2:
            return False, 'contencao'
    return True, float(theta)


def carrega(dirs):
    discos = []; recusados = {}; thetas = []
    for d in dirs:
        for f in sorted(Path(d).glob('grande-*.json')):
            r = json.loads(f.read_text())
            for t in r.get('tiles', []):
                if not t.get('ok'):
                    continue
                ok, info = reverifica(t, (r['Z'], r['Gp'], r['F']))
                if ok:
                    discos.append((Fr(t['sZ'][0]), Fr(t['sZ'][1]), Fr(t['rZ'])))
                    thetas.append(info)
                else:
                    recusados[info] = recusados.get(info, 0) + 1
    if recusados:
        print(f'tiles marcados ok mas RECUSADOS na reverificacao: {recusados}')
    if thetas:
        print(f'pior theta reverificado: {max(thetas):.7f}')
    return discos


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('dirs', nargs='+')
    ap.add_argument('--re-min', default='0'); ap.add_argument('--re-max', default='1.765')
    ap.add_argument('--im-max', default='0.25')
    ap.add_argument('--disco', default='0.401,0.125')
    ap.add_argument('--h0', default='0.05', help='lado das celulas iniciais')
    ap.add_argument('--prof', type=int, default=30,
                    help='a folga dos raios e 1e-6 r: celulas sobre vertices precisam descer a ~r 1e-6')
    ap.add_argument('--max-celulas', type=int, default=2_000_000)
    a = ap.parse_args()
    RE0, RE, IM, h0 = Fr(a.re_min), Fr(a.re_max), Fr(a.im_max), Fr(a.h0)
    cD, rD = (Fr(x) for x in a.disco.split(','))
    discos = carrega(a.dirs)
    print(f'{len(discos)} discos certificados')
    dentro = lambda x, y, cx, cy, r: (x - cx)**2 + (y - cy)**2 <= r*r

    def cantos(x0, y0, h):
        return [(x0, y0), (x0 + h, y0), (x0, y0 + h), (x0 + h, y0 + h)]

    def em_D(x0, y0, h):          # celula inteira no disco ABERTO D (cantos estritamente dentro)
        return all((x - cD)**2 + y*y < rD*rD for x, y in cantos(x0, y0, h))

    def coberta(x0, y0, h):
        cx_, cy_ = x0 + h/2, y0 + h/2
        for (dx, dy, r) in discos:
            # filtro barato: centro da celula perto do disco
            if (cx_ - dx)**2 + (cy_ - dy)**2 > (r + h)**2:
                continue
            if all(dentro(x, y, dx, dy, r) for x, y in cantos(x0, y0, h)):
                return True
        return False

    fila = []
    nx = -(-(RE - RE0) // h0); ny = -(-IM // h0)
    for i in range(int(nx)):
        for j in range(int(ny)):
            fila.append((RE0 + i*h0, j*h0, h0, 0))
    faltam = []
    resolvidas = 0
    if not discos:
        print('nenhum disco certificado: nada a verificar'); print('COBERTURA INCOMPLETA'); return
    visitadas = 0
    while fila:
        visitadas += 1
        if visitadas > a.max_celulas:
            print(f'limite de {a.max_celulas} celulas: a cobertura tem buracos grandes'); break
        x0, y0, h, p = fila.pop()
        # recorta a regiao: celulas que comecam fora de R sao descartadas; a regiao e fechada
        if x0 > RE or y0 > IM:
            continue
        if em_D(x0, y0, h) or coberta(x0, y0, h):
            resolvidas += 1
            continue
        if p < a.prof:
            h2 = h/2
            fila += [(x0, y0, h2, p + 1), (x0 + h2, y0, h2, p + 1),
                     (x0, y0 + h2, h2, p + 1), (x0 + h2, y0 + h2, h2, p + 1)]
        else:
            faltam.append((float(x0), float(y0), float(h)))
    print(f'celulas resolvidas: {resolvidas}; sem cobertura: {len(faltam)}')
    for c in faltam[:20]:
        print(f'  celula [{c[0]:.5f}, {c[0]+c[2]:.5f}] x [{c[1]:.5f}, {c[1]+c[2]:.5f}]')
    print('COBERTURA COMPLETA' if not faltam and not fila else 'COBERTURA INCOMPLETA')


if __name__ == '__main__':
    main()
