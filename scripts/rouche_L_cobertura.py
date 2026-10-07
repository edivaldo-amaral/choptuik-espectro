#!/usr/bin/env python3
"""Rouché de L: os discos certificados (saídas de rouche_L_disco com ok) cobrem o contorno inteiro
d([0, 3] x [-1/4, 1/4])? Aritmética racional; a meia-corda sqrt(r^2 - h^2) de cada disco numa aresta é
cotada por baixo. Imprime as lacunas (trechos de aresta sem disco), se houver."""
from __future__ import annotations

import argparse
from fractions import Fraction as Fr
import json
import math
from pathlib import Path

def arestas(re0, re1, im0, im1):
    """(nome, ponto inicial, direcao unitaria, comprimento) do bordo de [re0, re1] x [im0, im1] (racionais)."""
    return [(f'esquerda Re s = {re0}', (re0, im0), (0, 1), im1 - im0),
            (f'superior Im s = {im1}', (re0, im1), (1, 0), re1 - re0),
            (f'direita Re s = {re1}', (re1, im0), (0, 1), im1 - im0),
            (f'inferior Im s = {im0}', (re0, im0), (1, 0), re1 - re0)]


ARESTAS = arestas(Fr(0), Fr(3), Fr(-1, 4), Fr(1, 4))   # faixa A


def sqrt_inf(x: Fr) -> Fr:
    """Cota inferior racional de sqrt(x), x >= 0."""
    if x <= 0:
        return Fr(0)
    s = Fr(math.sqrt(float(x)))*(1 - Fr(1, 10**12))
    while s*s > x:
        s *= 1 - Fr(1, 10**9)
    return s


def intervalos(discos, a, u, L):
    out = []
    for (px, py), r in discos:
        # t da projeção e distância à reta (u é (1, 0) ou (0, 1))
        if u == (1, 0):
            t, h = px - a[0], abs(py - a[1])
        else:
            t, h = py - a[1], abs(px - a[0])
        if h > r:
            continue
        w = sqrt_inf(r*r - h*h)
        lo, hi = max(Fr(0), t - w), min(L, t + w)
        if lo < hi or (lo == hi and w > 0):
            out.append((lo, hi))
    return sorted(out)


def lacunas(iv, L):
    gaps = []; x = Fr(0)
    for lo, hi in iv:
        if lo > x:
            gaps.append((x, lo))
        x = max(x, hi)
    if x < L:
        gaps.append((x, L))
    return gaps


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('discos', type=Path, nargs='+', help='JSONs de rouche_L_disco --saida')
    ap.add_argument('--aceita-sem-certF', dest='aceita_sem_certF', action='store_true', help='diagnostico: aceita discos sem certF')
    ap.add_argument('--regiao', default='0,3,-1/4,1/4', help='re0,re1,im0,im1 (racionais; faixa B: 0,3,1/4,3/4)')
    a = ap.parse_args()
    re0, re1, im0, im1 = (Fr(x) for x in a.regiao.split(','))
    discos = []
    for f in a.discos:
        dd = json.load(open(f))
        if not dd.get('certF') and not a.aceita_sem_certF:
            print(f'  IGNORADO (sem certF, nao rigoroso): {f}')
            continue
        for x in dd['pontos']:
            if x['ok'] and x['r'] > 0:
                discos.append(((Fr(x['p'][0]), Fr(x['p'][1])), Fr(x['r'])))
    print(f'{len(discos)} discos certificados')
    tudo = True
    for nome, a0, u, L in arestas(re0, re1, im0, im1):
        g = lacunas(intervalos(discos, a0, u, L), L)
        if g:
            tudo = False
            for lo, hi in g:
                z0 = (a0[0] + lo*u[0], a0[1] + lo*u[1]); z1 = (a0[0] + hi*u[0], a0[1] + hi*u[1])
                print(f'  LACUNA na aresta {nome}: de {float(z0[0]):.6f}{float(z0[1]):+.6f}i a '
                      f'{float(z1[0]):.6f}{float(z1[1]):+.6f}i (comprimento {float(hi - lo):.2e})')
        else:
            print(f'  aresta {nome}: coberta')
    print('CONTORNO COBERTO' if tudo else 'contorno NAO coberto')


if __name__ == '__main__':
    main()
