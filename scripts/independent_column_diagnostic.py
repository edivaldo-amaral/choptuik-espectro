#!/usr/bin/env python3
"""DIAGNOSTICO EM PONTO FLUTUANTE (auditoria 16/09) -- NADA AQUI E PROVA.

Pergunta: o majorante uniforme d(G)=(beta+|s|) t(G) usado em GUARDED_INFINITE_COUPLING
e SHARP_* e quantas vezes maior que as normas de COLUNA reais de (B+s)Q no exterior?
Se a razao for grande, a estrategia de RT (EXTERIOR_ESTIMATE: normas de coluna
calculadas uma a uma numa faixa, bound analitico so muito longe) reduz as caixas;
se for ~1, o gargalo e estrutural e mais computacao nao resolve.

Norma: pesos RT fortes (65/64,5/4) x eta=(916,4096,243,1233), soma por coluna.
Operador: B = mu S Gamma1 + 2 S Xi Gamma2(RefA,.), Q = O_mu^-1 (formulas do TeX de RT,
eq. (22)-(23)); coeficientes de (23) vem de independent_component_beta.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import json

import numpy as np
from scipy.signal import fftconvolve

from independent_component_beta import parse_refa, selected_coefficients
from independent_free_inverse import exact_column, tail_bound, MU_REFA

ETA = np.array([916., 4096., 243., 1233.])
K1, K2 = 65/64, 5/4
MF, NF = 40, 100


def coef_arrays():
    fields, scale = parse_refa()
    tab = selected_coefficients(fields)
    out = {}
    for key, (c, vec) in tab.items():
        arrs = []
        for f in vec:
            a = np.zeros((2*MF+1, 2*NF+1), dtype=complex)
            for (m, n), (re, im) in f.items():
                z = complex(float(2*c*re*scale), float(2*c*im*scale))
                for nn in {n, -n}:
                    a[MF+m, NF+nn] = z
                    a[MF-m, NF-nn] = z.conjugate()
            arrs.append(a)
        out[key] = arrs
    return out


def column_ratio(coefs, d, m, n, s, mu=MU_REFA):
    """||(B+s) Q e||_eta / ||e||_eta para e real no endereco (d, m, n)."""
    kind = 'K' if d == 0 else 'J'
    col = [complex(float(z.re), float(z.im)) for z in exact_column(kind, m, n, mu)]
    Mo, No = m+MF, n+NF
    h = np.zeros((2*m+1, 2*n+1), dtype=complex)
    for j, z in enumerate(col):
        for jj in {j, -j}:
            h[m+m, n+jj] = z
            h[m-m, n-jj] = z.conjugate()
    par = (-1.0)**np.arange(-n, n+1)
    even, odd = h*(par == 1), h*(par == -1)
    outs = [np.zeros((2*Mo+1, 2*No+1), dtype=complex) for _ in range(4)]
    keys = [(('w1', None), h)] if d == 0 else [((f'w{d+1}', 'e'), even), ((f'w{d+1}', 'o'), odd)]
    for key, piece in keys:
        for i, cf in enumerate(coefs[key]):
            if np.any(cf) and np.any(piece):
                outs[i] += fftconvolve(cf, piece, mode='full')
    # restringe a m>=0, n>=0
    outs = [o[Mo:, No:].copy() for o in outs]
    # mu S Gamma1: (0, R(h1+(1-P)h2), -R h1, R(1-P)h4) na linha m
    hn = np.zeros(No+1, dtype=complex)
    hn[:n+1] = col
    Pn = (-1.0)**np.arange(No+1)

    def R(v):
        r = np.zeros_like(v)
        for j in range(len(v)):
            idx = np.arange(j+1, len(v), 2)
            r[j] = 2*np.sum(((-1.0)**np.arange(len(idx)))*v[idx])
        return r
    muf = float(mu)
    if d == 0:
        outs[1][m, :] += muf*R(hn)
        outs[2][m, :] += -muf*R(hn)
    elif d == 1:
        outs[1][m, :] += muf*R(hn-Pn*hn)
    elif d == 3:
        outs[3][m, :] += muf*R(hn-Pn*hn)
    outs[d][m, :] += s*hn
    mw = (2-(np.arange(Mo+1) == 0))[:, None]*K1**(np.arange(Mo+1)-m)[:, None]
    nw = (2-(np.arange(No+1) == 0))[None, :]*K2**(np.arange(No+1)-n)[None, :]
    wout = mw*nw
    win = (2-(m == 0))*(2-(n == 0))
    total = sum(ETA[i]*np.sum(wout*(np.abs(o.real)+np.abs(o.imag))) for i, o in enumerate(outs))
    return float(total/(ETA[d]*win))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--s', type=complex, default=1+0j)
    parser.add_argument('--out')
    args = parser.parse_args()
    coefs = coef_arrays()
    beta = 5191/500+abs(args.s)
    rows = []
    for n in (18, 36, 72, 144, 288, 576):
        for m in (0, 1, 2, 6, 12, 24, 48):
            for d in range(4):
                if m % 2 != (1 if d == 3 else 0):
                    continue
                nn = n+1 if (d == 0 and n % 2 == 0) else n
                r = column_ratio(coefs, d, m, nn, args.s)
                rows.append(dict(d=d, m=m, n=nn, column=r))
    by_n = {}
    for row in rows:
        by_n.setdefault(row['n'] - (row['n'] % 2 if row['d'] == 0 and row['n'] % 2 else 0), []).append(row)
    summary = []
    for n in (18, 36, 72, 144, 288, 576):
        cols = [r for r in rows if r['n'] in (n, n+1)]
        worst = max(cols, key=lambda r: r['column'])
        # termo radial do majorante uniforme: (beta+|s|) * 81/(n+1)
        radial_majorant = beta*float(Q(3, 2)*9/(Q(1, 6)*(n+1)))
        summary.append(dict(n=n, max_sampled_column=worst['column'], argmax=(worst['d'], worst['m']),
                            radial_uniform_majorant=radial_majorant,
                            ratio=radial_majorant/worst['column']))
    result = dict(s=str(args.s), note='DIAGNOSTICO float, amostra de colunas; nao e prova',
                  summary=summary, rows=rows)
    text = json.dumps(result, indent=2)
    if args.out:
        open(args.out, 'w').write(text+'\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
