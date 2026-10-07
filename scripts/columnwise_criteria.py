#!/usr/bin/env python3
"""Criterios 1 e 2 de docs/COLUMNWISE_EXTERIOR.md (amostra de colunas exteriores a G).

Criterio 1 (EXATO): para cada coluna amostrada fora de G=(107,384), um majorante
racional U(e) >= ||(H_true(s)-I)e|| (uniforme em |s|<=s_max, bola RT) e conferido
U(e) <= d*(G). Motores: nucleo C exato (n-D>=101) e implementacao Python direta
com 4 imagens (colunas de n pequeno na direcao de Fourier).
Criterio 2 (DIAGNOSTICO float, setor A): ||(B_ref+s_max)Qe|| em ponto flutuante,
para medir a folga d*/float e U/float. Nao e prova.
Saida: build/spectrum/columnwise-criteria.json (lida por test_columnwise_exterior.py).
"""
from __future__ import annotations

import json
from fractions import Fraction as Q

import columnwise_exterior as CE
import columnwise_phase2 as P2
import columnwise_validate as V

OUT = CE.ROOT/'build/spectrum/columnwise-criteria.json'


def sample():
    A_band = [(0, 1, 18), (0, 0, 0), (1, 0, 0), (1, 1, 6), (2, 1, 40), (3, 1, 1), (3, 0, 41), (3, 0, 105)]
    cols = []
    for n in (384, 385, 500, 501, 605):
        for d, p, m in A_band:
            if d == 0 and n % 2 == 0:
                continue
            cols.append(('A', d, p, m, n, 'c'))
    for n in (700, 1001):
        for d, p, m in ((0, 1, 20), (2, 0, 0), (3, 0, 57)):
            if d == 0 and n % 2 == 0:
                continue
            cols.append(('A', d, p, m, n, 'c'))
    cols += [('A', 2, 0, 108, 384, 'c'), ('A', 3, 0, 201, 701, 'c'), ('A', 0, 1, 110, 385, 'c')]
    cols += [('A', 1, 1, 108, 1, 'py'), ('A', 3, 0, 107, 3, 'py'), ('A', 0, 1, 150, 5, 'py')]
    for n in (385, 605):
        for d, p, m in ((0, 1, 4), (1, 0, 104), (1, 0, 0), (0, 0, 50)):
            if d == 0 and n % 2 == 0:
                continue
            cols.append(('sharp', d, p, m, n, 'c'))
    cols += [('sharp', 1, 1, 7, 700, 'c'), ('sharp', 0, 0, 30, 701, 'c'),
             ('sharp', 1, 0, 108, 2, 'py'), ('sharp', 0, 1, 107, 1, 'py')]
    return cols


def float_actual_A(d, p, m, n, s):
    """DIAGNOSTICO: ||(B_ref+s)Qe||/||e|| em float (coluna de Q por recorrencia float)."""
    import independent_column_diagnostic as D
    mu = float(CE.MU)
    k = 2 if d == 0 else 1
    b = m/2

    class Z:
        def __init__(self, z):
            self.re, self.im = z.real, z.imag

    def col(kind, m_, n_, mu_):
        Dj = lambda j: mu*(j+1)+1j*b
        c = [0j]*(n_+1)
        c[n_] = 1/Dj(n_)
        if n_-k >= 0:
            c[n_-k] = -2*mu*n_/(Dj(n_)*Dj(n_-k))
            j = n_-k
            while j-k >= 0:
                jj = j-k
                f = (mu*jj-1j*b)/Dj(jj) if k == 1 else (mu*(jj+1)-1j*b)/Dj(jj)
                c[jj] = -f*c[j]
                j = jj
        if p == 1:
            c = [1j*z for z in c]
        return [Z(z) for z in c]
    old = D.exact_column
    D.exact_column = col
    try:
        return D.column_ratio(_COEFS(), d, m, n, s)
    finally:
        D.exact_column = old


_C = []


def _COEFS():
    if not _C:
        import independent_column_diagnostic as D
        _C.append(D.coef_arrays())
    return _C[0]


def main():
    rep = json.loads((CE.ROOT/'build/spectrum/columnwise-exterior.json').read_text())['phase2']
    dstar = {s: Q(rep[s]['d_star']) for s in ('A', 'sharp')}
    results = []
    todo = sample()
    by_sys = {}
    for c in todo:
        if c[5] == 'c':
            by_sys.setdefault(c[0], []).append(c)
    for sysname, cs in by_sys.items():
        depth = 60 if sysname == 'A' else 100
        jobs = [(sysname, d, p, m, n, depth) for (_, d, p, m, n, _) in cs]
        for c, res in zip(cs, P2.run_row_c(sysname, jobs, 2)):
            results.append(dict(system=sysname, d=c[1], p=c[2], m=c[3], n=c[4], engine='C', depth=depth,
                                U=res['U'], U_decimal=res['U_decimal']))
    for (sysname, d, p, m, n, eng) in todo:
        if eng == 'py':
            res = V.upper((sysname, d, p, m, n, n+2))
            results.append(dict(system=sysname, d=d, p=p, m=m, n=n, engine='python-4-imagens', depth=n+2,
                                U=res['U'], U_decimal=res['U_decimal']))
    for r in results:
        r['d_star'] = str(dstar[r['system']])
        r['U_le_d_star'] = Q(r['U']) <= dstar[r['system']]
        if r['system'] == 'A':
            fl = float_actual_A(r['d'], r['p'], r['m'], r['n'], float(CE.SYSTEMS['A']['s_max']))
            r['float_actual_diagnostic'] = fl
            r['d_star_over_float'] = float(dstar['A'])/fl
            r['U_over_float'] = r['U_decimal']/fl
        print(json.dumps({k: r[k] for k in ('system', 'd', 'p', 'm', 'n', 'U_decimal', 'U_le_d_star')}
                         | ({'float': r['float_actual_diagnostic']} if 'float_actual_diagnostic' in r else {})),
              flush=True)
    A = [r for r in results if 'float_actual_diagnostic' in r]
    summary = dict(columns=len(results), criterion1_all=all(r['U_le_d_star'] for r in results),
                   has_m0=any(r['m'] == 0 for r in results), has_odd_m=any(r['m'] % 2 for r in results),
                   has_K=any(r['d'] == 0 for r in results), has_J=any(r['d'] > 0 for r in results),
                   criterion2_float_columns=len(A),
                   criterion2_max_d_star_over_float=max(r['d_star_over_float'] for r in A),
                   criterion2_max_U_over_float=max(r['U_over_float'] for r in A),
                   criterion2_count_d_star_le_2x=sum(r['d_star_over_float'] <= 2 for r in A),
                   note='criterio 1 exato (majorante U>=norma); criterio 2 so diagnostico float',
                   results=results)
    OUT.write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps({k: v for k, v in summary.items() if k != 'results'}, indent=2))


if __name__ == '__main__':
    main()
