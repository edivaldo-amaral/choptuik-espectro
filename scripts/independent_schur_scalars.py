#!/usr/bin/env python3
"""Recomputo INDEPENDENTE (auditoria 16/09) dos escalares de Schur.

A2: subsistema S+T do setor A (GUARDED_INFINITE_COUPLING.md).
A4: ponte sharp F+T e orcamento de Neumann da faixa U (SHARP_INTERMEDIATE_SHELL.md).

Entradas usadas e sua procedencia (todas recomputadas nesta auditoria):
  beta_eta <= 5191/500          independent_component_beta.py (recalculado 10.38108...)
  beta_sharp <= 462237/100000   independent_component_beta.py (recalculado 4.6223694...)
  ||Q F|| = 1/mu_RefA exato      independent_free_inverse.py
  cauda livre (rederivada)       independent_free_inverse.tail_bound
  parte RefA filtrada = 0        independent_component_beta.py (176x944, 160x960)
  kappa_crown                    relatorio OU independent_crown_tiles.py
  kappa_sharp                    relatorio OU independent_sharp_disk.py --mode right
  ||O^-1 D'|| <= 2(1+r)/((1-r) mu)  RT (37b) com kappa2 generico
Nenhum ponto flutuante nas decisoes; floats so para exibicao/busca inicial.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import json
from pathlib import Path

from independent_free_inverse import tail_bound, MU_REFA

EPS_MU = Q(43, 2**40)+Q(1, 2**277)
EPS_OM = Q(1, 2**25)+Q(1, 2**277)
KAPPA_CROWN_REPORT = Q(669414234798400, 18478486838353)
KAPPA_SHARP_REPORT = None  # lido do relatorio


def K2prime(kappa2: Q, mu: Q = MU_REFA) -> Q:
    r = 1/kappa2
    return 2*(1+r)/((1-r)*mu)


def crown_schur(kappa: Q, outer=(176, 944), beta=Q(5191, 500), radius=Q(5, 4)):
    e = K2prime(Q(5, 4))*EPS_MU
    d = (beta+radius)*tail_bound(*outer, Q(5, 4))+e
    c = 6*Q(4096, 243)*EPS_OM/MU_REFA
    schur = d+kappa*c*d
    return dict(outer=list(outer), d=float(d), c=float(c), kappa=float(kappa),
                schur=str(schur), schur_decimal=float(schur), schur_lt_0_997676=schur < Q(997676, 10**6),
                invertible=schur < 1, kappa_max_admissible=float((1-d)/(c*d)) if d < 1 else None)


def degree(delta: Q, ell: Q, r: Q, tol: Q) -> int:
    lo, hi = 0, 1
    f = lambda N: ell*r*delta**(N+1)/(1-delta) <= tol
    while not f(hi):
        lo, hi = hi, 2*hi
    while hi-lo > 1:
        mid = (lo+hi)//2
        if f(mid):
            hi = mid
        else:
            lo = mid
    return hi


def up(x: Q, den=10**6) -> Q:
    y = x*den
    return Q(-(-y.numerator//y.denominator), den)


def sharp_bridge(kappa: Q, outer, beta=Q(462237, 100000), radius=Q(4, 5)):
    k2 = Q(9, 8)
    e = K2prime(k2)*EPS_MU
    d = (beta+radius)*tail_bound(*outer, k2)+e
    qglob = max(1/MU_REFA, tail_bound(12, 36, k2))
    c = 3*EPS_OM*(1/MU_REFA)
    delta = d+c*kappa*d
    dF = (beta+radius)*tail_bound(12, 36, k2)+e
    ell = d+beta*qglob*kappa*d
    r = dF+c*kappa*dF
    dl, el, rr = up(delta), up(ell), up(r)
    M, N = outer
    dofs = (2*((M+1)//2)-1)*(N+N//2)-(2*((12+1)//2)-1)*(36+36//2)
    return dict(outer=list(outer), d=float(d), delta=str(delta), delta_decimal=float(delta),
                ell=float(ell), r=float(r), U_dofs=dofs,
                degree_1e12=degree(dl, el, rr, Q(1, 10**12)),
                degree_1e3=degree(dl, el, rr, Q(1, 10**3)),
                degree_1e1=degree(dl, el, rr, Q(1, 10)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--crown-json', type=Path, help='saida de independent_crown_tiles.py')
    parser.add_argument('--sharp-json', type=Path, help='saida de independent_sharp_disk.py --mode right')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    out = {}
    out['crown_schur_report_kappa'] = crown_schur(KAPPA_CROWN_REPORT)
    if args.crown_json:
        tiles = json.loads(args.crown_json.read_text())['tiles']
        kappa = max(Q(t['inverse_bound_exact']) for t in tiles)
        out['crown_schur_independent_kappa'] = crown_schur(kappa)
    rep = json.loads((root/'build/spectrum/sharp-guarded-bridge-160x960-preconditioned.json').read_text())
    ks = Q(rep['preconditioned_finite_inverse_bound'])
    out['sharp_report_kappa'] = float(ks)
    for outer in ((160, 960), (148, 840)):
        out[f'sharp_{outer[0]}x{outer[1]}_report_kappa'] = sharp_bridge(ks, outer)
    if args.sharp_json:
        mine = Q(json.loads(args.sharp_json.read_text())['preconditioned_inverse_bound'])
        out['sharp_independent_kappa'] = float(mine)
        for outer in ((160, 960), (148, 840)):
            out[f'sharp_{outer[0]}x{outer[1]}_independent_kappa'] = sharp_bridge(mine, outer)
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
