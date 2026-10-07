#!/usr/bin/env python3
"""Consolida a Fase 2 (docs/COLUMNWISE_EXTERIOR.md) e aplica a regra pre-registrada.

d*(G) = max( max_faixa U(e), d_R(N_far) [Fase 1], d_F(M) [Fase 1] ),  G=(M, N_lo).
Tudo em racionais exatos, lidos dos checkpoints e do relatorio da Fase 1.
Tambem refaz os escalares de Schur (GUARDED e ponte sharp) com d*(G).
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import json
from pathlib import Path

import columnwise_exterior as CE

ROOT = CE.ROOT
CK = ROOT/'build/spectrum/columnwise-phase2'


def phase1_cert(rep, sysname):
    per = {int(d): dict(k=(2 if CE.SYSTEMS[sysname]['kinds'][int(d)] == 'K' else 1),
                        beta_R=Q(v['beta_R']), beta_inf=Q(v['beta_inf']), bF=None)
           for d, v in rep[sysname]['per_component'].items()}
    return dict(sysname=sysname, e_B=Q(rep[sysname]['e_B']), e_mu=Q(rep[sysname]['e_mu']), per=per)


def band(sysname, M, depth, n_lo, n_hi):
    path = CK/f'{sysname}-M{M}-D{depth}.jsonl'
    rows = {}
    for line in path.read_text().splitlines():
        if line.strip():
            rec = json.loads(line)
            rows[rec['n']] = rec
    missing = [n for n in range(n_lo, n_hi) if n not in rows]
    worst, where = Q(0), None
    for n in range(n_lo, n_hi):
        if n in rows and Q(rows[n]['max_U']) > worst:
            worst, where = Q(rows[n]['max_U']), dict(n=n, **rows[n]['argmax'])
    columns = sum(rows[n]['columns'] for n in range(n_lo, n_hi) if n in rows)
    return dict(missing_rows=missing, max_U=worst, argmax=where, columns=columns,
                profile={n: rows[n]['max_U_decimal'] for n in range(n_lo, n_hi, 20) if n in rows})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=ROOT/'build/spectrum/columnwise-exterior-phase2.json')
    args = parser.parse_args()
    rep1 = json.loads((ROOT/'build/spectrum/columnwise-exterior.json').read_text())
    plan = {'A': dict(M=107, depth=60, n_lo=384, n_far=606), 'sharp': dict(M=107, depth=100, n_lo=384, n_far=606)}
    report = dict(phase=2, rule='N_G<=384 em A e sharp => rota de blocos continua; senao inviavel neste hardware',
                  conditional_on_published_RT_ball=True, physical_certificate=False)
    for sysname, pl in plan.items():
        cert = phase1_cert(rep1, sysname)
        dR = CE.radial_bound(cert, pl['n_far'])[0]
        dF = Q(rep1[sysname]['d_fourier_at_MG']) if rep1[sysname]['M_G_certified'] == pl['M'] else None
        b = band(sysname, pl['M'], pl['depth'], pl['n_lo'], pl['n_far'])
        complete = not b['missing_rows'] and dF is not None
        dstar = max(b['max_U'], dR, dF) if complete else None
        report[sysname] = dict(G=[pl['M'], pl['n_lo']], N_far=pl['n_far'], depth=pl['depth'],
                               band_columns=b['columns'], band_missing_rows=len(b['missing_rows']),
                               band_max_U=str(b['max_U']), band_max_U_decimal=float(b['max_U']),
                               band_argmax=b['argmax'], band_profile_every_20=b['profile'],
                               d_R_far=str(dR), d_R_far_decimal=float(dR),
                               d_F=str(dF) if dF is not None else None,
                               d_F_decimal=float(dF) if dF is not None else None,
                               d_star=str(dstar) if dstar is not None else None,
                               d_star_decimal=float(dstar) if dstar is not None else None,
                               N_G_certified=pl['n_lo'] if (dstar is not None and dstar < 1) else None,
                               box_dofs=CE.dofs(sysname, pl['M'], pl['n_lo']))
    NA, NS = report['A']['N_G_certified'], report['sharp']['N_G_certified']
    report['decision'] = ('rota de blocos continua (N_G<=384 em A e sharp)' if NA and NS and NA <= 384 and NS <= 384
                          else 'incompleto ou N_G>384')
    # Schur com d*(G)
    if report['A']['d_star']:
        dA = Q(report['A']['d_star'])
        cA = 6*Q(4096, 243)*CE.EPS_OM/CE.MU            # ||T H F||, G>=47x119
        for label, kappa in (('kappa_crown_independente', Q(314266490, 10**7)),
                             ('kappa_crown_relatorio', Q(669414234798400, 18478486838353))):
            s = dA+kappa*cA*dA
            report[f'guarded_schur_{label}'] = dict(schur=str(s), schur_decimal=float(s), invertible=s < 1,
                                                   margin_decimal=float(1-s),
                                                   U_dofs=CE.dofs('A', 107, 384)-CE.dofs('A', 6, 18))
    if report['sharp']['d_star']:
        dS = Q(report['sharp']['d_star'])
        cS = 3*CE.EPS_OM/CE.MU                          # ||T H F||, G>=53x137
        kappa = Q(1146011130, 10**8)
        s = dS+cS*kappa*dS
        report['sharp_bridge'] = dict(delta=str(s), delta_decimal=float(s), invertible=s < 1,
                                      margin_decimal=float(1-s),
                                      U_dofs=CE.dofs('sharp', 107, 384)-CE.dofs('sharp', 12, 36))
    args.out.write_text(json.dumps(report, indent=2, default=str)+'\n')
    print(json.dumps({k: (v if not isinstance(v, dict) else {kk: vv for kk, vv in v.items() if 'decimal' in kk or kk in ('N_G_certified', 'band_columns', 'band_missing_rows', 'G', 'invertible', 'U_dofs', 'band_argmax')})
                      for k, v in report.items()}, indent=2, default=str))


if __name__ == '__main__':
    main()
