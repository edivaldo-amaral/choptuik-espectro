#!/usr/bin/env python3
"""Montagem final do exterior (Marco 8 de docs/COLUMNWISE_EXTERIOR.md), cortes arbitrarios.

Le TODOS os checkpoints de build/spectrum/columnwise-phase2/, confere em racionais que as
regioes cobrem o complemento de G=(M_G,N_G) sem buracos, e monta
    d*(G) = max( max das regioes, d_R(N_far), d_F(M') ),
com d_R e d_F EXATOS da Fase 1 (d_F e exatamente Phi/b+e_mu, logo escala como 1/M).
Recalcula as margens do Schur GUARDED e da ponte sharp. Criterio pre-registrado: d*<=7/10.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import json
from pathlib import Path

import columnwise_exterior as CE

CK = CE.ROOT/'build/spectrum/columnwise-phase2'


def load_rows(sysname):
    """{n: lista de (m_lo, M, U)} a partir de todos os checkpoints do sistema."""
    rows = {}
    for path in sorted(CK.glob(f'{sysname}-*.jsonl')):
        for line in path.read_text().splitlines():
            if not line.strip():
                continue
            rec = json.loads(line)
            if rec['system'] != sysname:
                continue
            m_lo = rec.get('m_lo', 0)
            rows.setdefault(rec['n'], []).append((m_lo, rec['M'], Q(rec['max_U']), path.name))
    return rows


def covered(intervals, lo, hi):
    """intervalos [a,b) cobrem [lo,hi)?"""
    cur = lo
    for a, b in sorted(intervals):
        if a > cur:
            return False
        cur = max(cur, b)
        if cur >= hi:
            return True
    return cur >= hi


def assemble(sysname, M_G, N_G, M_far, N_far, rep):
    rows = load_rows(sysname)
    holes, worst, where = [], Q(0), None
    for n in range(0, N_far):
        need = (M_G, M_far) if n < N_G else (0, M_far)
        entries = rows.get(n, [])
        if not covered([(a, b) for a, b, _, _ in entries], *need):
            holes.append(n)
        for a, b, u, src in entries:
            if u > worst:
                worst, where = u, dict(n=n, m_lo=a, M=b, file=src)
    per = {int(d): dict(k=(2 if CE.SYSTEMS[sysname]['kinds'][int(d)] == 'K' else 1),
                        beta_R=Q(v['beta_R']), beta_inf=Q(v['beta_inf']), bF=None)
           for d, v in rep[sysname]['per_component'].items()}
    cert = dict(sysname=sysname, e_B=Q(rep[sysname]['e_B']), e_mu=Q(rep[sysname]['e_mu']), per=per)
    dR = CE.radial_bound(cert, N_far)[0]
    emu = Q(rep[sysname]['e_mu'])
    dF107 = Q(rep[sysname]['d_fourier_at_MG'])
    assert rep[sysname]['M_G_certified'] == 107
    dF = (dF107-emu)*Q(107, M_far)+emu
    dstar = max(worst, dR, dF) if not holes else None
    return dict(G=[M_G, N_G], N_far=N_far, M_far=M_far, holes=holes[:20], holes_count=len(holes),
                regions_max_U=str(worst), regions_max_U_decimal=float(worst), regions_argmax=where,
                d_R_far=str(dR), d_R_far_decimal=float(dR), d_F=str(dF), d_F_decimal=float(dF),
                d_star=str(dstar) if dstar else None,
                d_star_decimal=float(dstar) if dstar else None,
                meets_0_7=bool(dstar is not None and dstar <= Q(7, 10)),
                rows_seen=len(rows))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--M-G', type=int, default=107)
    parser.add_argument('--N-G', type=int, default=384)
    parser.add_argument('--M-far', type=int, default=160)
    parser.add_argument('--N-far', type=int, default=900)
    parser.add_argument('--out', type=Path, default=CE.ROOT/'build/spectrum/columnwise-exterior-phase3.json')
    args = parser.parse_args()
    repfile = CE.ROOT/'build/spectrum/columnwise-exterior.json'
    rep = json.loads(repfile.read_text())
    out = dict(marco=8, criterion='d*<=7/10 em A e sharp, cobertura sem buracos',
               conditional_on_published_RT_ball=True, physical_certificate=False)
    for sysname in ('A', 'sharp'):
        out[sysname] = assemble(sysname, args.M_G, args.N_G, args.M_far, args.N_far, rep)
    ok = all(out[s]['meets_0_7'] for s in ('A', 'sharp'))
    out['meta_atingida'] = ok
    if out['A']['d_star']:
        dA = Q(out['A']['d_star'])
        cA = 6*Q(4096, 243)*CE.EPS_OM/CE.MU
        for label, kappa in (('kappa_crown_independente', Q(314266490, 10**7)),
                             ('kappa_crown_relatorio', Q(669414234798400, 18478486838353))):
            s = dA+kappa*cA*dA
            out[f'guarded_schur_{label}'] = dict(schur=str(s), schur_decimal=float(s), invertible=s < 1,
                                                 margin_decimal=float(1-s))
    if out['sharp']['d_star']:
        dS = Q(out['sharp']['d_star'])
        s = dS+3*CE.EPS_OM/CE.MU*Q(1146011130, 10**8)*dS
        out['sharp_bridge'] = dict(delta=str(s), delta_decimal=float(s), invertible=s < 1,
                                   margin_decimal=float(1-s))
    args.out.write_text(json.dumps(out, indent=2)+'\n')
    rep['phase3_exterior'] = out
    repfile.write_text(json.dumps(rep, indent=2)+'\n')
    print(json.dumps({k: (v if not isinstance(v, dict) else
                          {kk: vv for kk, vv in v.items() if 'decimal' in kk or kk in
                           ('holes_count', 'meets_0_7', 'regions_argmax', 'invertible', 'rows_seen')})
                      for k, v in out.items()}, indent=2))


if __name__ == '__main__':
    main()
