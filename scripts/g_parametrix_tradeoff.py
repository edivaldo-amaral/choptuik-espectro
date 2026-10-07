#!/usr/bin/env python3
"""Curva de trade-off para caixas G'=(M',N') (docs/G_PARAMETRIX_FEASIBILITY.md, item 2).

DIAGNOSTICO em float (rotulado) + cotas EXATAS da Fase 1 onde baratas:
  * d_float(G'): max amostrado de ||(B_ref+5/4)Qe|| (setor A, norma eta) sobre colunas
    fora de G' (radial n=N', m amostrado; Fourier m=M', n amostrado), x1.5 (pre-registro);
  * d_R(N') exato (columnwise_exterior.radial_bound, betas da Fase 1);
  * d_F(M') exato por escalonamento EXATO d_F(M)=(d_F(107)-e_mu)*107/M+e_mu (a cota e Phi/b+e_mu).
"""
from __future__ import annotations

import json
from fractions import Fraction as Q

import columnwise_exterior as CE
from columnwise_criteria import float_actual_A

MS = (24, 48, 64, 107)
NS = (128, 192, 256, 384)


def main():
    rep = json.loads((CE.ROOT/'build/spectrum/columnwise-exterior.json').read_text())
    eA = rep['A']
    per = {int(d): dict(k=(2 if CE.SYSTEMS['A']['kinds'][int(d)] == 'K' else 1), beta_R=Q(v['beta_R']),
                        beta_inf=Q(v['beta_inf']), bF=None) for d, v in eA['per_component'].items()}
    cert = dict(sysname='A', e_B=Q(eA['e_B']), e_mu=Q(eA['e_mu']), per=per)
    dF107, emu = Q(eA['d_fourier_at_MG']), Q(eA['e_mu'])
    s = 1.25
    radial = {}
    for N in NS:
        best = (0, None)
        for d in range(4):
            n = N+1 if (d == 0 and N % 2 == 0) else N
            for m in (0, 1, 2, 3, 4, 6, 8, 12, 13, 16, 18, 20, 24, 28, 32, 40, 48, 56, 64, 80, 100):
                if m % 2 != CE.SYSTEMS['A']['mpar'][d]:
                    continue
                for p in ((0,) if m == 0 else (0, 1)):
                    v = float_actual_A(d, p, m, n, s)
                    if v > best[0]:
                        best = (v, (d, p, m, n))
        radial[N] = best
        print('radial', N, best, flush=True)
    fourier = {}
    for M in MS:
        best = (0, None)
        for d in range(4):
            m = M if M % 2 == CE.SYSTEMS['A']['mpar'][d] else M+1
            for n in (0, 1, 2, 3, 5, 9, 17, 25, 33, 49, 65, 97, 127):
                if d == 0 and n % 2 == 0:
                    continue
                for p in (0, 1):
                    v = float_actual_A(d, p, m, n, s)
                    if v > best[0]:
                        best = (v, (d, p, m, n))
        fourier[M] = best
        print('fourier', M, best, flush=True)
    rows = []
    for M in MS:
        for N in NS:
            dfl = max(radial[N][0], fourier[M][0])
            dR = CE.radial_bound(cert, N)[0]
            dF = (dF107-emu)*Q(107, M)+emu
            rows.append(dict(M=M, N=N, dofs_A=CE.dofs('A', M, N), dofs_sharp=CE.dofs('sharp', M, N),
                             d_float_max=dfl, d_float_x15=1.5*dfl,
                             d_R_exact=float(dR), d_F_exact=float(dF), d_phase1=float(max(dR, dF)),
                             rt_style_margin_float=1-1.5*dfl))
    out = dict(note='d_float e DIAGNOSTICO (x1.5 pre-registrado); d_R, d_F exatos (Fase 1)',
               s=s, radial=radial, fourier=fourier, rows=rows)
    (CE.ROOT/'build/spectrum/g-parametrix-tradeoff.json').write_text(json.dumps(out, indent=2)+'\n')
    for r in rows:
        print(r)


if __name__ == '__main__':
    main()
