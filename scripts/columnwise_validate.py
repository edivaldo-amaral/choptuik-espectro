#!/usr/bin/env python3
"""Criterios 1 e 2 de docs/COLUMNWISE_EXTERIOR.md: colunas exteriores amostradas.

Para cada coluna e FORA da caixa certificada calcula um ENCLAUSURAMENTO SUPERIOR
EXATO (inteiros/racionais, arredondamento diadico com raio propagado) de
    U(e) = ||B_ref Q e|| + (s_max+e_B)||Qe|| + e_mu  >=  ||(H_true(s)-I)e||,
e confere U(e) <= d*(G). O float (criterio 2) e so medida de folga.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as Q
import json
from pathlib import Path

import columnwise_exterior as CE

P_BITS = 96


def rnd(x: Q) -> int:
    return (2*x.numerator*(1 << P_BITS)+x.denominator)//(2*x.denominator)


def q_column(sysname, d, m, n, depth):
    """c_j para j em [max(0,n-depth), n], como (mid_re, mid_im) em 2^-P e raio (Fraction)."""
    kind = CE.SYSTEMS[sysname]['kinds'][d]
    k = 1 if kind == 'J' else 2
    mu, b = CE.MU, Q(m, 2)
    D = lambda j: (mu*(j+1), b)
    grid = Q(1, 1 << P_BITS)

    def inv(z):
        den = z[0]*z[0]+z[1]*z[1]
        return (z[0]/den, -z[1]/den)

    def mul(z, w):
        return (z[0]*w[0]-z[1]*w[1], z[0]*w[1]+z[1]*w[0])
    cols = {}
    cn = inv(D(n))
    cols[n] = ((rnd(cn[0]), rnd(cn[1])), grid)   # raio em MODULO
    if n-k < 0:
        return cols, k
    ck = mul(inv(D(n)), inv(D(n-k)))
    ck = (-2*mu*n*ck[0], -2*mu*n*ck[1])
    mid = (rnd(ck[0]), rnd(ck[1]))
    rad = 2*grid
    j = n-k
    while True:
        cols[j] = (mid, rad)
        if j-k < 0 or n-(j-k) > depth:
            break
        jj = j-k
        f = ((mu*jj, -b), D(jj)) if k == 1 else ((mu*(jj+1), -b), D(jj))
        num, den = f
        fi = mul(num, inv(den))
        exact = mul((-fi[0], -fi[1]), (Q(mid[0], 1 << P_BITS), Q(mid[1], 1 << P_BITS)))
        mid = (rnd(exact[0]), rnd(exact[1]))
        rad = rad+grid   # |f|<=1; raio em modulo, arredondamento <= grid/sqrt2
        j = jj
    return cols, k


def upper(args):
    sysname, d, p, m, n, depth = args
    cfg = CE.SYSTEMS[sysname]
    r = CE.r_of(sysname)
    cols, k = q_column(sysname, d, m, n, depth)
    tabs = CE.tables(sysname)
    (p1, q1), (p2, q2) = cfg['k1'], cfg['k2']
    out = [dict() for _ in range(cfg['ncomp'])]
    qnorm = Q(0)
    radpart = Q(0)
    field_abs = {}
    for j, ((xr, xi), rad) in cols.items():
        wrel = Q(q2, p2)**(n-j)*(Q(1, 2) if j == 0 and n > 0 else 1)
        qnorm += wrel*(Q(abs(xr)+abs(xi), 1 << P_BITS)+2*rad)
        # imagens de entrada: parte p, conj em -m
        base = (xr, xi) if p == 0 else (-xi, xr)
        imgs = {}
        for sg in (1, -1):
            v = base if sg == 1 else (base[0], -base[1])
            for tg in (1, -1):
                imgs[sg*m, tg*j] = v
        key = CE.key_for(sysname, d, j)
        for i, field in enumerate(tabs[key]):
            o = out[i]
            for (m0, n0), (vx, vy) in imgs.items():
                for (m1, n1), (a, bb) in field.items():
                    kk = (m0+m1, n0+n1)
                    x, y = o.get(kk, (0, 0))
                    o[kk] = (x+a*vx-bb*vy, y+a*vy+bb*vx)
            if rad:
                if (key, i) not in field_abs:
                    field_abs[key, i] = sum((abs(a)+abs(bb))*Q(p1, q1)**abs(m1)*Q(p2, q2)**abs(n1)
                                            for (m1, n1), (a, bb) in field.items())/(1 << 72)
                radpart += wrel*2*rad*Q(cfg['eta'][i], cfg['eta'][d])*field_abs[key, i]
        for comp, mult in CE.r_terms(sysname, d, j):
            o = out[comp]
            for sg in ({1, -1} if m else {1}):
                vx, vy = imgs[sg*m, j]
                L = 0
                while j-1-2*L >= 0:
                    amp = CE.MU_INT72*2*mult*(-1)**L
                    for nn in {j-1-2*L, -(j-1-2*L)}:
                        x, y = o.get((sg*m, nn), (0, 0))
                        o[sg*m, nn] = (x+amp*vx, y+amp*vy)
                    L += 1
            if rad:
                radpart += wrel*2*rad*Q(cfg['eta'][comp], cfg['eta'][d])*abs(mult)*2*CE.MU*r/(1-r*r)
    n_img = len({(sg*m, tg*n) for sg in (1, -1) for tg in (1, -1)})
    norm_in = n_img*cfg['eta'][d]*Q(p1, q1)**m*Q(p2, q2)**n*(1 << P_BITS)
    bq = CE.weighted_norm(sysname, out, m, n)/norm_in
    # colunas descartadas (profundidade): cota por grau com beta_glob
    jmin = min(cols)
    remainder = Q(0)
    if jmin-k >= 0:
        nu = Q(1) if k == 1 else Q(n, n-1)
        tail_c = CE.SQ2*2*nu/(CE.MU*(n+1))
        remainder = (cfg['beta_glob']+cfg['s_max'])*tail_c*r**(n-jmin+k)/(1-r**k)
    cert = CE.certify.__globals__  # noqa (so para acesso a constantes)
    if sysname == 'A':
        e_B = CE.EPS_MU*Q(43390, 2061)+6*Q(4096, 243)*CE.EPS_OM
    else:
        e_B = CE.EPS_MU*Q(144, 17)+3*CE.EPS_OM
    e_mu = CE.EPS_MU*2*(1+r)/((1-r)*CE.MU)
    U = bq+radpart+(cfg['s_max']+e_B)*qnorm+e_mu+remainder
    return dict(system=sysname, d=d, p=p, m=m, n=n, depth=depth, U=str(U), U_decimal=float(U),
                BQ_decimal=float(bq), q_decimal=float(qnorm), radius_part_decimal=float(radpart),
                remainder_decimal=float(remainder))


def schur_margins(rep):
    """Criterio 4: Schur GUARDED (S+T) e ponte sharp (F+T) com d*(G) no lugar de d(G).

    c_A = ||T H F|| para G >= 47x119 (parte RefA filtrada = 0): 6(4096/243) eps_om ||Q F||,
    ||Q F||=1/mu_ref (A2.1, A2.3). kappa_crown: 31.4266... (recertificado) e 36.2267 (relatorio).
    Sharp: c = 3 eps_om ||Q F|| (||Q F||=1/mu_ref), kappa_sharp = 11.46011129 (recertificado).
    """
    out = {}
    MG, NG = rep['A']['M_G_certified'], rep['A']['N_G_certified']
    if MG and NG:
        MG, NG = max(MG, 47), max(NG, 119)
        dA = max(Q(rep['A']['d_radial_at_NG']), Q(rep['A']['d_fourier_at_MG']))
        cA = 6*Q(4096, 243)*CE.EPS_OM/CE.MU
        for label, kappa in (('kappa_independente', Q(314266490, 10**7)),
                             ('kappa_relatorio', Q(669414234798400, 18478486838353))):
            s = dA+kappa*cA*dA
            out[f'guarded_{label}'] = dict(box=[MG, NG], d_star=float(dA), schur=float(s),
                                          schur_exact=str(s), invertible=s < 1,
                                          U_dofs=CE.dofs('A', MG, NG)-CE.dofs('A', 6, 18))
    MS, NS = rep['sharp']['M_G_certified'], rep['sharp']['N_G_certified']
    if MS and NS:
        MS, NS = max(MS, 53), max(NS, 137)   # parte RefA filtrada de T<-F e zero
        dS = max(Q(rep['sharp']['d_radial_at_NG']), Q(rep['sharp']['d_fourier_at_MG']))
        cS = 3*CE.EPS_OM/CE.MU
        kappa = Q(1146011130, 10**8)
        s = dS+cS*kappa*dS
        out['sharp_bridge'] = dict(box=[MS, NS], d_star=float(dS), delta=float(s), delta_exact=str(s),
                                   invertible=s < 1,
                                   U_dofs=CE.dofs('sharp', MS, NS)-CE.dofs('sharp', 12, 36))
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report', type=Path, default=CE.ROOT/'build/spectrum/columnwise-exterior.json')
    parser.add_argument('--workers', type=int, default=2)
    parser.add_argument('--depth', type=int, default=160)
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    rep = json.loads(args.report.read_text())
    jobs = []
    for sysname in ('A', 'sharp'):
        NG, MG = rep[sysname]['N_G_certified'], rep[sysname]['M_G_certified']
        cfg = CE.SYSTEMS[sysname]
        for d in range(cfg['ncomp']):
            mp = cfg['mpar'][d]
            K = cfg['kinds'][d] == 'K'
            ms = sorted({m for m in (0, 1, 5, 12, 13, 40, 41, 42, 140, 141) if m % 2 == mp})
            for m in ms[:4]:
                n = NG + (1 if K and NG % 2 == 0 else 0)
                for p in ((0,) if m == 0 else (0, 1)):
                    jobs.append((sysname, d, p, m, n, args.depth))
            mF = MG + (1 if MG % 2 != mp else 0)
            for n in (1, 3, 50, 101):
                if K and n % 2 == 0:
                    continue
                jobs.append((sysname, d, 1, mF, n, args.depth))
    results = []
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        for res in ex.map(upper, jobs):
            dstar = max(Q(rep[res['system']].get('d_radial_at_NG', '1')),
                        Q(rep[res['system']].get('d_fourier_at_MG', '1')))
            res['d_star'] = float(dstar)
            res['U_le_dstar'] = Q(res['U']) <= dstar
            print(json.dumps(res), flush=True)
            results.append(res)
    summary = dict(columns=len(results), all_pass=all(r['U_le_dstar'] for r in results),
                   max_U_over_dstar=max(r['U_decimal']/r['d_star'] for r in results),
                   schur=schur_margins(rep), results=results)
    if args.out:
        args.out.write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps({k: v for k, v in summary.items() if k != 'results'}, indent=2))


if __name__ == '__main__':
    main()
