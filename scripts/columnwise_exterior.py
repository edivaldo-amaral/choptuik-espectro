#!/usr/bin/env python3
"""Majorante exterior por colunas d*(G) (Fase 1, semianalitica, EXATA).

Pre-registro e prova: docs/COLUMNWISE_EXTERIOR.md. Nenhum ponto flutuante
entra em decisoes; floats aparecem so como campos *_decimal para leitura.

Sistemas:
  A      L_true(s) Q_ref, pesos (65/64,5/4) x eta=(916,4096,243,1233), |s|<=5/4.
  sharp  K_true(s) Qsharp_ref, pesos (129/128,9/8), sem eta, |s|<=4/5.

Para cada coluna unitaria e_{d,p,m,n}:
  ||(H_true-I)e|| <= sum_j (w_j/w_n)||c_j||_C betabar(d,m,j) + (|s|+e_B) q + e_mu,
c = Q e (RT (22): c_n=1/D_n, c_{n-k}=-2mu n/(D_n D_{n-k}), c_j=-f_j c_{j+k}),
beta(d,p,m,j)=||B_ref e_{d,p,m,j}||/||e|| EXATO (RefA suporte 41x101 => valores
independentes de j>=101 e de m>=41, a menos de uma cauda explicita do termo R).
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path

from independent_component_beta import parse_refa, selected_coefficients, sharp_coefficients

ROOT = Path(__file__).resolve().parents[1]
REFA_SHA = '5776e6038f994e06b4b82fc638d797154f630aeee37b251058512d3cb4cd450b'
MU = Q(722873400, 2**32)
MU_INT72 = 722873400 << 40            # mu_ref em unidades de 2^-72
EPS_MU = Q(43, 2**40)+Q(1, 2**277)    # RT secao 5 e (60)
EPS_OM = Q(1, 2**25)+Q(1, 2**277)     # RT (39b) e (60)
SQ2 = Q(99, 70)                        # >= sqrt(2)
G0 = Q(1207107, 1000000)               # >= (1+sqrt 2)/2
NU2 = Q(1155, 1000)                    # >= 2/sqrt(3)
SUPPORT_M, SUPPORT_N = 41, 101         # RefA: 0<=m<=40, 0<=n<=100

SYSTEMS = {
    'A': dict(k1=(65, 64), k2=(5, 4), eta=(916, 4096, 243, 1233), ncomp=4, s_max=Q(5, 4),
              beta_glob=Q(5191, 500), kinds=('K', 'J', 'J', 'J'), mpar=(0, 0, 0, 1)),
    'sharp': dict(k1=(129, 128), k2=(9, 8), eta=(1, 1), ncomp=2, s_max=Q(4, 5),
                  beta_glob=Q(462237, 100000), kinds=('K', 'J'), mpar=(0, 0)),
}


def r_of(sysname):
    a, b = SYSTEMS[sysname]['k2']
    return Q(b, a)


# ----------------------------------------------------------------- campos
_TABLES = {}


def tables(sysname):
    """dict key -> lista (por saida) de dicts Z^2 (m1,n1)->(re,im) em unidades 2^-72."""
    if sysname in _TABLES:
        return _TABLES[sysname]
    path = ROOT/'.cache/rt-1203.3766v1/sourcecode/RefA.dat'
    if hashlib.sha256(path.read_bytes()).hexdigest() != REFA_SHA:
        raise SystemExit('RefA diferente do auditado')
    fields, scale = parse_refa(path)
    assert scale == Q(1, 2**70)
    if sysname == 'A':
        raw, factor = selected_coefficients(fields), 2
    else:
        raw, factor = sharp_coefficients(fields), 1
    out = {}
    for key, (c, vec) in raw.items():
        F = 4*factor*c
        assert F.denominator == 1
        F = int(F)
        comps = []
        for f in vec:
            z2 = {}
            for (m, n), (a, b) in (f or {}).items():
                a, b = a*F, b*F
                for nn in {n, -n}:
                    z2[m, nn] = (a, b)
                    z2[-m, -nn] = (a, -b)
            comps.append({k: v for k, v in z2.items() if v != (0, 0)})
        out[key] = comps
    _TABLES[sysname] = out
    return out


def key_for(sysname, d, j):
    if d == 0:
        return ('w1', None)
    return (f'w{d+1}', 'e' if j % 2 == 0 else 'o')


def r_terms(sysname, d, j):
    """Saidas do termo mu*R (Gamma1): lista (comp, mult)."""
    if sysname == 'A':
        if d == 0:
            return [(1, 1), (2, -1)]
        if d in (1, 3) and j % 2 == 1:
            return [(d, 2)]
        return []
    return [(1, 1)] if d == 0 else []


def input_images(m, j, p):
    vals = {}
    for sg in (1, -1):
        for tg in (1, -1):
            vals[sg*m, tg*j] = (1, 0) if p == 0 else ((0, 1) if sg == 1 else (0, -1))
    return vals


def weighted_norm(sysname, out, center_m, center_n):
    """sum_i eta_i sum (|re|+|im|) k1^|M| k2^|N|, EXATO; devolve Fraction nas unidades 2^-72."""
    cfg = SYSTEMS[sysname]
    (p1, q1), (p2, q2) = cfg['k1'], cfg['k2']
    E1 = max([abs(M) for o in out for (M, N) in o] + [abs(center_m)])
    E2 = max([abs(N) for o in out for (M, N) in o] + [abs(center_n)])
    P1 = [p1**e*q1**(E1-e) for e in range(E1+1)]
    P2 = [p2**e*q2**(E2-e) for e in range(E2+1)]
    total = 0
    for i, o in enumerate(out):
        acc = 0
        for (M, N), (x, y) in o.items():
            if x or y:
                acc += (abs(x)+abs(y))*P1[abs(M)]*P2[abs(N)]
        total += cfg['eta'][i]*acc
    return Q(total, q1**E1*q2**E2*2**72)


def beta(sysname, d, p, m, j, r_window=None):
    """||B_ref e_{d,p,m,j}|| / ||e||  EXATO. r_window: max L do termo R (None = todos)."""
    cfg = SYSTEMS[sysname]
    if m == 0 and p == 1:
        raise ValueError('m=0 so parte real')
    tab = tables(sysname)[key_for(sysname, d, j)]
    img = input_images(m, j, p)
    out = [dict() for _ in range(cfg['ncomp'])]
    for i, field in enumerate(tab):
        o = out[i]
        for (m0, n0), (vx, vy) in img.items():
            for (m1, n1), (a, b) in field.items():
                k = (m0+m1, n0+n1)
                x, y = o.get(k, (0, 0))
                o[k] = (x+a*vx-b*vy, y+a*vy+b*vx)
    for comp, mult in r_terms(sysname, d, j):
        o = out[comp]
        for sg in ({1, -1} if m else {1}):
            vx, vy = img[sg*m, j]
            L = 0
            while j-1-2*L >= 0 and (r_window is None or L <= r_window):
                amp = MU_INT72*2*mult*(-1)**L
                for nn in {j-1-2*L, -(j-1-2*L)}:
                    x, y = o.get((sg*m, nn), (0, 0))
                    o[sg*m, nn] = (x+amp*vx, y+amp*vy)
                L += 1
    (p1, q1), (p2, q2) = cfg['k1'], cfg['k2']
    norm_in = len(img)*cfg['eta'][d]*Q(p1, q1)**m*Q(p2, q2)**j
    return weighted_norm(sysname, out, m, j)/norm_in


def r_tail(sysname, d):
    """Cota da parte do termo R com L>=50 (fora da janela), relativa, para todo j."""
    cfg, r = SYSTEMS[sysname], r_of(sysname)
    tail = 2*MU*r**101/(1-r*r)
    worst = Q(0)
    for par in (0, 1):
        s = sum((Q(cfg['eta'][c], cfg['eta'][d])*abs(mult) for c, mult in r_terms(sysname, d, 101+par)), Q(0))
        worst = max(worst, s)
    return worst*tail


def _job(args):
    sysname, d, p, m, j, window = args
    return args, beta(sysname, d, p, m, j, window)


def compute_betas(sysname, workers=2):
    cfg = SYSTEMS[sysname]
    jobs = []
    for d in range(cfg['ncomp']):
        mp = cfg['mpar'][d]
        mreps = [m for m in range(0, SUPPORT_M+2) if m % 2 == mp and m <= SUPPORT_M+1]
        # radial: j representativo >=101 por paridade, com janela L<=49
        jreps = [101] if cfg['kinds'][d] == 'K' else [101, 102]
        for m in mreps:
            for j in jreps:
                for p in ((0,) if m == 0 else (0, 1)):
                    jobs.append((sysname, d, p, m, j, 49))
        # Fourier: m representativo >=41 (paridade do componente), todos j<=100
        mrep = SUPPORT_M+1 if (SUPPORT_M+1) % 2 == mp else SUPPORT_M
        for j in range(0, SUPPORT_N):
            if cfg['kinds'][d] == 'K' and j % 2 == 0:
                continue
            for p in (0, 1):
                jobs.append((sysname, d, p, mrep, j, None))
    result = {}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for args, val in ex.map(_job, jobs, chunksize=4):
            result[args[1:]] = val
    return result


# ------------------------------------------------------------ certificado
def certify(sysname, betas):
    cfg, r = SYSTEMS[sysname], r_of(sysname)
    n_comp = cfg['ncomp']
    if sysname == 'A':
        g1 = Q(43390, 2061)
        e_B = EPS_MU*g1+6*Q(4096, 243)*EPS_OM
        e_mu = EPS_MU*2*(1+r)/((1-r)*MU)          # RT (37b): ||D'Q||
    else:
        e_B = EPS_MU*Q(144, 17)+3*EPS_OM
        e_mu = EPS_MU*2*(1+r)/((1-r)*MU)
    per = {}
    for d in range(n_comp):
        k = 1 if cfg['kinds'][d] == 'J' else 2
        mp = cfg['mpar'][d]
        rad = [v for (dd, p, m, j, w), v in betas.items() if dd == d and w == 49]
        beta_R = max(rad)+r_tail(sysname, d)
        mrep = SUPPORT_M+1 if (SUPPORT_M+1) % 2 == mp else SUPPORT_M
        bF = {}
        for j in range(SUPPORT_N):
            vals = [betas[(d, p, mrep, j, None)] for p in (0, 1) if (d, p, mrep, j, None) in betas]
            if vals:
                bF[j] = max(vals)
        beta_inf = max(betas[(d, p, mrep, jj, 49)] for p in (0, 1) for jj in ((101,) if k == 2 else (101, 102)))
        beta_inf += r_tail(sysname, d)
        per[d] = dict(k=k, beta_R=beta_R, beta_inf=beta_inf, bF=bF)
    return dict(sysname=sysname, e_B=e_B, e_mu=e_mu, per=per)


def radial_bound(cert, N):
    """d*_R(N): supremo sobre TODAS as colunas com n>=N (qualquer m)."""
    sysname = cert['sysname']
    cfg, r = SYSTEMS[sysname], r_of(sysname)
    worst, parts = Q(0), {}
    for d, info in cert['per'].items():
        k = info['k']
        if k == 1:
            C = 1+2*SQ2*r/(1-r)
            nu = Q(1)
        else:
            nu = Q(N, N-1)
            C = 1+2*SQ2*nu*r*r/(1-r*r)
        S_hat = (C+1/(4*C))/MU                      # >= (n+1) q para n>=N
        low = SQ2*2*nu/(MU*(N+1))*r**(N-(SUPPORT_N-1))/(1-r**k)
        val = ((info['beta_R']+cfg['s_max']+cert['e_B'])*S_hat/(N+1)
               + max(Q(0), cfg['beta_glob']-info['beta_R'])*low + cert['e_mu'])
        parts[d] = val
        worst = max(worst, val)
    return worst, parts


def fourier_bound(cert, M):
    """d*_F(M): supremo sobre colunas com m>=M (M>=42), qualquer n."""
    sysname = cert['sysname']
    cfg, r = SYSTEMS[sysname], r_of(sysname)
    if M < SUPPORT_M+1:
        raise ValueError('M>=42 exigido')
    worst, parts = Q(0), {}
    b = Q(M, 2)
    for d, info in cert['per'].items():
        k = info['k']
        nu = Q(1) if k == 1 else NU2
        bmax = max(max(info['bF'].values()), info['beta_inf'])
        B = lambda j: info['bF'][j] if j < SUPPORT_N else info['beta_inf']
        phi = Q(0)
        top = SUPPORT_N+220
        for n in range(0, top):
            if k == 2 and n % 2 == 0:
                continue
            acc = G0*B(n)
            i = 1
            while n-i*k >= 0 and i*k <= 220:
                acc += SQ2*nu*r**(i*k)*B(n-i*k)
                i += 1
            if n-i*k >= 0:
                acc += SQ2*nu*bmax*r**(i*k)/(1-r**k)
            phi = max(phi, acc)
        # n>=top: termos com j>=101 usam beta_inf; demais pesam r^(>=220)
        far = G0*info['beta_inf']+SQ2*nu*(info['beta_inf']*r**k/(1-r**k)+bmax*r**220/(1-r**k))
        phi = max(phi, far)
        qcoef = G0+SQ2*nu*r**k/(1-r**k)
        val = (phi+(cfg['s_max']+cert['e_B'])*qcoef)/b+cert['e_mu']
        parts[d] = val
        worst = max(worst, val)
    return worst, parts


def smallest(fn, lo, hi):
    if fn(hi)[0] >= 1:
        return None
    while hi-lo > 1:
        mid = (lo+hi)//2
        if fn(mid)[0] < 1:
            hi = mid
        else:
            lo = mid
    return hi


def dofs(sysname, M, N):
    cfg = SYSTEMS[sysname]
    total = 0
    for d in range(cfg['ncomp']):
        for m in range(M):
            if m % 2 != cfg['mpar'][d]:
                continue
            cnt = sum(1 for n in range(N) if not (cfg['kinds'][d] == 'K' and n % 2 == 0))
            total += cnt*(1 if m == 0 else 2)
    return total


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workers', type=int, default=2)
    parser.add_argument('--out', type=Path, default=ROOT/'build/spectrum/columnwise-exterior.json')
    args = parser.parse_args()
    report = dict(phase=1, method='semianalitico por grau; betas exatos; bounds analiticos de Q',
                  refa_sha256=REFA_SHA, physical_certificate=False,
                  conditional_on_published_RT_ball=True)
    for sysname in ('A', 'sharp'):
        betas = compute_betas(sysname, args.workers)
        cert = certify(sysname, betas)
        NG = smallest(lambda N: radial_bound(cert, N), 150, 4000)
        MG = smallest(lambda M: fourier_bound(cert, M), 42, 4000)
        entry = dict(
            e_B=str(cert['e_B']), e_mu=str(cert['e_mu']),
            per_component={d: dict(beta_R=str(i['beta_R']), beta_R_decimal=float(i['beta_R']),
                                   beta_inf=str(i['beta_inf']), beta_inf_decimal=float(i['beta_inf']),
                                   beta_F_max_decimal=float(max(i['bF'].values())))
                           for d, i in cert['per'].items()},
            N_G_certified=NG, M_G_certified=MG,
            betas_computed=len(betas))
        if NG:
            dR, partsR = radial_bound(cert, NG)
            entry.update(d_radial_at_NG=str(dR), d_radial_at_NG_decimal=float(dR),
                         d_radial_parts_decimal={d: float(v) for d, v in partsR.items()})
            d384 = radial_bound(cert, 384)[0]
            entry.update(d_radial_at_384=str(d384), d_radial_at_384_decimal=float(d384))
        if MG:
            dF, partsF = fourier_bound(cert, MG)
            entry.update(d_fourier_at_MG=str(dF), d_fourier_at_MG_decimal=float(dF),
                         d_fourier_parts_decimal={d: float(v) for d, v in partsF.items()})
        if NG and MG:
            entry['box_dofs'] = dofs(sysname, MG, NG)
        report[sysname] = entry
        print(json.dumps({sysname: {k: v for k, v in entry.items() if 'decimal' in k or k.endswith('certified') or k == 'box_dofs'}}, indent=2), flush=True)
    report['decision_rule'] = 'N_G<=384 em A e sharp => rota de blocos continua; senao inviavel neste hardware'
    NA, NS = report['A']['N_G_certified'], report['sharp']['N_G_certified']
    report['decision_phase1'] = ('continua' if NA and NS and NA <= 384 and NS <= 384
                                 else 'N_G(Fase 1) > 384')
    args.out.write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    main()
