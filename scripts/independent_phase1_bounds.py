"""Verificacao independente das cotas analiticas da Fase 1: d_R(N) e d_F(M).

Fisico, rodada 5b (17/09/2026). NAO importa scripts/columnwise_*.py: os beta
vem do motor proprio (independent_columnwise_crosscheck, TeX de RT, convolucao
direta). Os lemas de coluna sao rederivados em papel na secao 12 de
docs/AUDITORIA_FISICA_16SET.md; aqui eles sao aplicados com as MINHAS constantes.

  d_R(N) = max_d [ (beta_R + s_max + e_B) S(N)/(N+1)
                   + max(0, beta_glob - beta_R) low(N) + e_mu ],
  S(N) = (C + 1/(4C))/mu,  C = 1 + 2 sqrt2 nu R_k,  R_k = r^k/(1-r^k),
  nu = 1 (J), N/(N-1) (K),  low(N) = sqrt2 2 nu r^(N-100) / (mu (N+1) (1-r^k)).

  d_F(M) = max_d [ (phi_d + (s_max+e_B) qcoef_d) / (M/2) + e_mu ],
  phi_d = sup_n [ G0 B(n) + sum_{i>=1} sqrt2 nu_F r^(ik) B(n-ik) ],
  qcoef_d = G0 + sqrt2 nu_F r^k/(1-r^k),  G0 = (1+sqrt2)/2.

Constantes EXTERNAS (do Matematico, nao auditadas por mim): beta_glob, e_B, e_mu.
Os beta sao diagnostico em float (concordam com os exatos dele a ~1e-9).
"""
from fractions import Fraction as Fr
import argparse
import json
import math
import sys
from pathlib import Path

import independent_columnwise_crosscheck as CC

SQ2 = math.sqrt(2)
G0 = (1 + SQ2) / 2
NU2_THEIRS = 2 / math.sqrt(3)
MU = float(Fr(722873400, 2 ** 32))
EPS_MU = float(Fr(43, 2 ** 40) + Fr(1, 2 ** 277))
EPS_OM = float(Fr(1, 2 ** 25) + Fr(1, 2 ** 277))
PAR = {
    'A': dict(r=0.8, s_max=1.25, beta_glob=5191 / 500,
              e_B=EPS_MU * (43390 / 2061) + 6 * (4096 / 243) * EPS_OM,
              eta=(916, 4096, 243, 1233)),
    'sharp': dict(r=8 / 9, s_max=0.8, beta_glob=462237 / 100000,
                  e_B=EPS_MU * (144 / 17) + 3 * EPS_OM, eta=(1, 1)),
}
for s, p in PAR.items():
    p['e_mu'] = EPS_MU * 2 * (1 + p['r']) / ((1 - p['r']) * MU)
CACHE = Path('<scratch>'
             '/scratchpad/fisico/phase1_betas.json')


def r_tail(system, d):
    """Cauda do termo mu R_x Gamma1 com L>=50 (relativa), para todo j>=101."""
    cfg, r = CC.SYS[system], PAR[system]['r']
    tail = 2 * MU * r ** 101 / (1 - r * r)
    worst = 0.0
    for j in (101, 102):
        if system == 'A':
            if d == 0:
                terms = [(1, 1), (2, -1)]
            elif d in (1, 3) and j % 2 == 1:
                terms = [(d, 2)]
            else:
                terms = []
        else:
            terms = [(1, 1)] if d == 0 else []
        worst = max(worst, sum(cfg['eta'][c] / cfg['eta'][d] * abs(mult) for c, mult in terms))
    return worst * tail


def betas(system):
    key = system
    data = json.loads(CACHE.read_text()) if CACHE.exists() else {}
    if key in data:
        e = data[key]
        return e['beta_R'], {int(k): v for k, v in e['bF'].items()}, e['beta_inf']
    cfg = CC.SYS[system]
    beta_R, bF, beta_inf = {}, {}, {}
    for d in range(cfg['ncomp']):
        kind, mp = cfg['kinds'][d], cfg['mpar'][d]
        js = (101,) if kind == 'K' else (101, 102)
        best = 0.0
        for m in range(0, 43):
            if m % 2 != mp:
                continue
            for p in ((0,) if m == 0 else (0, 1)):
                for j in js:
                    best = max(best, CC.beta(system, d, p, m, j))
        beta_R[d] = best + r_tail(system, d)
        mrep = 42 if 42 % 2 == mp else 41
        tab = {}
        for j in range(0, 101):
            if kind == 'K' and j % 2 == 0:
                continue
            tab[j] = max(CC.beta(system, d, p, mrep, j) for p in (0, 1))
        bF[d] = tab
        beta_inf[d] = max(CC.beta(system, d, p, mrep, j) for p in (0, 1) for j in js) + r_tail(system, d)
        print(f'  {system} d={d}: beta_R={beta_R[d]:.9f} beta_inf={beta_inf[d]:.9f} '
              f'max bF={max(tab.values()):.9f}', flush=True)
    data[key] = dict(beta_R=beta_R, bF={d: {str(j): v for j, v in t.items()} for d, t in bF.items()},
                     beta_inf=beta_inf)
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    CACHE.write_text(json.dumps(data))
    return beta_R, {d: bF[d] for d in bF}, beta_inf


def d_radial(system, N, beta_R):
    cfg, par = CC.SYS[system], PAR[system]
    r = par['r']
    out = {}
    for d in range(cfg['ncomp']):
        k = 2 if cfg['kinds'][d] == 'K' else 1
        nu = N / (N - 1) if k == 2 else 1.0
        C = 1 + 2 * SQ2 * nu * r ** k / (1 - r ** k)
        S = (C + 1 / (4 * C)) / MU
        low = SQ2 * 2 * nu / (MU * (N + 1)) * r ** (N - 100) / (1 - r ** k)
        out[d] = ((beta_R[d] + par['s_max'] + par['e_B']) * S / (N + 1)
                  + max(0.0, par['beta_glob'] - beta_R[d]) * low + par['e_mu'])
    return max(out.values()), out


def d_fourier(system, M, bF, beta_inf, nu_K=NU2_THEIRS, top=400):
    cfg, par = CC.SYS[system], PAR[system]
    r, b = par['r'], M / 2
    out = {}
    for d in range(cfg['ncomp']):
        k = 2 if cfg['kinds'][d] == 'K' else 1
        nu = nu_K if k == 2 else 1.0
        tab = bF[d]
        bmax = max(max(tab.values()), beta_inf[d])
        B = lambda j: tab[j] if j <= 100 else beta_inf[d]
        phi = 0.0
        for n in range(0, top):
            if k == 2 and n % 2 == 0:
                continue
            acc, i = G0 * B(n), 1
            while n - i * k >= 0 and i * k <= 220:
                acc += SQ2 * nu * r ** (i * k) * B(n - i * k)
                i += 1
            if n - i * k >= 0:
                acc += SQ2 * nu * bmax * r ** (i * k) / (1 - r ** k)
            phi = max(phi, acc)
        far = G0 * beta_inf[d] + SQ2 * nu * (beta_inf[d] * r ** k / (1 - r ** k)
                                             + bmax * r ** 220 / (1 - r ** k))
        phi = max(phi, far)
        qcoef = G0 + SQ2 * nu * r ** k / (1 - r ** k)
        out[d] = (phi + (par['s_max'] + par['e_B']) * qcoef) / b + par['e_mu']
    return max(out.values()), out


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument('--N', type=int, nargs='*', default=[384, 606, 900])
    ap.add_argument('--M', type=int, nargs='*', default=[107, 160])
    a = ap.parse_args(argv)
    CC.background()
    for system in ('A', 'sharp'):
        print(f'== {system}')
        beta_R, bF, beta_inf = betas(system)
        beta_R = {int(k): v for k, v in beta_R.items()}
        bF = {int(d): {int(j): v for j, v in t.items()} for d, t in bF.items()}
        beta_inf = {int(k): v for k, v in beta_inf.items()}
        for N in a.N:
            val, per = d_radial(system, N, beta_R)
            print(f'  d_R({N}) = {val:.6f}   por componente: '
                  + ', '.join(f'{d}:{v:.5f}' for d, v in per.items()))
        for M in a.M:
            val, per = d_fourier(system, M, bF, beta_inf)
            v1, _ = d_fourier(system, M, bF, beta_inf, nu_K=1.0)
            print(f'  d_F({M}) = {val:.6f} (com nu_K=2/sqrt3 do autor); {v1:.6f} (com nu_K=1, provado)'
                  + '   por componente: ' + ', '.join(f'{d}:{v:.5f}' for d, v in per.items()))
        M0 = a.M[0]
        val0, _ = d_fourier(system, M0, bF, beta_inf)
        for M in a.M[1:]:
            val, _ = d_fourier(system, M, bF, beta_inf)
            lhs = (val - PAR[system]['e_mu']) * M
            rhs = (val0 - PAR[system]['e_mu']) * M0
            print(f'  forma 1/M: (d_F({M})-e_mu)*{M} = {lhs:.9f} vs (d_F({M0})-e_mu)*{M0} = {rhs:.9f}'
                  f'  [diferenca relativa {abs(lhs-rhs)/rhs:.1e}]')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
