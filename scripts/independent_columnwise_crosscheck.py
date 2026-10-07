"""Revisao cruzada (rodada 3, Fisico) do majorante exterior por colunas d*(G).

NAO importa scripts/columnwise_*.py nem independent_component_beta.py.
Usa o motor float64 de independent_numeric_audit.py (operadores reimplementados
do TeX de RT e validados contra as matrizes exportadas a 1e-15) e uma
resolucao triangular PROPRIA de O_mu c = e a partir da definicao do operador
(nao da formula (22) de RT). DIAGNOSTICO EM PONTO FLUTUANTE com margem declarada.

  beta  : invariancias beta(d,p,m,j) em j>=101 e m>=41 (pontos nao usados)
  betaR : max_m beta(j=101/102) por componente, para comparar com beta_R do relatorio
  U     : U_mine(e)=||B_ref Q e||+(s_max+e_B)||Q e||+e_mu nas colunas pedidas
"""
from fractions import Fraction as Fr
import argparse
import json
import sys

import numpy as np

import independent_numeric_audit as N

N.CAPS[0], N.CAPS[1] = 500, 2000


def conv_direct(u, v):
    """Convolucao por soma direta sobre as entradas nao nulas do fator mais esparso
    (sem FFT: o ruido absoluto da FFT e amplificado pelos pesos (5/4)^n)."""
    def sym(f):
        return np.concatenate([f.a[:, :0:-1], f.a], axis=1)
    us, vs = sym(u), sym(v)
    if np.count_nonzero(us) > np.count_nonzero(vs):
        us, vs = vs, us
    out = np.zeros((us.shape[0] + vs.shape[0] - 1, us.shape[1] + vs.shape[1] - 1), dtype=complex)
    ii, kk = np.nonzero(us)
    h0, h1 = vs.shape
    for i, k in zip(ii, kk):
        out[i:i + h0, k:k + h1] += us[i, k] * vs
    Ntot = u.N + v.N
    return N.Fld(out[:, Ntot:])


N.conv = conv_direct
EPS_MU = Fr(43, 2**40) + Fr(1, 2**277)
EPS_OM = Fr(1, 2**25) + Fr(1, 2**277)
SYS = {
    'A': dict(k1=65 / 64, k2=5 / 4, eta=(916, 4096, 243, 1233), kinds='KJJJ', mpar=(0, 0, 0, 1),
              s_max=1.25, ncomp=4),
    'sharp': dict(k1=129 / 128, k2=9 / 8, eta=(1, 1), kinds='KJ', mpar=(0, 0), s_max=0.8, ncomp=2),
}
_BG = {}


def background():
    if 'w' not in _BG:
        mu_fr, raw = N.parse_refa()
        _BG['mu_fr'] = mu_fr
        _BG['mu'], _BG['w'] = N.refa_fields(mu_fr, raw)
    return _BG['mu'], _BG['w']


def constants(system):
    mu = _BG['mu_fr']
    r = Fr(4, 5) if system == 'A' else Fr(8, 9)
    if system == 'A':
        e_B = EPS_MU * Fr(43390, 2061) + 6 * Fr(4096, 243) * EPS_OM
    else:
        e_B = EPS_MU * Fr(144, 17) + 3 * EPS_OM
    e_mu = EPS_MU * 2 * (1 + r) / ((1 - r) * mu)
    return float(e_B), float(e_mu)


def B_apply(system, fields):
    mu, w = background()
    if system == 'A':
        g1 = N.gamma1(fields)
        g2 = N.gamma2(w, fields)
        q = [mu * a + 2.0 * b.times_xi() for a, b in zip(g1, g2)]
        F, _ = N.selection_rt(q)
        return F
    g = N.gamma2_sharp(w, fields)
    return [g[0], mu * fields[0].div_xi() + g[1]]


def wnorm(system, fields, n_ref, m_ref):
    """norma l1 ponderada RT (soma sobre Z^2) dividida por k1^m_ref k2^n_ref."""
    cfg = SYS[system]
    tot = 0.0
    for i, f in enumerate(fields):
        if not np.any(f.a):
            continue
        ms = np.abs(np.arange(-f.M, f.M + 1))
        ns = np.arange(f.N + 1)
        wm = cfg['k1'] ** (ms - m_ref)
        wn = cfg['k2'] ** (ns - n_ref) * np.where(ns == 0, 1.0, 2.0)
        tot += cfg['eta'][i] * float(np.sum((np.abs(f.a.real) + np.abs(f.a.imag)) * wm[:, None] * wn[None, :]))
    return tot


def unit(system, d, p, m, n, values=None):
    cfg = SYS[system]
    comps = [N.Fld.zeros(max(m, 1), n) for _ in range(cfg['ncomp'])]
    z = 1.0 if p == 0 else 1.0j
    if values is None:
        values = {n: 1.0}
    M = comps[d].M
    for j, c in values.items():
        comps[d].a[M + m, j] = z * c
        if m:
            comps[d].a[M - m, j] = np.conj(z * c)
    return comps


def enorm(system, d, m, n):
    cfg = SYS[system]
    return cfg['eta'][d] * (2 if m else 1) * (2 if n else 1)


def beta(system, d, p, m, j):
    out = B_apply(system, unit(system, d, p, m, j))
    return wnorm(system, out, j, m) / enorm(system, d, m, j)


def q_column(system, d, m, n):
    """Resolve O c = e_n (componente d, modo m) pela definicao: O^A tem diagonal
    im/2+mu(j+1) e superdiagonal 2 mu k (k>j); O^B mantem so k-j par."""
    mu, _ = background()
    kind = SYS[system]['kinds'][d]
    size = n + 1
    O = np.zeros((size, size), dtype=complex)
    for j in range(size):
        O[j, j] = 0.5j * m + mu * (j + 1)
        ks = np.arange(j + 1, size)
        vals = 2 * mu * ks
        if kind == 'K':
            vals = np.where((ks - j) % 2 == 0, vals, 0.0)
        O[j, j + 1:] = vals
    rhs = np.zeros(size, dtype=complex)
    rhs[n] = 1.0
    from scipy.linalg import solve_triangular
    return solve_triangular(O, rhs, lower=False)


def U_mine(system, d, p, m, n):
    cfg = SYS[system]
    c = q_column(system, d, m, n)
    vals = {j: c[j] for j in range(n + 1) if c[j] != 0}
    Qe = unit(system, d, p, m, n, vals)
    qn = wnorm(system, Qe, n, m) / enorm(system, d, m, n)
    bq = wnorm(system, B_apply(system, Qe), n, m) / enorm(system, d, m, n)
    e_B, e_mu = constants(system)
    return dict(system=system, d=d, p=p, m=m, n=n, BQ=bq, q=qn, U=bq + (cfg['s_max'] + e_B) * qn + e_mu)


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', choices=('beta', 'betaR', 'U'))
    ap.add_argument('--cols', nargs='*', default=[])
    a = ap.parse_args(argv)
    background()
    if a.cmd == 'beta':
        tests = [('A', 3, 0, 41, [101, 137, 211]), ('A', 3, 1, 41, [101, 137]), ('A', 1, 0, 42, [102, 138, 210]),
                 ('A', 0, 0, 42, [101, 137, 211]), ('sharp', 0, 0, 42, [101, 137, 211]), ('sharp', 1, 1, 42, [102, 138])]
        for system, d, p, m, js in tests:
            vals = {j: beta(system, d, p, m, j) for j in js}
            print(system, 'd', d, 'p', p, 'm', m, {j: f'{v:.15g}' for j, v in vals.items()},
                  'max|dif| =', f'{max(vals.values()) - min(vals.values()):.3e}')
        for system, d, j, ms in (('A', 3, 101, [41, 53, 87]), ('A', 2, 102, [42, 88]), ('A', 0, 137, [42, 88]),
                                 ('sharp', 1, 102, [42, 88]), ('sharp', 0, 211, [42, 64])):
            vals = {(m, p): beta(system, d, p, m, j) for m in ms for p in (0, 1)}
            print(system, 'd', d, 'j', j, {k: f'{v:.15g}' for k, v in vals.items()},
                  'max|dif| =', f'{max(vals.values()) - min(vals.values()):.3e}')
        return
    if a.cmd == 'betaR':
        for system in ('A', 'sharp'):
            cfg = SYS[system]
            for d in range(cfg['ncomp']):
                best, arg = 0.0, None
                js = (101,) if cfg['kinds'][d] == 'K' else (101, 102)
                for m in range(0, 43):
                    if m % 2 != cfg['mpar'][d]:
                        continue
                    for p in ((0,) if m == 0 else (0, 1)):
                        for j in js:
                            v = beta(system, d, p, m, j)
                            if v > best:
                                best, arg = v, (m, p, j)
                print(system, 'd', d, 'max_m beta(j>=101) =', f'{best:.12f}', 'em (m,p,j) =', arg, flush=True)
        return
    for spec in a.cols:
        system, d, p, m, n = spec.split(':')
        res = U_mine(system, int(d), int(p), int(m), int(n))
        print(json.dumps(res), flush=True)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
