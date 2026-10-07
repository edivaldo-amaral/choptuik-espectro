"""H3 do testemunho de constraints: cota RACIONAL de ||C(s)|| com perda de raio.

Fisico, rodada 3 (17/09/2026). Ver docs/CONSTRAINT_WITNESS_NORMS.md.

Decomposicao (verificada universalmente em test_independent_constraint_operator_bound.py):
  C(s)h = D(s)h + A_w h,
  D(s)h = ( mu(xi dxi+2)h1 , -xi(dtau+s+mu)h2 + mu(xi dxi - xi^2 dxi)h2 + 2mu(h2-h3) ),
  A_w h = ( 2 O(xi Gamma2_1(w,h)) , -2 P(xi Gamma2_3(w,h)) ).
Entrada: Y+ = pesos (65/64,5/4), SEM eta (norma RT de SSPACE). Saida: Ysharp- =
(129/128,9/8). Normas l1 com soma por coluna (RT); complexificacao com modulo.

Parte EXATA (Fraction): c0 >= sup_col ||C_true(0) e||/||e||, uniforme na bola RT,
com projecao opcional na caixa de saida m<M_F; ||C1|| <= 9/8.
Parte DIAGNOSTICA (float): c_ref = ||Pi_F C_ref(s0) v|| / ||v||_{Y+} no 12x36.
"""
from fractions import Fraction as Fr
import argparse
import json
import math
import sys

import independent_numeric_audit as NA

A, B = Fr(65, 64), Fr(5, 4)          # pesos fortes
A2, B2 = Fr(129, 128), Fr(9, 8)      # pesos sharp menores
MU_REF = Fr(722873400, 2 ** 32)
EPS_MU = Fr(43, 2 ** 40) + Fr(1, 2 ** 277)
EPS_OM = Fr(1, 2 ** 25) + Fr(1, 2 ** 277)
MU_STAR = MU_REF + EPS_MU
Q_A, Q_B = A2 / A, B2 / B            # 129/130, 9/10
XI_NORM = B2                          # ||xi||_{Y-} = 9/8 (coluna n=0)
N0 = 200


def fourier_sup(M_F=None):
    """max_{0<=m<M_F} (m/2) Q_A^m; sem projecao, max em m<=130 (decrescente depois,
    pois (m+1)/m * 129/130 < 1 para m>=130)."""
    top = 131 if M_F is None else M_F
    best = Fr(0)
    for m in range(top):
        best = max(best, Fr(m, 2) * Q_A ** m)
    return best


# ---------------------------------------------- Chebyshev exato (normalizacao RT)
def xi_apply(v):
    out = {}
    for n, c in v.items():
        if n == 0:
            out[1] = out.get(1, 0) + c / 2
        else:
            out[n + 1] = out.get(n + 1, 0) + c / 2
            if n == 1:
                out[0] = out.get(0, 0) + c
            else:
                out[n - 1] = out.get(n - 1, 0) + c / 2
    return out


def xidxi_apply(v):
    out = {}
    for n, c in v.items():
        out[n] = out.get(n, 0) + n * c
        for j in range(n - 2, -1, -2):
            out[j] = out.get(j, 0) + 2 * n * c
    return out


def add(*terms):
    out = {}
    for coef, v in terms:
        for n, c in v.items():
            out[n] = out.get(n, 0) + coef * c
    return out


def nnorm(v, b):
    return sum(abs(c) * (1 if n == 0 else 2) * b ** n for n, c in v.items())


def D_column_ratios():
    """Por n<=N0: (k=1) ||(xi dxi+2)e||, (k=2) ||xi e|| e ||(xi dxi - xi^2 dxi + 2 - xi)e||,
    todos em Y- divididos por ||e||_{Y+} (m=0 e o pior caso; Q_A^m<=1)."""
    k1, x2, y2 = {}, {}, {}
    for n in range(N0 + 1):
        e = {n: Fr(1)}
        den = nnorm(e, B)
        if n % 2 == 1:
            k1[n] = nnorm(add((1, xidxi_apply(e)), (2, e)), B2) / den
        xe = xi_apply(e)
        xd = xidxi_apply(e)
        x2[n] = nnorm(xe, B2) / den
        y2[n] = nnorm(add((1, xd), (-1, xi_apply(xd)), (2, e), (-1, xe)), B2) / den
    return k1, x2, y2


def tail_bounds(n):
    """Majorantes decrescentes para n>=N0+1 (ver doc)."""
    q = Q_B
    S = q ** n * (n + 2 * n / (B2 ** 2 - 1))
    T1 = q ** n * ((n + 2) + 2 * n / (B2 ** 2 - 1))
    X = q ** n * (B2 + 1 / B2) / 2
    Y = S * (1 + XI_NORM) + 2 * q ** n + XI_NORM * q ** n
    return T1, X, Y


# ------------------------------------------------ fundo: normas das combinacoes
def background_combo_norms():
    mu, fields = NA.parse_refa()
    F = []
    for f in fields:
        F.append({(m, n): (Fr(r) * Fr(2) ** t, Fr(i) * Fr(2) ** t) for (m, n), (r, i, t) in f['exact'].items()})
    keys = set().union(*[set(f) for f in F])
    pw = lambda n: -1 if n % 2 else 1

    def combo(coefs):
        """coefs: lista (coef, comp(1..4), usa_P)."""
        tot = Fr(0)
        for (m, n) in keys:
            re = im = Fr(0)
            for c, comp, useP in coefs:
                val = F[comp - 1].get((m, n))
                if val:
                    s = c * (pw(n) if useP else 1)
                    re += s * val[0]
                    im += s * val[1]
            if re or im:
                tot += (abs(re) + abs(im)) * (1 if m == 0 else 2) * (1 if n == 0 else 2) * A2 ** m * B2 ** n
        return tot, sum(abs(c) for c, _, _ in coefs)

    C = {
        'k1': [(3, 2, True), (-2, 3, True), (1, 1, False), (-1, 2, False)],     # 1/4(...)*h1
        'k2e_r1': [(1, 2, True), (2, 3, True), (1, 1, False), (1, 2, False)],   # 1/2(...)*E h2
        'k2e_r2': [(-2, 2, True), (2, 3, True)],
        'k2o_r1': [(-2, 1, False)],
        'k2o_r2': [(2, 2, True), (-2, 3, True)],
        'k3e_r1': [(1, 2, True), (-1, 1, False), (1, 2, False)],
        'k3e_r2': [(2, 2, True)],
        'k3o_r1': [(1, 1, False), (-1, 2, True), (-1, 2, False)],
        'k3o_r2': [(-2, 2, True)],
        'k4_r1': [(2, 4, True)],
        'k4_r2': [(2, 4, True)],
    }
    return {name: combo(c) for name, c in C.items()}


def A_bounds(norms):
    """||A_w e||/||e|| (Y- <- Y-, logo <= Y- <- Y+), com a bola RT (eps_omega por componente)."""
    val = lambda name: norms[name][0] + norms[name][1] * EPS_OM
    xi = XI_NORM
    return {
        1: xi * Fr(1, 2) * val('k1'),
        2: xi * max(val('k2e_r1') + val('k2e_r2'), val('k2o_r1') + val('k2o_r2')),
        3: xi * max(val('k3e_r1') + val('k3e_r2'), val('k3o_r1') + val('k3o_r2')),
        4: xi * (val('k4_r1') + val('k4_r2')),
    }


def c0_bound(M_F=None, verbose=False):
    fs = fourier_sup(M_F)
    k1, x2, y2 = D_column_ratios()
    T1, X, Y = tail_bounds(N0 + 1)
    D1 = MU_STAR * max(max(k1.values()), T1)
    D2 = max(max(fs * x2[n] + MU_STAR * y2[n] for n in x2), fs * X + MU_STAR * Y)
    D3 = 2 * MU_STAR
    Ab = A_bounds(background_combo_norms())
    per = {1: D1 + Ab[1], 2: D2 + Ab[2], 3: D3 + Ab[3], 4: Ab[4]}
    c0 = max(per.values())
    if verbose:
        arg1 = max(k1, key=k1.get)
        arg2 = max(x2, key=lambda n: fs * x2[n] + MU_STAR * y2[n])
        print(f'  sup_m (m/2)(129/130)^m [m<{M_F or "inf"}] = {float(fs):.6f}')
        print(f'  D: k=1 {float(D1):.5f} (n={arg1}); k=2 {float(D2):.5f} (n={arg2}); k=3 {float(D3):.5f}')
        print('  A_w (com bola RT): ' + ', '.join(f'k={k}: {float(v):.5f}' for k, v in Ab.items()))
        print('  por componente de entrada: ' + ', '.join(f'k={k}: {float(v):.5f}' for k, v in per.items()))
    return c0, per


def rational_ceiling(x, den=10 ** 6):
    return Fr(math.ceil(x * den), den)


# --------------------------------------------------------- c_ref (diagnostico)
def c_ref(tag='12x36', targets=(0.401, 0.0414 + 0.5j), boxes=((12, 36), (4, 12))):
    import numpy as np
    SPEC = NA.SPEC
    L, addr = NA.read_rt_matrix(SPEC / f'rt-A-{tag}.dat', 'CHOPTUIK_RT_SPECTRAL_MATRIX_V1')
    C0, C1, caddr = NA.read_constraints(SPEC / f'rt-A-{tag}.dat.constraints', len(L))
    la, vr = np.linalg.eig(L)
    s_all = -la
    wv = np.array([(2 - (m == 0)) * (2 - (n == 0)) * float(A) ** m * float(B) ** n for _, _, m, n in addr])
    wc = np.array([(2 - (m == 0)) * (2 - (n == 0)) * float(A2) ** m * float(B2) ** n for _, _, m, n in caddr])
    out = []
    for t in targets:
        i = int(np.argmin(np.abs(s_all - t)))
        v = vr[:, i]
        v = v / np.sum(np.abs(v) * wv)
        c = (C0 + s_all[i] * C1) @ v
        row = dict(s0=[float(s_all[i].real), float(s_all[i].imag)])
        for (MF, NF) in boxes:
            mask = np.array([(m < MF and n < NF) for _, _, m, n in caddr])
            row[f'c_ref_F{MF}x{NF}'] = float(np.sum(np.abs(c[mask]) * wc[mask]))
        out.append(row)
    return out


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument('--no-cref', action='store_true')
    ap.add_argument('--out')
    a = ap.parse_args(argv)
    report = dict(norms='entrada (65/64,5/4) sem eta; saida (129/128,9/8)', exact=True)
    for M_F in (None, 12, 4):
        print(f'caixa de saida m<{M_F or "inf"}:')
        c0, per = c0_bound(M_F, verbose=True)
        c_op = c0 + Fr(5, 4) * XI_NORM
        c0r, copr = rational_ceiling(c0), rational_ceiling(c_op)
        assert c0r >= c0 and copr >= c_op
        print(f'  c0 <= {c0r} ; c_op(|s|<=5/4) = c0 + (5/4)(9/8) <= {copr} = {float(copr):.6f}')
        report[f'M_F={M_F}'] = dict(c0_upper=str(c0r), c_op_s_le_5_4_upper=str(copr),
                                    per_component={k: float(v) for k, v in per.items()})
    if not a.no_cref:
        rows = c_ref()
        for row in rows:
            s0 = complex(*row['s0'])
            for M_F, key in ((None, 'c_ref_F12x36'), (12, 'c_ref_F12x36'), (4, 'c_ref_F4x12')):
                c0 = Fr(report[f'M_F={M_F}']['c0_upper'])
                cop = float(c0) + abs(s0) * float(XI_NORM)
                row[f'ratio_{key}_M_F={M_F}'] = row[key] / cop
            print(json.dumps(row))
        report['c_ref_diagnostic_float'] = rows
    if a.out:
        with open(a.out, 'w') as fh:
            json.dump(report, fh, indent=2)
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
