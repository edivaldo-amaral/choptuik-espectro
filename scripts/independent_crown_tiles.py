#!/usr/bin/env python3
"""Recertificacao INDEPENDENTE (auditoria de 16/09) dos tiles da coroa 6x18.

Nao importa modulos do repositorio. Reaproveita do relatorio auditado apenas
os CENTROS e RAIOS dos tiles (e os segmentos, para conferir a geometria).
Tudo o mais e refeito:
  * leitura propria de rt-A-6x18.dat (matriz exportada antes da janela 14-16/09);
  * J_ref=O_{mu_RefA} construido da definicao de RT (22) e Q_ref=J_ref^-1 EXATO
    por substituicao regressiva em racionais complexos (sem residuo 108);
  * pesos RT fortes x eta=(916,4096,243,1233);
  * arredondamentos com erros de coluna calculados EXATAMENTE;
  * parametriz propria V=inv(H(s_j)) em float (so proposta), residuo exato;
  * orcamento do fundo derivado so de constantes RT:
      ||H_true-H_ref|| <= eps_mu K2' + K2 (eps_mu g1_eta + 6 (4096/243) eps_omega),
    K2=K2'=18/mu_ref (RT (37a),(37b)), eps_mu=43*2^-40+2^-277 (RT sec. 5, (60)),
    eps_omega=2^-25+2^-277 ((39b),(60)).
Para cada tile: defeito uniforme a+r*v, com a=||I-V Hhat(s_j)||+||V||(erros),
v=||V Q_SS||, e cota ||H_SS^-1||<=||V||/(1-a-r v).
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as Q
import hashlib
import json
from operator import mul
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
MATRIX = ROOT/'build/spectrum/rt-A-6x18.dat'
MATRIX_SHA = 'a3335d9550f9567db287c1ae105674a5914446f9d3772d3955e3628a747c5846'
REPORT = ROOT/'build/spectrum/crown-contour-6x18-background.json'
MU = Q(722873400, 2**32)
ETA = (916, 4096, 243, 1233)
G_BITS, V_BITS = 64, 56


def read_matrix(path=MATRIX):
    lines = path.read_text().splitlines()
    assert lines[0] == 'CHOPTUIK_RT_SPECTRAL_MATRIX_V1'
    dim = int(lines[1].split()[1])
    a0 = lines.index('BEGIN_ADDRESSES')+1
    addresses = []
    for k in range(dim):
        idx, d, part, m, n = map(int, lines[a0+k].split())
        assert idx == k
        addresses.append((d, part, m, n))
    e0 = lines.index('BEGIN_ENTRIES')+1
    A = [[Q(0)]*dim for _ in range(dim)]
    for line in lines[e0:]:
        f = line.split()
        if not f:
            continue
        r, c, coef = int(f[0]), int(f[1]), int(f[2])
        exp = int(f[3]) if len(f) > 3 else 0
        A[r][c] += Q(coef)*Q(2)**exp
    return addresses, A


def free_matrix(addresses, mu):
    """O_mu em DOFs reais: diag mu(n+1); d_tau = i m/2; 2 mu n abaixo (K: mesma paridade)."""
    look = {a: i for i, a in enumerate(addresses)}
    dim = len(addresses)
    J = [[Q(0)]*dim for _ in range(dim)]
    for col, (d, part, m, n) in enumerate(addresses):
        J[col][col] += mu*(n+1)
        if m:
            other = look[(d, 1-part, m, n)]
            # (i m/2)(x+iy) = -(m/2) y + i (m/2) x
            J[other][col] += Q(m, 2) if part == 0 else -Q(m, 2)
        step = 2 if d == 0 else 1
        for j in range(n-step, -1, -step):
            J[look[(d, part, m, j)]][col] += 2*mu*n
    return J


def exact_inverse_columns(addresses, mu):
    """Q=O_mu^-1 no box, por coluna complexa exata (preserva m, baixa n)."""
    look = {a: i for i, a in enumerate(addresses)}
    dim = len(addresses)
    Qm = [[Q(0)]*dim for _ in range(dim)]
    for col, (d, part, m, n) in enumerate(addresses):
        step = 2 if d == 0 else 1
        re = [Q(0)]*(n+1)
        im = [Q(0)]*(n+1)
        run = [(Q(0), Q(0)), (Q(0), Q(0))]
        for j in range(n, -1, -1):
            if (n-j) % step:
                continue
            dr, di = mu*(j+1), Q(m, 2)
            if step == 1:
                tr, ti = run[0][0]+run[1][0], run[0][1]+run[1][1]
            else:
                tr, ti = run[j % 2]
            rr, ri = (1 if j == n else 0)-tr, -ti
            den = dr*dr+di*di
            re[j], im[j] = (rr*dr+ri*di)/den, (ri*dr-rr*di)/den
            pr, pi = run[j % 2]
            run[j % 2] = (pr+2*mu*j*re[j], pi+2*mu*j*im[j])
        for j in range(n+1):
            if (n-j) % step:
                continue
            if part == 0:   # entrada real: c_j
                vals = ((0, re[j]), (1, im[j]))
            else:           # entrada imaginaria: i c_j
                vals = ((0, -im[j]), (1, re[j]))
            for p, x in vals:
                if x:
                    if m == 0 and p == 1:
                        raise AssertionError('m=0 deve ser real')
                    Qm[look[(d, p, m, j)]][col] = x
    return Qm


def colnorm(M, rows, cols):
    return max(sum(abs(M[i][j]) for i in rows) for j in cols)


def rnd(x: Q) -> int:
    return (2*x.numerator+x.denominator)//(2*x.denominator)


def prepare():
    raw = MATRIX.read_bytes()
    assert hashlib.sha256(raw).hexdigest() == MATRIX_SHA
    addresses, A = read_matrix()
    dim = len(addresses)
    J = free_matrix(addresses, MU)
    Qm = exact_inverse_columns(addresses, MU)
    # Checagem exata J Q = I no box.
    for j in range(dim):
        nz = [(k, Qm[k][j]) for k in range(dim) if Qm[k][j]]
        for i in range(dim):
            s = sum((J[i][k]*x for k, x in nz), Q(0))
            assert s == (1 if i == j else 0)
    w = [(2-(m == 0))*(2-(n == 0))*Q(65, 64)**m*Q(5, 4)**n*ETA[d] for d, part, m, n in addresses]
    Aw = [[A[i][j]*w[i]/w[j] for j in range(dim)] for i in range(dim)]
    Qw = [[Qm[i][j]*w[i]/w[j] for j in range(dim)] for i in range(dim)]
    Jw = [[J[i][j]*w[i]/w[j] for j in range(dim)] for i in range(dim)]
    B = [[Aw[i][j]-Jw[i][j] for j in range(dim)] for i in range(dim)]
    Bt = [[Aw[j][i]-Jw[i][j] for j in range(dim)] for i in range(dim)]
    allidx = range(dim)
    g = 1 << G_BITS
    Ar = [[rnd(g*x) for x in row] for row in Aw]
    Qr = [[rnd(g*x) for x in row] for row in Qw]
    EA = colnorm([[Aw[i][j]-Q(Ar[i][j], g) for j in allidx] for i in allidx], allidx, allidx)
    EQ = colnorm([[Qw[i][j]-Q(Qr[i][j], g) for j in allidx] for i in allidx], allidx, allidx)
    normA = colnorm(Aw, allidx, allidx)
    normQ = colnorm(Qw, allidx, allidx)
    normAr = Q(colnorm(Ar, allidx, allidx), g)
    shell = [i for i, (d, p, m, n) in enumerate(addresses) if m >= 4 or n >= 12]
    # Hr0 = (Ar Qr)_SS, denominador g^2.
    qcols = {j: [(k, Qr[k][j]) for k in allidx if Qr[k][j]] for j in shell}
    Hr0 = [[sum(Ar[i][k]*x for k, x in qcols[j]) for j in shell] for i in shell]
    QrS = [[Qr[i][j] for j in shell] for i in shell]
    errH0 = EA*normQ+normAr*EQ
    g1 = max(sum(ETA[i]*row[j] for i, row in enumerate(gamma1(Q(1)))) / ETA[j] for j in range(4))
    eps_mu = Q(43, 2**40)+Q(1, 2**277)
    eps_om = Q(1, 2**25)+Q(1, 2**277)
    K2 = 2/MU/(1-Q(4, 5))*(1+Q(4, 5))
    bg = eps_mu*K2+K2*(eps_mu*g1+6*Q(4096, 243)*eps_om)
    return dict(dim=dim, shell=len(shell), Hr0=Hr0, QrS=QrS, g=g, EQ=EQ, errH0=errH0, bg=bg,
                normA=normA, normQ=normQ, normB=colnorm(B, allidx, allidx),
                normB_transposed=colnorm(Bt, allidx, allidx), g1=g1, K2=K2)


def gamma1(mu):
    d = mu*Q(40, 9)
    return [[Q(0)]*4, [d, 2*d, Q(0), Q(0)], [d, Q(0), Q(0), Q(0)], [Q(0), Q(0), Q(0), 2*d]]


def matmul(X, Y):
    cols = list(zip(*Y))
    return [[sum(map(mul, row, col)) for col in cols] for row in X]


DATA = None


def init(data):
    global DATA
    DATA = data


def certify_tile(args):
    x, y, radius = args
    d = DATA
    g, n = d['g'], d['shell']
    Hr0, QrS = d['Hr0'], d['QrS']
    H = (np.array(Hr0, dtype=float)/g/g+(float(x)+1j*float(y))*np.array(QrS, dtype=float)/g)
    V = np.linalg.inv(H)
    gv = 1 << V_BITS
    Vr = [[round(float(z.real)*gv) for z in row] for row in V]
    Vi = [[round(float(z.imag)*gv) for z in row] for row in V]
    a = max(x.denominator, y.denominator)
    xn, yn = int(x*a), int(y*a)
    VrH, ViH = matmul(Vr, Hr0), matmul(Vi, Hr0)
    VrQ, ViQ = matmul(Vr, QrS), matmul(Vi, QrS)
    den = gv*g*g*a
    # V Hhat = [Vr Hr0 a + xn g Vr Q - yn g Vi Q] + i [Vi Hr0 a + xn g Vi Q + yn g Vr Q]  / den
    best = 0
    for j in range(n):
        s = 0
        for i in range(n):
            re = a*VrH[i][j]+xn*g*VrQ[i][j]-yn*g*ViQ[i][j]
            im = a*ViH[i][j]+xn*g*ViQ[i][j]+yn*g*VrQ[i][j]
            if i == j:
                re -= den
            s += abs(re)+abs(im)
        best = max(best, s)
    residual = Q(best, den)
    vnorm = Q(max(sum(abs(Vr[i][j])+abs(Vi[i][j]) for i in range(n)) for j in range(n)), gv)
    vq = Q(max(sum(abs(VrQ[i][j])+abs(ViQ[i][j]) for i in range(n)) for j in range(n)), gv*g)
    s_abs = abs(x)+abs(y)
    center = residual+vnorm*(d['errH0']+s_abs*d['EQ']+d['bg'])
    variation = vq+vnorm*d['EQ']
    defect = center+radius*variation
    return dict(center=[str(x), str(y)], radius=str(radius), center_defect=float(center),
                variation=float(variation), uniform_defect=float(defect),
                uniform_defect_exact=str(defect), accepted=defect < 1,
                inverse_bound=float(vnorm/(1-defect)) if defect < 1 else None,
                inverse_bound_exact=str(vnorm/(1-defect)) if defect < 1 else None)


def check_geometry(report):
    tiles = report['tiles']
    verts = [(Q(1, 8), Q(0)), (Q(1, 8), Q(1, 4)), (Q(1), Q(1, 4)), (Q(1), Q(0))]
    cur, leg = verts[0], 0
    for seg in report['upper_boundary_segments']:
        st, en = tuple(map(Q, seg['start'])), tuple(map(Q, seg['end']))
        t = tiles[seg['tile']]
        c, r = tuple(map(Q, t['center'])), Q(t['radius'])
        if st != cur or c != st:
            return False, 'lacuna'
        a, b = verts[leg], verts[leg+1]
        # en no segmento [a,b] e depois de st
        on_line = (a[0] == b[0] == en[0]) or (a[1] == b[1] == en[1])
        between = min(a[0], b[0]) <= en[0] <= max(a[0], b[0]) and min(a[1], b[1]) <= en[1] <= max(a[1], b[1])
        dist2 = (en[0]-st[0])**2+(en[1]-st[1])**2
        if not (on_line and between and 0 < dist2 <= r*r):
            return False, 'segmento fora'
        cur = en
        if cur == b:
            leg += 1
            if leg == 3:
                break
    return leg == 3 and cur == verts[3], 'ok'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tiles', default='all', help="'all' ou lista de indices separada por virgula")
    parser.add_argument('--workers', type=int, default=2)
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    report = json.loads(REPORT.read_text())
    data = prepare()
    head = dict(dimension=data['dim'], shell=data['shell'], errH0=float(data['errH0']),
                EQ=float(data['EQ']), background_budget=str(data['bg']),
                background_budget_decimal=float(data['bg']), g1_eta=str(data['g1']),
                exact_normQ_G=str(data['normQ']), normQ_G_decimal=float(data['normQ']),
                normA_G=float(data['normA']), normB_G=float(data['normB']),
                normB_G_if_transposed=float(data['normB_transposed']),
                K2=float(data['K2']))
    geometry_ok, why = check_geometry(report)
    head['geometry_upper_boundary_covered'] = geometry_ok
    head['geometry_note'] = why
    print(json.dumps(head, indent=2), flush=True)
    indices = range(len(report['tiles'])) if args.tiles == 'all' else [int(i) for i in args.tiles.split(',')]
    jobs = [(Q(report['tiles'][i]['center'][0]), Q(report['tiles'][i]['center'][1]),
             Q(report['tiles'][i]['radius'])) for i in indices]
    results = []
    with ProcessPoolExecutor(max_workers=args.workers, initializer=init, initargs=(data,)) as ex:
        for i, res in zip(indices, ex.map(certify_tile, jobs)):
            res['index'] = i
            res['report_inverse_bound'] = float(Q(report['tiles'][i]['parametrix_norm'])
                                                / (1-Q(report['tiles'][i]['uniform_defect'])))
            print(json.dumps(res), flush=True)
            results.append(res)
    summary = dict(head, tiles=results,
                   all_accepted=all(r['accepted'] for r in results),
                   max_uniform_defect=max(r['uniform_defect'] for r in results),
                   max_inverse_bound=max((r['inverse_bound'] or float('inf')) for r in results),
                   note='centros/raios do relatorio auditado; V, Q, pesos, erros e fundo independentes')
    if args.out:
        args.out.write_text(json.dumps(summary, indent=2)+'\n')


if __name__ == '__main__':
    main()
