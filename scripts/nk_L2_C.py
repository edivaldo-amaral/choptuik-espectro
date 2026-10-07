#!/usr/bin/env python3
"""S3b, L, modo rigoroso, ETAPA C: montagem do certificado NK bordejado de L (racionais exatos).

Le Z.json (etapa A: nivel Z') e T.json (etapa B: niveis de cauda e far). Norma do certificado:
||x||_v = max_l ||x_l||/v_l (niveis: Z' = (y_Z, lambda), D_J, R_J, far). Com A = diag(Mh^-1,
Hh_J^-1, I, I) e N as cotas por nivel de E = I - A DF(x0):

  PERRON: (N v)_l <= theta v_l nas linhas de cauda e far; na linha de Z', cuja entrada vem da
          norma CONJUNTA a sobre T, rhoZ v_0 + a sqrt(sum_{j em SZ} v_j^2) + epsZ v_far <= theta v_0,
          verificada por quadrados. Entao ||I - A DF(x0)||_v <= theta.
  2a ORDEM: DF(x) - DF(x0) = [(lambda - lambda0) Q0 dy + dlambda Q0 (y - y0); 0]. Q0 e diagonal em m
          e triangular superior em n: P_l Q0 P_c = 0 se l e c nao tem modos em comum (e P_far Q0
          P_c = 0 para c dentro de F). Com q_lc >= ||P_l Q0 P_c||, para ||x - x0||_v <= r:
          ||A (DF(x) - DF(x0))||_v <= Z2 r,  Z2 = 2 v_0 max_l ||A_l|| (sum_c q_lc v_c)/v_l.
  NK:     Y + theta r + Z2 r^2 <= r e theta + 2 Z2 r < 1, Y = max_l ||A_l P_l F(x0)||/v_l: existe um
          unico zero x* = (y*, lambda*) na bola, e DF(x*) e invertivel -- lambda* e autovalor
          ALGEBRICAMENTE SIMPLES de L (S3B_ROTA.md, 8), |lambda* - lambda0| <= v_0 r.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Fr
import json
import math
from pathlib import Path

import numpy as np

ARRED = 1 + Fr(1, 10**9)


def raiz_sup(x: Fr) -> Fr:
    """Racional >= sqrt(x)."""
    s = Fr(math.sqrt(float(x)))*(1 + Fr(1, 10**12))
    while s*s < x:
        s *= 1 + Fr(1, 10**9)
    return s


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--dir', type=Path, default=Path('build/s3b/rig'))
    a = ap.parse_args()
    z = json.loads((a.dir/'Z.json').read_text()); t = json.loads((a.dir/'T.json').read_text())
    nome = t['nome']; K = len(nome); iF = K - 1
    Nf = np.array(t['N'])
    Nf[0, 0] = z['rhoZ']; Nf[0, iF] = z['epsZ']
    SZ = t['SZ']
    # v: potencia nao linear em ponto flutuante (so propoe; a verificacao e exata)
    def Fv(v):
        w = Nf @ v
        w[0] = Nf[0, 0]*v[0] + z['aZ']*math.sqrt(sum(v[j]**2 for j in SZ)) + Nf[0, iF]*v[iF]
        return w
    v = np.ones(K)
    for _ in range(5000):
        w = Fv(v) + 1e-14; v = w/w.max()
    v = np.maximum(v, 1e-8)
    thf = float(max(Fv(v)/v))
    # verificacao exata
    N = [[Fr(float(x))*ARRED for x in linha] for linha in Nf.tolist()]
    aZ = Fr(z['aZ'])*ARRED
    vq = [Fr(float(x)) for x in v]
    theta = Fr(thf)*(1 + Fr(1, 10**7))
    for i in range(1, K):
        assert sum(N[i][j]*vq[j] for j in range(K)) <= theta*vq[i], f'linha {nome[i]}'
    base = (theta - N[0][0])*vq[0] - N[0][iF]*vq[iF] - sum(N[0][j]*vq[j] for j in range(1, K - 1) if j not in SZ)
    assert base >= 0 and base*base >= aZ*aZ*sum(vq[j]**2 for j in SZ), 'linha de Z\''
    assert theta < 1, f'theta = {float(theta):.6f} >= 1'
    print(f'PERRON verificado em racionais: theta = {float(theta):.6f}  ({K} niveis)')
    # 2a ordem
    niv = t['niveis']
    MZ = z['Z'][0]
    modos = [set(range(-(MZ - 1), MZ))] + [set(lv['modos']) for lv in niv] + [None]
    nA = [Fr(z['tV'])] + [Fr(lv['nV']) if lv['tipo'] == 'D' else Fr(1) for lv in niv] + [Fr(1)]
    nQ = [Fr(z['nQZ'])] + [Fr(lv['nQ']) for lv in niv] + [Fr(z['Qfar'])]
    def q(l, c):
        if c == iF:
            return nQ[c]                      # Q0 P_far desce para dentro de F
        if l == iF:
            return Fr(0)
        return nQ[c] if modos[l] & modos[c] else Fr(0)
    Z2 = 2*vq[0]*max(nA[l]*sum(q(l, c)*vq[c] for c in range(K))/vq[l] for l in range(K))*ARRED
    # residuo
    Yl = [Fr(z['YZ'])] + [Fr(float(y)) for y in t['Y']] + [Fr(0)]
    Y = max(Yl[l]/vq[l] for l in range(K))*ARRED
    # raio
    d = (1 - theta)**2 - 4*Z2*Y
    assert d > 0, 'o NK nao fecha: discriminante <= 0'
    rf = (float(1 - theta) - math.sqrt(float(d)))/(2*float(Z2))
    r = Fr(rf)*(1 + Fr(1, 10**6))
    assert Y + theta*r + Z2*r*r <= r, 'NK: T nao leva a bola nela mesma'
    assert theta + 2*Z2*r < 1, 'NK: T nao contrai na bola'
    lam_err = vq[0]*r
    nv2 = raiz_sup(sum(x*x for x in vq))
    print(f'Z2 = {float(Z2):.3e};  Y = {float(Y):.3e};  raio NK r = {float(r):.3e}')
    print(f'NK VERIFICADO: existe um unico zero (y*, lambda*) com ||x* - x0||_v <= {float(r):.3e};'
          f' |lambda* - {z["lam0"]}| <= {float(lam_err):.3e}; lambda* algebricamente simples.')
    out = dict(theta=str(theta), v=[str(x) for x in vq], nome=nome, Z2=str(Z2), Y=str(Y), r=str(r),
               lambda0=z['lam0'], erro_lambda=str(lam_err), norma_l2_v=str(nv2),
               erro_y_l2=str(nv2*r), certificado=True)
    (a.dir/'C.json').write_text(json.dumps(out, indent=1) + '\n')
    print(f'||y* - y0||_2 <= ||v||_2 r = {float(nv2*r):.3e}  (C.json)')


if __name__ == '__main__':
    main()
