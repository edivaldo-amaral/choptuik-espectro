#!/usr/bin/env python3
"""A folga de Lipschitz (Weyl) do certificado do simbolo, em aritmetica de bolas (Arb).

certify_symbol_numerical_range.folga calcula, por caixa da grade (theta, phi), uma cota da
variacao das formas S e M entre o centro e qualquer ponto da caixa:

    delta_j = L_theta,j dtheta/2 + L_phi,j dphi/2                 (Lipschitz exatos dos campos)
    dS = sum_j |S(e_j)| delta_j   (lados z e -z),  pS = |S(mu/xi = 1)|
    dxi = mu_max a/(b^2 + (a^2 - b^2) min cos^2 phi) dphi/2        (|d(1/xi)/dphi| na faixa)
    folga_lambda = max_faixa ||N^-1/2 Herm(dS + pS dxi) N^-1/2||_F,
    folga_norma  = max_faixa ||(dM + pM dxi) N^-1/2||_F,

em ponto flutuante com fator 1,001. Aqui o mesmo calculo e feito em bolas, e as cotas
finais (t + folga + correcao) sao refeitas com a folga em Arb.
"""
from __future__ import annotations

import argparse
import json
from fractions import Fraction as Q
from pathlib import Path
import sys

from flint import arb, acb, fmpq

sys.path.insert(0, str(Path(__file__).resolve().parent))
import certify_symbol_numerical_range as cs


def A(x) -> arb:
    x = Q(x)
    return arb(fmpq(x.numerator, x.denominator))


def folga_arb(forma, dim, lips, nt, nf, K2, mus):
    pi = arb.pi()
    dth, dph = 2*pi/nt, pi/nf
    delta = [A(Lt)*dth/2 + A(Lp)*dph/2 for Lt, Lp in lips]
    z4 = [acb(0)]*4
    zero = acb(0)
    dS = [[arb(0)]*dim for _ in range(dim)]
    dM = None
    for j in range(4):
        for lado in (0, 1):
            e = [acb(0)]*4; e[j] = acb(1)
            S, M = forma(e, z4, zero) if lado == 0 else forma(z4, e, zero)
            if dM is None:
                dM = [[arb(0)]*len(M[0]) for _ in range(len(M))]
            for r in range(dim):
                for c in range(dim):
                    dS[r][c] += abs(S[r][c])*delta[j]
            for r in range(len(M)):
                for c in range(len(M[0])):
                    dM[r][c] += abs(M[r][c])*delta[j]
    S1, M1 = forma(z4, z4, acb(1))                       # padrao do termo mu/xi
    pS = [[abs(x) for x in linha] for linha in S1]
    pM = [[abs(x) for x in linha] for linha in M1]
    kk = A(K2)
    aa, bb = (kk + 1/kk)/2, (kk - 1/kk)/2
    ns = [1/arb(2).sqrt()] + [arb(1)]*(dim - 1)
    mumax = max(A(m) for m in mus) if len(mus) > 1 else A(mus[0])
    fl = fn = arb(0)
    fl_up = fn_up = 0.0
    for i in range(nf):
        # min de cos^2 na faixa [i pi/nf, (i+1) pi/nf]: 0 se pi/2 esta nela (teste inteiro exato)
        if 2*i <= nf <= 2*(i + 1):
            cos2 = arb(0)
        else:
            c_lo = (pi*i/nf).cos()**2; c_hi = (pi*(i + 1)/nf).cos()**2
            cos2 = arb(min(c_lo.lower(), c_hi.lower()))
            cos2 = arb(max(float(cos2.lower()), 0.0)) if float(cos2.lower()) < 0 else cos2
        dxi = mumax*aa/(bb*bb + (aa*aa - bb*bb)*cos2)*dph/2
        # Frobenius de N^-1/2 Herm(Si) N^-1/2 e de Mi N^-1/2 (entradas >= 0)
        s2 = arb(0)
        for r in range(dim):
            for c in range(dim):
                x = ((dS[r][c] + pS[r][c]*dxi) + (dS[c][r] + pS[c][r]*dxi))/2*ns[r]*ns[c]
                s2 += x*x
        m2 = arb(0)
        for r in range(len(dM)):
            for c in range(dim):
                x = (dM[r][c] + pM[r][c]*dxi)*ns[c]
                m2 += x*x
        fl_up = max(fl_up, float(s2.sqrt().upper()))
        fn_up = max(fn_up, float(m2.sqrt().upper()))
    return fl_up, fn_up


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--operador', choices=('K', 'L'), default='K')
    ap.add_argument('--json', type=Path, required=True, help='certificado combinado ou de faixa')
    ap.add_argument('--nphi', type=int); ap.add_argument('--ntheta', type=int)
    a = ap.parse_args()
    cfg = cs.OPER[a.operador]
    K1, K2, forma, dim = cfg['k1'], cfg['k2'], cfg['forma'], cfg['dim']
    _, lips = cs.carrega(K1, K2)
    d = json.loads(a.json.read_text())
    if 'linhas' in d and list(d['linhas']) != [0, d['nphi']]:
        # certificado em faixas: a folga em Arb e o maximo sobre TODAS as faixas; compara com o
        # maximo gravado entre os JSON irmaos (symbol-X-faixa*.json)
        irmaos = [json.loads(q.read_text()) for q in sorted(a.json.parent.glob(a.json.name.split('faixa')[0] + 'faixa*.json'))]
        d = dict(d, folga_lambda=max(r['folga_lambda'] for r in irmaos), folga_norma=max(r['folga_norma'] for r in irmaos))
    nf = a.nphi or d.get('nphi'); nt = a.ntheta or d.get('ntheta')
    fl, fn = folga_arb(forma, dim, lips, nt, nf, K2, cs.MUS)
    print(f'{a.operador} ({nf}x{nt}): folga em Arb  lambda {fl:.8f}  norma {fn:.8f}')
    if 'folga_lambda' in d:
        print(f'   folga gravada (float x 1,001): lambda {d["folga_lambda"]:.8f}  norma {d["folga_norma"]:.8f}')
        ok = fl <= d['folga_lambda'] and fn <= d['folga_norma']
        tl, tn = A(d['t_lambda']), A(d['t_norma'])
        corr = d['correcao']
        R = float((tl + arb(fl) + arb(corr)).upper()); N = float((tn + arb(fn) + arb(corr)).upper())
        print(f'   refeito: max Re W(-B) <= {R:.6f}   ||B|| <= {N:.6f}   (gravado: {d.get("R_fundo")}, {d.get("norma_fundo")})')
        print('   a folga gravada cobre a de Arb: ' + ('SIM' if ok else 'NAO'))


if __name__ == '__main__':
    main()
