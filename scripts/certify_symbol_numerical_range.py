#!/usr/bin/env python3
"""Cota RIGOROSA, pelo simbolo pontual, da borda do campo de valores e da norma
da parte de fundo de K (constraints) ou de L (selecionado, F3).

O espaco. L^2 do toro de raios (kappa1, kappa2), coeficientes simetricos com
sinal: ||x||^2 = sum_{m,n em Z} kappa1^(2m) kappa2^(2n) |x_mn|^2. Nele l^1(w) de
RT esta contido com norma <= 1, e multiplicacao por campos, a reflexao xi -> -xi
(z -> -z) e a divisao por xi (1/xi(z)) sao todas PONTUAIS. Logo
    max Re W(-B) <= sup_toro lambda_max(N^-1/2 Herm(-S) N^-1/2),
    ||B||         <= sup_toro ||M N^-1/2||,
com S e M as formas de symbol_forms.py (valores no par z, -z; a componente
impar conta em dobro: N = diag(2, 1, ..., 1)).

Rigor. Centros das caixas em Arb (python-flint); o ponto flutuante so propoe o
limiar t, e em cada centro verifica-se t I - H > 0 (e t^2 I - M*M > 0) por LDL^T
em bolas. Folga por caixa: Weyl, com a norma de Frobenius das variacoes das
entradas, a partir de Lipschitz EXATOS dos campos (|re|+|im|) e de uma cota local
de |d(1/xi)/dphi| por faixa. mu em [1/6, 17/100] pelos extremos (convexidade).
A grade pode ser fatiada em faixas de phi (--linhas i0:i1) e as faixas combinadas
com --combinar.
"""
from __future__ import annotations

import argparse
import json
from fractions import Fraction as Q
import math
from pathlib import Path
import sys

import numpy as np
from flint import arb, acb, acb_mat, acb_poly, fmpq

sys.path.insert(0, str(Path(__file__).resolve().parent))
from background_couplings import RT_BACKGROUND_ERROR
from component_exterior import read_fields
from symbol_forms import forma_K, forma_L

REF = Path('.cache/rt-1203.3766v1/sourcecode/RefA.dat')
MUS = (Q(1, 6), Q(17, 100))
NM = 41
OPER = {'K': dict(k1=Q(129, 128), k2=Q(9, 8), forma=forma_K, dim=3),
        'L': dict(k1=Q(65, 64), k2=Q(5, 4), forma=forma_L, dim=7)}


def A(x: Q) -> arb:
    return arb(fmpq(x.numerator, x.denominator))


def carrega(K1, K2):
    campos = read_fields(REF)
    coefs, lips = [], []
    for f in campos:
        c = [[[Q(0), Q(0)] for _ in range(101)] for _ in range(NM)]
        for (m, n, part), v in f.items():
            c[m][n][part] += v
        coefs.append(c)
        Lt = Lp = Q(0)
        for m in range(NM):
            for n in range(101):
                re, im = c[m][n]
                if re or im:
                    a = abs(re) + abs(im)
                    cm = 1 if m == 0 else 2
                    Lt += cm*m*K1**m*(1 if n == 0 else K2**n + 1/K2**n)*a
                    Lp += cm*K1**m*n*(K2**n + 1/K2**n)*a
        lips.append((Lt, Lp))
    return coefs, lips


def polinomios(coefs):
    P, Pc = [], []
    for c in coefs:
        P.append([acb_poly([acb(A(re), A(im)) for re, im in c[m]]) for m in range(NM)])
        Pc.append([acb_poly([acb(A(re), -A(im)) for re, im in c[m]]) for m in range(NM)])
    return P, Pc


def linha_A(P, Pc, z: acb):
    """Matriz 8 x 81: linhas (sinal de z, campo j), colunas m = -40..40."""
    linhas = []
    for zz in (z, -z):
        iz = 1/zz
        for j in range(4):
            lin = [acb(0)]*(2*NM - 1)
            for m in range(NM):
                a0 = P[j][m][0] if P[j][m].degree() >= 0 else acb(0)
                lin[NM-1+m] = P[j][m](zz) + P[j][m](iz) - a0
                if m:
                    b0 = Pc[j][m][0] if Pc[j][m].degree() >= 0 else acb(0)
                    lin[NM-1-m] = Pc[j][m](zz) + Pc[j][m](iz) - b0
            linhas.append(lin)
    return acb_mat(linhas)


def pd(T):
    """LDL^T em bolas: todos os pivos com cota inferior > 0?"""
    n = len(T); T = [row[:] for row in T]
    for k in range(n):
        p = T[k][k].real
        if not p > 0:
            return False
        for i in range(k+1, n):
            f = T[i][k]/T[k][k]
            for j in range(k, n):
                T[i][j] = T[i][j] - f*T[k][j]
    return True


def folga(forma, dim, lips, nt, nf, linhas, K2, mus):
    """Folga de Weyl por caixa, maxima sobre as faixas pedidas (em float, com margem)."""
    dth, dph = 2*math.pi/nt, math.pi/nf
    delta = [float(Lt)*dth/2 + float(Lp)*dph/2 for Lt, Lp in lips]
    z4 = [0.0]*4
    dS = None; dM = None
    for j in range(4):
        for lado in (0, 1):
            e = [0.0]*4; e[j] = 1.0
            S, M = forma(e, z4, 0.0) if lado == 0 else forma(z4, e, 0.0)
            S, M = np.abs(np.array(S, complex)), np.abs(np.array(M, complex))
            dS = S*delta[j] if dS is None else dS + S*delta[j]
            dM = M*delta[j] if dM is None else dM + M*delta[j]
    S1, M1 = forma(z4, z4, 1.0)                      # padrao do termo mu/xi
    pS, pM = np.abs(np.array(S1, complex)), np.abs(np.array(M1, complex))
    kk = float(K2); aa, bb = (kk + 1/kk)/2, (kk - 1/kk)/2
    Nh = np.diag([1/math.sqrt(2)] + [1.0]*(dim-1))
    fl = fn = 0.0
    mumax = max(float(m) for m in mus)
    for i in range(*linhas):
        lo, hi = i*dph, (i+1)*dph
        cos2 = 0.0 if lo <= math.pi/2 <= hi else min(math.cos(lo)**2, math.cos(hi)**2)
        dxi = mumax*aa/(bb*bb + (aa*aa - bb*bb)*cos2)*dph/2
        Si, Mi = dS + pS*dxi, dM + pM*dxi
        fl = max(fl, np.linalg.norm(Nh @ ((Si + Si.T)/2) @ Nh, 'fro'))
        fn = max(fn, np.linalg.norm(Mi @ Nh, 'fro'))
    return fl*1.001, fn*1.001


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--operador', choices=('K', 'L'), default='K')
    ap.add_argument('--nphi', type=int, default=1024)
    ap.add_argument('--ntheta', type=int, default=2048)
    ap.add_argument('--linhas', default=None, help='faixa i0:i1 das linhas phi')
    ap.add_argument('--sem-rx', action='store_true', help='K: so Gamma2sharp (mu=0)')
    ap.add_argument('--t-lambda', required=False)
    ap.add_argument('--t-norma', required=False)
    ap.add_argument('--out', type=Path)
    ap.add_argument('--combinar', nargs='*', type=Path, help='junta relatorios de faixas')
    a = ap.parse_args()

    if a.combinar:
        rs = [json.loads(p.read_text()) for p in a.combinar]
        cob = sorted(tuple(r['linhas']) for r in rs)
        assert cob[0][0] == 0 and cob[-1][1] == rs[0]['nphi'] and all(cob[i][1] == cob[i+1][0] for i in range(len(cob)-1)), cob
        assert all(r['falhas'] == 0 and r['t_lambda'] == rs[0]['t_lambda'] and r['t_norma'] == rs[0]['t_norma'] for r in rs)
        fl, fn = max(r['folga_lambda'] for r in rs), max(r['folga_norma'] for r in rs)
        corr = rs[0]['correcao']
        Rb = float(Q(rs[0]['t_lambda'])) + fl + corr; Nb = float(Q(rs[0]['t_norma'])) + fn + corr
        print(f"{rs[0]['operador']}: {sum(r['pontos'] for r in rs)} pontos, 0 falhas; folga {fl:.5f}/{fn:.5f}")
        print(f'=> max Re W(-B) <= {Rb:.5f}     ||B|| <= {Nb:.5f}')
        if a.out:
            a.out.write_text(json.dumps(dict(operador=rs[0]['operador'], R_fundo=Rb, norma_fundo=Nb,
                                              faixas=[str(p) for p in a.combinar]), indent=2) + '\n')
        return

    cfg = OPER[a.operador]
    K1, K2, forma, dim = cfg['k1'], cfg['k2'], cfg['forma'], cfg['dim']
    coefs, lips = carrega(K1, K2)
    P, Pc = polinomios(coefs)
    k1, k2 = A(K1), A(K2)
    pi = arb.pi()
    nt, nf = a.ntheta, a.nphi
    i0, i1 = (map(int, a.linhas.split(':')) if a.linhas else (0, nf))
    W = acb_mat(2*NM-1, nt)
    for k in range(nt):
        w = k1*acb(0, 2*pi*arb(fmpq(2*k+1, 2*nt))).exp()
        wi = 1/w; pw = pm = acb(1)
        W[NM-1, k] = acb(1)
        for m in range(1, NM):
            pw *= w; pm *= wi
            W[NM-1+m, k] = pw; W[NM-1-m, k] = pm
    mus = (Q(0),) if a.sem_rx else MUS
    tl, tn = Q(a.t_lambda), Q(a.t_norma)
    t_l, t_n = A(tl), A(tn)
    ns = [arb(1)/arb(2).sqrt()] + [arb(1)]*(dim-1)
    falhas = pontos = 0; fmax_l, fmax_n = -1e9, 0.0
    for i in range(i0, i1):
        z = k2*acb(0, pi*arb(fmpq(2*i+1, 2*nf))).exp()
        V = linha_A(P, Pc, z) * W
        xi = (z + 1/z)/2
        for mu in mus:
            mx = A(mu)/xi
            for k in range(nt):
                S, M = forma([V[r, k] for r in range(4)], [V[r, k] for r in range(4, 8)], mx)
                H = [[-(S[r][c] + S[c][r].conjugate())/2*ns[r]*ns[c] for c in range(dim)] for r in range(dim)]
                T = [[(t_l if r == c else 0) - H[r][c] for c in range(dim)] for r in range(dim)]
                MN = [[M[r][c]*ns[c] for c in range(dim)] for r in range(len(M))]
                G = [[sum((MN[r][c1].conjugate()*MN[r][c2] for r in range(len(M))), acb(0)) for c2 in range(dim)] for c1 in range(dim)]
                T2 = [[(t_n*t_n if r == c else 0) - G[r][c] for c in range(dim)] for r in range(dim)]
                pontos += 1
                if not (pd(T) and pd(T2)):
                    falhas += 1
                if k % 11 == 0:
                    Hf = np.array([[complex(x.real.mid(), x.imag.mid()) for x in row] for row in H])
                    Mf = np.array([[complex(x.real.mid(), x.imag.mid()) for x in row] for row in MN])
                    fmax_l = max(fmax_l, np.linalg.eigvalsh(Hf)[-1]); fmax_n = max(fmax_n, np.linalg.norm(Mf, 2))
        if (i - i0) % max(1, (i1 - i0)//8) == 0:
            print(f'  linha {i} [{i0},{i1}): falhas {falhas}', flush=True)
    fl, fn = folga(forma, dim, lips, nt, nf, (i0, i1), K2, mus)
    corr = 40*float(RT_BACKGROUND_ERROR)
    print(f'{a.operador} linhas [{i0},{i1}): {pontos} pontos, falhas = {falhas}; amostra float {fmax_l:.5f}/{fmax_n:.5f};'
          f' folga {fl:.5f}/{fn:.5f}')
    print(f'=> (faixa) max Re W(-B) <= {float(tl)+fl+corr:.5f}   ||B|| <= {float(tn)+fn+corr:.5f}')
    if a.out:
        a.out.write_text(json.dumps(dict(operador=a.operador, nphi=nf, ntheta=nt, linhas=[i0, i1], pontos=pontos,
                                          falhas=falhas, t_lambda=str(tl), t_norma=str(tn), folga_lambda=fl,
                                          folga_norma=fn, correcao=corr, amostra=[fmax_l, fmax_n]), indent=2) + '\n')
    assert falhas == 0


if __name__ == '__main__':
    main()
