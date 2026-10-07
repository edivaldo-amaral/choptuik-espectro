#!/usr/bin/env python3
"""Cota RIGOROSA e fina da norma de multiplicacao por um campo em l^2(W) simetrico.

Convencao (ver certify_numerical_range_l2.py): coeficientes simetricos, peso
W(m,n) = kappa1^|m| kappa2^|n|. Em Chebyshev, para uma funcao radial g com
coeficientes simetricos g_u (u em Z, g_-u = g_u), nas coordenadas isometricas
y_n = sqrt(c_n) kappa^n x_n a multiplicacao e

    A = T + (H + E),
    T_nm = g_{n-m} kappa^{n-m}                       (Toeplitz em N)
    H_nm = g_{n+m} kappa^{n-m},   n,m >= 1           (Hankel)
    (H+E)_0m = (sqrt2-1) kappa^-m g_m,  (H+E)_n0 = (sqrt2-1) kappa^n g_n.

||T|| <= sup do simbolo = sup de g na elipse de Bernstein E_kappa; e
||H+E||_F^2 = sum_{u>=2} |g_u|^2 (k^2u - k^(4-2u))/(k^4-1)
            + (sqrt2-1)^2 sum_{u>=1} |g_u|^2 (k^2u + k^-2u).

Em Fourier, desigualdade TRIANGULAR (nao media): o shift ponderado por d tem
norma <= kappa1^|d|, logo ||M_f|| <= sum_d c_d kappa1^d ||M_rad[f_d]||.

O sup e calculado com bolas de Arb (python-flint) sobre INTERVALOS de angulo:
cada bola contem todos os valores do simbolo no intervalo. Refinamento
adaptativo ate a cota ficar a 1% do maximo nos centros. Nada de ponto
flutuante no que se afirma.
"""
from __future__ import annotations

from fractions import Fraction as Q
import math

import flint
from flint import arb, acb, acb_poly, fmpq

SQRT2M1_SQ_UP = Q(17158, 100000)          # (sqrt2-1)^2 = 0,171572875... <= 0,17158


def para_arb(x: Q) -> arb:
    return arb(fmpq(x.numerator, x.denominator))


def arb_para_fracao_acima(x: arb) -> Q:
    u = x.upper()
    man, exp = u.man_exp()
    return Q(int(man)) * (Q(2)**int(exp))


def raiz_acima(x: Q) -> Q:
    r = Q(math.isqrt(x.numerator * 10**24 // x.denominator) + 1, 10**12)
    while r*r < x:
        r += Q(1, 10**12)
    return r


def modos_fourier(campo: dict, paridade=None):
    """campo[(m,n,part)] -> {m: {n: (re, im)}} com filtro opcional de paridade em n."""
    out = {}
    for (m, n, part), v in campo.items():
        if paridade is not None and n % 2 != paridade:
            continue
        d = out.setdefault(m, {}).setdefault(n, [Q(0), Q(0)])
        d[part] += v
    return out


def sup_elipse(g: dict, kappa: Q, alvo_rel: Q = Q(1, 100), n0: int = 512, max_iter: int = 14):
    """Cota superior rigorosa de sup_{phi} |sum_u g_u kappa^u e^{i u phi}|, g simetrico."""
    if not g:
        return Q(0)
    U = max(g)
    coef = [acb(para_arb(g[u][0]), para_arb(g[u][1])) if u in g else acb(0) for u in range(U+1)]
    P = acb_poly(coef)
    g0 = coef[0]
    k = para_arb(kappa)
    dois_pi = 2*arb.pi()

    def simb(phi: arb) -> acb:
        z = k * acb(0, phi).exp()
        return P(z) + P(1/z) - g0

    # maximo nos centros (bolas pontuais: enclausuramento apertado) -> alvo
    centros = [simb(dois_pi*arb(fmpq(2*i+1, 2*n0))).abs_upper() for i in range(n0)]
    cmax = max(centros)
    limite = cmax * (1 + para_arb(alvo_rel))
    pend = [(fmpq(i, n0), fmpq(i+1, n0)) for i in range(n0)]
    cota = arb(0)
    for _ in range(max_iter):
        prox = []
        for a, b in pend:
            meio = (a + b)/2
            phi = dois_pi*arb(meio, (b - a)/2)          # bola cobrindo [a,b]*2pi
            v = simb(phi).abs_upper()
            if v < limite:
                cota = arb.max(cota, v) if hasattr(arb, 'max') else (v if v > cota else cota)
            else:
                prox += [(a, meio), (meio, b)]
        pend = prox
        if not pend:
            break
    for a, b in pend:                                     # restos: aceita a cota que houver
        phi = dois_pi*arb((a+b)/2, (b-a)/2)
        v = simb(phi).abs_upper()
        cota = v if v > cota else cota
    return arb_para_fracao_acima(cota)


def frobenius_resto(g: dict, kappa: Q) -> Q:
    """||H+E||_F, exato em racionais (cota superior da raiz)."""
    k2, k4 = kappa**2, kappa**4
    s = Q(0)
    for u, (re, im) in g.items():
        a2 = re*re + im*im
        if u >= 2:
            s += a2*(k2**u - k4/k2**u)/(k4 - 1)
        if u >= 1:
            s += SQRT2M1_SQ_UP*a2*(k2**u + 1/k2**u)
    return raiz_acima(s) if s else Q(0)


def norma_multiplicacao(campo: dict, kappa1: Q, kappa2: Q, paridade=None):
    """Cota rigorosa de ||M_campo|| em l^2(W) simetrico."""
    total = Q(0)
    for m, g in modos_fourier(campo, paridade).items():
        g = {n: v for n, v in g.items() if v[0] or v[1]}
        if not g:
            continue
        rad = sup_elipse(g, kappa2) + frobenius_resto(g, kappa2)
        total += (1 if m == 0 else 2) * kappa1**m * rad
    return total


if __name__ == '__main__':
    # autoteste: cota >= norma verdadeira numa truncacao; cota <= norma de campo
    import sys
    from pathlib import Path
    import numpy as np
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from component_exterior import read_fields, field_norm
    campos = read_fields(Path('.cache/rt-1203.3766v1/sourcecode/RefA.dat'))
    k1, k2 = Q(129, 128), Q(9, 8)
    for i, c in enumerate(campos):
        # so o modo de Fourier m=0, para comparar com a matriz radial verdadeira
        c0 = {key: v for key, v in c.items() if key[0] == 0}
        cota = norma_multiplicacao(c0, k1, k2)
        fn = field_norm(c0, fourier=k1, radial=k2)
        g = {n: complex(float(v[0]), float(v[1])) for n, v in modos_fourier(c0)[0].items()}
        N = 400
        A = np.zeros((N, N), complex)
        for n in range(N):
            for mm in range(N):
                cn = 1 if n == 0 else 2; cm = 1 if mm == 0 else 2
                val = g.get(abs(n-mm), 0) + (g.get(n+mm, 0) if mm > 0 else 0)
                A[n, mm] = math.sqrt(cn/cm) * float(k2)**(n-mm) * val
        verd = np.linalg.norm(A, 2)
        print(f'campo {i+1}, m=0: verdadeira(trunc 400) = {verd:.5f}  <=  cota = {float(cota):.5f}'
              f'  <=  norma de campo = {float(fn):.5f}   {verd <= float(cota) <= float(fn)}')
