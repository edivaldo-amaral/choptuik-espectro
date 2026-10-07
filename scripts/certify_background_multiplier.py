#!/usr/bin/env python3
"""Bound certificavel para a multiplicacao pelo fundo na realizacao l^2.

Ultimo elo da cadeia de F3 (ALTERNATIVE_ROUTES 5.9-5.20):

    R <= R_J + ||B||_2,   R_J < 0,08885 certificado (certify_free_numerical_range.py)

e ||B||_2 depende da norma do operador de multiplicacao pelo fundo. Na base de
Chebyshev com peso rho^n, a formula do produto
T_j T_n = (T_{j+n} + T_{|j-n|})/2 da a decomposicao EXATA

    M_f = Toeplitz + Resto,
    Toeplitz[i,j] = f_{|i-j|}/2 + delta_ij f_0/2      (simbolo = f na elipse)
    Resto[i,j]    = f_{i+j}/2   - delta_i0 f_j/2      (primeira linha NULA)

Logo  ||M_f|| <= sup_elipse|f| + ||Resto||, e para o Resto vale
||.||_2 <= ||.||_F com forma fechada

    ||Resto||_F^2 = (1/4) (rho^2/(rho^2-1)) sum_u f_u^2 (rho^u - rho^-u).

As duas parcelas sao somas finitas sobre os coeficientes publicados de RefA.

DIAGNOSTICO em ponto flutuante. A versao certificada roda a mesma formula em
racionais; nada aqui depende de convergencia ou extrapolacao.
"""
from __future__ import annotations
import argparse, json, re
from pathlib import Path
import numpy as np

K1_PADRAO, K2_PADRAO = 65/64, 5/4
REF = Path('.cache/rt-1203.3766v1/sourcecode/RefA.dat')


def le_campos(caminho: Path = REF):
    txt = caminho.read_text()
    pos, campos = 0, []
    while len(campos) < 4:
        i = txt.find('num_m', pos)
        if i < 0:
            break
        cab = txt[max(0, i-200):i+120]
        expo = int(re.findall(r'TwoExp\s+(-?\d+)', cab)[-1])
        nm = int(re.search(r'num_m\s+(\d+)', txt[i:i+40]).group(1))
        nn = int(re.search(r'num_n\s+(\d+)', txt[i:i+80]).group(1))
        j = txt.index('(((', i)
        d, k = 0, j
        while True:
            if txt[k] == '(':
                d += 1
            elif txt[k] == ')':
                d -= 1
                if d == 0:
                    break
            k += 1
        pr = re.findall(r'\((-?\d+),(-?\d+)\)', txt[j:k+1])
        campos.append(np.array([complex(int(a), int(b)) for a, b in pr]
                               ).reshape(nm, nn)*float(2.0**expo))
        pos = k
    return campos


def sup_dominio(campos, k1=K1_PADRAO, k2=K2_PADRAO, n=241):
    """sup de |campo| no anel |z|=k1 vezes a elipse de Bernstein de parametro k2."""
    th = np.linspace(0, 2*np.pi, n)
    ph = np.linspace(0, 2*np.pi, n)
    xi = (k2*np.exp(1j*ph) + np.exp(-1j*ph)/k2)/2
    ac = np.arccos(xi.astype(complex))
    total = 0.0
    for c in campos:
        M, N = c.shape
        T = np.cos(np.arange(N)[:, None]*ac[None, :])
        E = np.exp(1j*np.arange(M)[:, None]*th[None, :])*(k1**np.arange(M))[:, None]
        total += float(np.abs(np.einsum('mn,mt,np->tp', c, E, T)).max())
    return total


def resto_frobenius(campos, k1=K1_PADRAO, k2=K2_PADRAO):
    """Forma fechada de ||Resto||_F, somada em Fourier com peso l^1 (conservador)."""
    coef = (k2**2/(k2**2 - 1))/4
    total = 0.0
    for c in campos:
        M, N = c.shape
        s = 0.0
        for m in range(M):
            peso = (k1**m)*(2 if m else 1)
            s += peso**2*coef*sum(abs(c[m, u])**2*(k2**u - k2**(-u)) for u in range(N))
        total += float(np.sqrt(s))
    return total


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--kappa1', type=float, default=K1_PADRAO)
    ap.add_argument('--kappa2', type=float, default=K2_PADRAO)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    campos = le_campos()
    sup = sup_dominio(campos, a.kappa1, a.kappa2)
    res = resto_frobenius(campos, a.kappa1, a.kappa2)
    l1 = sum(float(np.sum(np.abs(c)*(a.kappa1**np.arange(c.shape[0]))[:, None]
                          * (a.kappa2**np.arange(c.shape[1]))[None, :]
                          * (2-(np.arange(c.shape[0]) == 0))[:, None]
                          * (2-(np.arange(c.shape[1]) == 0))[None, :])) for c in campos)
    print(f"kappa1 = {a.kappa1:.6f}   kappa2 = {a.kappa2:.6f}\n")
    print(f"  sup no dominio complexificado : {sup:.6f}")
    print(f"  ||Resto||_F (forma fechada)   : {res:.6f}")
    print(f"  ---------------------------------------")
    print(f"  ||M_fundo||_2 <=               : {sup+res:.6f}")
    print(f"  norma l^1 ponderada (hoje)    : {l1:.6f}")
    print(f"  ganho                          : {l1/(sup+res):.2f}x")
    rep = dict(kappa1=a.kappa1, kappa2=a.kappa2, sup=sup, resto_F=res,
               bound_l2=sup+res, norma_l1=l1, ganho=l1/(sup+res),
               nota='diagnostico; as duas parcelas sao somas finitas certificaveis')
    if a.out:
        Path(a.out).write_text(json.dumps(rep, indent=2)+'\n')
        print('\nescrito:', a.out)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
