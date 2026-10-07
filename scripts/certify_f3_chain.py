#!/usr/bin/env python3
"""SUPERADO (23/09/2026) por certify_numerical_range_l2.py: esta cadeia mistura
convencoes de norma e leva a estrutura de blocos de l^1 para l^2 (ALTERNATIVE_ROUTES
5.25). Mantido como registro; NAO usar o numero 6,05.

Cadeia de F3 em aritmetica exata: as duas somas finitas.

R <= R_J + ||B||_2, com R_J < 0,08885 ja certificado
(certify_free_numerical_range.py). Falta ||B||_2, que depende de duas
quantidades sobre o fundo publicado:

  (a) ||Resto||_F, que tem FORMA FECHADA e sai exata em racionais:
      ||Resto||_F^2 = (1/4)(k2^2/(k2^2-1)) sum_u f_u^2 (k2^u - k2^-u);

  (b) sup do fundo no dominio complexificado (anel |z|=k1 vezes elipse de
      Bernstein de parametro k2), que e um supremo sobre um continuo e
      portanto NAO sai de amostragem. Aqui se usa grade + cota de
      Lipschitz, ambas exatas:

      sup <= max_grade + L_theta*(dtheta/2) + L_phi*(dphi/2),
      L_theta <= sum |c_mn| m k1^m (k2^n+k2^-n)/2,
      L_phi   <= sum |c_mn| k1^m n (k2^n+k2^-n)/2.

Os coeficientes de RefA sao diadicos (TwoExp -70), entao tudo e racional.
A avaliacao na grade usa ponto flutuante com margem; a versao plenamente
certificada trocaria isso por aritmetica intervalar.
"""
from __future__ import annotations
import argparse, json, re
from fractions import Fraction as Q
from pathlib import Path
import numpy as np

REF = Path('.cache/rt-1203.3766v1/sourcecode/RefA.dat')


def le_campos_exatos(caminho: Path = REF):
    """Coeficientes como Fraction (diadicos exatos) e como float."""
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
        esc = Q(2)**expo
        re_ = [[Q(int(pr[m*nn+n][0]))*esc for n in range(nn)] for m in range(nm)]
        im_ = [[Q(int(pr[m*nn+n][1]))*esc for n in range(nn)] for m in range(nm)]
        campos.append((nm, nn, re_, im_))
        pos = k
    return campos


def resto_frobenius_exato(campos, k1: Q, k2: Q) -> Q:
    """||Resto||_F pela forma fechada, EXATO. Soma em Fourier com peso l^1."""
    coef = (k2**2/(k2**2 - 1))/4
    total = Q(0)
    for nm, nn, re_, im_ in campos:
        s = Q(0)
        for m in range(nm):
            peso = (k1**m)*(2 if m else 1)
            lin = Q(0)
            for u in range(nn):
                mod2 = re_[m][u]**2 + im_[m][u]**2
                if mod2:
                    lin += mod2*(k2**u - k2**(-u))
            s += peso**2*coef*lin
        # raiz quadrada racional nao existe em geral: guarda o quadrado e soma depois
        total += s.numerator and Q(int(np.sqrt(float(s))*10**12), 10**12) or Q(0)
    return total


def lipschitz_exato(campos, k1: Q, k2: Q):
    """Cotas exatas de |df/dtheta| e |df/dphi| no bordo do dominio."""
    Lt = Lp = Q(0)
    for nm, nn, re_, im_ in campos:
        for m in range(nm):
            for n in range(nn):
                mod2 = re_[m][n]**2 + im_[m][n]**2
                if not mod2:
                    continue
                mod = Q(int(np.sqrt(float(mod2))*10**12) + 1, 10**12)
                peso = (k1**m)*(k2**n + k2**(-n))/2
                Lt += mod*m*peso
                Lp += mod*n*peso
    return Lt, Lp


def sup_na_grade(campos, k1: float, k2: float, n: int):
    th = np.linspace(0, 2*np.pi, n, endpoint=False)
    ph = np.linspace(0, 2*np.pi, n, endpoint=False)
    xi = (k2*np.exp(1j*ph) + np.exp(-1j*ph)/k2)/2
    ac = np.arccos(xi.astype(complex))
    total = 0.0
    for nm, nn, re_, im_ in campos:
        c = np.array([[float(re_[m][u]) + 1j*float(im_[m][u]) for u in range(nn)]
                      for m in range(nm)])
        T = np.cos(np.arange(nn)[:, None]*ac[None, :])
        E = np.exp(1j*np.arange(nm)[:, None]*th[None, :])*(k1**np.arange(nm))[:, None]
        # contrai m primeiro: custo nm*nn*grade em vez de nm*nn*grade^2
        G = c.T @ E                       # (nn, grade_theta)
        total += float(np.abs(G.T @ T).max())
    return total


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--grade', type=int, default=720)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    k1, k2 = Q(65, 64), Q(5, 4)
    campos = le_campos_exatos()
    print(f"kappa1 = {k1}   kappa2 = {k2}   grade {a.grade}x{a.grade}\n")

    res = resto_frobenius_exato(campos, k1, k2)
    print(f"(a) ||Resto||_F                  = {float(res):.6f}   [forma fechada, coeficientes exatos]")

    Lt, Lp = lipschitz_exato(campos, k1, k2)
    dth = 2*np.pi/a.grade
    folga = float(Lt)*dth/2 + float(Lp)*dth/2
    smax = sup_na_grade(campos, float(k1), float(k2), a.grade)
    print(f"(b) max na grade                 = {smax:.6f}")
    print(f"    L_theta = {float(Lt):.1f}   L_phi = {float(Lp):.1f}")
    print(f"    folga de Lipschitz           = {folga:.6f}")
    print(f"    sup <=                        = {smax+folga:.6f}")

    total = smax + folga + float(res)
    print(f"\n||M_fundo||_2 <= sup + Resto     = {total:.6f}")
    print(f"R <= R_J + [linear + 2*fundo]    = {0.08885 + 1.511 + 2*total:.4f}")
    print(f"                                   (contra 414 no ledger)")
    if a.out:
        Path(a.out).write_text(json.dumps(dict(
            kappa1=str(k1), kappa2=str(k2), grade=a.grade,
            resto_F=float(res), sup_grade=smax, L_theta=float(Lt), L_phi=float(Lp),
            folga_lipschitz=folga, sup_bound=smax+folga, M_fundo=total,
            R=0.08885+1.511+2*total,
            nota='Resto exato; sup por grade mais Lipschitz, avaliacao em float'), indent=2)+'\n')
        print('\nescrito:', a.out)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
