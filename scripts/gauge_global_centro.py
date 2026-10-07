#!/usr/bin/env python3
"""Lema G (GAUGE_GLOBAL.md), fatos certificados:
  F1:  omega4*(tau, 0) nao e identicamente nula (centro; injetividade do gauge, Lema G);
  F1': omega4*(tau, -1) nao e identicamente nula (omega4^+ = -P omega4 no cone xi = 1; Lema G').

Na base de RT, omega(tau, xi) = sum_{m,n in Z} v_mn exp(i m tau/2 + i n theta), theta = arccos xi. Em xi = 0,
theta = pi/2, e o coeficiente de Fourier m = 1 de omega4(., 0) e
    c_1 = sum_{n in Z} v_{1n} i^n = v_10 + 2 sum_{n >= 1} v_1n cos(n pi/2)      (v_{1,-n} = v_{1n}).
Bola de RT: omega* = RefA + RefB + Corr, ||RefB|| <= 2^-25 [dkjfhk2], ||Corr|| <= 2^-277 [kdkkkskkskksksks], na
norma sum k1^|m| k2^|n| ||v_mn||_C com pesos >= 1 e ||z||_C >= |z|. Logo |c_1(omega*) - c_1(RefA)| <= 2^-25 + 2^-277.
Em xi = -1, theta = pi: c_1(-1) = sum_n v_1n (-1)^n = v_10 + 2 sum_{n >= 1} (-1)^n v_1n, com a mesma cota.
Confere tambem que o campo de indice 3 so tem m impar (e omega4, TAP). Tudo em racionais exatos."""
from __future__ import annotations

from fractions import Fraction as Fr
from pathlib import Path
import json
import math
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gauge_refB import ler

EPS = Fr(1, 2**25) + Fr(1, 2**277)


def main():
    campos, mu = ler('.cache/rt-1203.3766v1/sourcecode/RefA.dat')
    print(f'{len(campos)} campos; mu_RefA = {float(mu):.15f}')
    h, pares = campos[3]                                   # omega4 (SSPACE: omega1, omega2, omega3, omega4)
    esc = Fr(2)**h['TwoExp']
    print(f"omega4: off_m {h['off_m']}, num_m {h['num_m']}, off_n {h['off_n']}, num_n {h['num_n']}")
    for ii in range(h['num_m']):                           # omega4 e TAP: so m impar
        if (h['off_m'] + ii) % 2 == 0:
            assert all(int(a) == 0 and int(b) == 0 for a, b in pares[ii*h['num_n']:(ii + 1)*h['num_n']]), 'campo 3 com m par'
    print('campo 3: so m impar (omega4, TAP) - conferido')
    i = 1 - h['off_m']
    assert 0 <= i < h['num_m']
    cre = Fr(0); cim = Fr(0)
    cos_n = {0: 1, 1: 0, 2: -1, 3: 0}
    for k in range(h['num_n']):
        n = h['off_n'] + k
        re, im = pares[i*h['num_n'] + k]
        fat = cos_n[n % 4]*(1 if n == 0 else 2)
        cre += fat*int(re)*esc; cim += fat*int(im)*esc
    lim = max(abs(cre), abs(cim))                          # <= |c_1|
    print(f'c_1(RefA4(., 0)) = {float(cre):.12f} + {float(cim):.12f} i;  max(|Re|, |Im|) = {float(lim):.12f}')
    print(f'eps = 2^-25 + 2^-277 = {float(EPS):.3e}')
    ok = lim > EPS
    print(('F1 CERTIFICADO' if ok else 'F1 FALHA') + f': |c_1(omega4*(., 0))| >= {float(lim - EPS):.12f} > 0')
    dre = Fr(0); dim = Fr(0)
    for k in range(h['num_n']):
        n = h['off_n'] + k
        re, im = pares[i*h['num_n'] + k]
        fat = (-1)**n*(1 if n == 0 else 2)
        dre += fat*int(re)*esc; dim += fat*int(im)*esc
    lim2 = max(abs(dre), abs(dim))
    ok2 = lim2 > EPS
    print(f'c_1(RefA4(., -1)) = {float(dre):.12f} + {float(dim):.12f} i')
    print(("F1' CERTIFICADO" if ok2 else "F1' FALHA") + f': |c_1(omega4*(., -1))| >= {float(lim2 - EPS):.12f} > 0')
    # conferencia (nao faz parte da prova): omega4(tau, 0) de RefA numa malha, trocas de sinal em [0, 4 pi)
    vals = []
    for j in range(4000):
        t = 4*math.pi*j/4000; s = 0.0
        for ii in range(h['num_m']):
            m = h['off_m'] + ii; wm = 1 if m == 0 else 2
            cm = 0j
            for k in range(h['num_n']):
                n = h['off_n'] + k
                re, im = pares[ii*h['num_n'] + k]
                cm += cos_n[n % 4]*(1 if n == 0 else 2)*complex(int(re), int(im))*float(esc)
            s += wm*(cm*complex(math.cos(m*t/2), math.sin(m*t/2))).real
        vals.append(s)
    trocas = sum(1 for a, b in zip(vals, vals[1:] + vals[:1]) if a*b < 0)
    print(f'conferencia: omega4_RefA(tau, 0) em [0, 4pi): max |.| = {max(map(abs, vals)):.4f}, trocas de sinal = {trocas}')
    Path('build/t2').mkdir(parents=True, exist_ok=True)
    Path('build/t2/gauge_global_centro.json').write_text(json.dumps(dict(
        c1_re=str(cre), c1_im=str(cim), eps=str(EPS), cota_inferior=str(lim - EPS), certificado=ok,
        c1_cone_re=str(dre), c1_cone_im=str(dim), cota_inferior_cone=str(lim2 - EPS), certificado_cone=ok2,
        campo3_so_m_impar=True,
        conferencia_max=max(map(abs, vals)), conferencia_trocas=trocas), indent=1) + '\n')


if __name__ == '__main__':
    main()
