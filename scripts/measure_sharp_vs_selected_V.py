#!/usr/bin/env python3
"""DIAGNOSTICO EM PONTO FLUTUANTE -- NAO E PROVA.

Compara, no contorno Gamma_A = [1/8,1] x [-1/4,1/4], a norma da inversa
precondicionada que dimensiona a rota certificada:

    ||V(s)|| = || J (op + s I)^-1 ||_w      (norma de coluna l1 ponderada)

para o operador selecionado L (5 componentes, pesos RT) e para o operador de
constraints K (sistema sharp, 2 componentes, pesos (129/128, 9/8)).

Por que importa para S3: se K for injetivo em toda a regiao fora de um disco em
volta da sua raiz ~0,401, todo modo de L fora desse disco satisfaz as
constraints (identidade K C = M L). O custo de certificar isso escala com
max ||V_K|| no contorno, assim como o de L escala com max ||V_L||.
"""
import sys
from fractions import Fraction as Q
from pathlib import Path
import numpy as np
sys.path.insert(0, 'scripts')
from check_free_preconditioner import read_spectral_matrix, read_addresses, free_matrix as free_L
from build_tile import read_rt_weights
from independent_crown_tiles import free_matrix as free_K
from independent_sharp_matrix import read_sharp

MU = Q(722873400, 2**32)


def selecionado(malha):
    p = Path(f'build/spectrum/rt-A-{malha}.dat')
    dim, raw, expo = read_spectral_matrix(p)
    A = np.array([[float(x) for x in row] for row in raw]) * 2.0**expo
    J = np.array([[float(x) for x in row] for row in free_L(read_addresses(p, dim), MU)])
    w = np.array([float(x) for x in read_rt_weights(p, dim)])
    return A, J, w


def sharp(malha):
    addr, K = read_sharp(Path(f'build/spectrum/rt-sharp-{malha}.dat'))
    Km = np.array([[float(x) for x in row] for row in K])
    J = np.array([[float(x) for x in row] for row in free_K(addr, MU)])
    w = np.array([float((2-(m == 0))*(2-(n == 0))*Q(129, 128)**m*Q(9, 8)**n)
                  for d, p, m, n in addr])
    return Km, J, w


def norma_V(A, J, w, s):
    M = J @ np.linalg.inv(A + s*np.eye(len(A)))
    return float(np.max((np.abs(M) * w[:, None]).sum(axis=0) / w))


def contorno(k=12):
    esq = [complex(0.125, 0.25*t/k) for t in range(k+1)]
    topo = [complex(0.125 + 0.875*t/k, 0.25) for t in range(1, k+1)]
    dir_ = [complex(1.0, 0.25*t/k) for t in range(k)]
    return esq + topo + dir_


def main(malhas):
    pts = contorno()
    for malha in malhas:
        print(f'\n=== {malha} ===')
        for nome, op in (('L (selecionado)', selecionado), ('K (constraints)', sharp)):
            A, J, w = op(malha)
            raizes = sorted((r for r in -np.linalg.eigvals(A)
                             if r.real > 0 and abs(r.imag) <= 0.25 + 1e-9), key=lambda z: z.real)
            vals = [(norma_V(A, J, w, s), s) for s in pts]
            vmax, smax = max(vals)
            v1 = norma_V(A, J, w, 1.0)
            print(f'{nome:<17} dim={len(A):>5}  max ||V|| no contorno = {vmax:10.1f}'
                  f'  em s={smax.real:.3f}{smax.imag:+.3f}i   ||V(1)||={v1:8.1f}')
            print('   raizes em Re>0, |Im|<=1/4:',
                  ', '.join(f'{r.real:.4f}{r.imag:+.4f}i' for r in raizes[:8]) or 'nenhuma')


if __name__ == '__main__':
    main(sys.argv[1:] or ['8x24', '12x36'])
