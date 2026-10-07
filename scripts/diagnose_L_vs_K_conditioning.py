#!/usr/bin/env python3
"""DIAGNOSTICO EM PONTO FLUTUANTE -- NAO E PROVA.

Por que ||V_L|| e ~130x maior que ||V_K||? O operador ou a norma?
Mede V(s) = J (op + s)^-1 para L e K em varias normas comuns:
  propria   L: pesos RT do arquivo (com eta); K: (129/128, 9/8)
  forte     (65/64, 5/4), sem pesos de componente
  fraca     (129/128, 9/8), sem pesos de componente
  l1        sem peso algum
  l2        norma espectral
"""
import sys
from fractions import Fraction as Q
from pathlib import Path
import numpy as np
sys.path.insert(0, 'scripts')
from measure_sharp_vs_selected_V import selecionado, sharp
from weighted_norms import enderecos
from independent_sharp_matrix import read_sharp
from check_free_preconditioner import read_spectral_matrix


def pesos(addr, k1, k2):
    return np.array([float((2-(m == 0))*(2-(n == 0))*Q(k1)**m*Q(k2)**n) for _, _, m, n in addr])


def col_norm(M, w):
    return float(np.max((np.abs(M)*w[:, None]).sum(axis=0)/w))


def enderecos_L(malha):
    p = Path(f'build/spectrum/rt-A-{malha}.dat')
    dim, _, _ = read_spectral_matrix(p)
    return enderecos(p, dim)


def enderecos_K(malha):
    return read_sharp(Path(f'build/spectrum/rt-sharp-{malha}.dat'))[0]


def medir(malha, pontos):
    casos = (('L', selecionado, enderecos_L), ('K', sharp, enderecos_K))
    for nome, op, end in casos:
        A, J, w_prop = op(malha)
        addr = end(malha)
        normas = {'propria': w_prop,
                  'forte': pesos(addr, '65/64', '5/4'),
                  'fraca': pesos(addr, '129/128', '9/8'),
                  'l1': np.ones(len(A))}
        for s in pontos:
            M = J @ np.linalg.inv(A + s*np.eye(len(A)))
            cols = {k: col_norm(M, v) for k, v in normas.items()}
            cols['l2'] = float(np.linalg.norm(M, 2))
            print(f'{malha:>6} {nome} s={s:<6} ' +
                  '  '.join(f'{k}={v:9.2f}' for k, v in cols.items()))


if __name__ == '__main__':
    for malha in sys.argv[1:] or ['8x24', '12x36']:
        medir(malha, [1.0, 0.125])
