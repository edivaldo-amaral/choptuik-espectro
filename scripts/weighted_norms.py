#!/usr/bin/env python3
"""Norma ponderada de B e de V com pesos por DOF parametrizados."""
import sys
from fractions import Fraction as Q
from pathlib import Path
sys.path.insert(0, 'scripts')
from check_free_preconditioner import (read_spectral_matrix, read_addresses,
                                       free_matrix, exact_inverse, multiply,
                                       integer_matrix, one_norm, matrix_product)
from build_tile import approximate_inverse

MU = Q(722873400, 2**32)

def enderecos(path, dim):
    out = {}
    active = False
    with open(path, encoding='ascii') as f:
        for line in f:
            if line.strip() == 'BEGIN_ADDRESSES': active = True; continue
            if line.strip() == 'BEGIN_ENTRIES': break
            if active:
                i, comp, part, m, n = map(int, line.split())
                out[i] = (comp, part, m, n)
    return [out[i] for i in range(dim)]

def pesos(addr, k1, k2):
    return [(2-(m == 0))*(2-(n == 0))*k1**m*k2**n for (_, _, m, n) in addr]

def norma_B(malha, k1, k2):
    """Norma de coluna ponderada de B = A - J_mu, medida na matriz exportada."""
    p = Path(f'build/spectrum/rt-A-{malha}.dat')
    dim, raw, expo = read_spectral_matrix(p)
    addr_raw = read_addresses(p, dim)
    addr = enderecos(p, dim)
    w = pesos(addr, k1, k2)
    A = [[Q(x)*Q(2)**expo for x in row] for row in raw]
    J = free_matrix(addr_raw, MU)
    B = [[A[i][j]-J[i][j] for j in range(dim)] for i in range(dim)]
    return max(sum(abs(B[i][j])*w[i] for i in range(dim))/w[j] for j in range(dim))

def norma_V(malha, centro, k1, k2, bits=40):
    p = Path(f'build/spectrum/rt-A-{malha}.dat')
    dim, raw, expo = read_spectral_matrix(p)
    addr_raw = read_addresses(p, dim)
    w = pesos(enderecos(p, dim), k1, k2)
    principal = free_matrix(addr_raw, MU)
    inv_free = exact_inverse(principal)
    pencil = [[Q(x)*Q(2)**expo for x in row] for row in raw]
    for i in range(dim): pencil[i][i] += centro
    rel = multiply(pencil, inv_free)
    wt = [[x*w[i]/w[j] for j, x in enumerate(row)] for i, row in enumerate(rel)]
    num, den = integer_matrix(wt)
    guess = approximate_inverse(num, den, bits)
    prod = matrix_product(guess, num)
    scale = den*2**bits
    res = [[(scale if i == j else 0)-x for j, x in enumerate(row)] for i, row in enumerate(prod)]
    eta = one_norm(res, scale)
    return (one_norm(guess, 2**bits)/(1-eta), float(eta)) if eta < 1 else (None, float(eta))

if __name__ == '__main__':
    b = norma_B('6x18', Q(65,64), Q(5,4))
    print(f"||B_6x18|| com os pesos RT = {float(b):.6f}   (esperado 8.077805)")
