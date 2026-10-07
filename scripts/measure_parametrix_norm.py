#!/usr/bin/env python3
"""Mede ||V|| (finite_inverse_bound) aproximando da raiz, em varias malhas."""
import sys, time, json
from fractions import Fraction as Q
from pathlib import Path
sys.path.insert(0, 'scripts')
from check_free_preconditioner import (read_spectral_matrix, read_addresses,
                                       free_matrix, exact_inverse, multiply,
                                       integer_matrix, one_norm)
from build_tile import approximate_inverse, read_rt_weights

def matrix_product(a, b):
    from check_free_preconditioner import matrix_product as mp
    return mp(a, b)

MU = Q(722873400, 2**32)
BITS = 40

def norma_V(malha, centro, base=Q(0)):
    p = Path(f'build/spectrum/rt-A-{malha}.dat')
    dim, raw, expo = read_spectral_matrix(p)
    addr = read_addresses(p, dim)
    w = read_rt_weights(p, dim)
    principal = free_matrix(addr, MU)
    for i in range(dim): principal[i][i] += base
    inv_free = exact_inverse(principal)
    pencil = [[Q(x)*Q(2)**expo for x in row] for row in raw]
    for i in range(dim): pencil[i][i] += centro
    rel = multiply(pencil, inv_free)
    wt = [[x*w[i]/w[j] for j, x in enumerate(row)] for i, row in enumerate(rel)]
    num, den = integer_matrix(wt)
    guess = approximate_inverse(num, den, BITS)
    prod = matrix_product(guess, num)
    scale = den*2**BITS
    res = [[(scale if i == j else 0)-x for j, x in enumerate(row)] for i, row in enumerate(prod)]
    eta = one_norm(res, scale)
    if eta >= 1: return dim, None, float(eta)
    return dim, one_norm(guess, 2**BITS)/(1-eta), float(eta)

RAIZES = {'4x12': 0.231793, '6x18': 0.144011, '8x24': 0.166389, '10x30': 0.168311}
CENTROS = [('1', Q(1)), ('1/2', Q(1,2)), ('1/4', Q(1,4)),
           ('1/8 = borda', Q(1,8)), ('0,15', Q(3,20)), ('0,168311 = raiz 10x30', Q(168311,10**6))]

if __name__ == '__main__':
    for malha in sys.argv[1:]:
        print(f"\n=== malha {malha} (raiz mais proxima: {RAIZES.get(malha,'?')}) ===")
        print(f"{'centro s':<24} {'dim':>5} {'||V||':>16} {'dist a raiz':>12} {'defeito eta':>12}")
        for nome, c in CENTROS:
            t = time.time()
            dim, V, eta = norma_V(malha, c)
            d = abs(float(c)-RAIZES.get(malha, float('nan')))
            vs = f"{float(V):>16.2f}" if V is not None else f"{'DIVERGIU':>16}"
            print(f"{nome:<24} {dim:>5} {vs} {d:>12.4f} {eta:>12.2e}   ({time.time()-t:.0f}s)")
