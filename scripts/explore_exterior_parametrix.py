"""Diagnostico FINITO de correcoes de Neumann da coroa.

Usa floats; os pontos amostrados nao constituem um bound uniforme.
"""
import argparse
import json
from fractions import Fraction as Q
from pathlib import Path

import numpy as np
from scipy.linalg import eigvals

from analyze_spectrum import read_matrix
from check_free_preconditioner import free_matrix


def explore(path, core_m, core_n):
    a, addresses = read_matrix(path)
    principal = np.array(free_matrix(addresses, Q(722873400, 2**32)), dtype=float)
    q = np.linalg.inv(principal)
    eta = [916, 4096, 243, 1233]
    weights = np.array([eta[d]*(2-(m == 0))*(2-(n == 0))*(65/64)**m*(5/4)**n
                        for d, _, m, n in addresses])
    shell = [i for i, (_, _, m, n) in enumerate(addresses) if m >= core_m or n >= core_n]
    if not shell or len(shell) == len(addresses):
        raise ValueError('nucleo e coroa devem ser nao vazios')
    rows = []
    for s in (1+0j, 1+0.25j, 0.125+0.25j, 0.125+0j):
        h = (a+s*np.eye(len(a)))@q
        h = h*weights[:, None]/weights[None, :]
        hss = h[np.ix_(shell, shell)]
        qw = q*weights[:, None]/weights[None, :]
        qss = qw[np.ix_(shell, shell)]
        variants = []
        for damping in (1.0, 0.5, 0.25):
            e = np.eye(len(shell))-damping*hss
            power = np.eye(len(shell), dtype=complex)
            parametrix = np.zeros_like(power)
            norms = {}
            for k in range(1, 65):
                if k <= 32:
                    parametrix += damping*power
                power = power@e
                if k in (1, 2, 4, 8, 16, 32, 64):
                    norms[k] = float(np.linalg.norm(power, 1))
            variation = float(np.linalg.norm(parametrix@qss, 1))
            variants.append(dict(damping=damping,
                                 spectral_radius=float(max(abs(eigvals(e)))), power_norms=norms,
                                 parametrix32_norm=float(np.linalg.norm(parametrix, 1)),
                                 fixed_tile_variation=variation,
                                 proposed_tile_radius=(1-norms[32])/(2*variation)
                                 if norms[32] < 1 and variation else None))
        rows.append(dict(s=[s.real, s.imag], variants=variants))
    return dict(matrix=str(path), core=[core_m, core_n], shell_dimension=len(shell),
                samples=rows, infinite_tail_bounds=None, physical_certificate=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('matrix', type=Path)
    parser.add_argument('--core-m', type=int, default=4)
    parser.add_argument('--core-n', type=int, default=12)
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    output = json.dumps(explore(args.matrix, args.core_m, args.core_n), indent=2)+'\n'
    print(output, end='')
    if args.out:
        args.out.write_text(output)
