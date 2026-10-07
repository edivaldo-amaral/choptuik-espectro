"""Compara o pencil sharp com as constraints de modos selecionados finitos.

Localizador numerico, nao contagem certificada nem equivalencia de dominios.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eig

from analyze_spectrum import read_matrix, read_constraints


def weights(addresses, fourier=129/128, radial=9/8):
    return np.array([(2-(m == 0))*(2-(n == 0))*fourier**m*radial**n
                     for _, _, m, n in addresses])


def constraint_addresses(path):
    with path.open() as stream:
        for line in stream:
            if line.strip() == 'BEGIN_ADDRESSES':
                break
        result = []
        for line in stream:
            if line.strip() == 'BEGIN_ENTRIES':
                return result
            index, *address = map(int, line.split())
            if index != len(result):
                raise ValueError('enderecos de constraints fora de ordem')
            result.append(tuple(address))
    raise ValueError('enderecos de constraints incompletos')


def compare(selected_path, sharp_path):
    selected, aaddr = read_matrix(selected_path)
    sharp, kaddr = read_matrix(sharp_path, 'CHOPTUIK_RT_SHARP_MATRIX_V1')
    constraints_path = Path(str(selected_path)+'.constraints')
    if constraint_addresses(constraints_path) != kaddr:
        raise ValueError('sharp e constraints devem ter os MESMOS enderecos')
    c0, c1 = read_constraints(constraints_path, len(selected))
    ka, kv = eig(sharp)
    la, lv = eig(selected)
    wa, wk = weights(aaddr), weights(kaddr)
    wa_strong = weights(aaddr, 65/64, 5/4)
    norm = lambda x, w: float(np.sum(np.abs(x)*w))
    modes = []
    for index, s in enumerate(-la):
        if s.real <= 0 or abs(s.imag) > 1e-7:
            continue
        h = lv[:, index]
        c = (c0+s*c1)@h
        residual = (sharp+s*np.eye(len(sharp)))@c
        nearest = int(np.argmin(np.abs(-ka-s)))
        denom = norm(c, wk)
        modes.append(dict(s=float(s.real), constraint_ratio=denom/norm(h, wa),
                          constraint_ratio_strong_to_weak=denom/norm(h, wa_strong),
                          nearest_sharp_s=[float((-ka[nearest]).real), float((-ka[nearest]).imag)],
                          sharp_distance=float(abs(-ka[nearest]-s)),
                          propagated_relative_residual=norm(residual, wk)/denom if denom else None))
    return dict(selected=str(selected_path), sharp=str(sharp_path),
                selected_sha256=hashlib.sha256(selected_path.read_bytes()).hexdigest(),
                sharp_sha256=hashlib.sha256(sharp_path.read_bytes()).hexdigest(),
                selected_dimension=len(selected), sharp_dimension=len(sharp),
                weights=['129/128', '9/8'],
                positive_sharp_roots=[[float(s.real), float(s.imag)] for s in -ka if s.real > 0],
                selected_modes=sorted(modes, key=lambda row: row['s']),
                physical_certificate=False,
                limitation='finite truncation, approximate background, floating point eigenvectors')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('selected', type=Path)
    parser.add_argument('sharp', type=Path)
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    output = json.dumps(compare(args.selected, args.sharp), indent=2)+'\n'
    print(output, end='')
    if args.out:
        args.out.write_text(output)
