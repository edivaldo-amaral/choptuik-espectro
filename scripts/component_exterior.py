"""Majorante 4x4 da multiplicacao pelo fundo, preservando componentes.

Somente convolucoes sao majoradas; sinais dos coeficientes do fundo sao
combinados ANTES dos modulos. A norma reponderada e validada racionalmente.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import re

import numpy as np

from background_couplings import PAIR, RT_BACKGROUND_ERROR
from sharp_propagation import SOURCE, table
from structured_exterior import structured_input_tail


def read_fields(path):
    blocks = re.split(r'(?m)^Field\s*$', Path(path).read_text())[1:]
    if len(blocks) != 4:
        raise ValueError('esperados quatro campos')
    fields = []
    for block in blocks:
        header = {key: int(re.search(rf'(?m)^{key}\s+(-?\d+)\s*$', block)[1])
                  for key in ('TwoExp', 'off_m', 'off_n', 'num_m', 'num_n')}
        pairs = list(PAIR.finditer(block))
        if len(pairs) != header['num_m']*header['num_n']:
            raise ValueError('numero incorreto de coeficientes')
        field = {}
        for index, pair in enumerate(pairs):
            i, j = divmod(index, header['num_n'])
            m, n = i+header['off_m'], j+header['off_n']
            for part in range(2):
                value = Q(int(pair[part+1]))*Q(2)**header['TwoExp']
                if value:
                    field[m, n, part] = value
        fields.append(field)
    return fields


def field_expression(expr, fields):
    result = {}
    for field, coefficient in re.findall(r'A\(t,v([1-4]),([+-]?\d+)\)|P\(t\)', expr):
        if not field:
            result = {key: value*(-1)**key[1] for key, value in result.items()}
        else:
            for key, value in fields[int(field)-1].items():
                result[key] = result.get(key, Q(0))+int(coefficient)*value
    return result


def field_norm(field, parity=None, fourier=Q(65, 64), radial=Q(5, 4)):
    return sum((abs(v)*(2-(m == 0))*(2-(n == 0))*fourier**m*radial**n
                for (m, n, _), v in field.items() if parity is None or n % 2 == parity), Q(0))


def component_majorant(fields):
    rows = table('GAMMA2_macro.c', '_G2PAIR_')
    selected = (0, 1, 3, 4)
    result = [[Q(0) for _ in range(4)] for _ in range(4)]
    for out, row in enumerate(selected):
        for inp in range(4):
            columns = [(0, 1, Q(1, 2))] if inp == 0 else [
                (2*inp-1, 0, Q(1)), (2*inp, 1, Q(1))]
            bounds = []
            for k, parity, factor in columns:
                coefficient = field_expression(rows.get((row, k), ''), fields)
                # O(c*v) so retém a paridade complementar de c.
                bounds.append(factor*field_norm(coefficient, 1-parity if out == 0 else None))
            result[out][inp] = max(bounds)
    return result


def gamma1_majorant(mu_max=Q(17, 100)):
    d = Q(40, 9)*mu_max
    return [[Q(0)]*4, [d, 2*d, Q(0), Q(0)],
            [d, Q(0), Q(0), Q(0)], [Q(0), Q(0), Q(0), 2*d]]


def weighted_bound(matrix, weights):
    if len(weights) != len(matrix) or any(w <= 0 for w in weights):
        raise ValueError('pesos positivos exigidos')
    return max(sum(weights[i]*matrix[i][j] for i in range(len(weights)))/weights[j]
               for j in range(len(weights)))


def propose_weights(matrix):
    # Autovetor e apenas uma proposta; weighted_bound valida com racionais.
    array = np.array([[float(v) for v in row] for row in matrix])
    eigenvalues, vectors = np.linalg.eig(array.T)
    index = int(np.argmax(eigenvalues.real))
    v = np.abs(vectors[:, index].real)
    if not np.all(np.isfinite(v)) or np.min(v) <= 0:
        return [Q(1)]*len(matrix)
    v /= np.max(v)
    return [Q(max(1, round(float(x)*4096)), 4096) for x in v]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('background', type=Path)
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    nonlinear = component_majorant(read_fields(args.background))
    linear = gamma1_majorant()
    matrix = [[a+b for a, b in zip(x, y)] for x, y in zip(nonlinear, linear)]
    weights = propose_weights(matrix)
    correction = 6*RT_BACKGROUND_ERROR*max(weights)/min(weights)
    beta = weighted_bound(matrix, weights)+correction
    beta_upper = Q(-((-beta.numerator*1000)//beta.denominator), 1000)
    # Q0 e diagonal nas quatro componentes: seu bound nao muda com esses pesos.
    result = dict(background=str(args.background),
                  background_sha256=hashlib.sha256(args.background.read_bytes()).hexdigest(),
                  gamma2_source_sha256=hashlib.sha256((SOURCE/'GAMMA2_macro.c').read_bytes()).hexdigest(),
                  nonlinear_majorant=[[str(v) for v in row] for row in nonlinear],
                  nonlinear_decimal=[[float(v) for v in row] for row in nonlinear],
                  gamma1_majorant=[[str(v) for v in row] for row in linear],
                  weights=[str(w) for w in weights],
                  unweighted_reference_bound=str(weighted_bound(matrix, [Q(1)]*4)),
                  weighted_beta=str(beta), weighted_beta_decimal=float(beta),
                  weighted_beta_upper=str(beta_upper),
                  background_correction=str(correction),
                  conditional_on_RT_background_ball=True, physical_certificate=False,
                  TT_64x192_decimal=float((beta+Q(5, 4))*structured_input_tail(64, 192)),
                  TT_175x942_upper=str((beta_upper+Q(5, 4))*structured_input_tail(175, 942)))
    output = json.dumps(result, indent=2)+'\n'
    print(output, end='')
    if args.out:
        args.out.write_text(output)


if __name__ == '__main__':
    main()
