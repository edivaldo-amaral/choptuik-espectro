"""Auditoria racional de tiles da COROA FINITA de H=(A+sI)J^-1.

A inversa livre aproximada e validada por residuo exato; evita formar seu
denominador racional global gigante. Inclui opcionalmente o erro publicado
do fundo exato; nunca inclui o exterior infinito.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from build_contour import read_spectral_matrix
from build_tile import read_rt_weights
from check_free_preconditioner import free_matrix, integer_matrix, read_addresses


def ceil_grid(value, bits=40):
    scaled = value*(1 << bits)
    return Q(-(-scaled.numerator//scaled.denominator), 1 << bits)


def nearest_integer(value):
    return (2*value.numerator+value.denominator)//(2*value.denominator)


def norm_real(matrix, denominator):
    return Q(max(sum(abs(int(x)) for x in matrix[:, j]) for j in range(matrix.shape[1])), denominator)


def norm_complex(real, imag, denominator):
    return Q(max(sum(abs(int(a))+abs(int(b)) for a, b in zip(real[:, j], imag[:, j]))
                 for j in range(real.shape[1])), denominator)


def weighted_norm(matrix, denominator, weights):
    return max(sum(abs(int(matrix[i, j]))*weights[i] for i in range(len(weights)))
               / (denominator*weights[j]) for j in range(len(weights)))


def sparse_right_product(left, right):
    size = len(left)
    result = [[0]*size for _ in range(size)]
    for k, row in enumerate(right):
        nonzero = [(j, x) for j, x in enumerate(row) if x]
        for i in range(size):
            a = left[i][k]
            if a:
                for j, b in nonzero:
                    result[i][j] += a*b
    return np.array(result, dtype=object)


def exact_product(left, right):
    """Produto inteiro por limbs de 15 bits; cada dot int64 nao pode overflowar."""
    left, right = np.asarray(left, dtype=object), np.asarray(right, dtype=object)
    if left.ndim != 2 or right.ndim != 2 or left.shape[1] != right.shape[0]:
        raise ValueError('dimensoes incompativeis')
    base_bits, mask = 15, (1 << 15)-1
    if left.shape[1]*mask*mask >= 2**63:
        raise ValueError('acumulador int64 insuficiente')
    def limbs(matrix):
        count = (max((abs(int(v)).bit_length() for v in matrix.flat), default=0)+base_bits-1)//base_bits
        return [np.array([[(1 if int(v) >= 0 else -1)*((abs(int(v)) >> (base_bits*k)) & mask)
                           for v in row] for row in matrix], dtype=np.int64) for k in range(count)]
    result = np.zeros((left.shape[0], right.shape[1]), dtype=object)
    # A mesma garantia int64 vale para produtos esparsos: cada soma tem
    # no maximo inner_dimension termos, cada um <= mask^2 em modulo.
    from scipy.sparse import csr_matrix
    def prepared(matrix, transpose=False):
        values = []
        for a in limbs(matrix):
            sparse = csr_matrix(a.T if transpose else a) if np.count_nonzero(a)*8 <= a.size else None
            values.append((a, sparse))
        return values
    right_parts = prepared(right, transpose=True)
    for i, (a, sa) in enumerate(prepared(left)):
        for j, (b, sb_transpose) in enumerate(right_parts):
            product = sa@b if sa is not None else (
                (sb_transpose@a.T).T if sb_transpose is not None else a@b)
            result += product.astype(object)*(1 << (base_bits*(i+j)))
    return result


def prepare(path, core_m, core_n, bits=48):
    n, raw, exponent = read_spectral_matrix(path)
    addresses = read_addresses(path, n)
    shell = [i for i, (_, _, m, cn) in enumerate(addresses) if m >= core_m or cn >= core_n]
    if not shell or len(shell) == n:
        raise ValueError('nucleo e coroa devem ser nao vazios')
    eta = [916, 4096, 243, 1233]
    weights = [w*eta[a[0]] for w, a in zip(read_rt_weights(path, n), addresses)]
    principal = free_matrix(addresses, Q(722873400, 2**32))
    jn, jd = integer_matrix(principal)
    grid = 1 << bits
    proposed = np.linalg.inv(np.array(principal, dtype=float))
    qn = [[round(float(x)*grid) for x in row] for row in proposed]
    residual = np.eye(n, dtype=object)*(jd*grid)-sparse_right_product(jn, qn)
    residual_norm = weighted_norm(residual, jd*grid, weights)
    proposed_norm = weighted_norm(np.array(qn, dtype=object), grid, weights)
    # Verifica o bound 108 FINITAMENTE, sem depender do lema de EDP.
    if residual_norm >= 1 or proposed_norm/(1-residual_norm) > 108:
        raise ValueError('a auditoria racional nao justificou ||J^-1||<=108')
    qerror = ceil_grid(108*residual_norm)
    an, ad = (raw, 1 << -exponent) if exponent <= 0 else (
        [[v*(1 << exponent) for v in row] for row in raw], 1)
    anorm = ceil_grid(weighted_norm(np.array(an, dtype=object), ad, weights))
    aq = sparse_right_product(an, qn)
    hround, qround = [], []
    for i in shell:
        hround.append([nearest_integer(Q(int(aq[i, j]), ad)*weights[i]/weights[j]) for j in shell])
        qround.append([nearest_integer(Q(qn[i][j])*weights[i]/weights[j]) for j in shell])
    # Entradas acima representam operadores multiplicados por grid.
    rounding = Q(len(shell), 2*grid)
    return dict(h=np.array(hround, dtype=object), q=np.array(qround, dtype=object),
                grid=grid, qerror=qerror, anorm=anorm, rounding=rounding,
                inverse_norm_bound=ceil_grid(proposed_norm/(1-residual_norm)),
                size=len(shell), matrix_dimension=n)


def tile(data, real, imag):
    real, imag = Q(real), Q(imag)
    scale = math.lcm(real.denominator, imag.denominator)
    grid, h, q = data['grid'], data['h'], data['q']
    hr = scale*h+int(real*scale)*q
    hi = int(imag*scale)*q
    matrix = np.array(hr, dtype=float)/(scale*grid)+1j*np.array(hi, dtype=float)/(scale*grid)
    e = np.eye(data['size'])-matrix/2
    power = np.eye(data['size'], dtype=complex)
    v = np.zeros_like(power)
    for _ in range(32):
        v += power/2
        power = power@e
    vr = np.array([[round(float(x.real)*grid) for x in row] for row in v], dtype=object)
    vi = np.array([[round(float(x.imag)*grid) for x in row] for row in v], dtype=object)
    rr = np.eye(data['size'], dtype=object)*(scale*grid*grid)-(exact_product(vr, hr)-exact_product(vi, hi))
    ri = -(exact_product(vr, hi)+exact_product(vi, hr))
    vnorm = norm_complex(vr, vi, grid)
    # |s|<=|Re s|+|Im s|, rationnel et conservateur.
    distance = abs(real)+abs(imag)
    herror = ((data['anorm']+distance)*data['qerror']+(1+distance)*data['rounding']
              +data.get('background_H_error', Q(0)))
    a = ceil_grid(norm_complex(rr, ri, scale*grid*grid)+vnorm*herror)
    variation = ceil_grid(norm_complex(exact_product(vr, q), exact_product(vi, q), grid*grid)
                          +vnorm*(data['qerror']+data['rounding']))
    radius = Q(0)
    if a < 1:
        if not variation:
            raise ValueError('variacao zero inesperada')
        limit = (1-a)/(2*variation)
        radius = Q(1)
        while radius > limit:
            radius /= 2
    return dict(center=[str(real), str(imag)], center_defect=str(a),
                parametrix_norm=str(ceil_grid(vnorm)),
                variation=str(variation), radius=str(radius),
                uniform_defect=str(a+radius*variation),
                accepted=a+radius*variation < 1,
                physical_certificate=False)


def cover_upper_boundary(data):
    """Cobre tres lados da meia fronteira; a outra metade vem por conjugacao."""
    current = (Q(1, 8), Q(0))
    tiles, segments = [], []
    for target in ((Q(1, 8), Q(1, 4)), (Q(1), Q(1, 4)), (Q(1), Q(0))):
        while current != target:
            result = tile(data, *current)
            tiles.append(result)
            print(json.dumps(result), flush=True)
            if not result['accepted'] or Q(result['radius']) <= 0:
                return tiles, segments, False
            axis = 0 if target[0] != current[0] else 1
            remaining = target[axis]-current[axis]
            step = min(abs(remaining), Q(result['radius']))*(1 if remaining > 0 else -1)
            next_point = list(current)
            next_point[axis] += step
            segments.append(dict(start=list(map(str, current)), end=list(map(str, next_point)),
                                 tile=len(tiles)-1))
            current = tuple(next_point)
    return tiles, segments, True


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('matrix', type=Path)
    parser.add_argument('--out', type=Path)
    parser.add_argument('--cover-contour', action='store_true')
    parser.add_argument('--include-background', action='store_true',
                        help='inclui a bola RT e o erro de mu; usa os arquivos RefA e RefAplusB')
    args = parser.parse_args()
    if args.include_background and hashlib.sha256(args.matrix.read_bytes()).hexdigest() != \
            'a3335d9550f9567db287c1ae105674a5914446f9d3772d3955e3628a747c5846':
        parser.error('a ponte para a bola RT foi auditada somente para a matriz RefA 6x18 identificada por hash')
    data = prepare(args.matrix, 4, 12)
    if args.include_background:
        from reference_error_budget import read_mu, error_budget
        root = Path(__file__).resolve().parents[1]
        budget = error_budget(read_mu(root/'.cache/rt-1203.3766v1/sourcecode/RefA.dat'),
                              read_mu(root/'.cache/rt-1203.3766v1/RefAplusB.dat'))
        data['background_H_error'] = Q(budget['total_H_error'])
    tiles, segments, covered = [], [], False
    if args.cover_contour:
        tiles, segments, covered = cover_upper_boundary(data)
    else:
        for real, imag in ((Q(1), Q(0)), (Q(1), Q(1, 4)),
                           (Q(1, 8), Q(1, 4)), (Q(1, 8), Q(0))):
            result = tile(data, real, imag)
            tiles.append(result)
            print(json.dumps(result), flush=True)
    report = dict(matrix=str(args.matrix),
                  sha256=hashlib.sha256(args.matrix.read_bytes()).hexdigest(),
                  finite_crown_dimension=data['size'], core=[4, 12],
                  weights=[916, 4096, 243, 1233],
                  inverse_free_error=str(data['qerror']), tiles=tiles,
                  inverse_free_norm_bound=str(data['inverse_norm_bound']),
                  background_H_error=str(data.get('background_H_error', Q(0))),
                  conditional_on_RT_background_ball=args.include_background,
                  contour_coverage=covered, upper_boundary_segments=segments,
                  lower_boundary_by_real_conjugation=covered,
                  infinite_tail_bounds=None,
                  physical_certificate=False)
    if args.out:
        args.out.write_text(json.dumps(report, indent=2)+'\n')
