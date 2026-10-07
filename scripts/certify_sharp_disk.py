"""Exclusao racional de um disco para a matriz sharp FINITA de referencia."""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import numpy as np
from analyze_sharp import constraint_addresses
from build_contour import read_spectral_matrix
from certify_crown_tiles import exact_product, nearest_integer, norm_real


def check_compression(outer, outer_exponent, outer_addresses,
                      inner, inner_exponent, inner_addresses):
    lookup = {address: i for i, address in enumerate(outer_addresses)}
    if len(lookup) != len(outer_addresses) or len(set(inner_addresses)) != len(inner_addresses):
        raise ValueError('enderecos duplicados')
    if not all(address in lookup for address in inner_addresses):
        raise ValueError('caixas nao aninhadas')
    indices = [lookup[address] for address in inner_addresses]
    for i, oi in enumerate(indices):
        for j, oj in enumerate(indices):
            if Q(outer[oi][oj])*Q(2)**outer_exponent != Q(inner[i][j])*Q(2)**inner_exponent:
                raise ValueError('a compressao difere da matriz interna')
    return dict(exact_compression_match=True, inner_dimension=len(inner),
                shell_dimension=len(outer)-len(inner),
                all_shell_couplings_in_full_inverse=True)


def certify(raw, exponent, weights, center=Q(733, 1000), radius=Q(1, 32), bits=48,
            background_error=Q(0), principal=None):
    n, grid = len(raw), 1 << bits
    if (not n or radius < 0 or bits < 0 or background_error < 0 or len(weights) != n
            or min(weights) <= 0 or any(len(row) != n for row in raw)):
        raise ValueError('matriz quadrada, raio nao negativo e pesos positivos exigidos')
    rounded = np.array([[nearest_integer((Q(raw[i][j])*Q(2)**exponent
                         + (center if i == j else 0))*weights[i]/weights[j]*grid)
                         for j in range(n)] for i in range(n)], dtype=object)
    guess = np.linalg.inv(np.array(rounded, dtype=float)/grid)
    v = np.array([[round(float(x)*grid) for x in row] for row in guess], dtype=object)
    residual = np.eye(n, dtype=object)*grid**2-exact_product(v, rounded)
    vnorm = norm_real(v, grid)
    center_defect = norm_real(residual, grid**2)+vnorm*(Q(n, 2*grid)+background_error)
    uniform = center_defect+radius*vnorm
    result = dict(center=str(center), radius=str(radius), center_defect=str(center_defect),
                inverse_candidate_norm=str(vnorm), uniform_defect=str(uniform),
                uniform_defect_decimal=float(uniform), accepted=uniform < 1,
                inverse_bound=str(vnorm/(1-uniform)) if uniform < 1 else None,
                background_error_included=background_error > 0,
                background_operator_error=str(background_error), infinite_tail_included=False,
                physical_certificate=False)
    if principal is not None:
        if len(principal) != n or any(len(row) != n for row in principal):
            raise ValueError('principal com dimensao incompativel')
        right_residual = np.eye(n, dtype=object)*grid**2-exact_product(rounded, v)
        right_defect = norm_real(right_residual, grid**2)+vnorm*(Q(n,2*grid)+background_error+radius)
        jr = np.array([[nearest_integer(principal[i][j]*weights[i]/weights[j]*grid)
                        for j in range(n)] for i in range(n)], dtype=object)
        jvnorm = norm_real(exact_product(jr, v), grid**2)+Q(n,2*grid)*vnorm
        result['right_uniform_defect'] = str(right_defect)
        result['principal_times_candidate_norm'] = str(jvnorm)
        result['preconditioned_inverse_bound'] = str(jvnorm/(1-right_defect)) if right_defect < 1 else None
        result['preconditioned_inverse_bound_decimal'] = float(jvnorm/(1-right_defect)) if right_defect < 1 else None
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('matrix', type=Path)
    parser.add_argument('--out', type=Path)
    parser.add_argument('--include-background', action='store_true')
    parser.add_argument('--inner-matrix', type=Path)
    parser.add_argument('--audit-preconditioned', action='store_true')
    args = parser.parse_args()
    n, raw, exponent = read_spectral_matrix(args.matrix, 'CHOPTUIK_RT_SHARP_MATRIX_V1')
    addresses = constraint_addresses(args.matrix)
    if len(addresses) != n:
        parser.error('dimensao inconsistente')
    weights = [(2-(m == 0))*(2-(cn == 0))*Q(129, 128)**m*Q(9, 8)**cn
               for _, _, m, cn in addresses]
    background_error = Q(0)
    if args.include_background:
        if hashlib.sha256(args.matrix.read_bytes()).hexdigest() not in {
                'fd27e0221aafc94b5701ca1fb2f210dbfed5d20f43334436fcc7323ce1931e62',
                'aaaa6f933e4e12d3d5eb557d9e9130f22a50eddf1d6995bdfa55ca7ea4fbf338'}:
            parser.error('ponte de fundo validada somente para sharp RefA 8x24 e 12x36')
        from sharp_exterior import finite_background_error, mu_error
        background_error = finite_background_error(addresses, weights,
                            mu_error(Path(__file__).resolve().parents[1]))
    principal = None
    if args.audit_preconditioned:
        from check_free_preconditioner import free_matrix
        from reference_error_budget import read_mu
        principal = free_matrix(addresses, read_mu(Path(__file__).resolve().parents[1]/
                                  '.cache/rt-1203.3766v1/sourcecode/RefA.dat'))
    result = certify(raw, exponent, weights, background_error=background_error, principal=principal)
    if args.inner_matrix:
        _, inner, inner_exponent = read_spectral_matrix(args.inner_matrix, 'CHOPTUIK_RT_SHARP_MATRIX_V1')
        result['nested_shell'] = check_compression(raw, exponent, addresses,
                    inner, inner_exponent, constraint_addresses(args.inner_matrix))
        result['nested_shell']['inner_sha256'] = hashlib.sha256(args.inner_matrix.read_bytes()).hexdigest()
    result['conditional_on_published_RT_ball'] = args.include_background
    result.update(matrix=str(args.matrix), dimension=n, weights=['129/128', '9/8'],
                  sha256=hashlib.sha256(args.matrix.read_bytes()).hexdigest())
    output = json.dumps(result, indent=2)+'\n'
    if args.out:
        args.out.write_text(output)
    print(output, end='')
