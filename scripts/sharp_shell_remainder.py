"""Propaga o erro da inversa exterior ate o Schur da faixa intermediaria.

Nao constroi a parametriz da faixa. Os bounds analiticos dos blocos e a
bola RT continuam hipoteses; toda a aritmetica deste produtor e racional.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path

from sharp_exterior import free_tail


def schur_factors(kappa, ft, tf, uf, fu, ut, tu):
    values = list(map(Q, (kappa, ft, tf, uf, fu, ut, tu)))
    if min(values) < 0:
        raise ValueError('majorantes nao negativos exigidos')
    kappa, ft, tf, uf, fu, ut, tu = values
    return ut+uf*kappa*ft, tu+tf*kappa*fu


def remainder(delta, left, right, degree):
    delta, left, right = map(Q, (delta, left, right))
    if not 0 <= delta < 1 or min(left, right) < 0 or not isinstance(degree, int) or degree < 0:
        raise ValueError('contracao, normas e grau invalidos')
    return left*right*delta**(degree+1)/(1-delta)


def required_degree(delta, left, right, tolerance):
    tolerance = Q(tolerance)
    if tolerance <= 0:
        raise ValueError('tolerancia positiva exigida')
    if remainder(delta, left, right, 0) <= tolerance:
        return 0
    low, high = 0, 1
    while remainder(delta, left, right, high) > tolerance:
        low, high = high, 2*high
    while high-low > 1:
        mid = (low+high)//2
        if remainder(delta, left, right, mid) <= tolerance:
            high = mid
        else:
            low = mid
    return high


def close_shell(candidate_norm, finite_defect, operator_error):
    candidate_norm, finite_defect, operator_error = map(
        Q, (candidate_norm, finite_defect, operator_error))
    if min(candidate_norm, finite_defect, operator_error) < 0:
        raise ValueError('majorantes nao negativos exigidos')
    total = finite_defect+candidate_norm*operator_error
    return dict(total_defect=str(total), accepted=total < 1)


def ceil_bound(value, denominator=10**6):
    value = Q(value)*denominator
    return Q(-(-value.numerator//value.denominator), denominator)


def report_from_bridge(bridge, tolerance=Q(1,10**12)):
    if (bridge['inner'] != [12,36] or not bridge['finite_plus_far_tail_invertible']
            or not bridge['conditional_on_published_RT_ball']):
        raise ValueError('ponte sharp 12x36 com cauda contrativa exigida')
    kappa, c, d, delta, qbound = (Q(bridge[key]) for key in
        ('preconditioned_finite_inverse_bound', 'outward', 'tail_defect',
         'schur_defect', 'free_bound_used'))
    if (min(kappa, c, d, qbound) < 0 or not 0 <= delta < 1
            or delta != d+c*kappa*d or Q(bridge['inward']) != d):
        raise ValueError('identidade escalar de Schur inconsistente')
    outer = bridge['outer']
    if len(outer) != 2 or outer[0] <= 12 or outer[1] <= 36:
        raise ValueError('caixas nao aninhadas')
    beta_s = Q(462237,100000)+Q(4,5)
    principal_error = d-beta_s*free_tail(*outer)
    if principal_error < 0:
        raise ValueError('erro relativo do principal negativo')
    # U esta fora de F. Todos os acoplamentos U<->F,T sao mantidos.
    from_u = beta_s*free_tail(12,36)+principal_error
    uf = Q(462237,100000)*qbound
    left, right = schur_factors(kappa, d, c, uf, from_u, d, from_u)
    # Arredondamento para fora evita elevar denominadores gigantes a N.
    delta_up, left_up, right_up = map(ceil_bound, (delta, left, right))
    if delta_up >= 1:
        raise ValueError('grade insuficiente para preservar contracao estrita')
    degree = required_degree(delta_up, left_up, right_up, tolerance)
    error = remainder(delta_up, left_up, right_up, degree)
    return dict(bridge_outer=outer, remaining_shell_dimension=bridge['remaining_shell_dimension'],
                delta_upper=str(delta_up), left_factor_upper=str(left_up),
                right_factor_upper=str(right_up), target_error=str(tolerance),
                neumann_degree=degree, terms=degree+1,
                propagated_error_upper=str(ceil_bound(error, 10**18)),
                target_verified_exactly=error <= tolerance,
                finite_shell_parametrix_constructed=False,
                full_sharp_invertibility=False, physical_certificate=False,
                conditional_on_input_block_bounds_and_RT_ball=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('bridge', type=Path)
    parser.add_argument('--tolerance', type=Q, default=Q(1,10**12))
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    raw = args.bridge.read_bytes()
    result = report_from_bridge(json.loads(raw), args.tolerance)
    result['bridge_sha256'] = hashlib.sha256(raw).hexdigest()
    output = json.dumps(result, indent=2)+'\n'
    if args.out:
        args.out.write_text(output)
    print(output, end='')
