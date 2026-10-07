"""Bounds sharp nos pesos menores; nao valida os acoplamentos finito/cauda."""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
from background_couplings import RT_BACKGROUND_ERROR
from component_exterior import read_fields, field_expression, field_norm, weighted_bound
from sharp_propagation import table
from reference_error_budget import read_mu
from check_free_preconditioner import free_matrix


def free_tail(m, n, radial=Q(9, 8), mu_min=Q(1, 6), mu_max=Q(17, 100)):
    if m < 1 or n < 2 or radial <= 1 or not 0 < mu_min <= mu_max:
        raise ValueError('parametros de cauda invalidos')
    r, b = 1/radial, Q(m, 2)
    fourier = Q(3, 2)*max(1/(b*(1-r)),
                 1/(b*(1-r*r))+2*mu_max*r*r/(b*b*(1-r*r)))
    radial_factor = max(1+2*r/(1-r), 1+Q(2*n, n-1)*r*r/(1-r*r))
    return max(fourier, Q(3, 2)*radial_factor/(mu_min*(n+1)))


def sharp_majorant(fields, mu_max=Q(17, 100)):
    rows = table('GAMMA2_SHARP_macro.c', '_G2SHARPPAIR_')
    norm = lambda expr: field_norm(field_expression(expr, fields),
                                    fourier=Q(129, 128), radial=Q(9, 8))
    result = [[norm(rows[i, 1])/2,
               max(norm(rows[i, 4]), norm(rows[i, 5]))/2] for i in range(2)]
    result[1][0] += mu_max*Q(144, 17)
    return result


def mu_error(root):
    return abs(read_mu(root/'.cache/rt-1203.3766v1/RefAplusB.dat')
               -read_mu(root/'.cache/rt-1203.3766v1/sourcecode/RefA.dat'))+Q(1, 2**277)


def finite_background_error(addresses, weights, epsilon_mu, epsilon=RT_BACKGROUND_ERROR):
    if epsilon_mu < 0 or epsilon < 0:
        raise ValueError('erros nao negativos exigidos')
    principal = free_matrix(addresses, Q(1))
    # Remove partial_tau: seus unicos elementos trocam a parte real/imaginaria.
    jnorm = max(sum(abs(row[j])*weights[i] for i, row in enumerate(principal)
                    if addresses[i][1] == addresses[j][1])/weights[j]
                for j in range(len(weights)))
    return (jnorm+Q(144, 17))*epsilon_mu+3*epsilon


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('background', type=Path)
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    if hashlib.sha256(args.background.read_bytes()).hexdigest() != \
            '5776e6038f994e06b4b82fc638d797154f630aeee37b251058512d3cb4cd450b':
        parser.error('a bola usada corresponde ao RefA auditado')
    majorant = sharp_majorant(read_fields(args.background))
    beta = weighted_bound(majorant, [Q(1), Q(1)])+3*RT_BACKGROUND_ERROR
    beta = Q(-((-beta.numerator*1000000)//beta.denominator), 1000000)
    principal_error = 210*mu_error(Path(__file__).resolve().parents[1])
    rows = []
    # Todo o disco |s-.733|<=1/32 esta em |s|<4/5.
    for m, n in ((8, 24), (64, 192), (128, 768), (256, 1536)):
        d = (beta+Q(4, 5))*free_tail(m, n)+principal_error
        rows.append(dict(outer=[m, n], TT_defect=str(d),
                         TT_defect_decimal=float(d), diagonal_contraction=d < 1))
    result = dict(weights=['129/128', '9/8'], beta=str(beta), beta_decimal=float(beta),
                  principal_relative_error=str(principal_error), rows=rows,
                  conditional_on_published_RT_ball=True,
                  finite_tail_couplings_verified=False, physical_certificate=False)
    output = json.dumps(result, indent=2)+'\n'
    if args.out:
        args.out.write_text(output)
    print(output, end='')
