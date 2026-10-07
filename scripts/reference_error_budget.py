"""Erro do fundo RT no pencil precondicionado com mu de RefA FIXO."""
from fractions import Fraction as Q
import json
from pathlib import Path

from background_couplings import RT_BACKGROUND_ERROR
from component_exterior import gamma1_majorant, weighted_bound


def read_mu(path):
    with Path(path).open() as stream:
        lines = [next(stream).strip() for _ in range(4)]
    if lines[:2] != ['mu_MultiField', 'DyadicQ'] or not lines[2].startswith('TwoExp ') or not lines[3].startswith('coeff '):
        raise ValueError('cabecalho mu RT invalido')
    return Q(int(lines[3].split()[1]))*Q(2)**int(lines[2].split()[1])


def error_budget(mu_ref, mu_refined, omega_error=RT_BACKGROUND_ERROR):
    if not Q(1, 6) <= mu_ref <= Q(17, 100) or omega_error < 0:
        raise ValueError('parametros fora dos bounds usados')
    mu_error = abs(mu_refined-mu_ref)+Q(1, 2**277)
    eta = list(map(Q, (916, 4096, 243, 1233)))
    gamma1 = weighted_bound(gamma1_majorant(Q(1)), eta)
    principal_error = 114*mu_error
    bounded_error = 108*(gamma1*mu_error+6*max(eta)/min(eta)*omega_error)
    total = principal_error+bounded_error
    return dict(mu_reference=str(mu_ref), mu_refined=str(mu_refined),
                mu_error=str(mu_error), mu_error_decimal=float(mu_error),
                gamma1_weighted=str(gamma1), principal_relative_error=str(principal_error),
                total_H_error=str(total), total_H_error_decimal=float(total),
                conditional_on_published_RT_ball=True, physical_certificate=False)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    result = error_budget(read_mu(root/'.cache/rt-1203.3766v1/sourcecode/RefA.dat'),
                          read_mu(root/'.cache/rt-1203.3766v1/RefAplusB.dat'))
    output = json.dumps(result, indent=2)+'\n'
    print(output, end='')
    if args.out:
        args.out.write_text(output)
