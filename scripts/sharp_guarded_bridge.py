"""Schur sharp sobre F=(12x36) mais cauda infinita, mantendo a faixa omitida."""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
from analyze_sharp import constraint_addresses
from background_couplings import RT_BACKGROUND_ERROR
from component_exterior import read_fields, weighted_bound
from guarded_exterior import guarded_fields
from check_free_preconditioner import free_matrix
from reference_error_budget import read_mu
from sharp_exterior import sharp_majorant, free_tail, mu_error


def schur_bound(finite_inverse, outward, inward, tail_defect):
    if min(finite_inverse, outward, inward, tail_defect) < 0:
        raise ValueError('majorantes nao negativos exigidos')
    return tail_defect+outward*finite_inverse*inward


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path)
    parser.add_argument('--outer-m', type=int, default=256)
    parser.add_argument('--outer-n', type=int, default=1536)
    parser.add_argument('--finite-report', type=Path)
    parser.add_argument('--free-report', type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    matrix = root/'build/spectrum/rt-sharp-12x36.dat'
    report_path = args.finite_report or root/'build/spectrum/sharp-disk-12x36-background.json'
    background = root/'.cache/rt-1203.3766v1/sourcecode/RefA.dat'
    report = json.loads(report_path.read_text())
    expected = 'aaaa6f933e4e12d3d5eb557d9e9130f22a50eddf1d6995bdfa55ca7ea4fbf338'
    if (hashlib.sha256(matrix.read_bytes()).hexdigest() != expected or report['sha256'] != expected
        or not report['accepted'] or not report['conditional_on_published_RT_ball']
        or report['weights'] != ['129/128', '9/8']
        or Q(report['center']) != Q(733, 1000) or Q(report['radius']) != Q(1, 32)
        or hashlib.sha256(background.read_bytes()).hexdigest() !=
            '5776e6038f994e06b4b82fc638d797154f630aeee37b251058512d3cb4cd450b'):
        parser.error('dados diferentes da ponte auditada')
    defect, vnorm = Q(report['uniform_defect']), Q(report['inverse_candidate_norm'])
    if not 0 <= defect < 1 or Q(report['inverse_bound']) != vnorm/(1-defect):
        parser.error('bound de inversa inconsistente')
    addresses = constraint_addresses(matrix)
    weights = [(2-(m == 0))*(2-(n == 0))*Q(129,128)**m*Q(9,8)**n for _, _, m, n in addresses]
    principal = free_matrix(addresses, read_mu(background))
    jnorm = max(sum(abs(row[j])*weights[i] for i, row in enumerate(principal))/weights[j]
                for j in range(len(weights)))
    kappa = jnorm*Q(report['inverse_bound'])
    if report.get('preconditioned_inverse_bound') is not None:
        right = Q(report['right_uniform_defect'])
        direct = Q(report['preconditioned_inverse_bound'])
        if not 0 <= right < 1 or direct != Q(report['principal_times_candidate_norm'])/(1-right):
            parser.error('auditoria precondicionada inconsistente')
        kappa = min(kappa, direct)
    qbound = Q(204)
    if args.free_report:
        free = json.loads(args.free_report.read_text())
        if (free['operator'] != 'Jsharp_muRefA, NOT Ksharp' or free['address_source_sha256'] != expected
            or free['box'] != [12,36] or free['weights'] != ['129/128','9/8']
            or Q(free['center']) != 0 or Q(free['radius']) != 0
            or Q(free['mu']) != read_mu(background) or not free['accepted']):
            parser.error('relatorio livre incompativel')
        free_defect = Q(free['uniform_defect'])
        if not 0 <= free_defect < 1 or Q(free['inverse_bound']) != Q(free['inverse_candidate_norm'])/(1-free_defect):
            parser.error('bound livre inconsistente')
        qbound = min(qbound,max(Q(free['inverse_bound']),free_tail(12,36)))
    outer = (args.outer_m, args.outer_n)
    if outer[0] <= 12 or outer[1] <= 36:
        parser.error('a caixa externa deve conter estritamente 12x36')
    filtered = guarded_fields(read_fields(background), (12, 36), outer)
    beta_guard = weighted_bound(sharp_majorant(filtered, mu_max=Q(0)), [Q(1), Q(1)])
    outward = qbound*(beta_guard+3*RT_BACKGROUND_ERROR)
    d = (Q(462237,100000)+Q(4,5))*free_tail(*outer)+210*mu_error(root)
    delta = schur_bound(kappa, outward, d, d)
    result = dict(inner=[12,36], outer=outer, principal_finite_norm=str(jnorm),
                  free_bound_used=str(qbound),
                  remaining_shell_dimension=(2*((outer[0]+1)//2)-1)*(outer[1]+outer[1]//2)-len(addresses),
                  preconditioned_finite_inverse_bound=str(kappa),
                  preconditioned_finite_inverse_decimal=float(kappa),
                  outward=str(outward), outward_decimal=float(outward),
                  inward=str(d), tail_defect=str(d),
                  schur_defect=str(delta), schur_defect_decimal=float(delta),
                  finite_plus_far_tail_invertible=delta < 1,
                  full_sharp_invertibility=False, intermediate_shell_required=True,
                  conditional_on_published_RT_ball=True, physical_certificate=False,
                  matrix_sha256=expected,
                  free_report_sha256=hashlib.sha256(args.free_report.read_bytes()).hexdigest() if args.free_report else None,
                  finite_report_sha256=hashlib.sha256(report_path.read_bytes()).hexdigest())
    if args.out:
        args.out.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k,v in result.items() if 'decimal' in k or type(v) is bool}, indent=2))
