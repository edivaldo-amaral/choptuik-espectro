"""Acoplamento G->exterior distante por caudas de componentes, com erro RT.

Nao elimina a coroa intermediaria nem certifica o bloco exterior inverso.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path

from background_couplings import RT_BACKGROUND_ERROR
from component_exterior import read_fields, component_majorant, weighted_bound
from check_crown_cover import check
from reference_error_budget import read_mu, error_budget
from structured_exterior import structured_input_tail

ETA = list(map(Q, (916, 4096, 243, 1233)))


def guarded_fields(fields, inner, outer):
    if min(inner) < 1 or any(a > b for a, b in zip(inner, outer)):
        raise ValueError('caixas positivas aninhadas exigidas')
    dm, dn = (b-a+1 for a, b in zip(inner, outer))
    return [{key: value for key, value in field.items()
             if key[0] >= dm or key[1] >= dn} for field in fields]


def outward_bound(fields, inner, outer, free_bound, error=RT_BACKGROUND_ERROR):
    if free_bound <= 0 or error < 0:
        raise ValueError('bound livre positivo e erro nao negativo exigidos')
    majorant = component_majorant(guarded_fields(fields, inner, outer))
    beta = weighted_bound(majorant, ETA)+6*max(ETA)/min(ETA)*error
    return beta*free_bound


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('background', type=Path)
    parser.add_argument('crown_report', type=Path)
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    report = json.loads(args.crown_report.read_text())
    if (report['sha256'] != 'a3335d9550f9567db287c1ae105674a5914446f9d3772d3955e3628a747c5846'
        or report['weights'] != [916, 4096, 243, 1233]
        or hashlib.sha256(args.background.read_bytes()).hexdigest() !=
            '5776e6038f994e06b4b82fc638d797154f630aeee37b251058512d3cb4cd450b'):
        parser.error('dados diferentes da ponte RefA 6x18 auditada')
    qbound = Q(report['inverse_free_norm_bound'])
    # Norma l1: particionar as COLUNAS permite tomar maximo, nao soma.
    global_qbound = max(qbound, structured_input_tail(6, 18))
    if report.get('conditional_on_RT_background_ball') is not True:
        parser.error('a inversa da coroa deve incluir o erro do fundo')
    kappa = Q(check(report)['uniform_inverse_bound'])
    root = Path(__file__).resolve().parents[1]
    mu_budget = error_budget(read_mu(args.background),
                            read_mu(root/'.cache/rt-1203.3766v1/RefAplusB.dat'))
    principal_error = Q(mu_budget['principal_relative_error'])
    fields = read_fields(args.background)
    rows = []
    for outer in ((6, 18), (12, 36), (24, 72), (64, 192), (176, 944)):
        bound = outward_bound(fields, (6, 18), outer, qbound)
        # S=(6x18)-(4x12); T=I-outer. O restante intermediario NAO foi removido.
        d = (Q(5191, 500)+Q(5, 4))*structured_input_tail(*outer)+principal_error
        schur = d+kappa*bound*d
        rows.append(dict(inner=[6, 18], outer=outer, T_from_G=str(bound),
                         T_from_G_decimal=float(bound),
                         S_from_T=str(d), TT_defect=str(d),
                         S_plus_T_schur_defect=str(schur),
                         S_plus_T_schur_defect_decimal=float(schur),
                         S_plus_T_invertible=schur < 1,
                         intermediate_shell_required=outer != (6, 18)))
    improved_background_error = principal_error+global_qbound*(
        Q(mu_budget['gamma1_weighted'])*Q(mu_budget['mu_error'])
        +6*max(ETA)/min(ETA)*RT_BACKGROUND_ERROR)
    result = dict(free_bound=str(qbound), global_free_bound=str(global_qbound),
                  improved_background_H_error=str(improved_background_error),
                  improved_background_H_error_decimal=float(improved_background_error),
                  crown_inverse_bound=str(kappa), rows=rows,
                  background_sha256=hashlib.sha256(args.background.read_bytes()).hexdigest(),
                  crown_report_sha256=hashlib.sha256(args.crown_report.read_bytes()).hexdigest(),
                  uniform_in_s=True, conditional_on_published_RT_ball=True,
                  physical_certificate=False)
    output = json.dumps(result, indent=2)+'\n'
    if args.out:
        args.out.write_text(output)
    print(json.dumps([{k: row[k] for k in ('outer', 'T_from_G_decimal',
                     'S_plus_T_schur_defect_decimal', 'S_plus_T_invertible')} for row in rows], indent=2))
