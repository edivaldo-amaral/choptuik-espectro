"""Norma livre sharp global: colunas finitas auditadas mais cauda analitica."""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
from analyze_sharp import constraint_addresses
from certify_sharp_disk import certify
from check_free_preconditioner import free_matrix, integer_matrix
from sharp_exterior import free_tail


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    path = root/'build/spectrum/rt-sharp-12x36.dat'
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != 'aaaa6f933e4e12d3d5eb557d9e9130f22a50eddf1d6995bdfa55ca7ea4fbf338':
        parser.error('enderecos devem vir da matriz sharp auditada')
    addresses = constraint_addresses(path)
    weights = [(2-(m == 0))*(2-(n == 0))*Q(129,128)**m*Q(9,8)**n for _, _, m, n in addresses]
    mu = Q(722873400, 2**32)
    raw, denominator = integer_matrix(free_matrix(addresses, mu))
    if denominator & (denominator-1):
        raise ValueError('principal nao diadico')
    result = certify(raw, -(denominator.bit_length()-1), weights, center=Q(0), radius=Q(0))
    if not result['accepted']:
        raise ValueError('a inversa livre finita nao foi validada')
    tail = free_tail(12, 36)
    result.update(operator='Jsharp_muRefA, NOT Ksharp', mu=str(mu), box=[12,36],
                  weights=['129/128','9/8'], address_source_sha256=digest,
                  input_tail_bound=str(tail),
                  global_inverse_bound=str(max(Q(result['inverse_bound']),tail)),
                  global_inverse_bound_decimal=float(max(Q(result['inverse_bound']),tail)),
                  global_claim_requires_analytic_free_tail=True)
    if args.out:
        args.out.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('center_defect','input_tail_bound','global_inverse_bound_decimal')}, indent=2))
