"""Aritmetica da recuperacao de raio para raizes selecionadas; prova em papel."""
from fractions import Fraction as Q
import json
from sharp_exterior import free_tail
from structured_exterior import structured_input_tail


def transfer_bounds(m=1024, n=4096, spectral_radius=Q(5, 4)):
    if spectral_radius < 0:
        raise ValueError('raio espectral nao negativo exigido')
    strong_beta = Q(5191, 500)
    weak_beta = Q(162, 85)*strong_beta
    strong = (strong_beta+spectral_radius)*structured_input_tail(m, n)
    weak = (weak_beta+spectral_radius)*free_tail(m, n)
    return dict(box=[m, n], spectral_radius=str(spectral_radius),
                strong_beta=str(strong_beta), weak_beta=str(weak_beta),
                strong_defect=str(strong), weak_defect=str(weak),
                strong_defect_decimal=float(strong), weak_defect_decimal=float(weak),
                both_contract=strong < 1 and weak < 1,
                preconditioner='J_mu_true^-1, NOT J_mu_RefA^-1',
                conditional_on_published_RT_ball=True, physical_certificate=False)


if __name__ == '__main__':
    import argparse
    from pathlib import Path
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    output = json.dumps(transfer_bounds(), indent=2)+'\n'
    if args.out:
        args.out.write_text(output)
    print(output, end='')
