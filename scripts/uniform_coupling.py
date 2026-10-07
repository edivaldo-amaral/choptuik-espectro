"""Criterio uniforme de Schur/Neumann com dados de tile racionais.

Todos os dados devem se referir a MESMA caixa finita, norma e fundo.
O checker aritmetico nao verifica essa proveniencia nem a cobertura do contorno.
"""
from fractions import Fraction as Q

from check_frequency_blocks import two_block_test


def uniform_tile(a0, variation, radius, background_error, b, c, d):
    values = tuple(map(Q, (a0, variation, radius, background_error, b, c, d)))
    if min(values) < 0:
        raise ValueError('majorantes e raio devem ser nao negativos')
    a0, variation, radius, background_error, b, c, d = values
    a = a0+radius*variation+background_error
    result = dict(a=str(a), b=str(b), c=str(c), d=str(d),
                  **two_block_test(a, b, c, d))
    result['schur_defect'] = str(a+b*c/(1-d)) if d < 1 else None
    result['physical_certificate'] = False
    return result


def admissible_radius(a0, variation, background_error, b, c, d):
    """Supremo ESTRITO do raio; None se nem o centro passa ou limite infinito."""
    vals = tuple(map(Q, (a0, variation, background_error, b, c, d)))
    if min(vals) < 0:
        raise ValueError('majorantes nao negativos exigidos')
    a0, variation, error, b, c, d = vals
    if d >= 1:
        return dict(center_accepted=False, strict_radius_upper=None)
    margin = 1-a0-error-b*c/(1-d)
    return dict(center_accepted=margin > 0,
                strict_radius_upper=str(margin/variation) if margin > 0 and variation else None)
