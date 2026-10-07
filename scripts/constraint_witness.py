"""Testemunho racional condicional de constraint NAO nula num modo encerrado."""
from fractions import Fraction as Q


def violation_lower_bound(reference_lower, applied_error, operator_bound,
                          vector_error, spectral_error, vector_norm=Q(1), shift_bound=Q(9, 8)):
    values = tuple(map(Q, (reference_lower, applied_error, operator_bound,
                           vector_error, spectral_error, vector_norm, shift_bound)))
    if min(values) < 0:
        raise ValueError('bounds nao negativos exigidos')
    ref, applied, operator, vector, spectral, norm, shift = values
    lower = ref-applied-operator*vector-spectral*shift*(norm+vector)
    return dict(lower=str(lower), excludes_constraint_kernel=lower > 0,
                requires_certified_eigenvector_enclosure=True, physical_certificate=False)
