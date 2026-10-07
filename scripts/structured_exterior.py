"""Cauda de entrada melhorada pela cancelacao nas colunas da inversa livre.

Base REAL >=0; normas nos DOFs reais complexificados. Prova em
docs/STRUCTURED_EXTERIOR.md. Aritmetica racional nao certifica a PDE.
"""
from fractions import Fraction as Q
import json

from check_free_preconditioner import input_tail


def structured_input_tail(m: int, n: int, base: Q = Q(0),
                          mu_min: Q = Q(1, 6), mu_max: Q = Q(17, 100)) -> Q:
    if m < 1 or n < 2 or base < 0 or not 0 < mu_min <= mu_max:
        raise ValueError('exige M>=1, N>=2, base>=0 e 0<mu_min<=mu_max')
    # Cada coeficiente complexo corresponde a um bloco real 2x2:
    # ||bloco||_1 <= sqrt(2)|coeficiente| <= (3/2)|coeficiente|.
    b = Q(m, 2)
    r = Q(4, 5)
    fourier = max(Q(3, 2)/(b*(1-r)),
                  Q(3, 2)*(1/(b*(1-r*r))
                            + 2*mu_max*r*r/(b*b*(1-r*r))))
    radial = Q(27, 2)/(mu_min*(n+1))
    return min(input_tail(m, n, base, mu_min), max(fourier, radial))


def report(core=(12, 36), outer=(64, 192), distance=Q(5, 4)):
    if distance < 0 or any(p >= g for p, g in zip(core, outer)):
        raise ValueError('coroa estrita e distancia nao negativa exigidas')
    old = input_tail(*outer, Q(0))
    new = structured_input_tail(*outer)
    return dict(core=core, outer=outer, max_distance=str(distance),
                old_tail=str(old), structured_tail=str(new),
                T_from_S=str(23*structured_input_tail(*core)),
                S_from_T=str((23+distance)*new),
                TT_defect=str((23+distance)*new),
                TT_decimal=float((23+distance)*new),
                exterior_contraction=(23+distance)*new < 1,
                physical_certificate=False)


if __name__ == '__main__':
    print(json.dumps([report(outer=g) for g in
                      ((24, 72), (64, 192), (364, 1964))], indent=2))
