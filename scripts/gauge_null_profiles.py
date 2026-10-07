"""Geradores locais f(u)=u^n do gauge residual de coordenadas nulas.

Inclui a acao RT da translacao nula conjunta; nao prova equivalencia global.
"""
from fractions import Fraction as Q
from sharp_propagation import Poly, X


def power(p, n):
    result = Poly(1)
    for _ in range(n):
        result = result*p
    return result


def null_generator(n, mu=Q(1, 6)):
    if not isinstance(n, int) or n < 0 or mu <= 0:
        raise ValueError('n inteiro nao negativo e mu positivo exigidos')
    plus, minus = (-1)**n*power(1+X, n), (-1)**n*power(1-X, n)
    tau = (plus+minus)/(2*mu)
    xi = ((1+X)*minus-(1-X)*plus)/2
    return mu*(1-n), tau, xi


def phase_transversality(imaginary_reference, background_error):
    """b(h)=Re h4_(1,0); b(partial_tau w)=-Im w4_(1,0)/2."""
    if background_error < 0:
        raise ValueError('erro nao negativo exigido')
    return (abs(Q(imaginary_reference))-Q(background_error)/Q(65, 32))/2


def translation_profile(p, mu):
    """Perfil de Lie em s=mu: (1+partial_tau/mu+xi partial_xi)p."""
    if mu <= 0:
        raise ValueError('mu positivo exigido')
    return p+p.derivative(0)/mu+X*p.derivative(1)
