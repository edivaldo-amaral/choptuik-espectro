"""Testes algebricos do transporte de gauge; periodicidade global em papel."""
from fractions import Fraction as Q
from math import comb
from sharp_propagation import Poly, X
from geometric_reconstruction import inverse_tau_polynomial
from gauge_null_profiles import power


def transport(y, mu, shift, center):
    return y.derivative(0)+mu*(X-center)*y.derivative(1)+(shift-mu)*y


def inverse_polynomial(source, mu, shift, center):
    mu, shift, center = Q(mu), Q(shift), Q(center)
    if mu <= 0:
        raise ValueError('mu positivo exigido')
    coefficients = {}
    for (a, b, k), value in Poly(source).terms.items():
        for n in range(b+1):
            term = Poly({(a, 0, k): value*comb(b,n)*center**(b-n)})
            coefficients[n] = coefficients.get(n, Poly())+term
    result = Poly()
    for n, coefficient in coefficients.items():
        if coefficient == 0:
            continue
        denominator = shift+mu*(n-1)
        if denominator == 0:
            raise ValueError('ressonancia radial: compatibilidade temporal nao fornecida')
        result += inverse_tau_polynomial(coefficient, denominator)*power(X-center, n)
    return result


def coordinate_profiles(yminus, yplus, mu):
    return ((yplus+yminus)/(2*mu), ((1+X)*yminus-(1-X)*yplus)/2)


def inverse_periodic_mode(real, imag, mu, shift_real, frequency, center,
                          shift_imag=0):
    """Inverte um modo exp(i*m*tau/2), com coeficientes radiais racionais.

    Retorna as partes real/imaginaria do perfil. Na ressonancia compativel,
    escolhe zero para o coeficiente livre; nao remove uma obstrucao nao nula.
    Ao contrario de inverse_polynomial, este teste respeita periodicidade.
    """
    mu, sr, si, center = map(Q, (mu, shift_real, shift_imag, center))
    if mu <= 0 or sr <= 0 or not isinstance(frequency, int):
        raise ValueError('mu e Re s positivos e frequencia inteira exigidos')
    coefficients = []
    for source in (real, imag):
        terms = {}
        for (a, b, k), value in Poly(source).terms.items():
            if a or k:
                raise ValueError('o perfil do modo deve depender somente de xi')
            for n in range(b+1):
                terms[n] = terms.get(n, Q(0))+value*comb(b,n)*center**(b-n)
        coefficients.append(terms)
    yr, yi = Poly(), Poly()
    for n in coefficients[0].keys() | coefficients[1].keys():
        fr, fi = (terms.get(n, Q(0)) for terms in coefficients)
        d, b = sr+mu*(n-1), si+Q(frequency, 2)
        denominator = d*d+b*b
        if not denominator:
            if fr or fi:
                raise ValueError('obstrucao Floquet no traco do cone')
            continue
        basis = power(X-center, n)
        yr += ((d*fr+b*fi)/denominator)*basis
        yi += ((d*fi-b*fr)/denominator)*basis
    return yr, yi
