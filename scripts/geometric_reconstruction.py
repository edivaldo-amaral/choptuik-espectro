"""Identidades locais da reconstrucao geometrica linearizada, mu fixo.

Polinomios em tau sao usados somente para testar algebra diferencial;
a inversa global Floquet e demonstrada separadamente na nota matematica.
"""
from fractions import Fraction as Q
from sharp_propagation import Poly, X


def null_derivative(p, mu, sigma, shift=0):
    return sigma*(p.derivative(0)+shift*p)+mu*(1+sigma*X)*p.derivative(1)


def W_of(w, mu):
    return mu+X*(w[0]-w[1]+w[1].parity())/2


def temporal_form(h):
    """a para d_s z=a d tau+b d xi, com h^+=-P h^-."""
    return (-(1-X)*h.parity()-(1+X)*h)/2


def radial_form(h, mu):
    return (h-h.parity())/(2*mu)


def inverse_tau_polynomial(p, shift):
    shift = Q(shift)
    if shift == 0:
        raise ValueError('s nao nulo exigido')
    result, term, factor = Poly(), Poly(p), 1/shift
    while term != 0:
        result += factor*term
        term = term.derivative(0)
        factor /= -shift
    return result


def reconstruct_polynomial(h, mu, shift):
    numerator = h[0]-h[1]+h[1].parity()
    if numerator.parity() != -numerator:
        raise ValueError('numerador de delta Q deve ser impar')
    q = numerator.div_x()/2
    z = inverse_tau_polynomial(temporal_form(h[2]), shift)
    phi = inverse_tau_polynomial(temporal_form(h[3]), shift)
    return q, z, phi


def reconstruct_polynomial_chain(chain, mu, shift):
    result = []
    previous_z, previous_phi = Poly(), Poly()
    for h in chain:
        q, _, _ = reconstruct_polynomial(h, mu, shift)
        z = inverse_tau_polynomial(temporal_form(h[2])-previous_z, shift)
        phi = inverse_tau_polynomial(temporal_form(h[3])-previous_phi, shift)
        result.append((q, z, phi))
        previous_z, previous_phi = z, phi
    return result
