"""Identidade off-shell RT, em polinomios racionais locais (tau,xi,s).

Le as tabelas bilineares dos macros RT, sem executar codigo C. Os testes
verificam identidades locais, nao regularidade/Floquet no dominio infinito.
"""
from fractions import Fraction as Q
from pathlib import Path
import re

SOURCE = Path(__file__).resolve().parents[1] / '.cache/rt-1203.3766v1/sourcecode'


class Poly:
    def __init__(self, terms=0):
        if isinstance(terms, Poly):
            terms = terms.terms
        elif not isinstance(terms, dict):
            terms = {(0, 0, 0): Q(terms)}
        self.terms = {k: Q(v) for k, v in terms.items() if v}

    def __add__(self, other):
        result = dict(self.terms)
        for k, v in Poly(other).terms.items():
            result[k] = result.get(k, Q(0)) + v
        return Poly(result)

    __radd__ = __add__

    def __neg__(self):
        return Poly({k: -v for k, v in self.terms.items()})

    def __sub__(self, other):
        return self + -Poly(other)

    def __rsub__(self, other):
        return Poly(other) + -self

    def __mul__(self, other):
        result = {}
        for k, v in self.terms.items():
            for l, w in Poly(other).terms.items():
                key = tuple(a+b for a, b in zip(k, l))
                result[key] = result.get(key, Q(0)) + v*w
        return Poly(result)

    __rmul__ = __mul__

    def __truediv__(self, scalar):
        return self * (1/Q(scalar))

    def __eq__(self, other):
        return self.terms == Poly(other).terms

    def derivative(self, axis):
        return Poly({tuple(a-(i == axis) for i, a in enumerate(k)): v*k[axis]
                     for k, v in self.terms.items() if k[axis]})

    def parity(self):
        return Poly({k: v*(-1)**k[1] for k, v in self.terms.items()})

    def div_x(self):
        # (p(tau,xi,s)-p(tau,0,s))/xi, exatamente a regularizacao RT.
        return Poly({(a, b-1, c): v for (a, b, c), v in self.terms.items() if b})


TAU = Poly({(1, 0, 0): 1})
X = Poly({(0, 1, 0): 1})
SHIFT = Poly({(0, 0, 1): 1})


def odd(p):
    return (p-p.parity())/2


def even(p):
    return (p+p.parity())/2


def j(p, mu, shift=0):
    return p.derivative(0) + mu*((1+X)*p.derivative(1)+p) + shift*p


def table(name, pair):
    source = (SOURCE / name).read_text()
    pattern = rf'case {pair}\((\d+),(\d+)\): ([^\n]+?)break;'
    rows = {(int(n), int(k)): expr for n, k, expr in re.findall(pattern, source)}
    expected = 20 if pair == '_G2SHARPPAIR_' else 25
    if len(rows) != expected:
        raise ValueError(f'tabela RT inesperada: {len(rows)} != {expected}')
    return rows


def linear_expression(expr, fields):
    value = Poly()
    tokens = re.findall(r'A\(t,v([1-4]),([+-]?\d+)\)|P\(t\)', expr)
    for field, coefficient in tokens:
        value = value + int(coefficient)*fields[int(field)-1] if field else value.parity()
    return value


def gamma2(a, b):
    right = [b[0], b[1]+b[1].parity(), b[1]-b[1].parity(),
             b[2]+b[2].parity(), b[2]-b[2].parity(),
             b[3]+b[3].parity(), b[3]-b[3].parity()]
    rows = table('GAMMA2_macro.c', '_G2PAIR_')
    return [sum((linear_expression(expr, a)*right[k]/4
                 for (i, k), expr in rows.items() if i == n), Poly())
            for n in range(5)]


def gamma_sharp_full(w, residual):
    right = [term for p in residual for term in (p+p.parity(), p-p.parity())]
    rows = table('GAMMA2_SHARP_macro.c', '_G2SHARPPAIR_')
    return [sum((linear_expression(rows[n, k], w)*right[k]/4
                 for k in range(10)), Poly()) for n in range(2)]


def gamma1(w):
    a, b, c, d = w
    return [2*a, a+b-b.parity(), 2*(c-b).parity(), -a, d-d.parity()]


def omega(w, mu, shift=0):
    a, b, c, d = w
    principal = [j(a, mu, shift), j(b, mu, shift),
                 -j(b.parity(), mu, shift), j(c, mu, shift), j(d, mu, shift)]
    return [X*p+mu*g+X*q for p, g, q in zip(principal, gamma1(w), gamma2(w, w))]


def linearized_omega(w, h, mu, shift=SHIFT):
    # omega(h) contem o termo quadratico de h, que e retirado aqui.
    hh, wh, hw = gamma2(h, h), gamma2(w, h), gamma2(h, w)
    return [p-X*q+X*(a+b) for p, q, a, b in zip(omega(h, mu, shift), hh, wh, hw)]


def split_plus(q):
    selected = [even(q[0]).div_x(), q[1].div_x(), q[3].div_x(), q[4].div_x()]
    constraints = [odd(q[0]), -q[2].parity()]
    return selected, constraints


def inject_constraints(c):
    return [c[0], Poly(), c[1], Poly(), Poly()]


def lift_selected(f):
    return [X*f[0].parity(), X*f[1].parity(), Poly(),
            X*f[2].parity(), X*f[3].parity()]


def complete(w, e, mu, shift=0):
    a, b, c, d, _ = e
    nonlinear = gamma_sharp_full(w, e)
    return [X*(odd(j(a, mu, shift))+nonlinear[0])+mu*even(a),
            X*(j(c, mu, shift)-j(b, mu, shift).parity()+nonlinear[1])
            + mu*(a+b+b.parity()-2*d.parity())]


def propagation(w, c, mu, shift=SHIFT):
    return [p.div_x() for p in complete(w, inject_constraints(c), mu, shift)]


def forcing(w, f, mu, shift=SHIFT):
    return [-p.div_x() for p in complete(w, lift_selected(f), mu, shift)]
