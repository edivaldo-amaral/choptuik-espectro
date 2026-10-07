#!/usr/bin/env python3
"""Revisao V (03/10/2026), V4: conferencia NUMERICA (nao rigorosa) de W* = mu + xi(omega1 - omega2 + P omega2)/2 > 0
no cilindro real [0, 4 pi) x [-41/40, 41/40], com os dados RefA de RT, e das paridades do fundo.
RT afirmam (d) com |TF(omega - omega_RefA)| <= 2^-24 no cilindro; aqui so se mede a margem.
Rodar da raiz do repositorio.
"""
import sys
from pathlib import Path
import numpy as np
from numpy.polynomial import chebyshev as C

sys.path.insert(0, str(Path('scripts').resolve()))
from component_exterior import read_fields  # noqa: E402

REF = Path('.cache/rt-1203.3766v1/sourcecode/RefA.dat')
MU = 722873400/2**32
F = read_fields(REF)

ms = sorted({m for f in F for (m, n, p) in f})
assert min(ms) >= 0, 'RefA guarda so m >= 0'
coef = []
for f in F:
    M = max(m for (m, n, p) in f) + 1
    N = max(n for (m, n, p) in f) + 1
    a = np.zeros((M, N), complex)
    for (m, n, p), v in f.items():
        a[m, n] += float(v)*(1j if p else 1)
    coef.append(a)

# paridades: omega1 so n impar; omega1..3 so m par; omega4 so m impar
par = {
    'omega1 so n impar': not np.any(coef[0][:, 0::2]),
    'omega1..3 so m par': all(not np.any(coef[j][1::2, :]) for j in range(3)),
    'omega4 so m impar': not np.any(coef[3][0::2, :]),
}


def TF(a, tau, xi):
    """sum_m e^{i m tau/2} [a_m0 + 2 sum_{n>0} a_mn T_n(xi)], com a_{-m} = conj(a_m)."""
    M, N = a.shape
    w = np.full(N, 2.0); w[0] = 1.0
    fm = np.array([C.chebval(xi, a[m]*w) for m in range(M)])          # (M, len(xi))
    ph = np.exp(0.5j*np.outer(tau, np.arange(M)))                   # (len(tau), M)
    val = ph[:, :1]*fm[:1] + 2*np.real(ph[:, 1:] @ fm[1:])
    return np.real(val)


tau = np.linspace(0, 4*np.pi, 801)[:-1]
xi = np.linspace(-41/40, 41/40, 821)
w1 = TF(coef[0], tau, xi)
w2 = TF(coef[1], tau, xi)
w2P = TF(coef[1], tau, -xi)
W = MU + xi[None, :]*(w1 - w2 + w2P)/2
i, j = np.unravel_index(np.argmin(W), W.shape)
print('paridades do fundo:', par)
print(f'min W* na malha 800x821: {W.min():.6f} em tau={tau[i]:.4f}, xi={xi[j]:.4f}')
print(f'max W* na malha: {W.max():.6f}')
# no intervalo fisico 0 <= xi < 41/40
k = xi >= 0
print(f'min W* em 0 <= xi <= 41/40: {W[:, k].min():.6f}')
# perto do cone xi = 1
c = np.abs(xi - 1) < 0.01
print(f'min W* em |xi - 1| < 0.01: {W[:, c].min():.6f}')
print('erro de RefA (RT): |dW| <= (41/80)*(1+2)*2^-24 =', (41/80)*3*2**-24)
assert all(par.values())
assert W.min() > 0
print('W* > 0 na malha: OK (nao rigoroso; a afirmacao rigorosa e (d) de RT)')
