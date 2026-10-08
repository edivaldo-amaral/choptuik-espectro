#!/usr/bin/env python3
"""R5 (revisao de 08/10/2026): teste numerico do Lema 5.4 (deflacao de Brauer) e da formula de multiplicidade,
num modelo de matrizes; e conferencia da cobertura do circulo do Rouche de K (16 tiles).

Lema 5.4: L(s) = A + s, L(s) x0 = s x0 (A x0 = 0), L(s) x1 = (s - mu) x1 (A x1 = -mu x1),
L~(s) = L(s) + d0 x0 phi0^* + d1 x1 phi1^*. Afirma: det(I + L^-1 (L~ - L)) = q(s)/(s (s - mu)) com
q = det(diag(s, s - mu) + [d_j phi_i^* x_j]_ij) e mult_L~(s0) = mult_L(s0) + ord_s0 q/(s(s - mu)).
"""
import json
from fractions import Fraction as Fr
import math
from pathlib import Path
import numpy as np

rng = np.random.default_rng(7)
RAIZ = Path(__file__).resolve().parents[2]


def mults(M, tol=1e-6):
    """autovalores agrupados (multiplicidade algebrica)."""
    ev = np.linalg.eigvals(M)
    grupos = []
    for z in ev:
        for g in grupos:
            if abs(g[0] - z) < tol:
                g[1] += 1; break
        else:
            grupos.append([z, 1])
    return grupos


def teste(mu, outros, d0, d1, phi_aleat=True, jordan0=False):
    n = 2 + len(outros) + (1 if jordan0 else 0)
    Jm = np.zeros((n, n), complex)
    Jm[0, 0] = 0.0; Jm[1, 1] = -mu
    k = 2
    if jordan0:
        Jm[2, 2] = 0.0; Jm[0, 2] = 1.0; k = 3          # bloco de Jordan 2 em 0 (x0 = autovetor)
    for i, z in enumerate(outros):
        Jm[k + i, k + i] = z
    P = rng.normal(size=(n, n)) + 1j*rng.normal(size=(n, n))
    A = P @ Jm @ np.linalg.inv(P)
    x0, x1 = P[:, 0], P[:, 1]
    assert np.allclose(A @ x0, 0) and np.allclose(A @ x1, -mu*x1)
    phi = [rng.normal(size=n) + 1j*rng.normal(size=n) for _ in range(2)]
    At = A + d0*np.outer(x0, phi[0].conj()) + d1*np.outer(x1, phi[1].conj())
    # (1) identidade do determinante em pontos aleatorios
    err = 0
    for s in (0.3 + 0.2j, -1.7 + 0.1j, 2.2 - 0.9j):
        L = A + s*np.eye(n); Lt = At + s*np.eye(n)
        lhs = np.linalg.det(np.eye(n) + np.linalg.solve(L, Lt - L))
        Mq = np.diag([s, s - mu]) + np.array([[d0*phi[0].conj() @ x0, d1*phi[0].conj() @ x1],
                                              [d0*phi[1].conj() @ x0, d1*phi[1].conj() @ x1]])
        err = max(err, abs(lhs - np.linalg.det(Mq)/(s*(s - mu))))
    # (2) multiplicidades: raizes de L~ = raizes de L - {0, mu} (uma vez cada) + zeros de q
    M = np.array([[d0*phi[0].conj() @ x0, d1*phi[0].conj() @ x1], [d0*phi[1].conj() @ x0, d1*phi[1].conj() @ x1]])
    # q(s) = det(diag(s, s-mu) + M) = s^2 + (M00 + M11 - mu) s + (M00 (M11 - mu) - M01 M10)
    zq = np.roots([1, M[0, 0] + M[1, 1] - mu, M[0, 0]*(M[1, 1] - mu) - M[0, 1]*M[1, 0]])
    raizes_L = [-z for z in np.linalg.eigvals(A)]
    esperado = list(raizes_L)
    for r0 in (0.0, mu):
        i = int(np.argmin([abs(z - r0) for z in esperado])); esperado.pop(i)
    esperado += list(zq)
    obtido = [-z for z in np.linalg.eigvals(At)]
    esperado.sort(key=lambda z: (round(z.real, 5), round(z.imag, 5)))
    obtido.sort(key=lambda z: (round(z.real, 5), round(z.imag, 5)))
    dist = max(min(abs(a - b) for b in obtido) for a in esperado)
    return err, dist, zq


print('Lema 5.4 (modelo de matrizes):')
for args in ((0.1683, [-0.4, 0.73, 1.1 + 0.3j], 1.0, 1.1683),
             (0.1683, [-0.4, 0.73, 1.1 + 0.3j], 1.0, 1.1683, True, True),        # 0 com multiplicidade 2
             (0.1683, [-0.4, -1.0, 0.5], 1.0, 1.1683)):                          # outro autovalor de A em -1 (= zero de q?)
    err, dist, zq = teste(*args)
    print(f'  outros={args[1]} jordan0={len(args) > 5}: |det - q/(s(s-mu))| <= {err:.1e}; '
          f'raizes(L~) vs raizes(L) - {{0,mu}} + zeros(q): {dist:.1e}; zeros de q {np.round(zq, 4)}')
    assert err < 1e-9 and dist < 1e-6

# ---- Rouche de K: os 16 discos cobrem o circulo |s - 0.401| = 0.05?
d = json.loads((RAIZ/'build/s3b/rouche_K_circ.json').read_text())
rZ = Fr(d['rZ']); rho = Fr(d['rho'])
cent = [(Fr(t['dV'][0]), Fr(t['dV'][1])) for t in d['tiles']]
ang = sorted(math.atan2(float(y), float(x)) % (2*math.pi) for x, y in cent)
gaps = [(ang[(i + 1) % 16] - ang[i]) % (2*math.pi) for i in range(16)]
# pior ponto do circulo: no meio do maior arco; distancia corda = 2 rho sin(gap/4)... (meio do arco ate a ponta)
pior = max(2*float(rho)*math.sin(g/4) for g in gaps)
# conferencia racional grosseira: amostra 20000 pontos do circulo e mede a distancia ao centro de tile mais proximo
pts = [(math.cos(2*math.pi*k/20000)*0.05, math.sin(2*math.pi*k/20000)*0.05) for k in range(20000)]
dmax = max(min(math.hypot(px - float(x), py - float(y)) for x, y in cent) for px, py in pts)
print(f'Rouche de K: rZ = {float(rZ):.12f}; pior distancia (analitica, meio do arco) = {pior:.12f}; '
      f'amostrada = {dmax:.12f}; folga = {float(rZ) - pior:.2e}; raio dos centros = '
      f'{max(abs(math.hypot(float(x), float(y)) - 0.05) for x, y in cent):.1e} de 0.05')
assert pior < float(rZ)
print('OK')
