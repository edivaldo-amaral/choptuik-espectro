#!/usr/bin/env python3
"""R9 (revisao independente C1): Brauer em racionais.
L~(s) = L(s) + sum_k delta_k x_k phi_k^*, L(s) = L(0) + s, L(0) x_0 = 0, L(mu*) x_1 = 0 (x_k EXATOS).
Entao L(s)^-1 x_0 = x_0/s, L(s)^-1 x_1 = x_1/(s - mu*) e, com e_ij = phi_i^*(x_j - x_j^A), phi_i^* x_j^A = delta_ij:
  D(s) = det(I_2 + [delta_j phi_i^* x_j/(s - s_j)]) = q(s)/(s (s - mu*)),
  q(s) = (s + delta_0 (1 + e_00)) (s - mu* + delta_1 (1 + e_11)) - delta_0 delta_1 e_01 e_10.
Confere: delta = (1, 1 + mu_A) com mu_A = MU = 722873400/2^32 (o dado de schur.json), |e_ij| <= ||phi_i|| ||x_j - x_j^A||,
e min_{Re s >= 0} |q(s)| > 0 (cota inferior racional); e a identidade de D por conta simbolica com matrizes 2x2."""
import json, sys
from fractions import Fraction as Fr
sys.path.insert(0, 'scripts')
from tile_perron import DELTA_MU
import random

info = json.load(open('build/rouche_L_rig/z/lab2/schur.json'))['defl_refA']
MU = Fr(722873400, 2**32)
dl = [Fr(x) for x in info['delta']]
print(f'delta registrado = {info["delta"]};  1 + MU = {float(1 + MU)!r};  delta_1 == 1 + MU (float): {info["delta"][1] == float(1 + MU)}')
nphi = [Fr(x) for x in info['norma_phi']]
E = json.load(open('build/s4/rig/E.json'))
ex = [Fr(1.4563e-09) + Fr(1, 10**60), Fr(E['e'])]          # ||x_0 - x_0^A|| (d_tau RefB + Corr), ||x_1 - x_1^A|| (S4)
eb = [[nphi[i]*ex[j] for j in range(2)] for i in range(2)]  # |e_ij| <= ||phi_i|| ||x_j - x_j^A||
print('cotas |e_ij|:', [[f'{float(x):.2e}' for x in l] for l in eb])
dmu = Fr(DELTA_MU)                                         # |mu* - MU|
# para Re s >= 0: |s + d0(1 + e00)| >= d0 (1 - |e00|);  |s - mu* + d1 (1 + e11)| >= d1 (1 - |e11|) - mu* (Re s >= 0)
# mu* <= MU + dmu, d1 = 1 + MU  =>  >= 1 - d1 |e11| - dmu
a0 = dl[0]*(1 - eb[0][0]); a1 = 1 - dl[1]*eb[1][1] - dmu
cota_q = a0*a1 - dl[0]*dl[1]*eb[0][1]*eb[1][0]
print(f'DELTA_MU = {float(dmu):.2e};  min_(Re s >= 0) |q(s)| >= {float(cota_q):.15f} > 0: {cota_q > 0}')
# polos de D em Omega: s = mu* (dentro: 0 < mu* < 3, Im 0) e s = 0 (no bordo, excluido pela indentacao)
print(f'mu* em [{float(MU - dmu):.12f}, {float(MU + dmu):.12f}] contido em (0, 3): {MU - dmu > 0}')
# identidade de D: conferida EXATAMENTE (racionais) em 200 pontos aleatorios (polinomio de grau baixo)
random.seed(1); maxdif = Fr(0)
for _ in range(200):
    r_ = lambda: Fr(random.randint(-10**6, 10**6), random.randint(1, 10**6))
    s_, m_, d0, d1, e00, e01, e10, e11 = (r_() for _ in range(8))
    if s_ == 0 or s_ == m_:
        continue
    P = [[1 + e00, e01], [e10, 1 + e11]]; dd = [d0, d1]; sj = [Fr(0), m_]
    Mx = [[(1 if i == j else 0) + dd[j]*P[i][j]/(s_ - sj[j]) for j in range(2)] for i in range(2)]
    det = Mx[0][0]*Mx[1][1] - Mx[0][1]*Mx[1][0]
    q = (s_ + d0*(1 + e00))*(s_ - m_ + d1*(1 + e11)) - d0*d1*e01*e10
    maxdif = max(maxdif, abs(s_*(s_ - m_)*det - q))
print('max |s (s - mu) det(I + L^-1 E) - q(s)| em 200 pontos racionais =', maxdif)
print('Indentacao em 0 (Omega_eps exclui 0): M(L~) = M(L) + [zeros(D) - polos(D)] em Omega_eps = M(L) + (0 - 1)  =>  N(L) = N(L~) + 1')
