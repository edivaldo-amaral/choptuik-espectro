"""C3 (revisao CRH, 03/10/2026): H-reg' (HREC_IDA.md, R5 e par. 3), peso de Fourier kappa1' = 1.

Rodar da raiz do projeto:  OPENBLAS_NUM_THREADS=2 .venv/bin/python .codex-runs/2026-10-03-revisao-CRH/c3_hreg.py

  1. free_tail (sharp_exterior.py) nao tem parametro de peso de Fourier: assinatura e valores.
  2. Os majorantes de B (component_majorant / sharp_majorant, normas de campo de RefA) sao monotonos no
     peso de Fourier: recalculados com fourier = 1 e com o peso original; todos os indices m de RefA sao >= 0
     (senao fourier**m < 1 quebraria a monotonia).
  3. Convolucao em m com pesos kappa^|m| em l2: norma em kappa = 1 <= sum |b_k| <= sum kappa^|k| |b_k|.
  4. Um perfil C-infinito NAO analitico em tau, analitico em xi: h = exp(-1/sin^2(tau/4)) / (2 - xi).
     Coeficientes na base de RT e^{i m tau/2 + i n theta}, xi = cos theta. Mostra: decaimento em m mais rapido
     que qualquer potencia mas NAO geometrico (~exp(-C m^(2/3))) (logo h esta em Y_(1,k2') e NAO em Y_(k1,k2') para k1 > 1), e
     decaimento geometrico em n com razao 1/rho, rho = 2 + sqrt(3) (elipse de Bernstein ate o polo xi = 2).
     Tambem J h ~ (i m/2 + mu (n + 1)) h_mn continua somavel com peso (1, k2').
"""
from fractions import Fraction as Q
import inspect
import math
import sys
from pathlib import Path

import numpy as np

RAIZ = Path('<repo>')
sys.path.insert(0, str(RAIZ/'scripts'))
from sharp_exterior import free_tail                                      # noqa: E402
from component_exterior import read_fields, field_expression, weighted_bound  # noqa: E402
from sharp_propagation import table                                       # noqa: E402

# 1.
sig = inspect.signature(free_tail)
print('1. free_tail', sig, '-> parametros:', list(sig.parameters))
for (M, N, k2) in ((8, 24, Q(9, 8)), (64, 192, Q(9, 8)), (256, 1536, Q(1) + Q(1, 64))):
    print(f'   free_tail({M}, {N}, radial={k2}) = {float(free_tail(M, N, radial=k2)):.6f}')

# 2.
fields = read_fields(RAIZ/'.cache/rt-1203.3766v1/sourcecode/RefA.dat')
ms = [m for f in fields for (m, n, _) in f]
ns = [n for f in fields for (m, n, _) in f]
print(f'2. RefA: m em [{min(ms)}, {max(ms)}], n em [{min(ns)}, {max(ns)}]')


def fnorm(field, parity=None, fourier=Q(65, 64), radial=Q(5, 4)):
    # copia de component_exterior.field_norm com o peso de Fourier como parametro
    return sum((abs(v)*(2-(m == 0))*(2-(n == 0))*fourier**m*radial**n
                for (m, n, _), v in field.items() if parity is None or n % 2 == parity), Q(0))


def comp_major(fourier, radial):
    rows = table('GAMMA2_macro.c', '_G2PAIR_')
    selected = (0, 1, 3, 4)
    res = [[Q(0)]*4 for _ in range(4)]
    for out, row in enumerate(selected):
        for inp in range(4):
            cols = [(0, 1, Q(1, 2))] if inp == 0 else [(2*inp-1, 0, Q(1)), (2*inp, 1, Q(1))]
            res[out][inp] = max(fac*fnorm(field_expression(rows.get((row, k), ''), fields),
                                          1-par if out == 0 else None, fourier, radial)
                                for k, par, fac in cols)
    return res


def sharp_major(fourier, radial):
    rows = table('GAMMA2_SHARP_macro.c', '_G2SHARPPAIR_')
    nm = lambda e: fnorm(field_expression(e, fields), None, fourier, radial)
    return [[nm(rows[i, 1])/2, max(nm(rows[i, 4]), nm(rows[i, 5]))/2] for i in range(2)]


ok2 = True
for nome, f, (kf, kr) in (('component_majorant (65/64, 5/4)', comp_major, (Q(65, 64), Q(5, 4))),
                          ('sharp_majorant (129/128, 9/8)', sharp_major, (Q(129, 128), Q(9, 8)))):
    A = f(kf, kr)
    B1 = f(Q(1), kr)
    mono = all(B1[i][j] <= A[i][j] for i in range(len(A)) for j in range(len(A[0])))
    ok2 = ok2 and mono
    razao = max(float(B1[i][j]/A[i][j]) for i in range(len(A)) for j in range(len(A[0])) if A[i][j])
    wA = weighted_bound(A, [Q(1)]*len(A)); wB = weighted_bound(B1, [Q(1)]*len(A))
    print(f'   {nome}: entradas com fourier = 1 <= originais: {mono}; maior razao {razao:.6f}; '
          f'cota de coluna {float(wB):.6f} <= {float(wA):.6f}')
print('   monotonia no peso de Fourier:', 'OK' if ok2 else 'FALHA')

# 3.
rng = np.random.default_rng(3)
K = 12
b = (rng.normal(size=2*K+1) + 1j*rng.normal(size=2*K+1))*0.6**np.abs(np.arange(-K, K+1))
Mx = 200
ks = np.arange(-K, K+1)
idx = np.arange(-Mx, Mx+1)
C = np.zeros((len(idx), len(idx)), complex)
for a, k in enumerate(ks):
    C += b[a]*np.eye(len(idx), k=-k)
for kap in (1.0, 65/64, 1.2):
    W = np.diag(kap**np.abs(idx))
    nrm = np.linalg.norm(W @ C @ np.linalg.inv(W), 2)
    print(f'3. convolucao em l2 com peso {kap}^|m|: norma {nrm:.6f} <= sum kappa^|k| |b_k| = '
          f'{np.sum(kap**np.abs(ks)*np.abs(b)):.6f}')
print(f'   sum |b_k| = {np.sum(np.abs(b)):.6f}')

# 4.
Nt, Nth = 4096, 256
tau = 4*np.pi*np.arange(Nt)/Nt
th = 2*np.pi*np.arange(Nth)/Nth
with np.errstate(divide='ignore', over='ignore'):
    s2 = np.sin(tau/4)**2
    g = np.where(s2 > 1e-300, np.exp(-1/np.maximum(s2, 1e-300)), 0.0)
f = 1/(2 - np.cos(th))
H = np.outer(g, f)
c = np.fft.fft2(H)/(Nt*Nth)        # coeficiente de e^{i m tau/2} e^{i n theta}: indice m em [0, Nt)
cm = np.abs(np.fft.fft(g)/Nt)      # parte em tau
cn = np.abs(np.fft.fft(f)/Nth)     # parte em theta
print('4. h = exp(-1/sin^2(tau/4))/(2 - xi)  (precisao dupla: so coeficientes acima de ~1e-15 sao usados)')
for m in (10, 20, 40, 80, 100):
    print(f'   |g_m| m = {m:4d}: {cm[m]:.3e}   m^6 |g_m| = {m**6*cm[m]:.3e}   -log|g_m|/m^(2/3) = {-math.log(cm[m])/m**(2/3):.3f}'
          f'   -log|g_m|/m = {-math.log(cm[m])/m:.4f}')
C = min(-math.log(cm[m])/m**(2/3) for m in (40, 80, 100))
print(f'   -log|g_m|/m^(2/3) ~ {C:.2f}, quase constante, e -log|g_m|/m -> 0: decaimento superpolinomial e NAO geometrico')
print('   (Gevrey 3/2: o ponto de sela de exp(-16/tau^2 - i m tau/2) da exp(-C m^(2/3))). kappa1^m |g_m| cresce para')
print('   m > (C/log kappa1)^3:')
for k1 in (1.01, 65/64, 129/128):
    print(f'      kappa1 = {k1:.6f}: m* ~ {(C/math.log(k1))**3:.2e}; a serie sum kappa1^m |g_m| diverge. Com kappa1 = 1 converge.')
rho = 2 + math.sqrt(3)
for n in (5, 10, 20, 28):
    print(f'   |f_n| n = {n:3d}: {cn[n]:.3e}   rho^n |f_n| = {rho**n*cn[n]:.4f}  (rho = 2 + sqrt 3 = 3,732)')
mu = 0.168
k2 = 1.5
nn = np.arange(0, 29)
mm = np.arange(0, 101)
Sh = sum(cm[m]*cn[n]*k2**n for m in mm for n in nn)
SJ = sum(cm[m]*cn[n]*k2**n*abs(1j*m/2 + mu*(n+1)) for m in mm for n in nn)
print(f'   soma truncada (m <= 100, n <= 28) de |h_mn| k2^n, k2 = {k2}: {Sh:.6f};  com o fator de J |i m/2 + mu(n+1)|: {SJ:.6f}')
