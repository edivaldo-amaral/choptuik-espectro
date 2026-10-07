#!/usr/bin/env python3
"""Revisao V (03/10/2026), ilustracao NUMERICA (nao e prova) da volta nos autovetores de Galerkin.

Para cada raiz real s de L em Re s > 0 (mu*, s_K, 0,7332), com o autovetor truncado h de L(s) (fundo RefA):
  q = (h1 - h2 + P h2)/(2 xi),  z = T_s^-1 a3,  p = T_s^-1 a4   (como no documento)
e mede, numa malha (tau, xi) em [0, 4pi) x [-1, 1]:
  e3  = |D_s^- z - h3| e |D_s^+ z - h3^+|   (so precisa de L: dOmega_3^+ = 0 => compatibilidade (3,3))
  e2  = |D_s^- (z + xi^2 q/W*) - h2|        (precisa de C: Omega1 impar e Omega2# entram aqui)
  |C(s)h| relativo (constraint do autovetor truncado).
Esperado: em mu* e 0,7332 (constraints satisfeitas) e2 pequeno; em s_K (viola constraints) e2 = O(1).
Rodar da raiz do repositorio.
"""
import sys
from pathlib import Path
import numpy as np
from numpy.polynomial import chebyshev as Ch

sys.path.insert(0, str(Path('scripts').resolve()))
import signed_operator_L as sl  # noqa: E402
import signed_constraints as sc  # noqa: E402

MU = sl.MU
M, N = int(sys.argv[1]) if len(sys.argv) > 1 else 8, int(sys.argv[2]) if len(sys.argv) > 2 else 24
g = sl.campos()
L0, rs, offs = sl.operador(M, N, 0.0, g=g)
ev, V = np.linalg.eig(L0)
Ccal, sos, ris, oo, oi = sc.operador_C(M, N, 0.0, g=g)
_, _, _, _, _ = Ccal, sos, ris, oo, oi
C1, *_ = sc.operador_C(M, N, 1.0, g=g)
C1 = C1 - Ccal                                     # C(s) = C0 + s C1

Nt, Nx = 256, 64
tau = 4*np.pi*np.arange(Nt)/Nt
xk = np.cos(np.pi*(np.arange(Nx) + 0.5)/Nx)          # nos de Chebyshev (sem xi = 0)
TT, XX = np.meshgrid(tau, xk, indexing='ij')


def campo(coefs_por_m, xi):
    out = np.zeros((Nt, len(xi)), complex)
    for m, c in coefs_por_m.items():
        w = np.full(len(c), 2.0); w[0] = 1.0
        out += np.exp(0.5j*m*tau)[:, None]*Ch.chebval(xi, c*w)[None, :]
    return out


def fundo(j, xi):
    d = {}
    for m, a in g[j].items():
        d[m] = a
        if m > 0:
            d[-m] = a.conj()
    return campo(d, xi)


Wst = MU + XX*(fundo(0, xk) - fundo(1, xk) + fundo(1, -xk))/2
kfreq = np.fft.fftfreq(Nt, d=1.0/Nt)*0.5            # e^{i m tau/2}: m = fftfreq*Nt


def dtau(f):
    return np.fft.ifft(1j*kfreq[:, None]*np.fft.fft(f, axis=0), axis=0)


def dxi(f):
    out = np.zeros_like(f)
    for j in range(Nt):
        cr = Ch.chebfit(xk, f[j].real, Nx - 1); ci = Ch.chebfit(xk, f[j].imag, Nx - 1)
        out[j] = Ch.chebval(xk, Ch.chebder(cr)) + 1j*Ch.chebval(xk, Ch.chebder(ci))
    return out


def Ds(sig, f, s):
    return sig*(dtau(f) + s*f) + MU*(1 + sig*XX)*dxi(f)


def analisa(alvo):
    i = np.argmin(np.abs(ev + alvo))
    s = -ev[i]
    h = V[:, i]
    k = np.argmax(np.abs(h)); h = h/h[k]
    comps = {c: {} for c in range(4)}
    for r, o in zip(rs, offs[:-1]):
        for c, nn in r.comps:
            a0, a1 = r.off[c]
            full = np.zeros(N, complex); full[nn] = h[o + a0:o + a1]
            comps[c][r.m] = full
    H = [campo(comps[c], xk) for c in range(4)]
    HP = [campo(comps[c], -xk) for c in range(4)]
    q = (H[0] - H[1] + HP[1])/(2*XX)
    zp = {}
    for nome, c in (('z', 2), ('p', 3)):
        # T_s^-1 modo a modo: a = [(1-xi) h^+ - (1+xi) h^-]/2, h^+ = -P h
        a = ((1 - XX)*(-HP[c]) - (1 + XX)*H[c])/2
        A = np.fft.fft(a, axis=0)
        zp[nome] = np.fft.ifft(A/(s + 1j*kfreq[:, None]), axis=0)
    z, p = zp['z'], zp['p']
    rel = lambda e, ref: np.max(np.abs(e))/np.max(np.abs(ref))
    e3m = rel(Ds(-1, z, s) - H[2], H[2]); e3p = rel(Ds(1, z, s) - (-HP[2]), H[2])
    e4m = rel(Ds(-1, p, s) - H[3], H[3])
    e2 = rel(Ds(-1, z + XX**2*q/Wst, s) - H[1], H[1])
    Ch_ = (Ccal + s*C1) @ h
    cres = np.linalg.norm(Ch_)/np.linalg.norm(h)
    lres = np.linalg.norm((L0 + s*np.eye(len(L0))) @ h)/np.linalg.norm(h)
    print(f's = {s.real:.6f}{s.imag:+.1e}i  |L h|/|h| = {lres:.1e}  |C(s)h|/|h| = {cres:.2e}  '
          f'e3- = {e3m:.1e}  e3+ = {e3p:.1e}  e4- = {e4m:.1e}  e2 (log W) = {e2:.2e}')
    return dict(s=complex(s), C=cres, e3m=e3m, e3p=e3p, e4m=e4m, e2=e2)


print(f'Galerkin M = {M}, N = {N}, dim {len(L0)}; malha {Nt} x {Nx}')
res = {nome: analisa(alvo) for nome, alvo in (('mu*', 0.16831), ('s_K', 0.40102), ('0,7332', 0.73318))}
ok = (res['0,7332']['e2'] < 1e-2 and res['mu*']['e2'] < 1e-2 and res['s_K']['e2'] > 0.1)
print('PADRAO ESPERADO' if ok else 'PADRAO INESPERADO')
