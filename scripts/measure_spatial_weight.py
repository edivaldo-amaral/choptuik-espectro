#!/usr/bin/env python3
# DIAGNOSTICO EM PONTO FLUTUANTE -- NAO E PROVA. Ver S3_CONDICIONAMENTO_SHARP.md 10.4.
"""Teste em ponto flutuante: ||B Q0 P_T|| <= ||B W^-1|| * ||W Q0 P_T||, W = w(phi)."""
import sys; sys.path.insert(0, 'scripts')
import numpy as np
from pathlib import Path
from component_exterior import read_fields
from symbol_forms import forma_K
K1, K2 = 129/128, 9/8
MU = 722873400/2**32
F = read_fields(Path('.cache/rt-1203.3766v1/sourcecode/RefA.dat'))
C = []
for f in F:
    A = np.zeros((41, 101), complex)
    for (m, n, part), v in f.items(): A[m, n] += float(v)*(1j if part else 1)
    C.append(A)
def campo(A, w, z):
    n = np.arange(101); zn = z[..., None]**n + z[..., None]**(-n)
    Am = np.einsum('mn,...n->...m', A, zn) - A[:, 0]; Bm = np.einsum('mn,...n->...m', A.conj(), zn) - A[:, 0].conj()
    m = np.arange(41)
    return (Am*w[..., None]**m).sum(-1) + (Bm[..., 1:]*w[..., None]**(-m[1:])).sum(-1)
NT, NP = 128, 256                        # theta, phi em [0, 2pi)
th = np.linspace(0, 2*np.pi, NT, endpoint=False); ph = np.linspace(0, 2*np.pi, NP, endpoint=False)
TH, PH = np.meshgrid(th, ph, indexing='ij')
wv = K1*np.exp(1j*TH); z = K2*np.exp(1j*PH); xi = (z + 1/z)/2
vz = [campo(A, wv, z) for A in C]; vm = [campo(A, wv, -z) for A in C]
Ms = []
for mu in (1/6, 0.17):
    S, M = forma_K(vz, vm, mu/xi)
    Ms.append(np.moveaxis(np.array(M), (0, 1), (-2, -1)))     # (NT, NP, 4, 3)
def P1(wphi):
    """sup ||M N_w^-1/2||, N_w = diag(w(phi)^2 + w(phi+pi)^2, w(phi)^2, w(phi+pi)^2)."""
    wp = wphi; wq = np.roll(wphi, -NP//2)                        # w(phi), w(phi+pi)
    Nw = np.stack([np.sqrt(wp**2 + wq**2), wp, wq], -1)
    return max(np.linalg.norm(M/Nw[None, :, None, :], 2, axis=(-2, -1)).max() for M in Ms)
# perfil local do fundo: s(phi) = max_theta ||M(theta,phi)|| (sem peso, N padrao)
Nh = np.array([np.sqrt(2), 1, 1])
s_phi = np.maximum(*(np.linalg.norm(M/Nh, 2, axis=(-2, -1)).max(0) for M in Ms))
def livre(ns, passo, sigma):
    J = np.diag(MU*(ns + 1.0)).astype(complex) + sigma*np.eye(len(ns))
    for b in range(len(ns)): J[:b, b] += 2*MU*ns[b]
    return J
NR = 360; PG = 2048
phg = np.linspace(0, 2*np.pi, PG, endpoint=False)
def P2(wfun, ms, s0, radial_min=0):
    """sup_m ||W R_m^-1 P||: saida como funcao no circulo |z|=K2, norma L2 por Parseval na grade."""
    wg = wfun(phg)
    best = 0.0
    for m in ms:
        sig = s0 + 1j*m/2
        for ns, passo in ((np.arange(1, NR, 2), 2), (np.arange(NR), 1)):
            R = np.linalg.inv(livre(ns, passo, sig))
            cols = ns >= radial_min
            # sintese: Y(phi) = sum_{n in Z} y_|n| K2^n e^{i n phi}
            E = (K2**ns[:, None])*np.exp(1j*ns[:, None]*phg[None, :]) + \
                np.where(ns[:, None] > 0, (K2**(-ns[:, None]))*np.exp(-1j*ns[:, None]*phg[None, :]), 0)
            Ymap = (E.T @ R[:, cols])*wg[:, None]/np.sqrt(PG)          # (PG, ncols): valores ponderados
            wi = np.sqrt(np.where(ns[cols] == 0, 1.0, K2**(2.0*ns[cols]) + K2**(-2.0*ns[cols])))
            best = max(best, np.linalg.norm(Ymap/wi[None, :], 2))
    return best
wint = lambda arr: (lambda p: np.interp(p, ph, arr, period=2*np.pi))
print(f'perfil do fundo s(phi): max {s_phi.max():.3f} em phi={ph[s_phi.argmax()]:.2f}, min {s_phi.min():.3f} em phi={ph[s_phi.argmin()]:.2f}')
for M in (20, 30):
    ms = [M, M+4, M+12, M+30]
    for nome, warr in (('sem peso', np.ones(NP)), ('w = s^0.5', s_phi**0.5), ('w = s', s_phi), ('w = s^1.5', s_phi**1.5)):
        warr = warr/warr.max()
        a = P1(warr); b = P2(wint(warr), ms, 0.0)
        print(f'M={M} {nome:10s}: ||B W^-1|| = {a:.4f}  ||W Q0 P_F|| = {b:.4f}  produto = {a*b:.4f}')
print('referencias: soma modo a modo 0,879 (M=20) e 0,626 (M=30); verdadeiro (bloco 20x80) ~0,5')
