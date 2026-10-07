#!/usr/bin/env python3
"""Perturbacoes NAO esfericas axiais (impares), l >= 2, do espaco de Choptuik (fundo RefA de RT): exploracao
numerica (NAO rigorosa) do espectro de Floquet. Fase 1 da analise da simetria esferica (docs/NAO_ESFERICO.md).

Equacao de Gerlach-Sengupta (Martin-Garcia & Gundlach 1999, eq. 57, sem fonte para campo escalar):
    [r^-2 (r^4 Pi)^{|A}]_{|A} = (l - 1)(l + 2) Pi.
No fundo de RT (g^AB = e^{2 zeta} (D+ D- + D- D+)/2, sqrt|g| = e^{-2 zeta}/mu, D^s log r = mu/xi - w2^s), com G = r^2 Pi:
    (D+ + mu)(D-G + 2 v- G) + (D- - mu)(D+G + 2 v+ G) = 2 (l-1)(l+2) (W/xi)^2 G,   v^s = mu/xi - w2^s,
    W = mu + xi (w1 - w2- - w2+)/2.   Codificacao de RT: w2- = w2(tau, xi), w2+ = -w2(tau, -xi), xi em [-1, 1].
Floquet G = e^{s tau} g, g 4pi-periodica: D^s -> D^s + s s. Problema quadratico K0 + s K1 - 2 s^2 = 0 com
    K0 = D+D- + D-D+ + mu (D- - D+) + 2 v- D+ + 2 v+ D- + 2 (D+v- + D-v+) + 2 mu (v- - v+) - 2(l-1)(l+2)(W/xi)^2,
    K1 = 2 (D- - D+) - 2 mu + 2 (v- - v+),   D^t v^s = -mu^2 (1 + t xi)/xi^2 - D^t w2^s.
Colocacao: Fourier em tau (N_tau pontos em [0, 4pi)), Chebyshev-Gauss em xi (N_xi par: evita xi = 0).
Conversao para Martin-Garcia & Gundlach: lambda = 2 pi s / K, K = Delta/2 = 1,7227..."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys
import time

import numpy as np
from numpy.polynomial import chebyshev as C
import scipy.linalg as sla

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gauge_refB import ler

K_RT = 1.7227262011139107


def campos_refA(path='.cache/rt-1203.3766v1/sourcecode/RefA.dat'):
    campos, mu = ler(path)
    out = []
    for h, pares in campos:
        esc = 2.0**h['TwoExp']
        V = np.zeros((h['num_m'], h['num_n']), complex)
        for i in range(h['num_m']):
            for k in range(h['num_n']):
                re, im = pares[i*h['num_n'] + k]
                V[i, k] = complex(int(re), int(im))*esc
        assert h['off_m'] == 0 and h['off_n'] == 0
        out.append(V)
    return out, float(mu)


def avalia(V, tau, xi, dtau=0, dxi=0):
    """f(tau, xi) = sum_{m>=0} w_m Re(e^{i m tau/2} c_m(xi)), c_m = v_m0 + 2 sum_{n>=1} v_mn T_n(xi)."""
    nm, nn = V.shape
    a = V.copy(); a[:, 1:] *= 2
    if dxi:
        a = np.array([C.chebder(a[i], dxi) for i in range(nm)])
    cm = np.array([C.chebval(xi, a[i]) for i in range(nm)])          # (nm, len(xi))
    m = np.arange(nm)
    fase = np.exp(1j*np.outer(tau, m)/2)*((1j*m/2)**dtau)[None, :]      # (len(tau), nm)
    w = np.where(m == 0, 1.0, 2.0)
    return np.real((fase*w[None, :]) @ cm)                              # (len(tau), len(xi))


def matrizes(Nt, Nx, xs=1.0):
    tau = 4*np.pi*np.arange(Nt)/Nt
    k = np.fft.fftfreq(Nt, 1.0/Nt)*(2*np.pi/(4*np.pi))                  # frequencias m/2
    F = np.fft.fft(np.eye(Nt), axis=0)
    Dt = np.real_if_close(np.fft.ifft(1j*k[:, None]*F, axis=0))
    assert Nt % 2 == 1, 'N_tau impar (sem modo de Nyquist)'
    Dt = Dt.real
    j = np.arange(Nx)
    y = np.cos(np.pi*(2*j + 1)/(2*Nx))                                  # Chebyshev-Gauss em y = xi/xs, Nx par
    xi = xs*y
    Vd = C.chebvander(y, Nx - 1)
    Vinv = np.linalg.inv(Vd)
    Dc = np.array([C.chebval(y, C.chebder(np.eye(Nx)[i])) for i in range(Nx)]).T
    Dx = (Dc @ Vinv)/xs
    return tau, xi, Dt, Dx


def problema(l, Nt, Nx, campos, mu, xs=1.0, devolve_base=False):
    tau, xi, Dt, Dx = matrizes(Nt, Nx, xs)
    w1, w2 = campos[0], campos[1]
    T, X = np.meshgrid(tau, xi, indexing='ij')
    w2m = avalia(w2, tau, xi); w2p = -avalia(w2, tau, -xi)
    W = mu + X*(avalia(w1, tau, xi) - w2m - w2p)/2
    # D^+ w2^- e D^- w2^+ = (D^+ w2^-)(tau, -xi)
    Dpw2m = avalia(w2, tau, xi, dtau=1) + mu*(1 + X)*avalia(w2, tau, xi, dxi=1)
    Dmw2p = avalia(w2, tau, -xi, dtau=1) + mu*(1 - X)*avalia(w2, tau, -xi, dxi=1)
    vm = mu/X - w2m; vp = mu/X - w2p
    Dpvm = -mu**2*(1 + X)/X**2 - Dpw2m
    Dmvp = -mu**2*(1 - X)/X**2 - Dmw2p
    It, Ix = np.eye(Nt), np.eye(Nx)
    DT = np.kron(Dt, Ix); DX = np.kron(It, Dx)
    diag = lambda f: np.diag(f.ravel())
    Dp = DT + diag(mu*(1 + X)) @ DX
    Dm = -DT + diag(mu*(1 - X)) @ DX
    K0 = (Dp @ Dm + Dm @ Dp + mu*(Dm - Dp) + 2*diag(vm) @ Dp + 2*diag(vp) @ Dm
          + diag(2*(Dpvm + Dmvp) + 2*mu*(vm - vp) - 2*(l - 1)*(l + 2)*(W/X)**2))
    K1 = 2*(Dm - Dp) + diag(-2*mu + 2*(vm - vp))
    n = Nt*Nx
    A = np.block([[np.zeros((n, n)), np.eye(n)], [K0, K1]])
    B = np.block([[np.eye(n), np.zeros((n, n))], [np.zeros((n, n)), 2*np.eye(n)]])
    # regularidade no centro: G = xi^l * (polinomio PAR), base T_{2k}; projecao de minimos quadrados
    Pk = np.array([xi**l*C.chebval(xi/xs, np.eye(Nx)[2*k]) for k in range(Nx//2)]).T
    P1 = np.kron(np.eye(Nt), Pk)
    Z = np.zeros_like(P1)
    P = np.block([[P1, Z], [Z, P1]])
    Pinv = np.linalg.pinv(P)
    if devolve_base:
        return Pinv @ A @ P, Pinv @ B @ P, (Nt, Nx//2)
    return Pinv @ A @ P, Pinv @ B @ P


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--l', type=int, default=2)
    ap.add_argument('--res', default='9x24,13x32,17x40')
    ap.add_argument('--top', type=int, default=8)
    ap.add_argument('--xs', type=float, default=1.02, help='extremo do dominio em xi (> 1: cone interior; M de RT vai ate 41/40)')
    ap.add_argument('--suave', type=float, default=1e-3, help='filtro: |coef. de Chebyshev na metade superior| / |todos| <= este valor')
    ap.add_argument('--plano', action='store_true', help='teste: espaco de Minkowski (todos os omega = 0, W = mu)')
    a = ap.parse_args()
    campos, mu = campos_refA()
    if a.plano:
        campos = [0*V for V in campos]
    print(f'mu_RefA = {mu:.12f}, l = {a.l} (axial)')
    for r in a.res.split(','):
        Nt, Nx = map(int, r.split('x'))
        t0 = time.time()
        A, B, (nt, nk) = problema(a.l, Nt, Nx, campos, mu, a.xs, devolve_base=True)
        ev, vec = sla.eig(A, B)
        ok = np.isfinite(ev)
        ev, vec = ev[ok], vec[:, ok]
        # suavidade: coeficientes (tau, k) da componente g; fracao de energia em k >= nk/2
        g = vec[:nt*nk].reshape(nt, nk, -1)
        e_alto = np.sqrt(np.sum(np.abs(g[:, nk//2:])**2, axis=(0, 1)))
        e_tot = np.sqrt(np.sum(np.abs(g)**2, axis=(0, 1)))
        suave = e_alto/np.maximum(e_tot, 1e-300) <= a.suave
        ev = ev[suave]
        ev = ev[np.argsort(-ev.real)]
        # representante de Im s em (-1/4, 1/4] (aliases i/2 no periodo 4pi)
        rep = ev.real + 1j*((ev.imag + 0.25) % 0.5 - 0.25)
        print(f'{Nt}x{Nx}, xs = {a.xs} ({time.time() - t0:.0f}s), {len(ev)} autovalores suaves: maiores Re s (Im mod 1/2):')
        vistos = []
        for z in rep:
            if any(abs(z - v) < 1e-6 for v in vistos):
                continue
            vistos.append(z)
            lam = 2*np.pi*z/K_RT
            print(f'   s = {z.real:+.6f} {z.imag:+.6f}i   kappa*Delta = {lam.real*2*K_RT:+.4f}, '
                  f'omega*Delta/2pi = {lam.imag*2*K_RT/(2*np.pi):+.4f}')
            if len(vistos) >= a.top:
                break


if __name__ == '__main__':
    main()
