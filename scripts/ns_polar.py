#!/usr/bin/env python3
"""Perturbacoes polares (pares), l >= 2, do espaco de Choptuik (fundo RefA de RT): exploracao numerica (NAO rigorosa)
do espectro de Floquet. Coeficientes de build/naoesf/polar_coef.py (scripts/ns_polar_deriva.py).

Incognitas: Phi, k (paridade (-1)^l, = xi^l * polinomio par), A (= xi^l * polinomio), B = (-1)^l P A.
Equacoes: onda de Phi, traco de (inveqAB), componente nula EAm de (inveqAeven). Colocacao Fourier (tau, 4 pi, N_tau
impar) x Chebyshev-Gauss em xi em [-xs, xs] (xs > 1: cone interior), projecao de minimos quadrados nas bases
regulares, problema quadratico M0 + s M1 + s^2 M2 linearizado. Filtro de suavidade dos autovetores.
Fundo pelas variaveis de RT: zeta_x = (w3+ + w3-)/(2 mu), zeta_t = w3+ + mu - mu (1 + xi) zeta_x; idem log W com
w2 - w3 e phi com w4; W = mu + xi (w1 - w2- - w2+)/2; w^- = w(tau, xi), w^+ = -w(tau, -xi)."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys
import time

import numpy as np
from numpy.polynomial import chebyshev as C
import scipy.linalg as sla

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ns_axial as na

K_RT = na.K_RT


def derivs_campo(V, tau, xi):
    """valores e derivadas (t, x, tt, tx, xx) de w(tau, xi) e de w^+(tau, xi) = -w(tau, -xi), em grade (Nt, Nx)."""
    av = lambda xx, dt, dx: na.avalia(V, tau, xx, dtau=dt, dxi=dx)
    m = {k: av(xi, *d) for k, d in (('f', (0, 0)), ('t', (1, 0)), ('x', (0, 1)), ('tt', (2, 0)), ('tx', (1, 1)), ('xx', (0, 2)))}
    q = {k: av(-xi, *d) for k, d in (('f', (0, 0)), ('t', (1, 0)), ('x', (0, 1)), ('tt', (2, 0)), ('tx', (1, 1)), ('xx', (0, 2)))}
    p = {'f': -q['f'], 't': -q['t'], 'x': q['x'], 'tt': -q['tt'], 'tx': q['tx'], 'xx': -q['xx']}
    return m, p                                                   # sigma = -, sigma = +


def potencial(am, ap, mu, X):
    """F com D^s F = a^s + s mu (zeta) ou a^s (logW, phi: passar sem o termo): devolve F_x, F_t, F_xx, F_tx, F_tt, e F_xt
    pelo outro caminho (consistencia). am, ap: dicionarios de derivadas de a^- e a^+."""
    Fx = (ap['f'] + am['f'])/(2*mu)
    Fxx = (ap['x'] + am['x'])/(2*mu)
    Fxt = (ap['t'] + am['t'])/(2*mu)
    return Fx, Fxx, Fxt


def fundo(campos, mu, tau, xi):
    X = np.broadcast_to(xi[None, :], (len(tau), len(xi)))
    w1m, w1p = derivs_campo(campos[0], tau, xi)
    w2m, w2p = derivs_campo(campos[1], tau, xi)
    w3m, w3p = derivs_campo(campos[2], tau, xi)
    w4m, w4p = derivs_campo(campos[3], tau, xi)
    out = {}
    # zeta: D^s zeta = w3^s + s mu
    zx, zxx, zxt = potencial(w3m, w3p, mu, X)
    zt = w3p['f'] + mu - mu*(1 + X)*zx
    ztt = w3p['t'] - mu*(1 + X)*zxt
    ztx = w3p['x'] - mu*zx - mu*(1 + X)*zxx
    out.update(z_t=zt, z_x=zx, z_tt=ztt, z_tx=0.5*(ztx + zxt), z_xx=zxx)
    cons_z = np.max(np.abs(ztx - zxt))
    # log W: D^s logW = w2^s - w3^s
    bm = {k: w2m[k] - w3m[k] for k in w2m}; bp = {k: w2p[k] - w3p[k] for k in w2p}
    Lx, Lxx, Lxt = potencial(bm, bp, mu, X)
    Lt = bp['f'] - mu*(1 + X)*Lx
    Ltt = bp['t'] - mu*(1 + X)*Lxt
    Ltx = bp['x'] - mu*Lx - mu*(1 + X)*Lxx
    cons_L = np.max(np.abs(Ltx - Lxt))
    W = mu + X*(w1m['f'] - w2m['f'] - w2p['f'])/2
    Wx_alg = (w1m['f'] - w2m['f'] - w2p['f'])/2 + X*(w1m['x'] - w2m['x'] - w2p['x'])/2
    out.update(W=W, W_t=W*Lt, W_x=W*Lx, W_tt=W*(Ltt + Lt**2), W_tx=W*(0.5*(Ltx + Lxt) + Lt*Lx), W_xx=W*(Lxx + Lx**2))
    cons_W = np.max(np.abs(Wx_alg - W*Lx))
    # phi: D^s phi = w4^s
    px, pxx, pxt = potencial(w4m, w4p, mu, X)
    pt = w4p['f'] - mu*(1 + X)*px
    ptt = w4p['t'] - mu*(1 + X)*pxt
    ptx = w4p['x'] - mu*px - mu*(1 + X)*pxx
    cons_p = np.max(np.abs(ptx - pxt))
    out.update(p_t=pt, p_x=px, p_tt=ptt, p_tx=0.5*(ptx + pxt), p_xx=pxx)
    return out, dict(zeta=cons_z, logW=cons_L, W=cons_W, phi=cons_p)


def coeficientes():
    ns = {}
    exec(Path('build/naoesf/polar_coef.py').read_text(), ns)
    return ns['COEF']


def avalia_coef(expr, env):
    return eval(expr, {'sqrt': np.sqrt, 'pi': np.pi}, env)


def problema(l, Nt, Nx, campos, mu, xs, plano=False):
    tau, xi, Dt, Dx = na.matrizes(Nt, Nx, xs)
    if plano:
        T, Xg = np.meshgrid(tau, xi, indexing='ij')
        z = np.zeros_like(Xg)
        bg = dict(z_t=mu + z, z_x=z, z_tt=z, z_tx=z, z_xx=z, W=mu + z, W_t=z, W_x=z, W_tt=z, W_tx=z, W_xx=z,
                  p_t=z, p_x=z, p_tt=z, p_tx=z, p_xx=z)
        cons = {}
    else:
        bg, cons = fundo(campos, mu, tau, xi)
    Xg = np.broadcast_to(xi[None, :], (Nt, Nx))
    env = dict(bg); env.update(xi=Xg, mu=mu, l=l)
    COEF = coeficientes()
    It, Ix = np.eye(Nt), np.eye(Nx)
    DT = np.kron(Dt, Ix); DX = np.kron(It, Dx)
    R = np.kron(It, np.eye(Nx)[::-1])                        # reflexao xi -> -xi (pontos simetricos)
    n = Nt*Nx
    I = np.eye(n)
    sgn = (-1)**l
    eqs = ['wave', 'trace', 'EAm']
    unk = ['Phi', 'k', 'A']
    M = [np.zeros((3*n, 3*n)) for _ in range(3)]
    for ie, e in enumerate(eqs):
        for (u, i, j), expr in COEF[e].items():
            c = np.asarray(avalia_coef(expr, env), dtype=float)
            c = np.broadcast_to(c, (Nt, Nx)).ravel()
            Dj = np.linalg.matrix_power(DX, j) if j else I
            blocos = []                                     # (potencia de s, matriz)
            if i == 0:
                blocos = [(0, Dj)]
            elif i == 1:
                blocos = [(0, DT @ Dj), (1, Dj)]
            else:
                blocos = [(0, DT @ DT), (1, 2*DT), (2, I)]
            col = unk.index(u) if u != 'B' else 2
            for pw, Mat in blocos:
                op = c[:, None]*Mat
                if u == 'B':
                    op = sgn*(op @ R)
                M[pw][ie*n:(ie + 1)*n, col*n:(col + 1)*n] += op
    # bases regulares
    Peven = np.array([Xg[0]**l*C.chebval(xi/xs, np.eye(Nx)[2*k]) for k in range(Nx//2)]).T
    Pfull = np.array([Xg[0]**l*C.chebval(xi/xs, np.eye(Nx)[k]) for k in range(Nx)]).T
    P = sla.block_diag(np.kron(It, Peven), np.kron(It, Peven), np.kron(It, Pfull))
    Pi = np.linalg.pinv(P)
    Mr = [Pi @ Mk @ P for Mk in M]
    nr = Mr[0].shape[0]
    Aop = np.block([[np.zeros((nr, nr)), np.eye(nr)], [Mr[0], Mr[1]]])
    Bop = np.block([[np.eye(nr), np.zeros((nr, nr))], [np.zeros((nr, nr)), -Mr[2]]])
    dims = (Nt, Nx//2, Nx//2, Nx)
    return Aop, Bop, dims, cons


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--l', type=int, default=2)
    ap.add_argument('--res', default='9x24,13x32')
    ap.add_argument('--xs', type=float, default=1.02)
    ap.add_argument('--top', type=int, default=8)
    ap.add_argument('--suave', type=float, default=1e-3)
    ap.add_argument('--plano', action='store_true')
    a = ap.parse_args()
    campos, mu = na.campos_refA()
    print(f'mu_RefA = {mu:.12f}, l = {a.l} (polar)' + (' [PLANO]' if a.plano else ''))
    for r in a.res.split(','):
        Nt, Nx = map(int, r.split('x'))
        t0 = time.time()
        Aop, Bop, (nt, nh, nh2, nf), cons = problema(a.l, Nt, Nx, campos, mu, a.xs, a.plano)
        if cons:
            print('   consistencia do fundo (d_t d_x por dois caminhos):', {k: f'{v:.1e}' for k, v in cons.items()})
        ev, vec = sla.eig(Aop, Bop)
        ok = np.isfinite(ev) & (np.abs(ev) < 1e3)
        ev, vec = ev[ok], vec[:, ok]
        nr = nt*(nh + nh2 + nf)
        v = vec[:nr]
        partes = [v[:nt*nh].reshape(nt, nh, -1), v[nt*nh:nt*(nh + nh2)].reshape(nt, nh2, -1), v[nt*(nh + nh2):].reshape(nt, nf, -1)]
        e_alto = sum(np.sum(np.abs(p_[:, p_.shape[1]//2:])**2, axis=(0, 1)) for p_ in partes)
        e_tot = sum(np.sum(np.abs(p_)**2, axis=(0, 1)) for p_ in partes)
        suave = np.sqrt(e_alto/np.maximum(e_tot, 1e-300)) <= a.suave
        ev = ev[suave]
        ev = ev[np.argsort(-ev.real)]
        rep = ev.real + 1j*((ev.imag + 0.25) % 0.5 - 0.25)
        print(f'{Nt}x{Nx}, xs = {a.xs} ({time.time() - t0:.0f}s), {len(ev)} suaves: maiores Re s (Im mod 1/2):')
        vistos = []
        for z in rep:
            if any(abs(z - w) < 1e-5 for w in vistos):
                continue
            vistos.append(z)
            lam = 2*np.pi*z/K_RT
            print(f'   s = {z.real:+.6f} {z.imag:+.6f}i   kappa*Delta = {lam.real*2*K_RT:+.4f}, '
                  f'omega*Delta/2pi = {lam.imag*2*K_RT/(2*np.pi):+.4f}')
            if len(vistos) >= a.top:
                break


if __name__ == '__main__':
    main()
