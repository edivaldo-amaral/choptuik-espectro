#!/usr/bin/env python3
"""Revisao R2: conferencias simbolicas (sympy, funcoes genericas zeta, Q, phi de (tau, xi), fora das equacoes):

 (1) X0 = d_{u_-} + d_{u_+} = e^{mu tau}(mu^-1 d_tau + xi d_xi)  (regra da cadeia);
 (2) L_{X0} gbar tem componentes ++ e -- nulas (B = kappa L_{X0}(gbar, phi*) esta no ansatz);
 (3) a extracao de RT aplicada a L_{X0}(gbar, phi*) da delta omega_i^sigma = e^{mu tau}(Z + 1) omega_i^sigma,
     Z = mu^-1 d_tau + xi d_xi, para as 7 funcoes (omega1, omega2^+-, omega3^+-, omega4^+-): isto e,
     h0 = kappa e^{(mu - s0) tau} g*_perfil = kappa e^{i m0 tau/2} g*_perfil;
 (4) estrutura secular: para delta' = A + tau B no ansatz, a extracao e afim em tau:
     delta omega[A + tau B] = delta omega[A] + tau delta omega[B] + sigma * (variaveis de B), sem tau^2 nem
     termos com tau fora dessa forma; e o termo xi (D^sigma + sigma mu) em Omega^+ tem coeficiente xi em d_tau,
     de modo que, apos a divisao por xi, L'(s) = I (L-afim).
 (5) L_{tau X0} gbar = tau L_{X0} gbar + (dtau (x) X0^flat + X0^flat (x) dtau): a parte ++ do termo extra e
     2 g_{+-} d_+ tau, o que fecha a conta do transporte secular.
"""
import sys
import sympy as sp

t, x = sp.symbols('tau xi', real=True)
mu = sp.Symbol('mu', positive=True)
Z = sp.Function('zeta')(t, x)
Qf = sp.Function('Q')(t, x)
F = sp.Function('phi')(t, x)
W = mu + x**2*Qf
e2 = sp.exp(-2*Z)
r2 = sp.exp(-2*Z)*x**2/W**2          # b = coeficiente de g_{S^2}
g2 = sp.Matrix([[-e2*(1 - x**2), -e2*x/mu], [-e2*x/mu, e2/mu**2]])
crd = [t, x]


def lie_metric(X):
    """(L_X g_2D)_{ab} e X(b) para X = (X^tau, X^xi)."""
    L = sp.zeros(2, 2)
    for a in range(2):
        for b in range(2):
            L[a, b] = sum(X[c]*sp.diff(g2[a, b], crd[c]) for c in range(2)) + \
                sum(g2[c, b]*sp.diff(X[c], crd[a]) + g2[a, c]*sp.diff(X[c], crd[b]) for c in range(2))
    return L, sum(X[c]*sp.diff(r2, crd[c]) for c in range(2))


def D(sig, f):
    return sig*sp.diff(f, t) + mu*(1 + sig*x)*sp.diff(f, x)


def omegas(zeta, Qv, phi):
    Wv = mu + x**2*Qv
    L = zeta + sp.log(Wv)
    return {'w1': 2*x*Qv + D(-1, L) + D(1, L),
            'w2-': D(-1, L) + mu, 'w2+': D(1, L) - mu,
            'w3-': D(-1, zeta) + mu, 'w3+': D(1, zeta) - mu,
            'w4-': D(-1, phi), 'w4+': D(1, phi)}


def lin_omegas(dz, dQ, dphi):
    Wv = W
    dL = dz + x**2*dQ/Wv
    return {'w1': 2*x*dQ + D(-1, dL) + D(1, dL),
            'w2-': D(-1, dL), 'w2+': D(1, dL),
            'w3-': D(-1, dz), 'w3+': D(1, dz),
            'w4-': D(-1, dphi), 'w4+': D(1, dphi)}


def extrai(Lg, Lb, dphi):
    """Do ansatz: Lg = lambda * g2 (conferido), lambda = -2 dzeta; Lb/b = -2 dzeta - 2 dW/W."""
    lam = sp.simplify(Lg[1, 1]/g2[1, 1])
    for a in range(2):
        for b in range(2):
            assert sp.simplify(Lg[a, b] - lam*g2[a, b]) == 0, 'fora do ansatz'
    dz = -lam/2
    dW = sp.simplify(-W*(dz + Lb/(2*r2)))
    return dz, sp.simplify(dW/x**2), dphi


def main():
    ok = True
    # (1)
    up, um = -(1 + x)*sp.exp(-mu*t), -(1 - x)*sp.exp(-mu*t)
    Jm = sp.Matrix([[sp.diff(up, t), sp.diff(up, x)], [sp.diff(um, t), sp.diff(um, x)]])
    inv = Jm.inv()
    X0 = [sp.simplify(inv[0, 0] + inv[0, 1]), sp.simplify(inv[1, 0] + inv[1, 1])]
    c1 = [sp.simplify(X0[0] - sp.exp(mu*t)/mu), sp.simplify(X0[1] - sp.exp(mu*t)*x)]
    print('(1) X0 - e^{mu tau}(mu^-1, xi) =', c1)
    ok &= c1 == [0, 0]
    # (2)+(3)
    Lg, Lb = lie_metric(X0)
    dplus = [inv[0, 0], inv[1, 0]]
    dminus = [inv[0, 1], inv[1, 1]]
    gpp = sp.simplify(sum(dplus[a]*dplus[b]*Lg[a, b] for a in range(2) for b in range(2)))
    gmm = sp.simplify(sum(dminus[a]*dminus[b]*Lg[a, b] for a in range(2) for b in range(2)))
    print('(2) (L_X0 g)_{++} =', gpp, ', (L_X0 g)_{--} =', gmm)
    ok &= gpp == 0 and gmm == 0
    dphi = X0[0]*sp.diff(F, t) + X0[1]*sp.diff(F, x)
    dz, dQ, dph = extrai(Lg, Lb, dphi)
    lin = lin_omegas(dz, dQ, dph)
    om = omegas(Z, Qf, F)
    Zop = lambda f: sp.diff(f, t)/mu + x*sp.diff(f, x)
    for k in lin:
        diff = sp.simplify(lin[k] - sp.exp(mu*t)*(Zop(om[k]) + om[k]))
        print(f'(3) {k}: delta omega - e^(mu tau)(Z+1) omega = {diff}')
        ok &= diff == 0
    # (4) estrutura secular: variaveis (dz, dQ, dphi) = A + tau B com A, B genericos
    Az, AQ, Ap, Bz, BQ, Bp = [sp.Function(n)(t, x) for n in ('Az', 'AQ', 'Ap', 'Bz', 'BQ', 'Bp')]
    full = lin_omegas(Az + t*Bz, AQ + t*BQ, Ap + t*Bp)
    la, lb = lin_omegas(Az, AQ, Ap), lin_omegas(Bz, BQ, Bp)
    sig = {'w1': 0, 'w2-': -1, 'w2+': 1, 'w3-': -1, 'w3+': 1, 'w4-': -1, 'w4+': 1}
    extra = {'w1': 0, 'w2-': Bz + x**2*BQ/W, 'w2+': Bz + x**2*BQ/W, 'w3-': Bz, 'w3+': Bz, 'w4-': Bp, 'w4+': Bp}
    for k in full:
        d = sp.simplify(full[k] - la[k] - t*lb[k] - sig[k]*extra[k])
        print(f'(4) {k}: dw[A + tau B] - dw[A] - tau dw[B] - sigma*(B) = {d}')
        ok &= d == 0
    # (5)
    Xs = [t*X0[0], t*X0[1]]
    Ls, _ = lie_metric(Xs)
    spp = sp.simplify(sum(dplus[a]*dplus[b]*Ls[a, b] for a in range(2) for b in range(2)))
    gpm = sp.simplify(sum(dplus[a]*dminus[b]*g2[a, b] for a in range(2) for b in range(2)))
    dpt = sp.simplify(dplus[0])
    c5 = sp.simplify(spp - 2*gpm*dpt)
    print('(5) (L_{tau X0} g)_{++} - 2 g_{+-} d_+ tau =', c5)
    ok &= c5 == 0
    print('R2 SIMBOLICO: CONFERE' if ok else 'R2 SIMBOLICO: FALHA')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
