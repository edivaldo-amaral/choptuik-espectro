#!/usr/bin/env python3
"""G5, versao linear (sem fluxos): confere que as variacoes do ansatz produzidas por X_f satisfazem
   [delta_{X_f}](tau + 4 pi, xi) = [delta_{X_{f~}}](tau, xi),   f~(u) = Lambda^-1 f(Lambda u),
para delta zeta, delta W e delta phi, com fundo generico satisfazendo zeta(tau + 4pi) = zeta + 2K,
W e phi 4pi-periodicos.  Como as omega sao obtidas de (zeta, Q, phi) por operadores D^sigma de
coeficientes independentes de tau, e invariantes por zeta -> zeta + c, segue
   delta omega_{X_f} o Theta^2 = delta omega_{(Theta^2)^* X_f}.
As formulas de delta zeta e delta W vem do bloco [F] de g_simbolico.py (derivadas de L_X g)."""
import sympy as sp

tau, xi, mu, K = sp.symbols("tau xi mu K", real=True)
um, up, w = sp.symbols("u_m u_p w", real=True)
f = sp.Function("f")
Lam = sp.exp(-4 * sp.pi * mu)
ftil = sp.Lambda(w, f(Lam * w) / Lam)
U_p = -(1 + xi) * sp.exp(-mu * tau)
U_m = -(1 - xi) * sp.exp(-mu * tau)
# fundo: zeta = zeta0(tau, xi) com zeta0(tau+4pi) = zeta0 + 2K -> modelamos zeta = (K/(2pi)) tau + Zp(tau, xi),
# Zp 4pi-periodica; W, phi periodicas: usamos funcoes de (cos(tau/2), sin(tau/2), xi) genericas.
c, s = sp.cos(tau / 2), sp.sin(tau / 2)
Zp = sp.Function("Zp")(c, s, xi)
Wp = sp.Function("Wp")(c, s, xi)
Pp = sp.Function("Pp")(c, s, xi)
zeta = K * tau / (2 * sp.pi) + Zp


def campo(ff):
    dt = sp.exp(mu * tau) * (ff(U_p) + ff(U_m)) / (2 * mu)
    dx = sp.exp(mu * tau) * ((1 + xi) * ff(U_m) - (1 - xi) * ff(U_p)) / 2
    return dt, dx


def aplica(ff, F):
    dt, dx = campo(ff)
    return dt * sp.diff(F, tau) + dx * sp.diff(F, xi)


def dlinha(ff, u):
    return sp.diff(ff(w), w).subs(w, u)


def variacoes(ff):
    # delta zeta = X(zeta) + (f(u_-)+f(u_+))/(u_++u_-) - (f'(u_-)+f'(u_+))/2
    dz = aplica(ff, zeta) + (ff(U_m) + ff(U_p)) / (U_p + U_m) - (dlinha(ff, U_m) + dlinha(ff, U_p)) / 2
    # delta W = X(W) + W * extra, extra = [(u_- - u_+)(f'(u_-)+f'(u_+)) - 2(f(u_-)-f(u_+))]/(2(u_- - u_+))
    extra = ((U_m - U_p) * (dlinha(ff, U_m) + dlinha(ff, U_p)) - 2 * (ff(U_m) - ff(U_p))) / (2 * (U_m - U_p))
    dW = aplica(ff, Wp) + Wp * extra
    dphi = aplica(ff, Pp)
    return dz, dW, dphi


vf = variacoes(f)
vt = variacoes(ftil)
nomes = ("delta zeta", "delta W", "delta phi")
tudo_ok = True
for nome, a, b in zip(nomes, vf, vt):
    dif = sp.simplify(sp.expand(a.subs(tau, tau + 4 * sp.pi) - b))
    ok = dif == 0
    tudo_ok &= ok
    print(("OK   " if ok else "FALHA"), f"{nome}: [X_f](tau+4pi) = [X_f~](tau)", "" if ok else dif)
print("TUDO OK" if tudo_ok else "HA FALHAS")
