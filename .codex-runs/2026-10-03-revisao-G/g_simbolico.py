#!/usr/bin/env python3
"""Conferencias simbolicas (sympy) da revisao do Lema G.  Codigo proprio do revisor.

Blocos:
 [A] metrica de RT a partir do frame (ekhfkjfhkfhfkhf) e forma nos nulos u_+-; g_++ = g_-- = 0.
 [B] derivada de Lie em coordenadas nulas: (L_X g)_{++} = 2 g_{-+} d_+ X^-, (L_X g)_{--} = 2 g_{+-} d_- X^+,
     e componentes mistas/angulares para X radial.
 [C] delta tau, delta xi para X = f_-(u_-) d_- + f_+(u_+) d_+ (Passo 3).
 [D] equivariancia: (Theta^2)^* X tem f~(u) = Lambda^-1 f(Lambda u); em (tau, xi), X(tau + 4 pi, xi).
     Perfis f = u^n: expoentes mu(1 - n) e A_n, B_n de GAUGE_DOMAIN_CLASSIFICATION.
 [E] D^+-(u_-+) e a formula X(phi) = (e^{mu tau}/2mu)[f(u_-) omega4^+ - f(u_+) omega4^-] (Passo 5).
 [F] admissibilidade de X_f: delta zeta, delta W, delta Q a partir de L_X g; f = 1 da
     delta zeta = e^{mu tau}(Z zeta - 1), delta Q = e^{mu tau}(Z Q + 2Q); f = u^n: regularidade de delta Q no eixo.
 [G] delta omega_{X0} = e^{mu tau}(Z + 1) omega para as 7 funcoes (ekjehee), zeta, Q, phi genericas.
 [H] dominio M em (tau, xi) e curvas de nivel.
"""
import sympy as sp

ok = []


def confere(nome, expr):
    e = sp.simplify(expr)
    passou = (e == 0)
    ok.append((nome, passou))
    print(("OK   " if passou else "FALHA"), nome, "" if passou else f" resto = {e}")
    return passou


tau, xi, mu = sp.symbols("tau xi mu", real=True)
um, up = sp.symbols("u_m u_p", real=True)  # u_-, u_+
th, ph = sp.symbols("theta phi_ang", real=True)

# relacoes de coordenadas de RT: u_sigma = -(1 + sigma xi) e^{-mu tau}
U_p = -(1 + xi) * sp.exp(-mu * tau)
U_m = -(1 - xi) * sp.exp(-mu * tau)
TAU_u = -sp.log(-(up + um) / 2) / mu
XI_u = (up - um) / (up + um)

# ---------------------------------------------------------------- [A]
print("\n[A] metrica de RT a partir do frame, e forma nos nulos")
zeta = sp.Function("zeta")(tau, xi)
Q = sp.Function("Q")(tau, xi)
W = mu + xi**2 * Q
# Em coordenadas polares (tau, xi, theta, phi): e0 = d_tau + mu xi d_xi;
# sum_i e_i (x) e_i = mu^2 d_xi (x) d_xi + W^2 xi^-2 (g_S2)^-1   (deducao: e_i = mu (x^i/xi) d_xi + W T_i)
# Conferimos esta deducao em cartesianas, num ponto generico, em [A1].
ginv = sp.zeros(4, 4)
e0 = sp.Matrix([1, mu * xi, 0, 0])
ginv += -e0 * e0.T
ginv[1, 1] += mu**2
ginv[2, 2] += W**2 / xi**2
ginv[3, 3] += W**2 / (xi**2 * sp.sin(th) ** 2)
ginv = sp.exp(2 * zeta) * ginv
g = sp.simplify(ginv.inv())
# metrica do enunciado de RT em (u_-, u_+): g = e^{-2 zeta}( -2(du_- du_+ + du_+ du_-)/(mu^2 (u_+ + u_-)^2) + xi^2/W^2 g_S2 )
# Jacobiano (tau, xi) -> (u_-, u_+)
J = sp.Matrix([[sp.diff(U_m, tau), sp.diff(U_m, xi)], [sp.diff(U_p, tau), sp.diff(U_p, xi)]])  # d(u)/d(tau,xi)
g2 = g[:2, :2]
# g em (u): g_u = Jinv^T g2 Jinv
Jinv = sp.simplify(J.inv())
gu = sp.simplify(Jinv.T * g2 * Jinv)
confere("A: g_{--} = 0 (metrica do frame)", gu[0, 0])
confere("A: g_{++} = 0 (metrica do frame)", gu[1, 1])
gmp_RT = -2 * sp.exp(-2 * zeta) / (mu**2 * (U_p + U_m) ** 2)
confere("A: g_{-+} = -2 e^{-2zeta}/(mu^2 (u_++u_-)^2)", gu[0, 1] - gmp_RT)
confere("A: g_{theta theta} = e^{-2zeta} xi^2/W^2", g[2, 2] - sp.exp(-2 * zeta) * xi**2 / W**2)
confere("A: g_{tau theta} = g_{xi theta} = 0", g[0, 2] ** 2 + g[1, 2] ** 2 + g[0, 3] ** 2 + g[1, 3] ** 2)

# [A1] conferencia do frame em cartesianas: sum_i e_i (x) e_i aplicado a dx^a, dx^b
x1, x2, x3 = sp.symbols("x1 x2 x3", real=True)
X = [x1, x2, x3]
r = sp.sqrt(x1**2 + x2**2 + x3**2)
Qs = sp.Symbol("Qs")
# e_i = mu d_i + Q sum_k x^k (x^k d_i - x^i d_k): componentes cartesianas
E = []
for i in range(3):
    comp = [0, 0, 0]
    comp[i] += mu
    for k in range(3):
        comp[i] += Qs * X[k] * X[k]
        comp[k] += -Qs * X[k] * X[i]
    E.append(sp.Matrix(comp))
S = sum((E[i] * E[i].T for i in range(3)), sp.zeros(3, 3))
# esperado: mu^2 xhat xhat^T + (mu + Q r^2)^2 (I - xhat xhat^T)
xh = sp.Matrix(X) / r
esp = mu**2 * xh * xh.T + (mu + Qs * r**2) ** 2 * (sp.eye(3) - xh * xh.T)
confere("A1: sum e_i e_i = mu^2 rr + W^2 (I - rr) em cartesianas", sp.simplify((S - esp).subs({x1: sp.Rational(3, 7), x2: sp.Rational(-2, 5), x3: sp.Rational(1, 3)})).norm())
confere("A1 (simbolico completo)", sp.simplify(sp.expand(S - esp)).norm())
# D^sigma = sigma e0 + sum x^k e_k / xi em funcoes radiais = sigma d_tau + mu(1 + sigma xi) d_xi
for sg in (1, -1):
    radial = sum((X[k] * E[k] for k in range(3)), sp.zeros(3, 1)) / r  # componente cartesiana
    # para f(r): radial . grad f = (radial . xhat) f'
    coef = sp.simplify((radial.T * xh)[0])
    confere(f"A1: D^{sg:+d} radial = mu d_xi (coef)", coef - mu)

# ---------------------------------------------------------------- [B]
print("\n[B] derivada de Lie em coordenadas nulas")
# metrica generica esfericamente simetrica na forma dupla nula: G(u) (du_- du_+ + du_+ du_-) + R2(u) g_S2
Gf = sp.Function("G")(um, up)
R2 = sp.Function("R2")(um, up)
coords = [um, up, th, ph]
gN = sp.zeros(4, 4)
gN[0, 1] = gN[1, 0] = Gf
gN[2, 2] = R2
gN[3, 3] = R2 * sp.sin(th) ** 2
Xm = sp.Function("Xm")(um, up)
Xp = sp.Function("Xp")(um, up)
Xvec = [Xm, Xp, 0, 0]


def lie_g(gm, Xv, cs):
    n = len(cs)
    L = sp.zeros(n, n)
    for a in range(n):
        for b in range(n):
            s = 0
            for c in range(n):
                s += Xv[c] * sp.diff(gm[a, b], cs[c]) + gm[c, b] * sp.diff(Xv[c], cs[a]) + gm[a, c] * sp.diff(Xv[c], cs[b])
            L[a, b] = sp.simplify(s)
    return L


LX = lie_g(gN, Xvec, coords)
confere("B: (L_X g)_{--} = 2 g_{+-} d_- X^+", LX[0, 0] - 2 * Gf * sp.diff(Xp, um))
confere("B: (L_X g)_{++} = 2 g_{-+} d_+ X^-", LX[1, 1] - 2 * Gf * sp.diff(Xm, up))
confere("B: mistas (L_X g)_{+-theta}, _{+-phi} = 0 para X radial", sum(LX[a, b] ** 2 for a in range(2) for b in range(2, 4)))
confere("B: (L_X g)_{theta phi} = 0 e _{phi phi} = sin^2 * _{theta theta}", LX[2, 3] ** 2 + (LX[3, 3] - sp.sin(th) ** 2 * LX[2, 2]) ** 2)
# X com componente angular generica (nao radial): a condicao ++ / -- nao muda
Yth = sp.Function("Yth")(um, up, th, ph)
Yph = sp.Function("Yph")(um, up, th, ph)
LY = lie_g(gN, [Xm, Xp, Yth, Yph], coords)
confere("B: com componentes angulares, (L_X g)_{++} continua 2 g_{-+} d_+ X^-", LY[1, 1] - 2 * Gf * sp.diff(Xm, up))

# ---------------------------------------------------------------- [C]
print("\n[C] delta tau, delta xi (Passo 3)")
fm = sp.Function("f_m")
fp = sp.Function("f_p")
f = sp.Function("f")
# X = f_-(u_-) d_- + f_+(u_+) d_+ : X(tau) = f_- d tau/du_- + f_+ d tau/du_+
Xtau = fm(um) * sp.diff(TAU_u, um) + fp(up) * sp.diff(TAU_u, up)
Xxi = fm(um) * sp.diff(XI_u, um) + fp(up) * sp.diff(XI_u, up)
sub = {um: U_m, up: U_p}
Xtau_t = sp.simplify(Xtau.subs(sub))
Xxi_t = sp.simplify(Xxi.subs(sub))
confere("C: delta tau = e^{mu tau}[f(u_+) + f(u_-)]/(2 mu)", Xtau_t - sp.exp(mu * tau) * (fp(U_p) + fm(U_m)) / (2 * mu))
confere("C: delta xi = e^{mu tau}[(1+xi) f(u_-) - (1-xi) f(u_+)]/2", Xxi_t - sp.exp(mu * tau) * ((1 + xi) * fm(U_m) - (1 - xi) * fp(U_p)) / 2)
# eixo: X(u_+ - u_-) em xi = 0 e tangencia
axis = sp.simplify((Xxi_t).subs(xi, 0))
confere("C: delta xi em xi=0 = e^{mu tau}[f_-(u) - f_+(u)]/2, u = -e^{-mu tau}",
        axis - sp.exp(mu * tau) * (fm(-sp.exp(-mu * tau)) - fp(-sp.exp(-mu * tau))) / 2)
# f = 1: X0 = e^{mu tau}(mu^-1 d_tau + xi d_xi)
um1 = sp.Lambda(sp.Symbol("w"), 1)
X0tau = sp.simplify(Xtau_t.replace(fm, um1).replace(fp, um1).doit())
X0xi = sp.simplify(Xxi_t.replace(fm, um1).replace(fp, um1).doit())
confere("C: f=1 -> delta tau = e^{mu tau}/mu", X0tau - sp.exp(mu * tau) / mu)
confere("C: f=1 -> delta xi = e^{mu tau} xi", X0xi - sp.exp(mu * tau) * xi)

# ---------------------------------------------------------------- [D]
print("\n[D] equivariancia e perfis u^n")
Lam = sp.exp(-4 * sp.pi * mu)
# u_sigma(tau + 4pi, xi) = Lambda u_sigma(tau, xi)
confere("D: u_+(tau+4pi) = Lambda u_+", U_p.subs(tau, tau + 4 * sp.pi) - Lam * U_p)
confere("D: u_-(tau+4pi) = Lambda u_-", U_m.subs(tau, tau + 4 * sp.pi) - Lam * U_m)
# pullback por F = Theta^2: (F^*X)(p) = dF_p^{-1} X(F(p)).  Em (tau, xi), F = translacao, dF = id:
# (F^*X)(tau, xi) = X(tau + 4pi, xi).  Com X = X_f: componente tau
dtau_f = lambda ff: sp.exp(mu * tau) * (ff(U_p) + ff(U_m)) / (2 * mu)
dxi_f = lambda ff: sp.exp(mu * tau) * ((1 + xi) * ff(U_m) - (1 - xi) * ff(U_p)) / 2
ftil = sp.Lambda(sp.Symbol("w"), f(Lam * sp.Symbol("w")) / Lam)
lado_esq_tau = dtau_f(f).subs(tau, tau + 4 * sp.pi)
lado_dir_tau = dtau_f(ftil)
confere("D: (Theta^2)^*X_f = X_{f~}, f~(u) = Lambda^-1 f(Lambda u) [comp. tau]", sp.simplify(lado_esq_tau - lado_dir_tau))
confere("D: idem [comp. xi]", sp.simplify(dxi_f(f).subs(tau, tau + 4 * sp.pi) - dxi_f(ftil)))
# em coordenadas u: F(u) = Lambda u, dF = Lambda, (F^*X)^u(u) = Lambda^-1 X^u(Lambda u) -> mesmo f~ (conferido acima)
for n in range(0, 6):
    fn = sp.Lambda(sp.Symbol("w"), sp.Symbol("w") ** n)
    sn = mu * (1 - n)
    An = (-1) ** n * ((1 + xi) ** n + (1 - xi) ** n) / (2 * mu)
    Bn = (-1) ** n * ((1 + xi) * (1 - xi) ** n - (1 - xi) * (1 + xi) ** n) / 2
    confere(f"D: f=u^{n}: delta tau = e^(mu(1-n)tau) A_n", sp.simplify(sp.expand(dtau_f(fn) - sp.exp(sn * tau) * An)))
    confere(f"D: f=u^{n}: delta xi = e^(mu(1-n)tau) B_n", sp.simplify(sp.expand(dxi_f(fn) - sp.exp(sn * tau) * Bn)))
    # multiplicador: (Theta^2)^* X = e^{4 pi s_n} X  <=> f~ = e^{4pi s_n} f
    w = sp.Symbol("w")
    confere(f"D: f=u^{n}: f~ = e^(4 pi mu (1-n)) f", sp.simplify(ftil(w).subs(f(Lam * w), (Lam * w) ** n) - sp.exp(4 * sp.pi * sn) * w**n))

# ---------------------------------------------------------------- [E]
print("\n[E] Passo 5: D^sigma(u) e X(phi)")
Dp = lambda F: sp.diff(F, tau) + mu * (1 + xi) * sp.diff(F, xi)
Dm = lambda F: -sp.diff(F, tau) + mu * (1 - xi) * sp.diff(F, xi)
confere("E: D^+(u_+) = 0", Dp(U_p))
confere("E: D^-(u_-) = 0", Dm(U_m))
confere("E: D^+(u_-) = 2 mu e^{-mu tau}", Dp(U_m) - 2 * mu * sp.exp(-mu * tau))
confere("E: D^-(u_+) = -2 mu e^{-mu tau}", Dm(U_p) + 2 * mu * sp.exp(-mu * tau))
phi = sp.Function("varphi")(tau, xi)
Xphi = dtau_f(f) * sp.diff(phi, tau) + dxi_f(f) * sp.diff(phi, xi)
form = sp.exp(mu * tau) / (2 * mu) * (f(U_m) * Dp(phi) - f(U_p) * Dm(phi))
confere("E: X(phi) = (e^{mu tau}/2mu)[f(u_-) omega4^+ - f(u_+) omega4^-]", sp.simplify(sp.expand(Xphi - form)))
# comutador [D^sigma, X0] = e^{mu tau} D^sigma
X0 = lambda F: sp.exp(mu * tau) * (sp.diff(F, tau) / mu + xi * sp.diff(F, xi))
for nome, D in (("+", Dp), ("-", Dm)):
    confere(f"E: [D^{nome}, X0] = e^(mu tau) D^{nome}", sp.simplify(sp.expand(D(X0(phi)) - X0(D(phi)) - sp.exp(mu * tau) * D(phi))))

# ---------------------------------------------------------------- [F]
print("\n[F] admissibilidade de X_f: delta zeta, delta W, delta Q")
# metrica RT em (u_-, u_+): g_{-+} = -2 e^{-2 zeta}/(mu^2 (u_+ + u_-)^2), r^2 = e^{-2zeta} xi^2 / W^2
zu = sp.Function("zeta_u")(um, up)
Wu = sp.Function("W_u")(um, up)
gmp_u = -2 * sp.exp(-2 * zu) / (mu**2 * (up + um) ** 2)
R2u = sp.exp(-2 * zu) * XI_u**2 / Wu**2
gRT = sp.zeros(4, 4)
gRT[0, 1] = gRT[1, 0] = gmp_u
gRT[2, 2] = R2u
gRT[3, 3] = R2u * sp.sin(th) ** 2
Xf = [f(um), f(up), 0, 0]
LXf = lie_g(gRT, Xf, coords)
confere("F: (L_X g)_{++} = (L_X g)_{--} = 0 para X_f", LXf[0, 0] ** 2 + LXf[1, 1] ** 2)
# variacao do ansatz: delta g_{-+} = -2 delta zeta g_{-+};  delta r^2 = r^2(-2 delta zeta - 2 delta W/W)
dzeta = sp.simplify(-LXf[0, 1] / (2 * gmp_u))
dW = sp.simplify(Wu * (-LXf[2, 2] / R2u - 2 * dzeta) / 2)
# formula fechada: delta zeta = X(zeta) + (f(u_-)+f(u_+))/(u_++u_-) - (f'(u_-)+f'(u_+))/2
Xz = f(um) * sp.diff(zu, um) + f(up) * sp.diff(zu, up)
dzeta_esp = Xz + (f(um) + f(up)) / (up + um) - (sp.diff(f(um), um) + sp.diff(f(up), up)) / 2
confere("F: delta zeta = X(zeta) + (f_-+f_+)/(u_++u_-) - (f_-'+f_+')/2", dzeta - dzeta_esp)
XW = f(um) * sp.diff(Wu, um) + f(up) * sp.diff(Wu, up)
# delta W - X(W) = W * [ X(xi)/xi - X(zeta) - delta zeta ... ] -> calculamos o resto 'extra'
extra = sp.simplify(dW - XW)
print("     delta W - X(W) =", sp.simplify(extra / Wu), "* W")
# f = 1: extra deve ser 0 (delta W = X0(W)), e delta zeta = X0(zeta) - e^{mu tau}
extra1 = sp.simplify((extra / Wu).subs(f(um), 1).subs(f(up), 1).doit())
confere("F: f=1 -> delta W = X0(W) (logo delta Q = e^{mu tau}(ZQ + 2Q))", extra1)
dz1 = sp.simplify((dzeta - Xz).subs({f(um): 1, f(up): 1}).doit().subs(sub))
confere("F: f=1 -> delta zeta = X0(zeta) - e^{mu tau}", dz1 + sp.exp(mu * tau))
# delta Q a partir de delta W = X0(W): W = mu + xi^2 Q, X0 = e^{mu tau} Z, Z xi = xi
Zop = lambda F: sp.diff(F, tau) / mu + xi * sp.diff(F, xi)
confere("F: X0(mu + xi^2 Q) = xi^2 e^{mu tau}(ZQ + 2Q)", sp.simplify(sp.exp(mu * tau) * Zop(mu + xi**2 * Q) - xi**2 * sp.exp(mu * tau) * (Zop(Q) + 2 * Q)))
# f = u^n: o termo extra em (tau, xi) e sua paridade/regularidade no eixo
for n in range(0, 5):
    w = sp.Symbol("w")
    e_n = sp.simplify((extra / Wu).subs({f(um): um**n, f(up): up**n}).doit())
    e_n = sp.simplify(e_n.subs(sub))
    # delta Q = (X(W) + W*e_n)/xi^2: X(W) e regular (W regular e par); o termo W e_n / xi^2 precisa ser regular
    ser = sp.series(sp.simplify(e_n / xi**2), xi, 0, 3).removeO()
    print(f"     f=u^{n}: (delta W - X W)/W = {sp.factor(e_n)};  / xi^2 ~ {sp.simplify(ser)}")

# ---------------------------------------------------------------- [G]
print("\n[G] delta omega_{X0} = e^{mu tau}(Z + 1) omega  (ekjehee)")
F_ = zeta + sp.log(W)
D = {1: Dp, -1: Dm}
om = {}
om["1"] = 2 * xi * Q + Dm(F_) + Dp(F_)
for sg in (1, -1):
    om[f"2{sg:+d}"] = D[sg](F_) - sg * mu
    om[f"3{sg:+d}"] = D[sg](zeta) - sg * mu
    om[f"4{sg:+d}"] = D[sg](phi)
eps = sp.Symbol("eps")
dz = sp.exp(mu * tau) * (Zop(zeta) - 1)
dQ = sp.exp(mu * tau) * (Zop(Q) + 2 * Q)
dphi = sp.exp(mu * tau) * Zop(phi)
rep = {zeta: zeta + eps * dz, Q: Q + eps * dQ, phi: phi + eps * dphi}


def var(expr):
    e2 = expr.subs(rep, simultaneous=True)
    return sp.diff(e2, eps).subs(eps, 0)


for k, v in om.items():
    alvo = sp.exp(mu * tau) * (Zop(v) + v)
    confere(f"G: delta omega_{k} = e^(mu tau)(Z+1) omega_{k}", sp.simplify(sp.expand(var(v) - alvo)))
# conformal: zeta -> zeta + c nao muda omega
c = sp.Symbol("c")
for k, v in om.items():
    confere(f"G: omega_{k} invariante por zeta -> zeta + c", sp.simplify(v.subs(zeta, zeta + c).doit() - v))

# ---------------------------------------------------------------- [H]
print("\n[H] dominio M")
# M = {u_+ < 0, u_+ < u_- < -u_+/81} em (tau, xi)
e = sp.exp(-mu * tau)
print("     u_+ < 0  <=>", sp.simplify(sp.solve_univariate_inequality(-(1 + xi) < 0, xi)))
print("     u_+ < u_- <=>", sp.solve_univariate_inequality(-(1 + xi) < -(1 - xi), xi))
print("     u_- < -u_+/81 <=>", sp.solve_univariate_inequality(-(1 - xi) < (1 + xi) / 81, xi))

print("\nRESUMO:", sum(1 for _, p in ok if p), "de", len(ok), "conferencias OK")
falhas = [n for n, p in ok if not p]
if falhas:
    print("FALHAS:", falhas)
