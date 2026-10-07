"""C1 (revisao CRH, 03/10/2026): conferencia simbolica do Lema C (GAUGE_GLOBAL.md, par. 6).

Rodar com sympy no PYTHONPATH (instalado no scratchpad, nao no .venv):
    PYTHONPATH=<scratch>/pylib OPENBLAS_NUM_THREADS=2 .venv/bin/python c1_lema_c.py

Funcoes GENERICAS zeta(tau,xi), Q(tau,xi), phi(tau,xi); mu, eps simbolicos.
Confere:
  (a) u_sigma = -(1+sigma xi) e^{-mu tau}; psi_eps em (tau,xi) e' (T, Xi) com
      u_sigma(T,Xi) = u_sigma(tau,xi) + eps (preserva du_+-, eixo xi=0).
  (b) a metrica radial de RT, e^{-2 zeta}(-dtau^2 + (dxi - mu xi dtau)^2/mu^2), e'
      -2 e^{-2 zeta}(du_- du_+ + du_+ du_-)/(mu^2 (u_+ + u_-)^2)  (fator (u_+ + u_-)^2).
  (c) psi_eps^* g esta na forma de RT com zeta_eps, Q_eps (formulas do Lema C) e W_eps = W o psi_eps.
  (d) d/deps omega_eps |_{eps=0} = e^{mu tau}(Z+1) omega, nas SETE componentes
      omega1, omega2^+-, omega3^+-, omega4^+- (definicoes ekjehee de RT), Z = mu^{-1}d_tau + xi d_xi.
  (e) X0 = d_{u-} + d_{u+} = e^{mu tau}(mu^{-1} d_tau + xi d_xi).
"""
import sympy as sp

tau, xi = sp.symbols('tau xi', real=True)
mu = sp.symbols('mu', positive=True)
eps = sp.symbols('epsilon', real=True)
zeta = sp.Function('zeta')
Qf = sp.Function('Q')
phi = sp.Function('phi')

def u(sig, t, x):
    return -(1 + sig*x)*sp.exp(-mu*t)

# psi_eps em (tau, xi)
T = -sp.log(sp.exp(-mu*tau) - eps)/mu
Xi = xi*sp.exp(-mu*tau)/(sp.exp(-mu*tau) - eps)

ok = True
def check(nome, expr):
    global ok
    e = sp.simplify(expr)
    good = (e == 0)
    ok = ok and good
    print(f'[{"OK" if good else "FALHA"}] {nome}' + ('' if good else f'  resto = {e}'))

# (a)
for sig in (1, -1):
    check(f'(a) u_{"+" if sig>0 else "-"}(psi_eps) = u + eps', u(sig, T, Xi) - (u(sig, tau, xi) + eps))
check('(a) eixo: Xi(xi=0) = 0', Xi.subs(xi, 0))

# (b) metrica radial em (tau,xi): matriz G
def G_rt(z, q, t, x):
    e2 = sp.exp(-2*z)
    return sp.Matrix([[e2*(x**2 - 1), -e2*x/mu], [-e2*x/mu, e2/mu**2]])

def R_rt(z, q, t, x):
    return sp.exp(-z)*x/(mu + q*x**2)

z0, q0 = sp.symbols('z0 q0', real=True)
up, um = u(1, tau, xi), u(-1, tau, xi)
dup = sp.Matrix([sp.diff(up, tau), sp.diff(up, xi)])
dum = sp.Matrix([sp.diff(um, tau), sp.diff(um, xi)])
G_null = -2*sp.exp(-2*z0)*(dum*dup.T + dup*dum.T)/(mu**2*(up + um)**2)
check('(b) forma nula -2e^{-2z}(du_-du_+ + du_+du_-)/(mu^2(u_++u_-)^2) = forma RT',
      (G_null - G_rt(z0, q0, tau, xi)).applyfunc(sp.simplify).norm())

# (c) pullback por psi_eps
J = sp.Matrix([[sp.diff(T, tau), sp.diff(T, xi)], [sp.diff(Xi, tau), sp.diff(Xi, xi)]])
zT, qT = zeta(T, Xi), Qf(T, Xi)
G_pull = J.T*G_rt(zT, qT, T, Xi)*J
R_pull = R_rt(zT, qT, T, Xi)
S = up + um
# formulas do Lema C: e^{-2 zeta_eps} = e^{-2 zeta o psi}(S)^2/(S+2eps)^2 ; R o psi = e^{-zeta_eps} xi / W_eps
zeta_eps = zT - sp.log(S/(S + 2*eps))          # S<0 e S+2eps<0 para eps pequeno: razao positiva
W_eps = sp.exp(-zeta_eps)*xi/R_pull
Q_eps = sp.simplify((W_eps - mu)/xi**2)
check('(c) e^{-2zeta_eps} = e^{-2 zeta o psi} S^2/(S+2eps)^2',
      sp.exp(-2*zeta_eps) - sp.exp(-2*zT)*S**2/(S + 2*eps)**2)
check('(c) parte radial de psi^*g = forma RT com zeta_eps',
      (G_pull - G_rt(zeta_eps, Q_eps, tau, xi)).applyfunc(sp.simplify).norm())
check('(c) W_eps = W o psi_eps (o fator (u_++u_-)^2 cancela)',
      W_eps - (mu + qT*Xi**2))
check('(c) Q_eps = Q o psi * (Xi/xi)^2, com Xi/xi independente de xi',
      Q_eps - qT*(Xi/xi)**2)
check('(c) Xi/xi nao depende de xi', sp.diff(sp.simplify(Xi/xi), xi))

# (d) omega a partir de (zeta, Q, phi)
def Dsig(sig, f):
    return sig*sp.diff(f, tau) + mu*(1 + sig*xi)*sp.diff(f, xi)

def omegas(z, q, p):
    F = z + sp.log(mu + xi**2*q)
    return {
        'omega1': 2*xi*q + Dsig(-1, F) + Dsig(1, F),
        'omega2+': Dsig(1, F) - mu, 'omega2-': Dsig(-1, F) + mu,
        'omega3+': Dsig(1, z) - mu, 'omega3-': Dsig(-1, z) + mu,
        'omega4+': Dsig(1, p), 'omega4-': Dsig(-1, p),
    }

# Com funcoes genericas o sympy deixa objetos Subs com eps dentro (artefato, 1a tentativa).
# Como a identidade em (d) depende so do 2-jato de (zeta, Q, phi) no ponto, basta conferi-la
# para polinomios de grau 3 com coeficientes simbolicos (o 2-jato deles num ponto e' arbitrario).
def poly(nome):
    cs = sp.symbols(f'{nome}0:10')
    mons = [tau**i*xi**j for i in range(4) for j in range(4 - i)]
    return sp.Lambda((sp.Symbol('t_'), sp.Symbol('x_')),
                     sum(c*m for c, m in zip(cs, mons)).subs({tau: sp.Symbol('t_'), xi: sp.Symbol('x_')},
                                                             simultaneous=True))
zP, qP, pP = poly('a'), poly('b'), poly('c')
zeta_epsP = zP(T, Xi) - sp.log(S/(S + 2*eps))
Q_epsP = qP(T, Xi)*(Xi/xi)**2
w0 = omegas(zP(tau, xi), qP(tau, xi), pP(tau, xi))
we = omegas(zeta_epsP, Q_epsP, pP(T, Xi))
Zop = lambda f: sp.diff(f, tau)/mu + xi*sp.diff(f, xi)
for k in w0:
    dk = sp.diff(we[k], eps).subs(eps, 0).doit()
    alvo = sp.exp(mu*tau)*(Zop(w0[k]) + w0[k])
    check(f'(d) d/deps {k}_eps = e^(mu tau)(Z+1){k}', dk - alvo)

# variacoes de (zeta, Q, phi) separadas, como em GAUGE_RT_ACTION.md
dz = sp.diff(zeta_eps, eps).subs(eps, 0).doit()
dq = sp.diff(Q_eps, eps).subs(eps, 0).doit()
check('(d) delta zeta = e^{mu tau}(Z zeta - 1)', dz - sp.exp(mu*tau)*(Zop(zeta(tau, xi)) - 1))
check('(d) delta Q = e^{mu tau}(Z Q + 2Q)', dq - sp.exp(mu*tau)*(Zop(Qf(tau, xi)) + 2*Qf(tau, xi)))

# (e) gerador
check('(e) X0^tau = d/deps T = e^{mu tau}/mu', sp.diff(T, eps).subs(eps, 0) - sp.exp(mu*tau)/mu)
check('(e) X0^xi = d/deps Xi = xi e^{mu tau}', sp.diff(Xi, eps).subs(eps, 0) - xi*sp.exp(mu*tau))

print('TUDO OK' if ok else 'HA FALHAS')
