#!/usr/bin/env python3
"""Revisao V (03/10/2026): identidades de V1, V2 e V3 com funcoes GENERICAS (sympy).

Os Omega sao transcritos do TeX de RT (eezuuirzr, versao fora do bloco comment), NAO das
tabelas C nem de sharp_propagation.py. Tudo e off-shell: omega de fundo arbitrario, sem supor
Omega(omega) = 0, e perturbacao arbitraria e^{s tau} h.

Uso: PYTHONPATH=<sympy> python v123_identidades.py
"""
import sympy as sp

tau, xi, mu, s, eps = sp.symbols('tau xi mu s epsilon')
F = lambda n: sp.Function(n)(tau, xi)


def D(sig, f, shift=0):
    """D^sigma f = sigma d_tau f + mu(1+sigma xi) d_xi f; com shift: D_s^sigma (age em e^{s tau} f)."""
    return sig*(sp.diff(f, tau) + shift*f) + mu*(1 + sig*xi)*sp.diff(f, xi)


def Omegas(w, sig, Dop):
    """Os cinco Omega^sigma de RT (eezuuirzr). w: dict com 'w1', ('w2',+1), ('w2',-1), ...
    Dop(sig, f): derivada nula usada (D ou D_s)."""
    w1 = w['w1']
    w2 = {1: w[('w2', 1)], -1: w[('w2', -1)]}
    w3 = {1: w[('w3', 1)], -1: w[('w3', -1)]}
    w4 = {1: w[('w4', 1)], -1: w[('w4', -1)]}
    S, M = sig, -sig
    q = sp.Rational(1, 4)*xi
    O1 = (xi*(Dop(S, w1) + S*mu*w1) + mu*2*w1
          + q*(w1*w1 - 2*w2[M]*w1 - 6*w2[S]*w1 + 4*w3[S]*w1 + w2[M]*w2[M]
               - 2*w2[S]*w2[M] - 4*w3[S]*w2[M] + w2[S]*w2[S] + 4*w3[S]*w2[S] - 4*w4[S]*w4[S]))
    O2 = (xi*(Dop(S, w2[M]) + S*mu*w2[M]) + mu*(w1 + w2[M] + w2[S])
          + q*(w1*w1 - 2*w2[-1]*w1 - 2*w2[1]*w1 + w2[-1]*w2[-1] - 6*w2[1]*w2[-1] + w2[1]*w2[1]))
    O2s = (xi*(Dop(S, w2[S]) + S*mu*w2[S]) + mu*(2*w2[S] - 2*w3[S])
           + q*(-4*w2[S]*w2[S] + 8*w3[S]*w2[S] - 4*w4[S]*w4[S]))
    O3 = (xi*(Dop(S, w3[M]) + S*mu*w3[M]) + mu*(-w1)
          + q*(-w1*w1 + 2*w2[-1]*w1 + 2*w2[1]*w1 - w2[-1]*w2[-1] + 2*w2[1]*w2[-1]
               - w2[1]*w2[1] - 4*w4[1]*w4[-1]))
    O4 = (xi*(Dop(S, w4[M]) + S*mu*w4[M]) + mu*(w4[M] + w4[S])
          + q*(-4*w4[1]*w2[-1] - 4*w4[-1]*w2[1]))
    return {'1': O1, '2': O2, '2s': O2s, '3': O3, '4': O4}


def generic(prefix):
    w = {'w1': F(prefix + '1')}
    for i in (2, 3, 4):
        w[(f'w{i}', 1)] = F(f'{prefix}{i}p')
        w[(f'w{i}', -1)] = F(f'{prefix}{i}m')
    return w


def linearized(w, h, sig):
    """delta Omega^sigma = e^{-s tau} d/d eps Omega^sigma(w + eps e^{s tau} h)|_{eps=0}."""
    E = sp.exp(s*tau)
    wp = {k: w[k] + eps*E*h[k] for k in w}
    Om = Omegas(wp, sig, lambda S, f: D(S, f))
    return {k: sp.simplify(sp.expand(sp.diff(v, eps).subs(eps, 0)/E)) for k, v in Om.items()}


def zero(expr):
    return sp.simplify(sp.expand(expr)) == 0


resultados = {}
w = generic('w')
h = generic('h')

# (I1) identidade algebrica de log W, off-shell, sigma = +-, omega^+ e omega^- independentes.
W = mu + xi*(w['w1'] - w[('w2', -1)] - w[('w2', 1)])/2
for sig in (1, -1):
    Om = Omegas(w, sig, lambda S, f: D(S, f))
    lhs = D(sig, W) - W*(w[('w2', sig)] - w[('w3', sig)])
    rhs = (Om['1'] - Om['2'] - Om['2s'])/2
    resultados[f'I1 logW sigma={sig:+d}'] = zero(lhs - rhs)

# (I2) linearizacao (sem dividir por W) e versao dividida com o termo de fundo aproximado.
dW = xi*(h['w1'] - h[('w2', -1)] - h[('w2', 1)])/2          # = xi^2 q
for sig in (1, -1):
    dO = linearized(w, h, sig)
    Om = Omegas(w, sig, lambda S, f: D(S, f))
    lin = D(sig, dW, s) - dW*(w[('w2', sig)] - w[('w3', sig)]) - W*(h[('w2', sig)] - h[('w3', sig)])
    resultados[f'I2 logW linearizada sigma={sig:+d}'] = zero(lin - (dO['1'] - dO['2'] - dO['2s'])/2)
    resid = D(sig, W) - W*(w[('w2', sig)] - w[('w3', sig)])
    div = (D(sig, dW/W, s) - (h[('w2', sig)] - h[('w3', sig)])
           - (dO['1'] - dO['2'] - dO['2s'])/(2*W) + resid*dW/W**2)
    resultados[f'I2b dividida por W (termo -resid*dW/W^2) sigma={sig:+d}'] = zero(div)

# (I3) compatibilidade -2 mu xi (T_s b_i - d_xi a_i) = dOmega_j^- - dOmega_j^+, h^+- independentes.
dOm = linearized(w, h, -1)
dOp = linearized(w, h, 1)
for i, j in ((3, '3'), (4, '4'), (2, '2')):
    hp, hm = h[(f'w{i}', 1)], h[(f'w{i}', -1)]
    a = ((1 - xi)*hp - (1 + xi)*hm)/2
    b = (hp + hm)/(2*mu)
    comp = -2*mu*xi*(sp.diff(b, tau) + s*b - sp.diff(a, xi))
    resultados[f'I3 compat (i,j)=({i},{j})'] = zero(comp - (dOm[j] - dOp[j]))
    # controle: com o par trocado (i, j) a identidade deve FALHAR
outros = {'3': 4, '4': 2, '2': 3}
for j, i in outros.items():
    hp, hm = h[(f'w{i}', 1)], h[(f'w{i}', -1)]
    a = ((1 - xi)*hp - (1 + xi)*hm)/2
    b = (hp + hm)/(2*mu)
    comp = -2*mu*xi*(sp.diff(b, tau) + s*b - sp.diff(a, xi))
    resultados[f'I3 controle par errado (i,j)=({i},{j}) deve ser False'] = zero(comp - (dOm[j] - dOp[j]))
# controle: diferencas de Omega1 e Omega2sharp NAO sao exatas
resultados['I3 controle dOmega1^- - dOmega1^+ independente de a,b (deve ser False: nao e zero)'] = zero(dOm['1'] - dOp['1'])

# (I4) D_s^sigma z = h^sigma  <=>  T_s z = a, d_xi z = b  (h^sigma := D_s^sigma z).
z = F('z')
hp, hm = D(1, z, s), D(-1, z, s)
a = ((1 - xi)*hp - (1 + xi)*hm)/2
b = (hp + hm)/(2*mu)
resultados['I4 a(D_s z) = T_s z'] = zero(a - (sp.diff(z, tau) + s*z))
resultados['I4 b(D_s z) = d_xi z'] = zero(b - sp.diff(z, xi))
# reciproca: se T_s z = a e z_xi = b, entao D_s^sigma z = sigma a + mu(1+sigma xi) b = h^sigma
A, B, Hp, Hm = sp.symbols('A B Hp Hm')
aa = ((1 - xi)*Hp - (1 + xi)*Hm)/2
bb = (Hp + Hm)/(2*mu)
resultados['I4 reciproca sigma=+'] = zero((aa + mu*(1 + xi)*bb) - Hp)
resultados['I4 reciproca sigma=-'] = zero((-aa + mu*(1 - xi)*bb) - Hm)

# (I5) codificacao de RT: omega1 impar, omega_i^+ = -P omega_i; entao P Omega_i^s = -Omega_i^{-s}
#      (nao linear e linearizado). P implementado por xi -> -xi em funcoes genericas.
def P(expr):
    return expr.subs(xi, -xi).doit()


def encoded(prefix):
    U = sp.Function(prefix + '1')
    w1 = U(tau, xi) - U(tau, -xi)                     # impar generica
    e = {'w1': w1}
    for i in (2, 3, 4):
        f = F(f'{prefix}{i}')
        e[(f'w{i}', -1)] = f
        e[(f'w{i}', 1)] = -P(f)
    return e


we, he = encoded('W'), encoded('H')
for sig in (1, -1):
    Om_s = Omegas(we, sig, lambda S, f: D(S, f))
    Om_m = Omegas(we, -sig, lambda S, f: D(S, f))
    for k in Om_s:
        resultados[f'I5 P Omega_{k}^({sig:+d}) = -Omega_{k}^({-sig:+d})'] = zero(P(Om_s[k]) + Om_m[k])
    dS = linearized(we, he, sig)
    dM = linearized(we, he, -sig)
    for k in dS:
        resultados[f'I5 lin P dOmega_{k}^({sig:+d}) = -dOmega_{k}^({-sig:+d})'] = zero(P(dS[k]) + dM[k])

# (I6) V1: as componentes selecionadas de dOmega^+ (1 par, 2, 3, 4) se anulam em xi = 0;
#      as de constraint (1 impar, 2sharp) nao (controle).
dP = linearized(we, he, 1)
sel = {'1 parte par': (dP['1'] + P(dP['1']))/2, '2': dP['2'], '3': dP['3'], '4': dP['4']}
for k, v in sel.items():
    resultados[f'I6 dOmega^+_{k} em xi=0 e zero'] = zero(v.subs(xi, 0).doit())
nsel = {'1 parte impar (derivada em xi=0)': sp.diff((dP['1'] - P(dP['1']))/2, xi).subs(xi, 0).doit(),
        '2sharp': dP['2s'].subs(xi, 0).doit()}
for k, v in nsel.items():
    resultados[f'I6 controle dOmega^+_{k} em xi=0 nao e identicamente zero (deve ser False)'] = zero(v)

# (I7) o operador C de signed_constraints: parte linear livre de S#dOmega^+ (sem Gamma2):
#      (1/2)(1-P) dOmega_1^+ e -P dOmega_2sharp^+ = dOmega_2sharp^-, com w = 0.
w0 = {k: sp.Integer(0) for k in we}
dP0 = linearized(w0, he, 1)
c1 = sp.expand(((dP0['1'] - P(dP0['1']))/2))
h1 = he['w1']
resultados['I7 C linha 1 livre = mu(xi d_xi + 2) h1'] = zero(c1 - mu*(xi*sp.diff(h1, xi) + 2*h1))
c2 = sp.expand(-P(dP0['2s']))
h2, h3 = he[('w2', -1)], he[('w3', -1)]
alvo = (-xi*(sp.diff(h2, tau) + s*h2 + mu*h2) + mu*(xi - xi**2)*sp.diff(h2, xi) + 2*mu*(h2 - h3))
resultados['I7 C linha 2 livre = -xi(dt+s+mu)h2 + mu(xi-xi^2)dxi h2 + 2mu(h2-h3)'] = zero(c2 - alvo)

# Controle negativo de I3: alterar o termo quadratico de Omega3 para depender de sigma.
def Omegas_alterado(w, sig, Dop):
    O = Omegas(w, sig, Dop)
    O['3'] = O['3'] + sp.Rational(1, 4)*xi*(w[('w2', sig)]*w['w1'])
    return O


wp_ = {k: w[k] + eps*sp.exp(s*tau)*h[k] for k in w}
dalt = {sg: sp.expand(sp.diff(Omegas_alterado(wp_, sg, lambda S, f: D(S, f))['3'], eps).subs(eps, 0)
                      / sp.exp(s*tau)) for sg in (1, -1)}
hp3, hm3 = h[('w3', 1)], h[('w3', -1)]
a3 = ((1 - xi)*hp3 - (1 + xi)*hm3)/2
b3 = (hp3 + hm3)/(2*mu)
comp3 = -2*mu*xi*(sp.diff(b3, tau) + s*b3 - sp.diff(a3, xi))
resultados['CN I3 com Omega3 alterado (deve ser False)'] = zero(comp3 - (dalt[-1] - dalt[1]))

ok = True
for k, v in resultados.items():
    esperado = 'deve ser False' not in k
    marca = 'OK ' if v == esperado else 'FALHA'
    ok &= (v == esperado)
    print(f'{marca} {k}: {v}')
print('TODAS AS VERIFICACOES CONFEREM' if ok else 'HA DISCREPANCIAS')
