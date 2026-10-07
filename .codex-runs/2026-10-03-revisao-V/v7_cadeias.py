#!/usr/bin/env python3
"""Revisao V (03/10/2026), V5/V7 com funcoes genericas (sympy).

V7 (cadeias): para dzeta = e^{s tau}(z1 + tau z0), D^sigma dzeta = e^{s tau}(h1^sigma + tau h0^sigma) com
   h0^sigma = D_s^sigma z0,  h1^sigma = D_s^sigma z1 + sigma z0,
e a recorrencia do documento T_s z1 = a(h1) - z0, d_xi z1 = b(h1) vale; a compatibilidade da ordem 1 e
   T_s b(h1) - d_xi a(h1) + b(h0) = 0, que e o coeficiente de tau^0 de (dOmega_3^- - dOmega_3^+)/(-2 mu xi)
   para a perturbacao secular. Tambem: o termo secular da selecao/constraint e a derivada em s
   (L-afim), isto e, a constraint de jato e C(s0) h1 + C1 h0 = 0, nao C(s0) h1 = 0.
V5: delta g^{-1} pelo frame de RT em coordenadas cartesianas (formula do documento) = derivada
   de g^{-1}(Q + eps e^{s tau} q, zeta + eps e^{s tau} z); e (Theta^2)^* delta g = e^{4 pi s} e^{-4K} delta g
   quando zeta(tau + 4pi) = zeta + 2K e Q, q, z sao 4pi-periodicos.
"""
import sympy as sp
from pathlib import Path

src = Path(__file__).with_name('v123_identidades.py').read_text()
ns = {}
exec(src.split('resultados = {}')[0], ns)
Omegas, D, tau, xi, mu, s, eps, F = (ns[k] for k in ('Omegas', 'D', 'tau', 'xi', 'mu', 's', 'eps', 'F'))
zero = lambda e: sp.simplify(sp.expand(e)) == 0
res = {}

# ---- V7: jatos ----
z0, z1 = F('z0'), F('z1')
E = sp.exp(s*tau)
for sig in (1, -1):
    full = sp.expand(D(sig, E*(z1 + tau*z0))/E)
    h0 = D(sig, z0, s)
    h1 = D(sig, z1, s) + sig*z0
    res[f'D^{sig:+d} do jato = h1 + tau h0'] = zero(full - (h1 + tau*h0))
hp1, hm1 = D(1, z1, s) + z0, D(-1, z1, s) - z0
hp0, hm0 = D(1, z0, s), D(-1, z0, s)
a = lambda hp, hm: ((1 - xi)*hp - (1 + xi)*hm)/2
b = lambda hp, hm: (hp + hm)/(2*mu)
res['recorrencia T_s z1 = a(h1) - z0'] = zero(sp.diff(z1, tau) + s*z1 - (a(hp1, hm1) - z0))
res['d_xi z1 = b(h1)'] = zero(sp.diff(z1, xi) - b(hp1, hm1))
res['compat. ordem 1: T_s b(h1) - d_xi a(h1) + b(h0) = 0'] = zero(
    sp.diff(b(hp1, hm1), tau) + s*b(hp1, hm1) - sp.diff(a(hp1, hm1), xi) + b(hp0, hm0))
# a mesma compatibilidade sai de dOmega_3^- - dOmega_3^+ da perturbacao secular, off-shell
w = {'w1': F('w1')}
hh0 = {'w1': F('k1')}
hh1 = {'w1': F('l1')}
for i in (2, 3, 4):
    for sg, nm in ((1, 'p'), (-1, 'm')):
        w[(f'w{i}', sg)] = F(f'w{i}{nm}')
        hh0[(f'w{i}', sg)] = F(f'k{i}{nm}')
        hh1[(f'w{i}', sg)] = F(f'l{i}{nm}')
pert = {k: w[k] + eps*E*(hh1[k] + tau*hh0[k]) for k in w}
d3 = {sg: sp.expand(sp.diff(Omegas(pert, sg, lambda S, f: D(S, f))['3'], eps).subs(eps, 0)/E) for sg in (1, -1)}
dif = sp.expand(d3[-1] - d3[1])
A1 = a(hh1[('w3', 1)], hh1[('w3', -1)]); B1 = b(hh1[('w3', 1)], hh1[('w3', -1)])
B0 = b(hh0[('w3', 1)], hh0[('w3', -1)])
A0 = a(hh0[('w3', 1)], hh0[('w3', -1)])
# dif = -2 mu xi [(T_s b1 - d_xi a1 + b0) + tau (T_s b0 - d_xi a0)]  (sem usar coeff, pois as funcoes dependem de tau)
res['dOmega3^- - dOmega3^+ do jato = -2 mu xi [(T_s b1 - d_xi a1 + b0) + tau(T_s b0 - d_xi a0)]'] = zero(
    dif + 2*mu*xi*((sp.diff(B1, tau) + s*B1 - sp.diff(A1, xi) + B0) + tau*(sp.diff(B0, tau) + s*B0 - sp.diff(A0, xi))))
# termo secular de todos os Omega^+ = derivada em s do linearizado (L-afim): coef tau^0 = dO(h1) + d/ds dO(h0)
def lin(hv):
    p = {k: w[k] + eps*E*hv[k] for k in w}
    return {k: sp.expand(sp.diff(v, eps).subs(eps, 0)/E) for k, v in Omegas(p, 1, lambda S, f: D(S, f)).items()}


L1, L0 = lin(hh1), lin(hh0)
for k in L1:
    full = sp.expand(sp.diff(Omegas(pert, 1, lambda S, f: D(S, f))[k], eps).subs(eps, 0)/E)
    res[f'dOmega_{k}^+ do jato = [dO(h1) + d_s dO(h0)] + tau dO(h0)'] = zero(full - (L1[k] + sp.diff(L0[k], s)) - tau*L0[k])

# ---- V5: metrica ----
x1, x2, x3, K = sp.symbols('x1 x2 x3 K', real=True)
x = [x1, x2, x3]
r = sp.sqrt(x1**2 + x2**2 + x3**2)
Qf, Zf, qf, zf = (sp.Function(n) for n in ('Q', 'zeta', 'q', 'z'))


def frame(Qv):
    e0 = sp.Matrix([1, mu*x1, mu*x2, mu*x3])
    es = []
    for i in range(3):
        v = [0, 0, 0, 0]
        v[i + 1] += mu
        for k in range(3):
            v[i + 1] += Qv*x[k]*x[k]
            v[k + 1] -= Qv*x[k]*x[i]
        es.append(sp.Matrix(v))
    return e0, es


def ginv(Qv, Zv):
    e0, es = frame(Qv)
    return sp.exp(2*Zv)*(-e0*e0.T + sum((e*e.T for e in es), sp.zeros(4, 4)))


Qv, Zv, qv, zv = Qf(tau, r), Zf(tau, r), qf(tau, r), zf(tau, r)
deriv = sp.diff(ginv(Qv + eps*E*qv, Zv + eps*E*zv), eps).subs(eps, 0)
e0, es = frame(Qv)
de = []
for i in range(3):
    v = [0, 0, 0, 0]
    for k in range(3):
        v[i + 1] += qv*x[k]*x[k]
        v[k + 1] -= qv*x[k]*x[i]
    de.append(sp.Matrix(v))
formula = E*sp.exp(2*Zv)*(2*zv*(-e0*e0.T + sum((e*e.T for e in es), sp.zeros(4, 4)))
                          + sum((de[i]*es[i].T + es[i]*de[i].T for i in range(3)), sp.zeros(4, 4)))
res['V5 delta g^{-1} do documento = derivada de g^{-1}'] = all(zero(c) for c in (deriv - formula))
# Floquet: zeta = Zp + K tau/(2pi) com Zp 4pi-periodico; delta g = -g delta g^{-1} g.
Zp = sp.Function('Zp')
Zv2 = Zp(tau, r) + K*tau/(2*sp.pi)
pt = {x2: 0, x3: 0}
g = ginv(Qv, Zv2).subs(pt).inv()
dg = -g*formula.subs(Zv, Zv2).subs(pt)*g
per = {Qf(tau, r): Qf(tau + 4*sp.pi, r)}
shift = dg.subs(tau, tau + 4*sp.pi)
# sob a periodicidade de Q, q, z, Zp: substitui f(tau+4pi, r) -> f(tau, r)
rp = r.subs(pt)
rep = {fn(tau + 4*sp.pi, rp): fn(tau, rp) for fn in (Qf, qf, zf, Zp)}
shift = shift.xreplace(rep)
fator = sp.exp(4*sp.pi*s)*sp.exp(-4*K)
res['V5 (Theta^2)^* delta g = e^{4pi s} e^{-4K} delta g (no ponto (x1,0,0), x1 generico)'] = all(zero(sp.simplify(c)) for c in (shift - fator*dg))

ok = True
for k, v in res.items():
    ok &= bool(v)
    print(f"{'OK ' if v else 'FALHA'} {k}")
print('TODAS AS VERIFICACOES CONFEREM' if ok else 'HA DISCREPANCIAS')
