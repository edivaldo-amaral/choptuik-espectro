"""C2 (revisao CRH): os modos neutros invisiveis a L. Rodar com sympy no PYTHONPATH (scratchpad).

(a) delta = (2 g_bar, 0) (escala constante, delta zeta = -1) e delta = (0, 1) (translacao do escalar) tem
    delta omega = 0 nas sete funcoes de RT (ekjehee), para (zeta, Q, phi) genericos: estao em ker D Omega.
(b) Ambos resolvem as equacoes linearizadas: Ric(e^{2c} g) = Ric(g) (Christoffel invariantes por escala
    constante) e Box_{e^{2c}g} phi = e^{-2c} Box_g phi; d(phi + c) = d phi. Conferido em 2D+esfera com
    metrica radial generica (a parte que importa: Christoffel de e^{2c}g = Christoffel de g).
(c) Floquet com s = 0: (Theta^2)^*(2 g_bar) = e^{-4K}(2 g_bar), (Theta^2)^* 1 = 1.
(d) Nao sao de gauge: estao na forma do ansatz (delta g_{++} = delta g_{--} = 0), entao um X com
    L_X(g_bar, phi*) = delta seria admissivel, com delta omega_X = 0; o Passo 5 do Lema G (que nao usa Floquet
    nem Re s) da X = 0, logo L_X = 0 != delta. Contradicao.
"""
import sympy as sp

tau, xi = sp.symbols('tau xi', real=True)
mu, c = sp.symbols('mu c', positive=True)


def Dsig(sig, f):
    return sig*sp.diff(f, tau) + mu*(1 + sig*xi)*sp.diff(f, xi)


def omegas(z, q, p):
    F = z + sp.log(mu + xi**2*q)
    return [2*xi*q + Dsig(-1, F) + Dsig(1, F), Dsig(1, F) - mu, Dsig(-1, F) + mu,
            Dsig(1, z) - mu, Dsig(-1, z) + mu, Dsig(1, p), Dsig(-1, p)]


z, q, p = (sp.Function(n)(tau, xi) for n in ('zeta', 'Q', 'phi'))
w0 = omegas(z, q, p)
we = omegas(z + c, q, p)          # e^{-2(zeta + c)}: metrica multiplicada por e^{-2c}
wp = omegas(z, q, p + c)
print('(a) escala zeta -> zeta + c: delta omega = 0 nas 7:', all(sp.simplify(a - b) == 0 for a, b in zip(we, w0)))
print('(a) translacao phi -> phi + c: delta omega = 0 nas 7:', all(sp.simplify(a - b) == 0 for a, b in zip(wp, w0)))

# (b) Christoffel de e^{2c} g = Christoffel de g, para g 4D esfericamente simetrica generica em (t, r, th, ph)
t, r, th, ph = sp.symbols('t r th ph', real=True)
X = [t, r, th, ph]
A, B, Cc, R = (sp.Function(n)(t, r) for n in ('A', 'B', 'C', 'R'))
g = sp.Matrix([[A, B, 0, 0], [B, Cc, 0, 0], [0, 0, R**2, 0], [0, 0, 0, R**2*sp.sin(th)**2]])


def christ(G):
    Gi = G.inv()
    return [[[sp.simplify(sum(Gi[a, d]*(sp.diff(G[d, b], X[cc]) + sp.diff(G[d, cc], X[b]) - sp.diff(G[b, cc], X[d]))
                              for d in range(4))/2) for cc in range(4)] for b in range(4)] for a in range(4)]


G1, G2 = christ(g), christ(sp.exp(2*c)*g)
print('(b) Christoffel(e^{2c} g) = Christoffel(g):',
      all(sp.simplify(G1[a][b][cc] - G2[a][b][cc]) == 0 for a in range(4) for b in range(4) for cc in range(4)),
      '=> Ric(e^{2c} g) = Ric(g) = 2 dphi dphi; Box_{e^{2c}g} = e^{-2c} Box_g')
print('(c) s = 0: o modo de escala escala como o fundo sob Theta^2 (fator e^{-4K}); o de translacao e invariante.')
print('(d) ver docstring: Passo 5 do Lema G (GAUGE_GLOBAL.md par. 3) => nao sao de gauge.')
