"""Verificacao simbolica UNIVERSAL (exata) de identidades off-shell de RT.

Auditoria fisica de 16/09/2026. Transcricao feita diretamente de
.cache/rt-1203.3766v1 (choptuik.tex, secoes 2.2-2.3 e 2.4):
  eezuuirzr (Omega^sigma 2D), fkdkkdkkkeieiruus (Omega^+ Fourier-Chebyshev),
  Gamma1, Gamma2 (texto), a identidade diferencial homogenea antes de
  kdjskkkkksjsjhd, e Gamma2^sharp. NAO usa as tabelas dos macros C nem
  scripts/sharp_propagation.py (exceto em check_macro_tables, que as COMPARA).

Metodo. Campos sao polinomios em (tau, xi) com coeficientes SIMBOLICOS; s e mu
sao variaveis polinomiais. Uma identidade de grau d nos campos, com derivadas
de ordem <= k, envolve apenas os k-jatos dos campos nos pontos (tau0, +xi0) e
(tau0, -xi0) (a reflexao P e o unico termo nao local). Com a base
  k=2: tau^a xi^b, a<=2, b<=5-2a ;  k=1: b<=3-2a (a<=1) ;  k=0: b<=1, a=0
(somente b impares para a componente impar), interpolacao de Hermite realiza
QUALQUER k-jato nos dois pontos (xi0 != 0). Como as identidades nao tem tau
explicito, anular o polinomio para todos os coeficientes prova a identidade
para todos os campos C^k, todo s e todo mu. Nenhuma divisao por xi e usada
sem antes provar divisibilidade.
"""
from fractions import Fraction as Fr
import itertools
import sys
import time

CAP = [3]  # grau simbolico maximo mantido nas multiplicacoes


class Poly:
    """Chaves (t, x, s, m, syms): tau^t xi^x s^s mu^m prod c_syms."""
    __slots__ = ('d',)

    def __init__(self, d=None):
        self.d = {} if d is None else {k: v for k, v in d.items() if v}

    @staticmethod
    def const(c):
        return Poly({(0, 0, 0, 0, ()): Fr(c)})

    def __add__(self, o):
        o = o if isinstance(o, Poly) else Poly.const(o)
        r = dict(self.d)
        for k, v in o.d.items():
            nv = r.get(k, 0) + v
            if nv:
                r[k] = nv
            else:
                r.pop(k, None)
        p = Poly()
        p.d = r
        return p

    __radd__ = __add__

    def __neg__(self):
        p = Poly()
        p.d = {k: -v for k, v in self.d.items()}
        return p

    def __sub__(self, o):
        return self + (-(o if isinstance(o, Poly) else Poly.const(o)))

    def __rsub__(self, o):
        return Poly.const(o) - self

    def __mul__(self, o):
        if not isinstance(o, Poly):
            c = Fr(o)
            p = Poly()
            p.d = {k: v * c for k, v in self.d.items()} if c else {}
            return p
        cap = CAP[0]
        r = {}
        for k1, v1 in self.d.items():
            l1 = len(k1[4])
            for k2, v2 in o.d.items():
                if l1 + len(k2[4]) > cap:
                    continue
                syms = k1[4] + k2[4]
                if len(k1[4]) and len(k2[4]):
                    syms = tuple(sorted(syms))
                key = (k1[0] + k2[0], k1[1] + k2[1], k1[2] + k2[2], k1[3] + k2[3], syms)
                nv = r.get(key, 0) + v1 * v2
                if nv:
                    r[key] = nv
                else:
                    r.pop(key, None)
        p = Poly()
        p.d = r
        return p

    __rmul__ = __mul__

    def dt(self):
        return Poly({(k[0] - 1,) + k[1:]: v * k[0] for k, v in self.d.items() if k[0]})

    def dx(self):
        return Poly({(k[0], k[1] - 1) + k[2:]: v * k[1] for k, v in self.d.items() if k[1]})

    def P(self):
        return Poly({k: (-v if k[1] % 2 else v) for k, v in self.d.items()})

    def xi0_part(self):
        return Poly({k: v for k, v in self.d.items() if k[1] == 0})

    def div_x_exact(self):
        if any(k[1] == 0 for k in self.d):
            raise ArithmeticError('divisao por xi de polinomio nao divisivel')
        return Poly({(k[0], k[1] - 1) + k[2:]: v for k, v in self.d.items()})

    def is_zero(self):
        return not self.d

    def nterms(self):
        return len(self.d)


ONE = Poly.const(1)
X = Poly({(0, 1, 0, 0, ()): Fr(1)})
S = Poly({(0, 0, 1, 0, ()): Fr(1)})
MU = Poly({(0, 0, 0, 1, ()): Fr(1)})


def half(p):
    return p * Fr(1, 2)


def even(p):
    return half(p + p.P())


def odd(p):
    return half(p - p.P())


# ---------------------------------------------------------------- bases
def jet_basis(order, odd_only):
    table = {2: [(0, 5), (1, 3), (2, 1)], 1: [(0, 3), (1, 1)], 0: [(0, 1)]}[order]
    out = []
    for a, bmax in table:
        for b in range(bmax + 1):
            if odd_only and b % 2 == 0:
                continue
            out.append((a, b))
    return out


class SymbolSource:
    def __init__(self):
        self.n = 0

    def field(self, order, odd_only):
        p = Poly()
        for a, b in jet_basis(order, odd_only):
            p.d[(a, b, 0, 0, (self.n,))] = Fr(1)
            self.n += 1
        return p

    def multifield(self, order, first_odd=True):
        return [self.field(order, first_odd and i == 0) for i in range(4)]


# ------------------------------------------------------- RT, texto TeX
def D_sigma_plus_sigma_mu(f, sigma, shift):
    """(D^sigma + sigma mu) f, com d/dtau -> d/dtau + shift so na perturbacao."""
    tpart = f.dt() + (shift * f if shift is not None else Poly())
    return sigma * tpart + MU * (ONE + sigma * X) * f.dx() + sigma * MU * f


def omega_2d(sigma, w, shift=None, quad_only=False, lin_only=False):
    """eezuuirzr; w=(w1,w2,w3,w4), w_i^-=w_i, w_i^+=-P w_i (i=2,3,4)."""
    w1 = w[0]
    minus = {i: w[i - 1] for i in (2, 3, 4)}
    plus = {i: -w[i - 1].P() for i in (2, 3, 4)}
    om = lambda i, sg: (minus if sg < 0 else plus)[i]
    sg, ms = sigma, -sigma
    q = Fr(1, 4)
    lin = [
        X * D_sigma_plus_sigma_mu(w1, sg, shift) + 2 * MU * w1,
        X * D_sigma_plus_sigma_mu(om(2, ms), sg, shift) + MU * (w1 + om(2, ms) + om(2, sg)),
        X * D_sigma_plus_sigma_mu(om(2, sg), sg, shift) + MU * (2 * om(2, sg) - 2 * om(3, sg)),
        X * D_sigma_plus_sigma_mu(om(3, ms), sg, shift) + MU * (-w1),
        X * D_sigma_plus_sigma_mu(om(4, ms), sg, shift) + MU * (om(4, ms) + om(4, sg)),
    ]
    a1, a2m, a2p = w1, om(2, -1), om(2, 1)
    a2s, a2ms, a3s, a4s = om(2, sg), om(2, ms), om(3, sg), om(4, sg)
    a4m, a4p = om(4, -1), om(4, 1)
    quad = [
        X * q * (a1 * a1 - 2 * a2ms * a1 - 6 * a2s * a1 + 4 * a3s * a1 + a2ms * a2ms
                 - 2 * a2s * a2ms - 4 * a3s * a2ms + a2s * a2s + 4 * a3s * a2s - 4 * a4s * a4s),
        X * q * (a1 * a1 - 2 * a2m * a1 - 2 * a2p * a1 + a2m * a2m - 6 * a2p * a2m + a2p * a2p),
        X * q * (-4 * a2s * a2s + 8 * a3s * a2s - 4 * a4s * a4s),
        X * q * (-a1 * a1 + 2 * a2m * a1 + 2 * a2p * a1 - a2m * a2m + 2 * a2p * a2m
                 - a2p * a2p - 4 * a4p * a4m),
        X * q * (-4 * a4p * a2m - 4 * a4m * a2p),
    ]
    if quad_only:
        return quad
    if lin_only:
        return lin
    return [l + r for l, r in zip(lin, quad)]


def polarize(fun, w, h):
    """Parte bilinear B(w,h) de uma forma quadratica homogenea fun."""
    wh = [a + b for a, b in zip(w, h)]
    return [a - b - c for a, b, c in zip(fun(wh), fun(w), fun(h))]


def lin_omega_2d(sigma, w, h, shift):
    lin = omega_2d(sigma, h, shift, lin_only=True)
    bil = polarize(lambda v: omega_2d(sigma, v, quad_only=True), w, h)
    return [a + b for a, b in zip(lin, bil)]


def J(f, shift):
    tpart = f.dt() + (shift * f if shift is not None else Poly())
    return tpart + MU * ((ONE + X) * f.dx() + f)


def gamma1_tex(v):
    v1, v2, v3, v4 = v
    return [2 * v1, v1 + v2 - v2.P(), 2 * v3.P() - 2 * v2.P(), -v1, v4 - v4.P()]


def gamma2_tex(v, w):
    v1, v2, v3, v4 = v
    P = lambda f: f.P()
    cols = [
        (Fr(1, 4), [3 * P(v2) - 2 * P(v3) + v1 - v2, P(v2) + v1 - v2, Poly(), -P(v2) - v1 + v2, Poly()],
         w[0]),
        (Fr(1, 2), [P(v2) + 2 * P(v3) + v1 + v2, 2 * P(v2) + 2 * v2, -2 * P(v2) + 2 * P(v3),
                    -P(v2) - v2, P(v4) + v4], half(w[1] + P(w[1]))),
        (Fr(1, 2), [-2 * v1, P(v2) - v1 - v2, 2 * P(v2) - 2 * P(v3), v1, P(v4) - v4],
         half(w[1] - P(w[1]))),
        (Fr(1, 2), [P(v2) - v1 + v2, Poly(), 2 * P(v2), Poly(), Poly()], half(w[2] + P(w[2]))),
        (Fr(1, 2), [v1 - P(v2) - v2, Poly(), -2 * P(v2), Poly(), Poly()], half(w[2] - P(w[2]))),
        (Fr(1, 2), [-2 * P(v4), Poly(), -2 * P(v4), P(v4) + v4, P(v2) + v2], half(w[3] + P(w[3]))),
        (Fr(1, 2), [2 * P(v4), Poly(), 2 * P(v4), P(v4) - v4, P(v2) - v2], half(w[3] - P(w[3]))),
    ]
    out = [Poly() for _ in range(5)]
    for c, col, right in cols:
        for i in range(5):
            if col[i].d:
                out[i] = out[i] + c * col[i] * right
    return out


def omega_fc(v, shift=None, quad_only=False, lin_only=False):
    v1, v2, v3, v4 = v
    H = [J(v1, shift), J(v2, shift), -J(v2.P(), shift), J(v3, shift), J(v4, shift)]
    lin = [X * a + MU * b for a, b in zip(H, gamma1_tex(v))]
    quad = [X * g for g in gamma2_tex(v, v)]
    if quad_only:
        return quad
    if lin_only:
        return lin
    return [a + b for a, b in zip(lin, quad)]


def lin_omega_fc(w, h, shift):
    """D Omega^+(w)[h] com deslocamento: principal(h) + mu Gamma1 h + X(G2(w,h)+G2(h,w))."""
    lin = omega_fc(h, shift, lin_only=True)
    g1, g2 = gamma2_tex(w, h), gamma2_tex(h, w)
    return [a + X * (b + c) for a, b, c in zip(lin, g1, g2)]


def sharp_G_tex(w, Om):
    """(1/2){...} da identidade homogenea (secao 2.3), Om=(O1,O2,O2s,O3,O4)."""
    w1, w2, w3, w4 = w
    P = lambda f: f.P()
    O1, O2, O2s, O3, O4 = Om
    ev = lambda f: half(f + P(f))
    od = lambda f: half(f - P(f))
    terms = [
        ([2 * P(w2) - P(w3) + w1 - 2 * w2 + w3, P(w2) + w1 - w2], ev(O1)),
        ([P(w2) - P(w3) + w2 - w3, P(w2) + w1 - w2], od(O1)),
        ([-P(w3) - 3 * w1 + w3, -P(w2) - w1 - 7 * w2 + 4 * w3], ev(O2)),
        ([-P(w2) - P(w3) - w2 - w3, -P(w2) - w1 + w2 - 4 * w3], od(O2)),
        ([P(w3) - w1 - w3, 3 * P(w2) - w1 + w2], ev(O2s)),
        ([P(w2) + P(w3) + w2 + w3, 3 * P(w2) - w1 + w2], od(O2s)),
        ([2 * w1, 4 * w2], ev(O3)),
        ([-2 * P(w2) - 2 * w2, -4 * w2], od(O3)),
        ([2 * P(w4) - 2 * w4, -4 * w4], ev(O4)),
        ([2 * P(w4) + 2 * w4, 4 * w4], od(O4)),
    ]
    out = [Poly(), Poly()]
    for col, right in terms:
        for i in range(2):
            out[i] = out[i] + half(col[i] * right)
    return out


def sharp_T_tex(w, Om, shift=None):
    """Lado direito da identidade homogenea de RT (deve ser 0)."""
    O1, O2, O2s, O3, O4 = Om
    G = sharp_G_tex(w, Om)
    r1 = X * (odd(J(O1, shift)) + G[0]) + MU * even(O1)
    r2 = X * (J(O2s, shift) - J(O2, shift).P() + G[1]) + MU * (O1 + O2 + O2.P() - 2 * O3.P())
    return [r1, r2]


def gamma2_sharp_tex(v, c):
    v1, v2, v3, v4 = v
    P = lambda f: f.P()
    w1, w2s = c
    out0 = half((v2 + P(v2) - P(v3) - v3) * w1) + half((-v1 + P(v3) - v3) * half(w2s + P(w2s))) \
        + half((v2 + P(v2) + P(v3) + v3) * half(w2s - P(w2s)))
    out1 = half((v1 - v2 + P(v2)) * w1) + half((-v1 + v2 + 3 * P(v2)) * half(w2s + P(w2s))) \
        + half((-v1 + v2 + 3 * P(v2)) * half(w2s - P(w2s)))
    return [out0, out1]


def select(q):
    """S e S^sharp de RT; exige divisibilidade antes de dividir por xi."""
    F = [even(q[0]).div_x_exact() if q[0].d else Poly(),
         q[1].div_x_exact() if q[1].d else Poly(),
         q[3].div_x_exact() if q[3].d else Poly(),
         q[4].div_x_exact() if q[4].d else Poly()]
    C = [odd(q[0]), -q[2].P()]
    return F, C


# ------------------------------------------------------------- checks
def all_zero(polys):
    return all(p.is_zero() for p in polys)


def run(name, fn, degrees):
    ok = True
    for deg, order in degrees:
        CAP[0] = deg
        t0 = time.time()
        res = fn(order)
        good = all_zero(res)
        ok &= good
        size = sum(p.nterms() for p in res)
        print(f'  [{name}] grau<={deg} jatos={order}: {"OK" if good else "FALHOU"}'
              f' (termos residuais {size}, {time.time()-t0:.1f}s)')
    return ok


def check_gamma2_symmetry(order):
    src = SymbolSource()
    v, w = src.multifield(order), src.multifield(order)
    return [a - b for a, b in zip(gamma2_tex(v, w), gamma2_tex(w, v))]


def check_fc_equals_2d(order):
    src = SymbolSource()
    w = src.multifield(order)
    return [a - b for a, b in zip(omega_fc(w, S), omega_2d(+1, w, S))]


def check_reflection(order):
    src = SymbolSource()
    w = src.multifield(order)
    return [a.P() + b for a, b in zip(omega_2d(+1, w, S), omega_2d(-1, w, S))]


def check_rt_sharp_identity(order):
    src = SymbolSource()
    w = src.multifield(order)
    return sharp_T_tex(w, omega_2d(-1, w))


def check_sharp_restriction(order):
    """G(w, I c) == Gamma2^sharp(w, c) para c1 impar."""
    src = SymbolSource()
    w = src.multifield(order)
    c = [src.field(order, True), src.field(order, False)]
    full = sharp_G_tex(w, [c[0], Poly(), c[1], Poly(), Poly()])
    return [a - b for a, b in zip(full, gamma2_sharp_tex(w, c))]


def K_doc(w, c, shift):
    """SHARP_LINEARIZED_IDENTITY: (Kc)1=[dtau+s+mu(xi dxi+1)]c1+G_I1; (Kc)2=J_s c2+mu R_x c1+G_I2."""
    G = sharp_G_tex(w, [c[0], Poly(), c[1], Poly(), Poly()])
    k1 = c[0].dt() + shift * c[0] + MU * (X * c[0].dx() + c[0]) + G[0]
    k2 = J(c[1], shift) + MU * c[0].div_x_exact() + G[1] if c[0].d else J(c[1], shift) + G[1]
    return [k1, k2]


def M_doc(w, f, shift):
    """(Mf)1=-mu(xi dxi+2)P f1-(G_R)1; (Mf)2=P J_s(X P f2)-mu(Pf1+Pf2-f2+2f3)-(G_R)2."""
    f1, f2, f3, f4 = f
    Rf = [X * f1.P(), X * f2.P(), Poly(), X * f3.P(), X * f4.P()]
    G = sharp_G_tex(w, Rf)
    m1 = -MU * (X * f1.P().dx() + 2 * f1.P()) - G[0]
    m2 = J(X * f2.P(), shift).P() - MU * (f1.P() + f2.P() - f2 + 2 * f3) - G[1]
    return [m1, m2]


def check_linearized_KC_ML(order, with_background=True):
    src = SymbolSource()
    h = src.multifield(order)
    w = src.multifield(order) if with_background else [Poly() for _ in range(4)]
    dq = lin_omega_fc(w, h, S)
    F, C = select(dq)
    e = [-p.P() for p in omega_fc(w)]
    corr = sharp_G_tex(h, e)
    kc, ml = K_doc(w, C, S), M_doc(w, F, S)
    return [a - b + c for a, b, c in zip(kc, ml, corr)]


def check_C1_formula(order):
    """dC/ds = (0, -xi h2): coeficiente de s em C(s)h (M1=C1 segue da formula de M)."""
    src = SymbolSource()
    h = src.multifield(order)
    w = src.multifield(order)
    _, C = select(lin_omega_fc(w, h, S))
    coeff = lambda p: Poly({(k[0], k[1], k[2] - 1, k[3], k[4]): v for k, v in p.d.items() if k[2] == 1})
    higher = [Poly({k: v for k, v in p.d.items() if k[2] > 1}) for p in C]
    return [coeff(C[0]), coeff(C[1]) + X * h[1]] + higher


def muZ1(f):
    """mu (Z+1) f = dtau f + mu xi dxi f + mu f."""
    return f.dt() + MU * X * f.dx() + MU * f


def check_gauge_covariance(order):
    src = SymbolSource()
    w = src.multifield(order)
    h = [muZ1(p) for p in w]
    lhs = lin_omega_fc(w, h, MU)  # s = mu
    rhs = [muZ1(p) for p in omega_fc(w)]
    return [a - b for a, b in zip(lhs, rhs)]


def check_gauge_selected(order):
    src = SymbolSource()
    w = src.multifield(order)
    h = [muZ1(p) for p in w]
    Fl, Cl = select(lin_omega_fc(w, h, MU))
    F0, C0 = select(omega_fc(w))
    out = [a - (muZ1(b) + MU * b) for a, b in zip(Fl, F0)]
    out += [a - muZ1(b) for a, b in zip(Cl, C0)]
    return out


def check_phase_covariance(order):
    src = SymbolSource()
    w = src.multifield(order)
    h = [p.dt() for p in w]
    return [a - b.dt() for a, b in zip(lin_omega_fc(w, h, None), omega_fc(w))]


def check_reconstruction_compatibility(order):
    """-xi[T_s(h^+ + h^-) - 2mu d_xi a] = dOmega_j^- - dOmega_j^+, (i,j)=(2,2),(3,3),(4,4)."""
    src = SymbolSource()
    w, h = src.multifield(order), src.multifield(order)
    dm, dp = lin_omega_2d(-1, w, h, S), lin_omega_2d(+1, w, h, S)
    out = []
    for i, j in ((2, 1), (3, 3), (4, 4)):  # indices em (O1,O2,O2s,O3,O4)
        hm, hp = h[i - 1], -h[i - 1].P()
        a = half((ONE - X) * hp - (ONE + X) * hm)
        tb2mu = (hp + hm).dt() + S * (hp + hm)
        out.append(-X * (tb2mu - 2 * MU * a.dx()) - (dm[j] - dp[j]))
    return out


def check_W_identity(order):
    src = SymbolSource()
    w = src.multifield(order)
    W = MU + half(X * (w[0] - w[1] + w[1].P()))
    out = []
    for sg in (-1, 1):
        o = omega_2d(sg, w)
        w2 = w[1] if sg < 0 else -w[1].P()
        w3 = w[2] if sg < 0 else -w[2].P()
        DW = sg * W.dt() + MU * (ONE + sg * X) * W.dx()
        out.append(DW - W * (w2 - w3) - half(o[0] - o[1] - o[2]))
    return out


def check_macro_tables():
    """Compara sharp_propagation.gamma2/gamma_sharp_full (tabelas dos macros C)
    com a transcricao TeX, em todos os pares da base de 0-jatos (prova completa
    para formas bilineares algebricas)."""
    import sharp_propagation as sp
    basis_generic = [(0, 0), (0, 1)]
    basis_odd = [(0, 1)]

    def mono_sp(a, b):
        return sp.Poly({(a, b, 0): 1})

    def mono_me(a, b):
        return Poly({(a, b, 0, 0, ()): Fr(1)})

    def to_me(p):
        return Poly({(k[0], k[1], 0, 0, ()): Fr(v) for k, v in p.terms.items()})

    fields = []
    for comp in range(4):
        for (a, b) in (basis_odd if comp == 0 else basis_generic):
            fields.append((comp, a, b))
    bad = 0
    CAP[0] = 3
    for (ci, ai, bi), (cj, aj, bj) in itertools.product(fields, fields):
        vs = [sp.Poly() for _ in range(4)]
        ws = [sp.Poly() for _ in range(4)]
        vm = [Poly() for _ in range(4)]
        wm = [Poly() for _ in range(4)]
        vs[ci], vm[ci] = mono_sp(ai, bi), mono_me(ai, bi)
        ws[cj], wm[cj] = mono_sp(aj, bj), mono_me(aj, bj)
        g_sp = [to_me(p) for p in sp.gamma2(vs, ws)]
        g_me = gamma2_tex(vm, wm)
        bad += sum(not (x - y).is_zero() for x, y in zip(g_sp, g_me))
    # identidade sharp completa: residuo generico de cinco componentes
    for comp in range(4):
        for (a, b) in (basis_odd if comp == 0 else basis_generic):
            for rc in range(5):
                for (ra, rb) in basis_generic:
                    ws = [sp.Poly() for _ in range(4)]
                    wm = [Poly() for _ in range(4)]
                    ws[comp], wm[comp] = mono_sp(a, b), mono_me(a, b)
                    es = [sp.Poly() for _ in range(5)]
                    em = [Poly() for _ in range(5)]
                    es[rc], em[rc] = mono_sp(ra, rb), mono_me(ra, rb)
                    g_sp = [to_me(p) for p in sp.gamma_sharp_full(ws, es)]
                    g_me = sharp_G_tex(wm, em)
                    bad += sum(not (x - y).is_zero() for x, y in zip(g_sp, g_me))
    return bad


CHECKS = [
    ('Gamma2 (TeX) simetrico', check_gamma2_symmetry, [(2, 0)]),
    ('Omega^+ FC == Omega^+ 2D (f23ljf3lfj)', check_fc_equals_2d, [(1, 1), (2, 0)]),
    ('P Omega^+ = -Omega^-', check_reflection, [(1, 1), (2, 0)]),
    ('identidade sharp homogenea RT off-shell', check_rt_sharp_identity, [(1, 2), (2, 1), (3, 0)]),
    ('G(w,Ic) == Gamma2^sharp(w,c)', check_sharp_restriction, [(2, 0)]),
    ('K(s)C(s)h - M(s)L(s)h + G(h,e(w)) = 0', check_linearized_KC_ML, [(1, 2), (2, 1), (3, 0)]),
    ('dC/ds = (0,-xi h2)  (logo M1=C1)', check_C1_formula, [(1, 1), (2, 0)]),
    ('gauge: DOmega[mu(Z+1)w] (s=mu) = mu(Z+1)Omega', check_gauge_covariance, [(1, 2), (2, 1)]),
    ('gauge: L(mu)(Z+1)w=(Z+2)F, C(mu)(Z+1)w=(Z+1)c', check_gauge_selected, [(1, 2), (2, 1)]),
    ('fase: DOmega[dtau w] (s=0) = dtau Omega', check_phase_covariance, [(1, 2), (2, 1)]),
    ('reconstrucao: compatibilidade (2,2),(3,3),(4,4)', check_reconstruction_compatibility, [(1, 1), (2, 0)]),
    ('W: D^s W - W(w2-w3) = (O1-O2-O2s)/2', check_W_identity, [(1, 1), (2, 0)]),
]


def main(argv):
    quick = '--quick' in argv
    ok = True
    for name, fn, degrees in CHECKS:
        if quick:
            degrees = [(d, 0) for d, _ in degrees]
        ok &= run(name, fn, degrees)
    if not quick:
        sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent))
        bad = check_macro_tables()
        print(f'  [tabelas dos macros RT == TeX] divergencias: {bad}')
        ok &= bad == 0
    print('TODAS OK' if ok else 'HA FALHAS')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
