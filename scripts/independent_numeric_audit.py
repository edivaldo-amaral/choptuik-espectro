"""Auditoria fisica independente (16/09/2026) -- DIAGNOSTICO EM PONTO FLUTUANTE.

Nada aqui e prova. Reimplementa, a partir do TeX de RT (secao 2.4), os
operadores Fourier-Chebyshev (float64) e le os arquivos de matriz/RefA com
parsers proprios (nao usa analyze_spectrum/analyze_sharp). Subcomandos:

  validate  compara o operador reimplementado com as matrizes exportadas
            (selecionada, constraints C0/C1, sharp) coluna a coluna;
  gauge     modo s~mu: sobreposicao com (Z+1)omega_RefA, identidade exata
            s_N-mu = -y^H r/(y^H v) e decomposicao do residuo em defeito do
            fundo vs truncamento;
  sharp     raiz ~0.401: constraint C(s)h versus autovetor sharp;
  census    espectro de L_A nas faixas A e imagem de B (inclusive fora do
            retangulo 1/8<=Re<=1);
  exact     transversalidade |b(dtau omega*)|>117/1000 com Fraction.
"""
from fractions import Fraction as Fr
import argparse
import math
import re
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import eig
from scipy.signal import fftconvolve

ROOT = Path(__file__).resolve().parents[1]
REFA = ROOT / '.cache/rt-1203.3766v1/sourcecode/RefA.dat'
SPEC = ROOT / 'build/spectrum'


# ----------------------------------------------------------------- parsers
def parse_refa(path=REFA):
    text = Path(path).read_text()
    lines = text.splitlines()
    assert lines[0].strip() == 'mu_MultiField'
    mu_exp = int(lines[2].split()[1])
    mu_coeff = int(lines[3].split()[1])
    mu = Fr(mu_coeff, 1) * Fr(2) ** mu_exp
    fields = []
    i = 6
    for _ in range(4):
        assert lines[i].strip() == 'Field'
        two = int(lines[i + 1].split()[1])
        off_m, off_n = int(lines[i + 2].split()[1]), int(lines[i + 3].split()[1])
        num_m, num_n = int(lines[i + 4].split()[1]), int(lines[i + 5].split()[1])
        assert off_m == 0 and off_n == 0
        pairs = re.findall(r'\((-?\d+),(-?\d+)\)', lines[i + 6])
        assert len(pairs) == num_m * num_n
        exact = {}
        for idx, (re_, im_) in enumerate(pairs):
            m, n = divmod(idx, num_n)
            if int(re_) or int(im_):
                exact[m, n] = (int(re_), int(im_), two)
        fields.append(dict(num_m=num_m, num_n=num_n, exact=exact))
        i += 7
    return mu, fields


def read_rt_matrix(path, magic):
    with open(path) as fh:
        assert fh.readline().strip() == magic, path
        dim = int(fh.readline().split()[1])
        fh.readline()
        fh.readline()
        assert fh.readline().strip() == 'BEGIN_ADDRESSES'
        addr = []
        for j in range(dim):
            idx, d, k, m, n = map(int, fh.readline().split())
            assert idx == j
            addr.append((d, k, m, n))
        assert fh.readline().strip() == 'BEGIN_ENTRIES'
        A = np.zeros((dim, dim))
        for line in fh:
            f = line.split()
            coeff = int(f[2])
            e = int(f[3]) if len(f) == 4 else 0
            r, c = int(f[0]), int(f[1])
            assert A[r, c] == 0
            A[r, c] = math.ldexp(float(coeff), e) if abs(coeff) < 2 ** 1000 else float(Fr(coeff) * Fr(2) ** e)
    return A, addr


def read_constraints(path, columns):
    with open(path) as fh:
        assert fh.readline().strip() == 'CHOPTUIK_RT_SPECTRAL_CONSTRAINTS_V1'
        rows = int(fh.readline().split()[1])
        assert int(fh.readline().split()[1]) == columns
        assert fh.readline().strip() == 'BEGIN_ADDRESSES'
        addr = []
        for j in range(rows):
            idx, d, k, m, n = map(int, fh.readline().split())
            addr.append((d, k, m, n))
        assert fh.readline().strip() == 'BEGIN_ENTRIES'
        C0, C1 = np.zeros((rows, columns)), np.zeros((rows, columns))
        for line in fh:
            f = line.split()
            deg, r, c, coeff = int(f[0]), int(f[1]), int(f[2]), int(f[3])
            e = int(f[4]) if len(f) == 5 else 0
            (C0 if deg == 0 else C1)[r, c] = math.ldexp(float(coeff), e)
    return C0, C1, addr


# ------------------------------------------------ campos Fourier-Chebyshev
class Fld:
    """v_{mn}, m in [-M,M], n in [0,N]; TF(v)=sum_m e^{i m tau/2}[v_m0+2 sum_n v_mn T_n]."""

    def __init__(self, arr):
        self.a = arr

    @property
    def M(self):
        return (self.a.shape[0] - 1) // 2

    @property
    def N(self):
        return self.a.shape[1] - 1

    @staticmethod
    def zeros(M, N):
        return Fld(np.zeros((2 * M + 1, N + 1), dtype=complex))

    def resized(self, M, N):
        out = np.zeros((2 * M + 1, N + 1), dtype=complex)
        mm, nn = min(M, self.M), min(N, self.N)
        out[M - mm:M + mm + 1, :nn + 1] = self.a[self.M - mm:self.M + mm + 1, :nn + 1]
        return Fld(out)

    def __add__(self, o):
        if isinstance(o, (int, float)) and o == 0:
            return self
        M, N = max(self.M, o.M), max(self.N, o.N)
        return Fld(self.resized(M, N).a + o.resized(M, N).a)

    __radd__ = __add__

    def __sub__(self, o):
        return self + o * (-1)

    def __mul__(self, c):
        if isinstance(c, Fld):
            return conv(self, c)
        return Fld(self.a * c)

    __rmul__ = __mul__

    def P(self):
        return Fld(self.a * ((-1.0) ** np.arange(self.N + 1))[None, :])

    def dt(self):
        m = np.arange(-self.M, self.M + 1)
        return Fld(self.a * (0.5j * m)[:, None])

    def times_xi(self):
        v = np.concatenate([self.a, np.zeros((self.a.shape[0], 2), complex)], axis=1)
        out = np.zeros_like(v)
        out[:, 0] = v[:, 1]
        out[:, 1:-1] = 0.5 * (v[:, :-2] + v[:, 2:])
        return Fld(out[:, :-1])

    def one_plus_xi_dxi(self):
        n = np.arange(self.N + 1)
        kv = self.a * n[None, :]
        tail = np.cumsum(kv[:, ::-1], axis=1)[:, ::-1]  # sum_{k>=n} k v_k
        above = np.concatenate([tail[:, 1:], np.zeros((tail.shape[0], 1))], axis=1)
        return Fld(kv + 2 * above)

    def xi_dxi(self):
        n = np.arange(self.N + 1)
        kv = self.a * n[None, :]
        out = kv.copy()
        for parity in (0, 1):
            sel = kv[:, parity::2]
            tail = np.cumsum(sel[:, ::-1], axis=1)[:, ::-1]
            above = np.concatenate([tail[:, 1:], np.zeros((tail.shape[0], 1))], axis=1)
            out[:, parity::2] += 2 * above
        return Fld(out)

    def div_xi(self):
        """R_x: (v - v|_{xi=0})/xi; (R_x v)_n = 2 sum_L (-1)^L v_{n+1+2L}."""
        N = self.N
        out = np.zeros_like(self.a)
        for parity in (0, 1):
            idx = np.arange(1 + parity, N + 1, 2)  # n+1 com n de paridade `parity`
            if len(idx) == 0:
                continue
            sel = self.a[:, idx] * ((-1.0) ** np.arange(len(idx)))[None, :]
            tail = np.cumsum(sel[:, ::-1], axis=1)[:, ::-1]
            # para n = idx[j]-1: 2 sum_{L>=0} (-1)^L v_{idx[j+L]} = 2 (-1)^j tail_j
            out[:, idx - 1] = 2 * tail * ((-1.0) ** np.arange(len(idx)))[None, :]
        return Fld(out)

    def value(self, m, n):
        if abs(m) > self.M or n > self.N:
            return 0j
        return self.a[m + self.M, n]


CAPS = [120, 400]


def conv(u, v):
    def sym(f):
        return np.concatenate([f.a[:, :0:-1], f.a], axis=1)  # n = -N..N
    c = fftconvolve(sym(u), sym(v), mode='full')
    N = u.N + v.N
    out = Fld(c[:, N:])
    return out.resized(min(out.M, CAPS[0]), min(out.N, CAPS[1]))


def O_A(f, mu, shift=0.0):
    return f.dt() + f * shift + mu * (f.one_plus_xi_dxi() + f)


def O_B(f, mu, shift=0.0):
    return 0.5 * (O_A(f, mu, shift) + O_A(f.P(), mu, shift).P())


def gamma1(v):
    v1, v2, v3, v4 = v
    return [2 * v1, v1 + v2 - v2.P(), 2 * v3.P() - 2 * v2.P(), -1 * v1, v4 - v4.P()]


def gamma2(v, w):
    v1, v2, v3, v4 = v
    P = lambda f: f.P()
    Z = None
    cols = [
        (0.25, [3 * P(v2) - 2 * P(v3) + v1 - v2, P(v2) + v1 - v2, Z, -1 * P(v2) - v1 + v2, Z], w[0]),
        (0.5, [P(v2) + 2 * P(v3) + v1 + v2, 2 * P(v2) + 2 * v2, -2 * P(v2) + 2 * P(v3), -1 * P(v2) - v2,
               P(v4) + v4], 0.5 * (w[1] + P(w[1]))),
        (0.5, [-2 * v1, P(v2) - v1 - v2, 2 * P(v2) - 2 * P(v3), v1, P(v4) - v4], 0.5 * (w[1] - P(w[1]))),
        (0.5, [P(v2) - v1 + v2, Z, 2 * P(v2), Z, Z], 0.5 * (w[2] + P(w[2]))),
        (0.5, [v1 - P(v2) - v2, Z, -2 * P(v2), Z, Z], 0.5 * (w[2] - P(w[2]))),
        (0.5, [-2 * P(v4), Z, -2 * P(v4), P(v4) + v4, P(v2) + v2], 0.5 * (w[3] + P(w[3]))),
        (0.5, [2 * P(v4), Z, 2 * P(v4), P(v4) - v4, P(v2) - v2], 0.5 * (w[3] - P(w[3]))),
    ]
    out = [0] * 5
    for c, col, right in cols:
        for i in range(5):
            if col[i] is not None:
                out[i] = out[i] + c * conv(col[i], right)
    return [o if isinstance(o, Fld) else Fld.zeros(0, 0) for o in out]


def gamma2_sharp(v, c):
    v1, v2, v3, v4 = v
    P = lambda f: f.P()
    w1, w2 = c
    ev, od = 0.5 * (w2 + P(w2)), 0.5 * (w2 - P(w2))
    o0 = 0.5 * conv(v2 + P(v2) - P(v3) - v3, w1) + 0.5 * conv(-1 * v1 + P(v3) - v3, ev) \
        + 0.5 * conv(v2 + P(v2) + P(v3) + v3, od)
    o1 = 0.5 * conv(v1 - v2 + P(v2), w1) + 0.5 * conv(-1 * v1 + v2 + 3 * P(v2), ev) \
        + 0.5 * conv(-1 * v1 + v2 + 3 * P(v2), od)
    return [o0, o1]


def omega_plus_parts(mu, w, h, shift=0.0, background_quadratic=True):
    """X H_s h + mu Gamma1 h + X(2 Gamma2(w,h)) [linearizado] ou X Gamma2(w,w) se h is w."""
    h1, h2, h3, h4 = h
    H = [O_A(h1, mu, shift), O_A(h2, mu, shift), -1 * O_A(h2.P(), mu, shift),
         O_A(h3, mu, shift), O_A(h4, mu, shift)]
    g1 = gamma1(h)
    g2 = gamma2(w, h)
    factor = 2.0 if background_quadratic else 1.0
    return [a.times_xi() + mu * b + factor * c.times_xi() for a, b, c in zip(H, g1, g2)]


def selection_rt(q):
    """S = R_x o diag((1+P)/2,1,1,1) = diag((1-P)/2,1,1,1) o R_x (RT)."""
    F = [(0.5 * (q[0] + q[0].P())).div_xi(), q[1].div_xi(), q[3].div_xi(), q[4].div_xi()]
    C = [0.5 * (q[0] - q[0].P()), -1 * q[2].P()]
    return F, C


def refa_fields(mu_fr, fields, clip_m=None, clip_n=None):
    out = []
    for f in fields:
        M, N = f['num_m'] - 1, f['num_n'] - 1
        arr = np.zeros((2 * M + 1, N + 1), dtype=complex)
        for (m, n), (re_, im_, two) in f['exact'].items():
            if clip_m is not None and (m >= clip_m or n >= clip_n):
                continue
            z = complex(math.ldexp(float(re_), two), math.ldexp(float(im_), two))
            arr[M + m, n] = z
            arr[M - m, n] = z.conjugate()
        out.append(Fld(arr))
    return float(mu_fr), out


def basis_field(address, M, N):
    d, k, m, n = address
    comps = [Fld.zeros(M, N) for _ in range(4)]
    z = 1.0 if k == 0 else 1.0j
    comps[d].a[M + m, n] = z
    if m:
        comps[d].a[M - m, n] = np.conj(z)
    return comps


def dof_vector(fields, addresses):
    out = np.zeros(len(addresses))
    for j, (d, k, m, n) in enumerate(addresses):
        z = fields[d].value(m, n)
        out[j] = z.real if k == 0 else z.imag
    return out


def box_of(addresses):
    return 1 + max(a[2] for a in addresses), 1 + max(a[3] for a in addresses)


# -------------------------------------------------------------- comandos
def cmd_validate(args):
    mu_fr, raw = parse_refa()
    for tag in args.boxes:
        A, addr = read_rt_matrix(SPEC / f'rt-A-{tag}.dat', 'CHOPTUIK_RT_SPECTRAL_MATRIX_V1')
        C0, C1, caddr = read_constraints(SPEC / f'rt-A-{tag}.dat.constraints', len(A))
        Mb, Nb = box_of(addr)
        mu, w = refa_fields(mu_fr, raw, 2 * Mb, 2 * Nb)
        K = None
        sharp_path = SPEC / f'rt-sharp-{tag}.dat'
        if sharp_path.exists():
            K, kaddr = read_rt_matrix(sharp_path, 'CHOPTUIK_RT_SHARP_MATRIX_V1')
            assert kaddr == caddr
        cols = range(len(A)) if args.all_columns else np.random.default_rng(1).choice(len(A), min(40, len(A)), replace=False)
        errA = errA_trunc = errC0 = errC1 = 0.0
        scaleA = np.abs(A).max()
        for c in cols:
            h = basis_field(addr[c], Mb, Nb)
            q = omega_plus_parts(mu, w, h)
            F, C = selection_rt(q)
            colA = dof_vector(F, addr)
            errA = max(errA, np.abs(colA - A[:, c]).max())
            # variante: truncar produtos a caixa ANTES de dividir por xi
            qt = [p.resized(Mb - 1, Nb - 1) for p in q]
            Ft, _ = selection_rt(qt)
            errA_trunc = max(errA_trunc, np.abs(dof_vector(Ft, addr) - A[:, c]).max())
            errC0 = max(errC0, np.abs(dof_vector(C, caddr) - C0[:, c]).max())
            q1 = omega_plus_parts(mu, w, h, shift=1.0)
            _, Cs = selection_rt(q1)
            errC1 = max(errC1, np.abs(dof_vector(Cs, caddr) - dof_vector(C, caddr) - C1[:, c]).max())
        print(f'{tag}: dim={len(A)} colunas testadas={len(cols)} max|A|={scaleA:.3g}')
        print(f'  selecionada: erro (Pi L Pi, sem truncar intermediarios) = {errA:.3e};'
              f' variante truncada antes de R_x = {errA_trunc:.3e}')
        print(f'  constraints: erro C0 = {errC0:.3e}, erro C1 = {errC1:.3e}')
        if K is not None:
            errK = errKt = 0.0
            kc = range(len(K)) if args.all_columns else np.random.default_rng(2).choice(len(K), min(40, len(K)), replace=False)
            for c in kc:
                d, k, m, n = caddr[c]
                cf = [Fld.zeros(Mb, Nb), Fld.zeros(Mb, Nb)]
                z = 1.0 if k == 0 else 1.0j
                cf[d].a[Mb + m, n] = z
                if m:
                    cf[d].a[Mb - m, n] = np.conj(z)
                g = gamma2_sharp(w, cf)
                out = [O_B(cf[0], mu) + g[0], O_A(cf[1], mu) + mu * cf[0].div_xi() + g[1]]
                errK = max(errK, np.abs(dof_vector(out, caddr) - K[:, c]).max())
            print(f'  sharp: erro K(0) = {errK:.3e} (dim {len(K)})')


def eig_left_right(A):
    la, vl, vr = eig(A, left=True, right=True)
    return la, vl, vr


def weights(addresses, k1=129 / 128, k2=9 / 8):
    return np.array([(2 - (m == 0)) * (2 - (n == 0)) * k1 ** m * k2 ** n for _, _, m, n in addresses])


def cmd_gauge(args):
    mu_fr, raw = parse_refa()
    mu_ref = float(mu_fr)
    print(f'mu_RefA = {mu_fr} = {mu_ref:.17g}')
    for tag in args.boxes:
        A, addr = read_rt_matrix(SPEC / f'rt-A-{tag}.dat', 'CHOPTUIK_RT_SPECTRAL_MATRIX_V1')
        Mb, Nb = box_of(addr)
        mu, wc = refa_fields(mu_fr, raw, 2 * Mb, 2 * Nb)
        g_full = [f + (1.0 / mu) * f.dt() + f.xi_dxi() for f in wc]  # (Z+1) w_c, sem truncar
        v = dof_vector(g_full, addr)
        p = dof_vector([f.dt() for f in wc], addr)
        la, vl, vr = eig_left_right(A)
        s = -la
        i_mu = int(np.argmin(np.abs(s - mu)))
        i_0 = int(np.argmin(np.abs(s)))
        wts = weights(addr)
        print(f'\n== {tag} (dim {len(A)}), fundo recortado em {2*Mb}x{2*Nb}')
        for label, i, vec, target in (('s~mu', i_mu, v, mu), ('s~0 (fase)', i_0, p, 0.0)):
            x, y = vr[:, i], vl[:, i]
            cos = abs(np.vdot(vec, x)) / (np.linalg.norm(vec) * np.linalg.norm(x))
            shiftI = target
            r = A @ vec + shiftI * vec
            pred = -np.vdot(y, r) / np.vdot(y, vec)
            cond = np.linalg.norm(y) * np.linalg.norm(x) / abs(np.vdot(y, x))
            print(f'  {label}: s_N = {s[i].real:.12f}{s[i].imag:+.1e}i ; s_N - alvo = {(s[i]-target).real:+.3e}')
            print(f'     |cos|(gerador exato, autovetor) = {cos:.6f} ; cond autovalor = {cond:.1f}')
            print(f'     identidade -y^H r/y^H v = {pred.real:+.6e} (confere: {abs(pred-(s[i]-target)):.1e});'
                  f' ||r||/||v|| = {np.linalg.norm(r)/np.linalg.norm(vec):.3e}')
            if label == 's~mu' and args.decompose:
                # r = Pi (Z+2) F_{w_c} - Pi L (1-Pi)(Z+1)w_c  (se A = Pi L Pi)
                q = omega_plus_parts(mu, wc, wc, background_quadratic=False)
                F, _ = selection_rt(q)
                ZF = [f + (1.0 / mu) * f.dt() + f.xi_dxi() + f for f in F]  # (Z+2)F
                rdef = dof_vector(ZF, addr)
                outside = []
                for f in g_full:
                    o = Fld(f.a.copy())
                    o.a[f.M - (Mb - 1):f.M + Mb, :Nb] = 0.0
                    outside.append(o)
                qo = omega_plus_parts(mu, wc, outside)
                Fo, _ = selection_rt(qo)
                rtr = -dof_vector(Fo, addr)
                pd = -np.vdot(y, rdef) / np.vdot(y, v)
                pt = -np.vdot(y, rtr) / np.vdot(y, v)
                Fnorm = sum(np.abs(f.a).sum() for f in F)
                print(f'     ||r - (r_defeito + r_trunc)|| = {np.linalg.norm(r - rdef - rtr):.2e}')
                print(f'     parcela do DEFEITO do fundo: {pd.real:+.4e} ; parcela do TRUNCAMENTO: {pt.real:+.4e}')
                print(f'     ||S Omega^+(mu,RefA recortado)||_l1(coef) = {Fnorm:.3e}')


def cmd_sharp(args):
    for tag in args.boxes:
        A, addr = read_rt_matrix(SPEC / f'rt-A-{tag}.dat', 'CHOPTUIK_RT_SPECTRAL_MATRIX_V1')
        C0, C1, caddr = read_constraints(SPEC / f'rt-A-{tag}.dat.constraints', len(A))
        K, kaddr = read_rt_matrix(SPEC / f'rt-sharp-{tag}.dat', 'CHOPTUIK_RT_SHARP_MATRIX_V1')
        assert kaddr == caddr
        la, vl, vr = eig_left_right(A)
        ka, kl, kr = eig_left_right(K)
        s, sk = -la, -ka
        wk, wa = weights(kaddr), weights(addr)
        print(f'\n== {tag}: dim L={len(A)}, dim K={len(K)}')
        for target in [complex(t) for t in args.targets]:
            i = int(np.argmin(np.abs(s - target)))
            h = vr[:, i]
            c = (C0 + s[i] * C1) @ h
            j = int(np.argmin(np.abs(sk - s[i])))
            kvec = kr[:, j]
            cos2 = abs(np.vdot(kvec, c)) / (np.linalg.norm(kvec) * np.linalg.norm(c))
            cw = abs(np.vdot(kvec * np.sqrt(wk), c * np.sqrt(wk))) / (
                np.linalg.norm(kvec * np.sqrt(wk)) * np.linalg.norm(c * np.sqrt(wk)))
            res = (K + s[i] * np.eye(len(K))) @ c
            rel = np.sum(np.abs(res) * wk) / np.sum(np.abs(c) * wk)
            ratio_l2 = np.linalg.norm(c) / np.linalg.norm(h)
            ratio_w = np.sum(np.abs(c) * wk) / np.sum(np.abs(h) * wa)
            print(f'  L-raiz s={s[i].real:.9f}{s[i].imag:+.6f}i: ||C h||/||h|| l2={ratio_l2:.3e} pond={ratio_w:.3e};'
                  f' sharp mais proxima {sk[j].real:.9f}{sk[j].imag:+.1e}i (dist {abs(sk[j]-s[i]):.2e});'
                  f' |cos|(Ch, autovetor sharp) l2={cos2:.4f} pond={cw:.4f}; residuo rel. ||K(s)Ch||/||Ch||={rel:.3e}')
        pos = sorted([z for z in sk if z.real > 0], key=lambda z: (round(z.real, 3), z.imag))
        print('  raizes sharp com Re s>0 e |Im s|<=3/4:',
              ', '.join(f'{z.real:.6f}{z.imag:+.6f}i' for z in pos if abs(z.imag) <= 0.75))


def cmd_census(args):
    for tag in args.boxes:
        A, addr = read_rt_matrix(SPEC / f'rt-A-{tag}.dat', 'CHOPTUIK_RT_SPECTRAL_MATRIX_V1')
        s = -np.linalg.eigvals(A)
        pos = s[s.real > 0]
        stripA = [z for z in pos if -0.25 < z.imag <= 0.25]
        stripB = [z for z in pos if 0.25 < z.imag <= 0.75]
        print(f'\n== {tag} (dim {len(A)}): Re s>0 total {len(pos)}')
        print('  faixa A (-1/4,1/4]:', ', '.join(f'{z.real:.6f}{z.imag:+.4f}i' for z in sorted(stripA, key=lambda z: z.real)))
        print('  faixa B (1/4,3/4]: ', ', '.join(f'{z.real:.6f}{z.imag:+.4f}i' for z in sorted(stripB, key=lambda z: z.real)) or '(vazia)')
        near = sorted(pos, key=lambda z: abs(z.imag - 0.5))[:4]
        print('  4 raizes com Im mais proximo de 1/2:', ', '.join(f'{z.real:.6f}{z.imag:+.4f}i' for z in near))
        big = [z for z in pos if z.real > 1]
        print(f'  Re s>1: {len(big)} raizes; menores |Im-1/2| entre elas:',
              ', '.join(f'{z.real:.3f}{z.imag:+.3f}i' for z in sorted(big, key=lambda z: abs((z.imag % 1) - 0.5))[:3]))


def cmd_exact(args):
    mu, raw = parse_refa()
    re_, im_, two = raw[3]['exact'][1, 0]
    im = Fr(im_) * Fr(2) ** two
    print('mu_RefA =', mu, '=', float(mu))
    print('RefA4_(1,0) = (', Fr(re_) * Fr(2) ** two, ',', im, ')')
    eps = Fr(1, 2 ** 25) + Fr(1, 2 ** 277)
    weight = 2 * Fr(65, 64)
    lower = (abs(im) - eps / weight) / 2
    print('peso (m=1,n=0) =', weight, '; |b(dtau w)| >=', float(lower), '> 117/1000 ?', lower > Fr(117, 1000))


def M13_apply(mu, w, f1, f3):
    """Bloco (f1,f3) de M(s): independente de s (SHARP_LINEARIZED_IDENTITY, formulas de M).
    (Mf)1 = -mu(xi dxi+2)P f1 - G1(w,Rf) ; (Mf)2 = -mu(P f1 + 2 f3) - G2(w,Rf),
    Rf = (X P f1, 0, 0, X P f3, 0); G = (1/2){...} da identidade homogenea de RT."""
    Pf1, Pf3 = f1.P(), f3.P()
    zero = Fld.zeros(0, 0)
    e = [Pf1.times_xi(), zero, zero, Pf3.times_xi(), zero]
    G = sharp_G(w, e)
    r1 = -1 * mu * (Pf1.xi_dxi() + 2 * Pf1) - G[0]
    r2 = -1 * mu * (Pf1 + 2 * f3) - G[1]
    return r1, r2


def sharp_G(w, Om):
    w1, w2, w3, w4 = w
    P = lambda f: f.P()
    O1, O2, O2s, O3, O4 = Om
    ev = lambda f: 0.5 * (f + P(f))
    od = lambda f: 0.5 * (f - P(f))
    terms = [
        ([2 * P(w2) - P(w3) + w1 - 2 * w2 + w3, P(w2) + w1 - w2], ev(O1)),
        ([P(w2) - P(w3) + w2 - w3, P(w2) + w1 - w2], od(O1)),
        ([-1 * P(w3) - 3 * w1 + w3, -1 * P(w2) - w1 - 7 * w2 + 4 * w3], ev(O2)),
        ([-1 * P(w2) - P(w3) - w2 - w3, -1 * P(w2) - w1 + w2 - 4 * w3], od(O2)),
        ([P(w3) - w1 - w3, 3 * P(w2) - w1 + w2], ev(O2s)),
        ([P(w2) + P(w3) + w2 + w3, 3 * P(w2) - w1 + w2], od(O2s)),
        ([2 * w1, 4 * w2], ev(O3)),
        ([-2 * P(w2) - 2 * w2, -4 * w2], od(O3)),
        ([2 * P(w4) - 2 * w4, -4 * w4], ev(O4)),
        ([2 * P(w4) + 2 * w4, 4 * w4], od(O4)),
    ]
    out = [0, 0]
    for col, right in terms:
        if not np.any(right.a):
            continue
        for i in range(2):
            out[i] = out[i] + 0.5 * conv(col[i], right)
    return [o if isinstance(o, Fld) else Fld.zeros(0, 0) for o in out]


def cmd_m13(args):
    mu_fr, raw = parse_refa()
    for tag in args.boxes:
        _, _, caddr = read_constraints(SPEC / f'rt-A-{tag}.dat.constraints',
                                       int(open(SPEC / f'rt-A-{tag}.dat').readlines()[1].split()[1]))
        Mb, Nb = box_of(caddr)
        mu, w = refa_fields(mu_fr, raw, 2 * Mb, 2 * Nb)
        dim = len(caddr)
        Mat = np.zeros((dim, dim))
        leak = 0.0
        for c, (d, k, m, n) in enumerate(caddr):
            f = [Fld.zeros(Mb, Nb), Fld.zeros(Mb, Nb)]
            z = 1.0 if k == 0 else 1.0j
            f[d].a[Mb + m, n] = z
            if m:
                f[d].a[Mb - m, n] = np.conj(z)
            r1, r2 = M13_apply(mu, w, f[0], f[1])
            leak = max(leak, np.abs(0.5 * (r1 + r1.P()).a).max())
            Mat[:, c] = dof_vector([r1, r2], caddr)
        sv = np.linalg.svd(Mat, compute_uv=False)
        ev = np.linalg.eigvals(Mat)
        # principal (sem fundo) para comparacao
        print(f'{tag}: dim={dim}; sigma_min(M13_N)={sv[-1]:.4e}, sigma_max={sv[0]:.3e};'
              f' |autovalor| minimo={np.abs(ev).min():.4e}; parte par vazada na linha 1={leak:.1e}')


def M_apply(mu, w, F, shift):
    f1, f2, f3, f4 = F
    Rf = [f1.P().times_xi(), f2.P().times_xi(), Fld.zeros(0, 0), f3.P().times_xi(), f4.P().times_xi()]
    G = sharp_G(w, Rf)
    m1 = -1 * mu * (f1.P().xi_dxi() + 2 * f1.P()) - G[0]
    m2 = O_A(f2.P().times_xi(), mu, shift).P() - mu * (f1.P() + f2.P() - f2 + 2 * f3) - G[1]
    return [m1, m2]


def K_apply(mu, w, c, shift):
    g = gamma2_sharp(w, c)
    return [O_B(c[0], mu, shift) + g[0], O_A(c[1], mu, shift) + mu * c[0].div_xi() + g[1]]


def cmd_identity(args):
    """K(s)C(s)h - M(s)L(s)h + G(h,e(w)) = 0 com os operadores float (RefA, h suave aleatorio)."""
    mu_fr, raw = parse_refa()
    mu, w = refa_fields(mu_fr, raw, 12, 30)
    rng = np.random.default_rng(7)
    h = []
    for comp in range(4):
        f = Fld.zeros(6, 14)
        for m in range(0, 7):
            if (comp < 3 and m % 2) or (comp == 3 and m % 2 == 0):
                continue
            for n in range(15):
                if comp == 0 and n % 2 == 0:
                    continue
                z = (rng.normal() + 1j * rng.normal()) * 0.5 ** (m + n) if m else rng.normal() * 0.5 ** n
                f.a[6 + m, n] += z
                if m:
                    f.a[6 - m, n] += np.conj(z)
        h.append(f)
    CAPS[0], CAPS[1] = 200, 600
    for shift in (0.0, 0.3 + 0.2j):
        q = omega_plus_parts(mu, w, h, shift)
        F, C = selection_rt(q)
        kc = K_apply(mu, w, C, shift)
        ml = M_apply(mu, w, F, shift)
        e = [-1 * p.P() for p in omega_plus_parts(mu, w, w, 0.0, background_quadratic=False)]
        corr = sharp_G(h, e)
        res = [a - b + c for a, b, c in zip(kc, ml, corr)]
        scale = max(np.abs(x.a).max() for x in kc)
        print(f's={shift}: max|KC-ML+G(h,e)| = {max(np.abs(x.a).max() for x in res):.3e} (escala {scale:.3e});'
              f' sem correcao: {max(np.abs((a-b).a).max() for a, b in zip(kc, ml)):.3e}')


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest='cmd', required=True)
    for name in ('validate', 'gauge', 'sharp', 'census', 'm13'):
        p = sub.add_parser(name)
        p.add_argument('boxes', nargs='+')
        p.add_argument('--all-columns', action='store_true')
        p.add_argument('--decompose', action='store_true')
        p.add_argument('--targets', nargs='+', default=['0.7332', '0.401', '0.1683', '0'])
    sub.add_parser('exact')
    sub.add_parser('identity')
    args = ap.parse_args(argv)
    return dict(validate=cmd_validate, gauge=cmd_gauge, sharp=cmd_sharp,
                census=cmd_census, exact=cmd_exact, m13=cmd_m13, identity=cmd_identity)[args.cmd](args)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]) or 0)
