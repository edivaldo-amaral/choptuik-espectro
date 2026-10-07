#!/usr/bin/env python3
"""R2 (revisao S4): a identidade de gauge L_w(mu)(Z + 1)w = (Z + 2)F_w, testada nos operadores do
codigo (J_livre, bloco_B de signed_operator_L), com montagem propria de F_w e de Z.

F_w (selecionado, quadratico): F(w) = J(0)w + mu S Gamma1 w + S Xi Gamma2(w, w). Com
B_w = mu S Gamma1 + 2 S Xi Gamma2(w, .) (bloco_B com o fundo w) vale, SE Gamma2 e simetrica,
S Xi Gamma2(w, w) = (B_w - B_0) w / 2, logo F(w) = J(0)w + (B_w + B_0) w / 2. Testes:
 (T0) simetria de Gamma2: (B_v - B_0) w == (B_w - B_0) v, para v, w aleatorios;
 (T1) identidade para w polinomial aleatorio (caixa exata), com a convencao d_tau -> +i m/2 e,
      como controle, com -i m/2 e com (Z + 1) no lado direito;
 (T2) RefA: F(RefA) (defeito, deve ser pequeno se F nao tem termo constante e e o de RT) e
      L_A(mu_A)(Z_A + 1)RefA contra (Z_A + 2)F(RefA), na saida completa (sem truncar);
 (T3) ||L_A(mu_A) g_A|| / ||g_A|| com g_A de gauge_vector_L (o vetor do certificado) e com o meu.
 (T4) (Z + 1)w nao nulo e g_A contra o g proprio (mesmo vetor?)."""
from __future__ import annotations
import math, sys, time
from pathlib import Path
import numpy as np

RAIZ = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(RAIZ/'scripts'))
import signed_operator_L as sl

MU = sl.MU
T0 = time.time()
out = []
def P(s=''):
    print(s, flush=True); out.append(s)


def xidxi(a):
    ns = np.arange(len(a)); b = ns*a.astype(complex)
    for par in (0, 1):
        idx = ns[ns % 2 == par]
        v = (idx*a[idx])[::-1]
        cs = np.concatenate([[0], np.cumsum(v)[:-1]])[::-1]
        b[idx] += 2*cs
    return b


def campos_para_modos(g, M, N):
    """vetor por modo (layout Radial(m, N)), m com sinal, a partir de g[j][m >= 0]."""
    v = {}
    for m in sl.modos(M):
        r = sl.Radial(m, N); x = np.zeros(r.dim, complex)
        for c, ns in r.comps:
            a = sl.modo(g, c, m); full = np.zeros(max(N, len(a)), complex); full[:len(a)] = a
            o0, o1 = r.off[c]; x[o0:o1] = full[ns]
        v[m] = x
    return v


def aplica_Z(v, N, mu, k, sinal=1):
    """(Z + k) por modo: Z = mu^-1 d_tau + xi d_xi, d_tau -> sinal * i m/2."""
    w = {}
    for m, x in v.items():
        r = sl.Radial(m, N); y = np.zeros_like(x)
        for c, ns in r.comps:
            o0, o1 = r.off[c]
            full = np.zeros(N, complex); full[ns] = x[o0:o1]
            z = (k + sinal*1j*m/(2*mu))*full + xidxi(full)
            y[o0:o1] = z[ns]
        w[m] = y
    return w


def aplica_L(gbg, v, Nin, Nout, Mout, s, mu, livre=True):
    """(J(s) + B_gbg) v (B com fundo gbg; gbg=None: so mu S Gamma1), saida em |m| < Mout, n < Nout."""
    zero = [dict() for _ in range(4)]
    g = zero if gbg is None else gbg
    res = {}
    for mo in sl.modos(Mout):
        ro = sl.Radial(mo, Nout); acc = np.zeros(ro.dim, complex)
        for mi, x in v.items():
            if abs(mo - mi) <= sl.BAND_M and np.any(x):
                acc += sl.bloco_B(g, ro, sl.Radial(mi, Nin), mu) @ x
        if livre and mo in v:
            ri = sl.Radial(mo, Nin)
            pos = {}; o = 0
            for c, ns in ro.comps:
                for kk, n in enumerate(ns):
                    pos[(c, int(n))] = o + kk
                o += len(ns)
            idx = np.array([pos[(c, int(n))] for c, ns in ri.comps for n in ns])
            acc[idx] += sl.J_livre(ri, s, mu) @ v[mo]
        res[mo] = acc
    return res


def comb(*termos):
    ks = set().union(*[t[1].keys() for t in termos])
    return {m: sum(c*t[m] for c, t in termos if m in t) for m in ks}


def norma(v, pesos=False):
    import nk_L as nk
    t = 0.0
    for m, x in v.items():
        N = None
        t += float(np.sum(np.abs(x)**2)) if not pesos else 0.0
    return math.sqrt(t)


def aleatorio(M0, N0, rng, esc=0.3):
    g = [dict() for _ in range(4)]
    for j in range(4):
        for m in range(M0 + 1):
            if (j < 3 and m % 2) or (j == 3 and m % 2 == 0):
                continue
            a = np.zeros(101, complex)
            for n in range(N0 + 1):
                if j == 0 and n % 2 == 0:
                    continue
                a[n] = esc*(rng.standard_normal() + (1j*rng.standard_normal() if m > 0 else 0))/(1 + n + m)
            g[j][m] = a
    return g


def main():
    rng = np.random.default_rng(7)
    mu = 0.37                                     # mu generico: a identidade vale para todo mu
    M0, N0 = 4, 9
    Mo, No = 2*M0 + 2, 2*N0 + 6                   # caixa de saida exata
    gv, gw = aleatorio(M0, N0, rng), aleatorio(M0, N0, rng)
    v = campos_para_modos(gv, M0 + 1, N0 + 1); w = campos_para_modos(gw, M0 + 1, N0 + 1)
    # T0 simetria
    Bv_w = aplica_L(gv, w, N0 + 1, No, Mo, 0.0, mu, livre=False); B0_w = aplica_L(None, w, N0 + 1, No, Mo, 0.0, mu, livre=False)
    Bw_v = aplica_L(gw, v, N0 + 1, No, Mo, 0.0, mu, livre=False); B0_v = aplica_L(None, v, N0 + 1, No, Mo, 0.0, mu, livre=False)
    d1 = comb((1, Bv_w), (-1, B0_w)); d2 = comb((1, Bw_v), (-1, B0_v))
    num = math.sqrt(sum(float(np.sum(np.abs(d1[m] - d2[m])**2)) for m in d1)); den = math.sqrt(sum(float(np.sum(np.abs(d1[m])**2)) for m in d1))
    P(f'(T0) simetria de S Xi Gamma2: ||G(v,w) - G(w,v)|| / ||G(v,w)|| = {num/den:.2e}')
    # T1 identidade
    def F(gbg, x, Nin):
        Lw0 = aplica_L(gbg, x, Nin, No, Mo, 0.0, mu)             # J(0) + B_w
        B0 = aplica_L(None, x, Nin, No, Mo, 0.0, mu, livre=False)
        return comb((1, Lw0), (0.5, aplica_L(gbg, x, Nin, No, Mo, 0.0, mu, livre=False)), (-0.5, aplica_L(gbg, x, Nin, No, Mo, 0.0, mu, livre=False)),
                    (-0.5, aplica_L(gbg, x, Nin, No, Mo, 0.0, mu, livre=False)), (0.5, B0))
    for sinal, kF, nome in ((1, 2, 'convencao do autor: d_tau -> +i m/2, (Z + 2)F'),
                            (-1, 2, 'controle: d_tau -> -i m/2'),
                            (1, 1, 'controle: (Z + 1)F no lado direito')):
        Zw = aplica_Z(w, N0 + 1, mu, 1, sinal)
        lhs = aplica_L(gw, Zw, N0 + 1, No, Mo, mu, mu)
        Fw = F(gw, w, N0 + 1)
        rhs = aplica_Z(Fw, No, mu, kF, sinal)
        e = math.sqrt(sum(float(np.sum(np.abs(lhs[m] - rhs.get(m, 0))**2)) for m in lhs))
        n_ = math.sqrt(sum(float(np.sum(np.abs(lhs[m])**2)) for m in lhs))
        P(f'(T1) w aleatorio ({M0}x{N0}), mu = {mu}: {nome}: ||L_w(mu)(Z+1)w - (Z+k)F_w|| / ||L_w(mu)(Z+1)w|| = {e/n_:.2e}')
    # T2 RefA
    g = sl.campos()
    Ng = 101; Mg = 41
    wA = campos_para_modos(g, Mg, Ng)
    MoA, NoA = Mg + 41, Ng + 102
    No_s, Mo_s = NoA, MoA
    def FA(x):
        L0 = aplica_L(g, x, Ng, NoA, MoA, 0.0, MU)
        B0 = aplica_L(None, x, Ng, NoA, MoA, 0.0, MU, livre=False)
        BA = aplica_L(g, x, Ng, NoA, MoA, 0.0, MU, livre=False)
        return comb((1, L0), (-0.5, BA), (0.5, B0))            # J(0)w + (B_w + B_0) w/2
    FwA = FA(wA)
    nw = math.sqrt(sum(float(np.sum(np.abs(x)**2)) for x in wA.values()))
    nF = math.sqrt(sum(float(np.sum(np.abs(x)**2)) for x in FwA.values()))
    P(f'(T2) defeito do codigo em RefA: ||F(RefA)|| = {nF:.3e} (sem pesos; ||RefA|| = {nw:.3e})   [{time.time()-T0:.0f}s]')
    ZwA = aplica_Z(wA, Ng, MU, 1)
    lhs = aplica_L(g, ZwA, Ng, NoA, MoA, MU, MU)
    rhs = aplica_Z(FwA, NoA, MU, 2)
    nl = math.sqrt(sum(float(np.sum(np.abs(x)**2)) for x in lhs.values()))
    nr = math.sqrt(sum(float(np.sum(np.abs(x)**2)) for x in rhs.values()))
    e = math.sqrt(sum(float(np.sum(np.abs(lhs[m] - rhs.get(m, 0))**2)) for m in lhs))
    ng = math.sqrt(sum(float(np.sum(np.abs(x)**2)) for x in ZwA.values()))
    P(f'(T2) ||L_A(mu_A) g_A|| = {nl:.4e}, ||(Z_A + 2)F(RefA)|| = {nr:.4e}, diferenca {e:.2e}; ||g_A|| = {ng:.4e};'
      f' ||L_A g_A||/||g_A|| = {nl/ng:.2e}  (sem pesos, saida completa)   [{time.time()-T0:.0f}s]')
    # T3/T4: g_A do autor (gauge_vector_L.g_modo) contra o meu
    import gauge_vector_L as G
    dif = 0.0; nn2 = 0.0
    for m in sl.modos(Mg):
        ga, _ = G.g_modo(g, m, Ng)
        dif += float(np.sum(np.abs(ga - ZwA[m])**2)); nn2 += float(np.sum(np.abs(ga)**2))
    P(f'(T4) g_A (gauge_vector_L) contra (Z_A + 1)RefA proprio: diferenca relativa {math.sqrt(dif/nn2):.1e}')
    P(f'  [{time.time()-T0:.0f}s]')
    Path(__file__).with_name('r2_saida.txt').write_text('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
