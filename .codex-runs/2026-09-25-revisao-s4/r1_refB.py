#!/usr/bin/env python3
"""R1 (revisao S4): RefB exato com leitor PROPRIO, ||B_L(RefB)|| por Young com tabela PROPRIA
(transcrita do artigo de RT, eq. skdfhkfdjhdkhjdgfdjh33hk), ||RefB||_L, ||(Z_A + 1) RefB||_L,
conversao do Corr (bola de RT, norma l1 com pesos kappa^|m|, kappa^|n|) e DELTA_MU contra MU_80.

Nao usa gauge_refB.ler, fundo_L.matriz_M nem signed_operator_L.TABELA para as contas; so compara
a tabela propria com a do autor no fim."""
from __future__ import annotations
from fractions import Fraction as Fr
import math, re, sys, time, subprocess
from pathlib import Path
import numpy as np

RAIZ = Path(__file__).resolve().parents[2]
A_PATH = RAIZ/'.cache/rt-1203.3766v1/sourcecode/RefA.dat'
AB_PATH = RAIZ/'.cache/rt-1203.3766v1/RefAplusB.dat'
K1, K2 = Fr(65, 64), Fr(5, 4)
k1, k2 = 65/64, 5/4
T0 = time.time()

GRUPO = re.compile(r'\(((?:\(-?\d+,-?\d+\),?)+)\)')
PAR = re.compile(r'\((-?\d+),(-?\d+)\)')


def ler(path):
    """Leitor proprio: cabecalho do mu, num_comp; por Field: TwoExp, off_m, off_n, num_m, num_n e a
    lista aninhada (((par, ...), (par, ...)), ...) com GRUPOS explicitos. Devolve (mu exato, campos)
    com campos[f] = dict(exp, nm, nn, grupos) e grupos[i] = lista de (re, im) inteiros."""
    txt = path.read_text()
    cab = re.match(r'mu_MultiField\s+DyadicQ\s+TwoExp\s+(-?\d+)\s+coeff\s+(-?\d+)\s+MultiField\s+num_comp\s+(\d+)\s', txt)
    assert cab, 'cabecalho inesperado'
    mu = Fr(int(cab[2]))*Fr(2)**int(cab[1]); ncomp = int(cab[3])
    partes = txt.split('\nField\n')[1:]
    assert len(partes) == ncomp, (len(partes), ncomp)
    campos = []
    for p in partes:
        h = {}
        for k in ('TwoExp', 'off_m', 'off_n', 'num_m', 'num_n'):
            h[k] = int(re.search(rf'(?m)^{k}\s+(-?\d+)\s*$', p)[1])
        corpo = p[p.index('((('):].strip()
        assert corpo.startswith('(((') and corpo.endswith(')))')
        interno = corpo[1:-1]                        # tira o parentese externo
        grupos = []
        pos = 0
        for g in GRUPO.finditer(interno):
            assert interno[pos:g.start()] in ('', ','), 'lixo entre grupos'
            pos = g.end()
            pares = [(int(a), int(b)) for a, b in PAR.findall(g[1])]
            assert len(pares) == h['num_n'], (len(pares), h)
            grupos.append(pares)
        assert interno[pos:].strip() == '', 'sobra no fim'
        assert len(grupos) == h['num_m'], (len(grupos), h)
        assert h['off_m'] == 0 and h['off_n'] == 0
        campos.append(dict(exp=h['TwoExp'], nm=h['num_m'], nn=h['num_n'], grupos=grupos))
    return mu, campos


def main():
    out = []
    def P(s=''):
        print(s, flush=True); out.append(s)
    muA, A = ler(A_PATH); P(f'RefA lido: mu = {muA} = {float(muA):.17f}; campos: ' +
                           ', '.join(f"(2^{c['exp']}, {c['nm']}x{c['nn']})" for c in A) + f'  [{time.time()-T0:.0f}s]')
    muAB, AB = ler(AB_PATH); P(f'RefA+RefB lido: mu (TwoExp {math.log2(muAB.denominator):.0f}) = {float(muAB):.20f}; campos: ' +
                              ', '.join(f"(2^{c['exp']}, {c['nm']}x{c['nn']})" for c in AB) + f'  [{time.time()-T0:.0f}s]')
    # --- RefB EXATO: diferenca de diadicos em inteiros (unidade 2^eAB) ---
    B = []            # B[f] = matriz complexa (nm x nn) de floats (arredondados de valores exatos)
    maxrel = 0.0
    paridade = []
    for f, (ca, cb) in enumerate(zip(A, AB)):
        eA, eB = ca['exp'], cb['exp']; assert eB <= eA
        esc = 1 << (eA - eB)
        X = np.zeros((cb['nm'], cb['nn']), complex)
        nz_par = {'m_par': 0, 'm_impar': 0, 'n_par': 0, 'n_impar': 0}
        for m in range(cb['nm']):
            gb = cb['grupos'][m]
            ga = ca['grupos'][m] if m < ca['nm'] else None
            for n in range(cb['nn']):
                xr, xi = gb[n]
                if ga is not None and n < ca['nn']:
                    ar, ai = ga[n]; xr -= ar*esc; xi -= ai*esc
                if xr or xi:
                    X[m, n] = complex(math.ldexp(float(xr), eB) if xr else 0.0, math.ldexp(float(xi), eB) if xi else 0.0)
                    nz_par['m_par' if m % 2 == 0 else 'm_impar'] += 1
                    nz_par['n_par' if n % 2 == 0 else 'n_impar'] += 1
        B.append(X); paridade.append(nz_par)
        # conferencia de exatidao: um coeficiente grande refeito em Fraction
        m0, n0 = np.unravel_index(np.argmax(np.abs(X)), X.shape)
        ex = Fr(cb['grupos'][m0][n0][0])*Fr(2)**eB - (Fr(ca['grupos'][m0][n0][0])*Fr(2)**eA if m0 < ca['nm'] and n0 < ca['nn'] else 0)
        maxrel = max(maxrel, abs(float(ex) - X[m0, n0].real)/max(abs(float(ex)), 1e-300))
    P(f'RefB = (RefA+RefB) - RefA calculado em inteiros; conferencia Fraction do maior coeficiente: erro relativo {maxrel:.1e}')
    P('paridade dos coeficientes NAO nulos de RefB (conta por campo): ' + '; '.join(
        f'v{f+1}: {d}' for f, d in enumerate(paridade)))
    # layout: v1 so n impar e m par; v2, v3 m par; v4 m impar (setor A)
    P('paridade em RefA (nao nulos):')
    for f, c in enumerate(A):
        mm = [m for m in range(c['nm']) for n in range(c['nn']) if c['grupos'][m][n] != (0, 0)]
        nn = [n for m in range(c['nm']) for n in range(c['nn']) if c['grupos'][m][n] != (0, 0)]
        P(f'  v{f+1}: m pares {sum(1 for m in mm if m % 2 == 0)}, m impares {sum(1 for m in mm if m % 2)},'
          f' n pares {sum(1 for n in nn if n % 2 == 0)}, n impares {sum(1 for n in nn if n % 2)}')
    # --- normas de L (m com sinal: para m > 0 entra kappa1^2m + kappa1^-2m; omega2(0) = 1) ---
    def wn(nn):
        ns = np.arange(nn); return np.where(ns == 0, 1.0, k2**(2.0*ns) + k2**(-2.0*ns))
    def wm(m):
        return 1.0 if m == 0 else k1**(2.0*m) + k1**(-2.0*m)
    nB2 = sum(wm(m)*float(np.sum(wn(X.shape[1])*np.abs(X[m])**2)) for X in B for m in range(X.shape[0]))
    P(f'||RefB||_L = {math.sqrt(nB2):.6e}   (autor: 2,3506e-10)')
    # (Z_A + 1) = 1 + i m/(2 mu_A) + xi d_xi; xi d_xi em coeficientes simetricos (derivacao propria:
    # f = a_0 + 2 sum a_n T_n, xi T_n' = n T_n + 2n sum_{0 <= k < n, k = n mod 2} T_k (T_0 com peso 1/2)
    # => b_k = k a_k + 2 sum_{n > k, n = k mod 2} n a_n, k >= 0)
    muA_f = float(muA)
    def xidxi(a):
        nn = len(a); ns = np.arange(nn); b = ns*a
        for par in (0, 1):
            idx = ns[ns % 2 == par]
            v = (idx*a[idx])[::-1]
            cs = np.concatenate([[0], np.cumsum(v)[:-1]])[::-1]     # soma dos n > k de mesma paridade
            b[idx] += 2*cs
        return b
    # teste da formula de xi d_xi contra numpy.polynomial.chebyshev
    from numpy.polynomial import chebyshev as C
    rng = np.random.default_rng(0); a = rng.standard_normal(9)
    c = a.copy(); c[1:] *= 2                                   # coeficientes usuais de Chebyshev
    der = C.chebmulx(C.chebder(c)) if hasattr(C, 'chebmulx') else None
    x = np.linspace(-1, 1, 7)
    lhs = x*C.chebval(x, C.chebder(c)); bb = xidxi(a); cb_ = bb.copy(); cb_[1:] *= 2
    P(f'teste xi d_xi (coef. simetricos) contra chebder: max erro {np.max(np.abs(lhs - C.chebval(x, cb_))):.1e}')
    nZB2 = 0.0
    for X in B:
        for m in range(X.shape[0]):
            z = (1 + 1j*m/(2*muA_f))*X[m] + xidxi(X[m])
            nZB2 += wm(m)*float(np.sum(wn(X.shape[1])*np.abs(z)**2))
    P(f'||(Z_A + 1) RefB||_L = {math.sqrt(nZB2):.6e}   (autor: 2,4542e-08)')
    # --- tabela PROPRIA (transcrita do artigo, S Xi Gamma2(v, w)): entrada -> (fator, saidas[4]) ---
    # termos (coef, campo 0..3, refletido P)
    T = {
        ('w1', None): (Fr(1, 4), [[(1, 1, 1), (1, 1, 0), (-1, 2, 1), (-1, 2, 0)], [(1, 0, 0), (-1, 1, 0), (1, 1, 1)],
                                  [(-1, 0, 0), (1, 1, 0), (-1, 1, 1)], []]),
        ('w2', 'e'): (Fr(1, 2), [[(1, 0, 0), (1, 2, 1), (-1, 2, 0)], [(2, 1, 0), (2, 1, 1)], [(-1, 1, 0), (-1, 1, 1)],
                                 [(1, 3, 1), (1, 3, 0)]]),
        ('w2', 'o'): (Fr(1, 2), [[], [(-1, 0, 0), (1, 1, 1), (-1, 1, 0)], [(1, 0, 0)], [(1, 3, 1), (-1, 3, 0)]]),
        ('w3', 'e'): (Fr(1, 2), [[(-1, 0, 0)], [], [], []]),
        ('w3', 'o'): (Fr(1, 2), [[(-1, 1, 0), (-1, 1, 1)], [], [], []]),
        ('w4', 'e'): (Fr(1, 2), [[(1, 3, 0), (-1, 3, 1)], [], [(1, 3, 1), (1, 3, 0)], [(1, 1, 1), (1, 1, 0)]]),
        ('w4', 'o'): (Fr(1, 2), [[(1, 3, 0), (1, 3, 1)], [], [(1, 3, 1), (-1, 3, 0)], [(1, 1, 1), (-1, 1, 0)]]),
    }
    sys.path.insert(0, str(RAIZ/'scripts'))
    import signed_operator_L as sl
    iguais = all(float(T[k][0]) == sl.TABELA[k][0] and
                 [sorted(x) for x in T[k][1]] == [sorted(tuple(t) for t in x) for x in sl.TABELA[k][1]] for k in T)
    P(f'tabela propria (do artigo) == signed_operator_L.TABELA: {iguais}')
    cols = list(T)
    # Young: ||B_L(RefB)|| <= sum_d k1^|d| ||M(d)||_2, M(d)[cout, col] = 2 fator l1(sum |c| |v_j,d|)
    nn = B[0].shape[1]; ns = np.arange(nn); pl1 = np.where(ns == 0, 1.0, 2*k2**ns)
    absB = [np.abs(X) for X in B]
    Mmax = max(X.shape[0] for X in B) - 1
    tot = 0.0; por_d = []
    for d in range(0, Mmax + 1):
        M = np.zeros((4, len(cols)))
        for k, col in enumerate(cols):
            fator, saidas = T[col]
            for cout in range(4):
                if saidas[cout]:
                    v = sum(abs(c)*absB[j][d] for c, j, p in saidas[cout])
                    M[cout, k] = 2*float(fator)*float(np.sum(v*pl1))
        s = float(np.linalg.svd(M, compute_uv=False)[0])
        por_d.append(s)
        tot += (1 if d == 0 else 2)*k1**d*s                   # d e -d tem o mesmo |M|
    young = tot*(1 + 1e-9)
    P(f'||B_L(RefB)|| <= {young:.6e}  (Young, |d| <= {Mmax}; autor: 4,6748e-08 com fator 1,001)')
    dom = sorted(range(len(por_d)), key=lambda d: -k1**d*por_d[d])[:5]
    P('  maiores termos d: ' + ', '.join(f'd={d}: {k1**d*por_d[d]:.2e}' for d in dom))
    # quanto de RefB esta fora da caixa de RefA (m > 40 ou n > 100)
    fora = sum(wm(m)*float(np.sum(wn(nn)[ (np.arange(nn) > 100) | (m > 40)]*np.abs(X[m][(np.arange(nn) > 100) | (m > 40)])**2))
               for X in B for m in range(X.shape[0]))
    P(f'  parte de ||RefB||_L fora da caixa 41x101 de RefA: {math.sqrt(fora):.2e}')
    # --- Corr: ||B_L(Corr)|| <= max_j W_j ||Corr||_SSPACE, W_j = soma de 2 fator |c| nos termos do campo j ---
    W = [0.0]*4
    for col in cols:
        fator, saidas = T[col]
        for cout in range(4):
            for c, j, p in saidas[cout]:
                W[j] += 2*float(fator)*abs(c)
    R = 2.0**-277
    P(f'Corr: pesos W_j = {W}; ||B_L(Corr)|| <= max W_j * 2^-277 = {max(W)*R:.2e}  (< 1e-60: {max(W)*R < 1e-60})')
    # --- mu ---
    sys.path.insert(0, str(RAIZ/'scripts'))
    from tile_perron import MU_80, DELTA_MU
    d80 = abs(MU_80 - muAB)
    P(f'|MU_80 - mu_(A+B)| = {float(d80):.3e};  |mu* - MU_80| <= {float(d80 + Fr(1, 2**277)):.3e}  (<= 1e-80: {d80 + Fr(1, 2**277) <= Fr(1, 10**80)})')
    dA = abs(muA - muAB)
    P(f'|mu_A - mu_(A+B)| = {float(dA):.6e};  DELTA_MU = {DELTA_MU:.6e};  |mu_A - mu*| <= {float(dA + Fr(1, 2**277)):.6e} <= DELTA_MU: {float(dA + Fr(1, 2**277)) <= DELTA_MU}')
    P(f'  [{time.time()-T0:.0f}s]')
    Path(__file__).with_name('r1_saida.txt').write_text('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
