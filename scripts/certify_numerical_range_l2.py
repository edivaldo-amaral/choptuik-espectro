#!/usr/bin/env python3
"""Borda direita do campo de valores de L (F3) e de K (constraints), num UNICO
produto interno, com a estrutura de blocos de l^2 feita como l^2.

Substitui a cadeia de certify_f3_chain.py, que misturava convencoes
(ALTERNATIVE_ROUTES 5.25).

O espaco. Os coeficientes de RT sao os SIMETRICOS (indice n em Z, f_-n = f_n;
por isso o operador livre acopla 2*mu*n uniformemente, inclusive a linha n=0).
Usa-se l^2(W) simetrico, W(k) = kappa1^|k1| kappa2^|k2|, que em coordenadas reais
e D = diag(sqrt(c_m c_n) kappa1^m kappa2^n), c = 2-[.=0]. Nele:
  * l^1(w) de RT esta contido em l^2(W) com norma <= 1;
  * multiplicacao por f e convolucao em Z^2 com peso submultiplicativo, e Young
    da ||M_f|| <= norma l^1 de f -- exatamente a field_norm auditada;
  * a divisao regularizada por xi tem norma <= 2k^-1/(1-k^-2) (somas de linha e
    de coluna iguais);
  * pedacos (componente, paridade) sao ortogonais, entao a norma de um operador
    em blocos e limitada pela norma ESPECTRAL da matriz de normas de blocos --
    nao pelo maximo, que e a estrutura de l^1.

A parte livre. Herm(D J D^-1) >= lambda0 equivale, pelo pencil, a
Herm(Wt J) - lambda0 Wt definida positiva, Wt = D^2 = diag(c_n kappa2^(2n)),
racional. A matriz e semisseparavel de posto 1 mais diagonal, com recursao
escalar exata; a cauda fecha pelo argumento uniforme:
  gamma_{k+1} = rho_k [gamma_k + q_k(1-gamma_k)^2/(q_k(1-gamma_k)+1)] < rho_k
sempre que gamma_k < 1, e rho_k < 1 para k >= N0. Logo gamma_k < 1 e todo pivo
p_k = e_k (q_k (1-gamma_k) + 1) e positivo. Nada de intervalo amostrado.

A escala em mu. A parte radial livre e mu vezes uma matriz fixa, e a derivada
temporal e antissimetrica (some no Hermitiano). Entao lambda_min escala
linearmente com mu; certifica-se mu=1 e aplica-se mu em [1/6, 17/100] (RT).
"""
from __future__ import annotations

import argparse
import json
from fractions import Fraction as Q
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from background_couplings import RT_BACKGROUND_ERROR
from component_exterior import read_fields, field_expression, field_norm
from sharp_propagation import table

REF = Path('.cache/rt-1203.3766v1/sourcecode/RefA.dat')
MU_MIN, MU_MAX = Q(1, 6), Q(17, 100)


# ---------------------------------------------------------------- parte livre

ESPACO = 'simetrico'      # ou 'toro': L^2 do toro de raios (kappa1, kappa2)


def peso(n: int, kappa2: Q) -> Q:
    """D_n^2 do espaco escolhido (o fator de Fourier e constante por bloco e cancela)."""
    if ESPACO == 'toro':
        return Q(1) if n == 0 else kappa2**(2*n) + kappa2**(-2*n)
    return (1 if n == 0 else 2)*kappa2**(2*n)


def modelo(passo: int, kappa2: Q, N: int):
    """Indices n do modelo radial: passo 1 -> todos; passo 2 -> impares (v1, c1)."""
    ns = list(range(N)) if passo == 1 else list(range(1, 2*N, 2))
    Wt = [peso(n, kappa2) for n in ns]
    return ns, Wt


def recursao(passo: int, kappa2: Q, lam0: Q, N0: int):
    """Recursao escalar exata (mu=1) ate N0; devolve pivos, gamma_N0 e rho_N0."""
    ns, Wt = modelo(passo, kappa2, N0+1)
    alpha, pivos = Q(0), []
    for k in range(N0):
        u, v, e = Wt[k], Q(ns[k]), Wt[k]*(1-lam0)
        til = u - alpha*v
        p = til*v + e
        pivos.append(p)
        if p <= 0:
            return pivos, None, None
        alpha += til*til/p
    u, v = Wt[N0], Q(ns[N0])
    gamma = alpha*v/u
    n0, n1 = ns[N0], ns[N0]+passo
    rho = Q(n1, n0)*Wt[N0]/peso(n1, kappa2)            # decrescente em n nos dois espacos

    return pivos, gamma, rho


def ldlt_generico(passo: int, kappa2: Q, lam0: Q, N: int):
    """LDL^T exato de Herm(Wt J) - lam0 Wt (mu=1): conferencia da recursao."""
    ns, Wt = modelo(passo, kappa2, N)
    J = [[Q(0)]*N for _ in range(N)]
    for b, n in enumerate(ns):
        J[b][b] = Q(n+1)
        for a in range(b):
            J[a][b] = Q(2*n)
    M = [[(Wt[i]*J[i][j] + Wt[j]*J[j][i])/2 - (lam0*Wt[i] if i == j else 0)
          for j in range(N)] for i in range(N)]
    piv = []
    for k in range(N):
        p = M[k][k]
        piv.append(p)
        for i in range(k+1, N):
            f = M[i][k]/p
            for j in range(k, N):
                M[i][j] -= f*M[k][j]
    return piv


def lam_min_float(passo: int, kappa2: float, N: int = 60):
    ns, Wt = modelo(passo, Q(kappa2).limit_denominator(10**6), N)
    ns = np.array(ns, float); d = np.sqrt(np.array([float(x) for x in Wt]))
    J = np.diag(ns+1)
    for b in range(N):
        J[:b, b] = 2*ns[b]
    H = np.diag(d) @ J @ np.diag(1/d)
    return float(np.min(np.linalg.eigvalsh((H+H.T)/2)))


def certifica_livre(passo: int, kappa2: Q, N0: int = 40):
    lmin = lam_min_float(passo, float(kappa2))
    lam0 = Q(lmin - 1e-4).limit_denominator(10**5)
    while lam0 > Q(lmin) - Q(1, 10**5):
        lam0 -= Q(1, 10**5)
    pivos, gamma, rho = recursao(passo, kappa2, lam0, N0)
    ok = gamma is not None and gamma < 1 and rho < 1 and all(p > 0 for p in pivos)
    # R_free(mu) = -mu*lam0 no pior mu
    R = -(MU_MAX if lam0 < 0 else MU_MIN)*lam0
    return dict(passo=passo, kappa2=str(kappa2), lambda0_mu1=str(lam0),
                lambda_min_float_mu1=lmin, N0=N0, gamma_N0=float(gamma) if gamma else None,
                rho_N0=float(rho) if rho else None, menor_pivo=float(min(pivos)),
                certificado=ok, R_livre=str(R), R_livre_decimal=float(R))


# ------------------------------------------------------------- fundo e blocos

def raiz_racional_acima(x: Q) -> Q:
    """Racional r >= sqrt(x), verificado."""
    r = Q(np.sqrt(float(x))).limit_denominator(10**9) + Q(1, 10**9)
    while r*r < x:
        r += Q(1, 10**9)
    return r


NORMA = None      # substituida por --fino


def norma_campo(c, k1, k2, paridade=None):
    if NORMA is None:
        return field_norm(c, paridade, fourier=k1, radial=k2)
    return NORMA(c, k1, k2, paridade)


def majorante_L(campos, k2=Q(5, 4)):
    """4 saidas x 7 pedacos (v1; v2p,v2i; v3p,v3i; v4p,v4i), normas l^2(W)."""
    rows = table('GAMMA2_macro.c', '_G2PAIR_')
    sel = (0, 1, 3, 4)
    pedacos = [(0, 0, 1, Q(1, 2))] + [(i, 2*i-1, 0, Q(1)) for i in (1, 2, 3)] \
        + [(i, 2*i, 1, Q(1)) for i in (1, 2, 3)]
    ordem = [(0, 1), (1, 0), (1, 1), (2, 0), (2, 1), (3, 0), (3, 1)]   # (comp, paridade)
    col = {(i, p): j for j, (i, p) in enumerate(ordem)}
    M = [[Q(0)]*7 for _ in range(4)]
    for out, row in enumerate(sel):
        for comp, k, par, fac in pedacos:
            c = field_expression(rows.get((row, k), ''), campos)
            M[out][col[comp, par]] += fac*norma_campo(c, Q(65, 64), k2, 1-par if out == 0 else None)
    r = MU_MAX*2*(1/k2)/(1-1/k2**2)                  # mu * ||R_x|| = mu*40/9
    M[1][col[0, 1]] += r; M[1][col[1, 1]] += 2*r
    M[2][col[0, 1]] += r; M[3][col[3, 1]] += 2*r
    return M, ordem


def majorante_K(campos):
    """2 saidas x 3 pedacos (c1; c2 coluna 4; c2 coluna 5), pesos (129/128, 9/8)."""
    rows = table('GAMMA2_SHARP_macro.c', '_G2SHARPPAIR_')
    nrm = lambda e: norma_campo(field_expression(e, campos), Q(129, 128), Q(9, 8))
    M = [[nrm(rows[i, 1])/2, nrm(rows[i, 4])/2, nrm(rows[i, 5])/2] for i in range(2)]
    k = Q(9, 8)
    M[1][0] += MU_MAX*2*(1/k)/(1-1/k**2)             # mu * 144/17
    return M, [(0, 1), (1, None), (1, None)]


def nivel_componente(M, ordem, ncomp):
    """||B_ik|| <= sqrt(sum_p ||B_i,(k,p)||^2): pedacos ortogonais."""
    C = [[Q(0)]*ncomp for _ in range(len(M))]
    for i, linha in enumerate(M):
        for k in range(ncomp):
            s = sum((linha[j]**2 for j, (c, _) in enumerate(ordem) if c == k), Q(0))
            C[i][k] = raiz_racional_acima(s) if s else Q(0)
    return C


def definida_positiva(A):
    n = len(A); A = [r[:] for r in A]
    for k in range(n):
        if A[k][k] <= 0:
            return False
        for i in range(k+1, n):
            f = A[i][k]/A[k][k]
            for j in range(k, n):
                A[i][j] -= f*A[k][j]
    return True


def borda(r, C, correcao, eta=None):
    """Menor t racional com t I - (diag(r) + Sym(eta C eta^-1)) definida positiva."""
    n = len(r)
    eta = eta or [Q(1)]*n
    S = [[(eta[i]*C[i][j]/eta[j] + eta[j]*C[j][i]/eta[i])/2 + (r[i] if i == j else 0)
          for j in range(n)] for i in range(n)]
    t0 = float(np.max(np.linalg.eigvalsh(np.array([[float(x) for x in row] for row in S]))))
    t = Q(t0).limit_denominator(10**6) + Q(1, 10**6)
    while not definida_positiva([[(t if i == j else 0) - S[i][j] for j in range(n)] for i in range(n)]):
        t += Q(1, 10**6)
    return t + correcao, t0


def melhores_pesos(r, C):
    """Proposta em float (Nelder-Mead simples); a borda e sempre revalidada."""
    n = len(r)
    Cf = np.array([[float(x) for x in row] for row in C]); rf = np.array([float(x) for x in r])
    def f(x):
        e = np.exp(np.concatenate([[0.0], x]))
        S = (e[:, None]*Cf/e[None, :]); S = (S+S.T)/2 + np.diag(rf)
        return np.max(np.linalg.eigvalsh(S))
    x = np.zeros(n-1); passo = 0.5; melhor = f(x)
    for _ in range(400):
        melhorou = False
        for i in range(n-1):
            for s in (passo, -passo):
                y = x.copy(); y[i] += s
                v = f(y)
                if v < melhor:
                    melhor, x, melhorou = v, y, True
        if not melhorou:
            passo /= 2
            if passo < 1e-4:
                break
    e = np.exp(np.concatenate([[0.0], x]))
    return [Q(float(v)).limit_denominator(1000) for v in e]


# ---------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', type=Path)
    ap.add_argument('--fino', action='store_true',
                    help='sup na elipse (Arb) em vez de Young para o fundo')
    ap.add_argument('--so', choices=('L', 'K'))
    ap.add_argument('--toro', action='store_true', help='espaco L^2 do toro (para o limitante pelo simbolo)')
    ap.add_argument('--so-livre', action='store_true')
    a = ap.parse_args()
    global ESPACO
    if a.toro:
        ESPACO = 'toro'
    if a.fino:
        global NORMA
        from radial_multiplier_sup import norma_multiplicacao
        NORMA = lambda c, k1, k2, par: norma_multiplicacao(c, k1, k2, par)

    # 0. a recursao escalar reproduz o LDL^T generico
    for passo, k2 in ((1, Q(5, 4)), (2, Q(5, 4)), (1, Q(9, 8)), (2, Q(9, 8))):
        for lam in (Q(-1, 2), Q(-1, 10)):
            pr, _, _ = recursao(passo, k2, lam, 10)
            pg = ldlt_generico(passo, k2, lam, 10)
            assert pr == pg[:len(pr)], (passo, k2, lam)       # identicos ate o 1o pivo <= 0
            assert len(pr) == 10 or pr[-1] <= 0
    print('recursao escalar == LDL^T generico, exato (4 modelos x 2 lambda0, 10 pivos)\n')

    campos = read_fields(REF)
    rel = {}
    for nome, k2, modelos, maj, ncomp, corr in (
            ('L', Q(5, 4), {0: 2, 1: 1, 2: 1, 3: 1}, majorante_L, 4, 36*RT_BACKGROUND_ERROR),
            ('K', Q(9, 8), {0: 2, 1: 1}, majorante_K, 2, 20*RT_BACKGROUND_ERROR)):
        if a.so and nome != a.so:
            continue
        livres = {p: certifica_livre(p, k2) for p in set(modelos.values())}
        for p, d in livres.items():
            print(f"{nome} livre, passo {p}: lambda0(mu=1)={float(Q(d['lambda0_mu1'])):+.5f}  "
                  f"gamma_N0={d['gamma_N0']:.4f}  rho_N0={d['rho_N0']:.4f}  "
                  f"certificado={d['certificado']}  R_livre<={d['R_livre_decimal']:+.5f}")
        assert all(d['certificado'] for d in livres.values())
        if a.so_livre:
            rel[nome] = dict(livres=livres, espaco=ESPACO)
            continue
        r = [Q(livres[modelos[i]]['R_livre']) for i in range(ncomp)]
        M, ordem = maj(campos)
        C = nivel_componente(M, ordem, ncomp)
        R1, t1 = borda(r, C, corr)
        eta = melhores_pesos(r, C)
        R2, t2 = borda(r, C, corr, eta)
        R = min(R1, R2)
        print(f"{nome} fundo, nivel componente (decimal):")
        for linha in C:
            print('   ', '  '.join(f'{float(x):7.4f}' for x in linha))
        print(f"{nome}: R <= {float(R1):.4f} (sem pesos)   R <= {float(R2):.4f} (pesos {[float(e) for e in eta]})")
        print(f"   => R_{nome} <= {float(R):.4f}\n")
        rel[nome] = dict(livres=livres, majorante_componentes=[[str(x) for x in l] for l in C],
                         pesos=[str(e) for e in eta], R=str(R), R_decimal=float(R),
                         correcao_fundo=str(corr))
    if a.out:
        a.out.write_text(json.dumps(rel, indent=2) + '\n')
        print('escrito:', a.out)


if __name__ == '__main__':
    main()
