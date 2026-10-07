#!/usr/bin/env python3
"""Cotas SUPERIORES verificadas de normas espectrais em ponto flutuante IEEE (binario64,
arredondamento ao mais proximo), sem aritmetica intervalar.

Positividade verificada (Rump, "Verification of positive definiteness", BIT 46:433-452,
2006, Teorema 2.3, cota I da secao 3 e Corolario 2.4). Para A REAL SIMETRICA OU COMPLEXA
HERMITIANA n x n (o Teorema 2.3 cobre as duas, com as mesmas constantes), se o Cholesky em
ponto flutuante (arredondamento ao mais proximo, qualquer ordem de avaliacao) de A~ termina,
onde a~_ij = a_ij fora da diagonal e a~_ii <= a_ii - c, com

    c >= ||Delta(A)||_2,   ||Delta(A)||_2 <= gamma_{n+1}/(1 - gamma_{n+1}) tr(A) + n M eta,
    M = 3 (2n + max a_ii),   gamma_k = k eps/(1 - k eps),   eps = 2^-53,  eta = 2^-1074,

entao A e definida positiva. Aqui a hermitiana e realificada ([[Re, -Im], [Im, Re]], real
simetrica 2n x 2n, definida positiva se e so se A e): e so mais conservador (n -> 2n).

Hipoteses que o codigo garante: (i) A~ EXATAMENTE simetrica -- o Gram calculado em ponto
flutuante nao sai exatamente hermitiano, entao usa-se a parte hermitiana (G + G*)/2, que em
ponto flutuante e exatamente simetrica, com o erro de arredondamento dela somado a margem;
a parte hermitiana de um erro nao aumenta a norma espectral dele; (ii) c e inflado por
1 + 1e-6, o que cobre o arredondamento da propria soma tr(A) (~1e4 termos); (iii) a~_ii
arredondado ao mais proximo nao passa de a_ii - c porque c ja inclui 2 eps max a_ii.

Norma de X: t e cota se t^2 I - X^*X e definida positiva. O produto X^*X feito em
ponto flutuante tem erro |fl(G) - G| <= gamma'_k |X|^* |X| (k = linhas de X), cuja
norma espectral e <= gamma'_k ||X||_F^2; entra somado a c.
"""
from __future__ import annotations

import math

import numpy as np

U = 2.0**-53
ETA = 2.0**-1074


def gamma(n: int) -> float:
    """Constante de erro de produtos COMPLEXOS de comprimento n: gamma_{4(n+1)}."""
    nn = 4*(n + 1)
    assert nn*U < 0.01
    return nn*U/(1 - nn*U)


def gamma_real(n: int) -> float:
    assert n*U < 0.01
    return n*U/(1 - n*U)


def _margem_rump(d: np.ndarray, extra: float) -> float:
    """c do Corolario 2.4 de Rump para a real simetrica com diagonal d >= 0 (dimensao len(d)):
    cota I (gamma_{n+1}/(1 - gamma_{n+1}) tr + n M eta), mais extra, mais 2 eps max a_ii
    (arredondamento de a_ii - c ao mais proximo); o fator 1 + 1e-6 cobre o arredondamento
    das somas que formam c."""
    n = len(d)
    g = gamma_real(n + 1)
    dmax = float(d.max())
    c = g/(1 - g)*float(np.sum(d)) + n*3*(2*n + dmax)*ETA
    return (c + extra + 2*U*dmax)*(1 + 1e-6)


def _cholesky_no_lugar(Ar: np.ndarray) -> bool:
    """Cholesky (LAPACK dpotrf) sobrescrevendo Ar; Ar e simetrica, entao Ar.T (contigua em
    Fortran) e a mesma matriz e nao ha copia."""
    from scipy.linalg import lapack
    _, info = lapack.dpotrf(Ar.T, lower=0, clean=0, overwrite_a=1)
    return info == 0


def definida_positiva(A: np.ndarray, extra: float = 0.0) -> bool:
    """True so se A (hermitiana) - extra I for definida positiva, verificado pelo
    criterio de Rump na realificacao (real simetrica 2n x 2n)."""
    if np.iscomplexobj(A):
        Ar = np.block([[A.real, -A.imag], [A.imag, A.real]])
    else:
        Ar = np.array(A, float)
    # parte simetrica exata (ver a docstring do modulo); o arredondamento entra na margem
    fro_A = float(np.linalg.norm(Ar))*(1 + 1e-6)
    Ar = (Ar + Ar.T)*0.5
    extra = extra + U*fro_A
    d = np.diag(Ar).copy()
    if np.any(d < 0):
        return False
    Ar[np.diag_indices(len(d))] -= _margem_rump(d, extra)
    return _cholesky_no_lugar(Ar)


def _potencia(X, iters=60, seed=0):
    rng = np.random.default_rng(seed)
    v = rng.standard_normal(X.shape[1]) + 1j*rng.standard_normal(X.shape[1])
    v /= np.linalg.norm(v)
    for _ in range(iters):
        u = X.conj().T @ (X @ v)
        nu = np.linalg.norm(u)
        if nu == 0:
            return 0.0
        v = u/nu
    return float(np.linalg.norm(X @ v))


def cota_gram(G: np.ndarray, erro: float = 0.0, chute: float | None = None) -> float:
    """Cota superior verificada de sqrt(lambda_max(G_exata)), G hermitiana >= 0 calculada,
    com ||G_calculada - G_exata||_2 <= erro."""
    if G.size == 0:
        return math.sqrt(max(erro, 0.0))
    n = G.shape[0]
    if chute is None:
        chute = math.sqrt(max(float(np.linalg.eigvalsh(G)[-1]), 0.0)) if n <= 3000 else None
    if chute is None:
        # potencia sobre G (G = X^*X nao disponivel): maior autovalor de G
        rng = np.random.default_rng(0)
        v = rng.standard_normal(n) + 0j; v /= np.linalg.norm(v)
        lam = 0.0
        for _ in range(80):
            w = G @ v; lam = float(np.linalg.norm(w))
            if lam == 0:
                break
            v = w/lam
        chute = math.sqrt(lam)
    t = max(chute*(1 + 1e-6), math.sqrt(2*erro), 1e-150)
    escala = max(float(np.abs(np.diag(G)).max()), 1e-300)
    fro_G = float(np.linalg.norm(G))*(1 + 1e-6)
    complexa = np.iscomplexobj(G)
    m = 2*n if complexa else n
    Ar = np.empty((m, m))                     # um buffer so, reaproveitado entre tentativas
    for _ in range(200):
        # realificacao de t^2 I - G: [[t^2 - Re G, Im G], [-Im G, t^2 - Re G]]
        # parte hermitiana (G + G*)/2, exatamente simetrica em ponto flutuante (a + b = b + a,
        # a - b = -(b - a)); cada entrada arredonda uma vez: erro <= eps |G| entrada a entrada
        if complexa:
            np.add(G.real, G.real.T, out=Ar[:n, :n]); Ar[:n, :n] *= -0.5
            np.subtract(G.imag, G.imag.T, out=Ar[:n, n:]); Ar[:n, n:] *= 0.5
            np.negative(Ar[:n, n:], out=Ar[n:, :n]); Ar[n:, n:] = Ar[:n, :n]
        else:
            np.add(G, G.T, out=Ar); Ar *= -0.5
        Ar[np.diag_indices(m)] += t*t
        # t^2 - g_ii arredonda: erro <= eps (t^2 + |g_ii|); a simetrizacao: <= eps ||G||_F
        extra = erro + U*(t*t + escala) + U*fro_G
        d = np.diag(Ar).copy()
        if np.all(d >= 0):
            Ar[np.diag_indices(m)] -= _margem_rump(d, extra)
            if _cholesky_no_lugar(Ar):
                return t
        t *= 1.001
    raise RuntimeError('cota_gram nao convergiu')


def cota_norma(X: np.ndarray, erro_X: float = 0.0) -> float:
    """Cota superior verificada de ||X_exata||_2, com ||X_calculada - X_exata||_2 <= erro_X."""
    if X.size == 0:
        return erro_X
    m, n = X.shape
    Y = X if n <= m else X.conj().T              # Gram do lado menor
    k = Y.shape[0]
    G = Y.conj().T @ Y
    erro_G = gamma(k)*float(np.sum(np.abs(Y)**2))
    return cota_gram(G, erro_G, chute=_potencia(Y)) + erro_X


def _teste():
    rng = np.random.default_rng(1)
    for (m, n) in ((50, 30), (300, 300), (1200, 800)):
        X = rng.standard_normal((m, n)) + 1j*rng.standard_normal((m, n))
        verdade = np.linalg.norm(X, 2)
        c = cota_norma(X)
        print(f'{m}x{n}: verdadeira {verdade:.10f}  cota {c:.10f}  razao {c/verdade:.8f}')
        assert verdade <= c <= verdade*1.01
    # matriz singular e nula
    X = np.zeros((10, 5)); X[0, 0] = 3.0
    print('posto 1:', cota_norma(X))
    print('nula   :', cota_norma(np.zeros((4, 4))))
    # Gram acumulado com erro declarado
    Y = rng.standard_normal((400, 200))
    print('gram   :', cota_gram(Y.T @ Y, erro=1e-10), np.linalg.norm(Y, 2))


if __name__ == '__main__':
    _teste()
