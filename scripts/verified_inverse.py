#!/usr/bin/env python3
"""Cotas verificadas de ||H^-1 X|| e ||H^-1|| SEM formar o inverso (condicao de Loewner).

Para H inversivel e X quaisquer (matrizes de ponto flutuante, EXATAS como dados),
    ||H^-1 X||_2 <= t  <=>  X X^* <= t^2 H H^*  <=>  t^2 H H^* - X X^* >= 0,
por congruencia com H^-1. O inverso usado no certificado e o EXATO da matriz calculada H
(um operador fixo), e so se precisa verificar que uma hermitiana e semidefinida positiva:
criterio de Rump (BIT 46, 2006, Teorema 2.3 e Corolario 2.4, via verified_norms._margem_rump)
na realificacao [[Re G, -Im G], [Im G, Re G]] (real simetrica, com os mesmos autovalores,
duplicados), montada NUM buffer so (n = 14 mil: 6,5 GB).

Erros: P = fl(H H^*) com ||P - H H^*|| <= gamma_k ||H||_F^2 (k = dimensao interna); o mesmo
para S = fl(X X^*); a combinacao t^2 P - S arredonda por entrada (<= u (t^2 |P| + |S|)); a
simetrizacao (G + G^*)/2 e exata em ponto flutuante; o que sobra entra na margem de Rump.
"""
from __future__ import annotations

import math
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from verified_norms import U, gamma, _margem_rump, _cholesky_no_lugar


def pd_hermitiana(G, extra: float) -> bool:
    """True so se G (hermitiana calculada; ||G_calc - G_exata|| <= extra) e exatamente
    definida positiva: G_exata >= (margem) I > 0. Um buffer real 2n x 2n; G nao e alterada.
    Se G vier numa lista [G], ela e retirada e liberada depois da realificacao (memoria)."""
    caixa = G if isinstance(G, list) else None
    if caixa is not None:
        G = caixa.pop()
    n = G.shape[0]
    fro = float(np.linalg.norm(G))*(1 + 1e-6)
    if np.iscomplexobj(G):
        Ar = np.empty((2*n, 2*n))
        np.add(G.real, G.real.T, out=Ar[:n, :n]); Ar[:n, :n] *= 0.5
        np.subtract(G.imag.T, G.imag, out=Ar[:n, n:]); Ar[:n, n:] *= 0.5     # -Im(herm)
        np.negative(Ar[:n, n:], out=Ar[n:, :n]); Ar[n:, n:] = Ar[:n, :n]
    else:
        Ar = np.empty((n, n)); np.add(G, G.T, out=Ar); Ar *= 0.5
    if caixa is not None:
        del G                                   # so o buffer real fica para a Cholesky
    d = np.diag(Ar).copy()
    if np.any(d <= 0):
        return False
    Ar[np.diag_indices(len(d))] -= _margem_rump(d, extra + U*fro)
    ok = _cholesky_no_lugar(Ar)
    del Ar
    return ok


def gram_linhas(A: np.ndarray):
    """P = A A^* e cota do erro (gamma_k ||A||_F^2, k = numero de colunas)."""
    return A @ A.conj().T, gamma(A.shape[1])*float(np.sum(np.abs(A)**2))


def cota_loewner(P, erro_P, S, erro_S, t0, passo=1.002, tentativas=60):
    """Menor t (a partir de t0, em passos geometricos) com t^2 P_exata - S_exata >= 0
    verificado. S = None significa S = I (exata)."""
    fP = float(np.linalg.norm(P))
    fS = float(np.linalg.norm(S)) if S is not None else math.sqrt(P.shape[0])
    t = t0
    for _ in range(tentativas):
        G = (t*t)*P
        if S is None:
            G[np.diag_indices(len(G))] -= 1.0
        else:
            G -= S
        extra = t*t*erro_P + (erro_S if S is not None else 0.0) + 2*U*(t*t*fP + fS)
        cx = [G]; del G
        ok = pd_hermitiana(cx, extra)
        if ok:
            return t
        t *= passo
    raise RuntimeError('cota_loewner nao verificou')


def cota_inv(H, chute=None, P=None, erro_P=None):
    """Cota verificada de ||H^-1|| (H exata como dado)."""
    if P is None:
        P, erro_P = gram_linhas(H)
    if chute is None:
        chute = 1/float(np.linalg.svd(H, compute_uv=False)[-1]) if len(H) <= 4000 else None
    return cota_loewner(P, erro_P, None, 0.0, chute*(1 + 1e-6))


def cota_inv_X(H, X, erro_X=0.0, nHinv=None, chute=None, P=None, erro_P=None):
    """Cota verificada de ||H^-1 X_exata||, com ||X - X_exata|| <= erro_X (entao soma
    ||H^-1|| erro_X; nHinv e obrigatoria se erro_X > 0)."""
    if P is None:
        P, erro_P = gram_linhas(H)
    S, erro_S = gram_linhas(X)
    if chute is None:
        Y = np.linalg.solve(H, X)
        from tile_certificate import norma2
        chute = norma2(Y); del Y
    t = cota_loewner(P, erro_P, S, erro_S, chute*(1 + 1e-6))
    if erro_X:
        assert nHinv is not None
        t += nHinv*erro_X
    return t


def _teste():
    rng = np.random.default_rng(3)
    for n, k in ((200, 300), (600, 400)):
        H = np.eye(n) + 0.3*(rng.standard_normal((n, n)) + 1j*rng.standard_normal((n, n)))/math.sqrt(n)
        X = (rng.standard_normal((n, k)) + 1j*rng.standard_normal((n, k)))/math.sqrt(n)
        Hi = np.linalg.inv(H)
        v1 = np.linalg.norm(Hi, 2); v2 = np.linalg.norm(Hi @ X, 2)
        c1 = cota_inv(H); c2 = cota_inv_X(H, X)
        print(f'n = {n}: ||H^-1|| {v1:.6f} <= {c1:.6f} ({c1/v1:.5f});  ||H^-1 X|| {v2:.6f} <= {c2:.6f} ({c2/v2:.5f})')
        assert v1 <= c1 <= v1*1.003 and v2 <= c2 <= v2*1.003


if __name__ == '__main__':
    _teste()
