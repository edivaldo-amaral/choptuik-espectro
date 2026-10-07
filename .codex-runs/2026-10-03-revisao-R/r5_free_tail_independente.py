#!/usr/bin/env python3
"""Revisao R5: recalculo INDEPENDENTE da cauda livre ||Q (I - G)|| para Q = J_mu^-1 (componentes 2-4) e
K_mu^-1 (componente 1), na norma l1 de Chebyshev com peso radial k2' generico, e comparacao com
sharp_exterior.free_tail(M, N, k2').

Independente do formalismo de bandas de RT e de ALGEBRAIC_TAIL: os operadores sao montados na base T_n com
numpy.polynomial.chebyshev (chebder/chebmulx), J = i m/2 + mu((1 + xi) d_xi + 1), K = i m/2 + mu(xi d_xi + 1)
(RT: J_mu = (D^+ + mu), K_mu = (J + PJP)/2), e cada coluna da inversa e obtida por substituicao regressiva
(os operadores sao triangulares superiores em n, com diagonal mu(n + 1) + i m/2).
Norma: ||f|| = |c_0| + sum_{n>=1} k^n |c_n|  (f = sum c_n T_n), que e a norma simetrica de RT
sum_{n in Z} k^|n| |v_n|. A norma de operador l1 ponderada e o maximo das somas ponderadas de coluna.
Fator real: cada coeficiente complexo q vira o bloco real [[Re q, -Im q], [Im q, Re q]], de norma de coluna
|Re q| + |Im q|; reportamos esse valor exato e tambem (3/2)|q|, que e o que free_tail majora.
A cauda e amostrada: |m| >= M (n de 0 a NMAX) e n >= N (todos os m da amostra). O supremo na cauda infinita e
coberto pela prova em papel (relatorio); a amostra so testa se free_tail e de fato um majorante e quao justo.
"""
import sys
from fractions import Fraction as Fr
from pathlib import Path
import numpy as np
import numpy.polynomial.chebyshev as C

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from sharp_exterior import free_tail  # noqa: E402

MU = 0.16830707896344996951


def operador(nmax, tipo):
    """Matriz real (nmax+1)^2 de (1+xi)d_xi (tipo 'J') ou xi d_xi (tipo 'K') na base T_0..T_nmax."""
    A = np.zeros((nmax + 1, nmax + 1))
    for n in range(nmax + 1):
        e = np.zeros(n + 1); e[n] = 1.0
        d = C.chebder(e) if n > 0 else np.zeros(1)
        xd = C.chebmulx(d) if n > 0 else np.zeros(1)
        col = xd.copy()
        if tipo == 'J':
            col = C.chebadd(col, d)
        col = np.pad(col, (0, nmax + 1 - len(col)))[:nmax + 1]
        A[:, n] = col
    return A


def inversa(A, m, mu):
    """Inversa (triangular superior) de i m/2 + mu (A + 1), por solve_triangular (substituicao regressiva)."""
    from scipy.linalg import solve_triangular
    N1 = A.shape[0]
    Mz = mu*(A + np.eye(N1)) + 1j*(m/2)*np.eye(N1)
    return solve_triangular(Mz, np.eye(N1, dtype=complex), lower=False)


def somas_colunas(Qinv, k2):
    """Somas ponderadas de cada coluna: (real-bloco, complexa)."""
    N1 = Qinv.shape[0]
    w = np.array([1.0] + [k2**i for i in range(1, N1)])
    re = (w[:, None]*(np.abs(Qinv.real) + np.abs(Qinv.imag))).sum(axis=0)/w
    cp = (w[:, None]*np.abs(Qinv)).sum(axis=0)/w
    return re, cp


def main():
    casos = [(8, 24), (16, 64), (32, 128)]
    pesos = [Fr(9, 8), Fr(17, 16), Fr(33, 32), Fr(5, 4)]
    NMAX = 700
    ops = {t: operador(NMAX, t) for t in ('J', 'K')}
    # conferencia da diagonal: mu(n+1) + i m/2
    assert np.allclose(np.diag(ops['J']), np.arange(NMAX + 1))
    assert np.allclose(np.diag(ops['K']), np.arange(NMAX + 1))
    ms = sorted(set(list(range(0, 41)) + [48, 64, 96, 128, 200, 400]))
    falhou = False
    inv = {(t, m): inversa(ops[t], m, MU) for t in ('J', 'K') for m in ms}
    # conferencia: residuo da inversa
    for (t, m), Qi in list(inv.items())[:3]:
        Mz = MU*(ops[t] + np.eye(NMAX + 1)) + 1j*(m/2)*np.eye(NMAX + 1)
        assert np.abs(Mz @ Qi - np.eye(NMAX + 1)).max() < 1e-9
    for k2 in pesos:
        kf = float(k2)
        somas = {key: somas_colunas(Qi, kf) for key, Qi in inv.items()}
        for (M, N) in casos:
            sup_real, sup_c, arg = 0.0, 0.0, None
            for (t, m), (re, cp) in somas.items():
                for n in range(NMAX + 1):
                    if abs(m) >= M or n >= N:
                        if re[n] > sup_real:
                            sup_real, arg = float(re[n]), (t, m, n)
                        sup_c = max(sup_c, 1.5*float(cp[n]))
            ft = float(free_tail(M, N, radial=k2, mu_min=Fr(1, 6), mu_max=Fr(17, 100)))
            ok = sup_real <= ft and sup_c <= ft
            falhou |= not ok
            print(f"k2'={str(k2):6s} G={M}x{N}: sup amostrado real {sup_real:.5f} (em {arg}), "
                  f"(3/2)|q| {sup_c:.5f}; free_tail {ft:.5f}  {'OK' if ok else 'VIOLA'}")
    print('free_tail majora a amostra em todos os casos' if not falhou else 'HA VIOLACAO')
    return 1 if falhou else 0


if __name__ == '__main__':
    sys.exit(main())
