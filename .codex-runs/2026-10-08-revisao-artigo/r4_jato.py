#!/usr/bin/env python3
"""R4.2 (revisao de 08/10/2026): teste numerico do passo R7 do Lema 4.9 e da Prop. 4.12 com uma familia afim de
matrizes. Modelo de dimensao finita:
    L(s) = A + s,  C(s) = C0 + s C1,  K(s) = K0 + s,  M(s) = M0 + s M1,
com K(s) C(s) = M(s) L(s) para todo s. A "evolucao" e u' + A u = 0 (perfil e^{s tau} h <-> L(s) h = 0) e a
"constraint" e c(tau) = C1 u'(tau) + C0 u(tau) (a derivada em tau entra so com o coeficiente constante C1).

Confere:
 (a) a identidade KC = ML equivale a M1 = C1, M0 = C0 + K0 C1 - C1 A, K0 C0 = M0 A, e implica K0 D = D A com
     D = C0 - C1 A (como na Sec. 4.5);
 (b) u(tau) = e^{s0 tau} sum_j tau^j/j! h_{k-1-j} resolve u' + A u = 0 se e so se (h_0..h_{k-1}) e cadeia de Jordan
     L(s0) h_i = -h_{i-1};
 (c) c(tau) = e^{s0 tau} sum_j tau^j/j! [C(s0) h_{k-1-j} + C1 h_{k-2-j}]  (formula de jato);
 (d) c == 0  <=>  restricoes de jato  <=>  cadeia em ker(D|E);
 (e) dimensao fisica (por amostragem de c em varios tau) = dim E - rank(D|E);
 (f) exemplo 2x2: impor C(s0) h_j = 0 para cada j da a dimensao errada (nos dois sentidos).
"""
import numpy as np
from math import factorial
from scipy.linalg import expm, null_space

rng = np.random.default_rng(20261008)
TOL = 1e-9
saida = []


def p(*x):
    s = ' '.join(str(y) for y in x)
    print(s); saida.append(s)


def familia(A, D, C1, K0):
    """C0, M0, M1 a partir de A, D, C1, K0 (exige K0 D = D A)."""
    C0 = D + C1 @ A
    M1 = C1.copy()
    M0 = C0 + K0 @ C1 - C1 @ A
    return C0, M0, M1


def cadeia(A, s0, v, k):
    """h_{k-1} = v, h_{i-1} = -(A + s0) h_i."""
    n = len(A); N = A + s0*np.eye(n)
    h = [None]*k; h[k - 1] = v
    for i in range(k - 1, 0, -1):
        h[i - 1] = -N @ h[i]
    return h


def u_de(h, s0, tau):
    k = len(h)
    return np.exp(s0*tau)*sum(tau**j/factorial(j)*h[k - 1 - j] for j in range(k))


def du_de(h, s0, tau):
    k = len(h)
    v = s0*u_de(h, s0, tau)
    v = v + np.exp(s0*tau)*sum(tau**(j - 1)/factorial(j - 1)*h[k - 1 - j] for j in range(1, k))
    return v


# ------------------------------------------------------------------ caso 1: bloco de Jordan 3 + outros autovalores
s0 = 0.37 + 0.11j
n = 6
J = np.zeros((n, n), complex)
J[:3, :3] = -s0*np.eye(3) + np.diag([1, 1], 1)        # bloco de Jordan em -s0 (tamanho 3)
J[3, 3] = -1.3; J[4, 4] = 0.2 + 0.5j; J[5, 5] = -0.7
P = rng.normal(size=(n, n)) + 1j*rng.normal(size=(n, n))
A = P @ J @ np.linalg.inv(P)
Pinv = np.linalg.inv(P)
# D = polinomio em A projetado: D = W (A + s0)^1 restrito ao bloco + algo nos outros (K0 = A comuta)
for nome, D in (('D = (A+s0) (posto 1 em E)', A + s0*np.eye(n)),
                ('D = (A+s0)^2 (posto... ker maior)', np.linalg.matrix_power(A + s0*np.eye(n), 2)),
                ('D = 0 (tudo fisico)', np.zeros((n, n), complex)),
                ('D = I (nada fisico)', np.eye(n, dtype=complex))):
    K0 = A.copy()
    C1 = rng.normal(size=(n, n)) + 1j*rng.normal(size=(n, n))
    C0, M0, M1 = familia(A, D, C1, K0)
    # (a)
    err_a = max(np.abs((K0 + s*np.eye(n)) @ (C0 + s*C1) - (M0 + s*M1) @ (A + s*np.eye(n))).max()
                for s in (0, 1.0, -0.3 + 2j, 5j))
    err_KD = np.abs(K0 @ D - D @ A).max()
    E = P[:, :3]                                       # espaco de raizes em s0 (colunas do bloco)
    # (e) dimensao fisica por amostragem: v em E com D e^{-A tau} v = 0 para varios tau
    taus = np.linspace(-1, 1, 7)
    G = np.vstack([D @ expm(-A*t) @ E for t in taus])
    dim_fis = E.shape[1] - np.linalg.matrix_rank(G, tol=1e-8)
    dim_form = E.shape[1] - np.linalg.matrix_rank(D @ E, tol=1e-8)
    p(f'[{nome}] KC-ML: {err_a:.1e}; K0 D - D A: {err_KD:.1e}; dim fisica (amostragem) {dim_fis}; '
      f'formula dim E - posto(D|E) = {dim_form}')
    assert err_a < TOL and err_KD < TOL and dim_fis == dim_form
    # (b), (c), (d) para cadeias de comprimento k = 1, 2, 3 geradas por v aleatorio em E
    for k in (1, 2, 3):
        v = E @ (rng.normal(size=3) + 1j*rng.normal(size=3))
        # projeta v para que a cadeia tenha comprimento exatamente k: v em ker (A + s0)^k
        Nk = np.linalg.matrix_power(A + s0*np.eye(n), k)
        base = null_space(np.vstack([Nk, (np.eye(n) - E @ np.linalg.pinv(E))]))
        v = base @ (rng.normal(size=base.shape[1]) + 1j*rng.normal(size=base.shape[1]))
        h = cadeia(A, s0, v, k)
        res_cad = max(np.abs((A + s0*np.eye(n)) @ h[i] + (h[i - 1] if i > 0 else 0)).max() for i in range(k))
        res_ev = max(np.abs(du_de(h, s0, t) + A @ u_de(h, s0, t)).max() for t in (-0.7, 0.0, 0.4, 1.3))
        Cs = lambda x: (C0 + s0*C1) @ x
        jato = [Cs(h[0])] + [Cs(h[i]) + C1 @ h[i - 1] for i in range(1, k)]
        err_c = max(np.abs(C1 @ du_de(h, s0, t) + C0 @ u_de(h, s0, t)
                           - np.exp(s0*t)*sum(t**j/factorial(j)*jato[k - 1 - j] for j in range(k))).max()
                    for t in (-0.7, 0.0, 0.4, 1.3))
        Dh = [D @ x for x in h]
        err_Dj = max(np.abs(Dh[i] - jato[i]).max() for i in range(k))
        p(f'   k={k}: cadeia {res_cad:.1e}; u\'+Au {res_ev:.1e}; c(tau) - formula de jato {err_c:.1e}; '
          f'D h_i - jato_i {err_Dj:.1e}')
        assert res_cad < TOL and res_ev < TOL and err_c < 1e-8 and err_Dj < TOL

# ------------------------------------------------------------------ caso 2: contra-exemplos 2x2
s0 = 0.5
A = np.array([[-s0, 1.0], [0.0, -s0]])                # L(s0) = A + s0 = [[0,1],[0,0]]
h0 = np.array([1.0, 0.0]); h1 = np.array([0.0, -1.0])  # L(s0) h0 = 0, L(s0) h1 = -h0
assert np.allclose((A + s0*np.eye(2)) @ h1, -h0)
K0 = np.array([[-s0]])
for nome, D, C1 in (('D = e2^T, C1 h0 = D h1', np.array([[0.0, 1.0]]), np.array([[-1.0, 0.0]])),
                    ('D = 0, C1 h0 != 0', np.zeros((1, 2)), np.array([[1.0, 0.0]]))):
    C0, M0, M1 = familia(A, D, C1, K0)
    err_a = max(np.abs((K0 + s) @ (C0 + s*C1) - (M0 + s*M1) @ (A + s*np.eye(2))).max() for s in (0, 1.0, 2j))
    Cs = C0 + s0*C1
    taus = np.linspace(-1, 1, 7)
    G = np.vstack([D @ expm(-A*t) for t in taus])
    dim_fis = 2 - np.linalg.matrix_rank(G, tol=1e-10)
    # contagem "ingenua": vetores da cadeia com C(s0) h_j = 0 para cada j
    ingenua = int(abs((Cs @ h0)[0]) < 1e-12) + int(abs((Cs @ h0)[0]) < 1e-12 and abs((Cs @ h1)[0]) < 1e-12)
    jato = int(abs((Cs @ h0)[0]) < 1e-12) + int(abs((Cs @ h0)[0]) < 1e-12 and abs((Cs @ h1 + C1 @ h0)[0]) < 1e-12)
    p(f'[2x2, {nome}] KC-ML {err_a:.1e}; dim fisica {dim_fis}; contagem com jato {jato}; '
      f'contagem ingenua (C(s0)h_j = 0) {ingenua}')
    assert err_a < TOL and dim_fis == jato and ingenua != jato

p('TUDO CONFERE: formula de R7, restricoes de jato e Prop. 4.12 no modelo; a contagem ingenua erra nos 2x2.')
