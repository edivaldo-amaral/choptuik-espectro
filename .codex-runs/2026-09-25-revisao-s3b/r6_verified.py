#!/usr/bin/env python3
"""R6 (revisao independente de S3b): (a) verified_inverse contra referencias de alta precisao
(Arb, 400 bits, enclosure rigoroso) em matrizes pequenas, incluindo mal condicionadas; (b) fundo_L: a cota de
Young de cada bloco B_d contra a norma do bloco montado (secao finita); (c) a constante que liga
o erro de RT (norma l1 de RT, soma dos 4 campos) a cota de Young do operador: max_j ||C_j||_2."""
import math
import sys
from pathlib import Path

import numpy as np
import flint

sys.path.insert(0, str(Path('scripts').resolve()))
from verified_inverse import cota_inv, cota_inv_X, cota_loewner, gram_linhas
from verified_norms import cota_norma

flint.ctx.prec = 400
rng = np.random.default_rng(12345)


def _acb(A):
    return flint.acb_mat([[flint.acb(complex(x).real, complex(x).imag) for x in l] for l in np.atleast_2d(A)])


class Ref(float):
    """valor de referencia (ponto medio) com o intervalo rigoroso [lo, hi] (Arb, 400 bits)."""
    pass


def _norma_acb(Y):
    G = Y.conjugate().transpose()*Y if Y.ncols() <= Y.nrows() else Y*Y.conjugate().transpose()
    ev = G.eig(algorithm='rump')
    his = [float(e.real.upper()) for e in ev]; los = [float(e.real.lower()) for e in ev]
    k = int(np.argmax(his))
    r = Ref(math.sqrt(max(float(ev[k].real.mid()), 0.0)))
    r.lo = math.sqrt(max(max(los), 0.0)); r.hi = math.sqrt(max(his))
    return r


def ref_norma(A):
    """||A||_2 com enclosure rigoroso (Arb)."""
    return _norma_acb(_acb(A))


def ref_inv_X(H, X):
    return _norma_acb(_acb(H).solve(_acb(X)))


def cmat(n, k):
    return rng.standard_normal((n, k)) + 1j*rng.standard_normal((n, k))


def mal_cond(n, cond):
    U, _ = np.linalg.qr(cmat(n, n)); V, _ = np.linalg.qr(cmat(n, n))
    s = np.logspace(0, -math.log10(cond), n)
    return (U*s) @ V.conj().T


def caso(nome, H, X=None, erro_X=0.0):
    try:
        if X is None:
            c = cota_inv(H)
            ref = ref_inv_X(H, np.eye(len(H)))
        else:
            nH = cota_inv(H) if erro_X else None
            c = cota_inv_X(H, X, erro_X=erro_X, nHinv=nH)
            ref = ref_inv_X(H, X)
        ok = c >= ref.hi
        if c < ref.lo:
            ok = False
        print(f'  {nome:48s} verdadeiro {float(ref):.12e}  cota {c:.12e}  razao {c/float(ref):.6f}  '
              + ('OK' if ok else '*** COTA ABAIXO DO VERDADEIRO ***'))
        return ok
    except RuntimeError as e:
        print(f'  {nome:48s} nao verificou (RuntimeError: {e}) -- aceitavel (sem cota falsa)')
        return True


def parte_a():
    print('(a) verified_inverse contra Arb (400 bits; OK = cota >= limite superior do enclosure)')
    ok = True
    n = 30
    H = np.eye(n) + 0.3*cmat(n, n)/math.sqrt(n)
    ok &= caso('aleatoria bem condicionada ||H^-1||', H)
    ok &= caso('aleatoria bem condicionada ||H^-1 X|| (30x20)', H, cmat(n, 20))
    ok &= caso('aleatoria ||H^-1 x|| (vetor)', H, cmat(n, 1))
    for cond in (1e4, 1e8, 1e11, 1e13, 1e15):
        H = mal_cond(n, cond)
        ok &= caso(f'mal condicionada cond {cond:.0e}: ||H^-1||', H)
        ok &= caso(f'mal condicionada cond {cond:.0e}: ||H^-1 X||', H, cmat(n, 7))
    # triangular superior do tipo J (crescimento acima da diagonal)
    ns = np.arange(40.0); mu = 0.168
    J = np.diag(mu*(ns + 1) + 0.4 + 3j) + np.triu(np.tile(2*mu*ns, (40, 1)), 1)
    ok &= caso('triangular tipo J_mu (40)', J)
    w = np.sqrt((5/4)**(2*ns) + (5/4)**(-2*ns)); w[0] = 1
    Jw = (w[:, None]*J)/w[None, :]
    ok &= caso('triangular tipo J_mu pesada (40)', Jw)
    # linhas com escalas muito diferentes
    D = np.diag(np.logspace(-6, 6, n))
    ok &= caso('escalas 1e-6..1e6 nas linhas', D @ (np.eye(n) + 0.2*cmat(n, n)/math.sqrt(n)))
    # X quase nulo e X com erro declarado
    H = np.eye(n) + 0.3*cmat(n, n)/math.sqrt(n)
    ok &= caso('X minusculo (1e-200)', H, 1e-200*cmat(n, 3))
    ok &= caso('com erro_X = 1e-3', H, cmat(n, 5), erro_X=1e-3)
    # bordejada (estrutura de Mh: [[H, h], [l, 0]])
    h = cmat(n, 1)[:, 0]; l = cmat(n, 1)[:, 0]
    Mb = np.zeros((n + 1, n + 1), complex); Mb[:n, :n] = H; Mb[:n, n] = h; Mb[n, :n] = l.conj()
    ok &= caso('bordejada [[H,h],[l*,0]]', Mb)
    # cota_norma
    print('    cota_norma:')
    for (m, k, c) in ((20, 30, 1.0), (30, 30, 1e-10), (25, 10, 1e12)):
        U_, _ = np.linalg.qr(cmat(m, m)); V_, _ = np.linalg.qr(cmat(k, k)); r_ = min(m, k)
        A = c*(U_[:, :r_]*np.logspace(0, -12, r_)) @ V_[:, :r_].conj().T
        cn = cota_norma(A); r = ref_norma(A)
        okk = cn >= r.hi
        ok &= okk
        print(f'  cota_norma {m}x{k} escala {c:.0e} (valores singulares 1..1e-12)  verdadeiro {float(r):.12e}  cota {cn:.12e}  '
              + ('OK' if okk else '*** ABAIXO ***'))
    # Loewner com S = projecao (como nVhi) e chute muito abaixo do verdadeiro
    H = np.eye(n) + 0.5*cmat(n, n)/math.sqrt(n)
    P, eP = gram_linhas(H)
    proj = np.zeros((n, n)); proj[20:, 20:] = np.eye(n - 20)
    r = ref_inv_X(H, proj[:, 20:])
    t = cota_loewner(P, eP, proj, 0.0, 0.95*float(r))
    print(f'  Loewner S = projecao, chute 0,95 x verdadeiro: verdadeiro {float(r):.12e}  cota {t:.12e}  ' + ('OK' if t >= r else '*** ABAIXO ***'))
    ok &= t >= r.hi
    print(f'  (a): {"todas as cotas acima do verdadeiro" if ok else "FALHA"}')
    return ok


def parte_b():
    import signed_operator_L as sl
    import fundo_L
    print('(b) fundo_L: cota de Young de B_d (x kappa1^|d|) contra a norma do bloco pesado montado')
    g = sl.campos()
    K1, K2 = sl.K1, sl.K2
    N = 260
    pior = 0.0
    for d in list(range(-40, 41)):
        y = fundo_L.norma_Bd(g, d)
        for par in (0, 1):
            ri = sl.Radial(par, N); ro = sl.Radial(par + d, N + 101)
            X = sl.bloco_B(g, ro, ri)
            wi = np.sqrt(np.concatenate([fundo_L.cf.omega(ns, K2) for _, ns in ri.comps]))
            wo = np.sqrt(np.concatenate([fundo_L.cf.omega(ns, K2) for _, ns in ro.comps]))
            v = np.linalg.norm((wo[:, None]*X)/wi[None, :], 2)
            pior = max(pior, v/y)
            if v > y:
                print(f'  *** d = {d}, paridade {par}: bloco {v:.6e} > Young {y:.6e}')
        if d in (-40, -17, -1, 0, 1, 2, 16, 17, 30, 40):
            print(f'  d = {d:+3d}: Young {y:.4e}; bloco montado (n < {N}) max paridade ~ {v:.4e}')
    print(f'  maior razao bloco/Young: {pior:.4f} (tem de ser <= 1)')
    db = fundo_L.dB_L(g, 16)
    tot = sum(K1**abs(d)*fundo_L.norma_Bd(g, d) for d in range(-40, 41))
    print(f'  dB_L(16) = {db};  sum_d k1^|d| Young(d) = {tot:.4f} (simbolo usado: 3,964043)')
    # resto |d| > 16: por d
    print('  Young(d) por |d| (maior paridade de sinal): ' + ', '.join(
        f'{d}:{max(fundo_L.norma_Bd(g, d), fundo_L.norma_Bd(g, -d)):.1e}' for d in (8, 12, 16, 17, 18, 20, 21, 24, 30, 40)))


def parte_c():
    import signed_operator_L as sl
    print('(c) erro de RT: cota de Young de dB por unidade da norma de RT (soma dos 4 campos)')
    # campo j com um coeficiente unitario em (m, n): norma RT = (2 - d_m0)(2 - d_n0) k1^m k2^n;
    # Young: sum_{d = +-m} k1^|d| ||M(d)||_2 com M = C_j (2 - d_n0) k2^n |v|  => razao = ||C_j||_2
    C = [np.zeros((4, 7)) for _ in range(4)]
    ENT = [(0, None), (1, 'e'), (1, 'o'), (2, 'e'), (2, 'o'), (3, 'e'), (3, 'o')]
    for k, (cin, par) in enumerate(ENT):
        chave = ('w1', None) if cin == 0 else (f'w{cin+1}', par)
        fator, saidas = sl.TABELA[chave]
        for cout in range(4):
            for coef, j, _ in saidas[cout]:
                C[j][cout, k] += 2*fator*abs(coef)
    ns = [float(np.linalg.norm(c, 2)) for c in C]
    print('  ||C_j||_2 = ' + ', '.join(f'{x:.4f}' for x in ns) + f';  max = {max(ns):.4f}  (usado: 40)')
    print('  Frobenius: ' + ', '.join(f'{float(np.linalg.norm(c)):.4f}' for c in C))
    print('  => ||dB||_Young <= max_j ||C_j||_2 * ||d omega||_RT (Young e subaditiva; o max sobre a bola l1 de RT e')
    print('     atingido num coeficiente so). RT: ||RefB omega|| <= 2^-25 e a correcao <= 2^-277 (soma dos campos).')


if __name__ == '__main__':
    parte_a()
    parte_b()
    parte_c()
