#!/usr/bin/env python3
"""R6 (revisao T2): teste EXATO (racionais) da algebra de CONSTRAINT_ROOT_QUOTIENT.md com matrizes
aleatorias que satisfazem K C = M L, L(s) = A + s, K(s) = B + s, C = C0 + s C1, M = M0 + s M1.

Construcao geral: toda solucao de KC = ML (polinomial em s) tem M1 = C1, M0 = C0 + B C1 - C1 A e
B C0 = M0 A; com D := C0 - C1 A isso equivale a B D = D A. Logo sorteamos A, B, um entrelacador D
(nucleo de D -> B D - D A) e C1 livre, e pomos C0 = D + C1 A, M1 = C1, M0 = D + B C1.

Conferencias (por sorteio):
 1. KC = ML coeficiente a coeficiente (s^0, s^1, s^2);
 2. cadeias A h_j = -s0 h_j - h_{j-1}: D h_j = C(s0) h_j + C1 h_{j-1};
 3. dim do espaco FISICO, definido de forma independente como as solucoes u(t) = exp(-A t) u0, u0 em E,
    de L(d/dt) u = 0 com C(d/dt) u = 0 para todo t (isto e, D A^j u0 = 0 para todo j), igual a
    dim E - rank(D|E);
 4. se B + s0 e injetivo, D|E = 0 (toda a raiz satisfaz as constraints);
 5. o atalho errado (exigir C(s0) h = 0 em cada vetor de uma base de Jordan) pode dar outra dimensao.
"""
from fractions import Fraction as Fr
import json
import random
import sys


def zeros(r, c):
    return [[Fr(0)]*c for _ in range(r)]


def eye(n):
    M = zeros(n, n)
    for i in range(n):
        M[i][i] = Fr(1)
    return M


def mul(A, B):
    return [[sum((A[i][k]*B[k][j] for k in range(len(B))), Fr(0)) for j in range(len(B[0]))] for i in range(len(A))]


def add(A, B, a=1, b=1):
    return [[a*x + b*y for x, y in zip(ra, rb)] for ra, rb in zip(A, B)]


def rref(M):
    M = [row[:] for row in M]
    r = 0; piv = []
    rows, cols = len(M), len(M[0]) if M else 0
    for c in range(cols):
        p = next((i for i in range(r, rows) if M[i][c] != 0), None)
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        inv = 1/M[r][c]
        M[r] = [x*inv for x in M[r]]
        for i in range(rows):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [x - f*y for x, y in zip(M[i], M[r])]
        piv.append(c); r += 1
        if r == rows:
            break
    return M, piv


def rank(M):
    if not M or not M[0]:
        return 0
    return len(rref(M)[1])


def nullspace(M, ncols):
    if not M:
        return [[Fr(int(i == j)) for i in range(ncols)] for j in range(ncols)]
    R, piv = rref(M)
    livres = [c for c in range(ncols) if c not in piv]
    base = []
    for f in livres:
        v = [Fr(0)]*ncols; v[f] = Fr(1)
        for i, p in enumerate(piv):
            v[p] = -R[i][f]
        base.append(v)
    return base


def inv(M):
    n = len(M)
    R, piv = rref([row + e for row, e in zip(M, eye(n))])
    assert piv[:n] == list(range(n)), 'singular'
    return [row[n:] for row in R]


def matpow(A, k):
    R = eye(len(A))
    for _ in range(k):
        R = mul(R, A)
    return R


def jordan(blocos, n):
    """matriz de Jordan com blocos [(autovalor, tamanho)]."""
    J = zeros(n, n); i = 0
    for lam, t in blocos:
        for k in range(t):
            J[i + k][i + k] = Fr(lam)
            if k + 1 < t:
                J[i + k][i + k + 1] = Fr(1)
        i += t
    return J


def aleat(r, c, rng, a=3):
    return [[Fr(rng.randint(-a, a)) for _ in range(c)] for _ in range(r)]


def similar(J, rng):
    n = len(J)
    while True:
        P = aleat(n, n, rng, 2)
        if rank(P) == n:
            return mul(mul(P, J), inv(P))


def colunas(vs):
    """matriz cujas colunas sao os vetores vs."""
    return [[v[i] for v in vs] for i in range(len(vs[0]))]


def ensaio(rng, s0, blocosA, blocosB, nA, nB):
    A = similar(jordan(blocosA, nA), rng)
    B = similar(jordan(blocosB, nB), rng)
    # entrelacadores: D (nB x nA) com B D - D A = 0, como sistema linear em vec(D)
    eqs = []
    for i in range(nB):
        for j in range(nA):
            row = [Fr(0)]*(nB*nA)
            for k in range(nB):
                row[k*nA + j] += B[i][k]
            for k in range(nA):
                row[i*nA + k] -= A[k][j]
            eqs.append(row)
    base = nullspace(eqs, nB*nA)
    coef = [Fr(rng.randint(-2, 2)) for _ in base]
    vecD = [sum((c*b[t] for c, b in zip(coef, base)), Fr(0)) for t in range(nB*nA)]
    D = [[vecD[i*nA + j] for j in range(nA)] for i in range(nB)]
    C1 = aleat(nB, nA, rng)
    C0 = add(D, mul(C1, A))
    M1 = C1
    M0 = add(D, mul(B, C1))
    # 1. KC = ML
    ok1 = (add(mul(B, C0), mul(M0, A), 1, -1) == zeros(nB, nA)
           and add(add(C0, mul(B, C1)), add(M0, mul(M1, A)), 1, -1) == zeros(nB, nA)
           and C1 == M1)
    # E = ker (A + s0)^nA
    As0 = add(A, eye(nA), 1, s0)
    E = nullspace(matpow(As0, nA), nA)
    dimE = len(E)
    DE = mul(D, colunas(E)) if E else []
    rDE = rank(DE) if E else 0
    # 3. espaco fisico: u0 em E com D A^j u0 = 0, j = 0..nA-1
    obs = []
    for j in range(nA):
        obs += mul(mul(D, matpow(A, j)), colunas(E)) if E else []
    fis = dimE - (rank(obs) if obs else 0)
    ok3 = fis == dimE - rDE
    # 2. cadeia de Jordan: h0 em ker(A+s0), h1 com (A + s0) h1 = -h0 (se existir)
    ok2 = True
    Cs0 = add(C0, C1, 1, s0)
    N1 = nullspace(As0, nA)
    for h0 in N1:
        # resolver (A + s0) h1 = -h0
        aug = [row + [-h0[i]] for i, row in enumerate(As0)]
        R, piv = rref(aug)
        if nA in piv:
            continue
        h1 = [Fr(0)]*nA
        for i, p in enumerate(piv):
            h1[p] = R[i][nA]
        Dh1 = [sum((D[i][k]*h1[k] for k in range(nA)), Fr(0)) for i in range(nB)]
        alt = [sum((Cs0[i][k]*h1[k] + C1[i][k]*h0[k] for k in range(nA)), Fr(0)) for i in range(nB)]
        ok2 &= Dh1 == alt
    # 4. B + s0 injetivo => D|E = 0
    Bs0 = add(B, eye(nB), 1, s0)
    injetivo = rank(Bs0) == nB
    ok4 = (not injetivo) or rDE == 0
    # 5. atalho errado: base de Jordan com C(s0) h = 0 em cada vetor
    atalho = None
    if E:
        CE = mul(Cs0, colunas(E))
        atalho = dimE - rank(CE)
    return dict(ok1=ok1, ok2=ok2, ok3=ok3, ok4=ok4, dimE=dimE, rank_DE=rDE, fisico=fis,
                injetivo=injetivo, dimD=len(base), atalho_kerC=atalho)


def main():
    rng = random.Random(20261002)
    s0 = Fr(2)
    casos = [
        ('K(s0) nao injetivo, bloco 2', [(-2, 2), (-2, 1), (1, 1), (3, 1)], [(-2, 2), (5, 1), (-1, 1)], 5, 4),
        ('K(s0) injetivo', [(-2, 2), (-2, 1), (1, 1), (3, 1)], [(4, 1), (5, 2), (-1, 1)], 5, 4),
        ('K(s0) nao injetivo, simples', [(-2, 1), (1, 2), (3, 1)], [(-2, 1), (5, 1), (-1, 2)], 4, 4),
        ('L com Jordan 3, K Jordan 2', [(-2, 3), (0, 1)], [(-2, 2), (7, 1)], 4, 3),
    ]
    res = []
    for nome, bA, bB, nA, nB in casos:
        agg = dict(nome=nome, n=0, ok=True, exemplos=[])
        for t in range(25):
            r = ensaio(rng, s0, bA, bB, nA, nB)
            agg['n'] += 1
            agg['ok'] &= r['ok1'] and r['ok2'] and r['ok3'] and r['ok4']
            if t < 4:
                agg['exemplos'].append(r)
            difere = r['atalho_kerC'] is not None and r['atalho_kerC'] != r['fisico']
            agg['atalho_difere'] = agg.get('atalho_difere', 0) + int(difere)
        res.append(agg)
        print(f"{nome}: {agg['n']} sorteios, todas as conferencias {'OK' if agg['ok'] else 'FALHARAM'};"
              f" atalho ker C(s0) difere do fisico em {agg['atalho_difere']} sorteios;"
              f" exemplo {agg['exemplos'][0]}", flush=True)
    json.dump(res, open(sys.argv[1] if len(sys.argv) > 1 else 'r6_quociente.json', 'w'), indent=1, default=str)


if __name__ == '__main__':
    main()
