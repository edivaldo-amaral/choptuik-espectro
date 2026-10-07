#!/usr/bin/env python3
"""Certifica R_J = -lambda_min(Herm(J_mu)) na realizacao l^2 com pesos, exato.

R_J e a borda direita do campo de valores do operador LIVRE, e entra em
  spec(L) subset W(L)  =>  nao ha espectro para Re s > R_J + ||B||_2.
Ver ALTERNATIVE_ROUTES 5.10-5.12.

Duas observacoes tornam isto exato:

1. O minimo vive num problema 1-D radial da componente 3: diagonal mu(n+1),
   acoplamento 2*mu*n para TODO n'<n, pesos w_n = (2-[n=0]) kappa2^n.
   Estabiliza em N=18 (secao 5.11.3).

2. Conjugar por sqrt(w) introduz irracionais. Evita-se com o pencil
   generalizado: lambda_min(Herm(D J D^-1)) com D=diag(sqrt(w)) e o menor
   lambda de  Herm(W J) y = lambda W y,  W = diag(w) -- e W J e W sao
   RACIONAIS. Prova-se lambda_min >= lambda0 exibindo lambda0 racional com
   Herm(W J) - lambda0 W definida positiva, por LDL^T exato de pivos positivos.

DIAGNOSTICO do bloco finito. A cauda n>N ainda precisa de argumento separado.
"""
from __future__ import annotations
import argparse, json
from fractions import Fraction as Q

MU = Q(722873400, 2**32)
K2 = Q(5, 4)


def bloco(N: int):
    """W (diagonal) e Herm(W J), exatos, do modelo 1-D da componente 3."""
    w = [Q(2 if n else 1)*K2**n for n in range(N)]
    J = [[Q(0)]*N for _ in range(N)]
    for n in range(N):
        J[n][n] = MU*(n+1)
        for lo in range(n-1, -1, -1):
            J[lo][n] = 2*MU*n
    WJ = [[w[i]*J[i][j] for j in range(N)] for i in range(N)]
    H = [[(WJ[i][j] + WJ[j][i])/2 for j in range(N)] for i in range(N)]
    return w, H


def definida_positiva(M):
    """LDL^T exato sem pivotamento; devolve (ok, menor pivo)."""
    n = len(M)
    A = [row[:] for row in M]
    menor = None
    for k in range(n):
        p = A[k][k]
        if p <= 0:
            return False, p
        menor = p if menor is None else min(menor, p)
        for i in range(k+1, n):
            f = A[i][k]/p
            if f:
                for j in range(k, n):
                    A[i][j] -= f*A[k][j]
    return True, menor


def certifica(N: int, lam0: Q):
    w, H = bloco(N)
    M = [[H[i][j] - (lam0*w[i] if i == j else Q(0)) for j in range(N)] for i in range(N)]
    return definida_positiva(M)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--N', type=int, nargs='*', default=[18, 20, 30, 40])
    ap.add_argument('--lambda0', type=str, default='-1777/20000')
    ap.add_argument('--out', default=None)
    args = ap.parse_args()
    lam0 = Q(args.lambda0)
    print(f"lambda0 = {lam0} = {float(lam0):.9f}")
    print(f"se Herm(WJ) - lambda0 W for definida positiva, entao lambda_min > lambda0,")
    print(f"logo R_J < {float(-lam0):.9f}\n")
    print(f"{'N':>5} {'definida positiva?':>20} {'menor pivo':>16}")
    print('-'*44)
    linhas = []
    for N in args.N:
        ok, piv = certifica(N, lam0)
        print(f"{N:>5} {('SIM' if ok else 'nao'):>20} {float(piv):>16.6e}")
        linhas.append(dict(N=N, positiva=ok, menor_pivo=str(piv)))
    rep = dict(lambda0=str(lam0), R_J_bound=str(-lam0), blocos=linhas,
               nota='bloco finito; a cauda n>N exige argumento separado',
               mu=str(MU), kappa2=str(K2))
    if args.out:
        with open(args.out, 'w') as fh:
            json.dump(rep, fh, indent=2)
        print('\nescrito:', args.out)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
