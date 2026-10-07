#!/usr/bin/env python3
"""Revisao R5: (a) confere free_tail(M, N, k2') <= cota explicita de HREC_IDA.md em amostra grande de
(M, N, k2') racionais (a prova em papel esta no relatorio); (b) testa a extensao do Lema R5 a |s| <= S para
S > 3,1 (lacuna de escopo: a bijecao do par. 1 e a contagem de T2 usam Re s > 0 sem limite): para cada S ha
uma caixa G com contracao fraca e forte, porque (beta + S) * cota -> 0 e (beta+ + S) * structured_input_tail -> 0.
Aritmetica racional exata (fractions)."""
import random
import sys
from fractions import Fraction as Q
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from sharp_exterior import free_tail  # noqa: E402
from structured_exterior import structured_input_tail  # noqa: E402

BETA_P = Q(5191, 500)
MU_MIN, MU_MAX = Q(1, 6), Q(17, 100)


def D(k):
    return 2*k/(k*k - 1)


def cota(M, N, k2):
    r, b = 1/k2, Q(M, 2)
    return max(Q(3, 2)*(1 + 2*MU_MAX)/(b*(1 - r)), Q(3, 2)*(1 + 4/(1 - r))/(MU_MIN*(N + 1)))


def main():
    rnd = random.Random(20261003)
    viol = 0
    for _ in range(4000):
        M = rnd.choice([2, 3, 4, 7, 16, 100, 1024, 10**6]) + rnd.randint(0, 5)
        N = rnd.randint(2, 50) if rnd.random() < .5 else rnd.choice([100, 4096, 10**7])
        k2 = 1 + Q(1, rnd.randint(2, 10**6))
        if rnd.random() < .2:
            k2 = Q(rnd.randint(101, 300), 100)
        ft, ce = free_tail(M, N, radial=k2, mu_min=MU_MIN, mu_max=MU_MAX), cota(M, N, k2)
        viol += ft > ce
    print(f'(a) free_tail <= cota em 4000 amostras (M >= 2, N >= 2, k2 > 1): violacoes = {viol}')
    ok = viol == 0
    print('(b) extensao para |s| <= S, k2 = 9/8 e 1 + 2^-8:')
    for S in (Q(31, 10), Q(10), Q(50), Q(200)):
        for k2 in (Q(9, 8), 1 + Q(1, 256)):
            bw = BETA_P*D(k2)/D(Q(5, 4))
            M = 2
            while (bw + S)*cota(M, 2, k2) >= Q(9, 10) and (bw + S)*Q(3, 2)*(1 + 2*MU_MAX)/(Q(M, 2)*(1 - 1/k2)) >= Q(9, 10):
                M *= 2
            N = 4
            while (bw + S)*Q(3, 2)*(1 + 4/(1 - 1/k2))/(MU_MIN*(N + 1)) >= Q(9, 10):
                N *= 2
            fraco = (bw + S)*free_tail(M, N, radial=k2, mu_min=MU_MIN, mu_max=MU_MAX)
            forte = (BETA_P + S)*structured_input_tail(M, N)
            while forte >= 1:     # a caixa forte tambem tem de contrair
                M, N = 2*M, 2*N
                fraco = (bw + S)*free_tail(M, N, radial=k2, mu_min=MU_MIN, mu_max=MU_MAX)
                forte = (BETA_P + S)*structured_input_tail(M, N)
            good = fraco < 1 and forte < 1
            ok &= good
            print(f'   S = {float(S):6.1f}, k2 = {str(k2):8s}: G = {M} x {N}, fraco {float(fraco):.3f}, '
                  f'forte {float(forte):.3f}  {"OK" if good else "FALHA"}')
    print('CONFERE' if ok else 'FALHA')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
