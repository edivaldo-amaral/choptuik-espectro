#!/usr/bin/env python3
"""Recomputo INDEPENDENTE (auditoria de 16/09) da inversa livre de RT.

Nao importa nenhum modulo do repositorio. Tudo em fractions.Fraction.

O que este script faz:
  1. Constroi colunas EXATAS de Q(z0)=(O_mu+z0)^-1 por substituicao regressiva
     a partir da definicao de RT (arXiv:1203.3766v1, eq. (22) e a expressao
     de D'_A, D'_B em termos de Band/BandI, secao 3):
        (O_mu c)_j = (i m/2 + mu (j+1)) c_j + sum_{i>j} 2 mu i c_i
     (para K=O^B a soma e so sobre i = j mod 2).
     Q preserva m e so baixa n, logo a coluna (m,n) e EXATA em linhas 0..n:
     normas de coluna calculadas aqui nao sao truncadas.
  2. Norma de coluna nos DOFs reais complexificados, pesos RT
     (2-delta_n0) kappa2^n (o peso de Fourier cancela porque Q preserva m):
        sum_j w_j (|Re c_j|+|Im c_j|) / w_n.
  3. Rederivacao propria do majorante de cauda de ENTRADA ||Q(I-G)||
     (ver docs/AUDITORIA_MATEMATICA_16SET.md, secao A1) e comparacao com
     colunas exatas numa faixa exterior amostrada.
  4. Constantes de RT publicadas -> recuperacao de raio em 1024x4096 sem
     usar beta<10.382 (norma de componentes) do produtor auditado.
  5. Recomputo dos defeitos 117774/512125 e 8046603/10242500 a partir das
     constantes DECLARADAS em RADIUS_RECOVERY.md.

Ponto flutuante: nenhum.
"""
from __future__ import annotations

from fractions import Fraction as Q
import json

MU_REFA = Q(722873400, 2**32)          # RT secao 5
MU_MIN, MU_MAX = Q(1, 6), Q(17, 100)    # intervalo usado pelas notas
RT_L1 = Q(21, 2)                        # RT secao 5, bound de (39a)
RT_L2 = Q(1, 2**25)                     # RT (39b)
RT_R = Q(1, 2**277)                     # RT (60)
RT_REFB_MU = Q(43, 2**40)               # RT secao 5: |mu_refB| <= 43*2^-40


class C:
    """Complexo com partes racionais."""
    __slots__ = ('re', 'im')

    def __init__(self, re, im=Q(0)):
        self.re, self.im = Q(re), Q(im)

    def __add__(self, o):
        return C(self.re+o.re, self.im+o.im)

    def __sub__(self, o):
        return C(self.re-o.re, self.im-o.im)

    def __mul__(self, o):
        return C(self.re*o.re-self.im*o.im, self.re*o.im+self.im*o.re)

    def scale(self, x):
        return C(self.re*x, self.im*x)

    def inv(self):
        d = self.re*self.re+self.im*self.im
        return C(self.re/d, -self.im/d)

    def l1(self):
        return abs(self.re)+abs(self.im)


def exact_column(kind: str, m: int, n: int, mu: Q, z0: Q = Q(0)) -> list[C]:
    """Coluna n de (O_mu+z0)^-1 no modo de Fourier m. kind='J' ou 'K'."""
    if kind not in ('J', 'K') or n < 0 or mu <= 0 or z0 < 0:
        raise ValueError('parametros invalidos')
    step = 1 if kind == 'J' else 2
    c = [C(0) for _ in range(n+1)]
    running = [C(0), C(0)]  # soma de 2 mu i c_i por paridade de i
    for j in range(n, -1, -1):
        if (n-j) % step:
            continue
        diag = C(mu*(j+1)+z0, Q(m, 2))
        rhs = C(1 if j == n else 0)
        tail = running[0]+running[1] if step == 1 else running[j % 2]
        c[j] = (rhs-tail)*diag.inv()
        running[j % 2] = running[j % 2]+c[j].scale(2*mu*j)
    return c


def weight(n: int, kappa2: Q) -> Q:
    return (2-(n == 0))*kappa2**n


def column_norm(column: list[C], kappa2: Q) -> Q:
    n = len(column)-1
    return sum((weight(j, kappa2)*x.l1() for j, x in enumerate(column)), Q(0))/weight(n, kappa2)


def box_norm(M: int, N: int, mu: Q, kappa2: Q) -> tuple[Q, tuple]:
    """max de colunas exatas com 0<=m<M, 0<=n<N (J e K). Superconjunto dos DOFs reais."""
    best, where = Q(0), None
    for kind in ('J', 'K'):
        for m in range(M):
            for n in range(N):
                value = column_norm(exact_column(kind, m, n, mu), kappa2)
                if value > best:
                    best, where = value, (kind, m, n)
    return best, where


def tail_bound(M: int, N: int, kappa2: Q, mu_min: Q = MU_MIN, mu_max: Q = MU_MAX) -> Q:
    """Rederivacao propria: sup das colunas com m>=M OU n>=N, base z0=0 real.

    Coluna (m,n): e_n/D_n + q (e_{n-k} - f e_{n-2k} + ...), |f|<=1,
    q = -2 mu n/(D_n D_{n-k}); soma ponderada <= 1/|D_n| + |q| r^k/(1-r^k),
    r=1/kappa2; realificacao <= sqrt(2) <= 3/2.
    (a) |m|>=M, b=m/2>=M/2: k=1 -> 1/(b(1-r)); k=2 -> 1/(b(1-r^2))+2mu_max r^2/(b^2(1-r^2)).
    (b) n>=N>=2: |D_j|>=mu(j+1): k=1 -> (1+2r/(1-r))/(mu(n+1));
        k=2 -> (1+2n r^2/((n-1)(1-r^2)))/(mu(n+1)), decrescente em n.
    """
    if M < 1 or N < 2 or kappa2 <= 1:
        raise ValueError('M>=1, N>=2, kappa2>1')
    r, b = 1/kappa2, Q(M, 2)
    fourier = max(1/(b*(1-r)), 1/(b*(1-r*r))+2*mu_max*r*r/(b*b*(1-r*r)))
    radial = max(1+2*r/(1-r), 1+Q(2*N, N-1)*r*r/(1-r*r))/(mu_min*(N+1))
    return Q(3, 2)*max(fourier, radial)


def sampled_exterior_max(M: int, N: int, mu: Q, kappa2: Q, dm: int, dn: int) -> Q:
    """Colunas exatas na faixa exterior M<=m<M+dm (n<N+dn) ou N<=n<N+dn (m<M+dm)."""
    best = Q(0)
    for kind in ('J', 'K'):
        for m in range(M+dm):
            for n in range(N+dn):
                if m >= M or n >= N:
                    best = max(best, column_norm(exact_column(kind, m, n, mu), kappa2))
    return best


def rt_K1(kappa2: Q) -> Q:
    """RT (36): ||S Gamma1|| <= 4 kappa2^-1 (1-kappa2^-2)^-1."""
    return 4/kappa2/(1-1/kappa2**2)


def rt_K2(mu: Q, kappa2: Q) -> Q:
    """RT (37a): ||O_mu^-1|| <= 2 mu^-1 (1-kappa2^-1)^-1 (1+kappa2^-1)."""
    return 2/mu/(1-1/kappa2)*(1+1/kappa2)


def rt_only_beta(kappa2: Q) -> Q:
    """||B_true|| com B=mu S Gamma1 + 2 S Xi Gamma2(omega,.), so constantes RT.

    mu_true<=17/100; ||S Xi Gamma2(RefA,.)||<=L1=10.5 (vale em pesos menores:
    cada norma de campo em (33) decresce); correcao 2*3*(L2+R) por (34),(39b),(60).
    """
    return MU_MAX*rt_K1(kappa2)+2*RT_L1+6*(RT_L2+RT_R)


def radius_recovery(M=1024, N=4096, radius=Q(5, 4)):
    strong_t, weak_t = tail_bound(M, N, Q(5, 4)), tail_bound(M, N, Q(9, 8))
    declared_strong_beta = Q(5191, 500)
    declared_weak_beta = Q(162, 85)*declared_strong_beta
    rt_strong, rt_weak = rt_only_beta(Q(5, 4)), rt_only_beta(Q(9, 8))
    out = dict(
        box=[M, N],
        tail_strong=str(strong_t), tail_weak=str(weak_t),
        declared_strong_defect=str((declared_strong_beta+radius)*strong_t),
        declared_weak_defect=str((declared_weak_beta+radius)*weak_t),
        rt_only_strong_beta=str(rt_strong), rt_only_weak_beta=str(rt_weak),
        rt_only_strong_defect=str((rt_strong+radius)*strong_t),
        rt_only_weak_defect=str((rt_weak+radius)*weak_t),
        rt_only_strong_decimal=float((rt_strong+radius)*strong_t),
        rt_only_weak_decimal=float((rt_weak+radius)*weak_t),
    )
    out['declared_matches_report'] = (
        Q(out['declared_strong_defect']) == Q(117774, 512125)
        and Q(out['declared_weak_defect']) == Q(8046603, 10242500))
    out['rt_only_both_contract'] = (rt_strong+radius)*strong_t < 1 and (rt_weak+radius)*weak_t < 1
    return out


def main():
    result = {}
    # 1/mu e exatamente a norma da coluna (0,0) de J; comparar com os relatorios.
    for label, M, N, k2 in (('strong_6x18', 6, 18, Q(5, 4)), ('weak_12x36', 12, 36, Q(9, 8))):
        value, where = box_norm(M, N, MU_REFA, k2)
        result[f'exact_QF_norm_{label}'] = dict(value=str(value), decimal=float(value),
                                               argmax=where, equals_inverse_mu=value == 1/MU_REFA)
    crown_bound = Q(3266385568255, 549755813888)  # crown-contour-6x18-background.json
    result['crown_report_free_bound_is_upper'] = 1/MU_REFA <= crown_bound
    # Cauda: formula propria >= colunas exatas amostradas, tres valores de mu.
    checks = []
    for M, N, k2 in ((6, 18, Q(5, 4)), (12, 36, Q(9, 8)), (4, 12, Q(5, 4))):
        bound = tail_bound(M, N, k2)
        for mu in (MU_MIN, MU_REFA, MU_MAX):
            sample = sampled_exterior_max(M, N, mu, k2, 8, 24)
            checks.append(dict(box=[M, N], kappa2=str(k2), mu=str(mu), exact_sample_max=float(sample),
                               bound=float(bound), dominated=sample <= bound))
    result['tail_formula_vs_exact_columns'] = checks
    result['tail_176x944_strong'] = str(tail_bound(176, 944, Q(5, 4)))
    result['tail_6x18_strong'] = str(tail_bound(6, 18, Q(5, 4)))
    result['radius_recovery'] = radius_recovery()
    result['rt_constants'] = dict(K1_strong=str(rt_K1(Q(5, 4))), K1_weak=str(rt_K1(Q(9, 8))),
                                  K2_refA_strong=str(rt_K2(MU_REFA, Q(5, 4))),
                                  K2_refA_weak=str(rt_K2(MU_REFA, Q(9, 8))))
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
