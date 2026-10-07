#!/usr/bin/env python3
"""Analise de teto do defeito TT: existe piso que feche a rota do majorante?

O criterio de blocos de FREQUENCY_BLOCKS (existe r>0 com sum_i r_i e_ij < r_j
para cada coluna j) equivale, para E nao negativa, a rho(E)<1. A pergunta do
teto e portanto de Perron sobre a matriz majorante de blocos, nao sobre a
matriz infinita de DOFs.

Blocos: P (nucleo), S (coroa), T (exterior infinito). Com a parametriz
acoplada V sobre G=P+S, os quatro blocos finitos tem defeito a_fin, e os
acoplamentos com o exterior sao os de FREQUENCY_BLOCKS §7:

    e_TS = beta * t(P)              (nao depende de G)
    e_ST = v * (beta+rho) * t(G)
    e_TT =     (beta+rho) * t(G)

com t a cauda de ENTRADA de Q0 fora da caixa. Reduzindo ao criterio de dois
blocos com F=P+S, a condicao e

    e_TT < 1  e  e_TS*e_ST < (1-a_fin)*(1-e_TT).

Como e_TS nao depende de G e e_ST, e_TT sao proporcionais a t(G) -> 0, a
condicao SEMPRE passa para G grande: nao ha teto. O que existe e uma caixa
minima, que este script calcula exatamente em funcao de beta e do nucleo.

DIAGNOSTICO de dimensionamento. Nao certifica nada; v e a_fin sao parametros.
"""
from __future__ import annotations

import argparse
import json
from fractions import Fraction as Q

RHO = Q(5, 4)              # |s-z0| no contorno A atual, com z0=0
MU_MIN = Q(1, 6)

# niveis de majorante para ||B|| na norma correspondente
BETAS = {
    'uniforme (ALGEBRAIC_TAIL)': Q(23),
    'componentes (COMPONENT_EXTERIOR)': Q(5191, 500),      # 10,382
}


def t_structured(M: int, N: int) -> Q:
    """Cauda de entrada de Q0 fora de |m|<M, n<N (STRUCTURED_EXTERIOR)."""
    return max(Q(15, M), Q(81, N+1))


def t_algebraic(M: int, N: int) -> Q:
    """Cauda anterior, de ALGEBRAIC_TAIL, com z0=0 e mu_min=1/6."""
    return 18*max(Q(2, M), Q(1, 1)/(MU_MIN*(N-1)))


TAILS = {'estruturada': t_structured, 'algebrica': t_algebraic}


def min_box(beta: Q, core: tuple[int, int], tail, v: Q, a_fin: Q):
    """Menor caixa G que satisfaz o criterio de blocos, dado o nucleo P."""
    tP = tail(*core)
    eTS = beta*tP
    # x = (beta+rho)*t(G); condicao: x*(eTS*v + 1 - a_fin) < 1 - a_fin
    denom = eTS*v + 1 - a_fin
    x_max = (1 - a_fin)/denom                 # supremo aberto
    t_max = x_max/(beta + RHO)
    if t_max <= 0:
        return None
    # t(G)=max(15/M, 81/(N+1)) < t_max  =>  M > 15/t_max, N+1 > 81/t_max
    if tail is t_structured:
        M = int(Q(15)/t_max) + 1
        N = int(Q(81)/t_max) + 1
        while t_structured(M, N) >= t_max:
            M += 1
            N += 1
    else:
        M = int(Q(36)/t_max) + 1
        N = int(Q(108)/(MU_MIN*t_max)) + 2
        while t_algebraic(M, N) >= t_max:
            M += 1
            N += 1
    return dict(M=M, N=N, t_P=tP, e_TS=eTS, x_max=x_max, t_max=t_max)



def min_box_self_consistent(beta: Q, tail, v: Q, a_fin: Q, k: int = 2,
                            aspect: Q = Q(384, 107), cap: int = 40000):
    """Menor caixa com o NUCLEO otimizado junto, exigindo folga G >= k*P.

    Fixar P em 107x384 subestima o nucleo e infla a caixa. Aqui P cresce
    junto, e a caixa final e o MAXIMO entre o que o criterio exige e k*P:
    quando a separacao morde, e ela que dimensiona, nao o criterio.
    """
    M = 40
    while M < cap:
        N = int(M*aspect)
        r = min_box(beta, (M, N), tail, v, a_fin)
        if r['M'] <= k*M and r['N'] <= k*N:
            return dict(core=(M, N), box=(max(r['M'], k*M), max(r['N'], k*N)),
                        criterio=(r['M'], r['N']), separacao_morde=r['M'] <= k*M)
        M = int(M*1.04) + 1
    return None


def dofs(M: int, N: int) -> int:
    """Escala de graus de liberdade, normalizada pela caixa 107x384=143424."""
    return round(143424*(M*N)/(107*384))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--v', type=str, default='1',
                    help='majorante da norma da parametriz finita V (otimista: 1)')
    ap.add_argument('--a-fin', type=str, default='1/1000000000')
    ap.add_argument('--out', type=str, default=None)
    args = ap.parse_args()
    v, a_fin = Q(args.v), Q(args.a_fin)

    cores = [(12, 36), (24, 72), (64, 192), (107, 384), (214, 768)]
    linhas = []
    for nome_t, tail in TAILS.items():
        for nome_b, beta in BETAS.items():
            for core in cores:
                r = min_box(beta, core, tail, v, a_fin)
                if r is None:
                    continue
                d = dofs(r['M'], r['N'])
                linhas.append(dict(
                    cauda=nome_t, majorante=nome_b, beta=str(beta),
                    nucleo=f"{core[0]}x{core[1]}",
                    t_P=str(r['t_P']), e_TS=str(r['e_TS']),
                    caixa_G=f"{r['M']}x{r['N']}", M=r['M'], N=r['N'],
                    dofs_estimados=d, fator_sobre_107x384=round(d/143424, 1)))

    print(f"{'cauda':<12} {'majorante':<34} {'nucleo':>9} {'e_TS':>9} "
          f"{'caixa G minima':>16} {'x DOFs':>8}")
    print('-'*96)
    for l in linhas:
        print(f"{l['cauda']:<12} {l['majorante']:<34} {l['nucleo']:>9} "
              f"{float(Q(l['e_TS'])):>9.3f} {l['caixa_G']:>16} "
              f"{l['fator_sobre_107x384']:>8}")

    print()
    print('--- nucleo otimizado junto, com folga G >= 2P ---')
    print(f"{'majorante':<34} {'nucleo P':>12} {'caixa G':>14} {'x DOFs':>9}")
    print('-'*74)
    auto = []
    for nome_b, beta in BETAS.items():
        r = min_box_self_consistent(beta, t_structured, v, a_fin)
        if r is None:
            continue
        d = dofs(*r['box'])
        auto.append(dict(majorante=nome_b, beta=str(beta),
                         nucleo=f"{r['core'][0]}x{r['core'][1]}",
                         caixa_G=f"{r['box'][0]}x{r['box'][1]}",
                         fator_sobre_107x384=round(d/143424, 1)))
        print(f"{nome_b:<34} {r['core'][0]:>5}x{r['core'][1]:<6} "
              f"{r['box'][0]:>6}x{r['box'][1]:<7} {d/143424:>9.1f}")

    rep = dict(
        nota='DIAGNOSTICO de dimensionamento; v e a_fin sao parametros, nada e certificado',
        data='2026-09-20', v=str(v), a_fin=str(a_fin), rho=str(RHO),
        conclusao=('nao ha teto: e_TS nao depende de G e e_ST,e_TT sao proporcionais '
                   'a t(G)->0, logo o criterio passa para G grande. O que existe e '
                   'uma caixa minima, tabulada abaixo.'),
        linhas=linhas, nucleo_otimizado=auto)
    if args.out:
        with open(args.out, 'w') as fh:
            json.dump(rep, fh, indent=2)
        print('\nescrito:', args.out)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
