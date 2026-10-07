#!/usr/bin/env python3
"""Piso de similaridade GERAL contra piso de similaridade DIAGONAL.

O argumento que fecha as nove alavancas de reponderacao usa Perron-Frobenius:
    inf sobre pesos diagonais w>0 de ||D B D^-1||_1  =  rho(|B|).
Mas isso so cobre similaridades DIAGONAIS. Para similaridade geral,
    inf sobre S invertivel de ||S B S^-1||  =  rho(B),
e rho(B) <= rho(|B|) com desigualdade ESTRITA sempre que houver cancelamento.

A razao rho(|B|)/rho(B) mede quanto uma similaridade nao-diagonal poderia
comprar alem de TODA reponderacao possivel. Se for ~1, a direcao inteira esta
fechada e nao so a parte diagonal. Se for grande, ha uma lacuna real no mapa.
"""
import sys
from fractions import Fraction as Q
from pathlib import Path
import numpy as np
sys.path.insert(0, 'scripts')
from check_free_preconditioner import (read_spectral_matrix, read_addresses,
                                       free_matrix)

MU = Q(722873400, 2**32)


def matriz_B(malha):
    p = Path(f'build/spectrum/rt-A-{malha}.dat')
    dim, raw, expo = read_spectral_matrix(p)
    addr_raw = read_addresses(p, dim)
    esc = 2.0 ** expo
    J = free_matrix(addr_raw, MU)
    B = np.empty((dim, dim))
    for i in range(dim):
        linha = raw[i]
        Ji = J[i]
        for j in range(dim):
            B[i, j] = float(linha[j]) * esc - float(Ji[j])
    return B


def main(malhas):
    print('piso diagonal rho(|B|) contra piso geral rho(B)')
    print(f"{'malha':>7} {'dim':>6} {'rho(B)':>12} {'rho(|B|)':>12} {'razao':>8}")
    print('-' * 50)
    for malha in malhas:
        B = matriz_B(malha)
        rb = float(np.max(np.abs(np.linalg.eigvals(B))))
        ra = float(np.max(np.abs(np.linalg.eigvals(np.abs(B)))))
        print(f'{malha:>7} {B.shape[0]:>6} {rb:>12.6f} {ra:>12.6f} {ra/rb:>7.3f}x')


if __name__ == '__main__':
    main(sys.argv[1:] or ['4x12', '6x18', '8x24', '10x30', '12x36'])


def condicionamento(malhas):
    """O infimo rho(B) e aproximado por S mal-condicionada?

    Se B = X Lambda X^-1 entao S=X^-1 da ||S B S^-1||_2 = rho(B) exatamente.
    O preco e cond(X): todo OUTRO termo do argumento e multiplicado por ele.
    Este projeto ja foi mordido por esse modo de falha uma vez (peso balanceado
    por Perron com condicionamento 10^151), entao a medida e obrigatoria.
    """
    print()
    print('preco da similaridade nao-diagonal')
    print(f"{'malha':>7} {'rho(B)':>10} {'cond(X)':>12} {'defeito':>10}")
    print('-' * 44)
    for malha in malhas:
        B = matriz_B(malha)
        vals, X = np.linalg.eig(B)
        c = float(np.linalg.cond(X))
        # menor gap relativo entre autovalores: mede proximidade de Jordan
        v = np.sort_complex(vals)
        gap = float(np.min(np.abs(np.diff(v)))) if len(v) > 1 else float('nan')
        print(f'{malha:>7} {float(np.max(np.abs(vals))):>10.6f} {c:>12.4e} {gap:>10.2e}')
