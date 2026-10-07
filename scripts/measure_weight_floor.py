#!/usr/bin/env python3
"""Piso de Perron da reponderacao por DOF, e o que a familia produto captura.

A norma de coluna ponderada de B sob pesos diagonais w e ||D B D^-1||_1 com
D=diag(w); por Perron-Frobenius o infimo sobre TODOS os w>0 e rho(|B|). Isso
limita de uma vez kappa1, kappa2, pesos de componente e pesos livres por DOF.

DIAGNOSTICO: medido na truncagem 6x18. Para o operador infinito o piso e
MAIOR (acrescentar linhas e colunas nao negativas nao diminui rho).
"""
import sys, json
from fractions import Fraction as Q
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from weighted_norms import enderecos, pesos, MU
from check_free_preconditioner import read_spectral_matrix, read_addresses, free_matrix

ETA = [Q(229, 1024), Q(1), Q(243, 4096), Q(1233, 4096)]


def matriz_B(malha):
    p = Path(f'build/spectrum/rt-A-{malha}.dat')
    dim, raw, expo = read_spectral_matrix(p)
    addr_raw, addr = read_addresses(p, dim), enderecos(p, dim)
    A = [[Q(x)*Q(2)**expo for x in row] for row in raw]
    J = free_matrix(addr_raw, MU)
    return addr, np.array([[abs(float(A[i][j]-J[i][j])) for j in range(dim)]
                           for i in range(dim)])


def main():
    malha = sys.argv[1] if len(sys.argv) > 1 else '6x18'
    addr, B = matriz_B(malha)
    norma = lambda w: float(np.max((B*np.asarray(w, float)[:, None]).sum(0)/np.asarray(w, float)))
    w_eta = [float(pesos(addr, Q(65, 64), Q(5, 4))[i]*ETA[addr[i][0]]) for i in range(len(addr))]
    atual = norma(w_eta)
    rho = float(max(abs(np.linalg.eigvals(B))))
    melhor = min(((norma([(2-(m == 0))*(2-(n == 0))*np.exp(e1)**m*np.exp(e2)**n*float(ETA[c])
                          for (c, _, m, n) in addr]), np.exp(e1), np.exp(e2))
                  for e1 in np.linspace(-0.6, 0.6, 25)
                  for e2 in np.linspace(-1.2, 1.2, 49)), key=lambda t: t[0])
    rep = dict(malha=malha, norma_eta_atual=atual, piso_perron=rho,
               melhor_produto=melhor[0], kappa1=melhor[1], kappa2=melhor[2],
               ganho_produto=atual/melhor[0], ganho_maximo=atual/rho,
               fracao_capturada=float(np.log(atual/melhor[0])/np.log(atual/rho)),
               nota='truncagem; o piso do operador infinito e maior')
    print(json.dumps(rep, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
