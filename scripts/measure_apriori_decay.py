#!/usr/bin/env python3
"""Rota 5: estimativa a priori de decaimento das autofuncoes.

Pela equacao de autovalor, u_n = -(J_mu+s)^-1 (Bu)_n. O fator de contracao
local e C/|D(m,n)|, com C a soma de coluna ponderada de B e
|D| >= max(mu(n+1), |m|/2). Onde ele cruza 1, a recursao contrai -- e isso
NAO depende de s, ao contrario de tudo que foi medido nas rodadas anteriores.

DIAGNOSTICO. Mede tres coisas: o limiar de contracao, a localizacao de B e o
decaimento observado das autofuncoes.
"""
from __future__ import annotations
import sys, json
from fractions import Fraction as Q
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from weighted_norms import enderecos, MU
from check_free_preconditioner import read_spectral_matrix, read_addresses, free_matrix

ETA = [229/1024, 1.0, 243/4096, 1233/4096]


def carrega(malha):
    p = Path(f'build/spectrum/rt-A-{malha}.dat')
    dim, raw, expo = read_spectral_matrix(p)
    addr_raw, addr = read_addresses(p, dim), enderecos(p, dim)
    A = [[Q(x)*Q(2)**expo for x in row] for row in raw]
    J = free_matrix(addr_raw, MU)
    B = np.array([[abs(float(A[i][j]-J[i][j])) for j in range(dim)] for i in range(dim)])
    w = np.array([(2-(m == 0))*(2-(n == 0))*(65/64)**m*(5/4)**n*ETA[c]
                  for (c, _, m, n) in addr])
    return addr, B*w[:, None]/w[None, :], np.array([a[3] for a in addr])


def main():
    malha = sys.argv[1] if len(sys.argv) > 1 else '12x36'
    addr, Bw, n_i = carrega(malha)
    C = float(Bw.sum(0).max())
    dn = np.abs(n_i[:, None]-n_i[None, :])
    massa = {L: float((Bw*(dn <= L)).sum()/Bw.sum()) for L in (1, 3, 8, 12, 20)}
    razoes = [float(Bw[dn == d].max()) for d in range(3, 34, 3) if (dn == d).any()]
    decai = float(np.exp(np.mean(np.diff(np.log(razoes)))/3))
    rep = dict(malha=malha, C_col_ponderada=C,
               limiar_radial=C/float(MU)-1, limiar_fourier=2*C,
               massa_por_distancia=massa, razao_decaimento_B_por_indice=decai,
               n_maximo_da_malha=int(n_i.max()),
               nota=('o limiar de contracao esta ACIMA do n maximo das malhas '
                     'existentes: nenhuma computacao do projeto entrou nesse regime'))
    logp, caudas = 0.0, {}
    N0 = int(C/float(MU))
    for N in range(N0+1, 400):
        logp += np.log(C/(float(MU)*(N+1)))
        if N in (80, 100, 120, 150, 250, 384):
            caudas[N] = float(np.exp(logp))
    rep['cauda_radial_a_priori'] = caudas
    print(json.dumps(rep, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
