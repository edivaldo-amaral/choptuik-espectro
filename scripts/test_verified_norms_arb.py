#!/usr/bin/env python3
"""Auditoria: as normas verificadas em ponto flutuante (verified_norms, criterio de Rump)
contra uma cota INDEPENDENTE em aritmetica de bolas (Arb, python-flint).

Para G = X X^* (hermitiana >= 0) e W uma base de autovetores calculada em ponto flutuante:
S = W^* W e W^* G W = Lambda + E sao calculados em bolas; por congruencia (W invertivel),
    lambda_max(G) = max_y (y^* W^*GW y)/(y^* S y) <= (max Lambda + ||E||_F)/(1 - ||S - I||_F).
Nao usa Rump nem a analise de arredondamento do numpy. Matrizes: blocos reais do
certificado (configuracao pequena) e matrizes aleatorias mal condicionadas."""
from __future__ import annotations

from pathlib import Path
import sys
import time

import numpy as np
from flint import acb_mat, acb, arb, ctx

sys.path.insert(0, str(Path(__file__).resolve().parent))
from verified_norms import cota_norma

ctx.prec = 100


def para_acb(A):
    return acb_mat([[acb(float(z.real), float(z.imag)) for z in linha] for linha in np.asarray(A, complex).tolist()])


def fro_sup(M):
    s = arb(0)
    for linha in M.tolist():
        for z in linha:
            s += z.real.abs_upper()**2 + z.imag.abs_upper()**2      # >= 0 mesmo com bolas em 0
    return float(s.sqrt().upper())


def norma_arb_sup(X):
    """Cota superior de ||X||_2 por congruencia em bolas."""
    X = np.asarray(X, complex)
    if X.shape[0] > X.shape[1]:
        X = X.conj().T                          # ||X|| = ||X^*||; G = X X^* do lado menor
    G = X @ X.conj().T
    lam, W = np.linalg.eigh((G + G.conj().T)/2)
    Xa, Wa = para_acb(X), para_acb(W)
    Ga = Xa*Xa.transpose().conjugate()          # G exata em bolas
    WH = Wa.transpose().conjugate()
    S = WH*Wa
    n = S.nrows()
    for i in range(n):
        S[i, i] -= 1
    E = WH*Ga*Wa
    for i in range(n):
        E[i, i] -= float(lam[i])
    eS, eE = fro_sup(S), fro_sup(E)
    assert eS < 0.5
    return float(((float(lam.max()) + eE)/(1 - eS))**0.5)*(1 + 1e-12)


def main():
    import tile_perron as tp
    tp.RIG = True
    c = tp.Ctx((8, 24), (12, 36), (32, 77), 4)
    sil = lambda *a, **k: None
    qZ, V = tp.parte_Z(c, 0.3+0j, sil)
    q1, HZ1, QZ1 = tp.parte_T1(c, 0.3+0j, sil)
    rng = np.random.default_rng(7)
    U, _ = np.linalg.qr(rng.standard_normal((150, 150)) + 1j*rng.standard_normal((150, 150)))
    mal = U @ np.diag(np.logspace(0, -12, 150)) @ U.conj().T      # condicionamento 1e12
    casos = {'V (Z inverso, 252)': V, 'V H_Z,T1 (252 x 342)': V @ HZ1,
             'V Q_Z,T1': V @ QZ1, 'aleatoria mal condicionada (150)': mal,
             'posto 1 (150)': np.outer(rng.standard_normal(150), rng.standard_normal(150))}
    pior = 0.0
    for nome, X in casos.items():
        t0 = time.time()
        rump = cota_norma(X)
        arb_sup = norma_arb_sup(X)
        verdade = np.linalg.norm(X, 2)
        ok = rump >= arb_sup*(1 - 1e-12) or rump >= verdade
        pior = max(pior, arb_sup/rump)
        print(f'{nome:>34}: Rump {rump:.12f}  Arb(sup) {arb_sup:.12f}  (float {verdade:.12f})'
              f'  Rump >= Arb: {rump >= arb_sup}  [{time.time()-t0:.0f}s]', flush=True)
    print(f'maior Arb/Rump = {pior:.12f}  ->  ' + ('OK' if pior <= 1 else 'RUMP ABAIXO DA COTA EM BOLAS'))


if __name__ == '__main__':
    main()
