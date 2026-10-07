#!/usr/bin/env python3
"""R1 (revisao T2): ||B_L|| na norma do toro com sinal (pesos kappa1^m sqrt(omega2(n)), (65/64, 5/4)),
calculada em truncagens M x N do operador do codigo (signed_operator_L.bloco_B), por Lanczos.
Deve ficar ABAIXO da cota pontual de F3 (3,9641) e reproduzir o 'verdadeiro (12x36) 3,69' de §5.26,
o que mostra que F3 e o L-real usam o MESMO espaco. Testa tambem a norma 'errada' (pesos simetricos
kappa1^|m|, isto e, sem sinal) para ver se a distincao importa."""
import sys
import json
import numpy as np
from scipy.sparse.linalg import LinearOperator, svds

sys.path.insert(0, 'scripts')
import signed_operator_L as sl
import closed_form as cf

K1, K2 = 65/64, 5/4


def norma_B(M, N, g, sinal=True):
    rs = [sl.Radial(m, N) for m in sl.modos(M)]
    offs = np.cumsum([0] + [r.dim for r in rs])
    w = []
    for r in rs:
        rad = np.sqrt(np.concatenate([cf.omega(ns, K2) for _, ns in r.comps]))
        w.append((K1**r.m if sinal else K1**abs(r.m))*rad)

    cache = {}

    def bloco(i, j):
        # o bloco pesado so depende de (d, paridade de m_in): conferido em r5_iso.py (erro 0 exato)
        ro, ri = rs[i], rs[j]
        d = ro.m - ri.m
        if abs(d) > sl.BAND_M:
            return None
        chave = (d, ri.m % 2)
        if chave not in cache:
            X = sl.bloco_B(g, ro, ri)
            cache[chave] = (w[i][:, None]*X)/w[j][None, :]
        return cache[chave]

    blocos = {}
    for i in range(len(rs)):
        for j in range(len(rs)):
            b = bloco(i, j)
            if b is not None and np.any(b):
                blocos[(i, j)] = b
    n = int(offs[-1])

    def mv(x):
        x = np.asarray(x).ravel()
        y = np.zeros(n, complex)
        for (i, j), b in blocos.items():
            y[offs[i]:offs[i+1]] += b @ x[offs[j]:offs[j+1]]
        return y

    def rmv(x):
        x = np.asarray(x).ravel()
        y = np.zeros(n, complex)
        for (i, j), b in blocos.items():
            y[offs[j]:offs[j+1]] += b.conj().T @ x[offs[i]:offs[i+1]]
        return y

    A = LinearOperator((n, n), matvec=mv, rmatvec=rmv, dtype=complex)
    s = svds(A, k=1, return_singular_vectors=False, tol=1e-8, maxiter=3000)
    return float(s[0]), n


def main():
    g = sl.campos()
    res = []
    for M, N, sinal in ((12, 36, True), (24, 72, True), (40, 120, True), (60, 180, True), (80, 240, True)):
        v, n = norma_B(M, N, g, sinal)
        res.append(dict(M=M, N=N, sinal=sinal, dim=n, norma_B=v))
        print(f'{M}x{N} ({"com sinal" if sinal else "pesos |m|"}), dim {n}: ||B_L|| truncada = {v:.5f}   (F3: <= 3,9641)', flush=True)
    json.dump(res, open(sys.argv[1] if len(sys.argv) > 1 else 'r1_BL.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
