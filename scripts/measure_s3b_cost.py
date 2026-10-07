#!/usr/bin/env python3
"""DIAGNOSTICO EM PONTO FLUTUANTE -- NAO E PROVA.

Custo de S3b: certificar que, num disco D em volta de ~0,401, L tem uma raiz
algebricamente simples. Duas rotas, medidas pela grandeza que as dimensiona:
  (i)  winding 1 num circulo |s-s0|=r:  max ||V_L(s)|| no circulo
  (ii) Newton-Kantorovich no sistema aumentado F(x,s) = (H(s)x, l(x)-1),
       H(s) = (A+s)J^-1, x = J h.  DF invertivel <=> raiz algebricamente simples
       (Keldysh); o custo e ||DF^-1|| nos mesmos pesos.
Referencia: max ||V_L|| em Gamma_A (1846 em 12x36).
"""
import sys
import numpy as np
sys.path.insert(0, 'scripts')
from measure_sharp_vs_selected_V import selecionado, norma_V


def raiz_perto(A, alvo):
    ev = -np.linalg.eigvals(A)
    return ev[np.argmin(np.abs(ev - alvo))].real


def bordejado(A, J, w, s0):
    n = len(A)
    Jinv = np.linalg.inv(J)
    H = (A + s0*np.eye(n)) @ Jinv
    _, _, Vt = np.linalg.svd(H)
    x0 = Vt[-1].conj()
    x0 = x0 / np.sum(np.abs(x0)*w)                 # ||x0||_w = 1
    k = int(np.argmax(np.abs(x0)*w))
    DF = np.zeros((n+1, n+1), dtype=complex)
    DF[:n, :n] = H
    DF[:n, n] = Jinv @ x0                          # d/ds [(A+s)J^-1 x0] = J^-1 x0
    DF[n, k] = 1/x0[k]                             # l(x) = x_k / x0_k
    Dinv = np.linalg.inv(DF)
    W = np.append(w, 1.0)
    return float(np.max((np.abs(Dinv)*W[:, None]).sum(axis=0)/W))


def main(malhas):
    for malha in malhas:
        A, J, w = selecionado(malha)
        s0 = raiz_perto(A, 0.401)
        print(f'\n=== {malha}: raiz s0 = {s0:.6f}')
        for r in (0.05, 0.10, 0.15, 0.20):
            pts = [s0 + r*np.exp(1j*t) for t in np.linspace(0, np.pi, 9)]
            print(f'   circulo r={r:.2f}: max ||V_L|| = {max(norma_V(A, J, w, s) for s in pts):9.1f}')
        print(f'   Newton-Kantorovich: ||DF^-1|| = {bordejado(A, J, w, s0):9.1f}')


if __name__ == '__main__':
    main(sys.argv[1:] or ['8x24', '12x36'])
