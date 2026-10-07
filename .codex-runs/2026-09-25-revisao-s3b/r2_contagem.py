#!/usr/bin/env python3
"""R2 (revisao independente de S3b): autovalores de Kc_Z (Z = 20x80, fundo truncado |d| <= 20,
n <= 40) no disco |s - 0,401| < 0,05, com uma montagem INDEPENDENTE da do certificado:
sharp_operator_builder (base REAL (cos, sen), convolucao 2D nos campos Z^2, validada contra o
exportador de RT), com o truncamento do fundo feito aqui nos arrays Z^2. Kc_Z(s) = K + s, entao
as raizes sao s = -autovalores de K. Tambem: a mesma conta com o fundo completo e em Z maiores
(estabilidade da contagem), e a separacao do autovalor (condicao e distancia ao proximo)."""
import sys
from pathlib import Path

import numpy as np
import scipy.linalg as sla

sys.path.insert(0, str(Path('scripts').resolve()))
import sharp_operator_builder as sb

C, RHO = 0.401, 0.05
_orig = sb.campos_z2


def truncado(Dm=20, Dn=40):
    def f():
        out = []
        for G in _orig():
            H = G.copy()
            m = np.arange(-sb.MB, sb.MB + 1)[:, None]; n = np.arange(-sb.NB, sb.NB + 1)[None, :]
            H[(np.abs(m) > Dm) | (np.abs(n) > Dn)] = 0
            out.append(H)
        return out
    return f


def conta(M, N, trunc, rot):
    sb.campos_z2 = truncado() if trunc else _orig
    K, addr = sb.operador(M, N)
    w, vl, vr = sla.eig(K, left=True, right=True)
    s = -w
    d = np.abs(s - C)
    dentro = np.where(d < RHO)[0]
    ordem = np.argsort(d)
    txt = ', '.join(f'{s[i].real:.10f}{s[i].imag:+.1e}j' for i in dentro)
    cond = [1/abs(vl[:, i].conj() @ vr[:, i]) for i in dentro]
    prox = d[ordem[len(dentro)]] if len(ordem) > len(dentro) else float('nan')
    print(f'{rot:34s} dim {len(K):5d}: {len(dentro)} raiz(es) no disco: {txt};  cond. do autovalor '
          + ', '.join(f'{c:.2f}' for c in cond) + f';  proxima raiz a {prox:.4f} do centro (raio 0,05)', flush=True)
    return len(dentro)


if __name__ == '__main__':
    conta(20, 80, True, 'Kc_Z 20x80, fundo truncado (20,40)')
    conta(20, 80, False, 'K 20x80, fundo completo')
    conta(24, 96, True, 'Kc 24x96, fundo truncado')
    conta(16, 64, True, 'Kc 16x64, fundo truncado')
