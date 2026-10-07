#!/usr/bin/env python3
"""H-rec (ida), passo R5: recuperacao de raio a partir de QUALQUER peso fraco (k1', k2') com 1 <= k1' <= 65/64 (peso de Fourier 1 incluido),
1 < k2' <= 9/8, para |s| <= S com S QUALQUER (HREC_IDA.md; a caixa G cresce linearmente com S).

Mesmo argumento de RADIUS_RECOVERY.md (existencia forte, unicidade fraca, cadeias), com Q = J_mu^-1 e
H(s) = I + (B + s) Q, T = I - G, G = {|m| < M, n < N}:
  fraco:  ||T(H - I)T||_k' <= (beta(k2') + S) * free_tail(M, N, k2')
  forte:  ||T(H - I)T||_+  <= (beta+ + S) * structured_input_tail(M, N)          (Y+, pesos (65/64, 5/4))
beta(k2') = beta+ * D(k2')/D(5/4), D(k) = 2k/(k^2 - 1) a norma da divisao regularizada por xi (40/9 em 5/4,
144/17 em 9/8); as multiplicacoes nao aumentam a norma com pesos menores. free_tail nao depende de k1'
(Q preserva m). Cota explicita (papel, HREC_IDA.md): com r = 1/k2', b = M/2 >= 1, n >= 2,
  free_tail <= max( (3/2)(1 + 2 mu_max)/(b (1 - r)),  (3/2)(1 + 4/(1 - r))/(mu_min (N + 1)) ),
que tende a 0; logo para todo k2' > 1 existe G que contrai."""
from __future__ import annotations

from fractions import Fraction as Q
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from sharp_exterior import free_tail
from structured_exterior import structured_input_tail

BETA_P = Q(5191, 500)
MU_MIN, MU_MAX = Q(1, 6), Q(17, 100)


def D(k):
    return 2*k/(k*k - 1)


def beta(k2):
    return BETA_P*D(k2)/D(Q(5, 4))


def cota_explicita(M, N, k2):
    r = 1/k2; b = Q(M, 2)
    return max(Q(3, 2)*(1 + 2*MU_MAX)/(b*(1 - r)), Q(3, 2)*(1 + 4/(1 - r))/(MU_MIN*(N + 1)))


def caixa(k2, S):
    """Menores M, N (potencias de 2) pela cota explicita; confere free_tail <= cota e as duas contracoes."""
    bw = beta(k2)
    M = 2
    while (bw + S)*Q(3, 2)*(1 + 2*MU_MAX)/(Q(M, 2)*(1 - 1/k2)) >= Q(9, 10):
        M *= 2
    N = 4
    while (bw + S)*Q(3, 2)*(1 + 4/(1 - 1/k2))/(MU_MIN*(N + 1)) >= Q(9, 10):
        N *= 2
    ft = free_tail(M, N, radial=k2, mu_min=MU_MIN, mu_max=MU_MAX)
    assert ft <= cota_explicita(M, N, k2), (k2, M, N)
    fraco = (bw + S)*ft
    forte = (BETA_P + S)*structured_input_tail(M, N)
    return bw, M, N, fraco, forte


def main():
    linhas = []
    for S in (Q(31, 10), Q(10), Q(50), Q(200)):
        for k in ((3, 4, 5, 6, 8, 10, 12) if S == Q(31, 10) else (3, 6, 10)):
            k2 = 1 + Q(1, 2**k)
            bw, M, N, fraco, forte = caixa(k2, S)
            ok = fraco < 1 and forte < 1
            linhas.append(dict(S=str(S), k2=str(k2), beta=float(bw), M=M, N=N, fraco=float(fraco), forte=float(forte),
                               free_tail_le_cota=True, contrai=ok))
            print(f'S = {float(S):6.1f}, k2\' = 1 + 2^-{k:<2}: beta = {float(bw):9.2f}, G = {M} x {N}: '
                  f'fraco {float(fraco):.4f}, forte {float(forte):.4f}  ' + ('OK' if ok else 'FALHA'))
    Path('build/t2').mkdir(parents=True, exist_ok=True)
    Path('build/t2/radius_recovery_geral.json').write_text(json.dumps(dict(linhas=linhas), indent=1) + '\n')


if __name__ == '__main__':
    main()
