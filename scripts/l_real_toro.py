#!/usr/bin/env python3
"""L-real, lado do toro (T2_HIPOTESES.md §3.1): contração da cauda T(H - I)T no l2 do toro com sinal.

Recuperação de raio (como em RADIUS_RECOVERY.md) entre o toro (pesos kappa1^m sqrt(omega2(n)), o espaço de
todos os certificados de contagem) e a realização l1 de RT (Y+, pesos (65/64, 5/4)), com Q = J_mu^-1 (o mesmo
nos dois espaços, independente de s), H(s) = I + (B + s) Q e a MESMA caixa G = {|m| < M, n < N}:
  ||T (H - I) T||_toro <= (||B_L||_toro + |s|) ||J_mu^-1 T||_toro,   T = I - G (projeção ortogonal no toro).
J_mu^-1 é diagonal em m (o peso kappa1^m cancela dentro do modo):
  ||J^-1 T|| = max( max_{|m| < m0} cauda_Rinv(i m/2, N) [colunas n >= N do modo m],
                    max_{|m| >= m0} ||R_m^-1|| <= (1 + 2 c_w/(k2 - 1))/(|m|/2)   [Rinv_uniforme, s = 0] ),
com os modos m0 <= |m| < M cobertos pela cota do modo inteiro (que domina a das colunas n >= N).
||B_L||_toro <= 3,9641: símbolo de F3 em Arb (ALTERNATIVE_ROUTES.md §5.26); mais o erro do fundo de RT.
mu: a cota é calculada com mu_RefA -+ 1e-10 (|mu* - mu_RefA| <= 3,9e-11)."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent))
import closed_form as cf
import nk_L as nk

K2 = nk.K2


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--M', type=int, default=1024, help='so documenta a caixa: a cota de |m| >= m0 cobre M (revisao T2, item 9)')
    ap.add_argument('--N', type=int, default=4096)
    ap.add_argument('--m0', type=int, default=400)
    ap.add_argument('--S', type=float, default=3.1, help='cota de |s| na região (|3 + 0,75i| < 3,1)')
    ap.add_argument('--nB', type=float, default=3.9642, help='||B_L||_toro: a cota certificada do simbolo e 3,9641333 (symbol-L-2048x4096.json), arredondada para cima; antes 3,964101, abaixo dela (revisao do artigo, 08/10, I-3). A contracao vale para ||B|| ate ~15,9')
    ap.add_argument('--saida', type=Path, default=Path('build/t2/l_real_toro.json'))
    a = ap.parse_args()
    t0 = time.time()
    pior = 0.0; onde = None
    for mu in (nk.MU - 1e-10, nk.MU + 1e-10):
        for m in range(-(a.m0 - 1), a.m0):
            for passo in ((2, 1) if m % 2 == 0 else (1,)):
                v = cf.cauda_Rinv(1j*m/2, a.N, passo, mu, K2, a.N + 2500)
                if v > pior:
                    pior, onde = v, (m, passo, mu)
        print(f'mu {mu:.12f}: max cauda (|m| < {a.m0}) = {pior:.6f} em {onde}  [{time.time() - t0:.0f}s]', flush=True)
    cw = math.sqrt(1 + K2**-4)
    unif = (1 + 2*cw/(K2 - 1))/(a.m0/2)*1.0001          # |m| >= m0, s = 0
    qT = max(pior, unif)
    defeito = (a.nB + a.S)*qT
    print(f'||J^-1 T||_toro <= max({pior:.6f}, {unif:.6f}) = {qT:.6f};  ||T(H - I)T||_toro <= ({a.nB:.4f} + {a.S}) * {qT:.6f}'
          f' = {defeito:.4f}  ' + ('CONTRAI' if defeito < 1 else 'NAO contrai'))
    a.saida.parent.mkdir(parents=True, exist_ok=True)
    a.saida.write_text(json.dumps(dict(M=a.M, N=a.N, m0=a.m0, S=a.S, nB=a.nB, cauda_max=pior, onde=str(onde),
                                       uniforme=unif, qT=qT, defeito=defeito, contrai=defeito < 1), indent=1) + '\n')


if __name__ == '__main__':
    main()
