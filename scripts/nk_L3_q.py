#!/usr/bin/env python3
"""Normas RESTRITAS de Q0 para o termo de 2a ordem do NK de L (nk_L3_C.py) e o a_extra de Z'.

Q0 e diagonal em m e triangular superior em n. Em DF(x) - DF(x0) = [(lambda - lambda0) Q0 dy +
dlambda Q0 (y - y0); 0], o que entra por nivel e:
  q_ZT  = ||P_Z Q0 P_T|| = max sobre os modos de Z de ||R_m^-1 [n < NZ, n >= NZ]|| (a descida das
          colunas de cauda para dentro de Z),
  q_TT  = ||P_T Q0 P_T|| = max sobre os modos de F do bloco [linhas de T, colunas de T],
  q_Zf  = ||P_Z Q0 P_far||: descida de n >= Nc ate n < NZ (forma fechada),
  q_Tf  = ||P_T Q0 P_far|| <= ||Q0 P_far|| = Qfar.
E o maximo de ||R_m^-1|| nos modos MZ + banda <= |m| < MZ + 40 (a_extra de Z', lacuna L5 da
revisao). Todas as normas por cota_norma (Rump) com o erro de Q0 (Q_cert)."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import closed_form as cf
import nk_L as nk
import signed_operator_L as sl
from nk_L2_Z import Q_cert, FATOR_PESO
from verified_norms import cota_norma


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--dir', type=Path, default=Path('build/s3b/rig'))
    ap.add_argument('--banda', type=int, default=16)
    a = ap.parse_args()
    z = json.loads((a.dir/'Z.json').read_text())
    MZ, NZ = z['Z']; Mc, Nc = z['F']; lam0 = z['lam0']
    lam0 = complex(*lam0) if isinstance(lam0, list) else lam0      # faixa B: lambda0 complexo
    ctx = nk.Contexto(lam0, Nc)
    qZT = 0.0; qTT = 0.0; nQx = 0.0
    for m in range(-(Mc - 1), Mc):
        Q, dq, nq = Q_cert(m, Nc, lam0)
        nn = np.concatenate([ns for _, ns in ctx.rad(m, Nc).comps])
        colT = nk.colsT(ctx, m, MZ, NZ, Nc)
        if abs(m) < MZ:
            rZ = np.where(nn < NZ)[0]
            qZT = max(qZT, cota_norma(Q[np.ix_(rZ, colT)], dq))
        qTT = max(qTT, cota_norma(Q[np.ix_(colT, colT)], dq))
        if MZ + a.banda <= abs(m) < MZ + sl.BAND_M:
            nQx = max(nQx, nq + dq)
        ctx._rad.clear()
    qZf = cf.descida(NZ, Nc, sl.MU, sl.K2)
    out = dict(qZT=qZT*FATOR_PESO, qTT=qTT*FATOR_PESO, qZf=qZf, nQx_extra=nQx*FATOR_PESO, lam0=([lam0.real, lam0.imag] if isinstance(lam0, complex) else lam0))
    (a.dir/'Q.json').write_text(json.dumps(out, indent=1) + '\n')
    print(f'q_ZT = ||P_Z Q0 P_T|| <= {out["qZT"]:.5f};  q_TT = ||P_T Q0 P_T|| <= {out["qTT"]:.5f};'
          f'  q_Zf <= {qZf:.1e};  max ||Q0|| em {MZ + a.banda} <= |m| < {MZ + sl.BAND_M}: {out["nQx_extra"]:.5f}')


if __name__ == '__main__':
    main()
