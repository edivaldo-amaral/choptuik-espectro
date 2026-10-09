#!/usr/bin/env python3
"""N3: o resíduo B h0 é calculado com linhas truncadas em n < NZ + BAND_N + 1 = 229, e os blocos far <- janela
em n < Nc + BAND_N + 1 = 501. Isso só é exato se B leva a coluna n_in a linhas n_out <= n_in + 101.
Checa nos blocos montados (ctx._B) com linhas até n < N_in + 300: maior (n_out - n_in) com entrada não nula."""
import sys
sys.path.insert(0, 'scripts')
import numpy as np
import nk_L as nk
ctx = nk.Contexto(0.7331796596606924, 400)
pior = -10**9
for Ni in (128, 400):
    for d in (0, 1, 2, 5, 16, 17, 39, 40, -1, -16, -40):
        for par in (0, 1):
            No = Ni + 300
            X = ctx._B(d, par, No, Ni)
            ro = nk.sl.Radial(par + d, No); ri = nk.sl.Radial(par, Ni)
            no = np.concatenate([ns for _, ns in ro.comps]); ni = np.concatenate([ns for _, ns in ri.comps])
            nz = np.abs(X) > 0
            if nz.any():
                ii, jj = np.nonzero(nz)
                s = int((no[ii] - ni[jj]).max())
                pior = max(pior, s)
                print(f'Ni {Ni} d {d:4d} par {par}: max(n_out - n_in) = {s}; maior n_out = {no[ii].max()}', flush=True)
        ctx._B.cache_clear()
print('MAIOR SALTO n_out - n_in =', pior, ' (o código supõe <= BAND_N + 1 =', nk.BAND_N + 1, ')')
