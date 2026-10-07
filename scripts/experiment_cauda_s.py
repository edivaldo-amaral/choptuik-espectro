#!/usr/bin/env python3
"""Diagnóstico (ponto flutuante) para o desenho do Rouché de L: quanto o inverso de uma janela da
cauda num centro c0 aproxima o bloco em s = c0 + delta. rho(delta) = ||I - V_J(c0) Hh_J(s)||
= ||V_J(c0) (Hh_J(c0) - Hh_J(s))||, Hh_J(s) = I + P_J B Q0(s) P_J (janelas de nk_L3). Define o
espaçamento dos centros de cauda."""
from __future__ import annotations
import argparse, sys, time
from pathlib import Path
import numpy as np
from scipy.linalg import lu_factor
sys.path.insert(0, str(Path(__file__).resolve().parent))
import nk_L3
from nk_L3 import potencia_lu

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--c0', type=complex, default=0.3+0.25j)
    ap.add_argument('--deltas', default='0.05,0.1,0.2,0.4,0.1j,0.25j,-0.25j')
    ap.add_argument('--janelas', default='32..47,-15..0,144..159')
    ap.add_argument('--Z', default='32x128'); ap.add_argument('--F', default='200x400'); ap.add_argument('--W', type=int, default=16)
    a = ap.parse_args(); t0 = time.time()
    mk = lambda s: nk_L3.Cert3(argparse.Namespace(Z=a.Z, F=a.F, W=a.W, lam0=s))
    base = mk(a.c0)
    outros = {d: mk(a.c0 + complex(d)) for d in a.deltas.split(',')}
    for nome in a.janelas.split(','):
        i = [lv['nome'] for lv in base.niveis].index(nome)
        H0 = base.H(base.niveis[i]); lu = lu_factor(H0.copy(), check_finite=False)
        print(f'janela {nome} ({H0.shape[0]}): ||V_J(c0)|| {potencia_lu(lu, np.eye(H0.shape[0])):.4f}   [{time.time()-t0:.0f}s]', flush=True)
        for d, c in outros.items():
            X = H0 - c.H(c.niveis[i])
            print(f'   delta {d:>7}: rho = ||V_J(c0)(Hh(c0) - Hh(s))|| = {potencia_lu(lu, X):.4f}   [{time.time()-t0:.0f}s]', flush=True)
            del X; c.ctx.Q.cache_clear(); c.ctx._B.cache_clear()
        del H0, lu; base.ctx.Q.cache_clear(); base.ctx._B.cache_clear()

if __name__ == '__main__':
    main()
