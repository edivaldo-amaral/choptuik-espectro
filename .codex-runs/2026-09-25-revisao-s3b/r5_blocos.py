#!/usr/bin/env python3
"""R5 (revisao independente de S3b): normas VERDADEIRAS (ponto flutuante, montagem propria) de
blocos de E = I - A DF(x0) do certificado NK de L, para comparar com as cotas de T3.json.

Montagem propria: pesos kappa1^m sqrt(omega2(n)) calculados aqui; Q0 = inverso de (J_mu + lambda0)
por modo (J de signed_operator_L.J_livre, sem pesos, invertida aqui e pesada aqui); B de
signed_operator_L.bloco_B (sem pesos; banda COMPLETA |d| <= 40 para o operador 'verdadeiro' e
|d| <= 16 para Hh, o inverso aproximado do certificado); selecao das colunas de T (n >= NZ se
|m| < MZ; tudo se |m| >= MZ) escrita aqui. Normas por iteracao de potencia com a LU de Hh
(cotas INFERIORES da norma verdadeira; o certificado precisa estar ACIMA delas).

Uso: python r5_blocos.py I [J1 J2 ...] [--TZ] [--alto]   (I, J: nomes 'a..b' das janelas)"""
import argparse
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
from scipy.linalg import lu_factor, lu_solve

sys.path.insert(0, str(Path('scripts').resolve()))
import signed_operator_L as sl

K1, K2, MU = 65/64, 5/4, 722873400/2**32
LAM0 = 0.4010247330
MZ, NZ, MC, NC = 32, 128, 200, 400
T0 = time.time()


def log(s):
    print(f'{s}   [{time.time()-T0:.0f}s]', flush=True)


def omega(n):
    n = np.asarray(n, float)
    return np.where(n == 0, 1.0, K2**(2*n) + K2**(-2*n))


def layout(m, N):
    """(componente, n) por linha do modo m truncado em n < N (a ordem de sl.Radial)."""
    r = sl.Radial(m, N)
    comp = np.concatenate([np.full(len(ns), c) for c, ns in r.comps])
    nn = np.concatenate([ns for _, ns in r.comps])
    return r, comp, nn


def peso(m, N):
    _, _, nn = layout(m, N)
    return K1**m*np.sqrt(omega(nn))


def selT(m, N=NC):
    _, _, nn = layout(m, N)
    return np.where(nn >= NZ)[0] if abs(m) < MZ else np.arange(len(nn))


def selZ(m):
    _, _, nn = layout(m, NC)
    return np.where(nn < NZ)[0]


G = sl.campos()
_Qc, _Bc = {}, {}


def Q0(m):
    """Q0 pesado do modo m (n < NC)."""
    if m not in _Qc:
        if len(_Qc) > 60:
            _Qc.clear()
        r = sl.Radial(m, NC)
        J = sl.J_livre(r, LAM0)
        w = peso(m, NC)
        Qm = np.linalg.solve(J, np.eye(len(J)))
        _Qc[m] = (w[:, None]*Qm)/w[None, :]
    return _Qc[m]


def Bw(mo, mi, No=NC, Ni=NC):
    """Bloco pesado de B, modo mo (n < No) <- modo mi (n < Ni)."""
    k = (mo - mi, mi % 2, No, Ni)
    if k not in _Bc:
        if len(_Bc) > 200:
            _Bc.clear()
        ro, ri = sl.Radial(mo, No), sl.Radial(mi, Ni)
        X = sl.bloco_B(G, ro, ri)
        # pesos: kappa1^mo/kappa1^mi = kappa1^d; radiais
        wo = peso(mo, No)/K1**mo; wi = peso(mi, Ni)/K1**mi
        _Bc[k] = (K1**(mo - mi))*(wo[:, None]*X)/wi[None, :]
    return _Bc[k]


def janela(nome):
    a, b = nome.split('..')
    return list(range(int(a), int(b) + 1))


def bloco(lin_modos, lin_sel, col_modos, col_sel, banda):
    """P_lin B Q0 P_col (banda |d| <= banda), denso."""
    ol = np.cumsum([0] + [len(lin_sel(m)) for m in lin_modos])
    oc = np.cumsum([0] + [len(col_sel(m)) for m in col_modos])
    X = np.zeros((int(ol[-1]), int(oc[-1])), complex)
    for q, mi in enumerate(col_modos):
        cs = col_sel(mi)
        if not len(cs):
            continue
        Qi = Q0(mi)[:, cs]
        for p, mo in enumerate(lin_modos):
            if abs(mo - mi) <= banda:
                ls = lin_sel(mo)
                X[ol[p]:ol[p+1], oc[q]:oc[q+1]] = Bw(mo, mi)[ls] @ Qi
    return X


def potencia(apl, aplH, n, it=50, seed=0):
    rng = np.random.default_rng(seed)
    v = rng.standard_normal(n) + 1j*rng.standard_normal(n); v /= np.linalg.norm(v)
    s = 0.0
    for _ in range(it):
        u = apl(v); s = np.linalg.norm(u)
        v = aplH(u); nv = np.linalg.norm(v)
        if nv == 0:
            return 0.0
        v /= nv
    return float(np.linalg.norm(apl(v)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('I'); ap.add_argument('J', nargs='*')
    ap.add_argument('--TZ', action='store_true'); ap.add_argument('--alto', action='store_true')
    ap.add_argument('--saida', default=None)
    a = ap.parse_args()
    t3 = json.loads(Path('T3.json').read_text())
    nomes = t3['nome']; iI = nomes.index(a.I)
    mI = janela(a.I)
    out = dict(I=a.I)
    Hh = bloco(mI, selT, mI, selT, 16); Hh[np.diag_indices(len(Hh))] += 1.0
    n = len(Hh)
    log(f'janela {a.I}: dim {n}')
    # H verdadeiro - Hh = banda 17..40 (o resto)
    D = bloco(mI, selT, mI, selT, 40); D[np.diag_indices(n)] += 1.0; D -= Hh
    lu = lu_factor(Hh, overwrite_a=True, check_finite=False); del Hh
    sol = lambda x: lu_solve(lu, x, check_finite=False)
    solH = lambda x: lu_solve(lu, x, trans=2, check_finite=False)
    nV = potencia(sol, solH, n)
    out['nV'] = (nV, t3['nV'][iI])
    log(f'  ||Hh^-1||: verdadeiro ~ {nV:.5f}   cota T3 {t3["nV"][iI]:.5f}')
    rho = potencia(lambda v: sol(D @ v), lambda u: D.conj().T @ solH(u), n)
    out['rho'] = (rho, t3['N'][iI][iI])
    log(f'  ||I - Hh^-1 H_II|| (banda 17..40): ~ {rho:.3e}   cota T3 (diagonal) {t3["N"][iI][iI]:.3e};  ||resto_II|| ~ {np.linalg.norm(D, 2) if n < 3000 else float("nan"):.2e}')
    del D
    if a.alto:
        nn = np.concatenate([layout(m, NC)[2][selT(m)] for m in mI])
        alto = nn >= NC - 100 - 60
        P = lambda v: np.where(alto, v, 0)
        va = potencia(lambda v: sol(P(v)), lambda u: P(solH(u)), n)
        out['nVhi'] = (va, t3['nVhi'][iI])
        log(f'  ||Hh^-1 P_alto||: ~ {va:.5f}   cota T3 {t3["nVhi"][iI]:.5f}')
    for J in a.J:
        iJ = nomes.index(J); mJ = janela(J)
        X = bloco(mI, selT, mJ, selT, 40)
        v = potencia(lambda x: sol(X @ x), lambda u: X.conj().T @ solH(u), X.shape[1])
        out[f'N[{a.I}<-{J}]'] = (v, t3['N'][iI][iJ])
        log(f'  ||Hh_I^-1 P_I B Q0 P_J||, J = {J}: ~ {v:.5f}   cota T3 {t3["N"][iI][iJ]:.5f}')
        del X
    if a.TZ:
        mZ = list(range(-(MZ - 1), MZ))
        X = bloco(mI, selT, mZ, selZ, 40)
        v = potencia(lambda x: sol(X @ x), lambda u: X.conj().T @ solH(u), X.shape[1])
        out['TZ'] = (v, t3['TZ'][iI])
        log(f'  ||Hh_I^-1 P_I B Q0 P_Z||: ~ {v:.5f}   cota T3 {t3["TZ"][iI]:.5f}')
        del X
    if a.saida:
        Path(a.saida).write_text(json.dumps(out, indent=1) + '\n')


if __name__ == '__main__':
    main()
