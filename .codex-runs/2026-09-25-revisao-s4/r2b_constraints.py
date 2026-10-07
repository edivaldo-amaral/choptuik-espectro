#!/usr/bin/env python3
"""R2b (revisao S4): D g* = 0 via a identidade C_w(mu)(Z + 1)w = (Z + 1)c_w nos operadores do codigo
(signed_constraints.bloco_C), com c_w = (C_w(0) + C_0(0)) w / 2 (constraints quadraticas, sem termo
constante: conferido por c(RefA) ~ 0). Controles: (Z + 2) e d_tau -> -i m/2."""
from __future__ import annotations
import math, sys
from pathlib import Path
import numpy as np

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI.parents[1]/'scripts')); sys.path.insert(0, str(AQUI))
import signed_operator_L as sl
import signed_constraints as SC
import r2_identidade as R2

out = []
def P(s=''):
    print(s, flush=True); out.append(s)


def aplica_C(gbg, v, Nin, Nout, Mout, s, mu):
    g = [dict() for _ in range(4)] if gbg is None else gbg
    res = {}
    for mo in sl.modos(Mout):
        if mo % 2:
            continue
        so = SC.Saida(mo, Nout); acc = np.zeros(so.dim, complex)
        for mi, x in v.items():
            if abs(mo - mi) <= sl.BAND_M and np.any(x):
                acc += SC.bloco_C(g, so, sl.Radial(mi, Nin), s, mu) @ x
        res[mo] = acc
    return res


def Zsaida(v, N, mu, k, sinal=1):
    w = {}
    for m, x in v.items():
        so = SC.Saida(m, N); y = np.zeros_like(x)
        for c, ns in so.comps:
            o0, o1 = so.off[c]; full = np.zeros(N, complex); full[ns] = x[o0:o1]
            y[o0:o1] = ((k + sinal*1j*m/(2*mu))*full + R2.xidxi(full))[ns]
        w[m] = y
    return w


def nrm(v):
    return math.sqrt(sum(float(np.sum(np.abs(x)**2)) for x in v.values()))


def main():
    rng = np.random.default_rng(11)
    mu = 0.37; M0, N0 = 4, 9; Mo, No = 2*M0 + 2, 2*N0 + 8
    gw = R2.aleatorio(M0, N0, rng)
    w = R2.campos_para_modos(gw, M0 + 1, N0 + 1)
    cw = R2.comb((0.5, aplica_C(gw, w, N0 + 1, No, Mo, 0.0, mu)), (0.5, aplica_C(None, w, N0 + 1, No, Mo, 0.0, mu)))
    for sinal, k, nome in ((1, 1, 'd_tau -> +i m/2, (Z + 1)c_w'), (-1, 1, 'controle d_tau -> -i m/2'), (1, 2, 'controle (Z + 2)c_w')):
        lhs = aplica_C(gw, R2.aplica_Z(w, N0 + 1, mu, 1, sinal), N0 + 1, No, Mo, mu, mu)
        rhs = Zsaida(cw, No, mu, k, sinal)
        dif = math.sqrt(sum(float(np.sum(np.abs(lhs[m] - rhs[m])**2)) for m in lhs))
        P(f'(T5) w aleatorio, mu = {mu}: {nome}: ||C_w(mu)(Z+1)w - rhs|| / ||C_w(mu)(Z+1)w|| = {dif/nrm(lhs):.2e}')
    g = sl.campos(); Ng = 101
    wA = R2.campos_para_modos(g, 41, Ng)
    cA = R2.comb((0.5, aplica_C(g, wA, Ng, Ng + 104, 82, 0.0, sl.MU)), (0.5, aplica_C(None, wA, Ng, Ng + 104, 82, 0.0, sl.MU)))
    P(f'(T5) c(RefA) (constraints no fundo RefA, sem pesos) = {nrm(cA):.3e}; ||RefA|| = {nrm(wA):.3e}')
    lhs = aplica_C(g, R2.aplica_Z(wA, Ng, sl.MU, 1), Ng, Ng + 104, 82, sl.MU, sl.MU)
    P(f'(T5) ||C_A(mu_A) g_A|| = {nrm(lhs):.3e} (sem pesos); ||(Z_A + 1)c(RefA)|| = {nrm(Zsaida(cA, Ng + 104, sl.MU, 1)):.3e}')
    (AQUI/'r2b_saida.txt').write_text('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
