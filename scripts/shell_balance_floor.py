#!/usr/bin/env python3
"""O PISO exato da norma de coluna sobre qualquer familia de pesos diagonais.

A rota k=0 (cota analitica estilo RT, sem trabalho por coluna) e a unica que
tornaria a arquitetura (c) viavel nesta maquina. Ela exige colunas de
X=(B_true+s)Q_ref com norma <=~0,6 perto do nucleo; hoje chegam a ~58.

A norma usada no projeto e colnorm(M) = max_j sum_i (|Re M_ij| + |Im M_ij|),
ou seja a norma 1 da matriz NAO NEGATIVA  Bm = |Re X| + |Im X|.

Mudar os pesos eta_d, (65/64)^m, (5/4)^n -- ou qualquer outra familia de pesos
por DOF -- e exatamente uma similaridade diagonal positiva X -> D^-1 X D, que
leva Bm -> D^-1 Bm D entrada a entrada. Pelo teorema de Perron-Frobenius,

    min_{D>0 diagonal} || D^-1 Bm D ||_1  =  rho(Bm),

com o minimo atingido no autovetor de Perron. Logo rho(Bm) e um PISO que
nenhuma escolha de pesos pode furar. Se rho(Bm) > 0,6, a rota k=0 esta
fechada para esta discretizacao, em qualquer peso por DOF.

DIAGNOSTICO em float; o piso e um fato de algebra linear, nao uma prova
certificada (os autovalores sao calculados em ponto flutuante).
"""
from __future__ import annotations

import json

import numpy as np

import g_parametrix_prototype as G

CENTERS = (0.125+0.25j, 1.0, 0.733)


def perron(Bm):
    """raiz de Perron de Bm>=0, pelo espectro denso (1422x1422 e barato e confiavel;
    a iteracao de potencia nao converge aqui, o espectro tem topo quase degenerado).
    Devolve rho e o autovetor de Perron (nao negativo)."""
    lam, V = np.linalg.eig(Bm)
    i = int(np.argmax(np.abs(lam)))
    rho = float(np.abs(lam[i]))
    v = np.abs(np.real_if_close(V[:, i]))
    if v.sum() == 0:
        v = np.ones(len(Bm))
    v = v/v.sum()
    # certificado barato de Collatz-Wielandt no autovetor (so onde v>0)
    pos = v > 1e-14
    ratio = (Bm@v)[pos]/v[pos]
    return dict(rho=rho, cw_lower=float(np.min(ratio)), cw_upper=float(np.max(ratio)),
                support=int(pos.sum()), dim=len(Bm)), v


def main():
    out = dict(note='DIAGNOSTICO float; piso de Perron para a norma de coluna', centers={})
    A, addr = G.read(G.ROOT/'build/spectrum/rt-A-12x36.dat')
    n = len(addr)
    Qm = np.linalg.inv(G.free(addr))
    w = np.array([G.ETA[d]*(2-(m == 0))*(2-(nn == 0))*(65/64)**m*(5/4)**nn for d, p, m, nn in addr])
    for s0 in CENTERS:
        H = (w[:, None]*((A+s0*np.eye(n))@Qm)/w[None, :]).astype(complex)
        X = H-np.eye(n)
        Bm = np.abs(X.real)+np.abs(X.imag)
        atual = float(np.max(Bm.sum(axis=0)))
        cw, v = perron(Bm)
        rho = cw['rho']
        # o peso otimo é o autovetor de Perron; ele é da forma a^m b^n?
        lv = np.log(np.maximum(v, 1e-300))
        M = np.array([[1.0, a[2], a[3]] for a in addr])          # 1, m, n
        coef, *_ = np.linalg.lstsq(M, lv, rcond=None)
        resid = lv-M@coef
        ajuste = dict(log_a_por_m=float(coef[1]), log_b_por_n=float(coef[2]),
                      r2=float(1-np.var(resid)/np.var(lv)))
        # norma atingida por um peso PURAMENTE da forma a^m b^n (melhor ajuste)
        vfit = np.exp(M@coef)
        Bfit = (1/vfit)[:, None]*Bm*vfit[None, :]
        colfit = float(np.max(Bfit.sum(axis=0)))
        out['centers'][str(s0)] = dict(
            colnorm_atual=atual, piso_perron=rho, collatz_wielandt=cw, alvo_k0=0.6,
            razao_piso_alvo=rho/0.6,
            colnorm_com_peso_exponencial=colfit,
            ajuste_do_peso_otimo=ajuste,
            ganho_maximo_de_repesagem=atual/rho if rho > 0 else None)
        print(f's={s0}: colnorm={atual:.3f}  piso de Perron={rho:.4f} '
              f'[CW {cw["cw_lower"]:.3f}, {cw["cw_upper"]:.3f}; suporte {cw["support"]}/{cw["dim"]}]  '
              f'alvo k=0: 0,6  ->  ainda faltam {rho/0.6:.1f}x  (repesagem ganha {atual/rho:.1f}x)', flush=True)
    p = G.ROOT/'build/spectrum/shell-balance-floor.json'
    p.write_text(json.dumps(out, indent=2)+'\n')
    print('escrito:', p)


if __name__ == '__main__':
    main()
