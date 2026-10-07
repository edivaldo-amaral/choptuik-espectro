#!/usr/bin/env python3
"""DIAGNOSTICO em float (nao e prova): qual GRAU k cada pre-condicionador de casca exige.

Pergunta (rodada 6, 17/09/2026). O custo da arquitetura (c) e k=32 aplicacoes por
coluna da casca. A pergunta decisiva nao e "quanto vale o defeito", e sim:

    k_min(P) = menor k com ||(I - H_OO P)^k||_coluna <= alvo,

porque o custo por regiao e proporcional a k. Um pre-condicionador so serve se
baixar k_min, e o alvo tem de ser medido na NORMA DE COLUNA (o que se certifica),
nao no raio espectral (que nao se certifica e ignora a nao-normalidade).

Candidatos: Neumann amortecida (a rota atual), Jacobi, bloco-m, bloco-n, e
deflacao de posto r pelo projetor ESPECTRAL (autovetores a esquerda, inversa
exata no subespaco -- nao pinv, que so projeta ortogonalmente).

Cada candidato leva seu proprio amortecimento otimo, escolhido pelo raio
espectral (barato, via autovalores) e depois conferido na norma de coluna.
"""
from __future__ import annotations

import json

import numpy as np

import g_parametrix_prototype as G

CENTERS = (0.125+0.25j, 1.0, 0.733)
CORES = ((6, 18), (8, 24))
KMAX = 72
# 0,9 e o criterio pre-registrado; 1e-3 e 1e-5 sao o ponto de OPERACAO real,
# porque o raio do disco de validade uniforme escala com (1-alpha): baixar o alvo
# encolhe a regiao e multiplica o numero de regioes.
TARGETS = (0.9, 0.5, 1e-3, 1e-5)


def colnorm(M):
    return float(np.max(np.sum(np.abs(M.real)+np.abs(M.imag), axis=0)))


def build(s0):
    A, addr = G.read(G.ROOT/'build/spectrum/rt-A-12x36.dat')
    n = len(addr)
    Qm = np.linalg.inv(G.free(addr))
    w = np.array([G.ETA[d]*(2-(m == 0))*(2-(nn == 0))*(65/64)**m*(5/4)**nn for d, p, m, nn in addr])
    H = (w[:, None]*((A+s0*np.eye(n))@Qm)/w[None, :]).astype(complex)
    return H, addr


def curve(R):
    """||R^k|| ate KMAX, com o pico (a corcova da nao-normalidade) e os k_min."""
    d, Pw = {}, np.eye(len(R), dtype=complex)
    hit = {t: None for t in TARGETS}
    peak = 0.0
    for k in range(1, KMAX+1):
        Pw = Pw@R
        c = colnorm(Pw)
        peak = max(peak, c)
        d[k] = c
        for t in TARGETS:
            if hit[t] is None and c <= t:
                hit[t] = k
        if all(hit[t] is not None for t in TARGETS):
            break
    return dict(kmin_0_9=hit[0.9], kmin_0_5=hit[0.5], kmin_1e_3=hit[1e-3], kmin_1e_5=hit[1e-5],
                peak=peak,
                curve={k: d[k] for k in d if k in (1, 2, 4, 8, 16, 24, 32, 48, 64, 72)},
                last_k=max(d), last=d[max(d)])


def damp(M, grid=np.linspace(0.1, 1.8, 35)):
    """melhor amortecimento omega para I-omega*M, pelo raio espectral (barato)."""
    lam = np.linalg.eigvals(M)
    best = min(((float(np.max(np.abs(1-om*lam))), float(om)) for om in grid))
    return best[1], best[0]


def main():
    out = dict(note='DIAGNOSTICO float; nao e prova; alvo na norma de coluna', centers={})
    for s0 in CENTERS:
        H, addr = build(s0)
        rec = {}
        for cm, cnn in CORES:
            core = np.array([m < cm and nn < cnn for d, p, m, nn in addr])
            O = np.where(~core)[0]
            HOO, oaddr = H[np.ix_(O, O)], [addr[i] for i in O]
            N = len(O)
            I = np.eye(N, dtype=complex)
            lam, V = np.linalg.eig(HOO)
            sub = dict(shell_dofs=N, cand={})

            def add(name, P, cost):
                om, rho = damp(HOO@P)
                R = I-om*(HOO@P)
                c = curve(R)
                c.update(omega=om, rho=rho, custo_por_aplicacao=cost)
                sub['cand'][name] = c

            add('neumann_amortecida', I, '1 aplicacao de H')
            add('jacobi', np.diag(1/np.diag(HOO)), '1 aplicacao + 1 divisao diagonal')
            for name, grp in (('bloco_m', 2), ('bloco_n', 3)):
                P = np.zeros_like(HOO)
                for kv in sorted({a[grp] for a in oaddr}):
                    ii = np.array([i for i, a in enumerate(oaddr) if a[grp] == kv])
                    P[np.ix_(ii, ii)] = np.linalg.inv(HOO[np.ix_(ii, ii)])
                add(name, P, f'1 aplicacao + solve por bloco ({len({a[grp] for a in oaddr})} blocos)')
            # deflacao pelo projetor espectral exato (autovetores a esquerda = linhas de V^-1)
            Vi = np.linalg.inv(V)
            om0, _ = damp(HOO)
            dv = np.abs(1-om0*lam)
            ordv = np.argsort(-dv)
            for r in (8, 32):
                sel = ordv[:r]
                Pdef = (V[:, sel]*(1/lam[sel]))@Vi[sel, :]
                Prest = I-(V[:, sel])@Vi[sel, :]
                add(f'deflacao_espectral_r{r}', Pdef+om0*Prest, f'1 aplicacao + posto {r}')
            # onde moram os autovalores que sobram longe de 1, no omega bom
            sub['espectro_no_omega_bom'] = dict(
                omega=om0, rho=float(np.max(dv)),
                n_above_0_5=int(np.sum(dv > 0.5)), n_above_0_7=int(np.sum(dv > 0.7)),
                n_above_0_9=int(np.sum(dv > 0.9)),
                massa_m=_mass(V[:, ordv[:32]], oaddr, 2), massa_n=_mass(V[:, ordv[:32]], oaddr, 3))
            rec[f'core{cm}x{cnn}'] = sub
            print(f'  {cm}x{cnn} ok', flush=True)
        out['centers'][str(s0)] = rec
        print(f'--- s={s0} feito', flush=True)
    p = G.ROOT/'build/spectrum/shell-preconditioner-probe.json'
    p.write_text(json.dumps(out, indent=2)+'\n')
    print('escrito:', p)


def _mass(Vr, oaddr, grp):
    w = np.sum(np.abs(Vr)**2, axis=1)
    w = w/w.sum()
    agg = {}
    for i, a in enumerate(oaddr):
        agg[a[grp]] = agg.get(a[grp], 0.0)+float(w[i])
    return [[int(k), round(v, 4)] for k, v in sorted(agg.items(), key=lambda kv: -kv[1])[:6]]


if __name__ == '__main__':
    main()
