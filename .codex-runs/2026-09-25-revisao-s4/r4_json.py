#!/usr/bin/env python3
"""R4 (revisao S4): conferencias de consistencia dos JSON de S4 (lambda0, rt, eps_mu, juncao dos
checkpoints, origem na mesma estimativa) e comparacao do Q.json recalculado (rig_copia)."""
from __future__ import annotations
import json, math, sys
from pathlib import Path
import numpy as np

RAIZ = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(RAIZ/'scripts'))
import fundo_L, nk_L as nk, signed_operator_L as sl, closed_form as cf
from tile_perron import DELTA_MU

AQUI = Path(__file__).resolve().parent
rig = RAIZ/'build/s4/rig'
out = []
def P(s=''):
    print(s, flush=True); out.append(s)

z = json.loads((rig/'Z.json').read_text()); t = json.loads((rig/'T3.json').read_text())
c3 = json.loads((rig/'C3.json').read_text()); e = json.loads((rig/'E.json').read_text())
muA = sl.MU
P(f'mu_A = {muA!r}; Z.json lam0 = {z["lam0"]!r} ({z["lam0"] == muA}); C3.json lambda0 = {c3["lambda0"]!r} ({c3["lambda0"] == muA})')
P(f'EPS_FUNDO_L no ambiente desta conta = {fundo_L.EPS_FUNDO_L}; Z.json dB_L.rt = {z["dB_L"]["rt"]}; T3.json dB_L.rt = {t["dB_L"]["rt"]}; E.json rt = {e["rt"]}')
emu_mu = fundo_L.eps_mu_L(complex(muA)); emu_s3b = fundo_L.eps_mu_L(complex(0.4010247330))
P(f'eps_mu_L(mu_A) = {emu_mu!r}; eps_mu_L(0,401) = {emu_s3b!r}; Z.json eps_mu = {z["eps_mu"]!r} ({z["eps_mu"] == emu_mu}); T3.json eps_mu = {t["eps_mu"]!r} ({t["eps_mu"] == emu_mu})')
# dB_L recalculado com o fundo RefA e banda 16
fb = fundo_L.dB_L(sl.campos(), 16)
P(f'dB_L(banda 16) recalculado = {fb}; igual ao de Z.json: {fb == z["dB_L"]}; ao de T3.json: {fb == t["dB_L"]}')
# beta e pert de T3 recalculados (dependem de lambda0 via Qfar de Z.json e via eps_mu)
beta = (nk.NORMA_BL + fb['rt'])*z['Qfar'] + emu_mu
P(f'T3 beta = {t["beta"]!r}; recalculado (NORMA_BL + rt) Qfar(Z.json) + eps_mu(mu_A) = {beta!r}; igual: {beta == t["beta"]}')
P(f'T3 pert = {t["pert"]!r}; dB_L.total + eps_mu = {fb["total"] + emu_mu!r}; igual: {t["pert"] == fb["total"] + emu_mu}')
# Qfar recalculado em lambda0 = mu_A (e em 0,401 para contraste)
Mc, Nc = z['F']
def qfar(l0):
    q = max(max(cf.cauda_Rinv(l0 + 1j*m/2, Nc, p, sl.MU, sl.K2, 1500) for p in ((2, 1, 1) if m % 2 == 0 else (1,))) for m in range(-(Mc - 1), Mc))
    return max(q, cf.Rinv_uniforme(Mc, sl.MU, l0, sl.K2))
qf = qfar(muA)
P(f'Qfar(mu_A) recalculado = {qf!r}; Z.json Qfar = {z["Qfar"]!r}; igual: {qf == z["Qfar"]}; (em 0,401 seria {qfar(0.4010247330):.6f})')
s3b = RAIZ/'build/s3b/rig/T3.json'
if s3b.exists():
    ts = json.loads(s3b.read_text())
    P(f'contraste S3b (lambda0 = 0,401): T3 beta = {ts["beta"]:.6f}, eps_mu = {ts["eps_mu"]:.3e}, rt = {ts["dB_L"]["rt"]:.3e}')
# juncao dos checkpoints
a = json.loads((rig/'T3.lab.json').read_text()); b = json.loads((rig/'T3.lab2.json').read_text())
nomes = t['nome']
P(f'checkpoints: lab {len(a)} janelas, lab2 {len(b)}; intersecao {sorted(set(a) & set(b))}; uniao == nomes de T3: {sorted(set(a) | set(b)) == sorted(nomes)}; 26 janelas: {len(nomes) == 26}')
ok = True
for i, nm in enumerate(nomes):
    f = a.get(nm) or b.get(nm)
    for k in ('TZ', 'nV', 'nVhi', 'Y'):
        if f[k] != t[k][i]:
            ok = False; P(f'  DIFERE {nm} {k}: {f[k]} vs {t[k][i]}')
    Nrow = list(f['N']); Nrow[i] = f['rho']
    if Nrow != t['N'][i]:
        ok = False; P(f'  DIFERE {nm} linha N')
P(f'T3.json = juncao exata dos checkpoints (N com rho na diagonal): {ok}')
# mesma estimativa: cota = chute*(1 + 1e-4)*1.002^k*FATOR_PESO (cota_loewner a partir de est*(1+1e-4))
est = json.loads((RAIZ/'build/s4/nk3_mu.json').read_text())
P(f'nk3_mu.json: nomes iguais aos de T3: {est["nome"] == nomes}; beta estimado {est["beta"]:.8f} (T3 {t["beta"]:.8f})')
def passos(cota, chute):
    if chute <= 0:
        return None
    x = cota/(chute*(1 + 1e-4)*(1 + 1e-12))
    k = math.log(x)/math.log(1.002)
    return k
maus = []
for i, nm in enumerate(nomes):
    for k in ('nV',):
        kk = passos(t[k][i], est[k][i])
        if kk is None or abs(kk - round(kk)) > 1e-6 or kk < -1e-9:
            maus.append((nm, k, kk))
    j_viz = [j for j in range(len(nomes)) if est['N'][i][j] > 0 and j != i]
    for j in j_viz:
        # N certificado = (t + nV*erro)*FATOR: nao e exatamente chute*1.002^k; so exige cota >= chute
        if t['N'][i][j] < est['N'][i][j]:
            maus.append((nm, f'N<-{nomes[j]}', t['N'][i][j], est['N'][i][j]))
    if t['TZ'][i] and t['TZ'][i] < est['TZ'][i]:
        maus.append((nm, 'TZ', t['TZ'][i], est['TZ'][i]))
P('nV certificado = estimativa(nk3_mu) * 1,0001 * 1,002^k * (1 + 1e-12), k inteiro, em todas as 26 janelas: '
  + ('SIM' if not maus else f'NAO: {maus}'))
ks = [round(passos(t['nV'][i], est['nV'][i])) for i in range(len(nomes))]
P(f'  k por janela (nV): {ks}')
# Q.json recalculado
qo = json.loads((rig/'Q.json').read_text())
qn_p = AQUI/'rig_copia/Q.json'
if qn_p.exists():
    qn = json.loads(qn_p.read_text())
    P(f'Q.json original: {qo}')
    P(f'Q.json recalculado (nk_L3_q.py --dir rig_copia, lam0 de Z.json = mu_A): {qn}')
    P(f'identicos: {qo == qn}')
else:
    P('Q.json recalculado ainda nao existe')
(AQUI/'r4_saida.txt').write_text('\n'.join(out) + '\n')
