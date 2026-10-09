#!/usr/bin/env python3
"""N1: integridade e coerência dos dados de entrada de build/nk733 (código do revisor).

Confere: um só lambda0 em Z.json, Q.json, T3.json, C3.json; as caixas (Z 32x128, F 200x400, banda 16,
janelas W = 16, folga 60); que as 26 janelas de T3.json são exatamente as de nk_L3.janelas_L3(32, 200, 16)
(reconstruídas aqui SEM importar o código do autor); que T3.json é a junção EXATA de T3.parcial.json e
T3.joao.json (N, TZ, nV, nVhi, Y; rho na diagonal de N); que as constantes partilhadas (dB_L, eps_mu, nQZ,
Qfar) coincidem; h0w.npy (dimensão, norma); e que os números dos logs batem com os JSON."""
import json, math, re, sys, hashlib
from pathlib import Path
import numpy as np

D = Path('<repo>/build/nk733')
z = json.loads((D/'Z.json').read_text()); q = json.loads((D/'Q.json').read_text())
t = json.loads((D/'T3.json').read_text()); c = json.loads((D/'C3.json').read_text())
par = json.loads((D/'T3.parcial.json').read_text()); jo = json.loads((D/'T3.joao.json').read_text())
falhas = []
def chk(cond, msg):
    print(('OK    ' if cond else 'FALHA ') + msg)
    if not cond:
        falhas.append(msg)

# 1. lambda0
lz = z['lam0']; lq = q['lam0']; lt = t['lam0']; lc = c['lambda0']
print('lam0: Z', repr(lz), ' Q', repr(lq), ' T3', repr(lt), ' C3', repr(lc))
chk(lz == lq == lc == lt[0] and lt[1] == 0.0, 'um so lambda0 (Z, Q, T3, C3), parte imaginaria 0 em T3')
chk(lz == 0.7331796596606924, 'lambda0 = 0.7331796596606924 (o double do artigo)')
auto = (D/'autovalor.txt').read_text().split()
chk(float(auto[auto.index('LAMBDA') + 1]) == lz, 'autovalor.txt: LAMBDA igual a lam0')

# 2. caixas
chk(z['Z'] == [32, 128] and z['F'] == [200, 400] and z['banda'] == 16, 'Z.json: Z 32x128, F 200x400, banda 16')
# dimensão de Z: modos |m| < 32, n < 128; o layout radial tem 4 componentes; n = 14016 vem do log
chk(z['n'] == 14016, 'Z.json: n = 14016')
chk(t['folga_far'] == 60, 'T3.json: folga_far 60')
chk(t['eps_fundo_L'] == 4.7e-08 and z['dB_L']['rt'] == 4.7e-08, 'eps_fundo_L = 4,7e-8 em Z e T3')
chk(t['rouche'] is False, 'T3.json: rouche False (NK, com residuo)')

# janelas reconstruídas sem o código do autor: central |m| < c0 (c0 = 32 mod 16 = 0 -> 16) dividida
# em duas, depois [c0 + kW, c0 + (k+1)W) cortadas em Mc, positivas e depois negativas (reflexão)
MZ, Mc, W = 32, 200, 16
c0 = MZ % W or W
pos = [list(range(k, min(k + W, Mc))) for k in range(c0, Mc, W)]
neg = [[-m for m in reversed(j)] for j in pos]
cen = list(range(-(c0 - 1), c0)); meio = len(cen)//2
esperado = [cen[:meio + 1], cen[meio + 1:]] + pos + neg
chk(t['modos'] == esperado, f'T3.json: as {len(esperado)} janelas sao as de W = 16 (central dividida)')
chk(len(t['nome']) == 26 and len(set(t['nome'])) == 26, 'T3.json: 26 nomes distintos')
cob = sorted(m for j in t['modos'] for m in j)
chk(cob == list(range(-199, 200)), 'janelas cobrem |m| < 200 sem sobreposicao')
for k in ('N', 'TZ', 'nV', 'nVhi', 'Y', 'nQ', 'K', 'nQm'):
    chk(len(t[k]) == 26, f'T3.json: {k} com 26 entradas')
chk(len(t['Nf']) == 27 and all(len(r) == 26 for r in t['Nf']), 'T3.json: Nf 27 x 26 (25 grupos de F + 2 do far de Fourier)')
SZ_esp = [i for i, ms in enumerate(t['modos']) if min(abs(m) for m in ms) < MZ + 40]
chk(t['SZ'] == SZ_esp, f'T3.json: SZ = janelas com min|m| < 72: {SZ_esp}')

# 3. constantes partilhadas
chk(t['dB_L'] == z['dB_L'], 'dB_L (banda 16) igual em Z e T3')
chk(t['eps_mu'] == z['eps_mu'], 'eps_mu igual em Z e T3')
chk(t['nQZ'] == z['nQZ'] and t['Qfar'] == z['Qfar'], 'T3 leu nQZ e Qfar deste Z.json')
chk(abs(t['pert'] - (z['dB_L']['total'] + z['eps_mu'])) <= 1e-18, 'pert = dB_L.total + eps_mu')
beta = (3.964043 + z['dB_L']['rt'])*z['Qfar'] + z['eps_mu']
chk(beta == t['beta'], f'beta = (NORMA_BL + rt) Qfar + eps_mu = {beta!r} (T3 {t["beta"]!r})')
chk('tQ' in z, 'Z.json tem tQ (nk_L2_Z --so-tQ)')

# 4. junção
uni = dict(par); dup = set(par) & set(jo); uni.update(jo)
chk(not dup, f'parcial e joao disjuntos ({len(par)} + {len(jo)})')
chk(sorted(uni) == sorted(t['nome']), 'parcial U joao = as 26 janelas')
N = np.array(t['N'])
dif = 0
for i, nm in enumerate(t['nome']):
    f = uni[nm]
    linha = list(f['N']); linha[i] = f['rho']
    for k in ('TZ', 'nV', 'nVhi', 'Y'):
        if f[k] != t[k][i]:
            dif += 1; print('   difere', nm, k, f[k], t[k][i])
    if linha != t['N'][i]:
        dif += 1; print('   difere', nm, 'N')
    if f['N'][i] != 0.0:
        dif += 1; print('   diagonal do checkpoint nao nula', nm)
chk(dif == 0, 'T3.json = juncao EXATA (bit a bit) dos checkpoints, rho na diagonal')
orig = {nm: ('parcial' if nm in par else 'joao') for nm in t['nome']}
print('   origem:', ', '.join(f'{nm}:{o}' for nm, o in orig.items()))

# 5. h0w
h = np.load(D/'h0w.npy')
print(f'h0w.npy: shape {h.shape}, dtype {h.dtype}, ||h|| - 1 = {np.linalg.norm(h) - 1:.2e}, sha256 '
      + hashlib.sha256((D/'h0w.npy').read_bytes()).hexdigest()[:16])
chk(h.shape == (14016,) and abs(np.linalg.norm(h) - 1) < 1e-14, 'h0w: 14016 entradas, norma 1')

# 6. logs x JSON
eA = (D/'etapaA.txt').read_text()
mm = re.search(r"Z' RIGOROSO: \|\|Mh\^-1\|\| <= ([\d.]+), a <= ([\d.]+), rho <= ([\d.e+-]+), Y_Z <= ([\d.e+-]+), eps_far <= ([\d.e+-]+)", eA)
vals = list(map(float, mm.groups()))
print('etapaA:', vals)
chk(f"{z['tV']:.5f}" == mm.group(1) and f"{z['rhoZ']:.2e}" == mm.group(3) and f"{z['YZ']:.2e}" == mm.group(4),
    'etapaA.txt: t_V, rho, Y_Z batem com Z.json')
print(f"   a do log {mm.group(2)} vs Z.json aZ {z['aZ']:.5f}; eps_far do log {mm.group(5)} vs Z.json epsZ {z['epsZ']:.2e}")
tq = (D/'tQ.txt').read_text()
chk(f"{z['tQ']:.5f}" in tq, f"tQ.txt: tQ {z['tQ']:.5f} gravado em Z.json")
qq = (D/'q.txt').read_text()
chk(f"{q['qZT']:.5f}" in qq and f"{q['qTT']:.5f}" in qq and f"{q['nQx_extra']:.5f}" in qq, 'q.txt bate com Q.json')
chk('48 <= |m| < 72' in qq, 'q.txt: nQx_extra sobre 48 <= |m| < 72 (MZ + banda .. MZ + 40)')
eB = (D/'etapaB.txt').read_text() + (D/'etapaB_joao.txt').read_text()
nb = 0
for i, nm in enumerate(t['nome']):
    mm = re.search(re.escape(f'  {nm}: ||V|| <= ') + r'([\d.]+), \|\|V P_alto\|\| <= ([\d.]+), rho ([\d.e+-]+);.*; Y ([\d.e+-]+)', eB)
    if mm is None:
        nb += 1; print('   sem linha no log:', nm); continue
    ok = (f"{t['nV'][i]:.4f}" == mm.group(1) and f"{t['nVhi'][i]:.4f}" == mm.group(2)
          and f"{t['N'][i][i]:.1e}" == mm.group(3) and f"{t['Y'][i]:.1e}" == mm.group(4))
    if not ok:
        nb += 1; print('   log difere:', nm, mm.groups(), t['nV'][i], t['nVhi'][i], t['N'][i][i], t['Y'][i])
chk(nb == 0, 'logs etapaB*.txt batem com T3.json nas 26 janelas (nV, nVhi, rho, Y)')
print('\nFALHAS:', len(falhas))
for f in falhas:
    print('  -', f)
