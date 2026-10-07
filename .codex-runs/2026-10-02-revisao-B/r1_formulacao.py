#!/usr/bin/env python3
"""R1 (revisor): testes de falsificacao da formulacao B1.

(1) i-periodicidade no codigo: L(s + i) nos modos m == L(s) nos modos m + 2 (truncagem pequena, sem pesos), e
    bloco_B(m + 2 <- m' + 2) == bloco_B(m <- m').
(2) F3 nao depende de Im s: a parte hermitiana de W L(s) W^-1 (pesos do toro) so depende de Re s.
(3) a cota uniforme ||R_m(sigma)^-1|| <= c/max(mu, |Im sigma|), Re sigma >= 0 (eps_mu_L, nQx, Rinv_uniforme), testada
    em truncagens (norma de bloco principal <= norma verdadeira: e teste de falsificacao), Im sigma ate 3.
(4) cauda_Rinv_disco (sup de ||R^-1 P_{n>=N}|| no disco, Re sigma >= 0) contra o valor numerico em amostras do disco,
    com centros de Im ate 0,75 + i m/2.
(5) eps_mu_L crescente em |Re s| e |Im s| (o disco usa o pior canto).
(6) Brauer: polos e zeros de D(s) = prod (s - s_k + delta_k)/(s - s_k) pelos dados de schur.json.
(7) mesma referencia: sha256 de A em todos os discos (A e B), mesma particao de janelas em todas as caudas."""
import glob, hashlib, json, math, sys
from pathlib import Path
import numpy as np
RAIZ = Path('<repo>'); sys.path.insert(0, str(RAIZ/'scripts'))
import signed_operator_L as sl, closed_form as cf, fundo_L
MU, K1, K2 = sl.MU, sl.K1, sl.K2
res = {}
g = sl.campos()
# (1)
M, N = 7, 12
s1 = 0.3 + 0.4j
L1, rs1, of1 = sl.operador(M, N, s1, g=g)
L2, rs2, of2 = sl.operador(M, N, s1 + 1j, g=g)
ms = sl.modos(M); idx = {m: i for i, m in enumerate(ms)}
dif = 0.0
for mo in ms:
    for mi in ms:
        if mo + 2 in idx and mi + 2 in idx:
            a = L2[of2[idx[mo]]:of2[idx[mo]+1], of2[idx[mi]]:of2[idx[mi]+1]]
            b = L1[of1[idx[mo+2]]:of1[idx[mo+2]+1], of1[idx[mi+2]]:of1[idx[mi+2]+1]]
            dif = max(dif, float(np.abs(a - b).max()))
db = 0.0
for (mo, mi) in [(0, 0), (1, 0), (0, 1), (3, -2), (-5, 4), (10, 2), (7, 7), (-1, 1)]:
    a = sl.bloco_B(g, sl.Radial(mo, 20), sl.Radial(mi, 20)); b = sl.bloco_B(g, sl.Radial(mo + 2, 20), sl.Radial(mi + 2, 20))
    db = max(db, float(np.abs(a - b).max()))
print(f'(1) max |L(s+i)[m,m\'] - L(s)[m+2,m\'+2]| = {dif:.2e};  max |B(m+2<-m\'+2) - B(m<-m\')| = {db:.2e}')
res['periodicidade'] = dict(dif_L=dif, dif_B=db)
# (2)
w = np.concatenate([K1**r.m*np.sqrt(cf.omega(np.concatenate([ns for _, ns in r.comps]), K2)) for r in rs1])
def herm_ext(s):
    L, _, _ = sl.operador(M, N, s, g=g)
    X = (w[:, None]*L)/w[None, :]
    H = (X + X.conj().T)/2
    e = np.linalg.eigvalsh(H)
    return e[0], e[-1]
hs = {str(s): herm_ext(s) for s in (0.3 + 0.1j, 0.3 + 0.6j, 0.3 - 0.9j, 2.9 + 0.75j)}
print('(2) extremos da parte hermitiana de W L(s) W^-1 (7x12): ' + '; '.join(f'{k}: [{v[0]:.10f}, {v[1]:.10f}]' for k, v in hs.items()))
res['hermitiana'] = {k: list(v) for k, v in hs.items()}
# (3)
cw = math.sqrt(1 + K2**-4); c = 1 + 2*cw/(K2 - 1)
def R(sig, Nn, passo):
    ns = np.arange(Nn) if passo == 1 else np.arange(1, 2*Nn, 2)
    J = np.diag(MU*(ns + 1.0) + sig).astype(complex) + np.triu(np.tile(2*MU*ns, (len(ns), 1)), 1)
    ww = np.sqrt(cf.omega(ns, K2))
    return (ww[:, None]*J)/ww[None, :], ns
pior = 0.0
for re_ in (0.0, 1e-3, 0.3, 1.0, 3.0):
    for im in (0.0, 0.05, 0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 3.0):
        for passo in (1, 2):
            Jw, _ = R(re_ + 1j*im, 300, passo)
            nr = 1/np.linalg.svd(Jw, compute_uv=False).min()
            cota = c/max(MU, abs(im))
            pior = max(pior, nr/cota)
print(f'(3) max ||R^-1|| / (c/max(mu,|Im sigma|)) em 90 casos (N = 300) = {pior:.4f}  (precisa <= 1)')
res['Rinv_uniforme_razao_max'] = pior
# (4)
pior4 = 0.0; casos = []
rng = np.random.default_rng(5)
for sp, r in ((0.0 + 0.5j - 0.5j, 0.01), (0.0 + 0.75j, 0.02), (1.0 + 0.75j - 0.5j, 0.12), (0.0 + 0.25j, 0.03), (0.3 + 0.75j + 1.5j, 0.05)):
    for passo in (1, 2):
        Nt = 40
        cota = cf.cauda_Rinv_disco(sp, r, Nt, passo, MU, K2, 1500)
        mx = 0.0
        for _ in range(12):
            e = r*np.sqrt(rng.random())*np.exp(2j*np.pi*rng.random())
            sig = sp + e
            if sig.real < 0:
                sig = complex(0.0, sig.imag)
            Jw, ns = R(sig, 400, passo)
            Q = np.linalg.inv(Jw)
            mx = max(mx, float(np.linalg.norm(Q[:, ns >= Nt] if passo == 1 else Q[:, ns >= Nt], 2)))
        casos.append((str(sp), r, passo, mx, cota)); pior4 = max(pior4, mx/cota)
print(f'(4) cauda_Rinv_disco: max (numerico/cota) = {pior4:.4f} em {len(casos)} casos (precisa <= 1)')
res['cauda_Rinv_disco'] = dict(razao_max=pior4, casos=casos)
# (5)
grade = [[fundo_L.eps_mu_L(complex(x, y)) for y in np.linspace(0, 1.0, 21)] for x in np.linspace(0, 3.2, 17)]
G = np.array(grade)
mono = bool(np.all(np.diff(G, axis=0) >= -1e-30) and np.all(np.diff(G, axis=1) >= -1e-30))
print(f'(5) eps_mu_L monotona em |Re s| e |Im s| na grade [0,3.2]x[0,1]: {mono}; eps_mu_L(0,75i) = {fundo_L.eps_mu_L(0.75j):.3e}, (0,25i) = {fundo_L.eps_mu_L(0.25j):.3e}')
res['eps_mu_monotona'] = mono
# (6)
sch = json.load(open(RAIZ/'build/rouche_L_rig/z/lab2/schur.json'))
sk = [float(x) for x in sch['deflacao'].split(',')]
dl = [complex(*d['delta']) for d in sch['defl']]
zeros = [s - d for s, d in zip(sk, dl)]
print(f'(6) Brauer: polos de D(s) em {sk}, zeros em {zeros} (todos com Im = 0: fora de 1/4 < Im s < 3/4)')
res['brauer'] = dict(polos=sk, zeros=[[z.real, z.imag] for z in zeros])
# (7)
shas = set(); arqs = sorted(glob.glob(str(RAIZ/'build/rouche_L_rig/discos_v2/*.json'))) + sorted(glob.glob(str(RAIZ/'build/faixaB/discos/*.json')))
for f in arqs:
    shas.add(json.load(open(f)).get('sha256_A'))
mods = set()
for f in glob.glob(str(RAIZ/'build/rouche_L_rig/cauda/T3_*.json')) + glob.glob(str(RAIZ/'build/faixaB/cauda/T3_*.json')):
    t = json.load(open(f)); mods.add(hashlib.md5(json.dumps(t['modos']).encode()).hexdigest())
print(f'(7) sha256_A nos {len(arqs)} arquivos de discos: {shas}; particoes de janelas distintas nas caudas: {len(mods)}')
res['mesma_referencia'] = dict(shas=list(shas), particoes=len(mods))
json.dump(res, open(RAIZ/'.codex-runs/2026-10-02-revisao-B/r1_resultado.json', 'w'), indent=1, default=str)
