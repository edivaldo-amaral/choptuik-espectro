#!/usr/bin/env python3
"""N2: montagem de Perron de 3 níveis (Z', T, far) e NK, em racionais exatos, com código do revisor.

Não importa nada de scripts/. Cada float dos JSON entra como o racional EXATO do double (Fraction(float)).
Normas espectrais de matrizes pequenas não negativas por Collatz-Wielandt exato:
  ||A||_2^2 = rho(A^T A) <= max_i (A^T A u)_i / u_i  para qualquer u > 0  (A^T A >= 0),
com u = vetor de Perron aproximado (float, só propõe). Raiz quadrada para cima por racional verificado (s^2 >= x).

Duas montagens:
  (R) a do REVISOR: as entradas derivadas em N3, com a perturbação global (erro de RT + eps_mu) contada UMA vez
      por bloco, sem os fatores de folga ARRED = 1 + 1e-9 do autor;
  (A) a do AUTOR reescrita (mesmas fórmulas de nk_L3_C.py, incluindo as somas duplicadas e ARRED), para
      comparar com C3.json bit a bit onde der.
Em ambas: w = 1 (os pesos escolhidos no C3.json) e também os pesos do C3.json lidos do arquivo."""
import json, math
from fractions import Fraction as Fr
from pathlib import Path
import numpy as np

D = Path('<repo>/build/nk733')
z = json.loads((D/'Z.json').read_text()); qj = json.loads((D/'Q.json').read_text())
t = json.loads((D/'T3.json').read_text()); c3 = json.loads((D/'C3.json').read_text())
F = lambda x: Fr(float(x))                     # racional exato do double
ARRED = 1 + Fr(1, 10**9)
k = 26

def sqrt_sup(x: Fr) -> Fr:
    if x == 0:
        return Fr(0)
    s = Fr(math.sqrt(float(x)))
    while s*s < x:
        s *= 1 + Fr(1, 10**15)
    return s

def sqrt_inf(x: Fr) -> Fr:
    s = Fr(math.sqrt(float(x)))
    while s*s > x:
        s *= 1 - Fr(1, 10**15)
    return s

def norma2_sup(A):
    """cota superior exata de ||A||_2, A lista de listas de Fr >= 0."""
    m, n = len(A), len(A[0])
    G = [[sum(A[p][i]*A[p][j] for p in range(m)) for j in range(n)] for i in range(n)]
    Gf = np.array([[float(x) for x in l] for l in G])
    lam, V = np.linalg.eigh(Gf)
    u = np.abs(V[:, -1]); u = np.maximum(u/u.max(), 1e-12)
    uq = [Fr(float(x)) for x in u]
    rho = max(sum(G[i][j]*uq[j] for j in range(n))/uq[i] for i in range(n))
    return sqrt_sup(rho), math.sqrt(max(lam[-1], 0))

# ---------------- dados
N = [[F(x) for x in l] for l in t['N']]           # diagonal = rho_J
Nf = [[F(x) for x in l] for l in t['Nf']]
TZ = [F(x) for x in t['TZ']]; nV = [F(x) for x in t['nV']]; nVhi = [F(x) for x in t['nVhi']]
Y = [F(x) for x in t['Y']]; nQ = [F(x) for x in t['nQ']]
SZ = t['SZ']; modos = t['modos']
beta, desc = F(t['beta']), F(t['desc'])
rt, emu = F(z['dB_L']['rt']), F(z['eps_mu'])
resto, total = F(z['dB_L']['resto']), F(z['dB_L']['total'])
pertT = F(t['pert'])
tV, rhoZ, YZ, epsZ0 = F(z['tV']), F(z['rhoZ']), F(z['YZ']), F(z['epsZ'])
nQZ, Qf, nQT, erroY = F(z['nQZ']), F(z['Qfar']), F(z['nQT']), F(z['erro_Y'])
a_loew = F(z['a_loewner']); tQ = F(z['tQ'])
qZT, qTT, qZf, nQx = F(qj['qZT']), F(qj['qTT']), F(qj['qZf']), F(qj['nQx_extra'])
borda = [max(abs(m) for m in ms) >= 160 for ms in modos]
nQb = max(nQ[i] for i in range(k) if max(abs(m) for m in modos[i]) >= 176)
nQmax = max(nQ)
eg = rt + emu

def montagem(w, autor):
    wq = [F(x) for x in w]
    wmin = min(wq); wminSZ = min(wq[j] for j in SZ)
    WNW = [[wq[i]*N[i][j]/wq[j] for j in range(k)] for i in range(k)]
    nT, nT_f = norma2_sup(WNW)
    NfW = [[Nf[g][j]/wq[j] for j in range(k)] for g in range(len(Nf))]
    fT, fT_f = norma2_sup(NfW)
    M = [[None]*3 for _ in range(3)]
    if not autor:
        # REVISOR (derivação de N3): perturbação global (RT + mu) UMA vez por bloco, com nQmax; sem ARRED
        maxwV = max(wq[i]*nV[i] for i in range(k))
        aZ = sqrt_sup(a_loew**2 + (tV*resto*nQx)**2) + tV*erroY
        M[0][0] = rhoZ
        M[0][1] = aZ/wminSZ + tV*eg*nQmax/wmin
        M[0][2] = epsZ0 + tV*qZf
        TZw = sqrt_sup(sum((wq[i]*TZ[i])**2 for i in range(k)))
        M[1][0] = TZw + maxwV*eg*nQZ
        M[1][1] = nT + maxwV*eg*nQmax/wmin
        Tf = max(wq[i]*(nV[i] if borda[i] else nVhi[i]) for i in range(k))*beta \
            + sqrt_sup(sum((wq[i]*nV[i])**2 for i in range(k)))*desc
        M[1][2] = Tf + maxwV*eg*Qf
        M[2][0] = pertT*nQZ
        M[2][1] = fT + resto*nQb/wmin + eg*nQmax/wmin
        M[2][2] = beta
    else:
        # AUTOR: fórmulas de nk_L3_C.py (Q(x) = x (1 + 1e-9), aZ_rig e epsZ em float, somas duplicadas);
        # só as normas ||W N W^-1|| e ||Nf W^-1|| saem da cota do revisor (o autor usa Rump, em float)
        q = lambda x: F(x)*ARRED
        a_extra = z['tV']*z['dB_L']['resto']*qj['nQx_extra']
        aZ_rig = (math.sqrt(z['a_loewner']**2 + a_extra**2)
                  + z['tV']*(z['erro_Y'] + (z['dB_L']['rt'] + z['eps_mu'])*z['nQT']))*(1 + 1e-12)
        epsZ = z['epsZ'] + z['tV']*qj['qZf']*1.0001
        b_, d_, p_ = q(t['beta']), q(t['desc']), q(t['pert'])
        eg_ = q(z['dB_L']['rt']) + q(z['eps_mu']); nQmax_ = q(max(t['nQ']))
        nT = nT*ARRED; fT = fT*ARRED
        TZw = sqrt_sup(sum((wq[i]*q(t['TZ'][i]))**2 for i in range(k)))
        ZT = q(aZ_rig)/wminSZ
        Tf = max(wq[i]*q(t['nV'][i] if borda[i] else t['nVhi'][i]) for i in range(k))*b_ \
            + sqrt_sup(sum((wq[i]*q(t['nV'][i]))**2 for i in range(k)))*d_
        nQb_ = max(q(t['nQ'][i]) for i in range(k) if max(abs(m) for m in modos[i]) >= 176)
        fT += q(z['dB_L']['resto'])*nQb_/wmin
        maxwV = max(wq[i]*q(t['nV'][i]) for i in range(k))
        nT += maxwV*eg_*nQmax_/wmin; TZw += maxwV*eg_*q(z['nQZ']); Tf += maxwV*eg_*q(z['Qfar'])
        ZT += q(z['tV'])*eg_*nQmax_/wmin; fT += eg_*nQmax_/wmin
        M = [[q(z['rhoZ']), ZT, q(epsZ)], [TZw, nT, Tf], [p_*q(z['nQZ']), fT, b_]]
        aZ = q(aZ_rig)
    return M, dict(nT=nT, nT_float=nT_f, fT=fT, fT_float=fT_f, aZ=aZ, maxwV=maxwV, wmin=wmin, wq=wq)

def perron(M, v=None):
    Mf = np.array([[float(x) for x in l] for l in M])
    if v is None:
        lam, vec = np.linalg.eig(Mf); v = np.abs(vec[:, np.argmax(np.abs(lam))]); v = np.maximum(v/v.max(), 1e-9)
        v = [Fr(float(x)) for x in v]
    th = max(sum(M[i][j]*v[j] for j in range(3))/v[i] for i in range(3))
    return th, v, max(abs(np.linalg.eigvals(Mf)))

def nk(th, v, info, autor=False):
    v0, vT, vf = v
    wmin, maxwV, wq = info['wmin'], info['maxwV'], info['wq']
    if not autor:
        Z2 = 2*v0*max((tQ*v0 + tV*(qZT/wmin*vT + qZf*vf))/v0, maxwV*(qTT/wmin*vT + Qf*vf)/vT, Qf)
        YT = sqrt_sup(sum((wq[i]*Y[i])**2 for i in range(k))) + maxwV*rt
        Ymax = max(YZ/v0, YT/vT, rt/vf); rtv = rt
    else:
        q = lambda x: x*ARRED
        Z2 = 2*v0*max((q(tQ)*v0 + q(tV)*(q(qZT)/wmin*vT + q(qZf)*vf))/v0, maxwV*(q(qTT)/wmin*vT + q(Qf)*vf)/vT, q(Qf))*ARRED
        rtv = q(rt)*ARRED
        YT = sqrt_sup(sum((wq[i]*q(Y[i]))**2 for i in range(k))) + maxwV*rtv
        Ymax = max(q(YZ)/v0, YT/vT, rtv/vf)
    d = (1 - th)**2 - 4*Z2*Ymax
    assert d > 0
    r_sup = (1 - th - sqrt_inf(d))/(2*Z2)        # cota superior do menor raio
    r_inf = (1 - th - sqrt_sup(d))/(2*Z2)
    return dict(Z2=Z2, Y=Ymax, YZ_v0=(YZ if not autor else YZ*ARRED)/v0, YT_vT=YT/vT, Yf_vf=rtv/vf, disc=d,
                r_inf=r_inf, r_sup=r_sup, erro_lambda_sup=v0*r_sup)

def checa(th, Z2, Ym, r):
    return Ym + th*r + Z2*r*r <= r, th + 2*Z2*r < 1

out = {}
print('C3.json gravado: theta', float(Fr(c3['theta'])), ' Z2', float(Fr(c3['Z2'])), ' Y', float(Fr(c3['Y'])),
      ' r', float(Fr(c3['r'])), ' erro_lambda', float(Fr(c3['erro_lambda'])))
# matriz gravada
Mg = [[Fr(x) for x in l] for l in c3['M']]
vg = [Fr(x) for x in c3['v']]
thg = max(sum(Mg[i][j]*vg[j] for j in range(3))/vg[i] for i in range(3))
print('theta recalculado da M e v GRAVADAS:', float(thg), ' (gravado theta = th*(1+1e-9):', float(Fr(c3['theta'])), ')')
print('  igual a theta/(1+1e-9)?', thg*(1 + Fr(1, 10**9)) == Fr(c3['theta']))
out['theta_de_M_gravada'] = str(thg)
for autor in (False, True):
    nome = 'AUTOR' if autor else 'REVISOR'
    M, info = montagem([1.0]*k, autor)
    th, v, rho_f = perron(M)
    print(f'\n== montagem {nome} (w = 1)')
    for i in range(3):
        print('   ', '  '.join(f'{float(x):.6e}' for x in M[i]), '   gravado:', '  '.join(f'{float(x):.6e}' for x in Mg[i]))
    print(f'   ||N||_2 <= {float(info["nT"]):.8f} (float {info["nT_float"]:.8f});  ||Nf||_2 <= {float(info["fT"]):.8f} (float {info["fT_float"]:.8f})')
    print(f'   theta = {float(th):.9f} (raio espectral float {rho_f:.9f});  v = {[float(x) for x in v]}')
    if autor:
        th = th*(1 + Fr(1, 10**9))
        print(f'   theta*(1+1e-9) = {float(th):.12f} vs gravado {float(Fr(c3["theta"])):.12f}')
    R = nk(th, v, info, autor)
    print(f'   Z2 = {float(R["Z2"]):.6f};  Y = {float(R["Y"]):.6e} (Z\' {float(R["YZ_v0"]):.3e}, T {float(R["YT_vT"]):.4e}, far {float(R["Yf_vf"]):.2e})')
    print(f'   discriminante (1-theta)^2 - 4 Z2 Y = {float(R["disc"]):.8e} > 0')
    print(f'   r in [{float(R["r_inf"]):.10e}, {float(R["r_sup"]):.10e}];  |lambda* - lambda0| <= v0 r <= {float(R["erro_lambda_sup"]):.10e}')
    ok1, ok2 = checa(th, R['Z2'], R['Y'], R['r_sup'])
    print(f'   no r_sup: Y + theta r + Z2 r^2 <= r: {ok1};  theta + 2 Z2 r < 1: {ok2}')
    out[nome] = dict(M=[[str(x) for x in l] for l in M], theta=str(th), v=[str(x) for x in v], Z2=str(R['Z2']),
                     Y=str(R['Y']), disc=str(R['disc']), r_sup=str(R['r_sup']), erro_lambda_sup=str(R['erro_lambda_sup']),
                     floats=dict(theta=float(th), Z2=float(R['Z2']), Y=float(R['Y']), r=float(R['r_sup']),
                                 erro_lambda=float(R['erro_lambda_sup'])))
# condições NK com os números GRAVADOS (r gravado), em racionais
thG, Z2G, YG, rG = Fr(c3['theta']), Fr(c3['Z2']), Fr(c3['Y']), Fr(c3['r'])
ok1, ok2 = checa(thG, Z2G, YG, rG)
dG = (1 - thG)**2 - 4*Z2G*YG
rmin_sup = (1 - thG - sqrt_inf(dG))/(2*Z2G)
print(f'\n== gravado: discriminante {float(dG):.8e}; r gravado {float(rG):.12e} >= raiz menor {float(rmin_sup):.12e}: {rG >= rmin_sup};'
      f' invariante {ok1}; contrai {ok2} (theta + 2 Z2 r = {float(thG + 2*Z2G*rG):.6f})')
print(f'   erro_lambda gravado = v0 r? {Fr(c3["erro_lambda"]) == vg[0]*rG};  = {float(Fr(c3["erro_lambda"])):.12e}')
print(f'   raio de unicidade: theta + 2 Z2 rho < 1  <=>  rho < {float((1 - thG)/(2*Z2G)):.6e}')
out['gravado'] = dict(disc=str(dG), invariante=ok1, contrai=ok2, r_ge_raiz=rG >= rmin_sup,
                      erro_lambda_eq_v0r=Fr(c3['erro_lambda']) == vg[0]*rG)
Path(__file__).with_suffix('.json').write_text(json.dumps(out, indent=1) + '\n')

# comparação entrada a entrada, em racionais: a M gravada domina a do revisor?
Mr, _ = montagem([1.0]*k, False)
dom = all(Mg[i][j] >= Mr[i][j] for i in range(3) for j in range(3))
print('\nM gravada >= M do revisor, entrada a entrada (racionais):', dom)
for i in range(3):
    print('   folga relativa', '  '.join(f'{float((Mg[i][j] - Mr[i][j])/Mr[i][j]):+.2e}' for j in range(3)))
