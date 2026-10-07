"""C2 (revisao CRH, 03/10/2026): o eixo imaginario Re s = 0, Im s em [-1/4, 3/4], com os discos
certificados das faixas A (build/rouche_L_rig/discos_v2) e B (build/faixaB/discos). Codigo proprio,
racionais exatos (Fraction de cada float gravado, que e' o valor usado no Perron racional).

Confere:
  1. cada disco usado tem ok = True, r > 0 e theta (fracao exata gravada) < 1;
  2. a uniao das cordas dos discos na reta Re s = 0 cobre o segmento FECHADO [-1/4, 3/4]
     (meia-corda cotada por baixo em racionais: h_lo^2 <= r^2 - a^2);
  3. quais discos contem s = 0, i/4, i/2, 3i/4, -i/4 (cantos e juncao das faixas), com a folga;
  4. q(s) = (s + 1 + e00)(s - mu* + d1(1 + e11)) - d0 d1 e01 e10 sem zeros em Re s >= 0, com
     |e_ij| <= max||phi_i|| * max||P_Z(x_j - x_j^A)|| (schur.json, E.json, dtauB padrao do rouche_L_disco).
"""
from fractions import Fraction as F
import glob
import json
import math

RAIZ = '<repo>/'
arqs = sorted(glob.glob(RAIZ + 'build/rouche_L_rig/discos_v2/*.json')) + \
    sorted(glob.glob(RAIZ + 'build/faixaB/discos/*.json'))

def sqrt_lo(x):
    """racional y >= 0 com y^2 <= x (x racional >= 0)."""
    if x <= 0:
        return F(0)
    y = F(math.sqrt(float(x))) * F(999999, 1000000)
    while y*y > x:
        y *= F(999999, 1000000)
    return y

discos = []
ruins = []
tot = 0
for f in arqs:
    d = json.load(open(f))
    for p in d['pontos']:
        tot += 1
        a, b, r = F(p['p'][0]), F(p['p'][1]), F(p['r'])
        th = F(p['theta']) if '/' in p['theta'] or p['theta'].replace('.', '').isdigit() else None
        if not (p['ok'] and r > 0 and th is not None and th < 1):
            ruins.append((f.split('/')[-1], p['p'], p['r'], p['theta'][:20], p['ok']))
            continue
        discos.append(dict(arq=f.replace(RAIZ, ''), c=d['c'], a=a, b=b, r=r, th=th))
print(f'{tot} entradas lidas; {len(discos)} discos validos; descartadas: {len(ruins)}')
for x in ruins:
    print('   descartada:', x)

# 2. cordas no eixo
cordas = []
for D in discos:
    if D['a'] < 0:
        print('AVISO: centro com Re < 0', D)
    if D['r'] > abs(D['a']):
        h = sqrt_lo(D['r']**2 - D['a']**2)
        cordas.append((D['b'] - h, D['b'] + h, D))
cordas.sort(key=lambda t: t[0])
lo, hi = F(-1, 4), F(3, 4)
x = lo
usados = []
while x <= hi:
    cand = [c for c in cordas if c[0] <= x < c[1] or (c[0] <= x <= c[1] and x == hi)]
    if not cand:
        print(f'BURACO em Im s = {float(x)}')
        break
    best = max(cand, key=lambda c: c[1])
    usados.append(best)
    if x == hi:
        x = hi + 1
        break
    x = best[1]
    # a proxima corda deve comecar estritamente antes de x (sobreposicao); medir a sobreposicao
    prox = [c for c in cordas if c[0] < x < c[1]]
    if x < hi and not prox:
        print(f'BURACO logo apos Im s = {float(x)} (nenhuma corda contem o extremo no interior)')
        break
cobre = x > hi
print(f'2. cobertura do segmento fechado Re s = 0, Im s em [-1/4, 3/4]: {"SIM" if cobre else "NAO"}; '
      f'{len(usados)} cordas na cadeia gulosa; {len(cordas)} discos tocam o eixo')

# sobreposicao minima entre cordas consecutivas da cadeia
ovs = []
for c1, c2 in zip(usados, usados[1:]):
    ovs.append(c1[1] - c2[0])
print(f'   sobreposicao minima entre cordas consecutivas da cadeia: {float(min(ovs)):.3e}')
print(f'   corda extrema inferior vai ate {float(usados[0][0]):.6f}; superior ate {float(usados[-1][1]):.6f}')

# 3. pontos especiais
for nome, y in (('0', F(0)), ('-i/4', F(-1, 4)), ('i/4', F(1, 4)), ('i/2', F(1, 2)), ('3i/4', F(3, 4))):
    lst = []
    for D in discos:
        d2 = D['a']**2 + (D['b'] - y)**2
        if d2 < D['r']**2:
            folga = float(D['r']) - math.sqrt(float(d2))
            lst.append((folga, D))
    lst.sort(key=lambda t: -t[0])
    print(f'3. s = {nome}: contido no interior de {len(lst)} discos; melhor folga {lst[0][0]:.5f} '
          f'(p = {float(lst[0][1]["a"])}+{float(lst[0][1]["b"])}i, r = {float(lst[0][1]["r"]):.6f}, '
          f'theta = {float(lst[0][1]["th"]):.5f}, {lst[0][1]["arq"]})' if lst else f'3. s = {nome}: NAO COBERTO')
    if nome == '0':
        for fo, D in lst:
            print(f'      p = {float(D["a"])}{float(D["b"]):+.6f}i  r = {float(D["r"]):.6f}  |p| = {abs(float(D["b"])):.6f}'
                  f'  folga = {fo:.6f}  theta = {float(D["th"]):.5f}  {D["arq"]}')
        cent0 = [D for D in discos if D['a'] == 0 and D['b'] == 0]
        print(f'   discos com ponto-centro p EXATAMENTE em 0: {len(cent0)}  (o "disco_0.json" e o centro de CAUDA c = 0)')

# 4. q(s)
sch = json.load(open(RAIZ + 'build/rouche_L_rig/z/lab2/schur.json'))['defl_refA']
E = json.load(open(RAIZ + 'build/s4/rig/E.json'))['e']
dtauB = 1.4563e-09                      # padrao de rouche_L_disco.py (||d_tau RefB||)
nphi = max(F(x) for x in sch['norma_phi'])
ex = max(F(E), F(dtauB) + F(1, 10**60))
e = nphi*ex*F(1001, 1000)               # folga para o arredondamento das normas
d0, d1 = F(sch['delta'][0]), F(sch['delta'][1])
muA = d1 - 1
dmu = F(39, 10**12)
f1 = 1 - e                              # |s + 1 + e00| >= Re s + 1 - |e00| >= 1 - e
f2 = 1 - dmu - d1*e                     # s - mu* + d1(1 + e11) = s + 1 + (muA - mu*) + d1 e11
qlo = f1*f2 - d0*d1*e*e
print(f'4. |e_ij| <= {float(e):.3e}; delta = ({float(d0)}, {float(d1):.10f}); em Re s >= 0: |fator1| >= {float(f1):.9f}, '
      f'|fator2| >= {float(f2):.9f}, |q(s)| >= {float(qlo):.9f} > 0: {qlo > 0}')
print('   q(0) ~', float((1)*(-muA + d1)), '(com e = 0): o polo de det(I + L^-1 E) em 0 e simples')
