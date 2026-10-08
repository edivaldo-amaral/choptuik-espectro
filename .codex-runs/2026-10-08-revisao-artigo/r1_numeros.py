#!/usr/bin/env python3
"""R1 (revisao independente de 08/10/2026): auditoria dos numeros certificados impressos no manuscrito.

Le SO os arquivos de certificado (JSON/TXT em build/, dados de RT em .cache/) e compara com os numeros impressos
em paper/sec/*.tex. Para cada numero: impresso | valor exato do arquivo | tipo (sup = cota superior, inf = cota
inferior, ~ = valor aproximado/descritivo) | seguro? (sup: impresso >= exato; inf: impresso <= exato; ~: concorda
nos digitos impressos). Nao escreve nada fora desta pasta.

Uso: .venv/bin/python .codex-runs/2026-10-08-revisao-artigo/r1_numeros.py > r1_saida.txt
"""
from __future__ import annotations

from fractions import Fraction as Fr
import glob
import json
import math
import re
from pathlib import Path

import flint
from flint import arb, ctx

RAIZ = Path(__file__).resolve().parents[2]
ctx.prec = 400
J = lambda p: json.loads((RAIZ/p).read_text())
linhas = []


def reg(nome, impresso, exato, tipo, local, nota=''):
    """tipo: 'sup' (impresso deve ser >= exato), 'inf' (impresso <= exato), '~' (aproximacao),
    'eq' (deve coincidir), 'lt' (exato < limiar impresso: 'X < c')."""
    imp = Fr(impresso) if not isinstance(impresso, Fr) else impresso
    ex = Fr(exato) if not isinstance(exato, Fr) else exato
    if tipo == 'sup':
        ok = imp >= ex
    elif tipo == 'inf':
        ok = imp <= ex
    elif tipo == 'lt':
        ok = ex < imp
    elif tipo == 'eq':
        ok = imp == ex
    else:  # '~': concorda ate a ultima casa impressa (meia unidade)
        casas = len(str(impresso).split('.')[-1]) if '.' in str(impresso) else 0
        ok = abs(imp - ex) <= Fr(1, 2*10**casas) + Fr(1, 10**(casas + 3))
    linhas.append((nome, str(impresso), f'{float(ex):.10g}', tipo, 'SIM' if ok else 'NAO', local, nota))
    return ok


# ---------------------------------------------------------------- dados de RT (independente: parser proprio)
def le_rt(arq):
    txt = (RAIZ/arq).read_text()
    m = re.search(r'mu_MultiField\s+DyadicQ\s+TwoExp\s+(-?\d+)\s+coeff\s+(-?\d+)', txt)
    mu = Fr(int(m.group(2)))*Fr(2)**int(m.group(1))
    campos = []
    for bloco in re.split(r'(?m)^Field\s*$', txt)[1:]:
        h = {k: int(re.search(rf'(?m)^{k}\s+(-?\d+)\s*$', bloco).group(1))
             for k in ('TwoExp', 'off_m', 'off_n', 'num_m', 'num_n')}
        pares = re.findall(r'\((-?\d+),(-?\d+)\)', bloco)
        assert len(pares) == h['num_m']*h['num_n'], (arq, len(pares))
        campos.append((h, pares))
    return mu, campos


muA, camposA = le_rt('.cache/rt-1203.3766v1/sourcecode/RefA.dat')
with open(RAIZ/'.cache/rt-1203.3766v1/RefAplusB.dat') as f:
    cab = f.read(4000)
m_ = re.search(r'mu_MultiField\s+DyadicQ\s+TwoExp\s+(-?\d+)\s+coeff\s+(-?\d+)', cab)
muAB = Fr(int(m_.group(2)))*Fr(2)**int(m_.group(1))
# constantes de RT (choptuik.tex, eq. dnkdfkjhfkdhfkdhdfdk1/2): |K - K80| <= 1e-80, |mu - mu80| <= 1e-80
K80 = Fr(int('17227262011139106749857559727881918622210' '3458805781088634788403570540271455648855'), 10**80)
MU80 = Fr(int('1683070789634499695101349790428574207210' '0199080892966476395293134873313662587505'), 10**80)
reg('mu_A (RefA.dat) = 722873400*2^-32', '0.16830707900226116', muA, '~', 'roots.tex:57 (lambda0 de mu)')
reg('|mu_{A+B} - mu80| <= 1e-80 (consistencia RT)', '1e-80', abs(muAB - MU80), 'sup', 'setting.tex:18')
reg('mu = 0.168307078963... (enunciado)', '0.168307078963', MU80, '~', 'main_theorem.tex:41')
reg('mu ~ 0.16830708', '0.16830708', MU80, '~', 'setting.tex:18')
reg('mu ~ 0.16831', '0.16831', MU80, '~', 'main_theorem.tex:6')
reg('K ~ 1.72272620', '1.72272620', K80, '~', 'setting.tex:18')
reg('Delta = 2K = 3.44545240...', '3.44545240', 2*K80, '~', 'setting.tex:51 (truncado)')
reg('|mu - mu_A| (usado nos tiles: 3.9e-11)', '3.9e-11', abs(MU80 - muA) + Fr(1, 10**80), 'sup', 'scripts/check_cover.py (nao impresso)')

# ---------------------------------------------------------------- NK (Prop. 6.2)
NK = {'mu': ('build/s4/rig', dict(theta='0.880172', Z2='183.66', Y='6.60e-6', r='6.07e-5', erro='2.82e-5')),
      'sK': ('build/s3b/rig', dict(theta='0.839937', Z2='42.97', Y='5.45e-5', r='3.79e-4', erro='1.72e-4')),
      'sphys': ('build/nk733', dict(theta='0.805721', Z2='10.33', Y='4.19e-6', r='2.16e-5', erro='1.0163e-5')),
      'sB': ('build/faixaB/nk', dict(theta='0.91153', Z2='155.02', Y='5.97e-6', r='7.81e-5', erro='3.67e-5'))}
C3 = {}
for raiz, (pasta, imp) in NK.items():
    d = J(f'{pasta}/C3.json'); C3[raiz] = d
    for k, chave in (('theta', 'theta'), ('Z2', 'Z2'), ('Y', 'Y'), ('r', 'r'), ('erro', 'erro_lambda')):
        reg(f'NK {raiz}: {k}', imp[k], Fr(d[chave]), 'sup', 'roots.tex:57-60')
    th, Z2, Y, r = (Fr(d[k]) for k in ('theta', 'Z2', 'Y', 'r'))
    # condicoes do teorema NK como enunciado no artigo (com o Z2 gravado no papel de Z2)
    cond = (1 - th)**2 > 4*Z2*Y and Y + th*r + Z2*r*r <= r and th + 2*Z2*r < 1
    linhas.append((f'NK {raiz}: (1-th)^2>4 Z2 Y, Y+th r+Z2 r^2<=r, th+2Z2 r<1', '-', str(cond), 'cond', 'SIM' if cond else 'NAO',
                   'roots.tex:37-40', f'raio de unicidade (1-th)/(2Z2) = {float((1 - th)/(2*Z2)):.4e}'))
lam = {k: (complex(*C3[k]['lambda0']) if isinstance(C3[k]['lambda0'], list) else complex(C3[k]['lambda0'])) for k in C3}
reg('lambda0 sK', Fr(0.401024733), Fr(C3['sK']['lambda0']), 'eq', 'roots.tex:58')
reg('lambda0 sphys', Fr(0.7331796596606924), Fr(C3['sphys']['lambda0']), 'eq', 'roots.tex:59')
reg('lambda0 sB (Re)', Fr(0.0414194105346554), Fr(C3['sB']['lambda0'][0]), 'eq', 'roots.tex:60')
reg('lambda0 mu = mu_A', muA, Fr(C3['mu']['lambda0']), 'eq', 'roots.tex:57')
# enclausuramentos do Teorema Principal
eB = Fr(C3['sB']['erro_lambda']) + abs(Fr(C3['sB']['lambda0'][0]) - Fr('0.0414194105'))
reg('|sB - 0.0414194105| <= 3.7e-5', '3.7e-5', eB, 'sup', 'main_theorem.tex:42')
reg('|sK - 0.401024733| <= 1.72e-4', '1.72e-4', Fr(C3['sK']['erro_lambda']), 'sup', 'main_theorem.tex:43')
l0, e0 = Fr(C3['sphys']['lambda0']), Fr(C3['sphys']['erro_lambda'])
reg('sphys >= 0.73316949', '0.73316949', l0 - e0, 'inf', 'main_theorem.tex:44')
reg('sphys <= 0.73318983', '0.73318983', l0 + e0, 'sup', 'main_theorem.tex:44')
reg('raio dos discos de sK, sphys < 2e-4', '2e-4', max(Fr(C3['sK']['erro_lambda']), e0), 'lt', 'roots.tex:168')

# ---------------------------------------------------------------- lambda e 1/lambda em Arb (K de RT)
Kb = arb(flint.fmpq(K80.numerator, K80.denominator)) + arb(0, flint.fmpq(1, 10**80))
lo = arb(flint.fmpq((l0 - e0).numerator, (l0 - e0).denominator))
hi = arb(flint.fmpq((l0 + e0).numerator, (l0 + e0).denominator))
pi = arb.pi()
lam_lo = 2*pi*lo/Kb; lam_hi = 2*pi*hi/Kb
inv_lo = Kb/(2*pi*hi); inv_hi = Kb/(2*pi*lo)


def sup_arb(x):   # cota superior racional de uma bola
    m, r = x.mid(), x.rad()
    return Fr(str((m + r).upper().str(30, radius=False)).replace('e', 'e')) if False else Fr(float((m + r).upper()))


def lim(x):
    """extremos (inferior, superior) de uma bola arb como Fractions (via strings de 40 digitos, folga 1e-35)."""
    a = x.lower(); b = x.upper()
    fa = Fr(a.str(40, radius=False).replace(' ', '')); fb = Fr(b.str(40, radius=False).replace(' ', ''))
    return fa - Fr(1, 10**35), fb + Fr(1, 10**35)


reg('lambda >= 2.674040', '2.674040', lim(lam_lo)[0], 'inf', 'main_theorem.tex:69; intro:52; resumo')
reg('lambda <= 2.674115', '2.674115', lim(lam_hi)[1], 'sup', 'main_theorem.tex:69; intro:52; resumo')
reg('1/lambda >= 0.373955', '0.373955', lim(inv_lo)[0], 'inf', 'main_theorem.tex:70; intro:52; resumo')
reg('1/lambda <= 0.373966', '0.373966', lim(inv_hi)[1], 'sup', 'main_theorem.tex:70; intro:52; resumo')
# a partir do intervalo impresso [0.73316949, 0.73318983] (como um leitor faria)
lam_lo2 = 2*pi*arb('0.73316949')/Kb; lam_hi2 = 2*pi*arb('0.73318983')/Kb
linhas.append(('lambda a partir do intervalo impresso de s*', '-', f'[{lam_lo2.lower()}, {lam_hi2.upper()}]', '~', '-',
               'main_theorem.tex:44', 'so informativo'))

# ---------------------------------------------------------------- simbolos (Lema 5.2, Teorema 5.1(i))
sL = J('build/spectrum/symbol-L-2048x4096.json')
reg('max Re W(-B) <= 2.9492', '2.9492', Fr(sL['R_fundo']), 'sup', 'counting.tex:58')
reg('||B|| <= 3.9642', '3.9642', Fr(sL['norma_fundo']), 'sup', 'counting.tex:58; setting.tex:309; app A.3')
faixas = [J(f) for f in sL['faixas']]
folgaL = max(Fr(f['folga_lambda']) for f in faixas)
reg('folga de Lipschitz de L "0.0771" (o texto chama de minima; e a MAXIMA somada)', '0.0771', folgaL, 'sup',
    'counting.tex:68', 'a folga entra somada na cota; o seguro e arredondar para cima')
fp = J('build/spectrum/free-part-toro.json')
RJ = max(Fr(v['R_livre']) for v in fp['L']['livres'].values())
reg('max Re W(-J) <= -0.0232 (certificado gravado, mu >= 1/6)', '-0.0232', RJ, 'sup', 'counting.tex:59',
    'com mu verdadeiro: -mu*lambda0 = %.6f' % float(-MU80*Fr(fp['L']['livres']['1']['lambda0_mu1'])))
RJmu = max(-MU80*Fr(v['lambda0_mu1']) for v in fp['L']['livres'].values())
reg('max Re W(-J) <= -0.0232 (com mu verdadeiro)', '-0.0232', RJmu, 'sup', 'counting.tex:59')
reg('L(s) invertivel para Re s > 2.926 (certificado gravado)', '2.926', Fr(sL['R_fundo']) + RJ, 'sup', 'counting.tex:14,74')
reg('L(s) invertivel para Re s > 2.926 (mu verdadeiro)', '2.926', Fr(sL['R_fundo']) + RJmu, 'sup', 'counting.tex:14,74')
sK = J('build/spectrum/symbol-K-1024x2048.json')
reg('max Re W(-B#) <= 1.4396', '1.4396', Fr(sK['R_fundo']), 'sup', 'counting.tex:59')
reg('||B#|| <= 1.6497', '1.6497', Fr(sK['norma_fundo']), 'sup', 'counting.tex:60')
RJK = max(Fr(v['R_livre']) for v in fp['K']['livres'].values())
reg('max Re W(-K(0)) <= 1.765', '1.765', Fr(sK['R_fundo']) + RJK, 'sup', 'counting.tex:60; roots.tex:132')
reg('folga de Lipschitz de K "0.0535"', '0.0535', Fr(sK['folga_lambda']), 'sup', 'counting.tex:68')
linhas.append(('grade do simbolo de K', '2048x4096 (texto e Tabela 1)', f"{sK['nphi']}x{sK['ntheta']}", 'eq', 'NAO',
               'counting.tex:64; computation.tex:48', 'a grade de K gravada e 1024x2048'))

# ---------------------------------------------------------------- discos (Rouche)
def discos(pat):
    out = []
    for f in sorted(glob.glob(str(RAIZ/pat))):
        for x in J(f)['pontos']:
            out.append((x['ok'], float(x['r']), Fr(x['theta']), f))
    return out


dA = discos('build/rouche_L_rig/discos_v2/disco_*.json'); dB = discos('build/faixaB/discos/disco_*.json')
nA = sum(1 for o, r, t, f in dA if o and r > 0); nB = sum(1 for o, r, t, f in dB if o and r > 0)
reg('discos da faixa A', '96', nA, 'eq', 'counting.tex:206')
reg('discos adicionais da faixa B', '78', nB, 'eq', 'counting.tex:206')
reg('pontos nos arquivos de discos', '175', len(dA) + len(dB), 'eq', 'computation.tex:72')
reg('arquivos de discos', '22', len({f for *_, f in dA + dB}), 'eq', 'computation.tex:72')
reg('max theta faixa A', '0.99842', max(t for o, r, t, f in dA if o), 'sup', 'counting.tex:206')
reg('max theta faixa B (discos proprios)', '0.99908', max(t for o, r, t, f in dB if o), 'sup', 'counting.tex:207')
reg('pior Perron nos discos', '0.999074', max(t for o, r, t, f in dA + dB if o), 'sup', 'computation.tex:82')
rs = [r for o, r, t, f in dA + dB if o and r > 0]
reg('menor raio', '0.0028', Fr(min(rs)), 'inf', 'counting.tex:208', 'descritivo')
reg('maior raio', '0.4995', Fr(max(rs)), 'sup', 'counting.tex:208', 'descritivo')
dif = max(J(f).get('max_dtheta', 0) for f in glob.glob(str(RAIZ/'build/reproducao/discos/*_disco_*.json')))
reg('diferenca nivel 1 vs gravado', '7.8e-16', Fr(dif), 'sup', 'computation.tex:83')

# ---------------------------------------------------------------- janelas
jA = J('build/rouche_L_rig/janelas_v2.json'); jB = J('build/faixaB/janelas.json')
reg('max rho_J faixa A', '0.1488', max(Fr(v['max_rho']) for v in jA.values()), 'sup', 'counting.tex:229')
reg('max rho_J faixa B', '0.1517', max(Fr(v['max_rho']) for v in jB.values()), 'sup', 'counting.tex:229')
reg('retangulos por faixa', '10', len(jA), 'eq', 'counting.tex:226')
reg('retangulos por faixa (B)', '10', len(jB), 'eq', 'counting.tex:226')

# ---------------------------------------------------------------- contagem finita
ct = J('build/rouche_L_rig/contagem.json')
reg('||AU - UT||_F <= 5.97e-5', '5.97e-5', Fr(ct['eps_R']), 'sup', 'counting.tex:250')
reg('||I - U*U||_F <= 1.30e-6', '1.30e-6', Fr(ct['eps_G']), 'sup', 'counting.tex:251')
reg('eps_B gravado >= eps_R/sqrt(1-eps_G) (consistencia)', Fr(ct['eps_B']), Fr(ct['eps_R'])/Fr(math.sqrt(1 - ct['eps_G'])), 'sup', 'contagem.json')
txt = (RAIZ/'build/rouche_L_rig/contagem_rouche_v2.txt').read_text()
tmax = Fr(re.search(r'max t\(s\) no contorno <= ([0-9.]+)', txt).group(1))
teps = Fr(re.search(r't eps_B <= ([0-9.e+-]+)', txt).group(1))
reg('sup t <= 559.2', '559.2', tmax, 'sup', 'counting.tex:251', 'valor do TXT (3 casas)')
reg('t eps_B <= 0.0334', '0.0334', teps, 'sup', 'counting.tex:252', 'valor do TXT (4 alg.)')
reg('t eps_B recomputado = 559.151*eps_B', '0.0334', tmax*Fr(ct['eps_B']), 'sup', 'counting.tex:252')
emO = sorted(complex(*z).real for z in ct['em_Omega'])
reg('-T_ii em Omega_A: 0.401025...', '0.401025', Fr(emO[0]), '~', 'counting.tex:253')
reg('-T_ii em Omega_A: 0.733180...', '0.733180', Fr(emO[1]), '~', 'counting.tex:253')
reg('distancia da entrada mais proxima fora de Omega_A', '0.168', min(Fr(z[2]) for z in ct['fora_perto']), 'inf', 'counting.tex:254')

# ---------------------------------------------------------------- deflacao
sch = J('build/rouche_L_rig/z/lab2/schur.json')['defl_refA']
E = J('build/s4/rig/E.json')
reg('||P_Z(x0 - x0^A)|| <= 1.5e-9 (dtauB = 1.4563e-9 + arred.)', '1.5e-9', Fr('1.4563e-9') + Fr(1, 10**60) + Fr(2**-52)*10*Fr(sch['norma_x'][0]), 'sup', 'counting.tex:127')
reg('||x1 - x1^A|| <= 2.51e-8', '2.51e-8', Fr(E['e']), 'sup', 'counting.tex:128; roots.tex:92')
reg('||phi_i|| <= 5.5 (implicito em |e_ij| <= 5.5*2.51e-8)', '5.5', max(Fr(x) for x in sch['norma_phi']), 'sup', 'counting.tex:128')
eij = max(Fr(x) for x in sch['norma_phi'])*Fr(E['e'])
mod = 1 - eij - Fr(1, 10**9)
reg('fatores de q com modulo >= 1 - 1e-6', '0.999999', mod, 'inf', 'counting.tex:129')

# ---------------------------------------------------------------- identificacao de mu
reg('||L(lambda0) h0|| <= 3.76e-6', '3.76e-6', Fr(E['rL']), 'sup', 'roots.tex:94')
reg('||x_g - x0|| <= 1.82e-5', '1.82e-5', Fr(E['dist']), 'sup', 'roots.tex:94')
reg('r (E.json) <= 6.07e-5', '6.07e-5', Fr(E['r']), 'sup', 'roots.tex:94')
reg('dist < r (identificacao)', E['r'], Fr(E['dist']), 'lt', 'roots.tex:94')

# ---------------------------------------------------------------- testemunhos
for nome, pasta, imp in (('sK', 'build/s3b/rig', ('0.16975', '0.17173', '1.99e-3')),
                         ('sB', 'build/faixaB/nk', ('0.35469', '0.35518', '4.96e-4'))):
    D = J(f'{pasta}/D.json')
    reg(f'testemunho {nome}: ||Pi_F C h|| >= {imp[0]}', imp[0], Fr(D['folga']), 'inf', 'roots.tex:103')
    reg(f'testemunho {nome}: c_ref >= {imp[1]}', imp[1], Fr(D['c_ref']), 'inf', 'roots.tex:114-115')
    reg(f'testemunho {nome}: perdas < {imp[2]}', imp[2], Fr(D['perda']), 'lt', 'roots.tex:114-115')
    reg(f'testemunho {nome}: folga = c_ref - perda', D['folga'], Fr(D['c_ref']) - Fr(D['perda']), '~', 'D.json')
    linhas.append((f'caixa F do testemunho {nome}', '{|m|<12, n<48}', str(D['F']), 'eq', 'SIM' if D['F'] == [12, 48] else 'NAO', 'roots.tex:105', ''))

# ---------------------------------------------------------------- K: tiles e Rouche
rk = J('build/s3b/rouche_K_circ.json'); rkc = J('build/s3b/rouche_K_contagem_20x80.json')
reg('Rouche de K: raio do circulo', Fr(0.05), Fr(rk['rho']), 'eq', 'roots.tex:126,137')
reg('Rouche de K: Perron <= 0.8126', '0.8126', Fr(rk['pior_theta_rouche']), 'sup', 'roots.tex:139')
reg('Rouche de K: tiles no circulo', '16', len(rk['tiles']), 'eq', 'roots.tex:139')
linhas.append(('Rouche de K: caixa Z', '20x80', str(rk['Z']), 'eq', 'SIM' if rk['Z'] == [20, 80] else 'NAO', 'roots.tex:139', ''))
reg('Rouche de K: autovalor ~ 0.4010247311', '0.4010247311', Fr(complex(rkc['autovalores'][0].strip('()').replace('j', 'j')).real), '~', 'roots.tex:141')
cov = (RAIZ/'.codex-runs/2026-10-08-revisao-artigo/r6_cobertura_disco005.txt').read_text()
reg('tiles certificados', '154', int(re.search(r'(\d+) discos certificados', cov).group(1)), 'eq', 'roots.tex:133')
reg('pior Perron dos tiles <= 0.99612', '0.99612', Fr(re.search(r'pior theta reverificado: ([0-9.]+)', cov).group(1)), 'sup', 'roots.tex:135', 'reverificado por check_cover (7 casas)')

# ---------------------------------------------------------------- L-real e recuperacao de raio
lr = J('build/t2/l_real_toro.json'); rt_ = J('build/t2/radius_transfer_s3_1.json')
reg('||QT||_toro <= 0.0525', '0.0525', Fr(lr['qT']), 'sup', 'app_paper_lemmas.tex:69')
reg('defeito toro < 0.371 (gravado, com nB = %.7f)' % lr['nB'], '0.371', Fr(lr['defeito']), 'lt', 'setting.tex:308')
reg('defeito toro < 0.371 com ||B|| certificado 3.9641333', '0.371', (Fr(sL['norma_fundo']) + Fr(31, 10))*Fr(lr['qT']), 'lt', 'setting.tex:308')
reg('nB usado em l_real_toro >= ||B|| certificado', lr['nB'], Fr(sL['norma_fundo']), 'sup', 'build/t2/l_real_toro.json',
    'o certificado gravado usou ||B|| = 3.9641+1e-6, abaixo da cota do simbolo')
reg('defeito forte em Y+ <= 0.2665', '0.2665', Fr(rt_['strong_defect']), 'sup', 'setting.tex:309; app_paper_lemmas.tex:72')
reg('defeito fraco (radius_transfer) < 1', '1', Fr(rt_['weak_defect']), 'lt', 'radius_transfer_s3_1.json')
r5 = J('build/t2/radius_recovery_geral.json')['linhas']
ex = [l for l in r5 if l['S'] == '200' and l['k2'] == '65/64'][0]
linhas.append(('caixa de R5 para S=200, k2=1+2^-6', '2^17 x 2^20', f"{ex['M']} x {ex['N']}", 'eq',
               'SIM' if (ex['M'], ex['N']) == (2**17, 2**20) else 'NAO', 'physical_modes.tex:234; app A.3', ''))
linhas.append(('R5: todas as linhas contraem', '-', str(all(l['contrai'] for l in r5)), 'eq', 'SIM' if all(l['contrai'] for l in r5) else 'NAO', 'radius_recovery_geral.json', ''))
# derivadas do Lema 2.4 (dominios)
qa = Fr(129, 128)/Fr(65, 64); qb = Fr(9, 8)/Fr(5, 4)
reg('||d_tau|| <= 8385', '8385', qa/(2*(1 - qa)**2), 'sup', 'app_paper_lemmas.tex:112')
reg('||d_xi|| <= 1440', '1440', 2*qb/((Fr(9, 8) - 1)*(1 - qb)**2), 'sup', 'app_paper_lemmas.tex:112')

# ---------------------------------------------------------------- F1 (recalculo independente a partir de RefA.dat)
h, pares = camposA[3]
esc = Fr(2)**h['TwoExp']; i = 1 - h['off_m']
c0 = [Fr(0), Fr(0)]; c1 = [Fr(0), Fr(0)]
for k in range(h['num_n']):
    n = h['off_n'] + k
    a, b = (int(x) for x in pares[i*h['num_n'] + k])
    f0 = (1 if n == 0 else 2)*[1, 0, -1, 0][n % 4]          # i^n + i^-n
    f1 = (1 if n == 0 else 2)*(-1)**n
    c0[0] += f0*a*esc; c0[1] += f0*b*esc; c1[0] += f1*a*esc; c1[1] += f1*b*esc
eps = Fr(1, 2**25) + Fr(1, 2**277)
reg('F1: Re c1(omega_A,4(.,0)) = 0.072916...', '0.072916', c0[0], '~', 'physical_modes.tex:100', 'truncado')
reg('F1: Im c1(omega_A,4(.,0)) >= 0.493046', '0.493046', c0[1], 'inf', 'physical_modes.tex:93,100')
reg('F1: Re c1(.,-1) = -0.079473...', '-0.079473', c1[0], '~', 'physical_modes.tex:100')
reg('F1: Im c1(.,-1) >= 0.097753', '0.097753', c1[1], 'inf', 'physical_modes.tex:94,100')
reg('F1: 2^-25 + 2^-277 < 3e-8', '3e-8', eps, 'lt', 'physical_modes.tex:101')
g = J('build/t2/gauge_global_centro.json')
reg('F1: c1_re igual ao gravado', g['c1_re'], c0[0], 'eq', 'gauge_global_centro.json')
reg('F1: c1_im igual ao gravado', g['c1_im'], c0[1], 'eq', 'gauge_global_centro.json')
reg('F1: c1_cone_re igual ao gravado', g['c1_cone_re'], c1[0], 'eq', 'gauge_global_centro.json')
reg('F1: c1_cone_im igual ao gravado', g['c1_cone_im'], c1[1], 'eq', 'gauge_global_centro.json')

# ---------------------------------------------------------------- saida
w = [max(len(str(l[j])) for l in linhas) for j in range(7)]
print('| grandeza | impresso | exato (arquivo) | tipo | seguro? | local | nota |')
print('|---|---|---|---|---|---|---|')
for l in linhas:
    print('| ' + ' | '.join(str(x) for x in l) + ' |')
nao = [l for l in linhas if l[4] == 'NAO']
print(f'\n{len(linhas)} linhas; {len(nao)} com "NAO":')
for l in nao:
    print('  -', l[0], '| impresso', l[1], '| exato', l[2], '|', l[5])
