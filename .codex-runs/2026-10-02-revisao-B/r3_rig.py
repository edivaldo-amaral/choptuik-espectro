#!/usr/bin/env python3
"""R3 (revisor): conferencia dos niveis Z da faixa B (rig_*_B*.json e rig_*_Bteste.json).

Para cada rig: sha256 de A = 70791ce0...; sha256 de F = o do certF do centro; versao 2; pontos = os do tp_*.json
correspondente (copiado do laboratorio2 / PC 3) e os do disco que o usa; t_p do rig >= o t_p da FASE 1
impresso no log (linhas '(estimativa ...)', nao as '(do log anterior)'), com folga relativa <= 1e-6 (o rig
aplica ainda a correcao de deslocamento arredondado, que so aumenta). Os logs da fase 1: rig_B.txt (lab2 e Joao)
e build/caudaB/teste_reuso.txt (lab2), copiados para o scratchpad (caminhos em argv)."""
import glob, json, re, sys
from pathlib import Path

RAIZ = Path('<repo>')
SCR = Path(sys.argv[1])
SHA_A = '70791ce0c30f8065f93770c36ff89c1d8290978c5abf84ab0570313bc09fd002'
logs = [SCR/'lab2/rig_B.txt', SCR/'joao/rig_B.txt', SCR/'lab2/teste_reuso.txt']
fase1 = {}      # (prefixo F, p) -> [valores]
reuso = {}
pat1 = re.compile(r'^p = (.+?): t_p = \|\|\(A \+ p\)\^-1\|\| <= ([0-9.]+)\s+\((estimativa|do log anterior)')
for lg in logs:
    Fp = None
    for linha in lg.read_text().splitlines():
        m = re.match(r'sha256 A ([0-9a-f]+)\.\.\., F ([0-9a-f]+)\.\.\.', linha)
        if m:
            assert SHA_A.startswith(m.group(1)), lg
            Fp = m.group(2); continue
        m = pat1.match(linha)
        if m:
            p = complex(m.group(1).replace(' ', ''))
            dst = fase1 if m.group(3) == 'estimativa' else reuso
            dst.setdefault((Fp, p), []).append(float(m.group(2)))
certs = {}
for f in glob.glob(str(RAIZ/'build/faixaB/F/certF*.json')) + glob.glob(str(RAIZ/'build/rouche_L_rig/F/certF*.json')):
    d = json.load(open(f)); certs[complex(*d['c'])] = (f, d['sha256_F'])
# discos que usam cada rig
uso = {}
for f in glob.glob(str(RAIZ/'build/faixaB/discos/*.json')):
    d = json.load(open(f)); uso.setdefault(Path(d['zrig']).name, []).append((f, d))
tp_rem = {}
for f in glob.glob(str(SCR/'lab2/tp_*.json')) + glob.glob(str(SCR/'joao/tp_*.json')):
    tp_rem.setdefault(Path(f).name, []).append(json.load(open(f)))
problemas = []
saida = []
rigs = sorted(glob.glob(str(RAIZ/'build/faixaB/z/rig_*.json'))) + sorted(glob.glob(str(RAIZ/'build/faixaB/teste/rig_*.json')))
for f in rigs:
    d = json.load(open(f)); c = complex(*d['c']); nome = Path(f).name
    linha = dict(rig=nome, c=[c.real, c.imag], versao=d.get('versao'))
    ok_A = d['sha256_A'] == SHA_A
    fc, shaF = certs[c]
    ok_F = d['sha256_F'] == shaF
    Fp = d['sha256_F'][:16]
    pts = [complex(*x['p']) for x in d['pontos']]
    # tp: B2 e o _B do Joao (mesmo nome do tp)
    rot = nome.replace('rig_', 'tp_').replace('_B2.json', '_B.json')
    tps = tp_rem.get(rot, [])
    ok_tp_pts = any([complex(*q) for q in t['pontos']] == pts for t in tps)
    max_rel = 0.0; faltam = []; abaixo = []; dif_log = 0.0; tp_incoerente = []
    for x in d['pontos']:
        p = complex(*x['p'])
        v = fase1.get((Fp, p))
        if not v:
            faltam.append(p); continue
        v1 = max(v)
        # o log imprime 6 casas (arredondamento ao mais proximo): so e erro se o rig ficar abaixo de v1 - 5e-7
        if x['t_p'] < v1 - 5e-7*(1 + 1e-9):
            abaixo.append((p, x['t_p'], v1))
        dif_log = max(dif_log, abs(x['t_p'] - v1))
        # contra o tp_*.json (valores passados em --tps, precisao total): t_p(tp) <= t_p(rig) <= t_p(tp)(1 + 1e-6)
        for t in tps:
            for q, tv in zip(t['pontos'], t['t_p']):
                if complex(*q) == p and not (tv <= x['t_p'] <= tv*(1 + 1e-6)):
                    tp_incoerente.append((p, tv, x['t_p']))
        max_rel = max(max_rel, x['t_p']/v1 - 1)
    # discos
    dd = uso.get(nome, [])
    ok_disco = True; npts_disco = 0
    for fd, ds in dd:
        assert ds['sha256_F'] == d['sha256_F'] and ds['sha256_A'] == d['sha256_A']
        mp = {complex(*x['p']): x['t_p'] for x in d['pontos']}
        for x in ds['pontos']:
            npts_disco += 1
            if complex(*x['p']) not in mp or mp[complex(*x['p'])] != x['t_p']:
                ok_disco = False
    linha.update(sha_A=ok_A, sha_F_certF=ok_F, certF=Path(fc).name, n_pontos=len(pts), tp_json_mesmos_pontos=ok_tp_pts,
                 n_tp_json=len(tps), fase1_faltando=[str(p) for p in faltam], tp_abaixo_fase1=[str(a) for a in abaixo],
                 folga_rel_max=max_rel, dif_abs_log=dif_log, tp_incoerente=[str(a) for a in tp_incoerente], discos=[Path(fd).name for fd, _ in dd], pontos_discos=npts_disco, discos_coerentes=ok_disco)
    saida.append(linha)
    ruim = (not ok_A) or (not ok_F) or d.get('versao') != 2 or faltam or abaixo or max_rel > 1e-6 or not ok_tp_pts or not ok_disco or tp_incoerente
    if ruim:
        problemas.append(nome)
    print(f'{nome:40s} c={c}  shaA {"OK" if ok_A else "ERRO"}  shaF=certF {"OK" if ok_F else "ERRO"}  v{d.get("versao")}  '
          f'{len(pts)} pts  tp.json {"OK" if ok_tp_pts else "DIFERE/AUSENTE"}  fase1: faltam {len(faltam)}, abaixo {len(abaixo)}, '
          f'|rig - log| max {dif_log:.1e}, rig vs tp.json {"OK" if not tp_incoerente else "INCOERENTE"}  discos {[Path(fd).name for fd, _ in dd]} ({npts_disco} pts, {"coerentes" if ok_disco else "INCOERENTES"})')
# pontos dos discos sem rig
print('problemas:', problemas if problemas else 'nenhum')
json.dump(dict(rigs=saida, problemas=problemas), open(RAIZ/'.codex-runs/2026-10-02-revisao-B/r3_resultado.json', 'w'), indent=1)
