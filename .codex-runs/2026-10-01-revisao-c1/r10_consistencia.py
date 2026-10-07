#!/usr/bin/env python3
"""R10 (revisao independente C1): consistencia entre os artefatos.
 - T3: lam0 = centro do nome, eps_fundo_L = 4,7e-8, rouche, MESMA particao de janelas (nome, modos, SZ) em todos;
 - rig: c e pontos = plano (build/rouche_L_est/plano_rig.txt, 6 casas); nBhat igual; resto_F/qTT;
 - t_p do rig >= t_p dos JSONs tp_* das maquinas (quando existem);
 - disco: c, certF registrado, pontos = rig, t_p = rig, ok;
 - certF: c, sha256_F = sha256 do F local, (quando ausente: SEM SHA), posto;
 - E.json, schur.json (defl_refA) iguais nas tres maquinas."""
import hashlib, json, glob
from pathlib import Path

base = Path('build/rouche_L_rig')
centros = {'-0_25': complex(0, -0.25), '0': 0j, '0_25': 0.25j, '0_55p0_25': 0.55+0.25j, '0_55-0_25': 0.55-0.25j,
           '1_0': 1+0j, '1_5': 1.5+0j, '2_0': 2+0j, '2_5': 2.5+0j, '3_0': 3+0j}
lab = lambda c: f'{c.real + 0.0:+.4f}{c.imag + 0.0:+.4f}i'.replace('.', '_')
plano = {}
for l in open('build/rouche_L_est/plano_rig.txt'):
    if not l.strip():
        continue
    a, b = l.split(maxsplit=1)
    plano[complex(a.replace('j', 'j'))] = [complex(x) for x in b.strip().split(',')]
prob = []
ok = lambda cond, msg: (None if cond else prob.append(msg))
ref = None
shaF = {}
for f in sorted((base/'F').glob('F_*.npy')):
    h = hashlib.sha256()
    with open(f, 'rb') as fh:
        for bl in iter(lambda: fh.read(1 << 24), b''):
            h.update(bl)
    shaF[f.name] = h.hexdigest()
tpdir = Path('.codex-runs/2026-10-01-revisao-c1/r10_copia/lab2_tp')
for rot, c in centros.items():
    t3 = json.load(open(base/f'cauda/T3_{rot}.json'))
    ok(complex(*t3['lam0']) == c, f'T3_{rot}: lam0 {t3["lam0"]} != {c}')
    ok(t3['eps_fundo_L'] == 4.7e-08, f'T3_{rot}: eps_fundo_L {t3["eps_fundo_L"]}')
    ok(t3['rouche'], f'T3_{rot}: rouche falso')
    part = (t3['nome'], t3['modos'], t3['SZ'], t3['folga_far'], t3['dB_L'])
    if ref is None:
        ref = part
    ok(part[:4] == ref[:4], f'T3_{rot}: particao de janelas diferente')
    ok(part[4] == ref[4], f'T3_{rot}: dB_L diferente')
    rigs = glob.glob(str(base/f'z/*/rig_{lab(c)}_rig.json'))
    ok(len(rigs) == 1, f'{rot}: {len(rigs)} rigs')
    rg = json.load(open(rigs[0])); maq = Path(rigs[0]).parent.name
    ok(complex(*rg['c']) == c, f'rig {rot}: c')
    pts = [complex(*x['p']) for x in rg['pontos']]
    pl = plano.get(c)
    ok(pl is not None and len(pl) == len(pts) and all(abs(complex(round(p.real, 6), round(p.imag, 6)) - q) < 1e-9 for p, q in zip(pts, pl)),
       f'rig {rot}: pontos != plano')
    ok(len({x['nBhat'] for x in rg['pontos']}) == 1 and rg['pontos'][0]['nBhat'] == 3.850313347957818, f'rig {rot}: nBhat')
    tpf = tpdir/f'tp_{lab(c)}_rig.json'
    if tpf.exists():
        tp = json.load(open(tpf))
        ok(all(abs(complex(*a) - b) < 1e-12 for a, b in zip(tp['pontos'], pts)), f'tp {rot}: pontos')
        ok(all(x['t_p'] >= y for x, y in zip(rg['pontos'], tp['t_p'])), f'rig {rot}: t_p menor que o do tp_*.json')
        dif = max(x['t_p'] - y for x, y in zip(rg['pontos'], tp['t_p']))
        print(f'  {rot}: t_p do rig - t_p de tp_*.json (maquina): max {dif:.1e} (>= 0 exigido)')
    d = json.load(open(base/f'discos/disco_{rot}.json'))
    ok(complex(*d['c']) == c, f'disco {rot}: c')
    ok(d.get('certF') and Path(d['certF']).name == f'certF_{lab(c)}.json', f'disco {rot}: certF {d.get("certF")}')
    ok([complex(*x['p']) for x in d['pontos']] == pts, f'disco {rot}: pontos != rig')
    ok(all(x['t_p'] == y['t_p'] for x, y in zip(d['pontos'], rg['pontos'])), f'disco {rot}: t_p != rig')
    ok(all(x['ok'] for x in d['pontos']), f'disco {rot}: ponto nao ok')
    cF = json.load(open(base/f'F/certF_{lab(c)}.json'))
    ok(complex(*cF['c']) == c, f'certF {rot}: c')
    s = cF.get('sha256_F')
    estado = 'SEM sha256_F' if s is None else ('sha OK' if s == shaF[f'F_{lab(c)}.npy'] else 'sha DIFERENTE')
    if s is not None:
        ok(s == shaF[f'F_{lab(c)}.npy'], f'certF {rot}: sha256 difere do F local')
    fj = json.load(open(base/f'F/F_{lab(c)}.json'))
    print(f'{rot:10s} rig em {maq:5s} ({len(pts)} pontos), certF: {estado}, posto {cF["posto"]}, k {fj.get("k")}, rJ {fj.get("rJ")}, '
          f'qTT cert {cF["qTT"]:.6f} (fator {fj["qTT"]:.6f}), resto_F antigo {fj["resto_F"]:.3e}, theta max {max(x["theta_f"] for x in d["pontos"]):.6f}')
sch = [json.load(open(base/f'z/{m}/schur.json')) for m in ('lab', 'lab2', 'joao')]
ok(sch[0]['defl_refA'] == sch[1]['defl_refA'] == sch[2]['defl_refA'], 'schur.json: defl_refA difere entre maquinas')
gl = [json.load(open(base/f'z/{m}/globais.json')) for m in ('lab', 'lab2', 'joao')]
ok(gl[0] == gl[1] == gl[2], 'globais.json difere')
print('raizes de schur.json por maquina:', [s['raizes_em_Omega'][0][0] for s in sch])
print('PROBLEMAS:' if prob else 'nenhuma inconsistencia', *prob, sep='\n  ')
