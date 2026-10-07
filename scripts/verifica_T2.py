#!/usr/bin/env python3
"""Verificador de NIVEL 1 da cadeia de T2 (consolidacao, opcao C): refaz, a partir dos certificados gravados (arquivos
pequenos do repositorio), as verificacoes finais em aritmetica racional, sem recalcular as matrizes pesadas.

Nivel 1 (este script): para cada passo de docs/T2_ENUNCIADO.md,
  discos     : theta racional de CADA ponto de disco (faixas A e B) recalculado com o mesmo codigo do certificado
               (rouche_L_disco.contexto/theta_disco), igual ao gravado e < 1;
  cobertura  : a uniao dos discos cobre o contorno (rouche_L_cobertura.py, racional);
  contagem   : t * eps_B < 1 e a contagem de Schur (rouche_L_contagem.py rouche);
  nk         : o Perron de 3 niveis em racionais (nk_L3_C.py) refeito a partir de Z.json, Q.json e T3.json;
  s3a        : Perron de cada tile e cobertura em racionais (check_cover.py);
  lemas      : F1/F1', R5 e a transferencia de raio refeitos (gauge_global_centro, radius_recovery_geral,
               radius_transfer).
Nivel 2 (docs/REPRODUCAO.md): recalcular os certificados intermediarios a partir dos dados de RT (horas a dias).

Uso: .venv/bin/python scripts/verifica_T2.py [--so discos_A,discos_B,...] ; relatorio em build/reproducao/nivel1.txt."""
from __future__ import annotations

import argparse
from fractions import Fraction as Fr
import glob
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

RAIZ = Path(__file__).resolve().parent.parent
import os
# erro do fundo de RT em L: ||B_L(RefB)|| <= 4,67e-8 (S4_GAUGE.md, Young sobre RefAplusB - RefA), usado em TODOS os
# certificados de L (S4, S3b, faixas A e B, NK). Sem a variavel, fundo_L usa a cota generica 40 eps_omega, maior.
os.environ.setdefault('CHOPTUIK_EPS_FUNDO_L', '4.7e-08')
sys.path.insert(0, str(RAIZ/'scripts'))
PY = str(RAIZ/'.venv/bin/python')

DIR_SCHUR = 'build/rouche_L_rig/z/lab2'          # schur.json: dimensoes Z/F e dados da deflacao (defl_refA)


def nome_c(c, pasta):
    """arquivo de cauda T3 do centro c (convencao de nomes das pastas de cauda)."""
    re_, im_ = c
    def f(x):
        s = f'{x:g}'.replace('.', '_')
        return s
    if im_ == 0:
        cand = [f'T3_{f(re_)}.json', f'T3_{re_:.1f}'.replace('.', '_') + '.json']
    elif re_ == 0:
        cand = [f'T3_{f(im_)}.json']
    else:
        sep = 'p' if im_ > 0 else '-'
        cand = [f'T3_{f(re_)}{sep}{f(abs(im_))}.json']
    return [str(Path(pasta)/x) for x in cand]


_CACHE = {}


def _rd_com_cache():
    """rouche_L_disco com a cota do fundo (igual para todos os discos) e os dados centrais da cauda em cache."""
    import rouche_L_disco as rd
    if not _CACHE.get('ok'):
        orig_dB = rd.fundo_L.dB_L
        def dB_cache(campos, banda):
            k = ('dB', banda)
            if k not in _CACHE:
                _CACHE[k] = orig_dB(campos, banda)
            return _CACHE[k]
        rd.fundo_L.dB_L = dB_cache
        orig_cent = rd.dados_centrais
        def cent_cache(t3, c, MZ, NZ, Nc):
            import hashlib
            k = ('cent', hashlib.md5(json.dumps(t3, sort_keys=True).encode()).hexdigest(), c, MZ, NZ, Nc)
            if k not in _CACHE:
                _CACHE[k] = orig_cent(t3, c, MZ, NZ, Nc)
            return _CACHE[k]
        rd.dados_centrais = cent_cache
        _CACHE['ok'] = True
    return rd


def verifica_disco(arq, caudas, E='build/s4/rig/E.json', dtauB=1.4563e-09, tol=1e-12):
    """theta racional de cada ponto gravado, recalculado no (p, r) gravado. Criterio: theta < 1 em todos os pontos
    (rigoroso por si); a diferenca para o gravado deve estar na escala do arredondamento de ponto flutuante das
    entradas (BLAS/threads, ~1e-16), senao alguma entrada difere (alerta)."""
    rd = _rd_com_cache()
    d = json.loads(Path(arq).read_text())
    shaA = d.get('sha256_A') or '70791ce0c30f8065f93770c36ff89c1d8290978c5abf84ab0570313bc09fd002'
    erros = []
    for cauda in caudas:
        if not Path(cauda).exists():
            continue
        zr, t3, c, cent, g, MZ, NZ, Mc, Nc = rd.contexto(RAIZ/DIR_SCHUR, RAIZ/d['zrig'], cauda, RAIZ/E, dtauB,
                                                          RAIZ/d['certF'], shaA)
        zpts = {tuple(zp['p']): zp for zp in zr['pontos']}
        iguais = 0; total = 0; pior = Fr(0); difs = []
        for x in d['pontos']:
            if not x['ok']:
                continue
            total += 1
            p = complex(*x['p']); r = float(x['r'])
            zp = zpts[tuple(x['p'])]
            qz = rd.qz_ponto(cent, p, NZ)
            qf = rd.Qfar_disco(p, r, Mc, Nc)
            th, det = rd.theta_disco(p, r, zp, t3, c, cent, qz, g, exato=True, Qf_disco=qf)
            assert th < 1, f'{arq}: theta >= 1 em {x["p"]}'
            pior = max(pior, Fr(th))
            if str(th) == x['theta']:
                iguais += 1
            else:
                difs.append(abs(float(th) - x['theta_f']))
        dmax = max(difs) if difs else 0.0
        if dmax <= tol:
            return dict(ok=True, pontos=total, cauda=str(cauda).replace(str(RAIZ) + '/', ''), pior_theta=float(pior),
                        max_dtheta=dmax, exatos=iguais)
        erros.append((cauda, iguais, total, dmax, float(pior)))
    if erros:
        c_, ig, tot, dmax, pior = min(erros, key=lambda e: e[3])
        return dict(ok=False, pontos=tot, cauda=str(c_), pior_theta=pior, max_dtheta=dmax, exatos=ig,
                    alerta='theta < 1, mas a diferenca para o gravado excede o arredondamento: entrada diferente')
    return dict(ok=False, erro='nenhuma cauda encontrada', candidatos=caudas)


def candidatos_cauda(arq):
    d = json.loads(Path(arq).read_text())
    c = tuple(d['c'])
    pastas = ([RAIZ/'build/rouche_L_rig/cauda'] if 'rouche_L_rig' in str(arq)
              else [RAIZ/'build/faixaB/cauda', RAIZ/'build/rouche_L_rig/cauda'])
    cands = []
    for pasta in pastas:
        for nm in nome_c(c, pasta):
            cands += [nm, nm.replace('.json', '_lab2.json'), nm.replace('.json', '_labpc.json')]
    return cands


def resultado_disco(arq):
    """verifica um disco e grava build/reproducao/discos/<faixa>_<nome>.json (permite rodar discos em paralelo)."""
    rel = str(Path(arq).resolve().relative_to(RAIZ))
    faixa = 'A' if 'rouche_L_rig' in rel else 'B'
    out = RAIZ/'build/reproducao/discos'/f'{faixa}_{Path(arq).stem}.json'
    t0 = time.time()
    r = verifica_disco(RAIZ/rel, candidatos_cauda(RAIZ/rel))
    r['tempo_s'] = round(time.time() - t0); r['disco'] = rel
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(r, indent=1) + '\n')
    return r


def check_discos(faixa):
    if faixa == 'A':
        arqs = sorted(glob.glob(str(RAIZ/'build/rouche_L_rig/discos_v2/disco_*.json')))
        pastas = [RAIZ/'build/rouche_L_rig/cauda']
    else:
        arqs = sorted(glob.glob(str(RAIZ/'build/faixaB/discos/disco_*.json')))
        pastas = [RAIZ/'build/faixaB/cauda', RAIZ/'build/rouche_L_rig/cauda']
    res = {}
    for arq in arqs:
        pronto = RAIZ/'build/reproducao/discos'/f'{faixa}_{Path(arq).stem}.json'
        r = json.loads(pronto.read_text()) if pronto.exists() else resultado_disco(arq)    # reaproveita a rodada paralela
        res[str(Path(arq).relative_to(RAIZ))] = r
        print(f'   {Path(arq).name}: {r}', flush=True)
    ok = all(r['ok'] for r in res.values())
    return ok, res


def roda(cmd, timeout=7200):
    t0 = time.time()
    p = subprocess.run(cmd, cwd=RAIZ, capture_output=True, text=True, timeout=timeout)
    return p.returncode, p.stdout + p.stderr, round(time.time() - t0)


def check_cobertura(faixa):
    if faixa == 'A':
        discos = sorted(glob.glob('build/rouche_L_rig/discos_v2/disco_*.json', root_dir=RAIZ)); reg = '0,3,-1/4,1/4'
    else:
        discos = sorted(glob.glob('build/faixaB/discos/disco_*.json', root_dir=RAIZ)) + \
                 sorted(glob.glob('build/rouche_L_rig/discos_v2/disco_*.json', root_dir=RAIZ)); reg = '0,3,1/4,3/4'
    rc, out, t = roda([PY, '-W', 'ignore', 'scripts/rouche_L_cobertura.py', *discos, '--regiao', reg])
    ok = rc == 0 and 'CONTORNO COBERTO' in out
    return ok, dict(tempo_s=t, saida=out.strip().splitlines()[-6:])


def check_contagem(faixa):
    if faixa == 'A':
        discos = sorted(glob.glob('build/rouche_L_rig/discos_v2/disco_*.json', root_dir=RAIZ)); reg = '0,3,-0.25,0.25'
    else:
        discos = sorted(glob.glob('build/faixaB/discos/disco_*.json', root_dir=RAIZ)) + \
                 sorted(glob.glob('build/rouche_L_rig/discos_v2/disco_*.json', root_dir=RAIZ)); reg = '0,3,0.25,0.75'
    rc, out, t = roda([PY, '-W', 'ignore', 'scripts/rouche_L_contagem.py', 'rouche', '--contagem', 'build/rouche_L_rig/contagem.json',
                       '--discos', *discos, '--regiao', reg, '--diag', 'build/rouche_L_rig/diagT_lab2.npy'])
    linhas = out.strip().splitlines()
    ok = rc == 0 and any('OK' in l for l in linhas)
    zeros = [l for l in linhas if l.startswith('zeros de A + s')]
    return ok, dict(tempo_s=t, saida=linhas[-3:], zeros=zeros)


def check_nk(pasta, campos=('theta', 'r', 'erro_lambda')):
    """refaz a etapa C (Perron de 3 niveis em racionais) numa copia da pasta e compara com o C3.json gravado."""
    orig = RAIZ/pasta
    ref = json.loads((orig/'C3.json').read_text())
    with tempfile.TemporaryDirectory() as tmp:
        for f in orig.iterdir():
            if f.is_file() and f.suffix in ('.json', '.npy'):
                shutil.copy(f, tmp)
        (Path(tmp)/'C3.json').unlink(missing_ok=True)
        rc, out, t = roda([PY, '-W', 'ignore', 'scripts/nk_L3_C.py', '--dir', tmp])
        novo = json.loads((Path(tmp)/'C3.json').read_text()) if (Path(tmp)/'C3.json').exists() else {}
    iguais = {k: (str(novo.get(k)) == str(ref.get(k))) for k in campos if k in ref}
    ok = rc == 0 and 'NK VERIFICADO' in out and all(iguais.values())
    return ok, dict(tempo_s=t, iguais=iguais, saida=[l for l in out.splitlines() if 'NK' in l or 'PERRON' in l][-2:])


def check_etapa(pasta, script, saida, chaves, extra=()):
    """refaz a etapa D (testemunho de constraints) ou E (identificacao de S4) numa copia da pasta do NK e compara
    as chaves com o JSON gravado."""
    orig = RAIZ/pasta
    ref = json.loads((orig/saida).read_text())
    with tempfile.TemporaryDirectory() as tmp:
        for f in orig.iterdir():
            if f.is_file() and f.suffix in ('.json', '.npy') and f.name != saida:
                shutil.copy(f, tmp)
        rc, out, t = roda([PY, '-W', 'ignore', f'scripts/{script}', '--dir', tmp, *extra])
        novo = json.loads((Path(tmp)/saida).read_text()) if (Path(tmp)/saida).exists() else {}
    iguais = {k: (str(novo.get(k)) == str(ref.get(k))) for k in chaves}
    return rc == 0 and all(iguais.values()), dict(tempo_s=t, iguais=iguais, novo={k: novo.get(k) for k in chaves},
                                                  saida=out.strip().splitlines()[-1:] if out.strip() else [])


def check_script_json(cmd, saida_rel, chaves):
    """roda um script pequeno que grava um JSON e compara chaves com a versao do repositorio."""
    ref = json.loads((RAIZ/saida_rel).read_text())
    rc, out, t = roda(cmd)
    novo = json.loads((RAIZ/saida_rel).read_text())
    iguais = {k: novo.get(k) == ref.get(k) for k in chaves}
    if not all(iguais.values()):
        (RAIZ/saida_rel).write_text(json.dumps(ref, indent=1) + '\n')     # restaura o gravado
    return rc == 0 and all(iguais.values()), dict(tempo_s=t, iguais=iguais)


def check_s3a():
    rc, out, t = roda([PY, '-W', 'ignore', 'scripts/check_cover.py', 'build/cobertura_lab', 'build/cobertura_local'])
    linhas = out.strip().splitlines()
    ok = rc == 0 and any('COMPLETA' in l for l in linhas)
    return ok, dict(tempo_s=t, saida=linhas[-4:])


def check_janelas(faixa, tol=1e-12):
    """Janelas (blocos sem zeros em Omega, principio do maximo): refaz todas as regioes a partir das caudas gravadas
    e compara max rho_J com o gravado (mesmo valor em ponto flutuante; a fracao racional pode diferir no ultimo bit
    das entradas, como nos discos). Leve, mas ~4 min por regiao."""
    ref = json.loads((RAIZ/('build/rouche_L_rig/janelas_v2.json' if faixa == 'A' else 'build/faixaB/janelas.json')).read_text())
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp)/'janelas.json'
        rc, txt, t = roda([PY, '-W', 'ignore', 'scripts/rouche_L_janelas.py', '--faixa', faixa, '--saida', str(out)], timeout=6*3600)
        novo = json.loads(out.read_text()) if out.exists() else {}
    det = {}
    ok = rc == 0 and set(novo) == set(ref)
    for k in sorted(ref, key=int):
        a, b = novo.get(k, {}), ref[k]
        if not a:
            det[k] = 'ausente'; ok = False; continue
        ra, rb = float(Fr(a['max_rho'])), float(Fr(b['max_rho']))
        bom = ra < 1 and abs(ra - rb) <= tol*max(1, rb)
        det[k] = f'{ra:.6f} ({"=" if bom else "DIFERE de " + format(rb, ".6f")})'
        ok = ok and bom
    return ok, dict(tempo_s=t, regioes=det, saida=txt.strip().splitlines()[-2:])


def check_simbolo_L():
    """F3: recombina as faixas do simbolo de L (Arb, gravadas) e compara R_fundo e ||B|| com o certificado.
    O recalculo das faixas em Arb e nivel 2 (horas)."""
    ref = json.loads((RAIZ/'build/spectrum/symbol-L-2048x4096.json').read_text())
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp)/'L.json'
        rc, txt, t = roda([PY, '-W', 'ignore', 'scripts/certify_symbol_numerical_range.py', '--combinar', *ref['faixas'], '--out', str(out)])
        novo = json.loads(out.read_text()) if out.exists() else {}
    iguais = {k: novo.get(k) == ref.get(k) for k in ('R_fundo', 'norma_fundo')}
    ok = rc == 0 and all(iguais.values()) and ref['R_fundo'] <= 2.95 and ref['norma_fundo'] <= 3.9642
    return ok, dict(tempo_s=t, iguais=iguais, R_fundo=ref['R_fundo'], norma_fundo=ref['norma_fundo'])


def check_simbolo_K():
    """R#: certificado gravado (Arb, nivel 2 para recalcular); confere zero falhas e a cota."""
    d = json.loads((RAIZ/'build/spectrum/symbol-K-1024x2048.json').read_text())
    ok = d.get('R_fundo', 9) <= 1.4396 and d.get('norma_fundo', 9) <= 1.6497
    return ok, dict(R_fundo=d.get('R_fundo'), norma_fundo=d.get('norma_fundo'), nota='recalculo em Arb: nivel 2')


def check_rouche_K():
    """S3b: Rouche de K na circunferencia (gravado; o recalculo dos tiles e nivel 2)."""
    d = json.loads((RAIZ/'build/s3b/rouche_K_circ.json').read_text())
    c = json.loads((RAIZ/'build/s3b/rouche_K_contagem_20x80.json').read_text())
    ok = bool(d.get('certificado')) and all(t['ok'] for t in d['tiles'])
    return ok, dict(tiles=len(d['tiles']), pior_theta=d.get('pior_theta_rouche'), certificado=d.get('certificado'),
                    contagem={k: v for k, v in c.items() if not isinstance(v, (list, dict))}, nota='recalculo dos tiles: nivel 2')


def check_radius_transfer():
    ref = json.loads((RAIZ/'build/t2/radius_transfer_s3_1.json').read_text())
    rc, out, t = roda([PY, '-c', 'import sys, json; sys.path.insert(0, "scripts"); from fractions import Fraction as Q; '
                       'import radius_transfer as rt; print(json.dumps(rt.transfer_bounds(spectral_radius=Q(31, 10))))'])
    novo = json.loads(out.strip().splitlines()[-1]) if rc == 0 else {}
    chaves = ('strong_defect', 'weak_defect', 'both_contract')
    iguais = {k: novo.get(k) == ref.get(k) for k in chaves}
    return rc == 0 and all(iguais.values()) and ref.get('both_contract') is True, dict(tempo_s=t, iguais=iguais)


CHECKS = {
    'discos_A': lambda: check_discos('A'),
    'cobertura_A': lambda: check_cobertura('A'),
    'contagem_A': lambda: check_contagem('A'),
    'discos_B': lambda: check_discos('B'),
    'cobertura_B': lambda: check_cobertura('B'),
    'contagem_B': lambda: check_contagem('B'),
    'nk_s3b': lambda: check_nk('build/s3b/rig'),
    'nk_s4': lambda: check_nk('build/s4/rig'),
    'nk_raizB': lambda: check_nk('build/faixaB/nk'),
    'nk_0733': lambda: check_nk('build/nk733'),
    's3a': check_s3a,
    'F1': lambda: check_script_json([PY, 'scripts/gauge_global_centro.py'], 'build/t2/gauge_global_centro.json',
                                    ['c1_re', 'c1_im', 'cota_inferior', 'certificado', 'c1_cone_re', 'c1_cone_im', 'certificado_cone']),
    'R5': lambda: check_script_json([PY, 'scripts/radius_recovery_geral.py'], 'build/t2/radius_recovery_geral.json', ['linhas']),
    'testemunho_s3b': lambda: check_etapa('build/s3b/rig', 'nk_L3_D.py', 'D.json', ['c_ref', 'perda', 'folga', 'certificado']),
    'testemunho_raizB': lambda: check_etapa('build/faixaB/nk', 'nk_L3_D.py', 'D.json', ['c_ref', 'perda', 'folga', 'certificado']),
    'identificacao_s4': lambda: check_etapa('build/s4/rig', 'nk_L3_E.py', 'E.json', ['dist', 'r', 'identificado']),
    'janelas_A': lambda: check_janelas('A'),
    'janelas_B': lambda: check_janelas('B'),
    'L_real_toro': lambda: check_script_json([PY, 'scripts/l_real_toro.py'], 'build/t2/l_real_toro.json', ['cauda_max', 'qT', 'defeito', 'contrai']),
    'radius_transfer': check_radius_transfer,
    'simbolo_L': check_simbolo_L,
    'simbolo_K': check_simbolo_K,
    'rouche_K': check_rouche_K,
}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--so', default=None, help='lista de checks separados por virgula (padrao: todos)')
    ap.add_argument('--disco', default=None, help='verifica SO este disco (JSON) e grava o resultado (rodada paralela)')
    a = ap.parse_args()
    if a.disco:
        r = resultado_disco(a.disco)
        print(json.dumps(r)); return
    nomes = a.so.split(',') if a.so else list(CHECKS)
    rel = {}
    t0 = time.time()
    for n in nomes:
        print(f'== {n}', flush=True)
        try:
            ok, det = CHECKS[n]()
        except Exception as e:
            ok, det = False, dict(excecao=repr(e))
        rel[n] = dict(ok=ok, **det)
        print(f'   {"OK" if ok else "FALHA"}  {json.dumps(det, ensure_ascii=False)[:300]}', flush=True)
    Path(RAIZ/'build/reproducao').mkdir(parents=True, exist_ok=True)
    arq = RAIZ/'build/reproducao'/('nivel1.json' if not a.so else f'nivel1_{a.so.replace(",", "+")}.json')
    arq.write_text(json.dumps(rel, indent=1, ensure_ascii=False) + '\n')
    tudo = all(r['ok'] for r in rel.values())
    print(f'\nNIVEL 1: {"TUDO OK" if tudo else "HA FALHAS"} ({sum(r["ok"] for r in rel.values())}/{len(rel)}) [{time.time() - t0:.0f}s] -> {arq.relative_to(RAIZ)}')


if __name__ == '__main__':
    main()
