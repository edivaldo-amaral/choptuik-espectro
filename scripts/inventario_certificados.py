#!/usr/bin/env python3
"""Consolidacao (opcao C): inventario dos artefatos de cada alegacao da cadeia de T2.

Para cada passo de docs/T2_ENUNCIADO.md, lista os arquivos de que o certificado depende, seguindo as
dependencias gravadas nos proprios JSONs (zrig, certF, sha256_A, sha256_F dos discos), com tamanho,
sha256 e localizacao: notebook (este diretorio) e/ou laboratorio2 (indice remoto por ssh). O que nao
estiver em nenhum dos dois fica marcado 'A BUSCAR' (PC do laboratorio ou PC 3).

Saida: build/reproducao/inventario.json e docs/REPRODUCAO_MAPA.md (gerado; nao editar a mao)."""
from __future__ import annotations

import argparse
import glob
import hashlib
import json
from pathlib import Path
import subprocess
import time

RAIZ = Path(__file__).resolve().parent.parent
import os
LAB2 = os.environ.get('CHOPTUIK_HOST_REMOTO', '')    # <usuario>@host da maquina auxiliar (vazio: so local)
PESADO = 100*2**20                               # arquivos acima disto: sha so sob demanda/indice


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for bloco in iter(lambda: f.read(1 << 22), b''):
            h.update(bloco)
    return h.hexdigest()


def indice_remoto(host, sha_pesados=True):
    """{basename: [(caminho, tamanho, sha|None)]} de ~/Choptuik/build e .cache no host."""
    cmd = ("cd ~/Choptuik && find build .cache -type f -printf '%p\\t%s\\n' 2>/dev/null")
    out = subprocess.run(['ssh', '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=20', host, cmd],
                         capture_output=True, text=True, timeout=600).stdout
    idx = {}
    pesados = []
    for linha in out.splitlines():
        cam, tam = linha.rsplit('\t', 1)
        tam = int(tam)
        idx.setdefault(Path(cam).name, []).append([cam, tam, None])
        if sha_pesados and tam > PESADO and cam.endswith('.npy'):
            pesados.append(cam)
    if pesados:
        cmd = 'cd ~/Choptuik && sha256sum ' + ' '.join(pesados)
        out = subprocess.run(['ssh', '-o', 'BatchMode=yes', host, cmd], capture_output=True, text=True,
                             timeout=3600).stdout
        sha = {l.split()[1]: l.split()[0] for l in out.splitlines() if l.strip()}
        for v in idx.values():
            for e in v:
                if e[0] in sha:
                    e[2] = sha[e[0]]
    return idx


class Inventario:
    def __init__(self, remoto):
        self.remoto = remoto
        self.cache_sha = {}
        self.pesados_locais = None

    def local_pesados(self):
        """sha dos .npy pesados locais (para casar sha256_F / sha256_A)."""
        if self.pesados_locais is None:
            self.pesados_locais = {}
            for p in RAIZ.glob('build/**/*.npy'):
                if p.stat().st_size > PESADO:
                    self.pesados_locais[str(p.relative_to(RAIZ))] = sha256(p)
        return self.pesados_locais

    def item(self, cam, papel, sha_esperado=None):
        p = RAIZ/cam
        e = dict(caminho=cam, papel=papel)
        if p.exists():
            e['tamanho'] = p.stat().st_size
            if e['tamanho'] <= PESADO:
                e['sha256'] = self.cache_sha.setdefault(cam, sha256(p))
            else:
                e['sha256'] = self.local_pesados().get(cam) or sha256(p)
            e['notebook'] = True
        else:
            e['notebook'] = False
        rem = self.remoto.get(Path(cam).name, [])
        if sha_esperado:
            rem = [r for r in rem if r[2] is None or r[2].startswith(sha_esperado[:16])]
        e['laboratorio2'] = [r[0] for r in rem]
        if sha_esperado:
            e['sha256_esperado'] = sha_esperado
            if e.get('sha256') and not e['sha256'].startswith(sha_esperado[:16]):
                e['ALERTA'] = 'sha local difere do registrado'
        e['status'] = 'ok' if (e['notebook'] or e['laboratorio2']) else 'A BUSCAR'
        return e

    def por_sha(self, sha_pref, papel, nome_dica):
        """artefato pesado identificado so pelo sha (A.npy, F_*.npy)."""
        locais = [c for c, s in self.local_pesados().items() if s.startswith(sha_pref[:16])]
        remotos = [r[0] for v in self.remoto.values() for r in v if r[2] and r[2].startswith(sha_pref[:16])]
        return dict(caminho=locais[0] if locais else nome_dica, papel=papel, sha256_esperado=sha_pref,
                    notebook=bool(locais), laboratorio2=remotos,
                    status='ok' if (locais or remotos) else 'A BUSCAR')


def discos(inv, padrao):
    itens = []
    shaA = set(); shaF = set()
    for f in sorted(glob.glob(str(RAIZ/padrao))):
        rel = str(Path(f).relative_to(RAIZ))
        d = json.loads(Path(f).read_text())
        itens.append(inv.item(rel, 'disco (montagem Perron por ponto, theta racional)'))
        if d.get('zrig'):
            itens.append(inv.item(d['zrig'], 'nivel Z do disco (t_p por Loewner)'))
        if d.get('certF'):
            itens.append(inv.item(d['certF'], 'certificado do fator de posto baixo'))
        if d.get('sha256_A'):
            shaA.add(d['sha256_A'])
        if d.get('sha256_F'):
            shaF.add((d['sha256_F'], d.get('certF', '')))
    for s in sorted(shaA):
        itens.append(inv.por_sha(s, 'matriz de referencia A (Fourier-Chebyshev, Z = 32x128)', 'A.npy'))
    for s, cf in sorted(shaF):
        dica = Path(cf).name.replace('certF_', 'F_').replace('.json', '.npy')
        itens.append(inv.por_sha(s, 'fator F do disco (pesado)', dica))
    # sem duplicatas
    vistos = set(); out = []
    for e in itens:
        k = (e['caminho'], e.get('sha256_esperado'))
        if k not in vistos:
            vistos.add(k); out.append(e)
    return out


def arquivos(inv, padroes, papel):
    out = []
    for pad in padroes:
        achados = sorted(glob.glob(str(RAIZ/pad)))
        if not achados:
            out.append(inv.item(pad, papel))
        for f in achados:
            out.append(inv.item(str(Path(f).relative_to(RAIZ)), papel))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--sem-remoto', action='store_true')
    a = ap.parse_args()
    t0 = time.time()
    remoto = {} if (a.sem_remoto or not LAB2) else indice_remoto(LAB2)
    print(f'indice do laboratorio2: {sum(len(v) for v in remoto.values())} arquivos [{time.time() - t0:.0f}s]', flush=True)
    inv = Inventario(remoto)
    rt = '.cache/rt-1203.3766v1'
    passos = [
        ('0', 'existencia do fundo (RT): dados publicados', arquivos(inv, [f'{rt}/arXiv-1203.3766v1.tar.gz', f'{rt}/sourcecode/RefA.dat', f'{rt}/RefAplusB.dat'], 'dado externo (RT, arXiv 1203.3766v1)')),
        ('4', 'L-real: recuperacao de raio toro <-> RT', arquivos(inv, ['build/t2/l_real_toro.json', 'build/t2/l_real_toro.txt', 'build/t2/radius_transfer_s3_1.json', 'build/t2/beta_plus.json'], 'saida de l_real_toro.py / radius_transfer.py / independent_component_beta.py (beta+)')),
        ('5', 'F3: nada em Re s >= 2,93 (simbolo de L em Arb)', arquivos(inv, ['build/spectrum/symbol-L-2048x4096.json', 'build/spectrum/symbol-L-faixa*.json', 'build/spectrum/symbol-L-faixa*.log', 'build/spectrum/free-part-toro.json'], 'certificado do simbolo (Arb) / parte livre')),
        ('6', 'faixa A: N_A = 3 (Rouche de L)', discos(inv, 'build/rouche_L_rig/discos_v2/disco_*.json')
              + arquivos(inv, ['build/rouche_L_rig/cauda/T3_*.json', 'build/rouche_L_rig/janelas_v2.json', 'build/rouche_L_rig/janelas_v2.txt',
                                'build/rouche_L_rig/cobertura_v2.txt', 'build/rouche_L_rig/contagem.json', 'build/rouche_L_rig/contagem_rouche_v2.txt',
                                'build/rouche_L_rig/diagT_lab2.npy', 'build/rouche_L_rig/z/lab2/schur.json'], 'cauda / janelas / cobertura / contagem / deflacao')),
        ('7, 10', 'faixa B: N_B = 1 e a raiz B simples (NK)', discos(inv, 'build/faixaB/discos/disco_*.json')
              + arquivos(inv, ['build/faixaB/cauda/T3_*.json', 'build/faixaB/janelas.json', 'build/faixaB/janelas.txt', 'build/faixaB/cobertura.txt',
                                'build/faixaB/contagem_rouche.txt', 'build/faixaB/nk/*', 'build/faixaB/logs/*', 'build/rouche_L_rig/z/lab2/schur.json'], 'cauda / janelas / cobertura / contagem / NK / logs / deflacao')),
        ('8', 'S4: mu* simples, autoespaco span(g*), Dg* = 0', arquivos(inv, ['build/s4/rig/*', 'build/s4/corr_gauge.json', 'build/s4/nk3_mu.json', 'build/s4/*.txt'], 'NK bordejado de L em mu, identificacao (E), lema de Corr')),
        ('9', 's_K simples, viola as constraints (S3b)', arquivos(inv, ['build/s3b/rig/*', 'build/s3b/nk3_32x128.json'], 'NK de L, testemunho D')),
        ('11', 'K injetivo na faixa A fora de s_K: S3a (tiles), Rouche de K, R# (simbolo de K)', arquivos(inv, ['build/cobertura_lab/grande-*.json', 'build/cobertura_local/grande-*.json', 'build/s3b/rouche_K_circ.json', 'build/s3b/rouche_K_contagem_20x80.json', 'build/spectrum/symbol-K-1024x2048.json', 'build/spectrum/symbol-K-1024x2048.log', 'build/spectrum/free-part-toro.json'], 'tiles de Perron / Rouche de K no disco / simbolo de K / parte livre')),
        ('11a', '0,7332 simples (NK de L)', arquivos(inv, ['build/nk733/*'], 'NK bordejado de L em 0,7332')),
        ('11b', 'Lema G: fato F1/F1\' (omega4 nao nula no centro e no cone)', arquivos(inv, ['build/t2/gauge_global_centro.json'], 'racionais exatos')),
        ('11c', 'Lema R: R5 (recuperacao de raio de qualquer peso)', arquivos(inv, ['build/t2/radius_recovery_geral.json'], 'racionais exatos')),
    ]
    saida = []
    for passo, desc, itens in passos:
        n_buscar = sum(e['status'] != 'ok' for e in itens)
        print(f'passo {passo:6s} {len(itens):4d} artefatos, {n_buscar} a buscar  -- {desc}', flush=True)
        saida.append(dict(passo=passo, descricao=desc, artefatos=itens))
    Path(RAIZ/'build/reproducao').mkdir(parents=True, exist_ok=True)
    (RAIZ/'build/reproducao/inventario.json').write_text(json.dumps(saida, indent=1) + '\n')
    # resumo legivel
    L = ['# Mapa alegação → artefato (gerado por `scripts/inventario_certificados.py`; não editar)', '',
         f'Gerado em {time.strftime("%d/%m/%Y %H:%M")}. "NB": no notebook (repositório). "lab2": no laboratorio2. '
         '"A BUSCAR": em nenhum dos dois (PC do laboratório ou PC 3).', '']
    for s in saida:
        itens = s['artefatos']
        tot = sum(e.get('tamanho', 0) for e in itens)
        L += [f'## Passo {s["passo"]}: {s["descricao"]}', '',
              f'{len(itens)} artefatos, {tot/2**20:.1f} MB no notebook; '
              f'{sum(e["status"] != "ok" for e in itens)} a buscar.', '',
              '| artefato | papel | NB | lab2 | sha256 |', '|---|---|---|---|---|']
        for e in itens:
            sh = (e.get('sha256') or e.get('sha256_esperado') or '')[:12]
            alerta = (' **' + e['ALERTA'] + '**') if e.get('ALERTA') else ''
            L.append(f'| `{e["caminho"]}` | {e["papel"]} | {"sim" if e["notebook"] else "—"} | '
                     f'{"sim" if e["laboratorio2"] else "—"} | {sh}{alerta}{" **A BUSCAR**" if e["status"] != "ok" else ""} |')
        L.append('')
    (RAIZ/'docs/REPRODUCAO_MAPA.md').write_text('\n'.join(L) + '\n')
    print(f'gravados build/reproducao/inventario.json e docs/REPRODUCAO_MAPA.md [{time.time() - t0:.0f}s]')


if __name__ == '__main__':
    main()
