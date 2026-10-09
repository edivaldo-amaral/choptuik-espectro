#!/usr/bin/env python3
"""Exporta uma arvore LIMPA do repositorio para publicacao (docs/PUBLICACAO.md). Nao publica nada: so escreve num
diretorio local, que o autor revisa antes de qualquer envio.

- Copia apenas arquivos rastreados pelo git (git ls-files).
- Exclui os logs brutos das revisoes (.codex-runs/**/*.log, codex*.log) e os tarballs de snapshot
  (.snapshots/**/sources.tar.gz); eles vao para o arquivo de dados.
- Redige nos arquivos de texto: IPs das maquinas, usuarios de ssh, caminhos absolutos locais e "PC 3".
  Os rotulos de diretorio dentro dos certificados (z2/joao/...) sao mantidos, porque renomea-los mudaria o sha256.
- Confere ao final que nenhum padrao sensivel sobrou.

Uso: .venv/bin/python scripts/exporta_publico.py --destino <dir>"""
from __future__ import annotations

import argparse
from pathlib import Path
import re
import shutil
import subprocess

RAIZ = Path(__file__).resolve().parent.parent
# literais montados por concatenacao: o proprio script nao pode casar com os padroes que procura
HOME_ = '/home/'
TMPC = '/tmp/' + 'claude-'
USR = 'edivaldo'

EXCLUI = [
    re.compile(r'^\.codex-runs/.*\.log$'),
    re.compile(r'^\.snapshots/.*sources\.tar\.gz$'),
    re.compile(r'^prompt\.MD$'),                 # prompt interno do loop de setembro
    re.compile(r'^probe_result\.json$'),          # benchmark da maquina local
    re.compile(r'^\.codex-runs/.*\.pdf$'),          # PDFs de terceiros e rascunhos copiados para as revisoes
    re.compile(r'^paper/ESTRUTURA\.md$'),          # nota interna de planejamento do artigo
    re.compile(r'^submissao/'),                   # carta ao editor e notas da submissao (revisores sugeridos)
]
TROCAS = [
    (re.compile(r'\b100\.(?:\d{1,3})\.(?:\d{1,3})\.(?:\d{1,3})\b'), '<host>'),
    (re.compile(r'\b192\.168\.\d{1,3}\.\d{1,3}\b'), '<host-local>'),
    (re.compile(r'\b(?:usuario|joao|edivaldo)@'), '<usuario>@'),
    (re.compile(re.escape(HOME_ + USR) + r'/Choptuik'), '<repo>'),
    (re.compile(re.escape(TMPC) + r'\d+/[^\s\'"`)]*'), '<scratch>'),
    (re.compile(re.escape(HOME_) + r'(?:' + USR + r'|usuario|joao)\b'), '<home>'),
    (re.compile(r'PC do Jo(?:ã|a)o', re.I), 'PC 3'),
]
SENSIVEL = [
    re.compile(r'\b100\.(?:108|119|122)\.\d{1,3}\.\d{1,3}\b'),
    re.compile(r'192\.168\.\d{1,3}\.\d{1,3}'),
    re.compile(re.escape(HOME_ + USR)),
    re.compile(re.escape(TMPC)),
    re.compile(r'\b(?:usuario|joao)@'),
    re.compile(r'PC do Jo(?:ã|a)o', re.I),
]


def texto(p: Path) -> bool:
    try:
        b = p.read_bytes()[:4096]
    except OSError:
        return False
    return b'\x00' not in b


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--destino', type=Path, required=True)
    a = ap.parse_args()
    dest = a.destino.resolve()
    if dest == RAIZ or RAIZ in dest.parents:
        raise SystemExit('destino deve ficar fora do repositorio')
    # a exportacao copia a copia de trabalho: mudancas fora de commit iriam a publico (ja aconteceu com um
    # certificado reformatado pelo verificador, com sha256 diferente do inventario)
    sujos = subprocess.run(['git', 'diff', '--name-only', 'HEAD'], cwd=RAIZ, capture_output=True, text=True).stdout.split()
    if sujos:
        raise SystemExit('arquivos rastreados com mudancas fora de commit: ' + ', '.join(sujos[:10]))
    if dest.exists():
        shutil.rmtree(dest)
    arquivos = subprocess.run(['git', 'ls-files'], cwd=RAIZ, capture_output=True, text=True).stdout.splitlines()
    n = nex = nred = 0
    for rel in arquivos:
        if any(r.search(rel) for r in EXCLUI):
            nex += 1
            continue
        src = RAIZ/rel
        if not src.is_file():
            continue
        dst = dest/rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        if texto(src):
            s = src.read_text(errors='surrogateescape')
            s2 = s
            for padrao, troca in TROCAS:
                s2 = padrao.sub(troca, s2)
            if s2 != s:
                nred += 1
            dst.write_text(s2, errors='surrogateescape')
        else:
            shutil.copy2(src, dst)
        n += 1
    restos = []
    for p in dest.rglob('*'):
        if p.is_file() and texto(p):
            s = p.read_text(errors='surrogateescape')
            for r in SENSIVEL:
                m = r.search(s)
                if m:
                    restos.append((str(p.relative_to(dest)), m.group(0)))
    tam = sum(p.stat().st_size for p in dest.rglob('*') if p.is_file())
    print(f'exportados {n} arquivos ({tam/2**20:.1f} MB), excluidos {nex}, redigidos {nred}')
    print('padroes sensiveis restantes:', restos[:20] if restos else 'nenhum')
    if restos:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
