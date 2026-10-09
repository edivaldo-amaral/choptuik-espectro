#!/usr/bin/env python3
"""Digest sha256 por passo da cadeia de T2, para o apendice C do artigo.

Para cada passo de build/reproducao/inventario.json, monta as linhas "sha256  caminho" de todos os
artefatos (o sha dos arquivos pesados vem de build/reproducao/pacote_dados.json, pelo nome do arquivo),
ordena, e calcula o sha256 do texto. Quem tiver os arquivos pode refazer o inventario e conferir os
digests; quem so tiver o repositorio confere cada arquivo contra o inventario.

Saida: build/reproducao/digests.json, paper/sec/app_certificates_table.tex e paper/sec_pt/app_certificates_table.tex
(gerados; nao editar a mao).
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

RAIZ = Path(__file__).resolve().parent.parent

# passo do inventario -> (alegacao no artigo, verificacao do nivel 1)
ALEGACOES = {
    '0': (r'data of~\cite{RT2019} (\cref{sec:setting-background})', 'download script'),
    '4': (r'\cref{lem:Lreal}', r'\texttt{L\_real\_toro}, \texttt{radius\_transfer}, \texttt{beta\_plus}'),
    '5': (r'\cref{lem:F3}, \cref{thm:counts}~(i)', r'\texttt{simbolo\_L}'),
    '6': (r'\cref{thm:counts}, A band', r'\texttt{discos\_A}, \texttt{cobertura\_A}, \texttt{contagem\_A}, \texttt{janelas\_A}, \texttt{deflacao\_q}'),
    '7, 10': (r'\cref{thm:counts}, B band; \cref{prop:nk,prop:witness} at $\sB$',
              r'\texttt{discos\_B}, \texttt{cobertura\_B}, \texttt{contagem\_B}, \texttt{janelas\_B}, \texttt{nk\_raizB}'),
    '8': (r'\cref{prop:mu}', r'\texttt{nk\_s4}, \texttt{identificacao\_s4}'),
    '9': (r'\cref{prop:nk,prop:witness} at $\sK$', r'\texttt{nk\_s3b}, \texttt{testemunho\_s3b}'),
    '11': (r'\cref{prop:Kinj}', r'\texttt{s3a}, \texttt{rouche\_K}, \texttt{simbolo\_K}'),
    '11a': (r'\cref{prop:nk} at $\sphys$', r'\texttt{nk\_0733}'),
    '11b': (r'\cref{lem:F1}', r'\texttt{F1}'),
    '11c': (r'\cref{lem:recovery}', r'\texttt{R5}'),
}

# alegacoes na versao em portugues do artigo (paper/main_pt.tex); a coluna de verificacao e a mesma
ALEGACOES_PT = {
    '0': r'dados de~\cite{RT2019} (\cref{sec:setting-background})',
    '6': r'\cref{thm:counts}, faixa A',
    '7, 10': r'\cref{thm:counts}, faixa B; \cref{prop:nk,prop:witness} em $\sB$',
    '9': r'\cref{prop:nk,prop:witness} em $\sK$',
    '11a': r'\cref{prop:nk} em $\sphys$',
}
VERIF_PT = {'download script': 'script de download'}


def tabela(linhas_tex, cabecalho):
    return '\n'.join(['% gerado por scripts/digest_certificados.py; nao editar a mao',
                      r'\begin{tabular}{@{}p{0.33\textwidth}p{0.37\textwidth}rl@{}}', r'\toprule',
                      cabecalho, r'\midrule', *linhas_tex, r'\bottomrule', r'\end{tabular}']) + '\n'


def main():
    inv = json.loads((RAIZ / 'build/reproducao/inventario.json').read_text())
    pacote = json.loads((RAIZ / 'build/reproducao/pacote_dados.json').read_text())
    sha_pesado = {Path(it['arquivo']).name: it['sha256'] for it in pacote['itens']}
    saida, linhas_tex, linhas_pt = [], [], []
    for passo in inv:
        linhas = []
        for a in passo['artefatos']:
            sha = a.get('sha256') or sha_pesado.get(Path(a['caminho']).name)
            if sha is None:
                raise SystemExit(f"sem sha256: {a['caminho']} (passo {passo['passo']})")
            linhas.append(f"{sha}  {a['caminho']}")
        linhas.sort()
        digest = hashlib.sha256(('\n'.join(linhas) + '\n').encode()).hexdigest()
        saida.append(dict(passo=passo['passo'], descricao=passo['descricao'], arquivos=len(linhas), digest=digest))
        alegacao, verif = ALEGACOES[passo['passo']]
        linhas_tex.append(f"{alegacao} & {verif} & {len(linhas)} & \\texttt{{{digest[:16]}}} \\\\")
        linhas_pt.append(f"{ALEGACOES_PT.get(passo['passo'], alegacao)} & {VERIF_PT.get(verif, verif)} & {len(linhas)} & "
                         f"\\texttt{{{digest[:16]}}} \\\\")
    (RAIZ / 'build/reproducao/digests.json').write_text(json.dumps(saida, indent=1, ensure_ascii=False) + '\n')
    (RAIZ / 'paper/sec/app_certificates_table.tex').write_text(
        tabela(linhas_tex, r'claim & level-1 check & files & digest \\'))
    (RAIZ / 'paper/sec_pt/app_certificates_table.tex').write_text(
        tabela(linhas_pt, r'alegação & conferência do nível 1 & arquivos & digest \\'))
    for s in saida:
        print(f"{s['passo']:>6}  {s['arquivos']:3d}  {s['digest'][:16]}  {s['descricao']}")


if __name__ == '__main__':
    main()
