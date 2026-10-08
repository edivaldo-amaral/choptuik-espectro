#!/usr/bin/env python3
"""R7: refaz os digests do apendice C (copia de scripts/digest_certificados.py, escrevendo SO nesta pasta) e,
independentemente, recalcula o sha256 de cada arquivo do inventario que existe no disco, comparando com o
inventario (e com pacote_dados.json para os pesados). Compara os digests com paper/sec/app_certificates_table.tex
e com build/reproducao/digests.json."""
import hashlib, json, re
from pathlib import Path
RAIZ = Path(__file__).resolve().parents[2]
AQUI = Path(__file__).resolve().parent
inv = json.loads((RAIZ/'build/reproducao/inventario.json').read_text())
pacote = json.loads((RAIZ/'build/reproducao/pacote_dados.json').read_text())
sha_pesado = {Path(it['arquivo']).name: it['sha256'] for it in pacote['itens']}
tab = (RAIZ/'paper/sec/app_certificates_table.tex').read_text()
dig_tab = re.findall(r'& (\d+) & \\texttt\{([0-9a-f]{16})\}', tab)
dig_json = json.loads((RAIZ/'build/reproducao/digests.json').read_text())
out = []
def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()
divergentes = 0; ausentes = []
for k, passo in enumerate(inv):
    linhas = []
    for a in passo['artefatos']:
        s = a.get('sha256') or sha_pesado.get(Path(a['caminho']).name)
        linhas.append(f"{s}  {a['caminho']}")
        p = RAIZ/a['caminho']
        if p.exists():
            real = sha(p)
            esperado = a.get('sha256') or a.get('sha256_esperado') or sha_pesado.get(p.name)
            if real != esperado:
                divergentes += 1; out.append(f"  DIVERGE: {a['caminho']}: disco {real[:16]} inventario {str(esperado)[:16]}")
        else:
            ausentes.append(a['caminho'])
    linhas.sort()
    d = hashlib.sha256(('\n'.join(linhas) + '\n').encode()).hexdigest()
    nt, dt = dig_tab[k]
    dj = dig_json[k]
    ok = (d[:16] == dt and int(nt) == len(linhas) and dj['digest'] == d and dj['arquivos'] == len(linhas))
    out.append(f"passo {passo['passo']:>6}: {len(linhas):3d} arquivos, digest {d[:16]} | tabela {nt} {dt} | digests.json {dj['digest'][:16]} -> {'OK' if ok else 'DIFERE'}")
out.append(f"arquivos do inventario com sha no disco diferente do inventario: {divergentes}")
out.append(f"arquivos do inventario ausentes neste disco (pesados, no Zenodo/lab2): {len(ausentes)}: {sorted(set(ausentes))}")
print('\n'.join(out))
