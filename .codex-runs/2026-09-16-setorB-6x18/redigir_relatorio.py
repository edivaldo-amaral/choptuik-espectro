"""Escreve o relatório somente após a conclusão efetiva das verificações."""
import json, subprocess, hashlib
from pathlib import Path
pasta=Path('.codex-runs/2026-09-16-setorB-6x18')
estado=json.loads((pasta/'estado-final.json').read_text())
assert estado['codigo_verificacao']==0
origem=json.loads((pasta/'origem-polinomio.json').read_text())
resultados=json.loads((pasta/'resultados.json').read_text())
igualdades=json.loads((pasta/'igualdades-exatas.json').read_text())
assert all(p['iguais'] and p['bareiss']==p['polinomio'] for p in igualdades['pontos'])
with (pasta/'manifesto-final.log').open('x') as f:
    manifesto=subprocess.run(['sha256sum','-c','.snapshots/2026-09-16-pre-codex/MANIFEST.sha256'],stdout=f,stderr=subprocess.STDOUT)
assert manifesto.returncode==0
texto=['# Certificados de winding da matriz finita 6×18','',
'Os dois contornos têm winding **1**, com caixas racionais aceitas pelo kernel Lean. O resultado se limita à matriz finita indicada. O kernel confere a validade combinatória das caixas e a soma dos incrementos; as inclusões analíticas são fornecidas pelo estimador exato existente.','',
'## Origem e conferência independente','',
f"Matriz: `build/spectrum/rt-A-6x18.dat`, dimensão 333, sha256 `{origem['sha256']}`.",
'Cache calculado do zero nesta execução pela função `charpoly_exact` do pipeline inalterado (877 primos). A matriz foi conferida novamente antes do uso. Nenhum cache preexistente foi lido.',
'',
'Bareiss independente: leitura racional em `Fraction`, eliminação exata dos denominadores, Bareiss inteiro acelerado por GMP, com divisibilidade conferida em **cada** divisão e reconstrução do determinante como `Fraction`. A avaliação do polinômio também usa `Fraction`. Para `A=2⁻⁷²M`, foi comparado `det(A+zI)` com `2⁻²³⁹⁷⁶ q(2⁷²z)`.',
'',
'| Ponto racional | Igualdade exata Bareiss × polinômio | Valores racionais completos |',
'|---|---|---|',
'| `1/32` | `True` | [valor Bareiss](bareiss-ponto-1.json), [ambos os valores](igualdades-exatas.json) |',
'| `1/10` | `True` | [valor Bareiss](bareiss-ponto-2.json), [ambos os valores](igualdades-exatas.json) |',
'',
'Os numeradores e denominadores são extensos e estão registrados integralmente nos arquivos vinculados; a igualdade foi testada sobre os racionais, sem tolerância.',
'',
'## Contornos','',
'Grande: `∂([1/64,1] × [1/4,3/4])`. Pequeno: `∂([1/32,1/10] × [23/50,1/2])`. Orientação anti-horária, 64 bits no arredondamento para fora e subdivisões verticais ímpares.',
'',
'| Contorno | Nós | Arestas originais | Arestas após `--coarsen` | `none` | Winding | Lean |',
'|---|---:|---:|---:|---:|---:|---|']
for r in resultados:
    texto.append(f"| {'Grande' if 'sigma' in r['nome'] else 'Pequeno'} | {r['nos']} | {r['arestas_originais']} | {r['arestas_coarsen']} | {r['none']} | {r['winding']} | saída vazia, código 0 |")
texto += ['', 'Os arquivos Lean emitidos foram verificados também quanto à ausência dos quatro termos proibidos na tarefa. Nenhum arquivo de `formal/` foi editado.', '', '## Comandos exatos', '', 'Executados em `<repo>`. Ambiente dos processos de produção:', '', '```bash', 'OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 PYTHONDONTWRITEBYTECODE=1 nohup nice -n 10 .venv/bin/python -u .codex-runs/2026-09-16-setorB-6x18/calcular_polinomio.py > .codex-runs/2026-09-16-setorB-6x18/polinomio-execucao.log 2>&1 < /dev/null &', 'TMPDIR="$PWD/.codex-runs/2026-09-16-setorB-6x18" nice -n 10 cc -O2 -shared -fPIC .codex-runs/2026-09-16-setorB-6x18/bareiss_gmp.c -lgmp -o .codex-runs/2026-09-16-setorB-6x18/bareiss_gmp.so', 'OMP_NUM_THREADS=2 PYTHONDONTWRITEBYTECODE=1 nice -n 10 .venv/bin/python .codex-runs/2026-09-16-setorB-6x18/conferir_bareiss_gmp.py --testar', 'OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 PYTHONDONTWRITEBYTECODE=1 nohup nice -n 10 .venv/bin/python -u .codex-runs/2026-09-16-setorB-6x18/conferir_bareiss_gmp.py > .codex-runs/2026-09-16-setorB-6x18/bareiss-gmp.log 2>&1 < /dev/null &', 'OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 LEAN_NUM_THREADS=2 PYTHONDONTWRITEBYTECODE=1 PYTHONINTMAXSTRDIGITS=0 nohup nice -n 10 .venv/bin/python -u .codex-runs/2026-09-16-setorB-6x18/produzir_certificados.py > .codex-runs/2026-09-16-setorB-6x18/producao-execucao.log 2>&1 < /dev/null &', '```', '', 'Os comandos filhos abaixo herdaram o ambiente e a prioridade `nice 10` do orquestrador; incluem todas as tentativas de refinamento:', '', '```bash', (pasta/'comandos-pipeline.txt').read_text().rstrip(), '```', '', 'O script de reprodução emite o Lean usando `--coarsen --standalone --quiet`, confere os incrementos esperados e regenerados, e executa, para cada arquivo temporário:', '', '```bash', '(cd formal && nice -n 10 <home>/.elan/bin/lake env lean -j2 "$lean")', '```', '', '## Saída da reprodução', '', '`bash scripts/verify_reduced_b_6x18.sh` (código 0):', '', '```text', (pasta/'verificacao-final.log').read_text().rstrip(), '```', '', '## Tempo e integridade', '', f"Tempo desde o início do cálculo modular até o término das duas verificações: {estado['segundos_totais']:.1f} s ({estado['segundos_totais']/3600:.2f} h).", '', 'A tentativa inicial de Bareiss com frações em todas as operações foi interrompida por custo; foi substituída pela eliminação exata dos denominadores descrita acima. A primeira versão inteira Python também foi interrompida por custo; a versão concluída usou GMP e passou por testes independentes por expansão em permutações (incluindo troca de pivô, singularidade e entradas racionais). O primeiro lançamento desacoplado foi encerrado pelo ambiente; os processos efetivos usaram `nohup` com sessão supervisora mantida aberta. O orquestrador foi reiniciado antes de qualquer leitura de cache para habilitar inteiros decimais longos via `PYTHONINTMAXSTRDIGITS=0`.', '', 'Os 190 arquivos do manifesto sha256 permanecem inalterados ([checagem final](manifesto-final.log)). Pipeline e arquivos formais preservados. Nada ficou incompleto.', '']
with (pasta/'RELATORIO.md').open('x') as f:
    f.write('\n'.join(texto))
print('RELATORIO.md escrito; manifesto final integralmente preservado.')
