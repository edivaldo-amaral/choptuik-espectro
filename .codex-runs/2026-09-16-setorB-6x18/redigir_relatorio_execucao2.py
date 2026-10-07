"""Relatório final da retomada, somente após os dois certificados passarem."""
import json,time,subprocess
from pathlib import Path
p=Path('.codex-runs/2026-09-16-setorB-6x18')
estado=json.loads((p/'estado-final.json').read_text())
assert estado['codigo_verificacao']==0
origem=json.loads((p/'origem-polinomio-execucao2.json').read_text())
primeira=json.loads((p/'origem-polinomio.json').read_text())
conclusao=json.loads(sorted(p.glob('conclusao-polinomio-execucao2-*.json'))[-1].read_text())
resultados=json.loads((p/'resultados.json').read_text())
assert all(d['iguais'] and d['bareiss']==d['polinomio'] for d in conclusao['pontos'])
with (p/'manifesto-final-execucao2.log').open('x') as f:
    r=subprocess.run(['sha256sum','-c','.snapshots/2026-09-16-pre-codex/MANIFEST.sha256'],stdout=f,stderr=subprocess.STDOUT)
assert r.returncode==0
fim=time.time(); retomada=fim-origem['inicio']; total=fim-primeira['inicio']
texto=['# Winding da matriz finita 6×18','',
'Ambos os contornos têm winding **1**, com certificado de caixas aceito pelo kernel Lean. Este resultado se limita à matriz finita. O kernel verifica as caixas e os incrementos; as inclusões analíticas do determinante nas caixas são fornecidas pelo estimador exato existente.','',
'## Matriz, polinômio e Bareiss','',
f"Matriz `build/spectrum/rt-A-6x18.dat`, dimensão 333; SHA-256 `{origem['sha256']}`.",
f"Cache novo `polinomio-novo.txt`; SHA-256 `{conclusao['sha256_cache']}`.",
'','O cálculo foi reiniciado do zero em 17/09. A variante temporária reutiliza `_charpoly_mod`, `_coefficient_bound` e `_primes_above` de `scripts/build_contour.py`, sem alterar o arquivo. Usa os mesmos 877 primos e a mesma cota de Hadamard, com CRT incremental e checkpoints a cada 50 primos (último lote: 27). Cada checkpoint tem hash e vínculo com a matriz e o código. Não foi usado cache preexistente.','',
'Os dois determinantes Bareiss salvos e autorizados na retomada **não foram recalculados**. Foram obtidos com leitura e resultado em `Fraction`, eliminação exata de denominadores e Bareiss inteiro GMP com divisibilidade verificada em todas as operações. A comparação abaixo usa `Fraction`, sem tolerância, e confere também o SHA-256 da matriz. Para `A=2⁻⁷²M`, o valor polinomial comparado é `2⁻²³⁹⁷⁶ q(2⁷²z)`.','',
'| Ponto | Bareiss × polinômio | Valores racionais completos |',
'|---|---|---|',
'| `z=1/32` | igualdade exata: `True` | [Bareiss](bareiss-ponto-1.json), [ambos os valores](igualdades-exatas.json) |',
'| `z=1/10` | igualdade exata: `True` | [Bareiss](bareiss-ponto-2.json), [ambos os valores](igualdades-exatas.json) |',
'','Os numeradores e denominadores, com milhares de dígitos, estão integralmente nos arquivos vinculados. As comparações passaram antes de usar o cache nos contornos.','',
'## Certificados','',
'Grande: `∂([1/64,1] × [1/4,3/4])`. Pequeno: `∂([1/32,1/10] × [23/50,1/2])`. Orientação anti-horária; 64 bits de arredondamento para fora; subdivisões verticais ímpares.','',
'| Contorno | Nós | Arestas originais | Arestas após `--coarsen` | `none` | Winding | Lean |',
'|---|---:|---:|---:|---:|---:|---|']
for d in resultados:
    assert d['none']==0 and d['winding']==1
    texto.append(f"| {'Grande' if 'sigma' in d['nome'] else 'Pequeno'} | {d['nos']} | {d['arestas_originais']} | {d['arestas_coarsen']} | 0 | 1 | saída vazia; código 0 |")
texto+=['','## Comandos executados','',
'Diretório: `<repo>`. A sessão supervisora permaneceu aberta com `wait "$!"` após cada lançamento; os cálculos rodaram sob `nohup`, acompanhados pelos logs.','',
'```bash',
'OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 PYTHONDONTWRITEBYTECODE=1 PYTHONINTMAXSTRDIGITS=0 nohup nice -n 10 .venv/bin/python -u .codex-runs/2026-09-16-setorB-6x18/calcular_polinomio_checkpoint.py > .codex-runs/2026-09-16-setorB-6x18/polinomio-execucao2.log 2>&1 < /dev/null &',
'OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 LEAN_NUM_THREADS=2 PYTHONDONTWRITEBYTECODE=1 PYTHONINTMAXSTRDIGITS=0 PYTHONUNBUFFERED=1 nohup nice -n 10 .venv/bin/python -u .codex-runs/2026-09-16-setorB-6x18/continuar_execucao2.py > .codex-runs/2026-09-16-setorB-6x18/producao-execucao2.log 2>&1 < /dev/null &',
'```','',
'Comandos filhos do orquestrador, herdando o ambiente, `nice 10` e limite de memória de 1,5 GiB. Incluem todas as tentativas; arestas não classificadas causam refinamento, sem afrouxamento de critérios:','',
'```bash',(p/'comandos-pipeline.txt').read_text().rstrip(),'```','',
'O script de reprodução confere winding esperado 1, regenera os incrementos das caixas, compara-os aos arquivos entregues e emite Lean com `--coarsen --standalone --quiet`. Também rejeita os termos proibidos da tarefa. A checagem final de cada arquivo emitido é:','',
'```bash','(cd formal && nice -n 10 <home>/.elan/bin/lake env lean -j2 "$lean")','```','',
'## Saída final','',
'`bash scripts/verify_reduced_b_6x18.sh` terminou com código 0:','',
'```text',(p/'verificacao-final.log').read_text().rstrip(),'```','',
'## Tempo, memória e integridade','',
f'Retomada até o fechamento deste relatório: {retomada:.1f} s ({retomada/3600:.2f} h). Lapso desde o início da primeira execução: {total/3600:.2f} h, **incluindo a interrupção e a pausa entre sessões**. A primeira reconstrução morreu sem checkpoint após o marco 500/877 (2154 s), conforme registrado no log e informado na retomada.',
f"Pico RSS do cálculo modular retomado: {conclusao['pico_rss_kib']/1024:.1f} MiB; limite de memória virtual dos cálculos: 1536 MiB. Logs: [polinômio](polinomio-execucao2.log), [produção](producao-execucao2.log).",'',
'Os 190 arquivos do manifesto original continuam íntegros ([conferência final](manifesto-final-execucao2.log)); o pipeline e `formal/` não foram editados. Nada ficou incompleto.','']
with (p/'RELATORIO.md').open('x') as f: f.write('\n'.join(texto))
print('Relatório final escrito; manifesto preservado.')
