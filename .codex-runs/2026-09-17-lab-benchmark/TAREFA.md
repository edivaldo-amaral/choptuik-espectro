# Tarefa: medir o PC do laboratório e rederivar a viabilidade da arquitetura (c)

Você é o Codex, trabalhando no projeto Choptuik. Máquina alvo: um PC de mesa do
laboratório, acessível por `ssh <usuario>@<host>` (chave já instalada, sem
senha). O notebook local tem 2 núcleos físicos; o PC do laboratório tem mais, e é
exatamente isso que precisa ser medido.

## Contexto (leia antes de agir)

- `docs/G_PARAMETRIX_FEASIBILITY.md` §7: a arquitetura (c) precisa de ~228 h-núcleo
  por região certificada, e o total é 3,4–5,8·10⁴ h-núcleo. O critério
  pré-registrado é 2 h de relógio por região em 2 núcleos; aqui isso falha por ~57×.
- `docs/SHELL_PRECONDITIONER.md`: as duas rotas locais de barateamento estão
  fechadas (pré-condicionador: teto 29–32×, melhor real 1,93×; rota k=0: piso de
  Perron 9,5–12,3 contra alvo 0,6). Ou seja, **não tente baratear o algoritmo**.
  A única variável que sobrou é hardware.
- Os fatores de ganho dependem do caminho: o float (FFT) é limitado por banda de
  memória e ganhou só 1,54× quando o clock subiu 3,3×; a aritmética exata ganhou
  2,57×. **Não extrapole por clock nem por contagem de núcleos: meça.**

## O que fazer, nesta ordem

1. **Snapshot antes de tocar em nada**, no padrão do projeto:
   `.snapshots/2026-09-17-lab/MANIFEST.sha256` (com `sha256sum`) e `sources.tar.gz`.
2. **Specs reais do PC do laboratório** (não as da documentação): `nproc`,
   núcleos FÍSICOS (`lscpu`), clock sob carga total, `energy_performance_preference`
   e governor (o notebook rodava a 1/3 do clock por causa disso — confira lá),
   RAM, disco livre, versão do Python, presença de numpy/scipy/gmpy2.
3. **Levar os fontes** para lá (rsync do diretório do projeto, sem `.git`, sem
   `build/spectrum/*.dat`) e montar um venv com `requirements.txt`.
4. **Rodar os dois benchmarks do projeto**, iguais aos daqui:
   - `scripts/clock_recalibration.py` (FFT, caminho float);
   - núcleo C exato: 12 linhas da R1 do setor A, 1 worker, com um `--depth`
     NÃO usado em produção (61 ou 62) para não tocar nos checkpoints, e
     **apague o JSONL de bench depois**:
     `python columnwise_phase2.py A --n-lo 384 --n-hi 396 --M 160 --m-lo 107 --depth 61 --workers 1 --rows-per-call 4`
     Referência local (notebook, clock já corrigido): FFT 64×200 = 0,197 s;
     12 linhas da R1 = 17,5 s normalizadas para depth 60.
   - Meça também o **escalonamento com workers** (1, 2, 4, 8, ... até nproc),
     porque é isso que decide se os núcleos extras viram tempo de relógio.
5. **Rederivar**, com os números medidos, e SÓ com eles:
   - custo por região certificada em h-núcleo e em h de relógio usando todos os
     núcleos do PC do laboratório;
   - se o critério de 2 h por região é atingido, e por qual fator falha se não for;
   - quantos dias de relógio custaria o total de 3,4–5,8·10⁴ h-núcleo lá;
   - se a máquina aguenta a memória (matriz-livre, <1 GB por processo) com todos
     os workers simultâneos.

## Regras duras

- **Só crie arquivos novos.** Não edite nenhum arquivo existente, em especial
  `README.md`, `docs/PROOF_OBLIGATIONS.md`, `docs/RETOMADA.md`,
  `formal/ChoptuikFormal.lean` e o Makefile. O ledger é integrado pelo orquestrador.
- Arquivos novos permitidos, e só estes:
  - `.codex-runs/2026-09-17-lab-benchmark/RELATORIO.md`
  - `.codex-runs/2026-09-17-lab-benchmark/lab-specs.json`
  - `.codex-runs/2026-09-17-lab-benchmark/lab-bench.json`
  - `.snapshots/2026-09-17-lab/*`
- **Não toque nos checkpoints de produção** em `build/spectrum/columnwise-phase2/`
  (`*-D60.jsonl`, `*-D100.jsonl`). Se criar JSONL de bench, apague no fim.
- Float é diagnóstico, nunca prova. Nenhuma alegação matemática nova.
- Se algo não der certo, **diga que não deu**. Relatório parcial honesto vale mais
  que número inventado. Não reescreva conclusões existentes.
- Nada de `git`, nada de `sorry`/`axiom`/`native_decide`.

## Relatório

`RELATORIO.md`, em português, dizendo: specs medidas, os dois benchmarks com os
fatores contra o notebook, a curva de escalonamento por workers, a rederivação do
custo, e um veredito explícito de uma linha — a arquitetura (c) cabe no PC do
laboratório, ou não, e por qual fator.
