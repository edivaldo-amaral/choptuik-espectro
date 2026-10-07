# Execução do exterior (Marco 8), 17/09/2026 — sub-agente MATEMÁTICO

Meta: baixar d*(G) de 0,998 para ≤0,7, com G=107×384 inalterado.
Cortes novos: N_far=900, M'=160 (ver o Marco 8 em docs/COLUMNWISE_EXTERIOR.md).

## Comandos disparados (destacados, sobrevivem ao fim da sessão)

```bash
cd <repo>/scripts
setsid nohup ./columnwise_exterior_run.sh A 60     > <scratch>     2>&1 < /dev/null &
setsid nohup ./columnwise_exterior_run.sh sharp 100 > <scratch> 2>&1 < /dev/null &
```

Cada cadeia roda, em sequência, as três regiões do Marco 8:

| região | n | m | checkpoint |
|---|---|---|---|
| R1 | [384,606) | [107,160) | `{sis}-m107-M160-D{D}.jsonl` |
| R3 | [0,384) | [107,160) | o mesmo arquivo de R1 |
| R2 | [606,900) | [0,160) | `{sis}-M160-D{D}.jsonl` |

D=60 no setor A e D=100 no sharp. `--workers 1` por cadeia, duas cadeias
em paralelo: 2 núcleos no total, `nice -n 10`, `ulimit -v 3000000`.
Checkpoint por linha n, com `m_lo` gravado; reexecutar o mesmo comando
retoma de onde parou.

## PIDs e logs

- cadeia A: shell PID 41655, log `<scratch>`
- cadeia sharp: shell PID 41656, log `<scratch>`
- os PIDs dos processos Python mudam a cada região; acompanhe com
  `pgrep -fa "[c]olumnwise_phase2.py"`.

Se um log terminar com `FIM A` ou `FIM sharp`, aquela cadeia acabou.

## Estimativa

≈211 mil colunas no A e ≈91 mil no sharp. Nos primeiros minutos: ≈12 linhas
por 45 s na região R1. Estimativa total: **3 a 5 h de relógio** para as duas
cadeias em paralelo.

## COMANDO ÚNICO de montagem (rodar quando as duas cadeias terminarem)

```bash
cd <repo>/scripts && ../.venv/bin/python columnwise_exterior_assemble.py
```

Ele lê TODOS os checkpoints, confere em racionais que as regiões cobrem o
complemento de G sem buracos (`holes_count` deve ser 0), monta
d*=máx(regiões, d_R(900), d_F(160)) com as cotas EXATAS da Fase 1, aplica o
critério pré-registrado (`meets_0_7`), recalcula as margens do Schur
GUARDED e da ponte sharp, escreve
`build/spectrum/columnwise-exterior-phase3.json` e acrescenta o bloco
`phase3_exterior` a `build/spectrum/columnwise-exterior.json`.

Se `holes_count` não for 0, faltam linhas: reexecute a cadeia
correspondente (ela retoma pelo checkpoint) e rode a montagem de novo.

## Valores esperados (das cotas exatas já calculadas)

| sistema | d_R(900) | d_F(160) | máx. das regiões já certificadas |
|---|---:|---:|---:|
| A | 0,65116 | 0,66755 | 0,50493 |
| sharp | 0,67304 | 0,66784 | 0,46300 |

Se as regiões novas ficarem abaixo de 0,66, d* será ≈0,667 no A e ≈0,673 no
sharp, e a meta ≤0,7 é atingida. As primeiras linhas de R1 dão 0,42–0,45.

## Para a revisão do Físico

O modo "completo" do núcleo C (imagens em n dobradas, usado em R3 com
n<161) é código NOVO. Validei contra a implementação Python direta de 4
imagens em quatro colunas (A: m=108 n=1; m=107 n=3; m=110 n=25; sharp:
m=108 n=2), com o valor do C sempre maior, por ≤1,2e-4 relativo, e conferi
que numa coluna desdobrada (m=108, n=220) os dois motores coincidem em
6e-9. Vale reconferir esse modo de forma independente.

## Queda de energia, 17/09/2026 ~15:43 — e retomada

Falta de energia reiniciou a máquina (boot às 15:44). Estado dos checkpoints
verificado depois: nenhum registro corrompido (todos os JSONL fazem `parse`
limpo; o `record()` faz `fsync` por linha, então a interrupção não deixa
linha parcial).

| cadeia | região | cobertura no momento da queda | situação |
|---|---|---|---|
| sharp | R1+R3 (`sharp-m107-M160-D100`) | n=0..605, 606 linhas | **completa** |
| sharp | R2 (`sharp-M160-D100`) | n=606..899, 294 linhas | **completa** |
| A | R1+R3 (`A-m107-M160-D60`) | n=0..605, 606 linhas | **completa** |
| A | R2 (`A-M160-D60`) | n=606..793, 188 linhas | faltavam n=794..899 |

Ou seja, a cadeia `sharp` terminou inteira antes da queda; só a R2 do setor A
ficou incompleta, com 106 linhas pendentes.

Retomada (17/09 16:18), aproveitando que o núcleo do `sharp` ficou livre:

```bash
cd <repo>/scripts
ulimit -v 3000000
nice -n 10 ../.venv/bin/python columnwise_phase2.py A \
    --n-lo 606 --n-hi 900 --M 160 --m-lo 0 --depth 60 \
    --workers 3 --rows-per-call 4
```

O `main()` relê o JSONL e monta `done`, então as 188 linhas prontas são
puladas sem recomputar. Ritmo medido na retomada: ≈4 linhas/min com 3
workers. Os `max_U` da R2 do A seguem em 0,24–0,29, bem abaixo do 0,66 que a
meta exige.

Nota operacional: os logs da rodada original ficavam no scratchpad em `/tmp`,
que o reboot apagou. Os checkpoints em `build/spectrum/` são a única fonte
sobrevivente — e foram suficientes. Manter esse padrão.
