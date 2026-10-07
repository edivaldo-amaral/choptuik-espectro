# Relatório: o PC do laboratório medido, e a viabilidade de (c) rederivada

**Data:** 17/09/2026 · **Máquina alvo:** `<usuario>@<host>` (ComputadorHP)
**Executado por:** Claude, direto por ssh nesta sessão. **O Codex não foi lançado** —
a chave já estava instalada e a tarefa é toda medição, então rodei eu mesmo.
`TAREFA.md` e `LANCAR.sh` ficam no diretório como o enunciado que foi seguido.

Snapshot anterior a qualquer escrita: `.snapshots/2026-09-17-lab/MANIFEST.sha256`
(4323 arquivos) e `sources.tar.gz`.

---

## 1. Specs medidas (não as da documentação)

| | notebook (referência) | PC do laboratório |
|---|---|---|
| CPU | i5-7200U (Kaby Lake, 2017) | **i7-4770 (Haswell, 2013)** |
| núcleos físicos / lógicos | 2 / 4 | **4 / 8** |
| clock máximo nominal | 3100 MHz | 3900 MHz |
| **clock medido sob carga total** | 3100 MHz (EPP corrigido) | **3692 MHz nos 8 threads** |
| driver / governor | intel_pstate ativo / performance | intel_cpufreq passivo / **ondemand** |
| `energy_performance_preference` | performance | **não existe** (pstate passivo) |
| RAM | 7,6 GiB | **15 GiB** (+38 GiB swap) |
| disco livre | — | 336 GiB em `/`, 55 GiB em `/home` |
| SO / kernel | Manjaro / 6.6.144 | Ubuntu 20.04.6 / 5.15.0-139 |

**O clock não está estrangulado lá.** Era a primeira suspeita, pelo que aconteceu
aqui: subi 8 threads em carga cheia e os 8 CPUs foram a 3692 MHz e ficaram
(94,7 % do máximo nominal), voltando a 839 MHz em repouso. O `ondemand` do
`intel_cpufreq` está fazendo o trabalho; não há nada a corrigir no EPP porque em
modo passivo o arquivo nem existe. Detalhes em `lab-specs.json`.

## 2. Uma decisão de método: igualar a pilha de software

O Python do sistema lá é **3.8.10** (Ubuntu 20.04), com numpy 1.24.4 / scipy 1.10.1 —
o teto para 3.8. Aqui o projeto roda em **3.14.6** com numpy 2.5.3 / scipy 1.18.1.
Isso não é detalhe: o caminho exato é cerca de metade Python puro, e medi 13,00 s
em 3.8 contra 10,25 s em 3.14 para o mesmo trabalho — **1,27× só de versão do
interpretador**. Medir sem igualar teria comparado software, não hardware.

Instalei `uv` e o CPython 3.14.6 **como usuário, sem root e sem PPA**
(`~/.local/`), e montei `~/Choptuik/.venv` a partir do `requirements.txt`. O
`/usr/bin/python3` do sistema não foi tocado. Único uso de root na máquina: uma
consulta ao `apt` (o `libgmp-dev` já estava instalado; o cabeçalho fica em
`/usr/include/x86_64-linux-gnu/gmp.h`, por isso o `ls /usr/include/gmp.h` engana).
Tudo é reversível com `rm -rf ~/.local/bin/uv ~/.local/share/uv ~/Choptuik`.

O núcleo C foi **recompilado lá** (`gcc -O2 -fopenmp ... -lgmp`, gcc 9.4.0): o
binário trazido daqui não roda, pede GLIBC 2.34/2.38 que o focal não tem.

**Todos os números abaixo são das duas máquinas medidas hoje, com a mesma pilha.**

## 3. Caminho float (FFT): 1,42× por núcleo

Uma aplicação matriz-livre completa de B, 18 convoluções, melhor de 5, 1 núcleo:

| grade | notebook | laboratório | ganho |
|---|---:|---:|---:|
| 24×72 | 0,0616 s | 0,0509 s | 1,21× |
| 64×200 | 0,1970 s | **0,1384 s** | **1,42×** |
| 107×384 | 0,4542 s | 0,3032 s | 1,50× |

Continua valendo o que a §7 já dizia: o caminho float é limitado por banda de
memória. O laboratório tem 19 % mais clock e uma microarquitetura **quatro anos
mais velha**, e mesmo assim ganha 1,42× — o que ele tem de fato é DDR3 dual-channel
de desktop contra a memória de um ultrabook.

## 4. Caminho exato (núcleo C): 1,41× por núcleo

12 linhas da R1 do setor A, `--depth 61`, 1 worker, 3 repetições, melhor:

| | notebook | laboratório | ganho |
|---|---:|---:|---:|
| bruto (depth 61) | 14,44 s | **10,25 s** | 1,41× |
| normalizado a depth 60 | 14,20 s | 10,09 s | 1,41× |

Nota honesta sobre a referência: a §7 registra 17,5 s para esta medida no
notebook. Hoje o mesmo comando dá 14,2 s normalizados aqui — 1,23× mais rápido que
o registrado. Não sei explicar a diferença (candidatos: contenção na medida
original, estado térmico); por isso **não usei os 17,5 s** e comparei as duas
máquinas medidas hoje, lado a lado.

## 5. A curva de escalonamento — e o achado que importa

Esta era a pergunta que decidia tudo: **os núcleos extras viram tempo de relógio?**
A resposta depende de *como* se pede paralelismo, e a diferença é grande.

**`--workers` não é paralelismo de verdade.** No motor `c`, `--workers` vira
`OMP_NUM_THREADS` de *uma* chamada do núcleo C; o laço sobre blocos de linhas, e
todo o preparo e a remontagem em `Fraction`, são seriais no processo pai. 32 linhas
no laboratório:

| `--workers` | 1 | 2 | 4 | 8 |
|---|---:|---:|---:|---:|
| tempo | 25,12 s | 16,79 s | 12,88 s | 12,91 s |
| fator | 1,00× | 1,50× | **1,95×** | 1,95× |

Amdahl sobre esses pontos dá parte serial de ~35 % e teto de ~2,85×. Com 8 threads
a máquina entrega 1,95× — metade dos 4 núcleos físicos fica ociosa.

**Dividir as linhas entre processos independentes, sim.** As linhas `n` são
independentes; um processo por faixa, `--workers 1` em cada:

| processos | 1 | 2 | 4 | 8 |
|---|---:|---:|---:|---:|
| laboratório, 32 linhas | 25,13 s | 13,61 s | **8,55 s** | 9,22 s |
| fator | 1,00× | 1,85× | **2,94×** | 2,73× |
| notebook, 16 linhas | 19,02 s | 12,05 s | 11,33 s | — |
| fator | 1,00× | 1,58× | 1,68× | — |

Com 64 linhas (partida amortizada) os 4 processos sobem para 4,13 linhas/s, contra
1,27 linhas/s de um processo: **3,25×**.

Dois fatos práticos: o hyperthreading **atrapalha** (8 processos são piores que 4 —
o caminho é ALU inteira, e as duas threads brigam pela mesma unidade); e o teto
real da máquina é **2,94 núcleos efetivos de 4 físicos**, 74 % de eficiência.

O mesmo vale no caminho float, medido com processos concorrentes (vazão em
aplicações/s): laboratório 7,23 → 13,03 → 20,08 → 21,29 para 1, 2, 4 e 8 processos
(**2,94×**); notebook 5,01 → 8,44 → 10,55 (**2,11×**). Os 8 processos do laboratório
acrescentam só 6 % sobre 4 — de novo, hyperthreading não ajuda.

**Vazão de relógio, máquina contra máquina: 2,0× no float, 2,65× no exato.**

## 6. Memória

| | por processo | 4 processos |
|---|---:|---:|
| caminho exato | 58,3 MB | 239 MB |
| caminho float, grade 107×384 | 138,6 MB | ~0,55 GB |

Muito abaixo de 1 GB por processo. Oito processos do caminho float somam ~1,1 GB
dos 15 GB da máquina: **a memória não é restrição no laboratório**, com folga.

## 7. Rederivação do custo de (c), só com os números medidos

Base da §7: **228 h-núcleo** por região certificada (caminho float, ×3 do
certificado), critério pré-registrado de **2 h de relógio** por região, total de
**3,4–5,8·10⁴ h-núcleo**.

| quantidade | notebook (§7) | notebook (medido hoje) | **laboratório** |
|---|---:|---:|---:|
| por região, em h-núcleo | 228 | 228 | **160** |
| relógio por região, todos os núcleos | 114 h (2 núcleos ideais) | 108 h (2,11×) | **54 h** (2,94×) |
| violação do critério de 2 h | ~57× | 54× | **27×** |
| total, h-núcleo da própria máquina | 3,4–5,8·10⁴ | — | 2,4–4,1·10⁴ |
| total, dias de relógio | — | 676–1153 | **338–576 dias** |

Ou seja: o PC do laboratório corta o custo de relógio **por um fator 2,0** — e
nada mais. Uma região certificada sai de 108 h para 54 h; o critério de 2 h
continua violado por **27×**; o total sai de ~2 anos para **0,9 a 1,6 ano** de
relógio ininterrupto.

Para calibrar o que faltaria: atingir as 2 h por região exigiria **80 núcleos**
como os do laboratório com escalonamento ideal, ou **109 núcleos reais** com a
eficiência de 74 % que medi. E o total em 1000 núcleos desse porte, na mesma
eficiência, custaria **32–55 h — 1,4 a 2,3 dias**.

## 8. Veredito

**A arquitetura (c) não cabe no PC do laboratório: ele dá 2,0× de relógio sobre o
notebook e o critério de 2 h por região continua violado por 27×, com o total em
338–576 dias de relógio.** A conclusão da §5 não muda — inviável sem HPC —, mas
agora com o alvo dimensionado: da ordem de 100 núcleos para uma região dentro do
critério, ~1000 núcleos para o total em poucos dias.

## 9. Ressalvas

- Os 17,5 s da §7 não se reproduziram no notebook hoje (14,2 s normalizados). Usei
  medidas pareadas de hoje; a §7 não foi reescrita.
- O ganho de 1,42× (float) e 1,41× (exato) por núcleo é medida direta; o de relógio
  (2,0× e 2,65×) combina isso com a curva de escalonamento medida. Nada foi
  extrapolado por clock nem por contagem de núcleos.
- O caminho float foi medido pelo micro-benchmark de `clock_recalibration.py`, que é
  uma aplicação de B — não uma região inteira. A rederivação herda essa proporção da
  §7, não a remede.
- `ondemand` não foi trocado por `performance` no laboratório: sob carga total o
  clock já sobe a 3692 MHz e não havia o que ganhar.
- Float é diagnóstico, nunca prova. Nenhuma alegação matemática nova aqui.
- Nada da cópia do laboratório volta por rsync. Lá, o `clock_recalibration.py`
  reescreveu `build/spectrum/clock-recalibration-2026-09-17.json` da cópia; o arquivo
  deste repositório está intacto — a auditoria do manifesto do snapshot fecha sem
  nenhuma diferença.
- Nenhum checkpoint de produção foi tocado. Os benchmarks usaram `--depth 61` e o
  JSONL `A-m107-M160-D61.jsonl` foi apagado nas duas máquinas ao fim.
