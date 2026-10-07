# Viabilidade de uma parametriz certificada na caixa G

Autor: sub-agente MATEMÁTICO (Claude), rodada 4, 17/09/2026. Estudo de
viabilidade, não implementação.

> **Tempos recalibrados em 17/09/2026 (tarde), §7.** As estimativas de custo
> das §§3–5 foram medidas com a máquina a 1/3 do clock. Os valores válidos
> são os da §7; a decisão (INVIÁVEL SEM HPC) não mudou.

## 0. Pré-registro: regra de decisão

Escrita ANTES das curvas de trade-off, das estimativas de custo e de
qualquer protótipo (itens 2–4). Não será alterada depois; mudanças vão em
seções datadas no final.

**Pergunta.** Existe uma arquitetura que produza, nesta máquina (4 núcleos,
7,6 GB), um certificado de invertibilidade de H(s)=I+(B_true+s)Q_ref (setor
A, norma η) e de H♯(s) (sharp), uniforme em discos que cubram:
- o contorno Γ_A: 1/8≤Re s≤1, |Im s|≤1/4;
- o retângulo B: 1/8≤Re s≤1, 1/4≤Im s≤3/4;
- o disco sharp |s-0,733|≤1/32?

O certificado inclui o acoplamento com o exterior infinito e a bola RT.

**A rota é "VIÁVEL" se existir UMA arquitetura com, simultaneamente:**

1. memória de pico ≤5 GB por job;
2. tempo de relógio por centro ≤2 h em 2 núcleos;
3. número de centros ≤200 para Γ_A + retângulo B + disco sharp;
4. **defeito projetado ≤0,9**, uniforme em cada disco, incluindo o
   acoplamento com o exterior.

Senão, a rota é "**INVIÁVEL SEM HPC**", e registra-se o que HPC exigiria:
núcleos-hora, memória e número de centros.

**Como os números são projetados** (fixado antes de vê-los):

- **Defeito.** Uso normas de coluna em FLOAT (diagnóstico), multiplicadas
  pelo fator de segurança **1,5**. O fator vem da rodada 2: a cota
  certificada por coluna ficou ≤1,42× o float em 47 colunas. As cotas
  analíticas da Fase 1 entram sem fator, porque já são exatas. Se o
  defeito vier de uma série de Neumann amortecida, projeta-se a potência
  medida numericamente num centro representativo (o pior entre s=1/8+i/4,
  s=1 e s=0,733 sharp), também ×1,5.
- **Memória.** Contam os bytes das matrizes ou vetores guardados, com
  inteiros diádicos de 64 bits por entrada (mais raio, se houver). Não
  conto otimizações hipotéticas.
- **Tempo.** Medido em micro-benchmark na própria máquina e extrapolado
  linearmente no número de colunas e produtos. Na aritmética certificada
  aplica-se um fator ×3 sobre o float (inteiros e raio), a menos que o
  benchmark do núcleo C exato da rodada 2 dê número melhor.
- **Centros.** O raio do disco segue do pencil afim,
  r_j=(1-a_j)/(2‖V_j‖‖Q‖), com a_j o defeito no centro e ‖V_j‖ a norma
  projetada da parametriz. Conto os centros para cobrir o perímetro de Γ_A
  e do retângulo B e o disco sharp.
- **Tolerância de fronteira.** Se um número cair a menos de 10% de um
  limiar, a decisão é "**VIÁVEL COM RISCO**". Isso é registrado e tratado
  como viável, com o risco explícito.

Por que esses limiares:
- 5 GB deixa folga para o SO numa máquina de 7,6 GB dedicada. Com os jobs
  concorrentes de hoje (~3 GB), a execução seria sequencial.
- 200 centros × 2 h = 400 h em 2 núcleos, ou ≈8 dias com dois jobs em
  paralelo nos 4 núcleos. É o máximo razoável para uma rodada longa.
- 0,9 deixa margem de 0,1 para absorver os erros de implementação e o
  arredondamento da versão certificada.

## Estado

(seções datadas abaixo)

### 1. Precedente RT (lido em `choptuik.tex` §4–§5 e em `sourcecode/choptuik.c`; números conferidos nos arquivos)

**Caixas** (`SECTOR_init`):

| caixa | tamanho | papel |
|---|---|---|
| inner | 25×75 | inversa densa exata |
| outer K_S | 250×750 | onde a inversa é certificada |
| EXTERIOR_ESTIMATE | 1600×4800 | colunas estimadas uma a uma |

- DOFs reais em K_S: **654375** (setor selecionado, INV) e **280125**
  (sharp, SHARP_INV). Soma: 934500 linhas. Conferi a contagem: 93375 +
  2×186750 + 187500.
- DOFs de inner 25×75: ≈6475 (+μ).

**Construção de V** (`AUXINV_INNER_calculate_save`, `INV`, `INV_one_iteration`,
`J_UU_approx_inv_NoENL_ROUND_alloc_free`):
- Matriz DENSA de J_UU no inner, invertida com arredondamento diádico
  (`Matrix_Inverse_ROUND`, TwoExp=-75).
- Fora do inner, aproximação de 1ª ordem: J^-1 - J^-1 A J^-1.
- V é montada COLUNA A COLUNA. Para cada vetor unitário e∈K_S: refinamento
  iterativo r←r+M(Je-J_UU r), até 8 iterações, num setor LOCAL em torno de
  e (margem ±25 em m e ±75 em n). O setor cresce em (20,40) se não
  converge.
- Alvo: erro relativo ≤2^-16. Ao final, cota exata de ‖e-U f‖/‖e‖ e de
  ‖f‖/‖e‖.
- Aritmética: DyadicQ (mantissa GMP × 2^e), arredondada em TwoExp=-75 (e
  -40 a cada expansão), com cotas finais EXATAS. Executado em cluster HPC,
  paralelo por setores.

**Margens:**

| grandeza | valor |
|---|---|
| máx. ‖e-Uf‖/‖e‖ | 19596·2^-32 = 4,56e-6 (alvo 2^-16) |
| máx. ‖f‖/‖e‖ | 165,5, contra 𝔄(e)=256 |
| ‖1-UV‖ | ≤2^-16 |
| ‖V𝔄^-1‖ | ≤1 |
| 𝔄 | 256/16/4/2/1,5625/1 em cascas 50×150…250×750 |
| L5=‖U^-1𝔄^-1‖ | ≤1,004 |
| L6=‖U^-1‖ | ≤258 |
| L4=‖𝔄O^-1(A-KAK)‖ | ≤0,625 |

L4 vem de 26 868 000 colunas calculadas em 1600×4800 (máximo típico
550·2^-10≈0,537 na fronteira) mais a cota analítica além (M3≤0,027,
K3≤22,5). A contração final é S23=0,64, isto é, **margem 0,36**.

**Lição estrutural.**
1. RT NÃO fazem Schur com d≈1. Fazem Neumann única no ℓ¹ ponderado por 𝔄:
   ‖U^-1𝔄^-1‖·‖𝔄(H-U)‖ ≈1,004·0,625. A casca 𝔄=1,5625 absorve o fato de
   U-I ter colunas ~0,5 na fronteira de K_S.
2. A cauda analítica só entra muito longe (1600×4800), onde é desprezível.
   Nas nossas contas da rodada 2 é o contrário: a cauda analítica (0,998)
   domina.
3. A escala de RT (6,5·10^5 colunas locais com refinamento + 2,7·10^7
   colunas exteriores) exigiu HPC, e isso num ÚNICO ponto s. Não havia
   contorno.

### 2. Curva de trade-off (setor A)

Script `scripts/g_parametrix_tradeoff.py`, relatório
`build/spectrum/g-parametrix-tradeoff.json`. Conteúdo:
- d_float: DIAGNÓSTICO, máximo amostrado das normas reais de coluna de
  (B_ref+5/4)Qe fora de G';
- d_R e d_F: cotas EXATAS da Fase 1 (d_F escala exatamente como 1/M).

Normas reais máximas (float, s=5/4):

| direção | corte | máx. | onde |
|---|---|---:|---|
| radial | n=128 | 1,480 | K, m=6 |
| radial | n=192 | 0,992 | |
| radial | n=256 | 0,746 | |
| radial | n=384 | 0,498 | K, m=18 |
| Fourier | m=24 | 1,857 | d=2, n=65 |
| Fourier | m=48 | 0,951 | d=2, n=127 |
| Fourier | m=64 | 0,664 | |
| Fourier | m=107 | 0,294 | |

Na direção de Fourier o pior caso está em n moderado. A amostra da rodada
1 (n≤17) subestimava.

| G'=(M',N') | DOFs A | DOFs sharp | d_float×1,5 | d Fase 1 exata (máx(d_R,d_F)) |
|---|---:|---:|---:|---:|
| 24×128 | 10432 | 4416 | 2,79 | 4,55 |
| 48×192 | 31776 | 13536 | 1,49 | 3,04 |
| 48×384 | 63552 | 27072 | 1,43 | 2,23 |
| 64×256 | 56704 | 24192 | 1,12 | 2,28 |
| 64×384 | 85056 | 36288 | 1,00 | 1,67 |
| 107×256 | 95616 | 41088 | 1,12 | 2,28 |
| **107×384** | **143424** | **61632** | **0,75** | 1,52 (d_R); com a faixa exata da rodada 2: 0,998 |

**Melhor compromisso:** 107×384. É a ÚNICA caixa da grade com d_float×1,5<0,9.
Isso já era esperado das normas reais ~190/n e ~30/m, e o exponho sem
ajustar. A caixa 64×384 fica exatamente no limite (0,996).

Consequência para qualquer arquitetura no estilo de RT, com Neumann única
no ℓ¹ ponderado por 𝔄: o bloco a inverter tem ≥1,4·10^5 DOFs no A e
≥6·10^4 no sharp. Caixas menores só servem dentro de um esquema multinível
em que as cascas intermediárias também sejam invertidas numericamente.

### 3. Arquiteturas e custo (projeções com as regras da §0; tudo em float é DIAGNÓSTICO)

Medições em que as projeções se apoiam:
- `scripts/g_parametrix_prototype.py` → `build/spectrum/g-parametrix-prototype.json`;
- `scripts/g_parametrix_contour_norms.py` → `build/spectrum/g-parametrix-contour-norms.json`.

Ambos usam a matriz exata 12×36 (1422 DOFs) em float, na norma η × pesos RT.

**Micro-benchmark** (1 núcleo, float): uma aplicação matriz-livre completa
de B, com 18 convoluções FFT com RefA.

| grade | tempo |
|---|---:|
| 24×72 | 90 ms |
| 64×200 | 437 ms |
| 107×384 | 979 ms |

**Normas da resolvente pesada no protótipo 12×36** (float). Ao longo de
Γ_A e do retângulo B, ‖QH(s)^-1‖ vai de **84 a 6274**:

| ponto | ‖QH(s)^-1‖ | observação |
|---|---:|---|
| s=1/8 | 6274 | raiz de L_A em 0,168, a 0,043 |
| 1/8+i/4 | 879 | |
| 1+i/4 | 129 | |
| 1/8+i/2 | 1031 | candidato B em 0,041+0,5i |
| 1+i/2 | 84 | |

‖H(s)^-1‖ vai de 121 a 7722. Para comparação, a coroa sozinha tinha κ≤31
e RT tinham ‖U^-1‖≤258 num único ponto.

**(a) Inversa densa por blocos no estilo RT** (na caixa inteira).
- Memória: 143424² × 16 B ≈ **329 GB** (A) e 61632² × 16 B ≈ **61 GB**
  (sharp). Viola 5 GB por ≥12×.
- Montagem coluna a coluna como RT, sem guardar a matriz: o protótipo
  mostra que a aproximação de 1ª ordem na casca (P1, o esquema de RT) NÃO
  contrai nosso H(s0).

  | centro | defeito ℓ¹ de P1 | após 8 refinamentos |
  |---|---:|---:|
  | s=1 | 52 | 48 |
  | 1/8+i/4 | 312 | 3,8 |

  O motivo: nossas colunas de (B+s)Q na casca próxima têm norma até ≈58,
  contra as de RT, que ficavam abaixo de 1 fora do inner.
- Centros: parametriz por centro, com raio pela regra pré-registrada
  r=(1-a)/(2‖V‖‖Q‖), ‖Q‖=5,94, ‖V‖≈‖H^-1‖∈[121,7722]. Dá r≈10^-5–10^-3,
  isto é, **~5·10^4 centros**. Mesmo com a cota mais fina ‖QV‖, ~7·10^3.
  Viola 200 por 35–250×.

**(b) Núcleo denso + Neumann amortecida matriz-livre** (sugestão do
orquestrador), uma parametriz por centro. O protótipo passa no defeito:
Schur exato no núcleo 6×18 + Neumann amortecida de grau k na casca.

| centro | ω, k | defeito |
|---|---|---:|
| s=1 | ω=1/2, k=16 | 0,033 |
| s=1 | ω=1/2, k=32 | 9·10^-6 |
| 1/8+i/4 | ω=1, k=16 | 0,027 |
| 1/8+i/4 | ω=1, k=32 | 1,6·10^-5 |

Mas ‖P‖≈‖H^-1‖ herda os mesmos 10^2–10^3. Os centros são os mesmos de (a),
**~5·10^4**, e cada um exige a norma ℓ¹ de I-HP, ou seja, TODAS as ~1,4·10^5
colunas com k≈32 aplicações matriz-livres. Tempo por centro ≈
1,4·10^5 × 32 × 0,98 s ≈ 1,2·10^3 h em float, 3,7·10^3 h certificado.
Viola 2 h por ~10^3×.

**(c) Arquitetura proposta: Schur num núcleo + casca uniforme em s +
exterior no estilo RT.**
1. Núcleo C, por exemplo 8×24 (612 DOFs no A): o complemento de Schur
   S(s)=H_CC-H_CO V_O(s) H_OC é uma função analítica FINITA de s. O winding
   de det S(s) sai da maquinaria de caixas já existente, com erro de
   enclausuramento de V_O e do exterior. Não precisa de ‖H^-1‖, que é o que
   explode os centros.
2. Casca O=G'-C: V_O(s) = Neumann amortecida de grau 32, certificada
   UNIFORMEMENTE em regiões. O bloco de casca é bem condicionado
   (protótipo, float):

   | núcleo | ‖H_OO(s)^-1‖ | defeito k=32 | raio projetado (regra pré-registrada) |
   |---|---:|---:|---:|
   | 6×18 | 7–29 | ≤0,22 | 0,005–0,024 |
   | 8×24 | 5–13 | ≤8·10^-4 | 0,013–0,034 |

   Regiões: com núcleo 8×24, ≈50–90 para meio Γ_A e ≈60–100 para meio
   retângulo B (simetria s↦s̄+i), mais 1–2 no disco sharp: **≈110–190 ≤200**
   (menos de 10% abaixo do limiar → risco). Com núcleo 6×18 seriam ≈200–300.
3. Exterior além de G': colunas exatas no estilo RT, no núcleo C da rodada
   2 (≈0,04 s por coluna, exato), mais cauda analítica. É barato: ~3 h
   no total para empurrar as caudas a ≤0,7.

Custo da casca por REGIÃO, restrito à faixa em que a Neumann simples não
basta (m<64, n<200; ≈4,4·10^4 colunas no A): 4,4·10^4 × 32 × 0,44 s ≈
**170 h em float → 510 h certificado (×3), ≈255 h de relógio em 2 núcleos**.
Viola 2 h por ~130×. No sharp, ≈1/3 disso. Memória ✓ (matriz-livre, <1 GB).
Defeito projetado ✓ na casca. O acoplamento com o exterior (L4·L5 no
estilo RT) não foi medido; a projeção pela curva da §2 dá ≈0,75 em 107×384.

**Outras ideias consideradas e descartadas nesta rodada.**
- Pré-condicionador pela parte m-diagonal de B (blocos por m): a parte
  m≠0 de RefA sozinha já dá β_η=9,71 de 10,38 (majorante exato). Absorver
  m=0 quase não reduz as normas de coluna.
- Aplicar muitas colunas de uma vez: a matriz H_OO densa na faixa de
  4,4·10^4 DOFs ocuparia ≈31 GB, e esparsa ela não é, pelo suporte 81×201
  do fundo.

### 4. Protótipo: resultado

Feito em float, sobre a caixa reduzida 12×36 (núcleo 6×18, casca 1089
DOFs), nos centros s=1, s=1/8+i/4 e s=0,733. Tempo total <5 min de CPU.

**O defeito <1 É atingível:**

| arquitetura | centro | defeito |
|---|---|---:|
| (b)/(c) Schur + Neumann amortecida k=16 | s=1 | 0,033 |
| (b)/(c) Schur + Neumann amortecida k=16 | 1/8+i/4 | 0,027 |
| 1ª ordem estilo RT (P1) | s=1 | 52 |
| 1ª ordem estilo RT (P1) | 1/8+i/4 | 312 |

O centro s=0,733, no setor A, fica SOBRE a raiz física de L_A. Ali
‖H^-1‖≈8,9·10^5 e o defeito fica pequeno mas com P gigante: é o
comportamento esperado de um centro mal escolhido (o 0,733 pertence ao
disco SHARP, não ao contorno A). O gargalo não é o defeito: é
**tempo × número de centros**.

### 5. DECISÃO pela regra pré-registrada: **INVIÁVEL SEM HPC**

| arquitetura | memória ≤5 GB | tempo por centro ≤2 h | centros ≤200 | defeito ≤0,9 |
|---|---|---|---|---|
| (a) densa estilo RT | ✗ (329 GB) | – | ✗ (~5·10^4) | ✗ (1ª ordem não contrai) |
| (b) núcleo + Neumann por centro | ✓ | ✗ (~3,7·10^3 h) | ✗ (~5·10^4) | ✓ |
| (c) Schur + casca uniforme + exterior | ✓ | ✗ (~255 h por região) | ✓ com risco (110–190) | ✓ (casca); acoplamento ≈0,75 projetado |

Nenhuma arquitetura cumpre as quatro condições. A melhor, (c), falha só
no tempo, por ~130×.

**O que HPC exigiria para (c):**
- (110–190 regiões) × (≈510 h-núcleo no A + ≈170 h-núcleo no sharp, por
  região certificada) + exterior (~10^3 h-núcleo) ≈ **1–1,5·10^5
  h-núcleo**;
- memória por processo <2 GB, trivialmente paralelo por coluna e por
  região: ≈1000 núcleos por ~5 dias.

Para comparação, RT fizeram 6,5·10^5 colunas INV + 2,7·10^7 colunas
exteriores num único ponto s, também em HPC. O custo de (c) domina porque
cada coluna da casca pede 32 aplicações, não 1–8 como em RT.

**Ressalvas da projeção (não mudam a decisão, que tem folga de ~100×).**
- Os condicionamentos vêm da truncagem 12×36.
- O ×3 de custo certificado é otimista.
- O acoplamento L4·L5 com o exterior não foi medido numa caixa grande.
- Uma rota que reduzisse o grau k (melhor pré-condicionador de casca)
  ou certificasse a casca por uma cota analítica (k=0, estilo RT) baixaria
  o custo por região para ~0,5 h. Isso tornaria (c) viável, mas exige
  colunas de (B+s)Q com norma ≲0,6 já perto do núcleo. Hoje elas chegam a
  ≈58 (n~20) e ≈1,5 em n=128.

### 6. Implicações para os objetivos

- **Obj. 1** (‖C-C_N‖ uniforme e parametriz de U): a parametriz certificada
  da faixa U, necessária para GUARDED e para a ponte sharp, é inalcançável
  nesta máquina. A contribuição local possível é o exterior barato de (c),
  item 3, que reduziria as caudas analíticas da rodada 2 de 0,998 para
  ≲0,7.
- **Obj. 2 e 4** (testemunho S3 com H2, raiz 0,401 e raiz B): o
  enclausuramento do autopar exige a inversa na caixa com
  ‖QH^-1‖~10^3 perto das raízes. Inalcançável aqui, coerente com
  `CONSTRAINT_WITNESS_NORMS.md`.
- **Obj. 3** (winding certificado): a rota (c) é a arquitetura certa, porque
  troca ~5·10^4 centros por ~150 regiões e um determinante finito 612×612.
  Mas precisa de HPC (~10^5 h-núcleo) ou de um pré-condicionador de casca
  que dispense o grau 32. A pesquisa matemática de maior alavancagem passa
  a ser esse pré-condicionador. HPC fica como alternativa se houver
  alocação.

### Nota de reprodutibilidade (17/09)

`g_parametrix_contour_norms.py` foi reexecutado como script salvo. As normas
e os defeitos coincidem com a execução interativa anterior; ela deixou
`g-parametrix-resolvent.json` e `g-parametrix-shell.json`, agora
redundantes.

O micro-benchmark variou entre as duas execuções (máquina compartilhada):

| grade | 1ª execução | reexecução |
|---|---:|---:|
| 64×200 | 437 ms | 304 ms |
| 107×384 | 979 ms | 659 ms |

Com o valor mais favorável, o custo da casca em (c) cai de ≈255 h para
≈180 h de relógio por região em 2 núcleos. Ainda viola 2 h por ~90×. A
decisão não muda.

## 7. Recalibração de relógio — 17/09/2026, tarde (seção datada)

Todos os tempos das seções 3–5 foram medidos com a máquina a ~1/3 do clock:
`/etc/tlp.conf` fixava `CPU_ENERGY_PERF_POLICY_ON_AC=power`, que põe o HWP
`energy_performance_preference` em `power`, e o governor pedido
(`schedutil`) não existe no `intel_pstate` *active*. Corrigido e persistido
em `/etc/tlp.conf`; depois do reinício os 4 CPUs estão em `performance` e
em 3100 MHz (antes: 900–1000 MHz sob carga, 39 °C — não era térmico).

Os números abaixo são **rederivados por medida nova**, não reescalados pelo
clock nominal. Script: `scripts/clock_recalibration.py`; relatório
`build/spectrum/clock-recalibration-2026-09-17.json`. O JSON antigo
(`g-parametrix-contour-norms.json`) foi preservado como referência "antes".

**O clock subiu ~3,3×; o custo não caiu 3,3×.**

| caminho | benchmark | antes | agora | ganho |
|---|---|---:|---:|---:|
| float (FFT) | 1 aplicação de B, 64×200 | 0,304 s | 0,197 s | **1,54×** |
| float (FFT) | 1 aplicação de B, 107×384 | 0,659 s | 0,454 s | 1,45× |
| exato (núcleo C) | 12 linhas da R1 do A, 1 worker | 45 s | 17,5 s | **2,57×** |

O caminho float é limitado por banda de memória, não por relógio, e por isso
fica bem abaixo do clock; a aritmética exata, limitada pela ALU inteira,
chega perto. A medida do núcleo C foi repetida com a cadeia irmã competindo
(18,1 s contra 17,9 s sozinha), então o ganho não vem de ausência de
contenção. Normalizei por `depth`, já que a medida usou `depth` 61/62 para
não tocar nos checkpoints de produção.

**Custos rederivados** (base do texto principal: 0,44 s/aplicação, da 1ª
execução do benchmark; agora 0,197 s, melhor de 5 com a máquina ociosa):

| quantidade | §3–5 | rederivado |
|---|---:|---:|
| (c) casca por região, float | 170 h | **76 h** |
| (c) casca por região, certificado (×3) | 510 h | **228 h** |
| (c) relógio por região, 2 núcleos | 255 h | **114 h** |
| (c) violação do critério de 2 h | ~130× | **~57×** |
| (b) por centro, certificado | 3,7·10³ h | 1,7·10³ h |
| exterior estilo RT, exato | 0,04 s/coluna; ~3 h | 0,016 s/coluna; **~1,2 h** |
| HPC para (c) | 1–1,5·10⁵ h-núcleo | **3,4–5,8·10⁴ h-núcleo** |

**A decisão da §5 não muda: INVIÁVEL SEM HPC.** A arquitetura (c) continua
falhando só no tempo, agora por ~57× em vez de ~130×. Em 1000 núcleos, a
rota custaria 1,4–2,4 dias em vez de 4–6.

**Kaggle não resolve.** A conta gratuita dá 4 núcleos reais (contra 2 aqui)
e 30 GB, com teto de 12 h por sessão. Uma região certificada de (c) custa
228 h-núcleo → 57 h de relógio em 4 núcleos: ainda 29× acima do critério de
2 h, e acima do teto de sessão. O ganho por núcleo que havia sido medido
(≈2,9× em aritmética de `Fraction`) era contra a máquina estrangulada; com o
EPP corrigido ele cai para ≈1,1×, ou seja, paridade por núcleo. O que o
Kaggle oferece é 2× em núcleos, 4× em RAM e imunidade a queda de energia —
não uma mudança de ordem de grandeza.

**O que isto reforça.** A alavanca continua sendo matemática, não de
hardware: um pré-condicionador de casca que dispense a Neumann de grau 32
cortaria o custo por região para ~0,5 h e tornaria (c) viável aqui. Nenhum
ganho de máquina disponível chega perto de 57×.
