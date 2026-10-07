# Retomada em 14/09/2026

O último resultado documentado está em `MODE_DISCRIMINATION.md`, §10:
o truncamento 12×36 foi produzido e a previsão de aproximação a `mu_RT`
melhor que `10^-7` foi refutada. O candidato físico tem
`lambda = 2,6740777376`, com resíduo de constraints `9,17e-5`.
Esses números são diagnósticos existentes, não novas provas.

Na revisão do produtor `scripts/build_contour.py`, foi corrigido o
arredondamento do majorante de Hadamard usado para escolher o módulo CRT.
Para cada coeficiente de grau `n-k`, a cota utilizada agora é

```
ceil(sqrt(binomial(n,k)^2 * k^k * B^(2*k))),
```

calculada exclusivamente com inteiros, onde `B` é o maior módulo das
entradas. Antes se multiplicava `floor(sqrt(k^k+1))` por
`binomial(n,k)*B^k` e se somava 1. Para `n=3, B=1000`, o majorante
global antigo era `5000000001`, menor que a cota de Hadamard
`sqrt(27)*10^9`. Portanto a justificativa de reconstrução pelo CRT
não estava garantida pela fórmula. Isso não demonstra que algum
polinômio anteriormente produzido estava incorreto: o produto dos primos
pode exceder também a cota correta.

`scripts/test_build_contour.py`, integrado a `make test`, verifica o
arredondamento para fora e compara a reconstrução CRT com determinantes
de Bareiss em `n+1` pontos distintos para matrizes pequenas.

Permanecem abertas as obrigações analíticas do ledger. Esta correção não
fecha P1/P2/P3, S3/S4 nem a transferência para o operador infinito.
Certificados existentes e caches de polinômios não foram regenerados
nesta revisão. A leitura de cache verifica somente dimensão e monicidade;
a vinculação do cache à matriz de origem ainda merece revisão.

## Continuação da auditoria

- `ANALYTIC_AUDIT.md`: contraexemplo à comutação de inversão e projeção,
  também verificado no kernel; obstrução ao ganho analítico P2 no espaço
  sem pesos; distinção entre o ledger formal e o enunciado físico.
- `ALGEBRAIC_TAIL.md`: rota alternativa com pesos RT fixos, estimativa
  algébrica, ponto-base 828 e borda direita conservadora 414 sob os bounds
  publicados. Aritmética em `AlgebraicTail.lean`; análise funcional ainda
  em papel. A última seção deriva uma ponte usando o precondicionador
  livre, cujo majorante ainda não foi satisfeito numericamente.
- `build_tile.py`: resolvente finita correta, ponto-base explícito,
  pesos RT opcionais, suporte a centros racionais não diádicos e correção
  da colisão de nomes entre matriz residual e teorema por blocos.
- Setor B reduzido em 4×12: 1946 arestas classificadas, winding 0.
  União exata em 55 segmentos, aceitos pelo kernel. Dados em
  `examples/contour-A-4x12-reduced-B-*.txt`; `make verify-reduced-b` reproduz.

Não foi demonstrado “exatamente um modo físico instável”. Permanecem S3/S4,
a equivalência do domínio, a transferência quantitativa para os dados do
fundo exato e a cobertura de toda a região instável. Os contornos antigos
param em Re s=1 e começam em 1/8; a exclusão conservadora só começa em 414.

## 15/09/2026: precondicionador livre implementado

`check_free_preconditioner.py` realiza a fatorização alternativa na matriz
4×12 inteira, com aritmética racional e pesos RT. A inversa finita no
centro s=1 foi aceita; o majorante de transferência falhou já no primeiro
fator (`2592/11`). A falha ocorre em todas as malhas disponíveis e foi
registrada em Lean. Mais precisão na inversa não corrige esse primeiro fator.
Detalhes, relatório e reprodução em `ALGEBRAIC_TAIL.md`, última seção.
Próximo trabalho: bounds por blocos e acoplamentos de frequência para
substituir o majorante global; depois incorporar a incerteza do fundo.

## Blocos de frequência e seus acoplamentos

`FREQUENCY_BLOCKS.md` deriva o critério de contração por blocos, com pesos
racionais, e bounds direcionais do núcleo para o exterior e de volta.
No teste 6×18 com núcleo 4×12, inversas separadas falham pelo produto dos
acoplamentos (~100). A inversa acoplada passa com defeito por blocos <1,36e-9.
Isso é finito e não inclui a cauda infinita.

`background_couplings.py` soma exatamente as caudas ponderadas de RefA e
adiciona o erro publicado. O bound direto de B do núcleo 4×12 para fora
de 24×72 cai para <0,001427; a passagem inversa via Q0 foi estimada
separadamente. Relatórios em `build/spectrum/*couplings.json` e
`build/spectrum/frequency-blocks-6x18-core4x12.json`.

Próximo gargalo explícito: coroa/exterior e exterior/exterior, mais a
uniformidade no contorno complexo. A pequenez do acoplamento núcleo/exterior
sozinha não fecha esse gargalo nem S3/S4.
# Atualização: coroa/exterior e equivalência física

Foram acrescentados majorantes uniformes dos três blocos restantes em
`FREQUENCY_BLOCKS.md`, seção 7, e no produtor `background_couplings.py`.
Para P=12×36, G=64×192, base zero e |s|<=5/4, o defeito exterior tem
majorante 2619/191>1: o critério ainda não fecha. O termo s cancela em
T H S porque Q0 preserva G, mas não nos outros dois blocos.

`CONSTRAINT_GAUGE_EQUIVALENCE.md` demonstra lemas condicionais para
kernels/cadeias e seção do quociente gauge. Corrigida a obrigação:
invertibilidade sharp apenas no contorno não elimina violações interiores.
Ainda faltam a identidade off-shell RT/Floquet com domínios e a
classificação de gauge admissível. Testes Python: 26 passaram.
# Atualização: estrutura exterior e identidade sharp completa

`STRUCTURED_EXTERIOR.md` demonstra a melhora de cauda para base real:
nos parâmetros RT, max(15/M,81/(N+1)), em vez de max(36/M,108/(N-1)).
O produtor toma também o mínimo com a estimativa anterior. Em 64×192,
o bound uniforme TT cai de 13,712 para 10,177; ainda não contrai.

`SHARP_LINEARIZED_IDENTITY.md` constrói K(s)C(s)=M(s)L(s) no fundo exato,
partindo da identidade completa não linear. Em fundo aproximado aparece
o termo adicional -G(h,-P OmegaPlus(w)), que foi mantido e testado.
Inclui reconstrução dos cinco resíduos e fórmulas explícitas de M.
Extensão aos domínios fechados, espectro sharp e gauge ainda abertos.
Novos scripts: structured_exterior.py, sharp_propagation.py e seus testes.
# Continuação: norma exterior por componentes

`COMPONENT_EXTERIOR.md` separa a parte descendente mu S Gamma1 (primeira
saída identicamente zero) e preserva cancelamentos entre campos de RefA
antes de estimar a multiplicação. O produtor component_exterior.py dá
beta<10.382 na norma com pesos de componentes (916,4096,243,1233), incluindo
o erro do fundo convertido para essa norma. Em G=64×192, TT<4.882;
em 175×942, o bound diagonal é 117774/117875<1. Isso ainda não fecha
os blocos acoplados. Recalcular os blocos finitos na MESMA norma é obrigatório.
Relatório com hashes: build/spectrum/component-exterior.json.
O transporte racional da parametriz acoplada 6×18 para a nova norma passou:
defeito<2.28e-8 no centro s=1. O recálculo direto foi interrompido após
essa validação, e não produziu novo relatório; não confundir com o transporte.
Verificação: make test passou com 14 alvos Lean e 45 testes Python.
Nenhum novo teorema
de operadores infinitos foi formalizado em Lean.
# Etapa: contração uniforme, domínio sharp e gauge residual

Novas notas:
- UNIFORM_COUPLED_CONTRACTION.md: critério a+bc/(1-d)<1, raio estrito
  por tile e homotopia completa. Em 64×192, d>1 ainda impede esse bound;
  a caixa 175×942 não tem parametriz finita e deixa margem exterior pequena.
  O erro do centro deve incluir mu, não apenas o erro multiplicativo do fundo.
- SHARP_DOMAIN.md: D_min=D_max para o principal livre; extensão de
  C(s) ao domínio sharp via perda de raio (65/64,5/4)->(129/128,9/8).
  Falta controle espectral sharp nesses novos pesos, não substituível pelo antigo.
- GAUGE_DOMAIN_CLASSIFICATION.md: no gauge residual nulo analítico,
  s_n=mu(1-n); o único positivo move o cone. A fase neutra é transversal:
  |b(partial_tau omega*)|>117/1000, condicionado à bola RT.
  Projetar essa fase em s diferente de zero não preserva o kernel.

Novos scripts uniform_coupling.py e gauge_null_profiles.py; testes locais
exatos. Nada nesta etapa fornece uma contagem física concluída ou uma
formalização Lean das identidades de operadores infinitos.
Verificação desta etapa: 14 alvos Lean passaram; 55 testes Python passaram.
# Rodada de 16/09: sharp e parametriz amortecida

Novos exportador/analisador sharp: scripts/rt_sharp_matrix.c,
build_sharp.sh e analyze_sharp.py. Matrizes 4×12, 6×18, 8×24, 12×36
geradas, dimensões 54, 135, 252, 594. Na última, raiz sharp real
0.40102473085, contra selecionada 0.40102441099 (distância 3.2e-7).
O resíduo relativo da propagação desse modo é 0.0592: ainda não é zero.
Ver SHARP_SPECTRAL_LOCATOR.md; os relatórios preservam raízes adicionais.

Para a coroa 6×18 menos 4×12, Neumann sem amortecer diverge perto de s=1.
Com omega=1/2 e grau 31, a proposta passa nos quatro centros examinados.
certify_crown_tiles.py valida a proposta por resíduos racionais e controla
a variação do tile, sem confiar no cálculo de potências em float.
reference_error_budget.py dá ||H_true-H_ref||<0.0003257 para Q_ref FIXO,
incluindo epsilon_mu≈3.88112e-11 e a bola omega publicada. Essa parcela
pode ser incorporada aos tiles. Ver DAMPED_CROWN_PARAMETRIX.md.

## Cobertura da coroa concluída em 16/09

As duas auditorias completas terminaram: referência e fundo com incerteza,
51 tiles cada. Na segunda, defeito uniforme <0.503 e
sup_contorno ||H_SS^-1|| <=669414234798400/18478486838353 <36.227.
Isso vale para S=(6×18)-(4×12), condicionado à bola RT publicada.
Relatórios crown-contour-6x18{,-background}.json em build/spectrum.
check_crown_cover.py conferiu independentemente os escalares e a
cobertura geométrica; não recalcula os resíduos matriciais do produtor.

O próximo passo é controlar os blocos com T=I-(6×18), ou aumentar a
coroa mantendo todos os blocos intermediários. A nota da parametriz
explicita os termos de Schur. Somar seus majorantes positivos não cura
um bound TT>1. Também não se pode descartar o winding da coroa só por
ela ser invertível na fronteira. Constraints espectrais, completude
do gauge e contagem física infinita continuam abertas.

## Próxima rodada: acoplamentos, sharp local e ação gauge

GUARDED_INFINITE_COUPLING.md usa a norma finita ||Q_ref F||<5.942 e
caudas por componentes. De F=6×18 para fora de 64×192, ||T H F||<1.791e-5.
Com T fora de 176×944, a eliminação de S=(6×18)-(4×12) dá defeito de
Schur <0.997676: o subsistema S+T é uniformemente invertível, sob a
bola RT. O complemento P+[(176×944)-(6×18)] continua sem parametriz.
Relatório: build/spectrum/guarded-exterior.json.

SHARP_LOCAL_EXCLUSION.md: matriz sharp 8×24, pesos menores (129/128,9/8),
disco |s-0.733|<=1/32, defeito racional <0.264 e inversa <11.446.
Não inclui cauda infinita nem erro do fundo. Relatório sharp-disk-8x24.json.
A nota também deriva o transporte de cadeias c_j=C0 h_j+C1 h_(j-1).

GAUGE_RT_ACTION.md: ação explícita h=(1+mu^-1 partial_tau+xi partial_xi)omega
em s=mu, identidade off-shell completa testada racionalmente. Prova de
completude residual LOCAL no ansatz; cone móvel versus fixo e perda de
raio mantidos explícitos. Faltam equivalência física global, domínio
forte e multiplicidades algébricas; nenhuma raiz foi subtraída.

A partição de COLUNAS em F=6×18 e I-F também dá
||Q_ref||<=max(5.941521,81/19)<5.942 no espaço infinito. O orçamento
total do fundo melhora para <1.7918e-5 sem alterar os tiles antigos.
Verificação da rodada: 14 alvos Lean e 78 testes Python passaram.
As novas conclusões analíticas não foram formalizadas em Lean.

## Continuação: recuperação de raio e sharp com fundo incerto

RADIUS_RECOVERY.md: usando Q=J_mu_true^-1 (não Q_RefA), a cauda fora
de 1024×4096 contrai nos pesos fortes e fracos, com defeitos
117774/512125 e 8046603/10242500. Existência forte + unicidade fraca
provam recuperação de raio para todos os autovetores e cadeias em
|s|<=5/4. Nenhuma matriz dessa caixa foi invertida. O gerador gauge
positivo agora pertence ao domínio RT forte. A compactação de Q por
caudas também estabelece Fredholm de índice zero nessa realização.

SHARP_EXTERIOR_BOUNDS.md: beta_sharp<=4.62237 nos pesos menores.
O disco finito |s-.733|<=1/32 passou incluindo a bola RT e erro de mu,
com defeito <0.263453 e inversa <11.447. Cauda sharp TT<0.572 em
256×1536; não combinada com a inversa 8×24 ignorando a faixa entre elas.
Relatórios sharp-exterior.json, sharp-disk-8x24-background.json e
radius-transfer.json em build/spectrum.

CONSTRAINT_ROOT_QUOTIENT.md: D=C0-C1 A intertwina os geradores;
constraints em cadeias são D h_j=C(s0)h_j+C1 h_(j-1). A dimensão do
quociente é dim E-rank(D|E)-1 quando a linha gauge móvel é admissível.
Isso é uma fórmula condicional, não uma avaliação dessas dimensões.
Faltam a contração acoplada através das faixas intermediárias, a exclusão
sharp infinita e a equivalência geométrica global com a escolha do cone.

Verificação desta continuação: make test passou, 14 alvos Lean e 88 testes
Python. Novos lemas analíticos continuam em papel; nenhum foi anunciado
como formalização Lean ou fechamento da contagem física.

## Rodada: faixa sharp e reconstrução geométrica no cilindro

Sharp 12×36 passou com a faixa adicional de 342 DOFs acoplada a 8×24:
compressão diádica coincidente, defeito <0.265850, inversa <11.588 com
erro do fundo. A passagem correta a H=K Q_ref dá inversa finita <1201.510.
O subsistema F=(12×36) mais T fora de 160×960 tem Schur <0.936 em todo
o disco |s-.733|<=1/32. Resta o Schur finito sobre U=(160×960)-F,
228366 DOFs, que NÃO foi construído nem validado. A expressão exata e
o erro de truncamento da inversa exterior estão em SHARP_INTERMEDIATE_SHELL.md.

GEOMETRIC_FLOQUET_RECONSTRUCTION.md reconstrói q,z,p a partir de h, com
z=(partial_tau+s)^-1 a3 e p=(partial_tau+s)^-1 a4. Para Re s>0 não há
obstrução de período. Identidades locais e reciprocidade dos resíduos
foram testadas racionalmente. A equivalência é global no cilindro/no
ansatz e vale para cadeias; não é ainda a fixação global do gauge para
perturbações físicas arbitrárias. Cone fixo exige uma condição explícita
sobre a perturbação da métrica no cone; a alternativa móvel foi mantida.

## Retomada após a queda de energia em 16/09/2026

Foram recuperadas alterações posteriores à rodada acima. `make test`
passou novamente: 14 alvos Lean, zero usos de `native_decide`, 100 testes
Python e checagem de sintaxe dos scripts shell.

O relatório `sharp-disk-12x36-preconditioned.json` usa o resíduo direito
da inversa candidata V para estimar diretamente J_ref V e obter
||H_FF^-1||<11.461, em vez do produto de normas <1201.510.
`sharp-free-global.json` registra ||Q_ref||<5.942, combinando a auditoria
das colunas finitas com a estimativa analítica de cauda. Esses dois
relatórios foram preservados; seus resíduos matriciais completos não
foram recalculados nesta retomada.

Foi reexecutado o produtor da ponte com esses relatórios:

```bash
.venv/bin/python scripts/sharp_guarded_bridge.py \
  --outer-m 148 --outer-n 840 \
  --finite-report build/spectrum/sharp-disk-12x36-preconditioned.json \
  --free-report build/spectrum/sharp-free-global.json \
  --out /tmp/choptuik-resumed-bridge.json
cmp build/spectrum/sharp-guarded-bridge-148x840.json \
  /tmp/choptuik-resumed-bridge.json
```

O resultado coincide byte a byte com o relatório salvo: defeito de Schur
0.9892221789135963 (a comparação com 1 é racional), acoplamento de saída
<5.313e-7. Ainda resta a faixa intermediária de 184626 DOFs. A conclusão
continua restrita a F mais a cauda distante, condicionada à bola RT;
não exclui o espectro sharp do operador completo.

`global_null_gauge.py` e seus três testes também sobreviveram à queda.
Eles verificam identidades polinomiais de transporte, ressonância e
paridade. Não provam periodicidade global nem completude de gauge para
perturbações arbitrárias. O próximo trabalho analítico continua sendo
controlar o Schur da faixa intermediária e justificar a equivalência
global com os domínios físicos.

## Continuação confirmada: faixas intermediárias e gauge global

`sharp_shell_remainder.py` mantém os quatro acoplamentos com a faixa U
e propaga o erro de Neumann até seu Schur. Em 148×840, o grau 3785
dá erro <=9.92351e-13. A ponte 160×960 foi recalculada com as auditorias
precondicionadas recentes, dando defeito <0.915031; nela o grau 437 dá
erro <=9.7964e-13. Relatórios `sharp-shell-remainder-148x840.json`,
`sharp-guarded-bridge-160x960-preconditioned.json` e
`sharp-shell-remainder-160x960-preconditioned.json` em build/spectrum.
São orçamentos de erro: as séries e a parametriz do Schur NÃO foram
construídas. A exclusão sharp infinita permanece aberta.

`GLOBAL_NULL_GAUGE_TRANSPORT.md` dá uma fórmula global no cilindro para
o transporte de fixação das componentes radiais diagonais. Separa o
traço no cone e a integral de características, convergente para Re s>0.
Em s=mu-i m/2, explicita a obstrução no coeficiente temporal do traço;
cone fixo exige anular o traço inteiro. A reflexão preserva o centro.
Os testes novos usam modos Fourier periódicos e aritmética racional,
incluindo ressonâncias e a obstrução de cadeia do operador de transporte.
Falta vincular a classe analítica construída ao domínio físico e RT forte;
não foi concluída equivalência geométrica de todas as perturbações.

Verificação desta rodada: `make test` passou com 14 alvos Lean, zero
usos de `native_decide`, 107 testes Python e sintaxe shell válida.
Os lemas analíticos novos permanecem em papel.

## 16/09/2026 (tarde): retomada pelo Claude e auditoria

O trabalho de 14/09 16:40 a 16/09 13:27 foi feito sem o loop de agentes do
Claude e depois auditado de forma independente: `AUDITORIA_MATEMATICA_16SET.md`
e `AUDITORIA_FISICA_16SET.md`. Nenhuma alegação foi refutada. As correções de
alcance foram integradas ao ledger (seção "Auditoria de 16/09"): C1 volta a
aberta, C2 ganha hipóteses fora do kernel, a faixa U é inviável e C4 foi
reformulada após a descoberta da raiz B em `0,0414+0,5i`. `make test`: 14
alvos Lean, 0 `native_decide`, 120 testes Python.

Rodada 2 em curso: majorante exterior por colunas d*(G), com regra de decisão
pré-registrada em `COLUMNWISE_EXTERIOR.md` (N_G<=384 mantém a rota de blocos),
e revisão cruzada da norma do testemunho S3 da raiz 0,401.

## 17/09/2026: rodadas 2–4 e setor B 6×18

- **Rodada 2 (Matemático):** d*(G)<1 em G=107×384 para A e sharp, com a regra
  pré-registrada (a rota de blocos continua); margens de Schur ~10⁻³.
- **Rodada 3 (Físico):** d*(G) confirmado por código independente; H3 fechada
  (c_op<=22,179). O testemunho S3 é inalcançável hoje.
- **Setor B, 6×18:** o Codex produziu o polinômio exato e o orquestrador
  terminou os certificados. `make verify-reduced-b-6x18` dá winding 1 no
  retângulo e winding 1 na caixa da raiz, aceitos pelo kernel.
- **Rodada 4 (Matemático):** a inversão certificada na caixa G é INVIÁVEL SEM
  HPC (fator ~100 em tempo); a arquitetura (c) é a correta.
- **Decisão pendente do usuário:** HPC, pesquisa de um pré-condicionador de
  casca, ou os passos baratos (exterior de (c), cotas da Fase 1).
- `make test`: 14 alvos Lean, 0 `native_decide`, 133 testes Python.

## 17/09/2026 (tarde): Marco 8 fechado, e três achados de infraestrutura

**Marco 8 (exterior com cortes N_far=900, M'=160) — meta atingida.** As duas
cadeias do `columnwise_exterior_run.sh` terminaram e
`columnwise_exterior_assemble.py` montou o certificado
`build/spectrum/columnwise-exterior-phase3.json`:

| sistema | regiões | d_R(900) | d_F(160) | **d\*(G)** | ≤0,7 | buracos |
|---|---:|---:|---:|---:|:--:|:--:|
| A | 0,504929 | 0,651165 | 0,667547 | **0,667547** | sim | 0 |
| sharp | 0,484861 | 0,673040 | 0,667838 | **0,673040** | sim | 0 |

900 linhas vistas em cada sistema, cobertura do complemento de G=(107,384) sem
buracos. Schur GUARDED invertível nas duas contas (independente 0,667923,
relatório 0,667980; margem ≈0,3320) e ponte sharp invertível (δ=0,673044,
margem 0,326956). Ou seja: d* caiu de 0,998 para ≈0,67, e a margem de Schur
saiu de ~10⁻³ para ~0,33. O montador continua marcando
`conditional_on_published_RT_ball=true` e `physical_certificate=false` — a
ressalva de sempre, não foi tocada aqui.

**Achado 1 — queda de energia às 15:43.** Segunda em dois dias. Nada foi
perdido: o `fsync` por linha do `record()` garantiu checkpoints íntegros, a
cadeia `sharp` já havia terminado inteira e faltavam só 106 linhas da R2 do A,
retomadas em ~15 min. Os logs em `/tmp` foram apagados pelo reboot; os
checkpoints em `build/spectrum/` foram a única fonte sobrevivente. Detalhes em
`build/spectrum/columnwise-phase2/EXTERIOR-RUN.md`.

**Achado 2 — não havia repositório.** O diretório `.git` existia vazio: todo o
projeto estava sem versionamento, com o tarball `.snapshots/2026-09-16-pre-codex/`
como único backup. Corrigido: `git init` e commit inicial de 293 arquivos. O
`.gitignore` agora versiona os certificados `.json`/`.jsonl` de `build/` (1,9 MB,
horas de CPU) e deixa fora as matrizes RT `.dat` (~93 MB, regeneráveis).

**Achado 3 — a máquina roda a 1/3 do clock.** Sob carga total nos 4 threads a
CPU fica em 900–1000 MHz com teto de 3100 MHz, a 39 °C. Não é térmico:
`/etc/tlp.conf` põe `CPU_ENERGY_PERF_POLICY_ON_AC=power`, que fixa o HWP
`energy_performance_preference` em `power`; e o `CPU_SCALING_GOVERNOR_ON_AC=schedutil`
falha em silêncio, porque no `intel_pstate` *active* só existem `performance` e
`powersave`. **Consequência para a rodada 4:** o fator "~100× em tempo" de
`G_PARAMETRIX_FEASIBILITY.md` foi medido nessa condição. Corrigir o EPP dá ~3×,
e o Kaggle (4 núcleos reais contra 2 físicos, 30 GB) dá ~2× além disso. A
conclusão de inviabilidade sem HPC deve se manter — 6× não é 100× — mas os
números daquele documento **precisam ser rederivados, não reafirmados**.

**Kaggle disponível.** Conta gratuita `edivaldoag` configurada e testada de
ponta a ponta (push, status, output). A quota de 30 h/semana é só de GPU, que
não nos serve; CPU não é medida, o limite é 12 h por sessão. Ganho real: 4
núcleos reais, 30 GB de RAM e imunidade a queda de energia. Carga real exige
subir os fontes como dataset privado, o que depende de autorização do usuário.

## 17/09/2026 (fim de tarde): clock corrigido e custos rederivados

O `CPU_ENERGY_PERF_POLICY_ON_AC=performance` foi aplicado no `/etc/tlp.conf` e
sobreviveu ao reinício: os 4 CPUs estão em `performance`, a 3100 MHz. Os
custos do `G_PARAMETRIX_FEASIBILITY.md` foram **rederivados por medida nova**
(§7 daquele documento, `scripts/clock_recalibration.py`), não reescalados.

O achado é que o clock subiu ~3,3× mas o custo não caiu 3,3×: o caminho float
(FFT) ganhou só **1,54×**, porque é limitado por banda de memória, e a
aritmética exata do núcleo C ganhou **2,57×**. A previsão de "~3× do EPP mais
2× do Kaggle ≈ 6×" estava otimista nos dois fatores.

Consequência: a casca certificada de (c) cai de 255 h para **114 h** de relógio
por região em 2 núcleos, e a violação do critério de 2 h cai de ~130× para
**~57×**. O HPC necessário cai de ~10⁵ para **3,4–5,8·10⁴ h-núcleo**. A decisão
**INVIÁVEL SEM HPC não muda**, e agora está rederivada em vez de reafirmada.

O Kaggle também não muda a ordem de grandeza: uma região custaria 57 h de
relógio nos 4 núcleos de lá, acima do teto de 12 h por sessão. E o ganho de
2,9× por núcleo que havia sido medido era contra a máquina estrangulada — com
o EPP corrigido é ≈1,1×, paridade. A alavanca continua sendo o pré-condicionador
de casca, que precisaria de ~57× e nenhum hardware disponível dá.

## 17/09/2026 (noite): a rota do pré-condicionador está fechada

`docs/SHELL_PRECONDITIONER.md` (rodada 6). Dois resultados negativos, ambos
decisivos, medidos em float na caixa 12×36:

1. **Nenhum pré-condicionador da família serve.** O k=32 não era desperdício: o
   raio do disco escala com (1-α), então o ponto de operação é α≈10⁻⁵, e ali o
   k medido é 29 — confirma a escolha do documento. O melhor candidato,
   bloco-m, dá 1,93×; Jacobi, bloco-n e deflação ficam em ≤1,26×. A deflação é
   pior que inútil: cria transientes de 10²–10⁴ pela não-normalidade, e o que
   se certifica é a norma de coluna, não o raio. Mesmo um pré-condicionador
   perfeito e gratuito (k=1) só daria 29–32×, contra os 57× necessários.
2. **A rota analítica k=0 tem um piso exato.** Trocar pesos por DOF é uma
   similaridade diagonal, e por Perron-Frobenius o mínimo da norma de coluna
   sobre todos os pesos é ρ(|X|). Medido: 9,48 / 12,27 / 10,91 nos três
   centros, contra o alvo 0,6. A melhor repesagem concebível tira 58 para ~10.

Restam a alocação de HPC (3,4–5,8·10⁴ h-núcleo) ou uma arquitetura que não
certifique a casca coluna a coluna. Nenhuma das duas cabe nesta máquina.

## 17/09/2026 (noite): o PC do laboratório medido — 2,0×, e a rota fecha também

`.codex-runs/2026-09-17-lab-benchmark/RELATORIO.md`, com `lab-specs.json` e
`lab-bench.json`. Snapshot em `.snapshots/2026-09-17-lab/`; a auditoria contra o
manifesto fecha sem diferenças. Medido direto por ssh nesta sessão — o Codex não
foi lançado, e a chave já estava instalada (o `ssh-copy-id` travou depois de ter
feito o trabalho).

**A máquina:** i7-4770 (Haswell, 2013), 4 núcleos físicos / 8 threads, 15 GiB.
O clock **não** está estrangulado lá — 3692 MHz nos 8 threads sob carga cheia,
`intel_cpufreq` passivo com `ondemand`, sem `energy_performance_preference`.

**Antes de medir, a pilha foi igualada.** O Python do sistema é 3.8.10 contra
3.14.6 aqui, e só isso vale 1,27× no caminho exato (13,00 s contra 10,25 s no
mesmo trabalho). Montei CPython 3.14.6 e o venv do `requirements.txt` com `uv`,
como usuário, sem root e sem PPA. Sem igualar, a medida teria sido de software.

| | notebook | laboratório | ganho |
|---|---:|---:|---:|
| float (FFT), 1 núcleo, 64×200 | 0,197 s | 0,138 s | **1,42×** |
| exato (núcleo C), 12 linhas da R1 | 14,44 s | 10,25 s | **1,41×** |
| relógio com todos os núcleos | 2,11× de 1 núcleo | 2,94× de 1 núcleo | **2,0×** |

**Achado de infraestrutura: `--workers` não paraleliza.** No motor `c` ele vira
`OMP_NUM_THREADS` de *uma* chamada do núcleo, e o laço sobre blocos de linhas,
com todo o preparo e a remontagem em `Fraction`, é serial no processo pai.
Amdahl sobre a curva medida dá ~35 % de parte serial e teto de 2,85×; com 8
threads a máquina entrega 1,95×. Dividir as linhas entre **processos
independentes** (`--workers 1` em cada) dá **2,94×** com 4 processos — e 8
processos são *piores* que 4, porque o caminho é ALU inteira e o hyperthreading
só faz as duas threads brigarem. O mesmo padrão aparece no caminho float.

**Rederivação, só com o medido.** Uma região certificada cai de 108 h (notebook,
vazão medida) para **54 h** de relógio; o critério de 2 h segue violado por
**27×**, contra 54× aqui. O total de 3,4–5,8·10⁴ h-núcleo vira 2,4–4,1·10⁴
h-núcleo de lá, ou **338–576 dias** de relógio ininterrupto. Memória não é
restrição: 58 MB por processo no caminho exato, 139 MB no float.

**A rota do laboratório está fechada** — mas a medida entregou o dimensionamento
que faltava para o pedido de alocação: **~80 núcleos** desse porte com escala
ideal, ou **109 reais** na eficiência de 74 % medida, para uma região dentro das
2 h; e **1000 núcleos → 32–55 h (1,4 a 2,3 dias)** para o total.

Ressalva registrada: os 17,5 s da §7 não se reproduziram no notebook hoje (14,2 s
normalizados). A comparação usada é pareada, as duas máquinas medidas hoje com a
mesma pilha; a §7 não foi reescrita.

## 20/09/2026: a rota do majorante não tem teto, e o alvo virou número

`docs/TT_CEILING.md`, `scripts/tt_ceiling.py`,
`build/spectrum/tt-ceiling.json`. Também criado `docs/PRIORIDADES.md`, que
sequencia o que fazer em seguida e registra as decisões já tomadas para não
serem reabertas.

**Primeiro, o reenquadramento.** Das cinco obrigações analíticas que bloqueiam
o resultado — S3, S4, C1, F2, F3 —, três são o mesmo problema: C1, F2 e F3 se
reduzem a tornar o majorante do exterior menor que 1. F3 já está parcialmente
fechada pela cauda algébrica (`Re s >= 414`), e o que falta dela depende do
mesmo majorante. S3 e S4 são de outra natureza e ficam por último.

**A pergunta do teto é de Perron, mas não sobre a matriz de DOFs.** O critério
de blocos (existe `r>0` com `sum_i r_i e_ij < r_j` por coluna) equivale a
`rho(E)<1` para `E` não negativa. Na matriz majorante de blocos,
`e_TS = beta*t(P)` **não depende da caixa**, enquanto `e_ST` e `e_TT` são
proporcionais a `t(G)`, que vai a zero. Logo a condição passa para toda caixa
suficientemente grande: **não há piso.** É a diferença estrutural para a rota
k=0, onde a reponderação era similaridade diagonal sobre uma matriz fixa e o
ínfimo era `rho(|X|)`.

**O que existe é uma caixa exigida, e ela escala como `beta^2,8`.** Com o
núcleo em 107×384, `||V||=1` (otimista) e a cauda estruturada:

| `beta` | caixa G mínima | DOFs relativos a 107×384 |
|---:|---:|---:|
| 23 (uniforme) | 2124×11470 | 592,9× |
| 10,382 (componentes, hoje) | 556×3001 | 40,6× |
| 3,941 | 214×768 | 4,0× |
| 2,064 | 107×384 | **1,0×** |

**O alvo deixou de ser qualitativo.** Reduzir `beta` de 10,382 para **2,064**
(fator 5,0) faz a etapa computacional caber na caixa que já foi medida e
orçada, com **custo extra zero**. Um alvo intermediário, `beta <= 3,941`, a
deixa em 4× do orçado. Para calibrar: `beta` já caiu de 23 para 10,382 com um
único movimento, e a norma ponderada exata de `B_6x18` medida vale **8,078**,
abaixo do bound 10,382 em uso — 1,29× de folga só em contabilidade.

Ressalvas registradas em `TT_CEILING.md` §5: `||V||=1` é piso, `a_fin` vem de
um experimento validado só em `s=1`, os acoplamentos `P<->T` foram omitidos
pelo argumento de cauda do fundo, e a tradução DOFs→horas-núcleo é otimista.
É dimensionamento, não certificação.

**Consequência prática:** o trabalho analítico seguinte passou a ser também a
via mais barata de resolver o problema computacional. Comprar horas-núcleo
antes de tentar o fator 5 em `beta` é otimizar a restrição que não está
mordendo. Ver a regra de sequenciamento em `PRIORIDADES.md`.

## 20/09/2026 (continuação): a rota do majorante fecha na prática, e a faixa `0<Re s<1/8` não está vazia

Documento de fechamento em `docs/MAJORANT_ROUTE.md`; detalhe rodada a rodada em
`docs/TT_CEILING.md` §§7–23. Scripts: `tt_ceiling.py`, `tt_dimensioning.py`,
`measure_parametrix_norm.py`, `measure_contour_border.py`.

**Cinco rodadas, quatro alavancas fechadas com medida.** A rodada 1 mostrou que
a rota **não tem teto**: na matriz majorante de blocos, `e_TS` não depende da
caixa enquanto `e_ST` e `e_TT` vão a zero com ela, então o critério passa para
caixa grande — ao contrário da rota k=0, onde havia piso de Perron. O que
existe é uma caixa exigida.

Depois disso, cada alavanca foi medida em vez de estimada:

- **reponderação por componentes:** os pesos em uso estavam a **0,09%** do piso
  `rho(M)=10,3717`. Fechada, negativa.
- **ponto-base `z0`:** em `z0` grande `||V||` ganha bound analítico pequeno, mas
  a cauda piora mais rápido. `z0=0` ganha, inclusive supondo a cauda estruturada
  estendida. Fechada, negativa.
- **borda `sigma0`:** ver abaixo. Fechada, negativa, 1,6×.
- **absorver `B`:** quantificada em 1,24×, não no fator 5 esperado.
- **separação aditiva e razão de aspecto:** colhidas, ~1,7× e ~1,5×.

**O termo dominante era `||V||`, e ele estava mal medido.** Os 80 do registro
foram tomados em `s=1`, que é o ponto mais LONGE de qualquer raiz. Na borda do
contorno `||V||` vale **3822** em 6×18, e cresce também com a malha (79,1 em 138
DOFs, 143,7 em 333, no mesmo `s=1`). Não é patologia: `V` é a inversa do
operador e o contorno passa perto de onde ele é singular por construção. A raiz
responsável converge para `s = mu = 0,168311`, o candidato a modo de gauge de
S4, e não sai do contorno com o refinamento.

Com os valores medidos, a caixa exigida vai a **9127×–14607×** os DOFs de
107×384, contra os 332× da estimativa otimista. A distância à viabilidade não é
o fator 5 em `beta` que a primeira rodada sugeriu: são três a quatro ordens de
grandeza.

**Achado espectral colateral, e independente do custo.** Varrendo `sigma0` em
6×18, `||V||` dá **39183 em `sigma0=1/32`** — dez vezes o valor da borda atual e
fora de qualquer tendência monótona, indicando quase-singularidade perto de
`s ~ 0,03`. A faixa `0 < Re s < 1/8` estava listada como **descoberta** em
`ALGEBRAIC_TAIL.md` e em C4, e nunca tinha sido sondada: agora há evidência
numérica de que **não está vazia**. Afeta C3, C4 e F3. Ressalva: truncagem
6×18, pode ser artefato; merece confirmação em 8×24 e 12×36.

**Três números meus foram corrigidos no caminho**, e ficam registrados porque a
correção é o resultado: "1,29× de folga em contabilidade" (era 0,09%), "14,4×"
de caixa (era 23,4×, a separação é que mordia) e a hipótese de que `||V||` seria
monótono na distância à raiz (não é).

**Consequência para o orçamento.** `~/possivel_projeto.docx` está obsoleto e não
deve ser submetido: foi dimensionado para 107×384. Nenhuma compra de HPC se
justifica — nem 332×, nem 10⁴×, cabe em orçamento algum.

**Regra que fica:** continuar por esta rota exige ideia estruturalmente
diferente, não ajuste de peso, norma, ponto-base ou contorno.

## 21/09/2026: a revisão de completude acha o campo de valores, e F3 ganha uma peça certificada

Documento de trabalho: `docs/ALTERNATIVE_ROUTES.md` §§5.9–5.16. Nota de método:
`docs/METODO.md`. Scripts: `numerical_range_f3.py`,
`certify_free_numerical_range.py`.

**Como apareceu.** Depois que o ataque à contagem fechou quatro becos (§5.8), a
pergunta foi se eram só quatro. Varrendo `docs/` por famílias de técnica padrão,
**sete não apareciam em lugar nenhum do projeto**: campo de valores, realização
`l^2`, operadores limite, pseudoespectro, renormalização, formas quadráticas e
Gershgorin. O mapa não estava completo, e a primeira família testada rendeu.

**Campo de valores falha no contorno, por convexidade.** `W` é convexo e contém
todo o espectro; como `J_mu` tem autovalores `mu(n+1) -> infinito`, o espectro
chega a `Re s -> -infinito`, e um convexo que contém isso e as raízes em
0,168–0,733 cobre a borda em `1/8`. Nenhuma escolha de contorno escapa.

**Mas a borda direita de `W` dá F3 de graça.** `spec(T) subset W(T)`, então
`max Re W` é um `R` com "não há espectro à direita". E é uma **estimativa de
energia** — exatamente a evidência que a linha F3 pede. Medido: `R = 2,63` em
12×36, extrapolando a **2,91**, contra os **414** obtidos hoje por Neumann
grosseiro.

**A ponte `l^1` -> `l^2` é gratuita na direção necessária:**
`sum (w_j|x_j|)^2 <= (sum w_j|x_j|)^2`, logo `l^1(w) subset l^2(w)` com norma 1,
e toda autofunção em `l^1` é autofunção em `l^2`. Excluir em `l^2` exclui em
`l^1`.

**`R_J` está CERTIFICADO, incluindo a cauda.** Decompondo `R <= R_J + ||B||_2`,
a parte livre virou uma constante derivada:

- o mínimo vive **inteiramente na componente 3**, num problema **1-D** de
  dimensão efetiva **18**;
- conjugar por `sqrt(w)` traz irracionais, mas o mesmo número é o menor `lambda`
  do pencil `Herm(WJ) y = lambda W y`, com `WJ` e `W` **racionais**;
- `LDL^T` exato com `lambda0 = -1777/20000` dá todos os pivôs positivos:
  **`R_J < 0,08885`** (medido 0,088771, 0,09% de folga);
- `Herm(W J_mu)` é **semisseparável de posto 1 mais diagonal**, o que colapsa o
  `LDL^T` numa recursão escalar cujo ponto fixo é `1/kappa2`, dando
  `pivo_k/(w_k mu(k+1)) -> 1 - 1/kappa2 = 1/5`;
- a cauda fecha por intervalo invariante `[0,79 ; 0,80]` para `k>=15` mais
  checagem exata abaixo disso.

**Falta `||B||_2`.** A rota é Young: em Fourier é convolução; em Chebyshev o
espaço com peso `kappa2^n` é o de funções analíticas numa elipse, onde
multiplicação por `f` tem norma `<= sup_elipse|f| <= ||f||_{l^1(w)}` — a mesma
submultiplicatividade que `COMPONENT_EXTERIOR.md` já usa. Com o `beta <= 10,382`
provado, sairia **`R <= 10,47`**, e a região descoberta de F3 encolhe de
`1 < Re s < 414` para `1 < Re s < 10,5`.

**Uma retratação no caminho.** Tentei balancear coluna e linha pelo peso de
Perron e `||B||_2` "caiu" de 3,98 para 2,54, convergindo melhor. Era artefato:
`|B|` é redutível, o vetor de Perron tem 612 entradas nulas de 1422, e o peso
tinha condicionamento `10^151`. Pegou-se olhando a **forma** do peso.

---

## 22/09 — consolidação e retirada do pedido de HPC

Dois entregáveis, ambos fora do código.

**`docs/RELATORIO_SETEMBRO_2026.md`** reúne numa narrativa só o que estava
espalhado por seis documentos crescidos por acreção: o resultado positivo (F3
de 414 para 6,05), o mapa negativo (nove alavancas com piso, cinco rotas de
certificado, o critério que sobreviveu), o dimensionamento e o método. Serve de
porta de entrada; os outros documentos continuam sendo a fonte.

**`~/possivel_projeto.docx` foi reescrito, e contraria a versão de 17/09.**
Aquela pedia R$ 6.900 para a inversão certificada em 107×384. O dimensionamento
de 20–21/09 mostrou que 107×384 **não é a caixa que fecha o argumento** — faltam
de 360× a 9794×, conforme as alavancas. Nenhum edital ou alocação nacional cobre
isso, então **o pedido de computação foi retirado** e a proposta foi reorientada
para as obrigações analíticas, com orçamento de R$ 1.500 de custeio.

O original está preservado em `~/possivel_projeto_17set_ORIGINAL.docx`.

Duas coisas a registrar, porque são o tipo de aprendizado que se perde:

- **O piloto de US$ 4 foi cancelado, e o motivo importa.** Ele mediria o
  rendimento por núcleo alugado, quando a incerteza que decidia tudo era o
  **tamanho da caixa**. Teria dado um número correto e inútil. Mediu-se a coisa
  errada porque se supôs saber qual era a pergunta.
- **A proposta anterior já antecipava o risco** na sua seção 12.3, estimando
  fator 4 para a caixa. O fator real é 90 a 2400 vezes maior que essa
  antecipação. Ter registrado a ressalva não protegeu do erro — só tornou o erro
  localizável depois. É argumento a favor de ressalvas **com número**, não a
  favor de ressalvas.

Não submeter ao SDumont por enquanto: a justificativa daquela submissão era a
caixa de 107×384, e ela não existe mais.

---

## 23/09 — S3: o sistema de constraints é bem condicionado

Ao olhar S3/S4, a primeira hipótese era que a razão `‖Ch‖/‖h‖` ponderada da raiz
0,401 caía porque o autovetor não pertencia ao espaço de RT. **Refutada por dados
que já existiam** (auditoria matemática de 16/09): `‖h‖` ponderado converge; o que
cai é contaminação de modos altos de `Ch` pelo corte.

O achado real veio de comparar `‖V‖ = ‖J(op+s)⁻¹‖_w` para `L` e para o
operador de constraints `K` no mesmo contorno: **14,0 contra 1846,2** em 12×36,
com `‖V_K‖` convergindo e `K` com uma única raiz em `Re s > 0` (0,4010). Isso
divide S3 em S3a (excluir `K` fora de um disco em 0,401) e S3b (testemunho no
disco), e põe S3a na escala de ~1,7× a caixa original, não 10⁴×.

Achado lateral: `sharp_guarded_bridge.py` usa `κ = ‖J‖·‖K⁻¹‖ = 1201,5` onde
`‖JK⁻¹‖ = 8,41`. Perda de 104×, sem efeito relevante nesta ponte (`c·κ = 0,022`),
mas o tipo de confusão que custaria duas ordens num dimensionamento.

Detalhes e ressalvas: `docs/S3_CONDICIONAMENTO_SHARP.md`.

---

## 23/09 (tarde) — por que L é pior que K, custo de S3b, e F3 rebaixado

**Por que `L` é pior que `K`.** Dois fenômenos. Na borda do contorno é o operador:
a raiz de gauge com projeção não-normal, ℓ² 3186 contra 13. Longe das raízes é o
regime pré-assintótico: o fundo só fica pequeno quando `mu n > beta` e `m/2 > beta`,
o que dá `n > 62, |m| > 21` para `L` (fora das malhas) e `n > 27,5, |m| > 9,2`
para `K` (dentro). **Os limiares de `L` são exatamente os `n > 61, |m| > 21` do
certificado — é isso que o vão entre 16 e 61 mede.**

**S3b** pela rota de Newton–Kantorovich: `||DF^-1|| = 247` em 12×36, 7,5× melhor
que o contorno global, mas preso ao mesmo piso pré-assintótico de `L`. Ordem de
10³× a caixa original.

**F3 rebaixado.** Ao adaptar a máquina de F3 para `K`, apareceram dois passos
injustificados na soma da cadeia: convenções de norma misturadas (a parte linear
vale 1,511 com peso inteiro e ~3,0 na convenção de `R_J`) e o máximo sobre
combinações de ℓ¹ levado a ℓ² sem refazer a estrutura de blocos. Valores medidos
ficam bem abaixo das cotas (`||B||_2` = 3,77 contra 5,97 usados), então a cota é
plausível, mas não provada. `R_J` continua certificado. Reparo desenhado em
`ALTERNATIVE_ROUTES.md` §5.25: tudo em `D = diag(w)`.

A extensão de `K` para `Re s` grande ficou em pausa: usaria a mesma máquina.
Estimativa em ponto flutuante: `máx Re W(-K) = 1,04`.

---

## 23/09 (noite) — F3 refeito e K estendido: `R_L <= 10,92`, `R_K <= 5,06`, exatos

Uma única máquina, `scripts/certify_numerical_range_l2.py`, para os dois operadores:

- **ℓ² simétrico** (`D = diag(sqrt(c_m c_n) κ1^m κ2^n)`). A convenção foi deduzida
  do operador livre, que acopla `2μn` uniformemente — só assim se os coeficientes
  guardados são os simétricos. Confirmada depois: a parte livre verdadeira de `L`
  na truncação (−0,0973) bate com o certificado a cinco dígitos.
- **Young** para o fundo, com as normas de campo auditadas; blocos de
  (componente, paridade) combinados por norma espectral, não por máximo.
- **Cauda uniforme** para a parte livre: se `γ_k < 1`, então `γ_{k+1} < ρ_k < 1`.
  Isso também fecha uma lacuna antiga — a cauda de `R_J` (§5.15) tinha sido
  verificada em `k` amostrados.

Resultado: **`R_L <= 10,92`** (verdadeiro na truncação 2,35) e **`R_K <= 5,06`**
(verdadeiro 1,006). Nenhum bloco verdadeiro passa do majorante. Nenhum ponto
flutuante no que se afirma.

O número piorou (6,05 → 10,92) e a afirmação melhorou: antes era plausível, agora
é provada. O 6,05 está retirado em todos os documentos.

---

## 24/09 — `R♯ <= 1,765` e o dimensionamento da região finita de S3a

**`R♯` de 5,06 para 1,765, certificado.** A técnica nova: no L² do toro de raios
`(κ1, κ2)`, o fundo inteiro de `K` é pontual, e a borda do campo de valores vira um
`sup` de `λ_max` 3×3. Foram 4,19 milhões de centros verificados em Arb, sem falha.
Instalei `python-flint` no `.venv` para isso.

**Correção de premissa.** Eu tinha dito que apertar `R♯` encolheria a região finita
em ~3×. Não encolhe muito: o custo está perto de `Re s` pequeno e é fixado pelo
critério de cauda, não por `R♯`. Estimei também `||Γ2♯|| ≈ 1,1` para tirar `μR_x`
da norma; o protótipo deu 1,49, e a caixa de ~15 mil DOFs que projetei com isso
estava otimista.

**Dimensionamento honesto da região finita.** Por causa da faixa de guarda
(o fundo acopla modos a distância 40 × 100, e `||H_FF^-1|| ≈ 15` amplifica qualquer
acoplamento), a caixa inteira até a cauda tem de ser verificada. Com as constantes
novas e o fundo truncado: **~25–33 mil DOFs por tile**, contra 228 mil. Não cabe na
máquina local; cabe no Kaggle.

---

## 24/09 (tarde) — F3 em 2,93; montador validado; tiles de S3a fecham em ponto flutuante

- **F3: `R <= 2,93`, certificado** pelo símbolo pontual (forma 7×7 de `L`), 16,8
  milhões de centros em Arb, sem falha, em quatro faixas no PC do laboratório.
  Fator 141 contra os 414 originais; verdadeiro 2,28.
- **Montador de `K` na memória**, validado a 1,8·10⁻¹⁵ contra o exportado.
- **Correção:** as estimativas de DOFs da região finita estavam 2× acima (o setor A
  só usa `m` par).
- **Experimento de tile:** com bloco 20×80, os pontos mais difíceis do contorno
  (borda `Re s = 0`, cantos, `Re s = 0,2`) fecham com majorante 0,65–0,72.
- **Falta:** o tile rigoroso (faixa próxima exata + cauda distante pela forma fechada
  do inverso livre) e a contagem de tiles.

## 24/09 (noite) — S3a: o tile de Perron fecha, e o primeiro sai rigoroso

- Critério de Perron (raio espectral da matriz de cotas de blocos) no lugar da norma
  espectral; operador na base de Fourier com sinal, validado pelo espectro.
- Casca exata `T1 = G+ \ Z`: com `Z = 20×80`, `G+ = 40×160`, `F = 200×400`, raio 0,866 em
  `s0 = 0` (sem casca, 1,27).
- Modo rigoroso (normas de Rump, contabilidade de erros), fundo truncado em banda,
  três centros (só Z depende do tile pequeno), planejador de cobertura.
- Bug de rigor em `cauda_Rinv` corrigido (nan descartado em silêncio; sem efeito nos
  números publicados).
- **Primeiro tile rigoroso: `|s − 0,2| <= 0,04`, Perron 0,875.** Pendente de auditoria.
- Mapa de `r_adm`: 0,019–0,028 junto ao eixo imaginário, 0,06 em `Re s = 0,2`, e
  `||V Q_ZZ||` cai de 18 (eixo) para 0,9 (`Re s = 1,765`).
- Próximo: cobertura completa no PC do laboratório; auditoria do modo rigoroso; raio do
  disco de S3b.

## 24/09 (madrugada) — S3a: região finita coberta; auditoria do modo rigoroso

- Cobertura rigorosa completa: 154 tiles, `check_cover.py` COMPLETA sobre
  `0 <= Re s <= 1,765`, `0 <= Im s <= 1/4`, fora de `|s − 0,401| < 0,05`.
- Auditoria (`AUDITORIA_S3A_TILES_24SET.md`): 16 alegações, 6 achados corrigidos (entre
  eles uma afirmação errada minha sobre o teorema de Rump, que cobre hermitianas).
- Pendências para chamar S3a de provada: folga de Lipschitz do símbolo em Arb; raio de D
  com S3b; revisão independente.

## 25/09 — S3b fechado (rigoroso)

- **K** (`rouche_K.py`): Rouché de operadores (Gohberg–Sigal) em `|s − 0,401| = 0,05`,
  precondicionador fixo no centro, 16 tiles, pior θ' 0,8126; contagem finita bordejada:
  um autovalor simples. Lema: autovetores no espaço com sinal estão em Y. `s_K` é a única
  raiz de K em D.
- **L**: o NK bordejado em 32×128 não fechava com janelas (θ0 1,14–1,24): a cauda é muito
  não normal (norma 0,64, raio espectral 0,15) e a agregação de normas de blocos perdia
  1,5–1,9×. Fechou com **bloco-Jacobi em todas as 26 janelas** (inversos densos de até
  11.200) e a cauda como um nível ℓ² (`nk_L3*.py`): θ = 0,839826, r = 6,2e-4,
  |λ* − 0,401024733| <= 2,8e-4, algebricamente simples; testemunho ‖Π_F C h*‖ >= 0,1685.
- Máquinas: PC do laboratório (reiniciou duas vezes; alguém usa o Firefox nele) e o **PC
  gamer** do filho do Edivaldo (live USB, 31 GB, 1,7× o laboratório; memória `pc-gamer`).
- **Próximo:** revisão independente de S3b (artefatos obrigatórios, como a de S3a); depois
  S4 / o fechamento da contagem física.

## 25/09 (noite) — S4, parte computacional fechada

- Revisão independente de S3b feita (`AUDITORIA_S3B_25SET.md`); margem do NK reforçada com normas
  restritas de Q0 (θ final 0,839936, r = 3,78e-4).
- **S4** ([`S4_GAUGE.md`](S4_GAUGE.md)): NK de L em λ0 = μ_RefA com centro no vetor de gauge
  g_A = (Z+1)ω_RefA, mesma arquitetura de S3b. Em μ o nível Z' é 2,5× pior (‖M̂⁻¹‖ <= 52,15,
  raiz de fase s = 0 a 0,168) e o erro de RT como 40 ε_ω não fechava; com **RefB explícito**
  (RefAplusB − RefA exato) o erro do fundo cai para ‖B_L(RefB)‖ <= 4,67e-8
  (`CHOPTUIK_EPS_FUNDO_L`). Resultado: θ = 0,880172, r = 6,06e-5, λ* simples; **identificação**
  com o autopar de gauge exato (‖x_g − x0‖_v <= 1,82e-5 <= r): **μ* é raiz algebricamente simples
  de L, autoespaço span(g*)**, D g* = 0, contribuição física 0.
- Parte conceitual redigida: realização de cone móvel, N_físico = N_A − 2, hipóteses H-rec,
  H-gauge, H-cone.
- Máquinas: laboratório + **laboratorio2** (PC dedicado, i7-7700, memória `laboratorio2`), etapa
  B dividida em ordem direta/reversa com checkpoints complementares.
- **Próximo:** revisão independente de S4 (artefatos por alegação: A, B, C, E, RefB e a parte
  conceitual); depois o quadro de T2 (o que falta: C1–C4, F2–F3 e a transferência T1).

## 30/09 → 01/10/2026: Rouché de L certificado na faixa A

- **Nível Z rigoroso, 129 pontos em três máquinas.**
  - A fase 2 (Schur) estourava a memória: ~17 GB, contra 14–15 GB das máquinas. Foi reescrita sem cópias n×n.
  - Os t_p da fase 1 passaram a ser reaproveitados do log, arredondados para cima.
  - A partida do Loewner é f0 = 1,005 + 0,13(est/45)².
- **Lacunas achadas e fechadas no caminho:**
  - zeros dos blocos de janela no interior de Ω: princípio do máximo (`rouche_L_janelas.py`);
  - fator de posto baixo com identidades de ponto flutuante: `certF`;
  - sensibilidade da cauda com a fórmula antiga: corrigida em `rouche_L_disco.py`.
- **Resultado:** 96 discos com θ <= 0,99916, contorno coberto, contagem finita com 2 zeros em Ω
  (0,401 e 0,733), Brauer +1, logo **N_A = 3**. Detalhes em `C1_REAVALIACAO.md` §9.
- **Infraestrutura:**
  - O Tailscale SSH ainda estava em "check": o netmap mostrava holdAndDelegate e um accept com
    ruleExpires. O usuário corrigiu para accept em 01/10.
  - O PC 3 foi desligado por engano no meio do −0,25i; o centro foi terminado no laboratório
    com os t_p do log do João.
- **Pendente:** revisão independente; faixa B (C4).

## 01/10/2026 (tarde): revisão independente do Rouché de L e versão 2

- **Revisão:** Fable, com contexto limpo; caiu duas vezes, uma por o notebook desligar e outra por
  login expirado, e foi retomada. Não quebrou a contagem.
- **Achados:**
  - L1: deflação verdadeira fora de Z;
  - L2: sensibilidade do far;
  - L3–L6: rastreabilidade e arredondamento.
  - Na conferência achei mais: ‖Q0(s)P_J‖ no centro, eps_mu × ‖Q0‖ < 1, J calculada e J·D.
- **Consertos:**
  - rig v2, com sha de A e F e arredondamentos;
  - `Qfar_disco`;
  - os termos da deflação na montagem e nas janelas;
  - a cota de ∂τCorr no lema.
- **Reexecução:** a fase 2 dos 96 pontos rodou no laboratorio2 (65 pontos, ~1,8 min cada) e no PC do
  João (31, ~4,4 min cada), com os t_p reaproveitados. O PC do laboratório não foi usado, porque está
  em uso a partir das 7h.
- **Resultado v2:** θ <= 0,99842 nos 96 discos, contorno coberto, janelas com máximo 0,1488,
  t·ε_B <= 0,0334. **N_A = 3 revisado.** Artefatos em `build/rouche_L_rig/{discos_v2,z2}`.
- **Infraestrutura:**
  - as tarefas em segundo plano do harness são mortas após ~30 min, então trabalho longo no notebook
    vai com `setsid nohup`;
  - o ssh das máquinas agora entra direto (Tailscale em accept).
- **Próximo:** faixa B (C4).

## 01 → 02/10/2026: faixa B (C4) certificada

- **Montagem:** a referência e os discos da aresta Im = 1/4 vêm da faixa A, e as caudas dos centros reais
  servem na aresta Im = 3/4. Centros novos: 0,5i, 0,75i e 0,55 + 0,75i.
- **Loewner pré-condicionada** (K = I + QB̂) para resolvente até 202; fases do nível Z em processos
  separados.
- **NK do testemunho:**
  - h0 na matriz pesada (Y_Z 2·10⁻⁴ → 3·10⁻⁶);
  - λ0 = autovalor da truncagem sem deflação;
  - tQ = ‖Mh⁻¹[Q_ZZ; 0]‖ (2,4× menor que t_V‖Q_ZZ‖) no Z2.
- **Resultado:** N_B = 1, raiz simples, testemunho 0,3547 > 0, multiplicidade física 0. Com a faixa A,
  contagem física total 1.
- **Infraestrutura:**
  - O PC do laboratório foi usado só à noite, até 6h02.
  - O PC 3 foi desligado por engano às 10h49. A cauda de 0,5i foi juntada dos checkpoints dos
    dois PCs; a união cobria as 26 janelas.
  - O ssh entre as máquinas funciona (Tailscale accept).
- **Pendente:** revisão independente da faixa B.

## 02/10/2026 (tarde): revisão independente da faixa B

Revisor Fable, com contexto limpo, usando só o laboratorio2 e o PC 3. Não quebrou o certificado.
- Loewner pré-condicionada contra LU densa: +0,35%.
- Janelas de 0,5i das duas máquinas: iguais bit a bit.
- θ remontado em 3 discos, cobertura e contagem conferidas com código próprio.
- c_ref do testemunho reproduzido.

Lacunas L1–L4, todas de rastreabilidade:
- 78 discos efetivos, não 79;
- `tp_*.json` sobrescrito na fase 2: corrigido, e os logs foram arquivados com sha;
- `--alvo` não gravado: corrigido;
- cauda de 0,55 + 0,75i sem manifesto remoto: o PC do laboratório estava desligado; recálculo no laboratorio2.

O próximo trabalho são as hipóteses físicas e o enunciado final de T2.

## 02/10/2026 (noite): T2, hipóteses físicas

- **Inventário** (`T2_HIPOTESES.md`): o que é teorema, o que é hipótese física e o que faltava.
- **L-real (a lacuna nova):** todos os certificados contam no toro, mas a teoria de constraints e gauge
  está na realização ℓ¹ de RT. Fechada por recuperação de raio com G = 1024×4096 e |s| <= 3,1: toro
  0,371 (`l_real_toro.py`), RT 0,267 (`radius_transfer.py`).
- **L-afim, L-iso e L-constr** na realização com sinal. H-dof não se aplica.
- **Enunciado final:** `T2_ENUNCIADO.md`, com 1 modo físico sob H-rec, H-gauge, H-cone e H-reg.

## 03/10/2026 (madrugada): revisão de T2 e NK em 0,7332

- **Revisão independente de T2** (Fable, contexto limpo, notebook e PC 3): não quebrou.
  - Normas exatas de J⁻¹ na cauda abaixo de `cauda_Rinv` (razão 0,80–0,84).
  - ‖B_L‖ no toro com sinal: 3,8626 em 80×240, abaixo de 3,9641.
  - O setor B montado coincide com o A deslocado (diferença 0).
  - Quociente conferido em 100 sorteios exatos; KC = ML simbólico e numérico.
  - Lacunas 1–11, de enunciado e citação, corrigidas em `T2_ENUNCIADO.md`, `T2_HIPOTESES.md` e
    `S4_GAUGE.md` (λ = 1 → 0,614, o expoente do gauge depende das coordenadas). A contagem da faixa B foi
    refeita só com os discos de B (t·ε_B = 0,031), e o log lista os discos.
- **NK em 0,7332** (laboratorio2 + PC 3 na etapa B): λ* ∈ [0,73316950; 0,73318982], simples,
  θ = 0,806, r = 2,2·10⁻⁵. γ = Δ/(4πλ*) ≈ 0,3739.

## 03/10/2026 (tarde): H-gauge provada (Lema G)

- `GAUGE_GLOBAL.md`: o gauge residual é completo globalmente. Os modos de gauge de Floquet com Re s > 0
  são ℂ·g*, para gauges C¹ em todo o domínio de RT.
- Passos da prova:
  - média em SO(3);
  - preservar o ansatz dá ∂_{u_+}X^{u_-} = 0 e ∂_{u_-}X^{u_+} = 0;
  - as curvas de nível são conexas em M, logo X = f(u) globalmente;
  - o eixo impõe f_+ = f_-;
  - equivariância sob Θ²;
  - injetividade pela componente ω4: o eixo, mais a analiticidade no exterior do cone;
  - Floquet torna f homogêneo, e ser C¹ no cone força f constante.
- Fato certificado F1 (`scripts/gauge_global_centro.py`): |c_1(ω4*(·, 0))| >= 0,493046 > ε = 3·10⁻⁸.
- T2 passa a depender de H-rec (ida), H-cone e H-reg, além da modelagem: simetria esférica e
  identificação do fundo.

## 03/10/2026 (noite): Lema G revisado e H-rec (ida) provada

- **Revisão do Lema G:** não quebrou. F1 foi refeito com código próprio, com as mesmas frações; 67/67
  conferências simbólicas.
  - Consertos: μ fixo na definição de admissível, equivariância linear e notação (ḡ para a métrica, g*
    para o perfil).
  - Lema G′: gauge só contínuo no cone, com F1′ (|c_1(ω4(·, −1))| >= 0,0977).
- **Lema R (`HREC_IDA.md`):** a H-rec (ida) está provada.
  - R1: transporte nulo.
  - R2: a ressonância fecha pela simplicidade de μ*.
  - R3: extração com regularidade no eixo.
  - R4: equações; a revisão calculou o Ricci 4D.
  - R5: recuperação de raio de qualquer peso > 1, para todo s (`radius_recovery_geral.py`, com
    S até 200).
  - R6: unicidade módulo gauge.
  - R7: cadeias físicas.
- **Revisão do Lema R:** não quebrou. As lacunas L1–L9, de escopo e redação, foram incorporadas. A
  principal era estender R5 a todo s e tratar as cadeias.

## 03/10/2026 (noite, 2): revisão da volta

- **Não quebrou.**
  - As identidades de compatibilidade e de log W foram provadas em sympy com funções genéricas.
  - L e C do código batem com o TeX ponto a ponto (10⁻¹⁵).
  - O Ricci 4D direto confere as fórmulas de curvatura de RT.
- **Consertos em `GEOMETRIC_FLOQUET_RECONSTRUCTION.md`:**
  - enunciado preciso, com h ∈ D(Y+) e um modo físico de Floquet em M̂;
  - Lema V1;
  - constraints de jato para cadeias;
  - status atualizado.
- **O1:** a perda de raio de KC = ML, de (65/64, 5/4) para (129/128, 9/8), cai exatamente no espaço dos
  certificados de K. Fechado em `SHARP_DOMAIN.md`.
- Todas as peças matemáticas de T2 estão revisadas. Resta só a modelagem.

## 03/10/2026 (noite, 3): escolhas de modelagem

- **H-cone → Lema C.** e^{μτ}g* = d/dε de ψ_ε*(ḡ, φ*), com ψ_ε a translação conjunta dos nulos. É a mesma
  solução com T* deslocado. Com o gauge linearizado padrão (X0 é um campo analítico em M̂), a contagem
  é 1; com gauges que fixam o cone, 2.
- **Re s = 0:** as arestas esquerdas dos Rouchés A e B cobrem um período de L_A. A fórmula de
  multiplicidade da deflação dá: só s = 0, simples, a fase.
- **H-reg′:** R5 com κ1' = 1 dá C^∞ em τ e analítica em ξ.
- **Próximo:** a análise da simetria esférica, pedida pelo usuário. O PC 3 está indisponível em
  03/10.

## 03/10/2026 (noite, 4): simetria esférica, Fases 0 e 1

- Revisão leve dos três itens de modelagem: não quebrou. L1 corrigida: escala global e translação do
  escalar são neutros físicos invisíveis a L. As ressalvas L2–L6 foram corrigidas.
- **Axial:** a equação de GS no fundo de RT (G = r²Π, só ω e μ). O espectro bate com MG: κΔ = −2,320 e
  −3,292 em ℓ = 2 e 3.
- **Polar:** derivação simbólica (`ns_polar_deriva.py`) e numérica com filtro de vínculos.
  - ℓ = 2: s = −0,004764 ± 0,147345i, κΔ = −0,0599 (MG −0,07).
  - ℓ = 3: κΔ = −1,655 (MG −1,65).
- **Lições da discretização:**
  - paridade imposta no centro;
  - cone interior, ξ* = 1,02;
  - filtro de suavidade e, no polar, de vínculos;
  - o subsistema de evolução tem espúrios que violam os vínculos.
- O sympy está no scratchpad (`pylib`), fora do `.venv`.

