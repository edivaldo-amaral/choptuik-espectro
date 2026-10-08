# Ledger das obrigações de prova

Nenhuma linha marcada “aberta” pode ser substituída por uma diagonalização em
ponto flutuante.

## Auditoria de 16/09/2026 (tarde) do trabalho feito entre 14/09 e 16/09

Os parágrafos datados de 14/09 à noite a 16/09 manhã logo abaixo foram escritos
sem revisão independente. Dois auditores recalcularam e re-derivaram tudo de
forma independente: [`AUDITORIA_MATEMATICA_16SET.md`](AUDITORIA_MATEMATICA_16SET.md)
cobre cauda, Schur e coroa; [`AUDITORIA_FISICA_16SET.md`](AUDITORIA_FISICA_16SET.md)
cobre constraints, gauge e setor B. **Nenhuma alegação foi refutada.** O que mudou:

- **Confirmado com recomputação:** F1/recuperação de raio (defeitos exatos);
  Schur S+T <0,997676 (0,99759 com κ<31,43); coroa com 51 tiles (defeito
  <0,4885, κ≤31,43); discos e ponte sharp (δ 0,915031 e 0,989223); constantes
  18/(μ+a) e ‖B‖<23 contra RT (22), (33)–(39), (60). As identidades
  KC=ML−G(h,e(w)), a ação gauge, `s_n=mu(1−n)` e a reconstrução geométrica
  valem universalmente (`scripts/independent_rt_symbolic.py`, com mutações).
- **Sobrevendido:** (i) a faixa U "reduzida a um Schur finito de 228366 DOFs"
  é inviável nesta máquina (≈417 GB densos, amplificação ≈7,6·10⁴); (ii) C1
  e C2 soavam fechadas, mas P2 falha e C1-G está aberta; (iii) **C4:** o
  "winding 0 no setor B" é artefato do 4×12 (ver linha C4).
- **Fato novo:** de 6×18 a 12×36, `L_A` tem uma raiz convergente na faixa B,
  `0,041421+0,500008i` (real no setor B, `s_B≈0,0414`, `lambda≈0,151`), fora
  do retângulo certificado. Ela tem a assinatura de violação de constraints
  (‖Ch‖/‖h‖₂≈0,50 estável, colinear à raiz sharp a 7,5e-8). Confirmada pelo
  orquestrador com `scipy` direto sobre `rt-A-*.dat`.
- **Revisto:** a estagnação `|s−mu_RT|≈2,21e-6` é truncamento do gerador
  gauge exato (parcela do fundo ≈1,1e-10), não evidência contra `s=mu`.
- **Testes vazios:** `test_neutral_phase_projection_does_not_preserve_positive_kernel`
  não aplica operador; `test_geometric_residual_combinations_are_invertible`
  inverte uma matriz digitada à mão.
- **Todas as rotas finito→infinito** consomem o mesmo majorante exterior
  uniforme, que só contrai em caixas de ~5,8·10⁵ DOFs. A viabilidade de um
  majorante por colunas está em teste (rodada 2, regra de decisão pré-registrada
  em `COLUMNWISE_EXTERIOR.md`).

## Rodada 2 (17/09): majorante exterior por colunas — rota de blocos CONTINUA

[`COLUMNWISE_EXTERIOR.md`](COLUMNWISE_EXTERIOR.md), Marcos 1–7. Com a regra
pré-registrada (N_G<=384) aplicada sem mudança, a caixa G=107×384 dá
d*(G)<1 uniformemente em |s|<=5/4, com a bola RT:

| sistema | d*(G) | cota dominante | DOFs de G |
|---|---:|---|---:|
| A | 0,99820 | Fourier (Fase 1, m>=107) | 143424 |
| sharp | 0,99903 | radial (Fase 1, n>=606) | 61632 |

- **Faixa n∈[384,606), m<107, conferida coluna a coluna** em racionais
  exatos (núcleo C+GMP validado contra duas implementações Python): 57276
  colunas no setor A e 24642 no sharp, todas <1. Máximos 0,50493 e 0,46300.
- **Margens novas:** Schur GUARDED com defeito 0,998763 (margem 1,24e-3,
  U=143091 DOFs, antes ~5,8·10⁵); ponte sharp com 0,999033 (margem 9,7e-4,
  U=61038, antes 228366).
- **Hipóteses ainda em papel:** as invariâncias de β para j>=101 e m>=41
  (com cauda explícita) e as reduções por simetria de imagens, conferidas
  só em casos-teste. Nada disso está em Lean.
- **Critério 2 falhou como medida de folga** (d* chega a 17,3× o float),
  porque a folga está nas cotas analíticas além da faixa, não por coluna.
- **Testemunho S3** (seção (d) da auditoria matemática): a queda na norma com
  pesos é contaminação de modos altos, não C h→0. O testemunho precisa de
  norma com perda de raio e das hipóteses H1 (multiplicidade 1 do operador
  infinito), H2 (enclausuramento do autovetor verdadeiro), H3 (C limitado
  com perda de raio) e H4 (S4).
- **Próximo gargalo:** a parametriz da faixa U, ainda densa demais para 7 GB,
  e a redução de M_G conferindo coluna a coluna em m∈[48,107).

## Rodada 3 (17/09): revisão física do d*(G) e a hipótese H3

Revisão em [`AUDITORIA_FISICA_16SET.md`](AUDITORIA_FISICA_16SET.md) §11, com
código que não importa `columnwise_*`:

- **Confirmado:**
  - as invariâncias de β em j>=101 e m>=41, derivadas em papel e checadas
    em j=137, 211 e m=53, 87, 88, com diferenças sob a cauda R;
  - U(e) em 12 colunas, incluindo as duas argmax, com valores recomputados
    sempre ≤ os deles (+1e-6 a +1e-5);
  - a montagem {m>=107} ∪ faixa ∪ {n>=606}, sem buracos.
- **Não verificado de forma independente:** as cotas analíticas da Fase 1,
  d_R(606) e d_F(107), que FIXAM as margens de 10⁻³ (custo estimado de meio
  dia), e os escalares de Schur.
- **H3 fechada** ([`CONSTRAINT_WITNESS_NORMS.md`](CONSTRAINT_WITNESS_NORMS.md)):
  - C(s)=D(s)+A_ω, verificada exatamente para campos gerais;
  - entrada nos pesos fortes sem η, saída em (129/128,9/8), com a bola RT e
    uniforme em |s|<=5/4, que contém o retângulo B: c_op<=40,506 sem
    projeção, <=22,179 com a saída projetada em m<12;
  - ε_apl<3,04e-7.
- **Testemunho S3 inalcançável com a máquina atual** (diagnóstico em float):
  - c_ref/c_op≈1,24e-3 na raiz 0,401 e ≈5,25e-3 na raiz B;
  - o autovetor 12×36 difere do 10×30 em 3,8e-2 (A) e 5,6e-2 (B) nessa
    norma, 10–30× acima do ε_v tolerável;
  - com margens de Schur ~10⁻³, o enclausuramento amplificaria resíduos
    em ~10⁴.
  - Gargalo: H2, o enclausuramento do autopar infinito. H1 e H4 continuam
    abertas.

Retomada de 16/09: SHARP_INTERMEDIATE_SHELL.md agora propaga o resto de
Neumann por AMBOS os acoplamentos até o Schur intermediário. Para erro
<10^-12, os majorantes racionais pedem 3786 termos em 148×840, ou 438
em 160×960 com a ponte melhorada. A parametriz da faixa não foi construída;
S3 continua aberta. GLOBAL_NULL_GAUGE_TRANSPORT.md resolve o transporte
de fixação do bloco radial global numa classe analítica explícita, com
perda de domínio e condição de compatibilidade ressonante no cone.
Não identifica essa classe com o domínio físico/RT forte, e S4 continua
aberta. Os novos argumentos de operadores não foram formalizados em Lean.

Rodada seguinte: SHARP_INTERMEDIATE_SHELL.md valida a faixa 8×24->12×36
e o bloco F mais cauda infinita fora de 160×960 no disco físico candidato.
A exclusão sharp foi reduzida a um Schur finito não validado de 228366
DOFs; nenhum desses DOFs foi descartado. GEOMETRIC_FLOQUET_RECONSTRUCTION.md
dá equivalência linearizada global NO ANSATZ para Re s>0 e para jets,
incluindo centro e cone. Continua aberta a fixação de gauge global que
leve uma perturbação física geral ao ansatz com as condições desejadas.

Nova revisão: RADIUS_RECOVERY.md prova, sob os bounds RT, coincidência
dos kernels/cadeias selecionados nos pesos fortes e menores para |s|<=5/4.
Isso coloca o gerador gauge s=mu no domínio forte e dá Fredholm de índice
zero na realização RT (não a equivalência geométrica). A aritmética da
cauda usa defeitos <0.230 e <0.786, sem inversão de uma caixa gigante.
SHARP_EXTERIOR_BOUNDS.md inclui o erro do fundo no disco sharp finito:
defeito <0.263453. A cauda sharp isolada contrai em 256×1536, mas falta
a faixa intermediária e seus acoplamentos. CONSTRAINT_ROOT_QUOTIENT.md
separa dim E, rank das constraints e a linha gauge; essas dimensões
infinitas ainda não foram numericamente certificadas.

Continuação de 16/09: GUARDED_INFINITE_COUPLING.md controla o subsistema
S+T, S=(6×18)-(4×12), T fora de 176×944, com defeito de Schur <0.997676
uniformemente no contorno, condicionado à bola RT. A faixa intermediária
permanece no problema restante. SHARP_LOCAL_EXCLUSION.md dá exclusão
racional FINITA em |s-0.733|<=1/32, nos pesos menores, e explica o
transporte das constraints nas cadeias. GAUGE_RT_ACTION.md identifica
h=(1+mu^-1 partial_tau+xi partial_xi)omega como gerador nulo-translacional
em s=mu, com constraints nulas no fundo exato e perda de raio controlada.
Ainda não se subtrai sua multiplicidade da contagem nos pesos fortes.

Rodada de 16/09: `SHARP_SPECTRAL_LOCATOR.md` registra o novo exportador
sharp e a aproximação entre as raízes selecionada/sharp próximas de 0.401
(distância numérica 3.2e-7 em 12×36). Não é uma exclusão certificada.
`DAMPED_CROWN_PARAMETRIX.md` descreve a parametriz amortecida e a auditoria
racional de tiles da coroa finita, incluindo a mudança de mu e o erro do
fundo quando essa opção é usada. O bloco exterior infinito continua separado.
Cobertura concluída: 51 tiles, defeito <0.503 e inversa da coroa
uniformemente <36.227, incluindo o erro do fundo sob a hipótese RT.
Não implica contração acoplada nem winding zero da coroa.

Auditoria de 14/09/2026: veja [`ANALYTIC_AUDIT.md`](ANALYTIC_AUDIT.md).
A identificação da resolvente comprimida com a inversa finita foi refutada
por um contraexemplo verificado em Lean. P2 exige reformulação ou hipóteses
adicionais sobre os espaços; não é consequência da analiticidade do fundo.

Uma rota alternativa foi derivada em [`ALGEBRAIC_TAIL.md`](ALGEBRAIC_TAIL.md):
na realização por fechamento com pesos RT, a inversa livre deslocada dá
cauda algébrica uniforme e um ponto-base seguro em `s=828`, sob os bounds
publicados de RT. Os valores racionais estão verificados em Lean. A
equivalência desse domínio com o físico, a conversão de normas nos tiles e
C1-G continuam abertas; por isso as linhas históricas abaixo não devem ser
lidas como fechamento incondicional de F1/F3/C1.

Estimativa por blocos: [`FREQUENCY_BLOCKS.md`](FREQUENCY_BLOCKS.md) separa
núcleo, coroa finita e exterior infinito. A parametriz acoplada 6×18 passou
no centro real s=1; caudas direcionais do fundo foram calculadas exatamente,
com erro publicado separado. Agora existem majorantes uniformes de
coroa/exterior e exterior/exterior, mas a contração acoplada ainda falta.
`STRUCTURED_EXTERIOR.md` preserva cancelamentos da inversa livre;
`COMPONENT_EXTERIOR.md` dá beta<10.382 numa norma reponderada, com TT<4.882
em 64×192. Não combinar esses bounds com blocos finitos em outra norma.
`SHARP_LINEARIZED_IDENTITY.md` constrói a identidade diferencial completa
K(s)C(s)=M(s)L(s), com a correção do fundo aproximado explicitada.

Nova etapa: `UNIFORM_COUPLED_CONTRACTION.md` explicita o complemento de
Schur, o raio estrito por tile e a homotopia que preserva o bound.
`SHARP_DOMAIN.md` prova compatibilidade dos domínios com perda de raio;
o espectro sharp precisa ser controlado nesses pesos menores.
`GAUGE_DOMAIN_CLASSIFICATION.md` prova a transversalidade da fase e
classifica reparametrizações nulas analíticas que preservam o centro,
distinguindo cone fixo de cone móvel. Não prova completude do gauge físico.

## 01/10/2026: contagem na faixa A certificada pelo Rouché de L (revisada em 01/10; versão 2)

Rouché de operadores (Gohberg–Sigal, forma de homotopia) no contorno ∂([0, 3] × [−1/4, 1/4]):
- referência finita no nível Z, com deflação de Brauer conjunta biortogonal da fase e de μ;
- blocos de janela, bloco-Jacobi da cauda, como em S3b e S4.

Peças:
- 96 discos com Perron de 3 níveis em racionais, θ <= 0,99916, cobrindo o contorno;
- blocos de janela sem zeros em Ω, pelo princípio do máximo;
- certificado do fator de posto baixo (`certF`);
- contagem finita por Schur, com t·ε_B <= 0,033.

Resultado: **N(L; 0 < Re s < 3, |Im s| < 1/4) = 3** (μ*, s_K = 0,401 e 0,733), com multiplicidade.
Com F3 (R <= 2,93), vale para toda a faixa A, sem a fita 0 < Re s < 1/8. Detalhes, números e
artefatos estão em [`C1_REAVALIACAO.md`](C1_REAVALIACAO.md) §9. Isso substitui as rotas de C1, C1-G,
C2, C3 e T1 na faixa A. A faixa B (C4) continua fora.

**Revisão independente (01/10, `.codex-runs/2026-10-01-revisao-c1/`):** não quebrou a contagem. Achou duas
lacunas de derivação (L1 deflação verdadeira fora de Z; L2 sensibilidade do far) e quatro de rastreabilidade
e arredondamento (L3–L6). Todas foram consertadas e os 96 discos refeitos (versão 2, θ <= 0,99842, contorno
coberto, t·ε_B <= 0,0334). Detalhes em `C1_REAVALIACAO.md` §10.

## 02/10/2026: faixa B (C4) certificada e revisada

Rouché de L em ∂([0, 3] × [1/4, 3/4]) com a mesma referência deflacionada:
- 174 discos (θ <= 0,99907), Loewner pré-condicionada até t = 202;
- janelas sem zeros (máx 0,152);
- contagem finita com 1 zero (t·ε_B <= 0,033).

Logo **N(L; faixa B, Re s > 0) = 1**. NK bordejado em λ0 = 0,0414194105 + 0,5i: raiz simples (r = 7,8·10⁻⁵).
Testemunho: ‖Π_F C(λ*)h*‖ >= 0,3547, logo **multiplicidade física na faixa B = 0**. Com a faixa A, a contagem
física total em Re s > 0 é **1** (só a raiz 0,733), condicionada às hipóteses de S4 e da redução do setor B.
Detalhes em [`C1_REAVALIACAO.md`](C1_REAVALIACAO.md) §11. Revisão independente (§11.1): não quebrou o certificado; L1–L4 de rastreabilidade, todas fechadas (a cauda de 0,55 + 0,75i, refeita no laboratorio2, coincide a 5·10⁻¹⁶).

| ID | Obrigação | Estado | Evidência/saída exigida |
|---|---|---|---|
| E1 | Existência do fundo DSS no domínio RT | publicada | Teorema RT + dados dyádicos |
| C1-G | Comparação entre resolvente comprimida e inversa finita | **01/10: na faixa A, substituída pelo Rouché de L, certificado ([`C1_REAVALIACAO.md`](C1_REAVALIACAO.md) §9–10: N_A = 3); revisado em 01/10 (versão 2).** aberta, identidade de resíduo derivada | bound para `Π R0 E_N`, incluindo o resíduo fora do truncamento; ver `ANALYTIC_AUDIT.md` |
| E2 | Auditoria integral da CAP de existência | parcial | faltam dois arquivos `EXTERIOR_ESTIMATE` não depositados |
| S1 | Derivação de `L_A(s)=L_A(0)+sI` | derivada, implementação exploratória | equivalência geométrica formal ainda aberta |
| S2 | Pencil do setor oposto `L_B` | **reduzido a A**, sem implementação nativa | conjugação verificada em [`SECTOR_B.md`](SECTOR_B.md); `Sigma_B = Sigma_A - i/2` por Gohberg--Sigal. Nada precisa ser complexificado |
| S3 | Propagação/invertibilidade das constraints para `s` complexo | identidade e domínios RT com perda de raio demonstrados; parte espectral infinita aberta | **23/09:** `K` é ~130× mais bem condicionado que `L` no mesmo contorno (`||V_K||`≈14 contra 1846 em 12×36, convergindo) e tem uma só raiz em Re>0, 0,4010. S3 divide-se em S3a (excluir `K` fora de um disco em 0,401 ⇒ todo modo de `L` fora dele satisfaz as constraints) e S3b (testemunho no disco). S3a vive na escala de ~1,7× a caixa original, não 10⁴×. Ver [`S3_CONDICIONAMENTO_SHARP.md`](S3_CONDICIONAMENTO_SHARP.md). **24/09:** `R♯ <= 1,765` certificado pelo símbolo pontual (`K` sem espectro em `Re s > 1,765`); a região finita de S3a pede uma caixa verificada de ~25–33 mil DOFs por tile (§9). SHARP_DOMAIN e SHARP_LINEARIZED_IDENTITY dão KC=ML nos domínios; SHARP_EXTERIOR_BOUNDS inclui o fundo no disco finito. Falta controlar kernel sharp infinito e acoplamentos, não apenas a fronteira. **24/09 (noite):** tile de S3a por **critério de Perron** (raio espectral da matriz de cotas de blocos) com casca exata: `Z = 20×80` (inverso denso de 2280 DOFs), casca `40×160`, `F = 200×400` modo a modo, `far` pelo símbolo. Raio 0,57–0,87 na região; primeiro tile rigoroso (normas de Rump, erros a priori) certificado em `|s − 0,2| <= 0,04`, **pendente de auditoria**; **cobertura completa em 24/09: 154 tiles certificados cobrem toda a região finita de S3a fora do disco |s − 0,401| < 0,05** (verificador independente em racionais). Disco de S3b fixado em |s − 0,401| < 1/8 (`S3B_DISCO.md`). Revisão independente feita (agente de outro modelo, §6 da auditoria): achou duas lacunas (caminho de μR_x de far até G+; μ verdadeiro de RT a 3,9e-11 do diádico), ambas fechadas; 154 tiles continuam certificados (pior θ 0,99612) e a cobertura completa. **S3a: sem pendências conhecidas.** (`AUDITORIA_S3A_TILES_24SET.md`). **25/09 — S3b FECHADO (rigoroso; revisão independente feita, `AUDITORIA_S3B_25SET.md`: sem erro que derrube, L1–L5 corrigidas, margem do NK reforçada — θ = 0,839936, r = 3,78e-4, |λ* − 0,401024733| <= 1,71e-4, testemunho >= 0,16976):** K tem exatamente uma raiz em |s − 0,401| < 0,05, simples (Rouché de operadores, 16 tiles na circunferência, pior θ' 0,8126, mais a contagem finita bordejada), logo `s_K` é a única raiz de K em D; certificado NK bordejado de L (bloco-Jacobi em 26 janelas densas, cauda em ℓ², θ = 0,839826, r = 6,2e-4) dá λ* com \|λ* − 0,401024733\| <= 2,8e-4, algebricamente simples; testemunho ‖Π_F C(λ*)h*‖ >= 0,16849 > 0. **Exatamente uma raiz de L em D viola as constraints: `s_K`, simples.** Ver [`S3B_ROTA.md`](S3B_ROTA.md) §8 e §12. Ver §10.4–10.7 de [`S3_CONDICIONAMENTO_SHARP.md`](S3_CONDICIONAMENTO_SHARP.md). **Auditoria 16/09:** o exportador sharp calcula `Π♯K(0)Π♯` de RT, reproduzido a 1e-15. Nas raízes 0,401024 (A) e 0,041421+0,500008i (B), `C(s)h` é colinear ao autovetor sharp (1−\|cos\|≈7e-8), com ‖Ch‖/‖h‖₂≈0,32 e ≈0,50 estáveis; na norma (129/128,9/8) a razão da raiz 0,401 ainda CAI (0,172→0,056), em análise. Um testemunho exige enclosure certificado do autovetor E multiplicidade algébrica 1 |
| S4 | Classificação de gauge e equivalência ao quociente físico | **parte computacional certificada (25/09)**; equivalência global como hipóteses explícitas (H-rec, H-cone); H-gauge provada em 03/10 (Lema G, [`GAUGE_GLOBAL.md`](GAUGE_GLOBAL.md)) | **25/09 ([`S4_GAUGE.md`](S4_GAUGE.md)): μ* é raiz ALGEBRICAMENTE SIMPLES de L, com autoespaço span(g*), g* = (Z+1)ω*.** NK de L em λ0 = μ_RefA (mesma arquitetura de S3b, centro no vetor de gauge, erro do fundo por RefB explícito ‖B_L(RefB)‖ <= 4,67e-8): θ = 0,880172 (Perron em racionais), Z2 = 183,7, Y = 6,6e-6, r = 6,06e-5; identificação (`nk_L3_E.py`): ‖x_g − x0‖_v <= 1,82e-5 <= r, θ + 2 Z2 ρ = 0,902 < 1. Com D g* = 0, a contribuição física de μ é 0 (CONSTRAINT_ROOT_QUOTIENT) e N_físico = N_A − 2 na faixa A (−1/4 <= Im s <= 1/4; Σ_A é i-periódico), realização de cone móvel. **Revisão independente 26/09 ([`AUDITORIA_S4_26SET.md`](AUDITORIA_S4_26SET.md)):** nenhum erro; faixa fundamental e hipóteses (iv), (v), H-reg incorporadas; lema ‖(Z*+1)Corr‖_L <= 1,1e-75 escrito (`corr_gauge_bound.py`, constantes publicadas de RT). Antes: GAUGE_RT_ACTION e RADIUS_RECOVERY dão o gerador não nulo em s=mu no domínio forte. A proximidade da raiz finita não certifica sua identificação. CONSTRAINT_ROOT_QUOTIENT trata cadeias e a linha gauge; faltam fixação global, reconstrução física e escolha consistente do cone. **Auditoria 16/09:** `DΩ(w)[e^{μτ}(Z+1)w]=e^{μτ}(Z+1)Ω(w)` verificada universalmente; a distância 2,21e-6 da raiz a `mu` decompõe-se exatamente em truncamento (+2,2126e-6) e defeito do fundo (+1,1e-10), com \|cos\|((Z+1)ω, autovetor)→1. Isso fortalece a identificação numérica; a subtração continua dependendo de S4. Nenhum modo gauge nulo analítico cai na faixa fundamental de B |
| F1 | Família Fredholm analítica de índice zero | demonstrada em papel na realização RT, **recomputada de forma independente em 16/09**; ponte física aberta | Realização D:=J_μ^-1(Y), com Y o ℓ¹ ponderado de RT (65/64,5/4), complexificado. Q=J_μ^-1 compacto (‖Q(I−G)‖≤t(G)→0), B limitado por RT (33), (34), (36), H(s)=I+(B+s)Q Fredholm de índice 0, invertível para s≥414. Hipóteses: 1/6≤μ≤17/100 e a bola RT (39b), (60). RADIUS_RECOVERY: em \|s\|≤5/4, kernels e cadeias coincidem nos pesos (65/64,5/4) e (129/128,9/8). Não formalizado em Lean; não identifica D com o domínio geométrico |
| F2 | Coercividade de alta frequência | aberta | matrizes corretoras e certificado LDL* intervalar |
| F3 | Ausência de espectro para `Re(s)>=R` | **rota aberta com uma peça certificada** (21/09/2026) | **24/09: `R <= 2,93` certificado pelo símbolo pontual no L² do toro** (16,8 milhões de centros em Arb, sem falha; `ALTERNATIVE_ROUTES.md` §5.26). Antes: **23/09: `R <= 10,92` certificado em aritmética exata** (ℓ² simétrico, Young, cauda uniforme; `certify_numerical_range_l2.py`). A versão 6,05 de 22/09 misturava convenções e está retirada (`ALTERNATIVE_ROUTES.md` §5.25). estimativa de energia quantitativa — e o **campo de valores É uma**: `spec(T) subset W(T)`, logo a borda direita de `W` dá `R` direto. Realização `l^2` com pesos; a ponte é gratuita na direção necessária, porque `sum (w_j|x_j|)^2 <= (sum w_j|x_j|)^2` dá `l^1(w) subset l^2(w)` e toda autofunção em `l^1` é autofunção em `l^2`. Decompondo `R <= R_J + ||B||_2`: **`R_J < 0,08885` está CERTIFICADO** em aritmética racional, incluindo a cauda infinita — redução a um problema 1-D da componente 3, `LDL^T` exato do pencil `Herm(WJ)y=λWy`, e a cauda por recursão escalar com ponto fixo `1/kappa2` e intervalo invariante. Falta `||B||_2`, cuja rota é Young/elipse sobre a submultiplicatividade já usada em `COMPONENT_EXTERIOR`; com o `beta<=10,382` provado daria `R <= 10,47` contra os **414** atuais. Ver `ALTERNATIVE_ROUTES.md` §§5.9–5.16 e `scripts/certify_free_numerical_range.py` |
| C1 | Aproximação compacta `||C-C_N||<=epsilon_N` | **01/10: na faixa A, substituída pelo Rouché de L, certificado ([`C1_REAVALIACAO.md`](C1_REAVALIACAO.md) §9–10: N_A = 3); revisado em 01/10 (versão 2).** **aberta** (auditoria 16/09): a aritmética de composição em Lean está correta, mas as hipóteses P1–P3 NÃO são o caminho, porque P2 falha no modelo 1+∂τ (ANALYTIC_AUDIT) e C1-G está aberta. Rota substituta: cauda algébrica com pesos RT (ALGEBRAIC_TAIL, STRUCTURED_EXTERIOR), cujo majorante uniforme exige caixas ≳176×944. Histórico: composição fechada e ancorada ao tile; hipóteses PDE nomeadas | `TailBound.lean` compõe `epsilon_N` exatamente e `TailTileBridge.lean` exige `epsilon_N <= tailError`, de modo que o tile não pode declarar um controle de cauda que as constantes não sustentam. Os bounds viraram **P1/P2/P3** em [`TAIL_BOUND.md`](TAIL_BOUND.md). P2 é o difícil |
| C2 | Invertibilidade em todo o contorno | **01/10: na faixa A, substituída pelo Rouché de L, certificado ([`C1_REAVALIACAO.md`](C1_REAVALIACAO.md) §9–10: N_A = 3); revisado em 01/10 (versão 2).** infraestrutura em nível kernel. **Correção (auditoria 16/09):** além de `tailError`, ficam FORA do kernel C1-G, o erro do fundo, a cobertura do contorno e a identificação dos blocos com normas de operadores infinitos | `RationalCertificate.lean` verifica parametrices/Neumann por tile sem sair do kernel. `matrixVariation` **deixou de ser bound de EDP** e virou `def`: os campos são `tileRadius` e `matrixNorm`, e a não-negatividade da variação é agora **teorema**, não hipótese. Só `tailError` continua analítico |
| C3 | Winding 1 no setor A | **01/10: na faixa A, substituída pelo Rouché de L, certificado ([`C1_REAVALIACAO.md`](C1_REAVALIACAO.md) §9–10: N_A = 3); revisado em 01/10 (versão 2).** **maquinaria completa; o `3` é o limite, não artefato** | o pipeline roda da matriz até o kernel e `Gamma_A` dá winding **3**. Com quatro truncamentos, a terceira raiz **não** sai do contorno: instala-se em `s = 0,168311`, a `0,043` dentro da borda `sigma_0 = 1/8`. Minha leitura anterior, de que ela migraria para fora, está **refutada**. C3 continua exigindo uma cota inferior provada para o espectro físico não nulo |
| C4 | Multiplicidade **física** 0 no setor B (antes: "winding 0") | **02/10: certificada pelo Rouché de L em Ω_B = (0,3) × (1/4,3/4): N = 1 (0,0414 + 0,5i), raiz simples com testemunho de violação das constraints, logo multiplicidade física 0 ([`C1_REAVALIACAO.md`](C1_REAVALIACAO.md) §11); revisada em 02/10 (§11.1).** **alvo reformulado em 16/09**; regressão finita 4×12 apenas | `make verify-reduced-b` dá winding 0 em 4×12 sobre `∂([1/8,1]×[1/4,3/4])`, sem transferência. **O "winding 0 de `L_A` na faixa B" é numericamente FALSO:** de 6×18 a 12×36 há uma raiz convergente em `0,041421+0,500008i`, fora do retângulo, com assinatura de violação de constraints. Alvo correto: winding de `L_A` em `Γ_B(σ0)` menos as raízes com testemunho S3 igual a 0. Ainda descobertos: `0<Re s<1/8` (contém a raiz), `Re s=0`, `1<Re s<414`, as bordas `Im s=1/4, 3/4` fora de `[1/8,1]`, e toda a transferência infinita. Ver [`AUDITORIA_FISICA_16SET.md`](AUDITORIA_FISICA_16SET.md) §7. **17/09, regressão finita 6×18 (`make verify-reduced-b-6x18`), aceita pelo kernel:** winding **1** em `∂([1/64,1]×[1/4,3/4])` (1946 arestas, reduzidas a 87) e winding **1** na caixa `[1/32,1/10]×[23/50,1/2]` (450 arestas, reduzidas a 11). Na MATRIZ 6×18, a única raiz do retângulo é simples e fica na caixa. O polinômio característico exato foi conferido contra dois determinantes Bareiss exatos (`z=1/32`, `1/10`). Isso ainda não diz nada sobre o operador infinito nem sobre constraints |
| T1 | Transferência da contagem finita para o pencil | **01/10: na faixa A, substituída pelo Rouché de L, certificado ([`C1_REAVALIACAO.md`](C1_REAVALIACAO.md) §9–10: N_A = 3); revisado em 01/10 (versão 2).** formulada condicionalmente | formalização do lema Fredholm |
| T2 | Exatamente um modo físico instável | **03/10 (noite, 3): modelagem reduzida. H-cone virou o Lema C (o modo de μ* é tangente às translações de T*; contagens 1 e 2, ambas teoremas). Corolário do eixo imaginário: só a fase s = 0, simples (C1_REAVALIACAO §13). H-reg enfraquecida para H-reg′: C^∞ em τ, analítica em ξ. Restam: simetria esférica, H-reg′ (só em ξ) e identificação do fundo. Os três itens ainda sem revisão independente.** **03/10 (noite, 2): volta revisada (não quebrou; Lema V1: L e C juntos dão os dez δΩ nulos). A observação O1 da revisão (KC = ML com perda de raio) se fecha porque os certificados de K estão no toro♯(129/128, 9/8), exatamente o espaço de chegada de C (`SHARP_DOMAIN.md`, fechamento). Todas as peças matemáticas de T2 estão revisadas.** **03/10 (noite): H-rec (ida) virou o Lema R ([`HREC_IDA.md`](HREC_IDA.md)), revisado: todo modo físico de Floquet com Re s > 0 entra no ansatz por gauge analítico, com h ∈ D(Y+); a ressonância s ≡ μ fecha pela simplicidade de μ*; a recuperação de raio vale de qualquer peso > 1 e para todo s; as cadeias físicas têm comprimento 1 módulo gauge. Restam só a modelagem (simetria esférica, H-cone, H-reg, fundo) e a revisão própria da volta.** **03/10 (tarde): H-gauge virou o Lema G ([`GAUGE_GLOBAL.md`](GAUGE_GLOBAL.md)): os modos de gauge de Floquet com Re s > 0 são ℂ·g*, para gauges C¹ em todo o domínio de RT, sem analiticidade; fato certificado F1 (ω4*(·, 0) ≢ 0, |c_1| >= 0,493). Restam H-rec (ida), H-cone e H-reg.** **03/10: revisão independente de T2 feita (não quebrou; lacunas de enunciado e citação corrigidas) e NK em 0,7332: λ* ∈ [0,73316949; 0,73318983], simples ([`C1_REAVALIACAO.md`](C1_REAVALIACAO.md) §12). As quatro raízes estão localizadas rigorosamente.** **02/10 (noite): enunciado condicional em [`T2_ENUNCIADO.md`](T2_ENUNCIADO.md). As lacunas matemáticas L-real, L-afim, L-iso e L-constr estão fechadas ([`T2_HIPOTESES.md`](T2_HIPOTESES.md) §3.1); restam as hipóteses físicas H-rec, H-gauge, H-cone e H-reg.** **02/10: parte computacional fechada: N_A = 3 (μ gauge, s_K e 0,733) e N_B = 1 (viola as constraints), logo a contagem física é 1, condicionada às hipóteses de S4 (H-rec, H-gauge, H-cone) e de `SECTOR_B.md` §6. Faixas A e B revisadas; ficam as hipóteses físicas (lista em `C1_REAVALIACAO.md` §11.1).** aberta | depende de S3--S4, F1--F3 e C1--C4 (S2 saiu da lista: foi reduzida) |
| P2* | Regularidade analítica do **adjunto** | **eliminada** | majorante unilateral implementado: `maxDistance * tailError * (1 + maxDistance * inverseBound * operatorBound)`. Preço `N_F : 9 -> 13`. A regressão `rejectedUnilateralTile` exibe o preço: com os mesmos dados o bilateral passaria (`3/10`) e o unilateral recusa (`957/190`) |
| N-coh | Coerência de norma entre `tailError` e `matrixVariation` | **fechada** | `l^1` (soma por coluna) adotada; `inverseDefect` e `inverseBound` usam `oneNorm`. Único tile que discrimina é o 4x4, e ali `l^1` dá defeito **menor** (`545/16384` contra `279/8192`) |

## Sondagem de 20/09/2026, REFUTADA em 21/09: a faixa `0 < Re s < 1/8` segue descoberta, não povoada

Registrado aqui porque afeta C3, C4 e F3, e é independente da questão de custo
que motivou a medida. Varrendo a borda esquerda `sigma0` do contorno na malha
6×18 e medindo `||V||` (a norma da inversa diádica validada, `finite_inverse_bound`):

| `sigma0` | 1/8 | 1/10 | 1/16 | **1/32** | 1/64 | 1/256 |
|---|---:|---:|---:|---:|---:|---:|
| `||V||` | 3822 | 2399 | 3847 | **39183** | 5503 | 3413 |

O pico em `sigma0=1/32` é dez vezes o valor da borda atual e está fora de
qualquer tendência monótona na distância a `mu`, o que indica quase-singularidade
perto de `s ~ 0,03`. Essa faixa estava listada como **descoberta** em
`ALGEBRAIC_TAIL.md` e na linha C4; agora há evidência numérica de que é povoada.

**Consequência:** qualquer argumento que precise excluir espectro em
`0 < Re s < 1/8` — e C3 precisa, para uma cota inferior do espectro físico não
nulo — tem de lidar com isto, não pode supor a faixa vazia.

**A ressalva disparou. O refinamento refutou o achado.** Varrendo a mesma faixa
em três malhas:

| `Re s` | 4×12 | 6×18 | **8×24** |
|---|---:|---:|---:|
| 0,020 | 432 | 7229 | 5444 |
| **0,030** | 400 | **26295** | **3568** |
| 0,040 | 367 | 16453 | 2660 |
| 0,050 | 335 | 6466 | 2138 |
| 0,080 | 260 | 2699 | 1452 |

Em 4×12 não há feature; em 6×18 há um pico agudo; em 8×24 a curva volta a ser
**monótona decrescente**. O pico era artefato da malha 6×18, não raiz
convergente — o oposto do que aconteceu com a terceira raiz do setor A, que
apareceu e ficou.

**Conclusão corrigida:** a faixa `0 < Re s < 1/8` continua **descoberta** — não
há evidência de que seja povoada, e também não há prova de que seja vazia. O
estado da linha C3 não muda. Nada deve ser construído sobre o achado de 20/09.

Scripts: `scripts/measure_contour_border.py`, `scripts/scan_strip_sigma_min.py`.

## Lema de redução a certificar

Se `P(z)=P(z0)+(z-z0)J`, `P(z0)` é bijetivo, `R0=P(z0)^-1` é compacto e
`C=J R0`, escolha uma aproximação de posto finito `C_N` tal que

\[
  \|C-C_N\|\leq\varepsilon_N.
\]

Para um contorno `Gamma`, se

\[
 q_N=\sup_{z\in\Gamma}|z-z_0|\,\varepsilon_N
       \|[I+(z-z_0)C_N]^{-1}\|<1,
\]

a homotopia de Neumann liga `C` a `C_N` sem cruzar singularidades no
contorno. A multiplicidade total dentro de `Gamma` é então o winding do
determinante finito `det[I+(z-z0)C_N]`.

O ponto novo e difícil não é calcular o determinante: é produzir
`epsilon_N` uniforme para o fundo inteiro no ball RT e fechar todas as bordas
do domínio espectral.

[`formal/ChoptuikFormal.lean`](../formal/ChoptuikFormal.lean) formaliza a
contabilidade final e torna explícitas essas hipóteses. Ela não as postula como
axiomas e não finge que já foram demonstradas.

## O quadro com quatro truncamentos

Reproduzido de forma independente pelo orquestrador (raízes reais com
`Re s > 0,05`, por `scipy.linalg.eig` sobre as quatro matrizes):

| truncamento | dim | raízes |
|---|---:|---|
| 4x12 | 138 | 0,231793 · 0,438698 · 0,709627 |
| 6x18 | 333 | 0,144011 · 0,404385 · 0,731184 |
| 8x24 | 612 | 0,166389 · 0,401248 · 0,733091 |
| 10x30 | 975 | **0,168311** · 0,401024 · 0,733178 |

As três convergem em posição. O que as separa é **como** se comportam:

| raiz | posição | resíduo das constraints | sobreposição com gauge |
|---|---|---|---|
| `lambda≈2,67` | converge para 0,7332 | cai (fator ~5 por refinamento) | 0,117, estável |
| `lambda≈1,46` | converge para 0,4010 | **estagnado** (fator 1,00) | 0,029, nível genérico |
| `lambda≈0,61` | converge para **`mu_RT`** | cai | **0,916**, contra 0,032 genérico |

Ou seja, cada uma das três tem uma assinatura distinta, e as três assinaturas
são exatamente as que S3, S4 e o teorema-alvo preveem: uma raiz física, uma que
viola as constraints, uma alinhada com gauge. O quadro é **consistente** com
"exatamente um modo físico" — e consistência não é prova.

### Por que isto ainda não decide

O sub-agente que produziu o resultado registrou o veredito como **não decidiu**,
e a disciplina está certa:

- a coincidência com `mu_RT` é **pós-hoc** e não tem derivação neste repositório;
- os dois testes desenhados para decidir falharam cada um num controle. Em
  particular, a fixação de gauge de RT **move o autovalor físico** em `4,6e-2`,
  de forma convergente — então ela não serve como filtro;
- a quantificação "92 % dele é gauge" foi **corrigida pelo próprio agente**: o
  operador é não normal, a biortogonalidade dá `w_i·v_j = O(10^-15)` fora da
  diagonal, e os `0,916` são um ângulo de 24° entre autovetores de operador não
  normal, no par de pior condicionamento. A assimetria contra os 83° do modo
  físico sobrevive; a leitura percentual, não.

**Previsão registrada e REFUTADA.** O `12x36` (dimensão 1422) deu
`|s - mu_RT| = 2,21e-6`, não `< 10^-7`: a convergência estagnou (fator 1,7 no
último passo, depois de 504 no anterior). A identificação com `mu_RT` fica
**enfraquecida** — o modo está a poucos `10^-6`, mas não se pode mais afirmar que
converge para lá. Ver [`MODE_DISCRIMINATION.md`](MODE_DISCRIMINATION.md) §10.

## As constraints discriminam? Em parte, e só com dois truncamentos

Detalhes e números em [`MODE_DISCRIMINATION.md`](MODE_DISCRIMINATION.md).

**No `4x12`, não discriminam — e isso é demonstrável.** O argumento não é o valor
dos resíduos, é um **controle interno**: o modo de fase, que é gauge conhecido e
satisfaz as constraints por construção, tem o **pior** resíduo da tabela
(`4,43·10^-1`, quatro vezes o do candidato físico). Um teste que reprova o modo
bom conhecido não é um teste. As raízes estão isoladas, então não é contaminação
de autovetor: é erro de truncamento dominando.

**No `6x18`, discriminam — e para exatamente uma das duas raízes extras.**

| modo | `s`: 4x12 -> 6x18 | resíduo: 4x12 -> 6x18 | fator |
|---|---|---|---:|
| físico `lambda≈2,67` | 0,7096 -> 0,7312 | 1,14e-1 -> 1,94e-2 | **5,90** |
| segundo `lambda≈1,5` | 0,4387 -> 0,4044 | 3,23e-1 -> 3,19e-1 | **1,01** |
| terceiro `lambda≈0,5` | 0,2318 -> 0,1440 | 1,16e-1 -> 3,27e-2 | **3,56** |
| fase (gauge) | -0,0190 -> +0,0338 | 4,43e-1 -> 4,26e-2 | **10,41** |

**A distância entre 3 e 1 são dois problemas diferentes**, não um buraco só: uma
raiz que as constraints excluem (território de **S3**) e uma raiz que elas
**aceitam** e que ainda assim não converge de posição (território de **S4**, ou
controle de artefatos). Nenhum dos dois é do verificador de winding.

O que **não** está concluído, e é importante que não esteja: resíduo estagnado é
compatível tanto com uma raiz exata que viola as constraints (S3 falsa) quanto
com uma raiz espúria do truncamento, cujo autovetor não tem razão para satisfazer
nada. Dois truncamentos não determinam um limite. E a previsão de que o terceiro
modo seria o do tempo de acumulação (`lambda=1`) foi **refutada**: ele anda na
direção oposta.

## A primeira passagem pelos dados reais

`scripts/build_contour.py` leva a matriz `build/spectrum/rt-A-4x12.dat`
(dimensão 138) até um contorno consumível pelo produtor de certificados, em
aritmética exata. Dois resultados, ambos fechados pelo kernel:

| contorno | nós | classificam | winding |
|---|---:|---:|---:|
| `Gamma_A` completo | 1946 | 1946/1946 | **3** |
| entorno do candidato físico | 514 | 514/514 | **1** |

Três coisas que valem registro:

1. **Nenhuma aresta caiu em `none`**, o que valida na prática as duas regras de
   subdivisão. Mas não fechou na primeira tentativa: com 50 nós, **zero** de 50
   arestas classificavam e `K·L ≈ 7,3` contra o `< 1/2` exigido. A correção foi
   refinar o contorno, não afrouxar o critério.
2. O determinante relativo **telescopa** no pencil afim: `h(z_k) = p(z_k)/p(z_0)`
   com `p(z) = det(A + zI)`. Os dois caminhos dão os mesmos números; o custo
   difere por quatro ordens de grandeza (42 s por nó com Bareiss, contra 166 s
   uma vez por Hessenberg modular mais ~2 ms por nó). O caminho barato foi
   **verificado contra o caro** em dois pontos, com igualdade exata.
3. O `3` é o resultado que importa, e é honesto: a matriz truncada tem três
   modos com `Re s > 0` na faixa de A. A maquinaria de winding deixou de ser o
   gargalo de C3; o que falta é analítico.

### Uma previsão que não se confirmou

A §2.2 de [`WINDING_DATA.md`](WINDING_DATA.md) previa cruzamento **forçado** do
raio onde o contorno cruza `Im s = 0`. Não aconteceu: essas arestas deram
incremento 0 nos dois contornos. A previsão estava certa sobre `p`, e incompleta
sobre o que é de fato certificado: o certificado usa `h = p/p(z_0)` com `z_0`
complexo, e a normalização é uma rotação por constante não nula — não muda o
winding, mas **gira a reta de realidade para fora do raio**. Consequência dupla:
o risco de nó sobre o raio passa a ser controlável pela escolha do nó base, e em
troca a garantia de paridade de `Im p` se perde e passa a ser acidental.

## A costura entre a cauda e o tile

`TailData.epsilon` **calculava** `epsilon_N` e `Tile.tailError` **consumia** um
número, sem nada ligando os dois: um produtor podia declarar um `tailError`
menor que o `epsilon_N` que as suas próprias constantes produzem, e os dois
módulos aceitavam. `TailTileBridge.lean` exige `epsilon_N <= tailError` — o tile
pode ser conservador, nunca otimista. O sentido é o correto porque
`transferMajorant` é crescente em `tailError`: superestimar só dificulta
`transferMajorant < 1`; subestimar faz passar um certificado sem lastro.

Duas regressões documentam por que a checagem não é redundante:

- o tile mal ancorado (`tailError` uma ordem de grandeza **abaixo** de
  `epsilon_N`) tem `valid = true` e `transferMajorant = 147/190000`, que é
  **menor** que o do caso aceito (`147/19000`) — quem olhasse só o tile veria um
  certificado mais folgado, não um quebrado;
- um `TailData` inadmissível (`decayChebyshev = 1`) faz `epsilon_N` **diminuir**
  e `epsilon_N <= tailError` continuar valendo: a aritmética não denuncia, e só
  a checagem de admissibilidade pega.

## Construção do contorno: duas regras obrigatórias

Vêm da revisão cruzada e estão detalhadas em [`WINDING_DATA.md`](WINDING_DATA.md).
Como `C_N` é real e a caixa de truncamento é estável por `n -> -n`, vale
`f_N(conj z) = conj f_N(z)` **exatamente**. Daí:

- Em `Gamma_A` o eixo real cruza os dois lados verticais, `f_N` é real ali e,
  sendo a contagem alvo ímpar, os sinais em `sigma_0` e `R` são opostos. **O
  cruzamento do raio é forçado — e é dele que sai o `+1`.** Não é defeito.
- O único modo de falha real é um **nó** cair sobre uma reta de realidade onde
  `f_N < 0`: o valor fica *sobre* o raio, as duas arestas vizinhas devolvem
  `none` e **refinar não resolve**. Subdivisão uniforme com número **par** de
  segmentos no lado vertical produz exatamente isso.

Regras: (1) número **ímpar** de subdivisões nos lados verticais; (2) nós em
`±h`, aproveitando que `Im f_N(sigma_0 + iy)` é ímpar.

A economia de varrer meia faixa e dobrar por conjugação, sugerida antes em
[`SECTOR_B.md`](SECTOR_B.md) §7, é **incompatível** com o raio `R_{<0}`: ela
fecha sobre o eixo real, onde `f_N` é real num segmento inteiro. Se for adotada,
gire para o raio `i·R_{>0}`, o que não exige mudar o Lean — basta o produtor
emitir as caixas de `g = i·f`, e `i·([a,b]x[c,d]) = [-d,-c]x[a,b]` é exato em
`Rat`.

## Dois defeitos encontrados no que já existia

Ambos vieram da derivação de `epsilon_N` e nenhum estava registrado.

> **Resolvidos.** Os dois defeitos abaixo foram diagnosticados e **decididos**;
> falta aplicar as decisões ao código. Ficam registrados porque o raciocínio que
> levou a cada escolha é parte do ledger.

**P2\*, a hipótese tácita.** O truncamento **unilateral** `Pi C` tem o mesmo
determinante finito que o bilateral `Pi C Pi` (lema de Sylvester, demonstrado em
[`TAIL_BOUND.md`](TAIL_BOUND.md)), e só consome P2. Mas o majorante
`Tile.transferMajorant`, como está escrito em `RationalCertificate.lean`, é
majorante do truncamento **bilateral** — e esse exige regularidade analítica do
operador **adjunto**, que não se deduz de P2. Ou se troca o majorante pela versão
unilateral, ou P2\* entra no ledger como hipótese própria. Enquanto isso não for
decidido, o `tailError` que o tile consome não tem enunciado único.

**N-coh, a incoerência de norma.** Na realização de Wiener (`l^1` com pesos), a
norma de operador induzida é a máxima soma por **coluna**. `RatMatrix.infNorm`
calcula a máxima soma por **linha**, que é a norma induzida de `l^infinito`.
Enquanto as duas convivem sem conversão explícita, `tailError` e
`matrixVariation` podem estar medidos em normas diferentes, e a desigualdade de
Neumann que os soma não significa o que aparenta. A correção é barata (transpor,
ou fixar `l^infinito` em toda a cadeia e refazer as constantes de cauda), mas
tem de ser feita antes de qualquer certificado real.

## O que mudou de estado no winding

Até agora `WindingCertificate` somava incrementos `-1/0/+1` fornecidos por um
programa externo: a soma era conferida, mas a *origem* dos incrementos era um
oráculo, e `EdgeSignCertificate` não estava ligado a coisa alguma.

`formal/ChoptuikFormal/WindingInterval.lean` fecha esse buraco. Um
`IntervalWindingCertificate` carrega, por aresta, três caixas racionais
(`startBox`, `endBox`, `edgeBox`) e o kernel: (i) confere a forma e a inclusão
`edgeBox ⊇ startBox, endBox`; (ii) **calcula** o incremento de cada aresta pelo
cruzamento do raio `ℝ_{<0}`, em vez de recebê-lo; (iii) exige que a lista
encadeie e feche; (iv) prova `0 ∉ edgeBox` em todos os ramos, que é o
não-anulamento de `f` sobre `Γ`; (v) produz o `WindingCertificate` correspondente
com o mesmo total. `CertifiedSpectralObligations` em `ChoptuikFormal.lean`
consome esses certificados no lugar dos dois `Nat` asseverados.

A geometria deixou de ser hipótese. Cada aresta carrega agora `startPoint` e
`endPoint` exatos no plano `z`, e o kernel exige que o `endPoint` de uma aresta
seja **idêntico** ao `startPoint` da seguinte, que a última feche na primeira, e
que nenhuma aresta seja degenerada. Antes disso, um produtor com geometria
inconsistente podia emitir caixas que encadeavam perfeitamente sem descrever
curva alguma — as regressões `boxChainOnlyContour` e `degenerateContour` são
exatamente esses dois contornos, e são rejeitadas.

O que **continua hipótese analítica** do produtor externo, e está assim
documentado no módulo: que `startBox` contém `f(startPoint)`, `endBox` contém
`f(endPoint)` e `edgeBox` contém `f` sobre todo o segmento. Essa é a fronteira
com o estimador numérico, e é onde ela deve ficar declarada.

`make verify-winding` typechecka o Lean emitido pelo produtor contra o
verificador a cada execução. Foi a ausência dessa checagem que deixou os dois
divergirem silenciosamente uma vez.

A base de confiança do próprio verificador — quando o kernel checa e quando se
delega ao compilador — está registrada à parte em
[`KERNEL_TRUST.md`](KERNEL_TRUST.md).

## Rodada 4 (17/09): viabilidade da inversão certificada em G — INVIÁVEL SEM HPC

[`G_PARAMETRIX_FEASIBILITY.md`](G_PARAMETRIX_FEASIBILITY.md). A regra foi
pré-registrada na §0 antes dos números: memória ≤5 GB, ≤2 h por centro em 2
núcleos, ≤200 centros, defeito ≤0,9. Nenhuma arquitetura cumpre as quatro
condições. Os dados de matriz em float são diagnóstico.

- **Precedente RT.** Inversa densa exata em 25×75; inversa montada coluna a
  coluna em 250×750 (934500 colunas), com refinamento iterativo local em
  diádicos GMP; exterior coluna a coluna até 1600×4800. Contração final
  0,64, num ÚNICO ponto s, em HPC.
- **Trade-off.** 107×384 é a única caixa da grade com norma real de coluna
  ×1,5 abaixo de 0,9.
- **Arquiteturas.**
  - (a) Densa: 329 GB, ~5·10⁴ centros, e a 1ª ordem de RT não contrai.
  - (b) Núcleo denso com Neumann amortecida: ~3,7·10³ h por centro.
  - (c) Recomendada: Schur num núcleo 8×24, casca com Neumann uniforme em s
    por região e exterior coluna a coluna. Cabe em <1 GB e dá defeito
    ≤8·10⁻⁴ na casca, com ≈110–190 regiões, mas ≈180–255 h por região.
    Falha só no tempo, por ~100×. Em HPC custaria ~1–1,5·10⁵ núcleo-horas,
    trivialmente paralelas.
- **Protótipo em float** (12×36 com núcleo 6×18): defeito 0,033 (s=1) e
  0,027 (1/8+i/4). O defeito não é o gargalo; tempo × centros é.
- **Implicações.**
  - Obj. 1: a parametriz de U está fora de alcance nesta máquina. O
    exterior de (c), ~3 h de CPU, baixaria d* de 0,998 para ≲0,7.
  - Obj. 2 e 4: enclausurar os autopares exige ‖Q·H⁻¹‖~10³. Inalcançável
    aqui.
  - Obj. 3: (c) é a arquitetura certa. Ficaria viável nesta máquina com um
    pré-condicionador de casca que dispense o grau 32 (colunas de (B+s)Q
    ≲0,6 perto do núcleo; hoje chegam a ≈58). Esse é o problema de pesquisa
    de maior alavancagem.

## Rodada 5 (17/09): consolidação — exterior a d*<=0,7 e verificação das cotas

Direção escolhida pelo usuário depois da rodada 4: consolidar primeiro,
pesquisar o pré-condicionador de casca depois. HPC fica como alternativa a
solicitar.

- **Marco 8 do [`COLUMNWISE_EXTERIOR.md`](COLUMNWISE_EXTERIOR.md)
  (pré-registrado antes do cálculo):** para levar d* de 0,998 a <=0,7, os
  cortes passam a N_far=900 e M'=160, com G inalterada em 107×384. Cotas:
  d_R(900)=0,65116 (A) e 0,67304 (sharp); d_F(160)=0,66755 e 0,66784.
  Regiões novas R1, R2, R3 somam ≈211 mil colunas no A e ≈91 mil no sharp.
- **Verificação independente do Físico** ([`AUDITORIA_FISICA_16SET.md`](AUDITORIA_FISICA_16SET.md) §12):
  - **as quatro cotas CONFIRMADAS**, com os três lemas de coluna
    re-derivados de RT (22). Os valores do Físico ficam 5,1e-5 ABAIXO, que
    é a folga dos majorantes racionais (99/70>=√2 e afins);
  - **correção:** o que é exatamente proporcional a 1/M é `d_F − e_μ`, não
    `d_F`; e_μ vale 4,2e-9 (A) e 7,9e-9 (sharp);
  - **ν_K = 2/√3 é conservador em 15%** no sharp, sem efeito aqui;
  - **cobertura CONFIRMADA:** toda coluna fora de G cai em exatamente uma
    região ou numa cauda, com as fronteiras do lado certo;
  - **modo "completo" do núcleo C: confirmado com ressalva.** Os valores do
    produtor são sempre >= os independentes (direção certa), mas o viés é de
    8e-4 a 1,5e-3 relativo, cinco ordens acima do arredondamento, exceto em
    n=0, onde as imagens coincidem (9e-11). Suspeita: cancelamento perdido
    entre imagens em n para n<101. Não ameaça a meta (R3 dá ≈0,12 contra
    0,7), mas a validação anterior não é concordância no nível do
    arredondamento.
