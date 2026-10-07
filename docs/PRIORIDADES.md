# Sequência de prioridades

Documento de sequenciamento, criado em 20/09/2026. Não é um registro
cronológico (isso é `RETOMADA.md`) nem o ledger de obrigações (isso é
`PROOF_OBLIGATIONS.md`). É a resposta a uma pergunta só: **o que fazer em
seguida, e por quê.** Atualize o estado dos itens conforme forem fechando.

## DECISÃO PENDENTE (03/10/2026): o que fazer a seguir

Ver [`PANORAMA.md`](PANORAMA.md), que junta tudo o que foi produzido e lista os caminhos com custo e
valor:
- **A.** Artigo de T2.
- **B.** T3 não esférico.
- **C.** Consolidação e reprodutibilidade.
- **D.** Lean da cadeia final.
- **E.** Matemática além do espectral.
- **F.** Limpeza.

A sugestão é C (com F) → A → B. O usuário quer olhar o conjunto antes de decidir.

## ATUALIZAÇÃO 03/10/2026 (noite, 4): não esférico, Fases 0 e 1 feitas

[`NAO_ESFERICO.md`](NAO_ESFERICO.md): as equações de Gerlach–Sengupta, axiais e polares, foram reescritas
no fundo de RT, só com ω e μ. A numérica não rigorosa reproduz Martín-García e Gundlach (1999) em ℓ = 2 e
3, nos dois setores.
- **O ponto crítico é o modo polar ℓ = 2:** s = −0,004764 ± 0,147345i (κΔ = −0,0599). É amortecido, mas
  fica a 0,0048 do eixo.
- A fase rigorosa (T3) teria o tamanho do projeto esférico:
  - formulação bem condicionada (Galerkin com o inverso livre);
  - vínculos;
  - NK no modo de ℓ = 2;
  - cota uniforme em ℓ.

Também: revisão leve dos itens de modelagem feita. L1: o item 3 de T2 vale só para as raízes de L.

## ATUALIZAÇÃO 03/10/2026 (noite, 3): modelagem reduzida

1. **Lema C** (`GAUGE_GLOBAL.md` §6): o modo de μ* é tangente às translações de T*, que são soluções
   exatas isométricas ao fundo. H-cone deixa de ser hipótese, e T2 enuncia as duas contagens: 1 com o
   gauge linearizado padrão, 2 com gauges que fixam o cone.
2. **Eixo imaginário** (`C1_REAVALIACAO.md` §13): pelos Rouchés das faixas A e B, a única raiz em
   Re s = 0 é a fase s = 0, simples. Não há modos neutros físicos.
3. **H-reg′:** C^∞ em τ e analítica em ξ, porque R5 vale com peso de Fourier 1. Fica aberta só a
   regularidade C^∞ em ξ.

Pendências: a revisão independente dos três itens, que é leve, e a análise da simetria esférica.

## ATUALIZAÇÃO 03/10/2026 (noite, 2): volta revisada; todas as peças de T2 revisadas

- **Volta** (`GEOMETRIC_FLOQUET_RECONSTRUCTION.md`, revisão em `.codex-runs/2026-10-03-revisao-V/`): não
  quebrou.
  - O enunciado ficou preciso.
  - O Lema V1 diz que L e C juntos dão os dez δΩ nulos.
  - Para cadeias, valem as constraints de jato.
- **O1 da revisão:** KC = ML perde raio, e C(s)h cai em ℓ¹(129/128, 9/8). É exatamente o espaço em que S3a,
  S3b e R♯ certificam K injetivo, porque o toro♯(129/128, 9/8) contém esse ℓ¹ com norma <= 1. Fechado em
  `SHARP_DOMAIN.md`.

T2: "exatamente 1 modo instável físico", condicionado só à modelagem (simetria esférica, H-cone, H-reg e
identificação do fundo), com todas as peças matemáticas revisadas.

## ATUALIZAÇÃO 03/10/2026 (noite): H-rec (ida) provada (Lema R)

- **Lema G revisado** (`.codex-runs/2026-10-03-revisao-G/`): não quebrou. Ganhou a versão G′, com gauge só
  contínuo no cone e o fato F1′.
- **Lema R** ([`HREC_IDA.md`](HREC_IDA.md), revisado em `.codex-runs/2026-10-03-revisao-R/`): não quebrou.
  - Todo modo físico de Floquet com Re s > 0 entra no ansatz de RT por gauge analítico, com h ∈ D(Y+).
  - Peças novas:
    - a ressonância s ≡ μ, fechada pela simplicidade de μ* (o gauge secular daria uma cadeia de Jordan);
    - a recuperação de raio de qualquer peso > 1, para todo s;
    - as cadeias físicas (R7).

T2 fica condicionado só à modelagem: simetria esférica, H-cone, H-reg e identificação do fundo.

Próximos passos possíveis:
1. Revisão independente da volta (`GEOMETRIC_FLOQUET_RECONSTRUCTION.md`).
2. Formalização ou escrita do artigo.

## ATUALIZAÇÃO 03/10/2026 (tarde): H-gauge provada (Lema G)

[`GAUGE_GLOBAL.md`](GAUGE_GLOBAL.md): para todo campo de gauge C¹ em M × S² que preserva o ansatz, os modos
de gauge de Floquet com Re s > 0 são ℂ·g*, em s ≡ μ.
- A prova é global, sem analiticidade dos germes.
- Usa um único fato certificado, F1: ω4*(·, 0) ≢ 0.
- Ideia: as curvas de nível nulas são conexas em M, logo X = f(u); a equivariância sob Θ² torna f
  homogêneo; e ser C¹ no cone, que está dentro de M, força f constante.

Próximos passos:
1. Revisão independente do Lema G.
2. H-rec (ida): a existência da fixação de gauge.

## ATUALIZAÇÃO 03/10/2026: T2 revisado e a raiz física localizada

- **Revisão independente de T2** (`.codex-runs/2026-10-02-revisao-T2/`): não quebrou. O L-real, a redução
  B → A, o quociente e KC = ML conferem. As 11 lacunas eram de enunciado e citação, e foram corrigidas
  (classe de período completo de largura 1/2 com s_B real, hipóteses escondidas, H-rec e H-gauge como
  matemática não feita).
- **NK em 0,7332:** λ* ∈ [0,73316950; 0,73318982], simples (`C1_REAVALIACAO.md` §12).

Com isso a parte assistida por computador de T2 está completa. O que resta é matemática fora do
certificado: H-rec (ida) e H-gauge (global), e a formalização dos passos em papel 1–3 e 12.

## ATUALIZAÇÃO 02/10/2026 (noite): T2 enunciado, com as hipóteses físicas explícitas

[`T2_ENUNCIADO.md`](T2_ENUNCIADO.md) traz o teorema condicional: 1 modo instável físico, sob H-rec,
H-gauge, H-cone e H-reg. [`T2_HIPOTESES.md`](T2_HIPOTESES.md) traz o inventário e os lemas fechados:
- **L-real:** as contagens no toro valem na realização de RT, por recuperação de raio com contração 0,371
  no toro e 0,267 em RT;
- **L-afim;**
- **L-iso e L-constr.**

Próximos passos possíveis:
1. ~~Revisão independente de T2_HIPOTESES e T2_ENUNCIADO~~ **feita em 03/10**.
2. ~~NK em 0,7332~~ **feito em 03/10**.
3. As hipóteses físicas: H-rec (ida) e H-gauge (global) são problemas de análise geométrica fora do
   certificado.

## ATUALIZAÇÃO 01/10/2026: o Rouché de L está certificado

N(L; 0 < Re s < 3, |Im s| < 1/4) = 3, com 96 discos (θ <= 0,99916), cobertura do contorno, blocos de
janela sem zeros, certificado do fator e contagem finita por Schur ([`C1_REAVALIACAO.md`](C1_REAVALIACAO.md)
§9). Próximos passos, em ordem:

1. ~~Revisão independente do Rouché de L~~ **feita em 01/10** (`C1_REAVALIACAO.md` §10): L1–L6
   consertadas, 96 discos refeitos (versão 2), N_A = 3 mantido. Texto original do item:
   revisão independente, no protocolo usual: snapshot sha256 e artefatos por alegação.
   Pontos de maior risco:
   - as fórmulas de sensibilidade (disco e janelas);
   - o certificado do fator;
   - a forma de homotopia com referência s-dependente;
   - a contabilidade de Brauer.
2. ~~Faixa B (C4)~~ **certificada em 02/10** (`C1_REAVALIACAO.md` §11): N_B = 1, raiz simples que viola as
   constraints, logo multiplicidade física 0. **Revisada em 02/10** (§11.1; L1–L4 de rastreabilidade). Itens originais da revisão:
   - a Loewner pré-condicionada;
   - o reaproveitamento da referência e dos discos da faixa A;
   - a junção dos checkpoints;
   - o tQ no NK;
   - o testemunho em λ0 complexo.
3. A parte física de S3 (constraints) e T2.

## ATUALIZAÇÃO 26/09/2026: os itens 1–4 foram substituídos

A transferência finito → infinito (C1, C1-G, C2, C3 e T1) passa a ser feita por **Rouché de
operadores com referência finita**, a rota que certificou K em S3b, com a cauda de L tratada como
em S3b e S4 e a deflação de Brauer das raízes de gauge. O piloto deu θ0 = 0,81–0,89 no contorno
∂([0, 3] × [−1/4, 1/4]). O custo estimado é de ~300 tiles, 2–3 dias em três máquinas, sem HPC.
Com isso, a meta de β e a caixa de 556×3001 abaixo deixam de ser o caminho. Ver
[`C1_REAVALIACAO.md`](C1_REAVALIACAO.md). O texto abaixo fica como registro.

## O achado que organiza a lista

Das cinco obrigações analíticas em aberto que bloqueiam o resultado — S3, S4,
C1, F2, F3 —, **três são o mesmo problema**: C1, F2 e F3 se reduzem a *tornar
o majorante do exterior menor que 1*. F3 já está parcialmente fechada (a cauda
algébrica exclui `Re s >= 414`); o que falta dela, a faixa `1 < Re s < 414` e a
fita entre o eixo imaginário e `sigma0 = 1/8`, também depende de um majorante
melhor.

E esse problema tem métrica, com histórico de progresso real. Defeito TT em
G = 64×192:

| estado | defeito TT | ganho acumulado | fonte |
|---|---:|---:|---|
| uniforme, com `||B|| <= 23` | 13,712 | — | `FREQUENCY_BLOCKS.md` |
| cauda estruturada (cancelamento nas colunas) | 10,177 | 1,35× | `STRUCTURED_EXTERIOR.md` |
| + norma adaptada por componentes | **4,882** | 2,81× | `COMPONENT_EXTERIOR.md` |
| **alvo** | **< 1** | faltam ~4,9× | — |

Contraste com a rota que foi fechada: no pré-condicionador de casca faltavam
57× e o teto provado era 29–32×. Aqui faltam 4,9× e **não há teto provado** —
é exatamente isso que o item 1 abaixo vai determinar.

## Por que este é o item certo, e não computação

O tamanho da caixa é determinado pelo majorante, não pelo gosto: com
`||B|| <= 23` a condição de transferência exigia 865×2594; com a cauda
estruturada caiu para 364×1964. As horas-núcleo escalam com a caixa.

Portanto o majorante é **o único item da lista que ataca a obrigação
matemática e o orçamento ao mesmo tempo**. O pré-condicionador de casca só
atacava custo, e por isso seu fechamento negativo não doeu tanto quanto
parecia; este item, se ceder, resolve os dois.

## A sequência

### 1. Análise de teto do defeito TT — FECHADO POSITIVO (20/09/2026)

**Não há teto.** Ver `TT_CEILING.md`. O critério de blocos é `rho(E)<1` sobre
a matriz majorante de blocos, e nela `e_TS` não depende da caixa enquanto
`e_ST` e `e_TT` são proporcionais a `t(G) -> 0`: a condição passa para toda
caixa grande o bastante. O que existe é uma **caixa exigida**, com
`DOFs ~ beta^2,8`. Hoje, com `beta=10,382`, ela é 556×3001 — 40,6× os DOFs de
107×384. **Alvo: `beta <= 2,064` zera o custo extra**; `beta <= 3,941` o
deixa em 4×.

O item 2 abaixo está LIBERADO, e deixou de ser "melhorar o majorante" para
ser **"reduzir `beta` de 10,382 para 2,1"** — um fator 5,0.

<details>
<summary>Enunciado original do item</summary>

Antes de gastar semanas empurrando o majorante, determinar se ele tem piso,
pelo mesmo método que fechou a rota k=0 em um dia: reponderar por DOF é uma
similaridade diagonal, e por Perron–Frobenius o mínimo da norma de coluna
sobre todos os pesos positivos é `rho(|X|)`. Se algum subbloco finito do
exterior já tiver raio de Perron `>= 1`, **nenhuma reponderação** salva a rota
e o item 2 não deve ser tentado.

Critério de decisão: `rho >= 1` fecha a rota da reponderação pura;
`rho < 1` libera o item 2 e dá o alvo a perseguir.
</details>

### 2. Reduzir `beta` — PARCIALMENTE FECHADO (rodada 2, 20/09/2026)

Ver `TT_CEILING.md` §§7–10. **2c está morto** (reponderação já a 0,09% do
piso de Perron `rho(M)=10,3717`). **2b vale 1,24×** no melhor caso concebível,
não o fator 5. **O núcleo, que ninguém tinha otimizado, rendeu 1,7× de
graça**: a caixa exigida hoje é 23,4× os DOFs de 107×384, não 40,6×.

Para o alvo `beta<=2,064` seria preciso absorver ~93% do termo não linear.
A massa está em `h2` e `h4`; zerar o bloco de `h2` sozinho leva `rho` a 4,265.

Próximo passo antes de qualquer decisão de compra: fechar o dimensionamento
com os acoplamentos `P<->T` reais e `||V||>1`. Ambos empurram a caixa para
cima; 23,4× é piso.

- **2a.** Converter os blocos FF, FT e TF para a norma adaptada por
  componentes. `STRUCTURED_EXTERIOR.md` chama isso de obrigatório: hoje o TT
  vale 4,882 nessa norma e os demais blocos continuam na norma antiga, e o
  critério `b*c < (1-a)(1-d)` não pode ser aplicado com blocos em normas
  diferentes. É trabalho delimitado, não pesquisa.
- **2b.** Absorver parte de B numa parametriz exterior em vez de jogar tudo em
  `||B|| <= 23`. `COMPONENT_EXTERIOR.md` já exibe o principal triangular por
  componentes, `diag(K_mu, J_mu+mu R_x(1-P), J_mu, J_mu+mu R_x(1-P))`, e diz
  explicitamente que **a inversa desse novo principal ainda não foi
  certificada**. Há um fator ali contado como zero só porque ninguém foi
  buscá-lo. Pelo escalonamento `beta^2,8`, um fator 2 aqui vale 7× em DOFs.
- **2c.** Apertar a própria constante. A norma ponderada exata de `B_6x18`,
  medida, é **8,078**, contra o bound 10,382 em uso — há 1,29× de folga só em
  contabilidade, antes de qualquer ideia nova.

### 2bis. `||V||` — MEDIDO E FECHADO, NEGATIVO (rodada 4-5, 20/09/2026)

`TT_CEILING.md` §§17–23. Os 80 do registro eram `s=1`, o ponto mais longe de
qualquer raiz. Na borda do contorno `||V||` é **3822** (6×18), e cresce com a
malha. A caixa exigida vai a **9127× a 14607×** os DOFs de 107×384.

A alavanca `sigma0` foi medida e também está morta: `||V||` não é monótono na
distância a `mu`, e o pico em `sigma0=1/32` (39183) mostra que a faixa
`0 < Re s < 1/8` é povoada. Melhor ponto de toda a varredura: 1,6×.

**Todas as alavancas baratas estão exauridas.** A rota não tem teto matemático,
mas a distância à viabilidade é de três a quatro ordens de grandeza, não de um
fator 5. Continuar por aqui exige ideia estruturalmente diferente — não ajuste
de peso, norma, ponto-base ou contorno.

<details>
<summary>Enunciado anterior do item</summary>

Rodada 3 (`TT_CEILING.md` §§11–16): com `||V||=80` medido, o dimensionamento
fechado dá **332× os DOFs de 107×384**. `||V||` de 80 para 1 vale **30×**;
`beta` de 10,372 para 5 vale 3,7×. A prioridade trocou.

Dois experimentos baratos, antes de qualquer teoria:
- **medir `||V||` em 12×36**, que já existe, para saber se cresce com a caixa;
- **medir `||V||` perto da raiz `s=0,168311`**, a 0,043 da borda `sigma0=1/8`,
  porque é o `||V||` UNIFORME no contorno que entra, e perto de raiz ele
  cresce. Este é o risco que pode piorar tudo.

Fechados nesta rodada, para não reabrir: ponto-base `z0` (negativo, mesmo com
extensão hipotética da cauda estruturada), reponderação por componentes
(negativo, 0,09% de folga), separação aditiva (colhido) e razão de aspecto
5,4 (colhido).

</details>

### 3. Fechar C1-G — depende do item 2

A ponte até o determinante implementado. `ALGEBRAIC_TAIL.md` já dá o caminho
alternativo que evita a resolvente completa: usar a inversa LIVRE `Q0` como
pré-condicionador, para o qual `Q0 iota = iota Q0_N` é exato e Sylvester
reproduz o winding do mesmo pencil finito. Falta a condição quantitativa
`b*t_in*(1+b*q0*g) < 1`, que hoje falha já no primeiro fator em todas as
malhas disponíveis.

### 4. F3 e F2 — caem junto com o item 3

Com majorante bom, a faixa `1 < Re s < 414` e a fita até `sigma0` deixam de
exigir argumento novo. F2 (coercividade de alta frequência) é da mesma família.

### 5. S3 e S4 — FECHADOS (24–26/09/2026), com revisão independente

A expectativa abaixo estava calibrada para o pior caso; as duas saíram antes do previsto, pela rota
de enclosure certificado:

- **S3a** (24/09): 154 tiles de Perron; K injetivo em 0 <= Re s <= 1,765, |Im s| <= 1/4, fora do
  disco |s − 0,401| < 1/8.
- **S3b** (25/09, `AUDITORIA_S3B_25SET.md`): K tem uma raiz no disco, s_K, simples (Rouché); L tem
  s_K simples (NK com bloco-Jacobi); o testemunho viola as constraints.
- **S4** (25–26/09, `S4_GAUGE.md`, `AUDITORIA_S4_26SET.md`): μ* é raiz simples de L, com
  autoespaço span(g*) e D g* = 0 (NK com centro no gauge, RefB explícito, lema de Corr).

O que resta para T2 está na linha S4 de `PROOF_OBLIGATIONS.md`: C3/T1 (contagem N_A na faixa A,
com transferência finito → infinito), C4 (faixa B), F3 (1,765 < Re s < R) e as hipóteses físicas
H-rec, H-gauge, H-cone e H-reg.

Texto original do item, mantido para registro:


São de outra natureza: não são majorantes, são a identificação do quociente
físico. Muita maquinaria pronta — `DOmega(w)[e^{mu tau}(Z+1)w] = e^{mu tau}(Z+1)Omega(w)`
verificada universalmente, a distância 2,21e-6 da raiz a `mu` decomposta
exatamente em truncamento (+2,2126e-6) e defeito do fundo (+1,1e-10) — mas o
que falta é enclosure certificado de autovetor **e** multiplicidade algébrica
1, mais fixação global de gauge. Isso é pesquisa, não execução. Não começar
por aqui.

## Decisões já tomadas, para não reabrir

- **Pré-condicionador de casca e rota k=0: fechados, negativos.** Não reabrir
  sem uma ideia fora da família diagonal/blocos/deflação. Ver
  `SHELL_PRECONDITIONER.md`.
- **Máquinas locais: esgotadas.** Notebook e PC do laboratório medidos em
  17/09/2026; o melhor dos dois viola o critério de 2 h por região por 27× e
  levaria 338–576 dias. Ver `.codex-runs/2026-09-17-lab-benchmark/RELATORIO.md`.
- **Cluster PXE no laboratório de informática: rejeitado em 20/09/2026.** Não
  por dificuldade técnica, e sim por aritmética: 20 máquinas usadas 2 h por
  sessão, 3 vezes por semana, rendem ~348 h-núcleo por semana, menos que as
  ~494 de uma única máquina ligada 24/7. O ciclo de trabalho domina a contagem
  de núcleos. Só revisitar se as máquinas puderem ficar ligadas sem supervisão
  — aí o cenário vira 7 a 12 semanas e passa a valer. Nesse caso, preferir
  pendrive live com `toram` a PXE: elimina DHCP, TFTP, NFS e a dependência de
  rede, e PXE por WiFi praticamente não existe em hardware comum.
- **O orçamento de `~/possivel_projeto.docx` está OBSOLETO (20/09/2026).** Ele
  foi calculado para a caixa 107×384; o dimensionamento fechado mostra que
  fechar C1 exige ~332× isso. A conta de nuvem iria de US$ 626–1.066 para a
  ordem de US$ 200–350 mil, e o SDumont de 24–41 mil UA para 8–14 milhões.
  **Não submeter como está.** Revisar só depois que `||V||` for atacado.
- **Compra de HPC: parada até o item 2bis decidir.** As três opções levantadas em
  17/09/2026 estão em `~/possivel_projeto.docx`: SDumont/LNCC (grátis, fluxo
  contínuo, 24–41 mil UA), servidor dedicado (€250–426) e nuvem interrompível
  (US$ 626–1.066, 8 a 53 dias). A validação de uma região custa US$ 4 e pode
  ser feita a qualquer momento, independente de tudo.

## Frente aberta com progresso: F3 pelo campo de valores (21/09/2026)

A única frente do projeto que **melhorou** nesta semana. Ver
`ALTERNATIVE_ROUTES.md` §§5.9–5.17 e a linha F3 de `PROOF_OBLIGATIONS.md`.

| peça | estado |
|---|---|
| campo de valores é a estimativa de energia que F3 pede | **estabelecido** |
| ponte `l^1 -> l^2` | **gratuita**, `sum(w|x|)^2 <= (sum w|x|)^2` |
| `R_J < 0,08885` | **CERTIFICADO**, exato, incluindo a cauda |
| `||B||_2` pela rota `l^1` | estrutura do argumento pronta, falta redigir |
| `||B||_2` pela elipse | folga de 6,7× **medida** no fundo, transferência por escrever |
| transcrição para o Lean | não começada |

`R` sairia de **414** para **~10,5** pela rota `l^1`, ou **~4,7** pela elipse,
contra **~2,9** medido. A região descoberta de F3 encolhe na mesma proporção.

**Próximo passo:** redigir a estimativa bilinear de RT com o supremo na elipse
no lugar da soma `l^1`. É onde a folga de 6,7× ou se realiza ou se perde.

## Regra de sequenciamento

Não pular para computação enquanto o item 1 não fechar. O motivo está na seção
"Por que este é o item certo": se o majorante ceder, a caixa encolhe e o
orçamento encolhe junto. Comprar horas-núcleo antes disso é otimizar a
restrição que não está mordendo.
