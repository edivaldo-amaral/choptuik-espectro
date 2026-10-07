# Ledger de confiança do verificador formal

Este documento registra *em que* o projeto confia quando diz que um certificado
foi "verificado exatamente". Ele não trata da matemática do problema espectral;
trata da base de confiança do Lean. É uma obrigação distinta das de
[`PROOF_OBLIGATIONS.md`](PROOF_OBLIGATIONS.md) e estava até agora implícita.

## O problema

`scripts/check_formal.sh` rejeita `sorry`, `axiom` e `admit`. Isso é necessário
mas não suficiente: **`native_decide` não é nenhum dos três e ainda assim
acrescenta um axioma**. Ele delega a avaliação ao compilador Lean e ao runtime
compilado, fora do kernel, e sela o resultado num axioma opaco.

Medição direta neste toolchain (Lean 4.33.1):

```
theorem trustProbeNative : regressionTile.valid = true := by native_decide
#print axioms trustProbeNative
-- 'trustProbeNative' depends on axioms:
--   [propext, trustProbeNative._native.native_decide.ax_1_1]
```

O segundo axioma é gerado por declaração e não aparecia em nenhuma checagem do
projeto. Havia **10 ocorrências de `native_decide` em `formal/`**, todas nas
regressões de `RationalCertificate.lean` — ou seja, a frase "o kernel recomputa
cada produto de matriz, norma e desigualdade" valia para as *definições*, mas as
*regressões que as exercitam* eram fechadas fora do kernel.

## Os dois bloqueios, e por que só um deles é real

Trocar `native_decide` por `decide` falha no arquivo original. Foram isolados
dois bloqueios — e **um deles não é intransponível**, ao contrário do que a
primeira versão deste documento afirmava.

### Bloqueio 1: `Rat` — aparente, removível com `unseal`

A aritmética `Rat` não reduz durante a elaboração porque `Rat.mul`, `Rat.add`,
`Rat.sub` e `Rat.inv` são `@[irreducible]` no Lean 4.33. **O kernel nunca
respeitou irredutibilidade** — ela é uma diretiva de elaboração. Logo `unseal`
devolve a transparência sem acrescentar hipótese nenhuma: o kernel continua
checando o termo inteiro.

```
theorem a  : ((1/2 : Rat) < (2/3 : Rat)) := by decide                -- FALHA
unseal Rat.mul Rat.add Rat.sub Rat.inv in
theorem a' : ((1/2 : Rat) < (2/3 : Rat)) := by decide                -- PASSA
unseal Rat.mul Rat.add Rat.sub Rat.inv in
theorem c' : (mkRat 1 10 + mkRat 1 5 < mkRat 1 2) := by decide       -- PASSA
```

Sem `unseal`, apenas comparações de `mkRat n d` reduzem; com `unseal`, a
aritmética completa reduz — **inclusive a divisão**.

> **Correção registrada.** A primeira versão deste documento afirmava que
> "qualquer certificado cuja aritmética seja `Rat` está estruturalmente
> condenado ao `native_decide`", e concluía daí que seria preciso migrar para
> aritmética dyádica. As duas coisas são falsas. A afirmação foi refutada por
> medição e a conclusão de projeto que dela dependia foi descartada.

### Bloqueio 2: `Array` — real

```
def arr : Array Int := (Array.range 8).map (fun (i : Nat) => Int.ofNat i * 3)
theorem a1 : (arr[3]?.getD 0) = 9 := by decide
-- Tactic `decide` failed ... reduction got stuck
```

Este não sai com `unseal`: `Array.range` e `Array.map` são `semireducible`, não
`irreducible`, então não há atributo a remover — o que não reduz é a própria
definição. Literais `#[...]` e `.size` reduzem; `Array.map` não. Como
`RatMatrix.data : Array Rat` e `mul`/`sub`/`identity` passam por `Array.map`,
as regressões de `Tile` continuam presas mesmo com `unseal Rat.*` (verificado:
as cinco falham).

**A correção é trocar `Array` por `List`.** Com `List` + `Rat` + `unseal`, uma
réplica de `Tile` 4x4 — incluindo a divisão `‖V‖/(1−η)` — fecha por `decide` em
5,8 s com axiomas `[propext, Classical.choice, Quot.sound]`.

## Custo

Medições neste toolchain, produto de matriz cheia mais norma-infinito, tudo por
`decide` (requer `set_option maxRecDepth 100000`):

| representação | dimensão | tempo | axiomas |
|---|---:|---:|---|
| `List` + `Rat` + `unseal` | 4x4 | 5,8 s | `[propext, Classical.choice, Quot.sound]` |
| `RationalCertificate.lean` migrado | módulo inteiro, do zero | 46 s | padrão |

> **Atenção: as medições abaixo são SINTÉTICAS.** Foram feitas com matrizes de
> entradas pequenas (`(i mod 7 + 1)/128`, ~3 bits). Elas dizem o custo da
> *representação* — `List` + `Rat` + `unseal` contra `List` + dyádico — e **não**
> o custo sobre os dados reais, cujas entradas têm 54 a 75 bits. Para a parede
> real, veja a seção "A parede do tile" adiante.

Por regressão isolada, descontado o baseline de import (2,25 s), com a matriz
sintética:

| tile | tempo |
|---:|---:|
| 4x4 | 0,8 s |
| 8x8 | 11 s |
| 12x12 | 47 s |
| 13x13 | 68 s |
| 14x14 | 84 s |
| 16x16 | 138 s (só com `maxHeartbeats 0`) |

E, na variante dyádica sobre `List`, também sintética:

| tile | tempo | axiomas |
|---:|---:|---|
| 8x8 | 12 s | `[propext]` |
| 16x16 | 110 s | `[propext]` |

A razão 9,2x entre 8x8 e 16x16 é compatível com o O(n³) dos produtos mais o
crescimento em bits das mantissas.

**Extrapolação honesta.** As dimensões reais do projeto são 138 (4x12), 333
(6x18) e 612 (8x24). Em O(n³), 612 custa cerca de `(612/16)³ ≈ 5,6·10⁴` vezes o
caso 16x16 — da ordem de 10⁵ horas. Fechar a matriz cheia no kernel é inviável e
continuará inviável com qualquer margem realista de otimização. A consequência
não é desistir do nível K: é que **o certificado precisa ser decomposto em tiles
pequenos**, onde o kernel fecha. Essa é uma restrição de desenho, não um detalhe
de implementação.

### A parede do tile, medida sobre dados reais (obrigação K5)

Esta seção já esteve errada duas vezes. O registro de ambas fica, porque o modo
como cada uma caiu é informativo.

**Erro 1 — extrapolação.** O documento trazia "12x12 em 47 s", medido numa
matriz **sintética** com entradas de ~3 bits. As entradas reais têm 54 a 75 bits
(mediana 67, expoente diádico `-72`), e o número não se transporta.

**Erro 2 — a hipótese dos bits.** Daí concluí que a parede era o **comprimento
das entradas**, e que a saída seria reduzir bits nos dois lugares onde o
certificado tolera erro (`V`, que é um chute, e o `centerMatrix`, cujo
arredondamento pode ser absorvido). A primeira parte se implementou
(`build_tile.py --bits/--center-bits`, com o descarte emitido em
`centerResidual`); a segunda **não se confirmou**.

Medições com tiles construídos de `rt-A-4x12.dat`, ajustes padrão, uma de cada
vez para não contaminar tempos:

| dimensão | exato | 16 bits |
|---:|---|---|
| 8 | fecha | — |
| 20 | **fecha, 72 s** | fecha, 74 s |
| 40 | falha por `maxHeartbeats` em `whnf` | falha |
| 60 | falha | falha (179 s, levemente pior) |

Cortar de 67 para 16 bits não muda nada. A hipótese está refutada.

**O que de fato se sabe.** Todos esses tiles são **aritmeticamente válidos** em
`Fraction` — a matemática fecha em todas as dimensões testadas; o que não fecha é
o orçamento de elaboração. A falha é sempre `maxHeartbeats` em `whnf`, isto é,
**orçamento por declaração**, e não um teto do kernel.

Isso é exatamente o diagnóstico que apareceu em K6 para o certificado de winding,
onde a solução não foi inflar o orçamento nem reduzir dados: foi **decompor em
declarações menores**, o que resultou mais barato que o monolítico. A estrutura
do tile favorece a mesma saída, porque `infNorm` e `oneNorm` são um **máximo
sobre linhas** de somas de valores absolutos — decompor por blocos de linhas é
aritmeticamente trivial e daria uma declaração por bloco. Nada disso está
implementado, e é o próximo passo natural de K5.

### A saída de K5: decompor por blocos de colunas

`formal/ChoptuikFormal/TileBlocks.lean` faz para o tile o que `WindingArcs.lean`
fez para o contorno. A decomposição é por **colunas**, não por linhas, porque a
norma que o certificado usa desde a decisão de N-coh é a `l^1`, que é soma por
coluna: blocos de colunas são o corte natural e dispensam transposição. O lema
de composição é `maxList_append` — o máximo sobre a lista toda é o máximo dos
máximos por bloco —, ligado às normas por `infNorm_eq_maxList` e
`oneNorm_eq_maxList`.

Funciona, e o ganho é qualitativo: a dimensão 40, que **não fechava** de forma
monolítica, fecha decomposta.

| dimensão | colunas por bloco | trabalho por bloco (`n^2·b`) | tempo | por bloco |
|---:|---:|---:|---:|---:|
| 40 | 10 | 16 000 | 316 s | 79 s |
| 60 | 4 | 14 400 | 1280 s | 85 s |

**O tamanho de bloco não é uma constante.** Cada bloco de `b` colunas em
dimensão `n` exige `n^2·b` multiplicações racionais, porque cada entrada do
produto é um produto escalar de comprimento `n`. Para manter o custo por bloco
no patamar que fecha, é preciso `b` proporcional a `1/n^2` — e as duas medições
acima confirmam isso: com `n^2·b` praticamente igual, o tempo por bloco é
praticamente igual (79 s contra 85 s). Com blocos de 10 fixos, as dimensões 60 e
138 estouram.

**A dimensão real não é alcançada, e isso está medido até o piso da abordagem.**

| dimensão | colunas por bloco | resultado |
|---:|---:|---|
| 40 | 10 | fecha, 316 s |
| 60 | 4 | fecha, 1280 s |
| 138 | 2 | falha: 3979 s, 71 estouros (de 69 blocos) |
| 138 | **1** | falha: 7688 s, **101 declarações estouram** |

`b = 1` — uma declaração por coluna — é o **piso absoluto**: não existe bloco
menor. Como ele falha, a decomposição por colunas não chega à dimensão 138, e
nenhum ajuste do tamanho de bloco resolve.

Um detalhe do modo de falha que orienta o próximo passo: com `b = 1` as
primeiras ~37 colunas **fecham** (as falhas começam só na linha 4781 do arquivo
gerado). O custo não é uniforme entre colunas, cresce ao longo do arquivo. Isso
é compatível com o crescimento em bits se acumulando na elaboração, e sugere que
o gargalo não é o trabalho aritmético de uma coluna isolada mas o estado que a
elaboração carrega.

**Consequência para C2.** Verificar o tile em nível kernel na dimensão real exige
**estrutura**, não partição. Partição já foi levada ao limite. As direções que
sobram, nenhuma explorada:

- explorar a estrutura de **produto tensorial** Fourier x Chebyshev da matriz, que
  pode fatorar o produto em fatores muito menores;
- certificar um **problema reduzido** de dimensão menor, em vez do truncamento
  cheio, com um argumento que transfira a conclusão;
- aceitar nível N (`native_decide`) para o volume e reservar o nível K para os
  pontos de fecho — o que custaria o axioma opaco e teria de ser declarado.

### O limite do contorno (obrigação K6)

O mesmo teto vale para o certificado de winding, e agora foi medido sobre
contornos gerados a partir da matriz espectral real (`build/spectrum/rt-A-4x12.dat`,
via `scripts/build_contour.py`):

| contorno | arestas | ajustes padrão | com `maxHeartbeats 0` |
|---|---:|---|---:|
| em torno do candidato físico | 514 | **fecha**, ~2,5 a 3 min | — |
| `Gamma_A` completo | 1946 | **não fecha** | fecha, 673 s |

O custo é essencialmente **linear** no número de arestas (144 s para 514, 673 s
para 1946: razão 4,7 para 3,8 vezes mais arestas). Não há explosão — o problema
é só o orçamento por declaração.

O caso de 1946 arestas falha por `maxHeartbeats` em **dois** pontos distintos, e
a distinção importa:

```
:10:4:     (deterministic) timeout at `«LCNF compiler»`
:11691:47: (deterministic) timeout at `whnf`
```

O segundo é a checagem da proposição, esperado. O primeiro é o **compilador**
tentando gerar código executável para a `def` do certificado — trabalho que não
serve para nada aqui, já que o objeto só existe para ser reduzido pelo kernel.

Consequência de desenho, análoga à de K5 para os tiles: o certificado de winding
precisa ser **decomposto em arcos**. Isso é natural na estrutura, porque o total
é uma soma sobre a lista.

`WindingArcs.lean` faz isso, e a medição sobre o contorno real fecha a questão:

| forma | orçamento | tempo |
|---|---|---:|
| monolítico, 1946 arestas | exige `maxHeartbeats 0` | 673 s |
| 4 arcos de ~487 arestas | **ajustes padrão** | **362 s** |

Decompor não é só mais auditável — é mais barato. O orçamento por declaração
não era apenas um teto arbitrário: a elaboração de uma única declaração enorme
custa mais que a soma das partes. A propriedade que isso preserva é a que dá
sentido ao projeto: um auditor externo reproduz a verificação sem afrouxar
limite nenhum.

### Cuidado com `maxHeartbeats` como remédio

Subir `maxHeartbeats` faz um módulo lento compilar, e por isso mascara a
diferença entre duas causas muito distintas:

- **custo real de kernel**, como o tile 16x16, onde a aritmética exata é
  genuinamente cara e só resta escolher uma dimensão menor;
- **prova mal escrita**, onde o elaborador normaliza um tipo caro sem
  necessidade.

Caso concreto registrado neste repositório: `TailTileBridge.lean` compilava em
**161 s** e o `set_option maxHeartbeats 2000000` "resolvia". Medindo cada
regressão isoladamente, todas custavam menos de 1 s — o gasto inteiro vinha de
um `have` **sem anotação de tipo**. Com a anotação, o módulo passou a **4,4 s** e
o `set_option` foi removido.

A regra prática: antes de subir o orçamento, meça as regressões isoladamente. Se
elas forem baratas e o módulo caro, o problema é a prova, não o kernel.

## Política de dois níveis

- **Nível K (kernel).** `decide`, com `unseal` onde a irredutibilidade atrapalha.
  Base de confiança: o kernel do Lean e seus axiomas padrão. `unseal` **não**
  rebaixa o nível.
- **Nível N (nativo).** `native_decide`. Base de confiança: kernel + compilador
  Lean + runtime compilado + um axioma opaco por declaração. Aceitável para
  varredura em massa, **desde que contado**. Não pode ser descrito como
  "verificado pelo kernel".

## Estado

| ID | Obrigação de confiança | Estado |
|---|---|---|
| K1 | Nenhum `sorry`/`axiom`/`admit` | satisfeita, checada por `check_formal.sh` |
| K2 | Contagem explícita de `native_decide` | satisfeita, `scripts/count_native_decide.py` |
| K3 | Aritmética do certificado redutível no kernel | **satisfeita**: `RationalCertificate.lean` migrado para `List` + `unseal`; `native_decide` = **0** no repositório inteiro |
| K4 | Winding em nível K | **satisfeita**: `WindingInterval.lean`, todas as regressões por `decide` |
| K5 | Tiles pequenos o bastante para o kernel fechar | **medida até o piso; negativa para a dimensão real**. A decomposição por colunas (`TileBlocks.lean`) faz fechar dim 40 e 60, que não fechavam. Dim 138 falha mesmo com uma declaração por coluna, que é o piso. C2 na dimensão real exige **estrutura**, não partição |
| K6 | Certificado de winding dentro do orçamento | **satisfeita**: `WindingArcs.lean` decompõe o contorno em arcos. As 1946 arestas de `Gamma_A`, que não fechavam nos ajustes padrão, fecham como 4 arcos de ~487 sem inflar orçamento — e mais rápido que o monolítico |

## Migração já aplicada

| o que | antes | depois | motivo |
|---|---|---|---|
| regressões de `WindingCertificate` | `native_decide` | `decide` | aritmética inteira pura |
| regressões de `EdgeSignCertificate` | `native_decide` | `decide` | só comparam; literais com `mkRat` |
| `WindingInterval.lean` | — | `decide` em tudo | desenhado sobre comparações |
| `TailBound.lean` | — | `decide` em tudo (19 regressões) | `List`-free, `unseal Rat.*` |
| `DeterminantBound.lean` | — | `decide` em tudo | majorantes de `det` por Jacobi |
| regressões de `Tile`/`BoundaryCertificate` | `native_decide` | `decide` | `RatMatrix.data` migrado de `Array Rat` para `List Rat`, mais `unseal` |

Contagem de usos reais de `native_decide` em `formal/`: **10 -> 0**.

A migração do `Tile` foi verificada contra o estado anterior: as definições
matemáticas (`inverseDefect`, `inverseBound`, `transferMajorant`,
`inequalities`, `shapeCorrect`, `infNorm`, `rowAbsSum`, `dot`,
`BoundaryCertificate.valid`) são idênticas linha a linha — a troca
`Array -> List` exigiu apenas `data.size -> data.length`, `Array.range ->
List.range` e `#[x] -> [x]`. `RatMatrix.entry` não mudou uma letra. Os três
valores exatos das regressões originais (`1/20`, `10/19`, `3/100`) continuam
válidos, agora recomputados pelo kernel.

A contagem é feita por `scripts/count_native_decide.py`, que remove comentários
de linha e de bloco antes de contar — menções em prosa não inflam o número.

## O que o nível K já sustenta

`ChoptuikFormal.lean` tem `CertifiedSpectralObligations`, em que os dois windings
deixaram de ser `Nat` asseverados e passaram a ser certificados intervalares
aceitos. O teorema final

```
theorem exactlyOnePhysicalUnstableModeCertified
    (obligations : CertifiedSpectralObligations physicalMultiplicity) :
    physicalMultiplicity = 1
```

depende de `[propext, Quot.sound]` — axiomas padrão, sem avaliação nativa. As
obrigações analíticas continuam hipóteses; a contagem discreta não é mais uma
delas.
