import Std
import ChoptuikFormal.RationalCertificate

/-!
# Blocos de colunas: um tile grande partido em declarações pequenas

Este módulo é o análogo, para `Tile`, do que `ChoptuikFormal.WindingArcs` é
para o certificado de winding — e existe pela mesma razão medida.

Um `Tile` de dimensão `20` da matriz espectral real fecha por `decide` nos
ajustes padrão em cerca de 72 s. O de dimensão `40` **não** fecha: estoura
`maxHeartbeats` em `whnf`. A falha não é aritmética — todos esses tiles são
válidos em `Fraction` — nem é um teto do kernel: é o **orçamento por
declaração**. Reduzir o comprimento das entradas não muda nada (medido: 67 bits
contra 16 bits dá 153 s contra 154 s em dimensão 40).

## O que se decompõe, e por quê

A norma que o certificado usa é `oneNorm`, a máxima soma de valores absolutos
por **coluna**:

```
oneNorm M = infNorm (transpose M) = max_j (soma_i |M_ij|) .
```

Um máximo sobre todas as colunas é o máximo dos máximos sobre blocos de
colunas. Logo, se as colunas `[0, n)` do resíduo `I - V·F₀` forem partidas em
blocos consecutivos, cada bloco pode ser uma declaração separada e o total é
`qmax` dos parciais. É o análogo de `sumIncrements_append` de `WindingArcs`,
trocando soma por máximo — e é o `theorem blockNorm_split` abaixo.

O ponto decisivo é que **a coluna `j` do resíduo só precisa da coluna `j` de
`F₀`**: `(I - V·F₀)_{·j} = e_j - V·(F₀)_{·j}`. Nenhum bloco calcula o produto
inteiro. A soma dos custos é a mesma do monolítico; o que muda é que ela passa
a ser paga em `k` orçamentos, e não em um só.

## O que a decomposição não faz

Vale dizer explicitamente, porque é fácil vender isto errado: partir em blocos
**não reduz o trabalho aritmético**. O produto `V·F₀` continua custando
`O(n³)` operações exatas em `Rat`, distribuídas entre os blocos. O que a
decomposição compra é **viabilidade** dentro do orçamento padrão — nenhum
`set_option maxHeartbeats` — e nada além disso. Em `WindingArcs` a
decomposição saiu mais barata que o monolítico porque o contorno grande
disparava também o compilador LCNF; aqui o custo é aritmética pura, e ele não
desaparece. Os números medidos estão em `docs/KERNEL_TRUST.md`, obrigação K5.

## Como o resultado final é o mesmo do tile monolítico

Nada aqui enfraquece o que se afirma. `tileValidOfBlockNorms` termina em
`tile.valid = true`, exatamente a proposição que `by decide` no tile inteiro
produziria, e por `decide_eq_true` aplicado a uma prova de `tile.inequalities`
— a mesma `Prop` de `RationalCertificate`, sem variante nem enfraquecimento. O
kernel recomputa cada produto, cada norma e cada desigualdade; só que em
pedaços.

Os dois valores que atravessam a fronteira entre declarações são racionais
exatos: a norma do resíduo e `‖V‖₁`. Não são estimativas, são os valores que o
kernel calculou no bloco correspondente.
-/

namespace Choptuik.TileBlocks

open Choptuik.RationalCertificate

/- `RationalCertificate` mantém `Tile.inequalitiesDecidable` `private`, então,
   como em `WindingArcs`, reconstruo o utilitário que permite provar
   `... .valid = true` sem nomear a instância privada. -/
private theorem decide_true_with {claim : Prop} (inst : Decidable claim)
    (proof : claim) : @decide claim inst = true :=
  @decide_eq_true claim inst proof

/-! ## `qmax` e o máximo de uma lista

`RationalCertificate.absSum` é `private`; a definição abaixo tem o mesmo corpo,
letra por letra, de modo que `infNorm_eq_maxList` vale por `rfl`. -/

/-- Soma dos valores absolutos de uma linha. -/
def absSum (row : List Rat) : Rat :=
  row.foldl (fun total x => total + qabs x) 0

/-- Máximo de uma lista de racionais, a partir de `0` — a mesma dobra que
`infNorm` faz sobre as linhas. -/
def maxList (values : List Rat) : Rat :=
  values.foldl qmax 0

private theorem qmax_of_le {x y : Rat} (le : y ≤ x) : qmax x y = x := by
  unfold qmax
  exact if_neg (Rat.not_lt.mpr le)

private theorem qmax_nonneg {x : Rat} (nonneg : 0 ≤ x) (y : Rat) : 0 ≤ qmax x y := by
  unfold qmax
  split
  · exact Rat.le_trans nonneg (Rat.le_of_lt (by assumption))
  · exact nonneg

private theorem qmax_assoc (x y z : Rat) : qmax (qmax x y) z = qmax x (qmax y z) := by
  by_cases lowMid : x < y
  · by_cases midHigh : y < z
    · have lowHigh : x < z := Std.lt_trans lowMid midHigh
      simp [qmax, lowMid, midHigh, lowHigh]
    · simp [qmax, lowMid, midHigh]
  · by_cases midHigh : y < z
    · simp [qmax, lowMid, midHigh]
    · have midLow : y ≤ x := Rat.not_lt.mp lowMid
      have highMid : z ≤ y := Rat.not_lt.mp midHigh
      have highLow : ¬ x < z := Rat.not_lt.mpr (Rat.le_trans highMid midLow)
      simp [qmax, lowMid, midHigh, highLow]

private theorem foldl_qmax_nonneg : ∀ (values : List Rat) (start : Rat), 0 ≤ start →
    0 ≤ values.foldl qmax start
  | [], _, nonneg => nonneg
  | value :: rest, start, nonneg =>
      foldl_qmax_nonneg rest (qmax start value) (qmax_nonneg nonneg value)

theorem maxList_nonneg (values : List Rat) : 0 ≤ maxList values :=
  foldl_qmax_nonneg values 0 (Rat.le_refl)

private theorem foldl_qmax_start : ∀ (values : List Rat) (start seed : Rat),
    values.foldl qmax (qmax start seed) = qmax start (values.foldl qmax seed)
  | [], _, _ => rfl
  | value :: rest, start, seed => by
      show rest.foldl qmax (qmax (qmax start seed) value)
        = qmax start (rest.foldl qmax (qmax seed value))
      rw [qmax_assoc, foldl_qmax_start rest start (qmax seed value)]

/-- **O lema de composição, na forma de listas.** O máximo de uma concatenação
é o máximo dos máximos. É este fato — e só ele — que permite que cada bloco de
colunas vire uma declaração separada. -/
theorem maxList_append (left right : List Rat) :
    maxList (left ++ right) = qmax (maxList left) (maxList right) := by
  have step := foldl_qmax_start right (maxList left) 0
  rw [qmax_of_le (maxList_nonneg left)] at step
  show (left ++ right).foldl qmax 0 = _
  rw [List.foldl_append]
  exact step

/-! ## Fatos estruturais sobre as listas de `RatMatrix`

Tudo abaixo é aritmética de listas: nenhum `Rat` é somado ou multiplicado.
Serve para relacionar o objeto que o `Tile` define com os blocos que as
declarações separadas calculam. -/

theorem length_chunkRows : ∀ (count width : Nat) (data : List Rat),
    (chunkRows count width data).length = count
  | 0, _, _ => rfl
  | count + 1, width, data => by
      show (chunkRows count width (data.drop width)).length + 1 = count + 1
      rw [length_chunkRows count width (data.drop width)]

theorem length_rowList (matrix : RatMatrix) : matrix.rowList.length = matrix.rows :=
  length_chunkRows _ _ _

theorem chunkRows_rect : ∀ (count width : Nat) (data : List Rat),
    count * width ≤ data.length → ∀ row ∈ chunkRows count width data, row.length = width
  | 0, _, _, _, _, member => absurd member (by simp [chunkRows])
  | count + 1, width, data, bound, row, member => by
      have headBound : width ≤ data.length := by
        have : width ≤ (count + 1) * width := by
          rw [Nat.succ_mul]; exact Nat.le_add_left _ _
        exact Nat.le_trans this bound
      have tailBound : count * width ≤ (data.drop width).length := by
        rw [List.length_drop]
        rw [Nat.succ_mul] at bound
        omega
      rcases List.mem_cons.mp member with head | tail
      · subst head
        rw [List.length_take]
        omega
      · exact chunkRows_rect count width (data.drop width) tailBound row tail

theorem rect_rowList (matrix : RatMatrix) (wf : matrix.wellFormed) :
    ∀ row ∈ matrix.rowList, row.length = matrix.cols :=
  chunkRows_rect matrix.rows matrix.cols matrix.data (Nat.le_of_eq wf.symm)

theorem length_concatRows : ∀ (rows : List (List Rat)) (width : Nat),
    (∀ row ∈ rows, row.length = width) → (concatRows rows).length = rows.length * width
  | [], _, _ => by simp [concatRows]
  | row :: rest, width, rect => by
      show (row ++ concatRows rest).length = (rest.length + 1) * width
      rw [List.length_append, rect row (by simp),
          length_concatRows rest width (fun r member => rect r (by simp [member])),
          Nat.succ_mul]
      omega

/-- Quebrar em linhas o que foi concatenado a partir de linhas devolve as
linhas. É o lema que liga `mul` e `transpose`, que produzem `data` por
`concatRows`, de volta à lista de linhas. -/
theorem chunkRows_concatRows : ∀ (rows : List (List Rat)) (width : Nat),
    (∀ row ∈ rows, row.length = width) →
    chunkRows rows.length width (concatRows rows) = rows
  | [], _, _ => rfl
  | row :: rest, width, rect => by
      have rowLength : row.length = width := rect row (by simp)
      have restRect : ∀ r ∈ rest, r.length = width := fun r member => rect r (by simp [member])
      have takeRow : row.take width = row := List.take_of_length_le (Nat.le_of_eq rowLength)
      have dropRow : row.drop width = [] := List.drop_of_length_le (Nat.le_of_eq rowLength)
      show (row ++ concatRows rest).take width
          :: chunkRows rest.length width ((row ++ concatRows rest).drop width) = _
      rw [List.take_append, List.drop_append, rowLength, Nat.sub_self, takeRow, dropRow,
          List.take_zero, List.drop_zero, List.append_nil, List.nil_append,
          chunkRows_concatRows rest width restRect]

/-! ### `transposeRows` -/

theorem length_transposeRows_le : ∀ (width : Nat) (rows : List (List Rat)),
    (transposeRows width rows).length ≤ width
  | width, [] => by simp [transposeRows]
  | width, row :: rest => by
      show (List.zipWith (fun x column => x :: column) row (transposeRows width rest)).length ≤ width
      rw [List.length_zipWith]
      exact Nat.le_trans (Nat.min_le_right _ _) (length_transposeRows_le width rest)

theorem length_transposeRows : ∀ (width : Nat) (rows : List (List Rat)),
    (∀ row ∈ rows, row.length = width) → (transposeRows width rows).length = width
  | width, [], _ => by simp [transposeRows]
  | width, row :: rest, rect => by
      show (List.zipWith (fun x column => x :: column) row (transposeRows width rest)).length = width
      rw [List.length_zipWith, rect row (by simp),
          length_transposeRows width rest (fun r member => rect r (by simp [member]))]
      omega

private theorem mem_zipWith_cons : ∀ {row : List Rat} {columns : List (List Rat)}
    {item : List Rat}, item ∈ List.zipWith (fun x column => x :: column) row columns →
    ∃ head tail, item = head :: tail ∧ tail ∈ columns
  | _ :: _, [], _, member => absurd member (by simp)
  | [], _, _, member => absurd member (by simp)
  | x :: row, column :: columns, item, member => by
      rcases List.mem_cons.mp member with head | tail
      · exact ⟨x, column, head, by simp⟩
      · obtain ⟨h, t, itemEq, memTail⟩ := mem_zipWith_cons tail
        exact ⟨h, t, itemEq, by simp [memTail]⟩

theorem rect_transposeRows : ∀ (width : Nat) (rows : List (List Rat)),
    (∀ row ∈ rows, row.length = width) →
    ∀ column ∈ transposeRows width rows, column.length = rows.length
  | width, [], _, column, member => by
      have : column = [] := List.eq_of_mem_replicate member
      simp [this]
  | width, row :: rest, rect, column, member => by
      obtain ⟨head, tail, columnEq, memTail⟩ := mem_zipWith_cons member
      have tailLength : tail.length = rest.length :=
        rect_transposeRows width rest (fun r m => rect r (by simp [m])) tail memTail
      subst columnEq
      simp [tailLength]

theorem length_colList (matrix : RatMatrix) (wf : matrix.wellFormed) :
    matrix.colList.length = matrix.cols :=
  length_transposeRows matrix.cols matrix.rowList (rect_rowList matrix wf)

theorem rect_colList (matrix : RatMatrix) (wf : matrix.wellFormed) :
    ∀ column ∈ matrix.colList, column.length = matrix.rows := by
  intro column member
  have := rect_transposeRows matrix.cols matrix.rowList (rect_rowList matrix wf) column member
  rw [this, length_rowList]

/-- As linhas da transposta são as colunas. -/
theorem rowList_transpose (matrix : RatMatrix) (wf : matrix.wellFormed) :
    matrix.transpose.rowList = matrix.colList := by
  show chunkRows matrix.cols matrix.rows (concatRows matrix.colList) = matrix.colList
  rw [← length_colList matrix wf]
  exact chunkRows_concatRows matrix.colList matrix.rows (rect_colList matrix wf)

/-! ## As normas em termos de `maxList` -/

theorem infNorm_eq_maxList (matrix : RatMatrix) :
    matrix.infNorm = maxList (matrix.rowList.map absSum) := by
  unfold maxList
  rw [List.foldl_map]
  rfl

theorem oneNorm_eq_maxList (matrix : RatMatrix) (wf : matrix.wellFormed) :
    matrix.oneNorm = maxList (matrix.colList.map absSum) := by
  show matrix.transpose.infNorm = _
  rw [infNorm_eq_maxList, rowList_transpose matrix wf]

/-! ## Blocos de colunas

`blockNorm width rows start count` é a `oneNorm` restrita às colunas
`[start, start + count)` da matriz cujas linhas são `rows`. Note que `rows`
pode ser — e no uso real é — uma lista **não avaliada**: o `drop`/`take` sobre
a transposta força apenas o esqueleto das listas, de modo que só as colunas
selecionadas têm as suas entradas de fato calculadas. É isso que faz de cada
bloco uma fração do custo, e não uma cópia dele. -/

/-- As colunas `[start, start + count)`. -/
def colBlock (width : Nat) (rows : List (List Rat)) (start count : Nat) : List (List Rat) :=
  ((transposeRows width rows).drop start).take count

/-- A norma `l¹` do bloco de colunas: máximo, sobre as colunas do bloco, da
soma de valores absolutos. -/
def blockNorm (width : Nat) (rows : List (List Rat)) (start count : Nat) : Rat :=
  maxList ((colBlock width rows start count).map absSum)

/-- **O lema de composição.** Um bloco de `a + b` colunas é o `qmax` do bloco
das primeiras `a` com o bloco das `b` seguintes. Aplicado `k-1` vezes, ele
troca uma declaração impossível por `k` declarações possíveis mais álgebra
grátis: a recomposição não reexamina entrada nenhuma.

Os parâmetros `count` e `next` são redundantes (`count = a + b`,
`next = start + a`) e estão aqui para que o gerador possa escrever os numerais
já somados, casando sintaticamente com os lemas dos blocos. -/
theorem blockNorm_split (width : Nat) (rows : List (List Rat))
    (start count a b next : Nat) (countEq : count = a + b) (nextEq : next = start + a) :
    blockNorm width rows start count
      = qmax (blockNorm width rows start a) (blockNorm width rows next b) := by
  subst countEq
  subst nextEq
  unfold blockNorm colBlock
  rw [List.take_add, List.drop_drop, List.map_append, maxList_append]

/-- O bloco que cobre todas as colunas é a norma inteira. -/
theorem oneNorm_eq_blockNorm (matrix : RatMatrix) (wf : matrix.wellFormed) :
    matrix.oneNorm = blockNorm matrix.cols matrix.rowList 0 matrix.cols := by
  rw [oneNorm_eq_maxList matrix wf]
  unfold blockNorm colBlock
  show _ = maxList ((((matrix.colList).drop 0).take matrix.cols).map absSum)
  rw [List.drop_zero, List.take_of_length_le (Nat.le_of_eq (length_colList matrix wf))]

/-! ## O resíduo `I - V·F₀` linha a linha

`Tile.inverseDefect` calcula `((identity n).sub (V.mul F₀)).oneNorm`. A
definição abaixo descreve as **linhas** desse resíduo diretamente, sem passar
pelo layout row-major, e `residualMatrix_rowList` prova que as duas coincidem.
A vantagem é que `residualRows` é uma `zipWith`: fatiar as suas colunas força
apenas o esqueleto, nunca as entradas descartadas. -/

/-- Linha `i` do resíduo: `e_i - (linha i de V) · F₀`, com `cols` as colunas de
`F₀`. -/
def residualRows (idRows vRows cols : List (List Rat)) : List (List Rat) :=
  List.zipWith
    (fun idRow vRow => List.zipWith (fun a b => a - b) idRow (cols.map (dotRow vRow)))
    idRows vRows

/-- O resíduo do tile, exatamente a matriz que `Tile.inverseDefect` normaliza. -/
def residualMatrix (tile : Tile) : RatMatrix :=
  (RatMatrix.identity tile.centerMatrix.rows).sub
    (tile.approxInverse.mul tile.centerMatrix)

theorem inverseDefect_eq_residual (tile : Tile) :
    tile.inverseDefect
      = (residualMatrix tile).oneNorm + tile.approxInverse.oneNorm * tile.matrixVariation :=
  rfl

theorem chunkRows_zipWith (f : Rat → Rat → Rat) : ∀ (count width : Nat) (left right : List Rat),
    chunkRows count width (List.zipWith f left right)
      = List.zipWith (List.zipWith f) (chunkRows count width left) (chunkRows count width right)
  | 0, _, _, _ => rfl
  | count + 1, width, left, right => by
      show (List.zipWith f left right).take width
          :: chunkRows count width ((List.zipWith f left right).drop width) = _
      rw [List.take_zipWith, List.drop_zipWith,
          chunkRows_zipWith f count width (left.drop width) (right.drop width)]
      rfl

private theorem chunkRows_congr {count₁ count₂ width₁ width₂ : Nat} (countEq : count₁ = count₂)
    (widthEq : width₁ = width₂) (data : List Rat) :
    chunkRows count₁ width₁ data = chunkRows count₂ width₂ data := by
  rw [countEq, widthEq]

/-- As linhas de um produto: linha `i` do produto é a linha `i` da esquerda
contra todas as colunas da direita. Nenhuma outra linha da esquerda participa —
é o que torna a fatia de linhas barata, e por transposição a fatia de colunas. -/
theorem rowList_mul (left right : RatMatrix) (rightWf : right.wellFormed) :
    (left.mul right).rowList
      = left.rowList.map (fun row => right.colList.map (dotRow row)) := by
  have rect : ∀ row ∈ left.rowList.map (fun row => right.colList.map (dotRow row)),
      row.length = right.cols := by
    intro row member
    obtain ⟨source, _, rowEq⟩ := List.mem_map.mp member
    subst rowEq
    rw [List.length_map, length_colList right rightWf]
  have counted : (left.rowList.map (fun row => right.colList.map (dotRow row))).length
      = left.rows := by
    rw [List.length_map, length_rowList]
  show chunkRows left.rows right.cols
      (concatRows (left.rowList.map (fun row => right.colList.map (dotRow row)))) = _
  rw [← counted]
  exact chunkRows_concatRows _ right.cols rect

/-- **As linhas do resíduo são `residualRows`.** Aqui a matriz que o `Tile`
define e o objeto que os blocos calculam passam a ser o mesmo. -/
theorem residualMatrix_rowList (tile : Tile) (shape : tile.shapeCorrect) :
    (residualMatrix tile).rowList
      = residualRows (RatMatrix.identity tile.centerMatrix.rows).rowList
          tile.approxInverse.rowList tile.centerMatrix.colList := by
  obtain ⟨centerWf, _, square, inverseRows, _, _⟩ := shape
  have productRows : tile.approxInverse.rows = tile.centerMatrix.rows := inverseRows
  have productCols : tile.centerMatrix.cols = tile.centerMatrix.rows := square.symm
  show chunkRows tile.centerMatrix.rows tile.centerMatrix.rows
      (List.zipWith (fun x y => x - y)
        (RatMatrix.identity tile.centerMatrix.rows).data
        (tile.approxInverse.mul tile.centerMatrix).data) = _
  rw [chunkRows_zipWith]
  rw [chunkRows_congr productRows.symm productCols.symm
        (tile.approxInverse.mul tile.centerMatrix).data]
  show List.zipWith (List.zipWith (fun x y => x - y))
      (RatMatrix.identity tile.centerMatrix.rows).rowList
      (tile.approxInverse.mul tile.centerMatrix).rowList = _
  rw [rowList_mul tile.approxInverse tile.centerMatrix centerWf, List.zipWith_map_right]
  rfl

theorem residualMatrix_wellFormed (tile : Tile) (shape : tile.shapeCorrect) :
    (residualMatrix tile).wellFormed := by
  obtain ⟨centerWf, _, square, inverseRows, _, _⟩ := shape
  have identityLength :
      (RatMatrix.identity tile.centerMatrix.rows).data.length
        = tile.centerMatrix.rows * tile.centerMatrix.rows := by
    show ((List.range (tile.centerMatrix.rows * tile.centerMatrix.rows)).map _).length = _
    rw [List.length_map, List.length_range]
  have productLength : (tile.approxInverse.mul tile.centerMatrix).data.length
      = tile.centerMatrix.rows * tile.centerMatrix.rows := by
    show (concatRows (tile.approxInverse.rowList.map
      (fun row => tile.centerMatrix.colList.map (dotRow row)))).length = _
    rw [length_concatRows _ tile.centerMatrix.cols, List.length_map, length_rowList,
        inverseRows, square]
    intro row member
    obtain ⟨_, _, rowEq⟩ := List.mem_map.mp member
    subst rowEq
    rw [List.length_map, length_colList tile.centerMatrix centerWf]
  show (List.zipWith (fun x y => x - y)
      (RatMatrix.identity tile.centerMatrix.rows).data
      (tile.approxInverse.mul tile.centerMatrix).data).length
    = tile.centerMatrix.rows * tile.centerMatrix.rows
  rw [List.length_zipWith, identityLength, productLength]
  simp

/-! ## O bloco de um tile -/

/-- As linhas do resíduo do tile, na forma que os blocos fatiam. -/
def residualRowList (tile : Tile) : List (List Rat) :=
  residualRows (RatMatrix.identity tile.centerMatrix.rows).rowList
    tile.approxInverse.rowList tile.centerMatrix.colList

/-- **A declaração cara, uma por bloco.** A norma `l¹` do resíduo restrita às
colunas `[start, start + count)`. -/
def residualBlockNorm (tile : Tile) (start count : Nat) : Rat :=
  blockNorm tile.centerMatrix.rows (residualRowList tile) start count

/-- A composição dos blocos do tile: aplicação direta de `blockNorm_split`,
sem nenhuma hipótese sobre o tile, porque é álgebra de listas. -/
theorem residualBlockNorm_split (tile : Tile) (start count a b next : Nat)
    (countEq : count = a + b) (nextEq : next = start + a) :
    residualBlockNorm tile start count
      = qmax (residualBlockNorm tile start a) (residualBlockNorm tile next b) :=
  blockNorm_split _ _ start count a b next countEq nextEq

/-- O bloco que cobre todas as colunas é a norma do resíduo inteiro — a mesma
que `Tile.inverseDefect` usa. -/
theorem oneNorm_residualMatrix (tile : Tile) (shape : tile.shapeCorrect) :
    (residualMatrix tile).oneNorm = residualBlockNorm tile 0 tile.centerMatrix.rows := by
  rw [oneNorm_eq_blockNorm (residualMatrix tile) (residualMatrix_wellFormed tile shape)]
  show blockNorm tile.centerMatrix.rows (residualMatrix tile).rowList 0 tile.centerMatrix.rows = _
  rw [residualMatrix_rowList tile shape]
  rfl

/-! ## Do bloco ao tile

As três declarações abaixo são o que o arquivo gerado usa. Cada uma consome
racionais **exatos** vindos das declarações dos blocos; nenhuma delas recalcula
o produto de matrizes. -/

/-- `inverseDefect` a partir da norma do resíduo e de `‖V‖₁`. -/
theorem inverseDefect_eq (tile : Tile) (shape : tile.shapeCorrect)
    (residual inverse : Rat)
    (residualEq : residualBlockNorm tile 0 tile.centerMatrix.rows = residual)
    (inverseEq : tile.approxInverse.oneNorm = inverse) :
    tile.inverseDefect = residual + inverse * tile.matrixVariation := by
  rw [inverseDefect_eq_residual tile, oneNorm_residualMatrix tile shape, residualEq, inverseEq]

/-- `transferMajorant` a partir de `inverseDefect` e de `‖V‖₁`. -/
theorem transferMajorant_eq (tile : Tile) (defect inverse : Rat)
    (defectEq : tile.inverseDefect = defect)
    (inverseEq : tile.approxInverse.oneNorm = inverse) :
    tile.transferMajorant
      = tile.maxDistance * tile.tailError *
          (1 + tile.maxDistance * (inverse / (1 - defect)) * tile.operatorBound) := by
  show tile.maxDistance * tile.tailError *
      (1 + tile.maxDistance * (tile.approxInverse.oneNorm / (1 - tile.inverseDefect))
        * tile.operatorBound) = _
  rw [defectEq, inverseEq]

/-- `Tile.shapeCorrect` é decidível, mas a instância não se sintetiza sozinha:
`RatMatrix.wellFormed` é uma `def` semirredutível e a síntese não a desdobra.
Com esta instância, as obrigações de forma do arquivo gerado fecham por
`decide` — custo `O(n²)` de percurso de lista, sem uma única operação
aritmética. -/
instance decidableShapeCorrect (tile : Tile) : Decidable tile.shapeCorrect := by
  unfold Tile.shapeCorrect RatMatrix.wellFormed
  infer_instance

/-- **O tile aceito por partes.** Mesma conclusão que `by decide` no tile
inteiro produziria, e pela mesma `Prop`: `Tile.inequalities`. -/
theorem tileValidOfParts (tile : Tile)
    (shape : tile.shapeCorrect)
    (radius : 0 ≤ tile.tileRadius) (matrixNorm : 0 ≤ tile.matrixNorm)
    (distance : 0 ≤ tile.maxDistance) (tailError : 0 ≤ tile.tailError)
    (operator : 0 ≤ tile.operatorBound)
    (defect : tile.inverseDefect < 1)
    (transfer : tile.transferMajorant < 1) :
    tile.valid = true :=
  decide_true_with _
    ⟨shape, radius, matrixNorm, distance, tailError, operator, defect, transfer⟩

/-- **O certificado de tile por blocos.** `residual` é o `qmax` das normas dos
blocos e `inverse` é `‖V‖₁`; as duas desigualdades finais são aritmética sobre
racionais literais, e fecham por `decide` em tempo desprezível. O resultado é
literalmente `tile.valid = true`. -/
theorem tileValidOfBlockNorms (tile : Tile) (residual inverse : Rat)
    (shape : tile.shapeCorrect)
    (residualEq : residualBlockNorm tile 0 tile.centerMatrix.rows = residual)
    (inverseEq : tile.approxInverse.oneNorm = inverse)
    (radius : 0 ≤ tile.tileRadius) (matrixNorm : 0 ≤ tile.matrixNorm)
    (distance : 0 ≤ tile.maxDistance) (tailError : 0 ≤ tile.tailError)
    (operator : 0 ≤ tile.operatorBound)
    (defect : residual + inverse * tile.matrixVariation < 1)
    (transfer : tile.maxDistance * tile.tailError *
        (1 + tile.maxDistance
          * (inverse / (1 - (residual + inverse * tile.matrixVariation)))
          * tile.operatorBound) < 1) :
    tile.valid = true := by
  have defectEq : tile.inverseDefect = residual + inverse * tile.matrixVariation :=
    inverseDefect_eq tile shape residual inverse residualEq inverseEq
  refine tileValidOfParts tile shape radius matrixNorm distance tailError operator ?_ ?_
  · rw [defectEq]; exact defect
  · rw [transferMajorant_eq tile _ inverse defectEq inverseEq]; exact transfer

/-! ## Regressões

`unseal` abaixo remove o atributo `@[irreducible]` de quatro funções de `Rat`.
**Isso não aumenta a base de confiança**: irredutibilidade é diretiva de
*elaboração*, e o kernel nunca a respeitou. Ver `docs/KERNEL_TRUST.md`.

As regressões usam `tile4` e `tile8` de `RationalCertificate`, que fecham
monoliticamente por `decide` — é justamente por isso que servem: aqui as duas
rotas são exercitadas lado a lado e o kernel exige que deem o **mesmo**
resultado. Nos tamanhos em que só a rota dos blocos fecha não existe termo de
comparação, e a regressão perderia o sentido. -/

unseal Rat.mul Rat.add Rat.sub Rat.inv

set_option maxRecDepth 100000

/-! ### `tile4`: dois blocos de duas colunas -/

theorem tile4Shape : tile4.shapeCorrect := by decide

theorem tile4BlockLow : residualBlockNorm tile4 0 2 = 265 / 16384 := by decide
theorem tile4BlockHigh : residualBlockNorm tile4 2 2 = 15 / 1024 := by decide

theorem tile4Residual : residualBlockNorm tile4 0 4 = 265 / 16384 := by
  rw [residualBlockNorm_split tile4 0 4 2 2 2 rfl rfl, tile4BlockLow, tile4BlockHigh]
  decide

theorem tile4InverseNorm : tile4.approxInverse.oneNorm = 35 / 32 := by decide

/-- O valor exato que a regressão monolítica de `RationalCertificate` afirma
(`545/16384`), reobtido pela rota dos blocos. -/
theorem tile4DefectByBlocks : tile4.inverseDefect = 545 / 16384 := by
  rw [inverseDefect_eq tile4 tile4Shape _ _ tile4Residual tile4InverseNorm]
  decide

-- E o monolítico, recomputado aqui pelo kernel, dá o mesmo número.
example : tile4.inverseDefect = 545 / 16384 := by decide

/-! ### `tile8`: dois blocos de quatro colunas, e quatro de duas

O bloco baixo sozinho vale `499/8192` e o alto vale `553/8192`: o máximo está
no segundo bloco, de modo que a composição não é decorativa — sem ela o
certificado afirmaria um número menor que o verdadeiro. -/

theorem tile8Shape : tile8.shapeCorrect := by decide

theorem tile8BlockLow : residualBlockNorm tile8 0 4 = 499 / 8192 := by decide
theorem tile8BlockHigh : residualBlockNorm tile8 4 4 = 553 / 8192 := by decide

example : residualBlockNorm tile8 0 4 ≠ residualBlockNorm tile8 0 8 := by
  rw [tile8BlockLow]
  decide

theorem tile8Residual : residualBlockNorm tile8 0 8 = 553 / 8192 := by
  rw [residualBlockNorm_split tile8 0 8 4 4 4 rfl rfl, tile8BlockLow, tile8BlockHigh]
  decide

theorem tile8InverseNorm : tile8.approxInverse.oneNorm = 157 / 128 := by decide

/-- A mesma norma do resíduo, agora com quatro blocos de duas colunas: a
composição é associativa e o resultado não depende de como se corta. -/
theorem tile8ResidualByFour : residualBlockNorm tile8 0 8 = 553 / 8192 := by
  rw [residualBlockNorm_split tile8 0 8 2 6 2 rfl rfl,
      residualBlockNorm_split tile8 2 6 2 4 4 rfl rfl,
      residualBlockNorm_split tile8 4 4 2 2 6 rfl rfl]
  rw [show residualBlockNorm tile8 0 2 = 961 / 16384 by decide,
      show residualBlockNorm tile8 2 2 = 499 / 8192 by decide,
      show residualBlockNorm tile8 4 2 = 1063 / 16384 by decide,
      show residualBlockNorm tile8 6 2 = 553 / 8192 by decide]
  decide

/-- `inverseDefect` de `tile8` por blocos: `355/4096`, o valor que a regressão
monolítica de `RationalCertificate` afirma. -/
theorem tile8DefectByBlocks : tile8.inverseDefect = 355 / 4096 := by
  rw [inverseDefect_eq tile8 tile8Shape _ _ tile8Residual tile8InverseNorm]
  decide

example : tile8.inverseDefect = 355 / 4096 := by decide

/-- **A regressão que importa.** `tile8.valid = true` obtido sem nenhum `decide`
sobre o produto de matrizes inteiro: só os dois blocos, a composição por `qmax`
e aritmética sobre literais. A proposição é exatamente a mesma que
`example : tile8.valid = true := by decide` afirma em `RationalCertificate`. -/
theorem tile8ValidByBlocks : tile8.valid = true :=
  tileValidOfBlockNorms tile8 (553 / 8192) (157 / 128) tile8Shape tile8Residual
    tile8InverseNorm (by decide) (by decide) (by decide) (by decide) (by decide)
    (by decide) (by decide)

-- O monolítico, para comparação, recomputado pelo kernel neste mesmo módulo.
example : tile8.valid = true := by decide

/-- Um tile recusado continua recusado: a rota dos blocos não pode produzir
`valid = true` quando o defeito não cabe. Aqui `‖I - V F‖₁ = 0` e todo o
defeito vem da variação anunciada, que é grande demais. -/
example : residualBlockNorm rejectedDefectTile 0 1 = 0 := by decide

example : rejectedDefectTile.valid = false := by decide

#print axioms tile8ValidByBlocks
#print axioms tileValidOfBlockNorms
#print axioms blockNorm_split

end Choptuik.TileBlocks
