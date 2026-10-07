import Std
import Lean.Elab.Tactic.Decide

/-!
Exact checker for the rational inequalities in a Neumann/Fredholm boundary
certificate.

This module deliberately has no floating-point input.  A numerical program may
suggest the matrices and bounds, but Lean recomputes every matrix product, norm
and strict inequality with `Rat`.

The checker does *not* manufacture the analytic estimate `tailError`: a
complete Choptuik proof must derive it from the PDE.  What is proved here is
that, once that upper bound is supplied, the finite Neumann inequalities used on
every contour tile really hold.

Nota: `matrixVariation` saiu dessa lista e **nao** e mais um campo. Para o
pencil afim `F_N(z) = I + (z-z0) C_N` vale `F_N(z) - F_N(z_j) = (z - z_j) C_N`
exatamente, logo a variacao no tile e `raio_do_tile * norma(C_N)` e aritmetica
sobre a matriz racional, nao estimativa. Agora ela e um `def` derivado de tres
dados -- `tileRadius`, `matrixNorm` e a norma de `centerResidual` --, e nao um
numero que o produtor entrega pronto. Ver `docs/WINDING_DATA.md` secao 3.6.

O terceiro ingrediente permite que `centerMatrix` seja um **arredondado** de
`F_N(z_j)`: pela desigualdade triangular basta somar `||F_N(z_j) - centro||` a
variacao, e esse termo e a norma de uma matriz escrita no certificado,
recomputada pelo kernel, nao um numero estimado. Ver o campo `centerResidual`.

## Nivel K: tudo aqui e recomputado pelo kernel

Todas as regressoes deste arquivo fecham por `decide`, isto e, sao recomputadas
pelo kernel do Lean. Nao ha avaliacao nativa em lugar nenhum, e portanto nenhum
axioma opaco gerado por declaracao. Duas medidas foram necessarias para chegar
a isso, ambas documentadas em `docs/KERNEL_TRUST.md`:

* **`RatMatrix.data` e uma `List Rat`, nao um `Array Rat`.** `Array.map` e
  `Array.range` nao reduzem no kernel -- nem o caso trivial `arr[3]?.getD 0 = 9`
  fecha por `decide`. Como `mul`, `sub` e `identity` passam por `map`, o
  certificado inteiro ficava preso. `List` reduz. A matematica e identica,
  entrada por entrada.
* **`unseal Rat.mul Rat.add Rat.sub Rat.inv`** antes das regressoes. Ver o
  comentario na secao de regressoes: e diretiva de elaboracao e nao mexe na
  base de confianca.

Uma terceira medida veio depois, e e de custo, nao de possibilidade: as
operacoes de matriz deixaram de **indexar** a lista (`RatMatrix.entry`) e
passaram a percorre-la estruturalmente. O produto caiu de `O(n^5)` passos de
reducao para `O(n^3)`, e com isso um tile de dimensao `20` da matriz espectral
real, que antes era morto pelo sistema depois de mais de uma hora, passou a
fechar em `73 s` no orcamento **padrao**. Ver a secao
"Representacao estrutural" logo abaixo.
-/

namespace Choptuik.RationalCertificate

/-- Absolute value on exact rationals. -/
def qabs (x : Rat) : Rat :=
  if x < 0 then -x else x

/-- Maximum on exact rationals. -/
def qmax (x y : Rat) : Rat :=
  if x < y then y else x

/-- A row-major rational matrix.  `wellFormed` rejects malformed external
certificates before any entry is used. -/
structure RatMatrix where
  rows : Nat
  cols : Nat
  data : List Rat
deriving Repr, DecidableEq

def RatMatrix.wellFormed (matrix : RatMatrix) : Prop :=
  matrix.data.length = matrix.rows * matrix.cols

def RatMatrix.entry (matrix : RatMatrix) (row col : Nat) : Rat :=
  matrix.data[row * matrix.cols + col]?.getD 0

/-! ### Representacao estrutural: `entry` saiu do caminho quente

`RatMatrix.entry` indexa uma `List`, e indexar uma lista custa `O(indice)`
passos de reducao no kernel -- nao existe acesso aleatorio a percorrer. A
versao anterior de `mul` chamava `entry` duas vezes por parcela de cada produto
interno, o que da `O(n^5)` passos de lista num produto `n x n`. Esse custo nao
aparece na contagem de operacoes aritmeticas e era o que dominava o tempo de
kernel: medido sobre a matriz espectral real, um tile `8x8` com entradas de 73
bits custa o **mesmo** que um tile `8x8` com entradas de 12 bits (13,3 s contra
14,3 s). Nao eram os bits.

As operacoes abaixo percorrem as listas **estruturalmente** (`chunkRows`,
`zipWith`, `foldl`, `map`), sem indexar: `O(n^3)` no produto e `O(n^2)` nas
demais. A matematica e a mesma, entrada por entrada, para matrizes
`wellFormed`. Isso nao fica por conta de revisao visual: as definicoes
indexadas antigas continuam no arquivo com o sufixo `Spec`, e as regressoes
`*Spec_agrees*` exigem que o **kernel** recompute as duas e as ache iguais. -/

/-- Quebra `data` em `count` blocos de `width` entradas. A recursao e
estrutural no contador: nao ha obrigacao de terminacao, e o resultado tem
exatamente `count` linhas. -/
def chunkRows : Nat → Nat → List Rat → List (List Rat)
  | 0, _, _ => []
  | count + 1, width, data =>
      data.take width :: chunkRows count width (data.drop width)

/-- As linhas da matriz. Para uma matriz `wellFormed` esta e exatamente a
decomposicao de `data` em `rows` linhas de `cols` entradas. -/
def RatMatrix.rowList (matrix : RatMatrix) : List (List Rat) :=
  chunkRows matrix.rows matrix.cols matrix.data

/-- Transposicao de uma lista de linhas, dada a largura: `foldr` com
`zipWith (· :: ·)` sobre `List.replicate width []`. -/
def transposeRows (width : Nat) (rows : List (List Rat)) : List (List Rat) :=
  rows.foldr
    (fun row columns => List.zipWith (fun x column => x :: column) row columns)
    (List.replicate width [])

/-- As colunas da matriz. -/
def RatMatrix.colList (matrix : RatMatrix) : List (List Rat) :=
  transposeRows matrix.cols matrix.rowList

/-- Concatenacao de linhas de volta ao layout row-major. -/
def concatRows : List (List Rat) → List Rat
  | [] => []
  | row :: rest => row ++ concatRows rest

/-- Produto interno de duas linhas, entrada a entrada. -/
def dotRow : List Rat → List Rat → Rat
  | x :: xs, y :: ys => x * y + dotRow xs ys
  | _, _ => 0

/-- Exact dense matrix multiplication: linha `i` de `left` contra coluna `j` de
`right`, em ordem row-major. -/
def RatMatrix.mul (left right : RatMatrix) : RatMatrix where
  rows := left.rows
  cols := right.cols
  data := concatRows (left.rowList.map fun row => right.colList.map (dotRow row))

/-- Exact matrix subtraction.  Shape compatibility is checked separately by
the certificate predicate; `zipWith` trunca no menor, de modo que uma matriz
mal formada nunca fabrica entrada. -/
def RatMatrix.sub (left right : RatMatrix) : RatMatrix where
  rows := left.rows
  cols := left.cols
  data := List.zipWith (fun x y => x - y) left.data right.data

def RatMatrix.identity (size : Nat) : RatMatrix where
  rows := size
  cols := size
  data := (List.range (size * size)).map fun index =>
    if index / size = index % size then 1 else 0

private def absSum (row : List Rat) : Rat :=
  row.foldl (fun total x => total + qabs x) 0

/-- Induced infinity norm, computed exactly as the maximum absolute row sum. -/
def RatMatrix.infNorm (matrix : RatMatrix) : Rat :=
  matrix.rowList.foldl (fun largest row => qmax largest (absSum row)) 0

/-- Soma exata. Nao existia no arquivo original; serve so para montar as
matrizes das regressoes maiores (`I + E`) sem escrever `n^2` literais. -/
def RatMatrix.add (left right : RatMatrix) : RatMatrix where
  rows := left.rows
  cols := left.cols
  data := List.zipWith (fun x y => x + y) left.data right.data

/-- Transposta. Serve para exprimir a norma por coluna sem duplicar codigo. -/
def RatMatrix.transpose (matrix : RatMatrix) : RatMatrix where
  rows := matrix.cols
  cols := matrix.rows
  data := concatRows matrix.colList

/-- Norma induzida por `l^1`: maxima soma de valores absolutos por COLUNA.

**Esta e a norma usada pelo certificado** (`inverseDefect` e `inverseBound`),
por decisao analitica registrada em `docs/WINDING_DATA.md`: `l^1` nos
coeficientes controla a funcao, `l^infty` nao controla nada, e `l^1` e algebra
de Banach com constante `1` -- que e exatamente o que o termo bilinear e a
multiplicacao pelo fundo exigem. `infNorm` (por LINHA) continua no arquivo
porque `oneNorm` e definida atraves dela, e porque o contraste em
`asymmetricMatrix` mostra que as duas sao mesmo diferentes. -/
def RatMatrix.oneNorm (matrix : RatMatrix) : Rat :=
  matrix.transpose.infNorm

/-- A norma por coluna e literalmente a norma por linha da transposta. Vale por
`rfl`: nao ha nada a conferir por inspecao. -/
theorem RatMatrix.oneNorm_eq_transpose_infNorm (matrix : RatMatrix) :
    matrix.oneNorm = matrix.transpose.infNorm := rfl

/-! ### Nao negatividade das normas

`infNorm` e um `foldl` de `qmax` a partir de `0`, logo nunca desce abaixo de
zero. Isso deixa de ser hipotese e vira teorema -- e e o que sustenta
`Tile.matrixVariation_nonneg` depois que a variacao passou a incluir o erro de
arredondamento do centro. -/

private theorem foldl_qmax_nonneg {α : Type} (f : α → Rat) :
    ∀ (items : List α) (start : Rat), 0 ≤ start →
      0 ≤ items.foldl (fun largest item => qmax largest (f item)) start
  | [], _, nonneg => nonneg
  | item :: rest, start, nonneg =>
      foldl_qmax_nonneg f rest (qmax start (f item)) <| by
        show 0 ≤ qmax start (f item)
        unfold qmax
        split
        · exact Rat.le_trans nonneg (Rat.le_of_lt (by assumption))
        · exact nonneg

theorem RatMatrix.infNorm_nonneg (matrix : RatMatrix) : 0 ≤ matrix.infNorm :=
  foldl_qmax_nonneg _ _ _ (Rat.le_refl)

theorem RatMatrix.oneNorm_nonneg (matrix : RatMatrix) : 0 ≤ matrix.oneNorm :=
  matrix.transpose.infNorm_nonneg

/-! ### As definicoes indexadas antigas, mantidas como especificacao

Estas sao exatamente as definicoes que este arquivo usava antes da troca para
percurso estrutural, letra por letra. Elas nao entram em `Tile`: existem para
que a troca de implementacao seja **conferida pelo kernel** e nao apenas
revisada a olho. As regressoes que as comparam estao na secao de regressoes,
depois do `unseal`. -/

private def dot (size : Nat) (left right : Nat → Rat) : Rat :=
  (List.range size).foldl (fun total index => total + left index * right index) 0

/-- Especificacao indexada de `RatMatrix.mul`. -/
def RatMatrix.mulSpec (left right : RatMatrix) : RatMatrix where
  rows := left.rows
  cols := right.cols
  data := (List.range (left.rows * right.cols)).map fun index =>
    let row := index / right.cols
    let col := index % right.cols
    dot left.cols (fun k => left.entry row k) (fun k => right.entry k col)

/-- Especificacao indexada de `RatMatrix.sub`. -/
def RatMatrix.subSpec (left right : RatMatrix) : RatMatrix where
  rows := left.rows
  cols := left.cols
  data := (List.range (left.rows * left.cols)).map fun index =>
    let row := index / left.cols
    let col := index % left.cols
    left.entry row col - right.entry row col

/-- Especificacao indexada de `RatMatrix.add`. -/
def RatMatrix.addSpec (left right : RatMatrix) : RatMatrix where
  rows := left.rows
  cols := left.cols
  data := (List.range (left.rows * left.cols)).map fun index =>
    let row := index / left.cols
    let col := index % left.cols
    left.entry row col + right.entry row col

private def RatMatrix.rowAbsSum (matrix : RatMatrix) (row : Nat) : Rat :=
  (List.range matrix.cols).foldl
    (fun total col => total + qabs (matrix.entry row col)) 0

/-- Especificacao indexada de `RatMatrix.infNorm`. -/
def RatMatrix.infNormSpec (matrix : RatMatrix) : Rat :=
  (List.range matrix.rows).foldl
    (fun largest row => qmax largest (matrix.rowAbsSum row)) 0

/-- Especificacao indexada de `RatMatrix.transpose`. -/
def RatMatrix.transposeSpec (matrix : RatMatrix) : RatMatrix where
  rows := matrix.cols
  cols := matrix.rows
  data := (List.range (matrix.cols * matrix.rows)).map fun index =>
    matrix.entry (index % matrix.rows) (index / matrix.rows)

/-- Data for one tile of the spectral counting contour.

`centerMatrix` is the finite matrix at the tile centre -- exactly `F_N(z_j)`
when `centerResidual` is empty, and a rounding of it otherwise, the discarded
part being `centerResidual` -- and `approxInverse` is a proposed inverse `V_j`.
`maxDistance` must bound `|z-z₀|`, and `tailError` must bound `||C-C_N||`.

`tileRadius`, `matrixNorm` e `centerResidual` substituem o antigo campo
`matrixVariation`: `raio * ||C_N||` **e** a variacao da matriz no tile,
exatamente, e `||centerResidual||_1` e o que separa o centro escrito do centro
verdadeiro (ver o cabecalho do modulo e `Tile.matrixVariation`).
`operatorBound` majora `norma(C) <= M_J * M_0` e e o preco do truncamento
unilateral adotado em `transferMajorant`.
-/
structure Tile where
  centerMatrix : RatMatrix
  approxInverse : RatMatrix
  /-- Raio do tile em torno de `z_j`. -/
  tileRadius : Rat
  /-- Majorante de `norma(C_N)`, calculado sobre a matriz racional. -/
  matrixNorm : Rat
  maxDistance : Rat
  tailError : Rat
  /-- Majorante de `norma(C) = norma(J R_0) <= M_J * M_0`. -/
  operatorBound : Rat
  /-- Residuo de arredondamento do centro.

  **Obrigacao do produtor, entrada por entrada:**
  `centerMatrix.entry i j + centerResidual.entry i j = F_N(z_j) i j`
  para todo `i, j`, com `entry` devolvendo `0` fora da forma declarada.

  Esta obrigacao e da **mesma especie** que a que o arquivo ja carregava para
  `centerMatrix` sozinho -- "esta matriz e `F_N(z_j)`" --, e o caso antigo e o
  caso particular `centerResidual = 0x0`, que e o valor padrao do campo. O que
  muda e que agora `centerMatrix` pode ser uma versao **arredondada** de
  `F_N(z_j)`, com o residuo exato guardado a parte, e o kernel calcula
  `Tile.centerError = ||centerResidual||_1` em vez de receber um numero
  estimado. -/
  centerResidual : RatMatrix := { rows := 0, cols := 0, data := [] }
deriving Repr, DecidableEq

/-- `||F_N(z_j) - centerMatrix||_1`, recomputado pelo kernel a partir do
residuo. Vale `0` quando o centro e exato (`centerResidual` vazio). -/
def Tile.centerError (tile : Tile) : Rat :=
  tile.centerResidual.oneNorm

/-- Majorante de `||F_N(z) - centerMatrix||` sobre o tile inteiro.

Pela desigualdade triangular,
`||F_N(z) - centerMatrix|| <= ||F_N(z) - F_N(z_j)|| + ||F_N(z_j) - centerMatrix||`.
A primeira parcela e `raio * ||C_N||`, com igualdade no pior caso para o pencil
afim, e a segunda e exatamente `Tile.centerError`. Nao ha hipotese nova: a
primeira e o produto de dois dados que o produtor ja tem e a segunda e a norma
de uma matriz escrita no certificado, recomputada pelo kernel. -/
def Tile.matrixVariation (tile : Tile) : Rat :=
  tile.tileRadius * tile.matrixNorm + tile.centerError

def Tile.shapeCorrect (tile : Tile) : Prop :=
  tile.centerMatrix.wellFormed ∧
  tile.approxInverse.wellFormed ∧
  tile.centerMatrix.rows = tile.centerMatrix.cols ∧
  tile.approxInverse.rows = tile.centerMatrix.rows ∧
  tile.approxInverse.cols = tile.centerMatrix.cols ∧
  tile.centerResidual.wellFormed

/-- `||I - V F₀|| + ||V|| δ_F`. -/
def Tile.inverseDefect (tile : Tile) : Rat :=
  let size := tile.centerMatrix.rows
  let residual := (RatMatrix.identity size).sub
    (tile.approxInverse.mul tile.centerMatrix)
  residual.oneNorm + tile.approxInverse.oneNorm * tile.matrixVariation

/-- Neumann upper bound `||V|| / (1-η)` for the inverse on a tile. -/
def Tile.inverseBound (tile : Tile) : Rat :=
  tile.approxInverse.oneNorm / (1 - tile.inverseDefect)

/-- Majorante de transferencia de cauda, na forma **unilateral**
`|z-z0| * eps_N * (1 + |z-z0| * beta_N * ||C||)`.

Esta e a forma correta para `C_N = Pi C` (truncamento so a esquerda), que tem o
mesmo determinante finito pelo lema de Sylvester; ela aparece via Woodbury. A
forma antiga, `|z-z0| * eps_N * max(1, beta_N)`, e a do truncamento **bilateral**
`Pi C Pi` e consome a hipotese `P2*` -- regularidade analitica do adjunto --,
que estava tacita e nao demonstrada. Trocar elimina essa hipotese; o preco,
medido em `DeterminantBound`, e um fator `30,02` no majorante e `N_F : 9 -> 13`.
Ver `docs/WINDING_DATA.md`. -/
def Tile.transferMajorant (tile : Tile) : Rat :=
  tile.maxDistance * tile.tailError *
    (1 + tile.maxDistance * tile.inverseBound * tile.operatorBound)

/-- Mathematical content checked for a contour tile. -/
def Tile.inequalities (tile : Tile) : Prop :=
  tile.shapeCorrect ∧
  0 ≤ tile.tileRadius ∧
  0 ≤ tile.matrixNorm ∧
  0 ≤ tile.maxDistance ∧
  0 ≤ tile.tailError ∧
  0 ≤ tile.operatorBound ∧
  tile.inverseDefect < 1 ∧
  tile.transferMajorant < 1

private def Tile.inequalitiesDecidable (tile : Tile) :
    Decidable tile.inequalities := by
  unfold Tile.inequalities Tile.shapeCorrect RatMatrix.wellFormed
  infer_instance

/-- Executable, exact-rational tile checker. -/
def Tile.valid (tile : Tile) : Bool :=
  @decide tile.inequalities tile.inequalitiesDecidable

/-- The executable checker cannot accept a tile unless all advertised
Neumann and tail inequalities are true in `Rat`. -/
theorem Tile.valid_sound (tile : Tile) (accepted : tile.valid = true) :
    tile.inequalities := by
  exact @of_decide_eq_true tile.inequalities tile.inequalitiesDecidable accepted

theorem Tile.inverseDefect_lt_one (tile : Tile) (accepted : tile.valid = true) :
    tile.inverseDefect < 1 :=
  (tile.valid_sound accepted).2.2.2.2.2.2.1

theorem Tile.transferMajorant_lt_one (tile : Tile) (accepted : tile.valid = true) :
    tile.transferMajorant < 1 :=
  (tile.valid_sound accepted).2.2.2.2.2.2.2

/-- A variacao da matriz nao precisa ser exigida como hipotese: ela e
**derivada**, e sua nao negatividade e consequencia da dos dois fatores mais a
de uma norma. -/
theorem Tile.matrixVariation_nonneg (tile : Tile) (accepted : tile.valid = true) :
    0 ≤ tile.matrixVariation :=
  Rat.add_nonneg
    (Rat.mul_nonneg (tile.valid_sound accepted).2.1 (tile.valid_sound accepted).2.2.1)
    tile.centerResidual.oneNorm_nonneg

/-- O erro de arredondamento do centro e nao negativo em qualquer tile, sem
hipotese nenhuma: e uma norma. -/
theorem Tile.centerError_nonneg (tile : Tile) : 0 ≤ tile.centerError :=
  tile.centerResidual.oneNorm_nonneg

structure BoundaryCertificate where
  tiles : List Tile
  nonempty : tiles ≠ []

/-- Every contour tile must pass; `List.all` short-circuits on failure. -/
def BoundaryCertificate.valid (certificate : BoundaryCertificate) : Bool :=
  certificate.tiles.all Tile.valid

/-- Sound aggregation over the entire contour. -/
theorem BoundaryCertificate.valid_sound
    (certificate : BoundaryCertificate)
    (accepted : certificate.valid = true) :
    ∀ tile ∈ certificate.tiles, tile.inequalities := by
  intro tile member
  apply Tile.valid_sound tile
  exact List.all_eq_true.mp accepted tile member

/-! ## Regressoes, todas fechadas pelo kernel

`unseal` abaixo remove o atributo `@[irreducible]` de quatro funcoes de `Rat`.
**Isso nao aumenta a base de confianca.** Irredutibilidade e uma diretiva de
*elaboracao*: o kernel nunca a respeitou. O que `unseal` faz e devolver ao
elaborador a permissao de desdobrar `Rat.mul`, `Rat.add`, `Rat.sub` e `Rat.inv`
ao montar o termo de prova; o termo resultante continua sendo checado inteiro
pelo kernel, exatamente como qualquer outro. A diferenca para a avaliacao
nativa e categorica: aquela delega a computacao ao compilador e sela o
resultado num axioma opaco gerado por declaracao, que aparece no
`#print axioms`. Os `#print axioms` no fim deste arquivo listam apenas
`propext`, `Classical.choice` e `Quot.sound`.

`Rat.div` nao precisa de `unseal`: ja e `semireducible` e reduz sozinha assim
que `Rat.inv` esta liberada -- tentar libera-la da erro de reducibilidade. -/

unseal Rat.mul Rat.add Rat.sub Rat.inv

set_option maxRecDepth 100000

/-! Um certificado exato pequeno, usado como regressao.

`F0 = [2]`, `V = [1/2]`, `raio = 1/10` e `||C_N|| = 1`, logo
`matrixVariation = 1/10`, `eta = 1/20` e o majorante do inverso e `10/19`.
Numa matriz `1x1` a soma por linha e a soma por coluna coincidem, de modo que a
troca de `l^infty` por `l^1` **nao muda** `inverseDefect` nem `inverseBound`
aqui -- os dois continuam `1/20` e `10/19`, como antes.

O que muda e o majorante de transferencia. Com `|z-z0| = 3`, `eps_N = 1/100` e
`||C|| = 1`, a forma unilateral da
`3 * (1/100) * (1 + 3 * (10/19) * 1) = 147/1900 ~ 0,077`,
contra `3/100 = 0,03` da forma bilateral antiga. O fator entre as duas e o
preco de nao supor `P2*`.
-/

def scalarMatrix (x : Rat) : RatMatrix where
  rows := 1
  cols := 1
  data := [x]

def regressionTile : Tile where
  centerMatrix := scalarMatrix 2
  approxInverse := scalarMatrix (1 / 2)
  tileRadius := 1 / 10
  matrixNorm := 1
  maxDistance := 3
  tailError := 1 / 100
  operatorBound := 1

example : regressionTile.matrixVariation = 1 / 10 := by decide
example : regressionTile.inverseDefect = 1 / 20 := by decide
example : regressionTile.inverseBound = 10 / 19 := by decide
example : regressionTile.transferMajorant = 147 / 1900 := by decide
example : regressionTile.valid = true := by decide

def regressionBoundary : BoundaryCertificate where
  tiles := [regressionTile, regressionTile]
  nonempty := by decide

example : regressionBoundary.valid = true := by decide

/-! ### Tiles maiores

As matrizes usam a estrutura de Neumann `F0 = I + E`, `V = I - E`, de modo que
`V F0 = I - E^2` e o residuo e genuinamente pequeno. As entradas de `E` variam
com o indice (`(indice mod 7) + 1`) para que o produto nao tenha entradas
repetidas: com entradas todas iguais o kernel reaproveitaria subtermos e o
tempo medido seria otimista. Os valores exatos foram conferidos por fora, em
`fractions.Fraction`, antes de serem escritos aqui. -/

/-- `E` com entradas `((indice mod 7) + 1) / denominator`. -/
def rampMatrix (size denominator : Nat) : RatMatrix where
  rows := size
  cols := size
  data := (List.range (size * size)).map fun index =>
    mkRat (Int.ofNat (index % 7 + 1)) denominator

def neumannCenter (size denominator : Nat) : RatMatrix :=
  (RatMatrix.identity size).add (rampMatrix size denominator)

def neumannInverse (size denominator : Nat) : RatMatrix :=
  (RatMatrix.identity size).sub (rampMatrix size denominator)

def neumannTile (size denominator : Nat) : Tile where
  centerMatrix := neumannCenter size denominator
  approxInverse := neumannInverse size denominator
  tileRadius := 1 / 512
  matrixNorm := 8
  maxDistance := 2
  tailError := 1 / 256
  operatorBound := 10

/-- Tile `4x4`. -/
def tile4 : Tile := neumannTile 4 128

-- Sob `l^1` (coluna) o defeito cai: `279/8192 ~ 0,03406` virou
-- `545/16384 ~ 0,03326`. A norma por coluna de `V` e a por linha coincidem
-- aqui (`35/32`), mas as do residuo nao: `265/16384` contra `139/8192`.
example : tile4.matrixVariation = 1 / 64 := by decide
example : tile4.inverseDefect = 545 / 16384 := by decide
example : tile4.inverseBound = 17920 / 15839 := by decide
example : tile4.transferMajorant = 374239 / 2027392 := by decide
example : tile4.valid = true := by decide

/-- Tile `8x8`. -/
def tile8 : Tile := neumannTile 8 128

-- Neste `8x8` as duas normas coincidem por acidente da matriz-rampa, de modo
-- que `inverseDefect` e `inverseBound` nao mudaram com a troca para `l^1`. So
-- o majorante de transferencia mudou, e mudou por causa da forma unilateral.
example : tile8.inverseDefect = 355 / 4096 := by decide
example : tile8.inverseBound = 5024 / 3741 := by decide
example : tile8.transferMajorant = 104221 / 478848 := by decide
example : tile8.valid = true := by decide

/-! ### A norma por coluna e mesmo diferente da norma por linha -/

/-- `[[1, 2], [0, 0]]`: somas por linha `3` e `0`; somas por coluna `1` e `2`.
A norma `l^infty` e `3`, a norma `l^1` e `2`. -/
def asymmetricMatrix : RatMatrix where
  rows := 2
  cols := 2
  data := [1, 2, 0, 0]

example : asymmetricMatrix.infNorm = 3 := by decide
example : asymmetricMatrix.oneNorm = 2 := by decide
example : asymmetricMatrix.transpose.data = [1, 0, 2, 0] := by decide

/-! ### A implementacao estrutural concorda com a especificacao indexada

A troca de `entry` (indexacao de lista, `O(n^5)` no produto) por percurso
estrutural (`O(n^3)`) nao e uma refatoracao que se aceita por leitura: as
regressoes abaixo mandam o **kernel** recomputar as duas definicoes e exigem
igualdade dos objetos inteiros -- forma e dados. Ha um caso **retangular** de
proposito, porque e onde a aritmetica de indices `row * cols + col` erra sem
que um caso quadrado denuncie. -/

/-- `2x3`, entradas distintas. -/
def rect23 : RatMatrix where
  rows := 2
  cols := 3
  data := [1, 2, 3, 4, 5, 6]

/-- `3x2`, entradas distintas. -/
def rect32 : RatMatrix where
  rows := 3
  cols := 2
  data := [7, 8, 9, 10, 11, 12]

example : rect23.mul rect32 = rect23.mulSpec rect32 := by decide
example : rect32.mul rect23 = rect32.mulSpec rect23 := by decide
example : rect23.transpose = rect23.transposeSpec := by decide
example : rect32.transpose = rect32.transposeSpec := by decide
example : rect23.infNorm = rect23.infNormSpec := by decide
example : rect23.transpose.infNorm = rect23.transposeSpec.infNormSpec := by decide
example : rect23.sub rect23 = rect23.subSpec rect23 := by decide
example : rect23.add rect23 = rect23.addSpec rect23 := by decide

example : asymmetricMatrix.transpose = asymmetricMatrix.transposeSpec := by decide
example : asymmetricMatrix.infNorm = asymmetricMatrix.infNormSpec := by decide

-- E nas matrizes que o certificado usa de fato, `4x4` e `6x6`.
example :
    (neumannInverse 4 128).mul (neumannCenter 4 128)
      = (neumannInverse 4 128).mulSpec (neumannCenter 4 128) := by decide
example :
    (neumannInverse 6 128).mul (neumannCenter 6 128)
      = (neumannInverse 6 128).mulSpec (neumannCenter 6 128) := by decide
example : (neumannCenter 6 128).transpose = (neumannCenter 6 128).transposeSpec := by decide
example : (neumannCenter 6 128).oneNorm = (neumannCenter 6 128).transposeSpec.infNormSpec := by
  decide

/-! ### Centro arredondado, com o arredondamento absorvido

O centro exato e `F_N(z_j) = [2001/1000]`. O tile guarda o centro
**arredondado** `[2]` e o residuo exato `[1/1000]`, cuja soma e o centro
exato -- que e a obrigacao declarada no campo `centerResidual`. O kernel
calcula `centerError = 1/1000` e o soma a `raio * ||C_N|| = 1/10`, de modo que
`matrixVariation = 101/1000` e `inverseDefect = 0 + (1/2)(101/1000) = 101/2000`.
Compare com `regressionTile`, de centro exato: la o defeito e `1/20 = 50/1000`.
O arredondamento e pago, e aparece no numero. -/
def roundedCenterTile : Tile where
  centerMatrix := scalarMatrix 2
  approxInverse := scalarMatrix (1 / 2)
  tileRadius := 1 / 10
  matrixNorm := 1
  maxDistance := 3
  tailError := 1 / 100
  operatorBound := 1
  centerResidual := scalarMatrix (1 / 1000)

example : roundedCenterTile.centerError = 1 / 1000 := by decide
example : roundedCenterTile.matrixVariation = 101 / 1000 := by decide
example : roundedCenterTile.inverseDefect = 101 / 2000 := by decide
example : roundedCenterTile.valid = true := by decide

/-- **Rejeicao: arredondamento grande demais.** Mesmo tile, residuo `4`: o
centro escrito esta tao longe do centro verdadeiro que a serie de Neumann nao
converge mais. `matrixVariation = 1/10 + 4` e `inverseDefect = 41/20 >= 1`.
E o mecanismo que impede trocar bits por um certificado vazio. -/
def rejectedRoundingTile : Tile :=
  { roundedCenterTile with centerResidual := scalarMatrix 4 }

example : rejectedRoundingTile.inverseDefect = 41 / 20 := by decide
example : rejectedRoundingTile.valid = false := by decide

/-- **Rejeicao: residuo mal formado.** `data` tem uma entrada e a forma
anunciada pede quatro. `shapeCorrect` recusa antes de qualquer norma. -/
def rejectedResidualShapeTile : Tile :=
  { roundedCenterTile with
      centerResidual := { rows := 2, cols := 2, data := [1] } }

example : rejectedResidualShapeTile.valid = false := by decide

/-! ### Rejeicoes -/

/-- **Rejeicao: defeito de inverso `>= 1`.** `V = [1/2]` continua sendo o
inverso exato de `F0 = [2]`, mas a variacao anunciada da matriz no tile e tao
grande que a serie de Neumann nao converge: `eta = 0 + (1/2)*4 = 2 >= 1`. -/
def rejectedDefectTile : Tile where
  centerMatrix := scalarMatrix 2
  approxInverse := scalarMatrix (1 / 2)
  tileRadius := 4
  matrixNorm := 1
  maxDistance := 3
  tailError := 1 / 100
  operatorBound := 1

example : rejectedDefectTile.inverseDefect = 2 := by decide
example : rejectedDefectTile.valid = false := by decide

/-- **Rejeicao: majorante de transferencia `>= 1`.** Aqui `eta = 1/20 < 1` e o
tile passaria na primeira condicao, mas a cauda anunciada e grande demais:
`transferMajorant = 3 * 1 * (1 + 3 * (10/19) * 1) = 147/19 >= 1`. -/
def rejectedTransferTile : Tile where
  centerMatrix := scalarMatrix 2
  approxInverse := scalarMatrix (1 / 2)
  tileRadius := 1 / 10
  matrixNorm := 1
  maxDistance := 3
  tailError := 1
  operatorBound := 1

example : rejectedTransferTile.inverseDefect = 1 / 20 := by decide
example : rejectedTransferTile.transferMajorant = 147 / 19 := by decide
example : rejectedTransferTile.valid = false := by decide

/-- **Rejeicao: o preco exato de abandonar `P2*`.** Este tile passaria no
majorante **bilateral** antigo e e recusado pelo **unilateral** novo, com os
mesmos dados. E a regressao que torna o custo da decisao visivel em vez de
apenas documentado: com `||C|| = 10` e `|z-z0| = 3`, o fator
`1 + 3 * (10/19) * 10 = 319/19 ~ 16,8` substitui `max(1, 10/19) = 1`. -/
def rejectedUnilateralTile : Tile where
  centerMatrix := scalarMatrix 2
  approxInverse := scalarMatrix (1 / 2)
  tileRadius := 1 / 10
  matrixNorm := 1
  maxDistance := 3
  tailError := 1 / 10
  operatorBound := 10

-- O majorante bilateral antigo, escrito a mao com os mesmos dados: passaria.
example :
    rejectedUnilateralTile.maxDistance * rejectedUnilateralTile.tailError *
      qmax 1 rejectedUnilateralTile.inverseBound = 3 / 10 := by decide
-- O unilateral, que e o que o certificado agora exige: nao passa.
example : rejectedUnilateralTile.transferMajorant = 957 / 190 := by decide
example : rejectedUnilateralTile.inverseDefect = 1 / 20 := by decide
example : rejectedUnilateralTile.valid = false := by decide

/-- **Rejeicao: forma errada.** `V` e `1x2` e nao casa com `F0`. -/
def rejectedShapeTile : Tile where
  centerMatrix := scalarMatrix 2
  approxInverse := { rows := 1, cols := 2, data := [1, 1] }
  tileRadius := 1 / 10
  matrixNorm := 1
  maxDistance := 3
  tailError := 1 / 100
  operatorBound := 1

example : rejectedShapeTile.valid = false := by decide

/-- **Rejeicao: um tile ruim entre tiles bons derruba o contorno inteiro.** -/
def rejectedBoundary : BoundaryCertificate where
  tiles := [regressionTile, rejectedTransferTile, regressionTile]
  nonempty := by decide

example : rejectedBoundary.valid = false := by decide

/-! Winding discreto: o programa externo fornece incrementos de cruzamento
orientado já verificados em cada aresta do contorno (por exemplo, via um
certificado intervalar de sinais). Cada aresta pode contribuir `-1`, `0` ou
`1`; Lean confere essa condição, rejeita o contorno vazio, soma exatamente os
incrementos e compara a soma ao alvo. -/

structure WindingCertificate where
  increments : List Int
  target : Int
deriving Repr, DecidableEq

def WindingCertificate.total (certificate : WindingCertificate) : Int :=
  certificate.increments.foldl (fun total step => total + step) 0

def WindingCertificate.incrementAdmissible (increment : Int) : Prop :=
  increment = -1 ∨ increment = 0 ∨ increment = 1

def WindingCertificate.claims (certificate : WindingCertificate) : Prop :=
  certificate.increments ≠ [] ∧
  (∀ increment ∈ certificate.increments,
    WindingCertificate.incrementAdmissible increment) ∧
  certificate.total = certificate.target

private def WindingCertificate.claimsDecidable
    (certificate : WindingCertificate) : Decidable certificate.claims := by
  unfold WindingCertificate.claims WindingCertificate.incrementAdmissible
  infer_instance

def WindingCertificate.valid (certificate : WindingCertificate) : Bool :=
  @decide certificate.claims certificate.claimsDecidable

theorem WindingCertificate.valid_sound
    (certificate : WindingCertificate)
    (accepted : certificate.valid = true) :
    certificate.claims := by
  exact @of_decide_eq_true certificate.claims certificate.claimsDecidable accepted

theorem WindingCertificate.total_eq_target_of_valid
    (certificate : WindingCertificate)
    (accepted : certificate.valid = true) :
    certificate.total = certificate.target :=
  (certificate.valid_sound accepted).2.2

def regressionWinding : WindingCertificate where
  increments := [1, -1, 1]
  target := 1

-- Nivel K: aritmetica inteira pura, reduz no kernel (ver docs/KERNEL_TRUST.md).
example : regressionWinding.valid = true := by decide
example : regressionWinding.total = 1 := by decide

def rejectedWinding : WindingCertificate where
  increments := [2]
  target := 2

example : rejectedWinding.valid = false := by decide

/-! ### `EdgeSignCertificate` (SUPERADO)

Certificado local de sinal: um intervalo racional fechado que nao cruza zero
fixa o sinal de uma quantidade real.

**Este certificado foi superado por
`Choptuik.WindingInterval.IntervalWindingCertificate` e nao deve ser usado em
contagem de winding.** Dois motivos, e o segundo e o grave:

1. Ele esta desligado de tudo: nada em `WindingCertificate` exige que um
   incremento venha de um `EdgeSignCertificate`. Ele certifica um fato que
   ninguem consome.
2. Ele e intrinsecamente fraco para a tarefa. O que ele certifica e o sinal de
   um intervalo **real**; o que a contagem por principio do argumento precisa e
   o cruzamento de um **raio no plano complexo** por uma aresta do contorno.
   Sao coisas diferentes: `sign = 1` com `0 < lower` nao tem relacao nenhuma
   com o incremento de winding de uma aresta. `IntervalWindingCertificate` faz
   a coisa certa, com caixas complexas e classificacao do cruzamento de
   `R_{<0}`, e liga o resultado a `WindingCertificate` pela ponte
   `IntervalWindingCertificate.toWinding`.

Mantido no arquivo porque remover API publica e decisao do dono do
repositorio. -/
structure EdgeSignCertificate where
  lower : Rat
  upper : Rat
  sign : Int
deriving Repr, DecidableEq

def EdgeSignCertificate.claims (edge : EdgeSignCertificate) : Prop :=
  edge.lower ≤ edge.upper ∧
  ((edge.sign = 1 ∧ 0 < edge.lower) ∨
   (edge.sign = -1 ∧ edge.upper < 0))

private def EdgeSignCertificate.claimsDecidable
    (edge : EdgeSignCertificate) : Decidable edge.claims := by
  unfold EdgeSignCertificate.claims
  infer_instance

def EdgeSignCertificate.valid (edge : EdgeSignCertificate) : Bool :=
  @decide edge.claims edge.claimsDecidable

theorem EdgeSignCertificate.valid_sound
    (edge : EdgeSignCertificate) (accepted : edge.valid = true) :
    edge.claims := by
  exact @of_decide_eq_true edge.claims edge.claimsDecidable accepted

-- Nivel K: `EdgeSignCertificate` so COMPARA racionais, nunca os soma nem os
-- divide. Comparacoes de `mkRat n d` reduzem no kernel; literais escritos com
-- `/` nao, porque passam por `Rat.div`. Escrevendo com `mkRat`, estas
-- regressoes dispensam avaliacao nativa. Ver docs/KERNEL_TRUST.md.
def regressionPositiveEdge : EdgeSignCertificate where
  lower := mkRat 1 10
  upper := mkRat 1 5
  sign := 1

example : regressionPositiveEdge.valid = true := by decide

def regressionCrossingEdge : EdgeSignCertificate where
  lower := mkRat (-1) 10
  upper := mkRat 1 10
  sign := 1

example : regressionCrossingEdge.valid = false := by decide

#print axioms RatMatrix.oneNorm_eq_transpose_infNorm
#print axioms Tile.valid_sound
#print axioms Tile.inverseDefect_lt_one
#print axioms Tile.matrixVariation_nonneg
#print axioms Tile.centerError_nonneg
#print axioms RatMatrix.oneNorm_nonneg
#print axioms Tile.transferMajorant_lt_one
#print axioms BoundaryCertificate.valid_sound
#print axioms WindingCertificate.valid_sound
#print axioms WindingCertificate.total_eq_target_of_valid
#print axioms EdgeSignCertificate.valid_sound

end Choptuik.RationalCertificate
