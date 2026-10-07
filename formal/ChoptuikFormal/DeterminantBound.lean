import Std
import ChoptuikFormal.TailBound

/-!
# O elo que falta entre o operador e o certificado de winding

`ChoptuikFormal.WindingInterval` consome, por aresta, três caixas racionais e
exige que `edgeBox` contenha `f(z) = det[I + (z-z₀)C_N]` para **todo** `z` do
segmento. Hoje nada no repositório produz esse enclosure. Este módulo compõe,
em `Rat` exato, as duas estimativas que o produzem a partir de dados que a
maquinaria de tiles já calcula:

1. **Derivada logarítmica (fórmula de Jacobi).** `f'(z) = f(z)·tr(F(z)⁻¹C_N)`,
   logo `|f'/f| ≤ K` com `K` composto em `LogDerivativeData.bound`. O ponto
   importante é que `tr` **não** é majorado por `n·‖·‖`: o termo dominante
   `tr(V C_N)` é calculado, e só o resto é majorado (por Cauchy--Schwarz de
   Frobenius). Sem isso o fator `n = 612` tornaria a malha de arestas
   inviável.
2. **Raio de inflação da aresta.** De `|f'| ≤ K|f|` sai
   `sup_aresta |f| ≤ m/(1-K·L)` e `|f(z) - f(z_k)| ≤ K·L·m/(1-K·L) =: δ`,
   com `m ≥ |f(z_k)|` e `L` o comprimento da aresta. `edgeBox` é então
   `startBox` inflado por `δ`.

Além disso o módulo fixa duas decisões de desenho analisadas em
`docs/TAIL_BOUND.md` e `docs/WINDING_DATA.md`:

3. **Majorante de transferência unilateral.** Com `C_N = Π C` (que tem o mesmo
   determinante finito, pelo lema de Sylvester), o majorante correto é
   `|λ|·ε_N·(1 + |λ|·β_N·‖C‖)`, e **não** `|λ|·ε_N·max(1, β_N)`, que é o do
   truncamento bilateral `Π C Π` e consome a hipótese `P2*` (regularidade
   analítica do adjunto). Aqui os dois são computados lado a lado para que o
   preço da troca fique explícito e verificado.
4. **`matrixVariation` não é um bound de EDP.** Para o pencil afim
   `F_N(z) = I + (z-z₀)C_N` vale `F_N(z) - F_N(z_j) = (z-z_j)C_N` exatamente,
   logo `‖F_N(z) - F_N(z_j)‖ ≤ raio_do_tile · ‖C_N‖`. É aritmética sobre a
   matriz racional, não uma estimativa a demonstrar.

Como em `TailBound`, nada aqui é um teorema sobre o operador de Choptuik: o que
o kernel aceita é a composição. Que `traceTerm`, `residualFrobenius`,
`inverseBound` e `matrixFrobenius` sejam majorantes verdadeiros continua sendo
obrigação do produtor, e está especificado em `docs/WINDING_DATA.md`.
-/

namespace Choptuik.DeterminantBound

/- Mesma observação de `TailBound`: `Rat.mul`, `Rat.add`, `Rat.sub` e `Rat.inv`
   são `@[irreducible]` no Lean 4.33 e `unseal` devolve a transparência normal,
   sem acrescentar hipótese nenhuma ao kernel. Todas as regressões fecham por
   `decide`. -/
unseal Rat.mul Rat.add Rat.sub Rat.inv

/-! ## 1. Derivada logarítmica pela fórmula de Jacobi -/

/-- Dados para majorar `|f'/f| = |tr(F(z)⁻¹ C_N)|` num tile.

Decomposição usada, com `V` o inverso aproximado do tile:

```
tr(F⁻¹C_N) = tr(V C_N) + tr((F⁻¹ - V) C_N) ,
|tr((F⁻¹-V)C_N)| ≤ ‖F⁻¹-V‖_F · ‖C_N‖_F ≤ ‖I - V F‖_F · ‖F⁻¹‖ · ‖C_N‖_F .
```

`traceTerm` é `|tr(V C_N)|`, calculado exatamente pelo produtor (uma soma
`Σ_{ij} V_{ij} (C_N)_{ji}`, barata). `residualFrobenius` majora `‖I - V F(z)‖_F`
**em todo o tile**, `inverseBound` majora `‖F(z)⁻¹‖` e `matrixFrobenius` majora
`‖C_N‖_F`. Frobenius envolve raiz quadrada: o produtor entrega um majorante
racional, e é ele que entra aqui. -/
structure LogDerivativeData where
  /-- `|tr(V C_N)|`, calculado, não majorado. -/
  traceTerm : Rat
  /-- Majorante racional de `‖I - V F(z)‖_F` no tile. -/
  residualFrobenius : Rat
  /-- Majorante de `‖F(z)⁻¹‖` no tile (o `inverseBound` do tile de Neumann). -/
  inverseBound : Rat
  /-- Majorante racional de `‖C_N‖_F`. -/
  matrixFrobenius : Rat
deriving Repr, DecidableEq

/-- `K ≥ sup |f'/f|` sobre o tile. -/
def LogDerivativeData.bound (data : LogDerivativeData) : Rat :=
  data.traceTerm + data.residualFrobenius * data.inverseBound * data.matrixFrobenius

def LogDerivativeData.admissible (data : LogDerivativeData) : Prop :=
  0 ≤ data.traceTerm ∧ 0 ≤ data.residualFrobenius ∧
  0 ≤ data.inverseBound ∧ 0 ≤ data.matrixFrobenius

/-! ## 2. Raio de inflação da aresta -/

/-- Dados de uma aresta do contorno.

`nodeModulus` é um majorante racional de `|f(z_k)|`; a partir de uma caixa
racional basta tomar `|Re|_max + |Im|_max`, que majora o módulo sem raiz
quadrada. `edgeLength` é o comprimento do segmento — racional, porque o
contorno de `SECTOR_B.md` §8 é um retângulo de lados paralelos aos eixos. -/
structure EdgeEnclosureData where
  /-- `K` da seção 1. -/
  logDerivative : Rat
  /-- `L`, comprimento da aresta. -/
  edgeLength : Rat
  /-- `m ≥ |f(z_k)|`. -/
  nodeModulus : Rat
deriving Repr, DecidableEq

/-- `K·L`, a quantidade que controla tudo. -/
def EdgeEnclosureData.product (data : EdgeEnclosureData) : Rat :=
  data.logDerivative * data.edgeLength

/-- `sup_aresta |f| ≤ m/(1-K·L)`, de `S ≤ m + K·L·S`. -/
def EdgeEnclosureData.supModulus (data : EdgeEnclosureData) : Rat :=
  data.nodeModulus / (1 - data.product)

/-- `δ = K·L·m/(1-K·L) ≥ |f(z) - f(z_k)|` na aresta. É por este `δ` que o
produtor infla `startBox` para obter `edgeBox`. -/
def EdgeEnclosureData.radius (data : EdgeEnclosureData) : Rat :=
  data.product * data.supModulus

/-- Condições sem as quais o raio não significa nada. `product < 1/2` é o que
garante `δ < m`, isto é, que a caixa inflada ainda pode excluir a origem; sem
isso nenhuma aresta classifica e o contorno é inútil. -/
def EdgeEnclosureData.admissible (data : EdgeEnclosureData) : Prop :=
  0 ≤ data.logDerivative ∧ 0 < data.edgeLength ∧ 0 < data.nodeModulus ∧
  data.product < 1 / 2 ∧ data.radius < data.nodeModulus

private def EdgeEnclosureData.admissibleDecidable (data : EdgeEnclosureData) :
    Decidable data.admissible := by
  unfold EdgeEnclosureData.admissible
  infer_instance

def EdgeEnclosureData.valid (data : EdgeEnclosureData) : Bool :=
  @decide data.admissible data.admissibleDecidable

theorem EdgeEnclosureData.valid_sound (data : EdgeEnclosureData)
    (accepted : data.valid = true) : data.admissible :=
  @of_decide_eq_true data.admissible data.admissibleDecidable accepted

/-- A inflação não engole o valor no nó: é o que permite que a aresta ainda
seja classificável pelo cruzamento do raio. -/
theorem EdgeEnclosureData.radius_lt_nodeModulus (data : EdgeEnclosureData)
    (accepted : data.valid = true) : data.radius < data.nodeModulus :=
  (data.valid_sound accepted).2.2.2.2

/-- Composição das seções 1 e 2: os dados de aresta produzidos por um tile. -/
def edgeFromTile (logData : LogDerivativeData) (edgeLength nodeModulus : Rat) :
    EdgeEnclosureData where
  logDerivative := logData.bound
  edgeLength := edgeLength
  nodeModulus := nodeModulus

/-! ## 3. Majorante de transferência: unilateral contra bilateral -/

/-- Dados escalares do majorante de transferência num tile.

`operatorBound` majora `‖C‖ = ‖J R₀‖ ≤ M_J · M₀`, ou seja, é o produto das
constantes de **P3** e **P1** que `TailBound.TailData` já carrega. É o preço do
truncamento unilateral: ele aparece via Woodbury,
`(I + λ Π C)⁻¹ = I - λ ι (I + λ Π C Π|_{V_N})⁻¹ Π C`. -/
structure TransferData where
  /-- `sup_Γ |z - z₀|`. -/
  maxDistance : Rat
  /-- `ε_N`. -/
  tailError : Rat
  /-- `β_N ≥ ‖(I + λ M_N)⁻¹‖`, a quantidade finita do tile. -/
  inverseBound : Rat
  /-- `‖C‖ ≤ M_J · M₀`. -/
  operatorBound : Rat
deriving Repr, DecidableEq

/-- **Majorante unilateral** (`C_N = Π C`), que só consome `P2`:
`|λ|·ε_N·(1 + |λ|·β_N·‖C‖)`. -/
def TransferData.unilateralMajorant (data : TransferData) : Rat :=
  data.maxDistance * data.tailError *
    (1 + data.maxDistance * data.inverseBound * data.operatorBound)

/-- **Majorante bilateral** (`C_N = Π C Π`), que é o implementado hoje em
`Tile.transferMajorant` e consome `P2*`: `|λ|·ε_N·max(1, β_N)`. -/
def TransferData.bilateralMajorant (data : TransferData) : Rat :=
  data.maxDistance * data.tailError *
    (if data.inverseBound < 1 then 1 else data.inverseBound)

def TransferData.admissible (data : TransferData) : Prop :=
  0 ≤ data.maxDistance ∧ 0 ≤ data.tailError ∧
  0 ≤ data.inverseBound ∧ 0 ≤ data.operatorBound

/-- Certificado de um tile na versão unilateral. -/
structure UnilateralTile where
  data : TransferData
deriving Repr, DecidableEq

def UnilateralTile.claims (tile : UnilateralTile) : Prop :=
  tile.data.admissible ∧ tile.data.unilateralMajorant < 1

private def UnilateralTile.claimsDecidable (tile : UnilateralTile) :
    Decidable tile.claims := by
  unfold UnilateralTile.claims TransferData.admissible
  infer_instance

def UnilateralTile.valid (tile : UnilateralTile) : Bool :=
  @decide tile.claims tile.claimsDecidable

theorem UnilateralTile.valid_sound (tile : UnilateralTile)
    (accepted : tile.valid = true) : tile.claims :=
  @of_decide_eq_true tile.claims tile.claimsDecidable accepted

theorem UnilateralTile.majorant_lt_one (tile : UnilateralTile)
    (accepted : tile.valid = true) : tile.data.unilateralMajorant < 1 :=
  (tile.valid_sound accepted).2

/-! ## 4. `matrixVariation` do pencil afim é aritmética, não hipótese -/

/-- Para `F_N(z) = I + (z-z₀)C_N` vale `F_N(z) - F_N(z_j) = (z - z_j)C_N`, logo
`‖F_N(z) - F_N(z_j)‖ ≤ raio · ‖C_N‖` **com igualdade no pior caso**. O campo
`matrixVariation` do tile é portanto calculável a partir da matriz racional; o
ledger o listava junto com os bounds de EDP, e não é um deles. -/
def affineMatrixVariation (tileRadius matrixNorm : Rat) : Rat :=
  tileRadius * matrixNorm

/-! ## Regressões (todas por `decide`, nenhuma por `native_decide`) -/

/-- Dados plausíveis, **não certificados**, de um tile na malha 8×24
(`n = 612`): o termo de traço calculado vale `12`, o resíduo de Frobenius do
inverso aproximado vale `10⁻⁵`, `β = 50` e `‖C_N‖_F ≤ 200`. -/
def sampleLogDerivative : LogDerivativeData where
  traceTerm := 12
  residualFrobenius := 1 / 100000
  inverseBound := 50
  matrixFrobenius := 200

example : sampleLogDerivative.bound = 121 / 10 := by decide

/-- Com `K = 12,1`, aresta de comprimento `1/500` e nó normalizado (`m = 1`,
que é o que a cadeia relativa de `docs/WINDING_DATA.md` entrega), a inflação da
aresta é `δ = 121/4879 ≈ 2,5 %` do módulo do nó. -/
def sampleEdge : EdgeEnclosureData :=
  edgeFromTile sampleLogDerivative (1 / 500) 1

example : sampleEdge.product = 121 / 5000 := by decide
example : sampleEdge.radius = 121 / 4879 := by decide
example : sampleEdge.valid = true := by decide

/-- Aresta longa demais: `K·L = 121/50 > 1`, a estimativa de Gronwall não
fecha e o dado é recusado. É o sinal de que o contorno precisa de mais nós. -/
def sampleEdgeTooLong : EdgeEnclosureData :=
  edgeFromTile sampleLogDerivative (1 / 5) 1

example : sampleEdgeTooLong.valid = false := by decide

/-- Aresta no limite: `K·L` pouco abaixo de `1/2` ainda passa, mas a inflação
já é quase do tamanho do valor no nó (`δ ≈ 0,96·m`), o que na prática impede a
classificação. O verificador aceita; o produtor não deveria chegar perto. -/
def sampleEdgeMarginal : EdgeEnclosureData :=
  edgeFromTile sampleLogDerivative (49 / 1210) 1

example : sampleEdgeMarginal.product = 49 / 100 := by decide
example : sampleEdgeMarginal.radius = 49 / 51 := by decide
example : sampleEdgeMarginal.valid = true := by decide

/-- Dados de transferência: `|λ| ≤ 3`, `β = 50`, `‖C‖ ≤ M_J·M₀ = 10`. São as
mesmas constantes com que `TailBound.neumannTarget = 1/150` foi calibrado. -/
def sampleTransfer (tailError : Rat) : TransferData where
  maxDistance := 3
  tailError := tailError
  inverseBound := 50
  operatorBound := 10

-- O preço exato da troca: o majorante unilateral é `1501/50 = 30,02` vezes o
-- bilateral, com os mesmos dados.
example : (sampleTransfer (1 / 10000)).unilateralMajorant = 4503 / 10000 := by decide
example : (sampleTransfer (1 / 10000)).bilateralMajorant = 3 / 200 := by decide

-- `ε_N = 10⁻³` basta para o bilateral e **não** basta para o unilateral.
example : (sampleTransfer (1 / 1000)).bilateralMajorant < 1 := by decide
example : ¬ ((sampleTransfer (1 / 1000)).unilateralMajorant < 1) := by decide
example : { data := sampleTransfer (1 / 1000) : UnilateralTile }.valid = false := by decide
example : { data := sampleTransfer (1 / 10000) : UnilateralTile }.valid = true := by decide

/-! ### Quanto custa, em modos de Fourier, abandonar `P2*`

Com as constantes `H-opt` de `TailBound` (`M₀ = 10`, `A₀ = 4`, `M_J = 1`,
`q_F = q_C = 2/5`, realização de Wiener) e `N_C = 24`, o bilateral fecha em
`N_F = 9` e o unilateral em `N_F = 13`. **Quatro modos de Fourier é o preço
inteiro de não precisar demonstrar `P2*`.** -/

def unilateralTarget : Rat := 1 / 4503

example : unilateralTarget = 1 / (3 * (1 + 3 * 50 * 10)) := by decide

example : (TailBound.optimisticGrid 9 24).epsilon ≤ TailBound.neumannTarget := by decide
example : ¬ ((TailBound.optimisticGrid 12 24).epsilon ≤ unilateralTarget) := by decide
example : (TailBound.optimisticGrid 13 24).epsilon ≤ unilateralTarget := by decide

example :
    { data := { sampleTransfer (TailBound.optimisticGrid 13 24).epsilon with } : UnilateralTile
    }.valid = true := by decide

example :
    { data := { sampleTransfer (TailBound.optimisticGrid 12 24).epsilon with } : UnilateralTile
    }.valid = false := by decide

/-- `matrixVariation` de um tile de raio `1/64` com `‖C_N‖ = 20`. -/
example : affineMatrixVariation (1 / 64) 20 = 5 / 16 := by decide

#print axioms EdgeEnclosureData.valid_sound
#print axioms EdgeEnclosureData.radius_lt_nodeModulus
#print axioms UnilateralTile.valid_sound
#print axioms UnilateralTile.majorant_lt_one

end Choptuik.DeterminantBound
