import Std
import ChoptuikFormal.RationalCertificate
import ChoptuikFormal.TailBound

/-!
# A costura entre o bound de cauda e o tile de Neumann

Dois módulos deste repositório deveriam se encaixar e, até aqui, não se
falavam:

* `Choptuik.TailBound.TailData.epsilon` **calcula** `ε_N` a partir das
  constantes nomeadas das três hipóteses de EDP (`P1`, `P2`, `P3`);
* `Choptuik.RationalCertificate.Tile.tailError` **consome** um `ε_N`.

Nada no tipo ligava um ao outro. Um produtor podia declarar em `tailError` um
número **menor** que o `ε_N` que as suas próprias constantes produzem, e os dois
módulos aceitavam sem reclamar — o `Tile` porque `tailError` menor só torna
`transferMajorant < 1` mais fácil, e o `TailBound` porque ele nunca olha o tile.
É o mesmo defeito que `WindingCertificate` tinha antes de `WindingInterval`: um
número que vem de lugar nenhum. Este módulo fecha a costura.

## O sentido da desigualdade, e por que é este

A cadeia é:

```
‖C - C_N‖  ≤  ε_N = data.epsilon        (conteúdo de P1/P2/P3, não provado aqui)
           ≤  tile.tailError            (VERIFICADO por este módulo)
```

`tile.tailError` precisa **majorar** `‖C - C_N‖`, e é isso que a composição
acima entrega. O sentido não é convenção: `Tile.transferMajorant` é

```
maxDistance * tailError * (1 + maxDistance * inverseBound * operatorBound)
```

que é **monótona crescente** em `tailError` quando os demais fatores são não
negativos — e `Tile.inequalities` exige `0 ≤ maxDistance` e
`0 ≤ operatorBound`. Logo superestimar `tailError` só torna a condição
`transferMajorant < 1` mais difícil de satisfazer: é o erro seguro. Subestimar
é que é fatal, porque faz o certificado passar com um controle de cauda que as
constantes de EDP não sustentam. Daí a exigência ser `data.epsilon ≤
tile.tailError`, e não a igualdade nem o contrário: **o tile pode ser
conservador, nunca otimista**.

O teorema `JustifiedTile.tailBound_transfers` transforma esse parágrafo em
dedução: dado o conteúdo de EDP `‖C - C_N‖ ≤ ε_N` como hipótese explícita, ele
conclui `‖C - C_N‖ ≤ tile.tailError` por transitividade.

## Unilateral ou bilateral: qual `ε_N` é o certo aqui

Esta ponte usa `TailData.epsilon`, o truncamento **unilateral** `C_N = Π C`,
que só consome `P2`. Isso é uma mudança em relação ao que o docstring de
`TailBound.twoSidedEpsilon` recomenda, e a razão é que a recomendação de lá
ficou obsoleta: ela descreve o `Tile.transferMajorant` antigo, de forma
`max(1, β_N)`, que era o do truncamento bilateral `Π C Π` e consumia a hipótese
`P2*` (regularidade analítica do adjunto). `Tile.transferMajorant` foi desde
então reescrito na forma unilateral `(1 + |λ| β_N ‖C‖)` exatamente para
dispensar `P2*`, e é com `TailData.epsilon` que ele casa. Ver
`docs/WINDING_DATA.md` e o §3 de `DeterminantBound`. O ponto está registrado no
relatório; não editei `TailBound.lean`, que é de outro autor.

Justificar contra `twoSidedEpsilon` continuaria **sólido**, apenas mais
pessimista: como `twoSidedEpsilon r a = r.epsilon + a.epsilon`, um `tailError`
que domina a soma domina a fortiori cada parcela. A regressão
`twoSidedIsStronger` registra isso.

## O que continua sendo hipótese

Tudo o que já era. Este módulo **não** demonstra `‖C - C_N‖ ≤ ε_N`: isso é o
conteúdo conjunto de `P1`, `P2` e `P3`, que nenhum arquivo deste repositório
prova. O que ele acrescenta é que o número usado pelo tile não pode mais ser
menor do que aquele que as constantes declaradas produzem — a hipótese fica
declarada num lugar só, em vez de aparecer duas vezes com valores diferentes.
-/

namespace Choptuik.TailTileBridge

open Choptuik.RationalCertificate
open Choptuik.TailBound

/- Mesma observação de `TailBound` e `RationalCertificate`: `unseal` é diretiva
   de elaboração, o kernel nunca respeitou irredutibilidade, e todas as
   regressões abaixo fecham por `decide`. -/
unseal Rat.mul Rat.add Rat.sub Rat.inv

/- `TailBound` mantem a sua instancia de decidibilidade `private`, mas
   `TailData.admissible` em si e publica, de modo que a instancia pode ser
   reconstruida aqui sem tocar naquele arquivo. Expo-la como `instance` deixa
   `decide` usa-la diretamente nas regressoes. -/
instance instDecidableAdmissible (data : TailData) : Decidable data.admissible := by
  unfold TailData.admissible
  infer_instance

/-! ## Um tile justificado -/

/-- Um tile do contorno **junto com** as constantes de EDP que justificam o seu
`tailError`. -/
structure JustifiedTile where
  /-- Constantes `P1`/`P2`/`P3` e os truncamentos `N_F`, `N_C`. -/
  tailData : TailData
  /-- O tile de Neumann que consome `ε_N`. -/
  tile : Tile
deriving Repr, DecidableEq

/-- Conteúdo verificado: as constantes são admissíveis, o tile passa em todas
as suas próprias desigualdades, e o `tailError` do tile domina o `ε_N` composto
a partir das constantes. -/
def JustifiedTile.claims (justified : JustifiedTile) : Prop :=
  justified.tailData.admissible ∧
  justified.tile.valid = true ∧
  justified.tailData.epsilon ≤ justified.tile.tailError

private def JustifiedTile.claimsDecidable
    (justified : JustifiedTile) : Decidable justified.claims := by
  unfold JustifiedTile.claims
  infer_instance

/-- Verificador executável e exato em `Rat`. -/
def JustifiedTile.valid (justified : JustifiedTile) : Bool :=
  @decide justified.claims justified.claimsDecidable

theorem JustifiedTile.valid_sound
    (justified : JustifiedTile) (accepted : justified.valid = true) :
    justified.claims :=
  @of_decide_eq_true justified.claims justified.claimsDecidable accepted

/-- O tile embutido satisfaz todas as desigualdades de Neumann e de cauda. -/
theorem JustifiedTile.tile_inequalities
    (justified : JustifiedTile) (accepted : justified.valid = true) :
    justified.tile.inequalities :=
  Tile.valid_sound justified.tile (justified.valid_sound accepted).2.1

/-- O `tailError` declarado é dominado pelo `ε_N` composto. -/
theorem JustifiedTile.epsilon_le_tailError
    (justified : JustifiedTile) (accepted : justified.valid = true) :
    justified.tailData.epsilon ≤ justified.tile.tailError :=
  (justified.valid_sound accepted).2.2

/-- **O que a ponte serve para concluir.** Dado o conteúdo analítico de
`P1`/`P2`/`P3` na forma `‖C - C_N‖ ≤ ε_N` — hipótese explícita, que este
repositório não demonstra —, o `tailError` usado pelo tile é de fato um
majorante de `‖C - C_N‖`.

É este teorema que torna o sentido da desigualdade uma dedução em vez de uma
convenção: ele só fecha porque a exigência é `epsilon ≤ tailError`. -/
theorem JustifiedTile.tailBound_transfers
    (justified : JustifiedTile) (accepted : justified.valid = true)
    {truncationError : Rat}
    (pde : truncationError ≤ justified.tailData.epsilon) :
    truncationError ≤ justified.tile.tailError :=
  Rat.le_trans pde (justified.epsilon_le_tailError accepted)

/-- Empacotamento: um certificado aceito entrega um tile válido cujo controle
de cauda está ancorado em constantes nomeadas, e não asseverado. -/
theorem JustifiedTile.justified_summary
    (justified : JustifiedTile) (accepted : justified.valid = true) :
    justified.tile.inequalities ∧
    justified.tailData.admissible ∧
    justified.tailData.epsilon ≤ justified.tile.tailError :=
  ⟨justified.tile_inequalities accepted,
   (justified.valid_sound accepted).1,
   justified.epsilon_le_tailError accepted⟩

/-! ## O contorno inteiro sob as mesmas constantes

O `ε_N` não é por tile: ele depende só do truncamento `(N_F, N_C)` e das
constantes de EDP, que são globais. Logo o contorno inteiro deve ser
justificado por **um único** `TailData`; é isso que a estrutura abaixo impõe, e
é a diferença entre "cada tile tem algum `ε_N`" e "todos os tiles têm o mesmo,
composto das mesmas constantes". -/

structure JustifiedBoundary where
  tailData : TailData
  boundary : BoundaryCertificate

def JustifiedBoundary.claims (justified : JustifiedBoundary) : Prop :=
  justified.tailData.admissible ∧
  justified.boundary.valid = true ∧
  ∀ tile ∈ justified.boundary.tiles,
    justified.tailData.epsilon ≤ tile.tailError

private def JustifiedBoundary.claimsDecidable
    (justified : JustifiedBoundary) : Decidable justified.claims := by
  unfold JustifiedBoundary.claims
  infer_instance

def JustifiedBoundary.valid (justified : JustifiedBoundary) : Bool :=
  @decide justified.claims justified.claimsDecidable

theorem JustifiedBoundary.valid_sound
    (justified : JustifiedBoundary) (accepted : justified.valid = true) :
    justified.claims :=
  @of_decide_eq_true justified.claims justified.claimsDecidable accepted

/-- Cada tile do contorno, visto como tile justificado pelas constantes
globais. -/
def JustifiedBoundary.tileAt
    (justified : JustifiedBoundary) (tile : Tile) : JustifiedTile where
  tailData := justified.tailData
  tile := tile

/-- **A ponte, estendida ao contorno.** Se o contorno é aceito, então cada um
dos seus tiles, emparelhado com as constantes globais, forma um `JustifiedTile`
aceito. Nenhum tile escapa da ancoragem. -/
theorem JustifiedBoundary.tileAt_valid
    (justified : JustifiedBoundary) (accepted : justified.valid = true) :
    ∀ tile ∈ justified.boundary.tiles, (justified.tileAt tile).valid = true := by
  intro tile member
  have verified := justified.valid_sound accepted
  have tileValid : tile.valid = true :=
    List.all_eq_true.mp verified.2.1 tile member
  have dominates : justified.tailData.epsilon ≤ tile.tailError :=
    verified.2.2 tile member
  exact @decide_eq_true _ (JustifiedTile.claimsDecidable _)
    ⟨verified.1, tileValid, dominates⟩

/-- Consequência direta: todo tile do contorno satisfaz as suas desigualdades
**e** tem o controle de cauda ancorado nas mesmas constantes. -/
theorem JustifiedBoundary.every_tile_justified
    (justified : JustifiedBoundary) (accepted : justified.valid = true) :
    ∀ tile ∈ justified.boundary.tiles,
      tile.inequalities ∧ justified.tailData.epsilon ≤ tile.tailError := by
  intro tile member
  have verified := justified.valid_sound accepted
  exact ⟨Tile.valid_sound tile (List.all_eq_true.mp verified.2.1 tile member),
         verified.2.2 tile member⟩

/-! ## Regressões, todas fechadas pelo kernel

As constantes são as de `H-opt` de `docs/TAIL_BOUND.md` na malha `12×24`, que é
a menor malha que fecha o alvo de Neumann naquele módulo. O valor composto é

```
ε_N = 40 · ((2/5)^13 + (2/5)^25) = 16000268435456/59604644775390625 ≈ 2,684·10⁻⁴ .
```

O tile é o escalar de `RationalCertificate`: `F₀ = [2]`, `V = [1/2]`, raio
`1/10`, `‖C_N‖ = 1`, `|z-z₀| = 3`, `‖C‖ = 1`. Com `inverseBound = 10/19`, a
condição `transferMajorant < 1` só exige `tailError < 19/147 ≈ 0,129`, o que é
**muito** mais frouxo que `ε_N`: é precisamente essa folga que deixava passar
um `tailError` fisicamente injustificado. -/

def bridgeGrid : TailData := optimisticGrid 12 24

example : bridgeGrid.epsilon = 16000268435456 / 59604644775390625 := by decide

/-- Tile cujo `tailError = 10⁻³` domina o `ε_N ≈ 2,684·10⁻⁴` composto. -/
def justifiedScalarTile : Tile where
  centerMatrix := scalarMatrix 2
  approxInverse := scalarMatrix (1 / 2)
  tileRadius := 1 / 10
  matrixNorm := 1
  maxDistance := 3
  tailError := 1 / 1000
  operatorBound := 1

def acceptedBridge : JustifiedTile where
  tailData := bridgeGrid
  tile := justifiedScalarTile

example : justifiedScalarTile.valid = true := by decide
example : justifiedScalarTile.transferMajorant = 147 / 19000 := by decide
example : acceptedBridge.valid = true := by decide

/-- **Rejeição 1: o defeito que a ponte existe para pegar.**

Idêntico ao tile aceito, exceto por `tailError = 10⁻⁴`, que é **menor** que o
`ε_N ≈ 2,684·10⁻⁴` produzido pelas próprias constantes declaradas. O tile
continua perfeitamente válido isoladamente — `transferMajorant` fica até
*menor*, `147/190000` —, e era exatamente assim que o número passava
despercebido: quem olha só o tile vê um certificado mais folgado, não um
certificado quebrado. A ponte recusa. -/
def underJustifiedTile : Tile :=
  { justifiedScalarTile with tailError := 1 / 10000 }

def rejectedTooSmall : JustifiedTile where
  tailData := bridgeGrid
  tile := underJustifiedTile

-- O tile sozinho passa, e passa com folga maior que o aceito:
example : underJustifiedTile.valid = true := by decide
example : underJustifiedTile.transferMajorant = 147 / 190000 := by decide
-- Mas o `tailError` declarado não é sustentado pelas constantes:
example : ¬ (bridgeGrid.epsilon ≤ underJustifiedTile.tailError) := by decide
example : rejectedTooSmall.valid = false := by decide

/-- **Rejeição 2: `TailData` inadmissível, com o defeito escondido da
aritmética.**

`decayChebyshev = 1` não é decaimento: a série de Chebyshev não converge e
`q^(N_C+1)` não é cauda de nada. O detalhe que torna esta regressão útil é que
o produtor zerou `chebyshevTailFactor`, de modo que a **aritmética** não
denuncia nada — `ε_N` fica até menor que o do caso aceito, e a desigualdade
`ε_N ≤ tailError` continua valendo. Só a checagem de admissibilidade pega. É a
razão de `admissible` ser exigida separadamente em vez de deduzida do valor
composto. -/
def inadmissibleGrid : TailData :=
  { bridgeGrid with decayChebyshev := 1, chebyshevTailFactor := 0 }

def rejectedInadmissible : JustifiedTile where
  tailData := inadmissibleGrid
  tile := justifiedScalarTile

-- A aritmética não denuncia: o `ε_N` até diminui e a desigualdade vale.
example : inadmissibleGrid.epsilon = 65536 / 244140625 := by decide
example : inadmissibleGrid.epsilon ≤ justifiedScalarTile.tailError := by decide
example : justifiedScalarTile.valid = true := by decide
-- E ainda assim é recusado, pela admissibilidade:
example : ¬ inadmissibleGrid.admissible := by decide
example : rejectedInadmissible.valid = false := by decide

/-- **Rejeição 3: constantes boas, tile ruim.** As duas metades são exigidas:
aqui o `ε_N` é dominado com folga, mas a cauda declarada é grande demais para o
majorante de transferência do próprio tile. -/
def rejectedBadTile : JustifiedTile where
  tailData := bridgeGrid
  tile := { justifiedScalarTile with tailError := 1 }

example : bridgeGrid.epsilon ≤ (1 : Rat) := by decide
example : ({ justifiedScalarTile with tailError := 1 } : Tile).valid = false := by decide
example : rejectedBadTile.valid = false := by decide

/-- Justificar contra o `ε_N` bilateral é estritamente mais forte: a soma
domina cada parcela. Um `tailError` que passe pelo bilateral passa a fortiori
por esta ponte. -/
example : bridgeGrid.epsilon ≤ twoSidedEpsilon bridgeGrid bridgeGrid := by decide

example : twoSidedEpsilon bridgeGrid bridgeGrid ≤ justifiedScalarTile.tailError := by
  decide

/-! ### O contorno inteiro -/

def acceptedBoundary : JustifiedBoundary where
  tailData := bridgeGrid
  boundary :=
    { tiles := [justifiedScalarTile, justifiedScalarTile, justifiedScalarTile]
      nonempty := by decide }

example : acceptedBoundary.valid = true := by decide

/-- Um único tile mal ancorado no meio de tiles bons derruba o contorno
inteiro — que é o comportamento certo, porque `ε_N` é global. -/
def rejectedBoundary : JustifiedBoundary where
  tailData := bridgeGrid
  boundary :=
    { tiles := [justifiedScalarTile, underJustifiedTile, justifiedScalarTile]
      nonempty := by decide }

-- O contorno passa no verificador antigo, que só olha os tiles:
example : rejectedBoundary.boundary.valid = true := by decide
-- E é recusado pela ponte, por causa do tile do meio:
example : rejectedBoundary.valid = false := by decide

#print axioms JustifiedTile.valid_sound
#print axioms JustifiedTile.tile_inequalities
#print axioms JustifiedTile.epsilon_le_tailError
#print axioms JustifiedTile.tailBound_transfers
#print axioms JustifiedTile.justified_summary
#print axioms JustifiedBoundary.valid_sound
#print axioms JustifiedBoundary.tileAt_valid
#print axioms JustifiedBoundary.every_tile_justified

end Choptuik.TailTileBridge
