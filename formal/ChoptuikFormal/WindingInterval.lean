import Std
import ChoptuikFormal.RationalCertificate

/-!
# Princípio do argumento discreto sobre caixas racionais

Este módulo fecha o buraco registrado no ledger (obrigações `C3`/`C4`). Até
aqui `Choptuik.RationalCertificate.WindingCertificate` apenas somava
incrementos `-1/0/+1` entregues por um programa externo: nada obrigava um
incremento a vir de um enclosure verificado, de modo que a contagem final
confiava num oráculo. Aqui o incremento de cada aresta é *calculado* pelo
kernel do Lean a partir de caixas racionais exatas, e a ponte
`IntervalWindingCertificate.toWinding` mostra que o certificado inteiro antigo
é consequência do certificado intervalar.

Não há ponto flutuante em lugar nenhum: todas as coordenadas são `Rat` e todas
as comparações são decididas pelo kernel.

## O que continua sendo hipótese do produtor externo

Seja `Γ` o contorno poligonal fechado `z₀ → z₁ → ⋯ → z_{n-1} → z_n = z₀` e
`f : Γ → ℂ` contínua (no caso de interesse `f(z) = det[I + (z-z₀)C_N]`). Para a
aresta `k` o produtor declara três caixas e afirma:

* `startBox` contém `f(z_k)`;
* `endBox`   contém `f(z_{k+1})`;
* `edgeBox`  contém `f(z)` para **todo** `z` no segmento `[z_k, z_{k+1}]`.

Essas três inclusões são a parte analítica e **não** são verificadas aqui (nem
poderiam ser: `f` não aparece neste módulo). Verificado aqui está todo o resto:
que as caixas são bem formadas; que `edgeBox` contém `startBox` e `endBox`
(consequência necessária da afirmação acima, logo um teste de consistência
barato mas real); que a lista encadeia e fecha; que cada aresta admite
classificação; e que a soma dos incrementos é o alvo declarado.

## Por que os três ramos determinam o incremento

Levante o argumento de `f` continuamente ao longo de `Γ`, obtendo `θ(t)` com
`w = (θ(1) - θ(0)) / 2π` igual ao winding em torno de `0`. Ponha
`m_k = ⌊(θ_k - π) / 2π⌋` em cada nó. Como `θ_n = θ₀ + 2πw` e `⌊x + w⌋ = ⌊x⌋ + w`
para `w` inteiro, a soma telescópica `Σ_k (m_{k+1} - m_k)` vale exatamente `w`.
O incremento da aresta `k` é portanto `m_{k+1} - m_k`, e cada ramo de
`edgeCrossing` o determina:

* `edgeBox ⊆ {Re > 0}`: ao longo da aresta `θ mod 2π ∈ (-π/2, π/2)`, logo
  `(θ - π)/2π` nunca é inteiro; por continuidade `⌊(θ-π)/2π⌋` é constante na
  aresta e o incremento é `0`.
* `edgeBox ⊆ {Im > 0}` ou `edgeBox ⊆ {Im < 0}`: então `θ mod 2π ∈ (0, π)` ou
  `θ mod 2π ∈ (-π, 0)`, e de novo `(θ-π)/2π` nunca é inteiro: incremento `0`.
* `edgeBox ⊆ {Re < 0}`: então `θ mod 2π ∈ (π/2, 3π/2)`; como esse conjunto é
  união disjunta de intervalos abertos de comprimento `π`, existe um único
  inteiro `m` com `θ(t) ∈ (π/2 + 2πm, 3π/2 + 2πm)` em toda a aresta. Dentro
  desse intervalo vale `⌊(θ-π)/2π⌋ = m - 1` exatamente quando `Im f > 0`
  (segundo quadrante) e `⌊(θ-π)/2π⌋ = m` quando `Im f ≤ 0` (terceiro
  quadrante). Portanto `Im: + → -` dá `m - (m-1) = +1`, `Im: - → +` dá `-1`, e
  **os dois nós com o mesmo sinal estrito de `Im` dão `0`**: `(m-1) - (m-1)` ou
  `m - m`. Como `m` é constante na aresta inteira, esse `0` vale mesmo que
  `Im f` oscile várias vezes dentro dela — o critério só depende dos extremos.

Orientação, conferida de forma independente: para `f(z) = z` no círculo
`z = e^{iθ}` com `θ : 0 → 2π` (sentido anti-horário, winding `+1`), o único
cruzamento do eixo real negativo acontece em `θ = π`, onde `Im` passa de
positivo para negativo. Logo `Im: + → -` contribui `+1`, exatamente como no
terceiro ramo acima.

Fora desses três ramos a aresta é recusada (`none`) e o certificado inteiro é
rejeitado: o contorno precisa ser refinado. Cada um dos três ramos implica
`0 ∉ edgeBox` (teorema `WindingEdge.zero_not_mem_edgeBox`), que é o que garante
que `f` não se anula sobre `Γ` — sem isso a contagem não significaria nada.

## Encadeamento e fechamento

O telescopamento só vale se as arestas formarem um ciclo. Lean confere a
igualdade das *caixas* nos nós compartilhados: `endBox` da aresta `k` igual a
`startBox` da aresta `k+1`, e `endBox` da última igual a `startBox` da primeira.
Isso é necessário (caixas distintas exibiriam uma lista que não descreve um
ciclo) e é o bastante para o argumento combinatório, pois garante que os dois
lados de cada nó usam exatamente os mesmos sinais ao decidir `⌊(θ-π)/2π⌋`. A
identidade dos pontos `z_k` em si continua a cargo do produtor, junto com as
três inclusões analíticas.

## Nota de confiança

Todas as regressões deste arquivo usam `decide`, isto é, são recomputadas pelo
kernel. Nenhuma usa avaliação nativa, que acrescentaria `Lean.ofReduceBool` e
passaria a confiar no compilador fora do kernel.

**Correção de uma nota anterior.** Uma versão deste cabeçalho afirmava que a
forma `n / d` "trava a redução no kernel" e que os racionais não inteiros
*precisam* ser escritos como `mkRat n d`. Isso está errado, e a origem do erro
era ter medido sem `unseal`. `Rat.mul`, `Rat.add`, `Rat.sub` e `Rat.inv` são
`@[irreducible]` no Lean 4.33; isso é uma diretiva de **elaboração**, e o
kernel nunca a respeitou. Com `unseal Rat.mul Rat.add Rat.sub Rat.inv` — como
já se faz em `RationalCertificate`, `TailBound` e `DeterminantBound` — as
formas `n / d`, `+`, `-`, `*` e `^` reduzem normalmente. A restrição não é do
Lean e não deve ser propagada como se fosse.

O que continua verdadeiro: `mkRat` é inofensivo e dispensa `unseal`, porque
este módulo só **compara** racionais e nunca os soma. O produtor
`scripts/make_winding_certificate.py` pode seguir emitindo `mkRat`, mas por
conveniência, não por necessidade. A seção de regressões traz `unseal` e uma
regressão escrita com `/` e com aritmética `centro ± raio`, que fecha por
`decide` — a correção acima está verificada no arquivo, não apenas afirmada.
-/

namespace Choptuik.WindingInterval

/-- Auxiliar puramente logístico: `decide` com uma instância explícita. É
preciso porque a instância usada por `WindingCertificate.valid` é `private` no
módulo `RationalCertificate` e não pode ser nomeada aqui. -/
private theorem decide_true_with {claim : Prop} (instance_ : Decidable claim)
    (proof : claim) : @decide claim instance_ = true :=
  @decide_eq_true claim instance_ proof

/-! ## Caixas racionais complexas -/

/-- Retângulo fechado `[reLo, reHi] × [imLo, imHi] ⊆ ℂ` com cantos racionais
exatos. -/
structure ComplexBox where
  reLo : Rat
  reHi : Rat
  imLo : Rat
  imHi : Rat
deriving Repr, DecidableEq

/-- Pertinência do ponto `x + i y` à caixa. Proposição, usada só nos
enunciados de solidez. -/
def ComplexBox.mem (box : ComplexBox) (x y : Rat) : Prop :=
  box.reLo ≤ x ∧ x ≤ box.reHi ∧ box.imLo ≤ y ∧ y ≤ box.imHi

/-- Uma caixa é bem formada quando os cantos estão na ordem certa. -/
def ComplexBox.wellFormed (box : ComplexBox) : Bool :=
  box.reLo ≤ box.reHi && box.imLo ≤ box.imHi

/-- `box ⊆ {Re > 0}`. -/
def ComplexBox.strictlyPosRe (box : ComplexBox) : Bool := 0 < box.reLo

/-- `box ⊆ {Re < 0}`. -/
def ComplexBox.strictlyNegRe (box : ComplexBox) : Bool := box.reHi < 0

/-- `box ⊆ {Im > 0}`. -/
def ComplexBox.strictlyPosIm (box : ComplexBox) : Bool := 0 < box.imLo

/-- `box ⊆ {Im < 0}`. -/
def ComplexBox.strictlyNegIm (box : ComplexBox) : Bool := box.imHi < 0

/-- `outer.contains inner` decide a inclusão de retângulos `inner ⊆ outer`. -/
def ComplexBox.contains (outer inner : ComplexBox) : Bool :=
  outer.reLo ≤ inner.reLo && inner.reHi ≤ outer.reHi &&
  outer.imLo ≤ inner.imLo && inner.imHi ≤ outer.imHi

theorem ComplexBox.le_of_wellFormed (box : ComplexBox)
    (accepted : box.wellFormed = true) :
    box.reLo ≤ box.reHi ∧ box.imLo ≤ box.imHi := by
  simp only [ComplexBox.wellFormed, Bool.and_eq_true, decide_eq_true_eq] at accepted
  exact accepted

/-- A inclusão decidida por `contains` é mesmo a inclusão de conjuntos. -/
theorem ComplexBox.mem_of_contains {outer inner : ComplexBox}
    (inclusion : outer.contains inner = true) {x y : Rat}
    (inside : inner.mem x y) : outer.mem x y := by
  simp only [ComplexBox.contains, Bool.and_eq_true, decide_eq_true_eq] at inclusion
  obtain ⟨⟨⟨reLo, reHi⟩, imLo⟩, imHi⟩ := inclusion
  obtain ⟨xLo, xHi, yLo, yHi⟩ := inside
  exact ⟨Rat.le_trans reLo xLo, Rat.le_trans xHi reHi,
         Rat.le_trans imLo yLo, Rat.le_trans yHi imHi⟩

theorem ComplexBox.not_mem_zero_of_strictlyPosRe (box : ComplexBox)
    (strict : box.strictlyPosRe = true) : ¬ box.mem 0 0 := by
  intro inside
  have positive : (0 : Rat) < box.reLo := by
    unfold ComplexBox.strictlyPosRe at strict
    exact of_decide_eq_true strict
  exact absurd positive (Rat.not_lt.mpr inside.1)

theorem ComplexBox.not_mem_zero_of_strictlyNegRe (box : ComplexBox)
    (strict : box.strictlyNegRe = true) : ¬ box.mem 0 0 := by
  intro inside
  have negative : box.reHi < 0 := by
    unfold ComplexBox.strictlyNegRe at strict
    exact of_decide_eq_true strict
  exact absurd inside.2.1 (Rat.not_le.mpr negative)

theorem ComplexBox.not_mem_zero_of_strictlyPosIm (box : ComplexBox)
    (strict : box.strictlyPosIm = true) : ¬ box.mem 0 0 := by
  intro inside
  have positive : (0 : Rat) < box.imLo := by
    unfold ComplexBox.strictlyPosIm at strict
    exact of_decide_eq_true strict
  exact absurd positive (Rat.not_lt.mpr inside.2.2.1)

theorem ComplexBox.not_mem_zero_of_strictlyNegIm (box : ComplexBox)
    (strict : box.strictlyNegIm = true) : ¬ box.mem 0 0 := by
  intro inside
  have negative : box.imHi < 0 := by
    unfold ComplexBox.strictlyNegIm at strict
    exact of_decide_eq_true strict
  exact absurd inside.2.2.2 (Rat.not_le.mpr negative)

/-! ## Pontos do plano `z` -/

/-- Um ponto racional exato do plano `z`, isto é, um **vértice do contorno**
`Γ`. Não confundir com `ComplexBox`: as caixas vivem no plano de *valores* de
`f`, os pontos vivem no plano do *argumento* `z`. -/
structure ComplexPoint where
  re : Rat
  im : Rat
deriving Repr, DecidableEq

/-! ## Arestas do contorno -/

/-- Uma aresta do contorno: o segmento `[z_k, z_{k+1}]` no plano `z`, com o
enclosure dos valores de `f` sobre ele.

`startPoint` e `endPoint` são os dois vértices, **exatos**; é sobre eles que o
kernel verifica que a lista de arestas descreve mesmo uma curva fechada (ver
`chained`, `closed` e `WindingEdge.nondegenerate`).

Semântica pretendida das três caixas, que é a hipótese analítica a cargo do
produtor externo e que o kernel **não** pode verificar (ver o cabeçalho do
módulo): `startBox ∋ f(startPoint)`, `endBox ∋ f(endPoint)` e
`edgeBox ∋ f(z)` para todo `z` no segmento `[startPoint, endPoint]`. Em
particular `edgeBox` tem de conter `startBox` e `endBox`; essa consequência é
verificada aqui por `WindingEdge.shapeOk`. -/
structure WindingEdge where
  startPoint : ComplexPoint
  endPoint : ComplexPoint
  startBox : ComplexBox
  endBox : ComplexBox
  edgeBox : ComplexBox
deriving Repr, DecidableEq

/-- Checagens de forma: as três caixas são bem formadas e `edgeBox` contém as
caixas dos dois nós. -/
def WindingEdge.shapeOk (edge : WindingEdge) : Bool :=
  edge.startBox.wellFormed && edge.endBox.wellFormed && edge.edgeBox.wellFormed &&
  edge.edgeBox.contains edge.startBox && edge.edgeBox.contains edge.endBox

/-- Uma aresta de comprimento zero não é aresta: ela não contribui para o
winding e ainda assim ocuparia uma posição na cadeia, de modo que um produtor
poderia usá-la para satisfazer o encadeamento sem descrever curva nenhuma. -/
def WindingEdge.nondegenerate (edge : WindingEdge) : Bool :=
  edge.startPoint != edge.endPoint

/-- Classificação exata do cruzamento do raio `ℝ_{<0}` pela aresta.

`some 0` quando a aresta comprovadamente não tem cruzamento líquido do eixo
real negativo; `some 1` / `some (-1)` quando ela o cruza uma vez em cada
sentido, com o sinal fixado pela orientação (`Im : + → -` vale `+1`, ver o
cabeçalho do módulo); `none` quando as caixas fornecidas não decidem nada e o
contorno precisa ser refinado.

O quarto caso do semiplano esquerdo — os **dois** nós com o mesmo sinal estrito
de `Im` — também dá `0`, e isso é demonstrável pelo mesmo argumento do piso: em
`edgeBox ⊆ {Re < 0}` o inteiro `m` é constante na aresta (conexidade), e o piso
vale `m-1` em todo ponto com `Im f > 0` e `m` em todo ponto com `Im f < 0`.
Dois nós com `Im > 0` dão `(m-1) - (m-1) = 0`; dois com `Im < 0` dão
`m - m = 0`. Note que isso vale mesmo que `Im f` oscile **várias vezes** dentro
da aresta: o critério só olha os dois extremos, e `m` não muda no meio. É o
caso que aparece na vizinhança do cruzamento forçado do eixo real
(`docs/WINDING_DATA.md` §2.2), onde as arestas ficam inteiras no semiplano
esquerdo; devolver `none` ali forçava refinamento sem necessidade. -/
def WindingEdge.edgeCrossing (edge : WindingEdge) : Option Int :=
  if edge.edgeBox.strictlyPosRe then
    -- Semiplano direito: não encontra o eixo real negativo.
    some 0
  else if edge.edgeBox.strictlyPosIm || edge.edgeBox.strictlyNegIm then
    -- Semiplano superior ou inferior aberto: não encontra o eixo real.
    some 0
  else if edge.edgeBox.strictlyNegRe then
    -- Semiplano esquerdo: o cruzamento é decidido pelos sinais de Im nos nós.
    -- Os dois casos de MESMO sinal vêm primeiro de propósito: assim
    -- `edgeCrossing_eq_zero_of_sameSignIm` vale sem hipótese de boa formação, e
    -- uma caixa malformada (que teria `imLo > 0` e `imHi < 0` ao mesmo tempo)
    -- produz `0` em vez de `±1` — o erro conservador, que não infla a contagem.
    -- Caixas malformadas são rejeitadas por `shapeOk` de todo modo.
    if edge.startBox.strictlyPosIm && edge.endBox.strictlyPosIm then
      some 0
    else if edge.startBox.strictlyNegIm && edge.endBox.strictlyNegIm then
      some 0
    else if edge.startBox.strictlyPosIm && edge.endBox.strictlyNegIm then
      some 1
    else if edge.startBox.strictlyNegIm && edge.endBox.strictlyPosIm then
      some (-1)
    else
      none
  else
    none

/-- A aresta foi classificada por algum dos três ramos. -/
def WindingEdge.classified (edge : WindingEdge) : Bool :=
  edge.edgeCrossing.isSome

/-- Incremento da aresta. Vale `0` também no caso não classificado; por isso
todo certificado exige `classified` separadamente. -/
def WindingEdge.increment (edge : WindingEdge) : Int :=
  edge.edgeCrossing.getD 0

/-- O incremento é sempre um dos três valores admitidos pelo certificado
inteiro de `RationalCertificate`. -/
theorem WindingEdge.increment_admissible (edge : WindingEdge) :
    RationalCertificate.WindingCertificate.incrementAdmissible edge.increment := by
  unfold RationalCertificate.WindingCertificate.incrementAdmissible
  unfold WindingEdge.increment WindingEdge.edgeCrossing
  repeat' split
  all_goals decide

/-- **Quarto caso do semiplano esquerdo.** Com `edgeBox ⊆ {Re < 0}` e os dois
nós com o mesmo sinal estrito de `Im`, o incremento é `0`.

Este é o conteúdo do argumento do piso recolhido no cabeçalho do módulo: no
semiplano esquerdo o inteiro `m` é constante ao longo da aresta, e o piso vale
`m-1` onde `Im f > 0` e `m` onde `Im f < 0`; dois nós do mesmo lado dão
diferença nula. O enunciado não supõe boa formação das caixas, e é por isso que
os dois casos de mesmo sinal são testados **antes** dos casos de sinal oposto na
definição de `edgeCrossing`. -/
theorem WindingEdge.edgeCrossing_eq_zero_of_sameSignIm (edge : WindingEdge)
    (leftHalf : edge.edgeBox.strictlyNegRe = true)
    (sameSign :
      (edge.startBox.strictlyPosIm && edge.endBox.strictlyPosIm) = true ∨
      (edge.startBox.strictlyNegIm && edge.endBox.strictlyNegIm) = true) :
    edge.edgeCrossing = some 0 := by
  cases positiveRe : edge.edgeBox.strictlyPosRe with
  | true => simp [WindingEdge.edgeCrossing, positiveRe]
  | false =>
    cases horizontal : (edge.edgeBox.strictlyPosIm || edge.edgeBox.strictlyNegIm) with
    | true => simp [WindingEdge.edgeCrossing, positiveRe, horizontal]
    | false =>
      cases bothPositive :
          (edge.startBox.strictlyPosIm && edge.endBox.strictlyPosIm) with
      | true =>
        simp [WindingEdge.edgeCrossing, positiveRe, horizontal, bothPositive,
          leftHalf]
      | false =>
        -- Se não são os dois com `Im > 0`, a hipótese força os dois com
        -- `Im < 0`, e esse é o segundo caso testado.
        have bothNegative :
            (edge.startBox.strictlyNegIm && edge.endBox.strictlyNegIm) = true := by
          rcases sameSign with positives | negatives
          · simp [bothPositive] at positives
          · exact negatives
        simp [WindingEdge.edgeCrossing, positiveRe, horizontal, bothPositive,
          bothNegative, leftHalf]

/-- Corolário imediato: nesse caso o incremento somado é `0`. -/
theorem WindingEdge.increment_eq_zero_of_sameSignIm (edge : WindingEdge)
    (leftHalf : edge.edgeBox.strictlyNegRe = true)
    (sameSign :
      (edge.startBox.strictlyPosIm && edge.endBox.strictlyPosIm) = true ∨
      (edge.startBox.strictlyNegIm && edge.endBox.strictlyNegIm) = true) :
    edge.increment = 0 := by
  unfold WindingEdge.increment
  rw [edge.edgeCrossing_eq_zero_of_sameSignIm leftHalf sameSign]
  rfl

/-- **Não anulamento sobre o contorno.** Qualquer um dos três ramos de
`edgeCrossing` exclui a origem de `edgeBox`. Como `edgeBox` deve conter todos
os valores de `f` sobre a aresta, é este teorema que garante que `f` não se
anula sobre `Γ`. -/
theorem WindingEdge.zero_not_mem_edgeBox (edge : WindingEdge) {increment : Int}
    (classified : edge.edgeCrossing = some increment) : ¬ edge.edgeBox.mem 0 0 := by
  unfold WindingEdge.edgeCrossing at classified
  split at classified
  · next strict => exact ComplexBox.not_mem_zero_of_strictlyPosRe _ strict
  · split at classified
    · next strict =>
        rcases (Bool.or_eq_true _ _).mp strict with positive | negative
        · exact ComplexBox.not_mem_zero_of_strictlyPosIm _ positive
        · exact ComplexBox.not_mem_zero_of_strictlyNegIm _ negative
    · split at classified
      · next strict => exact ComplexBox.not_mem_zero_of_strictlyNegRe _ strict
      · exact absurd classified (by simp)

theorem WindingEdge.zero_not_mem_startBox (edge : WindingEdge) {increment : Int}
    (shape : edge.shapeOk = true) (classified : edge.edgeCrossing = some increment) :
    ¬ edge.startBox.mem 0 0 := by
  intro inside
  simp only [WindingEdge.shapeOk, Bool.and_eq_true] at shape
  exact edge.zero_not_mem_edgeBox classified
    (ComplexBox.mem_of_contains shape.1.2 inside)

theorem WindingEdge.zero_not_mem_endBox (edge : WindingEdge) {increment : Int}
    (shape : edge.shapeOk = true) (classified : edge.edgeCrossing = some increment) :
    ¬ edge.endBox.mem 0 0 := by
  intro inside
  simp only [WindingEdge.shapeOk, Bool.and_eq_true] at shape
  exact edge.zero_not_mem_edgeBox classified
    (ComplexBox.mem_of_contains shape.2 inside)

/-! ## Certificado intervalar do winding -/

/-- Arestas consecutivas casam, **no plano `z` e no plano de valores**.

As duas igualdades têm papéis diferentes e nenhuma das duas dispensa a outra:

* `endPoint = startPoint` é o que faz a lista descrever uma **curva**. Sem
  ela, um produtor com geometria inconsistente (nós repetidos, arestas fora de
  ordem, contorno que não fecha no plano `z`) emitiria caixas que encadeiam
  perfeitamente sem que exista contorno nenhum, e a soma dos incrementos não
  significaria coisa alguma.
* `endBox = startBox` é o que garante que os dois lados de cada nó leem os
  **mesmos sinais** ao decidir `⌊(θ-π)/2π⌋`. Ela é o que torna a soma
  telescópica consistente, e é mais forte do que parece: impede, por exemplo,
  um ciclo de duas arestas com incrementos `+1, +1`, que exigiria que a mesma
  caixa tivesse `imLo > 0` e `imHi < 0`. -/
def chained : List WindingEdge → Bool
  | [] => true
  | [_] => true
  | first :: second :: rest =>
      (first.endPoint == second.startPoint) &&
      (first.endBox == second.startBox) &&
      chained (second :: rest)

/-- O contorno fecha: a última aresta volta ao primeiro vértice, e a caixa
correspondente também coincide. A lista vazia não fecha nada e é rejeitada. -/
def closed (edges : List WindingEdge) : Bool :=
  match edges.head?, edges.getLast? with
  | some first, some last =>
      (last.endPoint == first.startPoint) && (last.endBox == first.startBox)
  | _, _ => false

/-- Certificado intervalar completo de um winding number. -/
structure IntervalWindingCertificate where
  edges : List WindingEdge
  target : Int
deriving Repr, DecidableEq

def IntervalWindingCertificate.increments
    (certificate : IntervalWindingCertificate) : List Int :=
  certificate.edges.map WindingEdge.increment

def IntervalWindingCertificate.total
    (certificate : IntervalWindingCertificate) : Int :=
  certificate.increments.foldl (fun total step => total + step) 0

/-- Conteúdo matemático verificado para um contorno inteiro. -/
def IntervalWindingCertificate.claims
    (certificate : IntervalWindingCertificate) : Prop :=
  certificate.edges ≠ [] ∧
  (∀ edge ∈ certificate.edges, edge.shapeOk = true) ∧
  (∀ edge ∈ certificate.edges, edge.nondegenerate = true) ∧
  (∀ edge ∈ certificate.edges, edge.classified = true) ∧
  chained certificate.edges = true ∧
  closed certificate.edges = true ∧
  certificate.total = certificate.target

private def IntervalWindingCertificate.claimsDecidable
    (certificate : IntervalWindingCertificate) : Decidable certificate.claims := by
  unfold IntervalWindingCertificate.claims
  infer_instance

/-- Verificador executável e exato do certificado intervalar. -/
def IntervalWindingCertificate.valid
    (certificate : IntervalWindingCertificate) : Bool :=
  @decide certificate.claims certificate.claimsDecidable

/-- O verificador não aceita um certificado sem que todas as afirmações sejam
verdadeiras. -/
theorem IntervalWindingCertificate.valid_sound
    (certificate : IntervalWindingCertificate)
    (accepted : certificate.valid = true) :
    certificate.claims :=
  @of_decide_eq_true certificate.claims certificate.claimsDecidable accepted

theorem IntervalWindingCertificate.total_eq_target_of_valid
    (certificate : IntervalWindingCertificate)
    (accepted : certificate.valid = true) :
    certificate.total = certificate.target :=
  (certificate.valid_sound accepted).2.2.2.2.2.2

/-- **A lista de arestas descreve mesmo um ciclo fechado no plano `z`.**
Esta é a afirmação que o certificado antigo não fazia: ali só as caixas eram
comparadas, e uma geometria inconsistente passava despercebida. -/
theorem IntervalWindingCertificate.geometryClosed
    (certificate : IntervalWindingCertificate)
    (accepted : certificate.valid = true) :
    chained certificate.edges = true ∧ closed certificate.edges = true :=
  ⟨(certificate.valid_sound accepted).2.2.2.2.1,
   (certificate.valid_sound accepted).2.2.2.2.2.1⟩

/-- **Nenhuma aresta é degenerada.** -/
theorem IntervalWindingCertificate.nondegenerate_edges
    (certificate : IntervalWindingCertificate)
    (accepted : certificate.valid = true) :
    ∀ edge ∈ certificate.edges, edge.startPoint ≠ edge.endPoint := by
  intro edge member
  have distinct : edge.nondegenerate = true :=
    (certificate.valid_sound accepted).2.2.1 edge member
  unfold WindingEdge.nondegenerate at distinct
  exact bne_iff_ne.mp distinct

/-- **`f` não se anula em nenhuma aresta do contorno aceito.** -/
theorem IntervalWindingCertificate.zero_not_mem_edgeBoxes
    (certificate : IntervalWindingCertificate)
    (accepted : certificate.valid = true) :
    ∀ edge ∈ certificate.edges, ¬ edge.edgeBox.mem 0 0 := by
  intro edge member
  have isClassified : edge.classified = true :=
    (certificate.valid_sound accepted).2.2.2.1 edge member
  cases crossing : edge.edgeCrossing with
  | none =>
      unfold WindingEdge.classified at isClassified
      rw [crossing] at isClassified
      exact absurd isClassified (by decide)
  | some increment => exact edge.zero_not_mem_edgeBox crossing

/-- Versão nos nós: nem `f(z_k)` nem `f(z_{k+1})` se anulam. -/
theorem IntervalWindingCertificate.zero_not_mem_nodeBoxes
    (certificate : IntervalWindingCertificate)
    (accepted : certificate.valid = true) :
    ∀ edge ∈ certificate.edges,
      ¬ edge.startBox.mem 0 0 ∧ ¬ edge.endBox.mem 0 0 := by
  intro edge member
  have shape : edge.shapeOk = true :=
    (certificate.valid_sound accepted).2.1 edge member
  have isClassified : edge.classified = true :=
    (certificate.valid_sound accepted).2.2.2.1 edge member
  cases crossing : edge.edgeCrossing with
  | none =>
      unfold WindingEdge.classified at isClassified
      rw [crossing] at isClassified
      exact absurd isClassified (by decide)
  | some increment =>
      exact ⟨edge.zero_not_mem_startBox shape crossing,
             edge.zero_not_mem_endBox shape crossing⟩

/-! ## A ponte para o certificado inteiro existente

Esta é a parte que fecha o buraco: o `WindingCertificate` de
`RationalCertificate` deixa de receber incrementos de um oráculo e passa a
recebê-los de enclosures verificados. -/

/-- Extrai o certificado inteiro correspondente. Os incrementos não são
fornecidos: são computados a partir das caixas. -/
def IntervalWindingCertificate.toWinding
    (certificate : IntervalWindingCertificate) :
    RationalCertificate.WindingCertificate where
  increments := certificate.increments
  target := certificate.target

theorem IntervalWindingCertificate.toWinding_increments
    (certificate : IntervalWindingCertificate) :
    certificate.toWinding.increments = certificate.increments := rfl

theorem IntervalWindingCertificate.toWinding_total
    (certificate : IntervalWindingCertificate) :
    certificate.toWinding.total = certificate.total := rfl

/-- **Ponte.** Um certificado intervalar válido produz um `WindingCertificate`
válido. -/
theorem IntervalWindingCertificate.toWinding_valid
    (certificate : IntervalWindingCertificate)
    (accepted : certificate.valid = true) :
    certificate.toWinding.valid = true := by
  have verified := certificate.valid_sound accepted
  have nonempty : certificate.toWinding.increments ≠ [] := by
    show certificate.edges.map WindingEdge.increment ≠ []
    cases listing : certificate.edges with
    | nil => exact absurd listing verified.1
    | cons _ _ => simp
  have admissible : ∀ increment ∈ certificate.toWinding.increments,
      RationalCertificate.WindingCertificate.incrementAdmissible increment := by
    intro increment member
    have listed : increment ∈ certificate.edges.map WindingEdge.increment := member
    obtain ⟨edge, _, image⟩ := List.mem_map.mp listed
    rw [← image]
    exact edge.increment_admissible
  have summed : certificate.toWinding.total = certificate.toWinding.target :=
    certificate.total_eq_target_of_valid accepted
  have full : certificate.toWinding.claims := ⟨nonempty, admissible, summed⟩
  unfold RationalCertificate.WindingCertificate.valid
  exact decide_true_with _ full

/-- **Ponte, forma usada na contagem.** Um certificado intervalar válido
entrega o mesmo total que o certificado inteiro extraído, e esse total é o
alvo. -/
theorem IntervalWindingCertificate.toWinding_total_eq_target
    (certificate : IntervalWindingCertificate)
    (accepted : certificate.valid = true) :
    certificate.toWinding.valid = true ∧
    certificate.toWinding.total = certificate.total ∧
    certificate.toWinding.total = certificate.target :=
  ⟨certificate.toWinding_valid accepted,
   certificate.toWinding_total,
   certificate.toWinding_total ▸ certificate.total_eq_target_of_valid accepted⟩

/-! ## Regressões

Todas decididas pelo kernel com `decide`; nenhuma usa `native_decide`. -/

/-! **Aritmética racional no kernel.** `unseal` abaixo devolve ao elaborador a
permissão de desdobrar as quatro funções `@[irreducible]` de `Rat`. Não é
hipótese nenhuma: irredutibilidade é diretiva de elaboração e o kernel nunca a
respeitou; o termo continua sendo checado inteiro. Com isso, `+`, `-`, `*`, `/`
e `^` sobre `Rat` reduzem por `decide`, e o helper `nodeBox` abaixo pode
escrever `centro ± raio` como se escreve no papel.

As caixas das primeiras regressões continuam com os quatro cantos já avaliados,
porque foram escritas assim e não há motivo para mexer; `dipContour` usa
`nodeBox` e literais com `/` justamente para deixar registrado, verificado pelo
kernel, que a restrição anunciada numa versão anterior deste arquivo não
existe. -/

unseal Rat.mul Rat.add Rat.sub Rat.inv

/-- Caixa quadrada de raio `radius` em torno de `(centerRe, centerIm)`. -/
def nodeBox (centerRe centerIm radius : Rat) : ComplexBox where
  reLo := centerRe - radius
  reHi := centerRe + radius
  imLo := centerIm - radius
  imHi := centerIm + radius

/-- Menor caixa que contém as duas caixas dadas. Só usa comparações, logo
reduz no kernel. -/
def hull (left right : ComplexBox) : ComplexBox where
  reLo := if left.reLo ≤ right.reLo then left.reLo else right.reLo
  reHi := if left.reHi ≤ right.reHi then right.reHi else left.reHi
  imLo := if left.imLo ≤ right.imLo then left.imLo else right.imLo
  imHi := if left.imHi ≤ right.imHi then right.imHi else left.imHi

theorem ComplexBox.contains_hull_left (left right : ComplexBox) :
    (hull left right).contains left = true := by
  unfold ComplexBox.contains hull
  simp only [Bool.and_eq_true, decide_eq_true_eq]
  refine ⟨⟨⟨?_, ?_⟩, ?_⟩, ?_⟩ <;> split <;>
    first
      | exact Rat.le_refl
      | assumption
      | (rename_i notLess; exact Rat.le_of_lt (Rat.not_le.mp notLess))

theorem ComplexBox.contains_hull_right (left right : ComplexBox) :
    (hull left right).contains right = true := by
  unfold ComplexBox.contains hull
  simp only [Bool.and_eq_true, decide_eq_true_eq]
  refine ⟨⟨⟨?_, ?_⟩, ?_⟩, ?_⟩ <;> split <;>
    first
      | exact Rat.le_refl
      | assumption
      | (rename_i notLess; exact Rat.le_of_lt (Rat.not_le.mp notLess))

/-! ### Geometria do contorno no plano `z`

Os quatro vértices abaixo formam o retângulo `[1, 2] × [-1/4, 1/4]`, que é a
forma do contorno `Γ_A` de `docs/SECTOR_B.md` §8. Percorridos na ordem
`SW → SE → NE → NW → SW` eles dão o sentido anti-horário.

As coordenadas usam `/`, isto é, divisão de `Rat`, recomputada pelo kernel. -/

def vertexSW : ComplexPoint := { re := 1, im := -1 / 4 }
def vertexSE : ComplexPoint := { re := 2, im := -1 / 4 }
def vertexNE : ComplexPoint := { re := 2, im := 1 / 4 }
def vertexNW : ComplexPoint := { re := 1, im := 1 / 4 }

/-- Aresta cujo `edgeBox` é o hull das caixas dos dois nós. Só é legítimo
quando o produtor já sabe que a imagem do segmento cabe nesse hull. -/
def straightEdge (startPoint endPoint : ComplexPoint)
    (startBox endBox : ComplexBox) : WindingEdge where
  startPoint := startPoint
  endPoint := endPoint
  startBox := startBox
  endBox := endBox
  edgeBox := hull startBox endBox

/-- Quatro caixas de valores de `f`, um quadrado de lado `20` centrado na
origem, cada nó já inflado por `1`: `(±10, ±10)` vira
`[±10-1, ±10+1] × [±10-1, ±10+1]`. -/
def cornerPP : ComplexBox := { reLo := 9, reHi := 11, imLo := 9, imHi := 11 }
def cornerMP : ComplexBox := { reLo := -11, reHi := -9, imLo := 9, imHi := 11 }
def cornerMM : ComplexBox := { reLo := -11, reHi := -9, imLo := -11, imHi := -9 }
def cornerPM : ComplexBox := { reLo := 9, reHi := 11, imLo := -11, imHi := -9 }

/-- Contorno quadrado em torno da origem, sentido anti-horário: winding `+1`.
O único cruzamento do eixo real negativo é o lado esquerdo, percorrido de
`Im > 0` para `Im < 0`. -/
def squareCounterClockwise : IntervalWindingCertificate where
  edges :=
    [ straightEdge vertexSW vertexSE cornerPP cornerMP,
      straightEdge vertexSE vertexNE cornerMP cornerMM,
      straightEdge vertexNE vertexNW cornerMM cornerPM,
      straightEdge vertexNW vertexSW cornerPM cornerPP ]
  target := 1

example : squareCounterClockwise.valid = true := by decide
example : squareCounterClockwise.total = 1 := by decide
example : squareCounterClockwise.increments = [0, 1, 0, 0] := by decide
example : squareCounterClockwise.toWinding.valid = true := by decide

/-- O mesmo quadrado percorrido no sentido horário: winding `-1`. O contorno no
plano `z` também é percorrido ao contrário. -/
def squareClockwise : IntervalWindingCertificate where
  edges :=
    [ straightEdge vertexSW vertexNW cornerPP cornerPM,
      straightEdge vertexNW vertexNE cornerPM cornerMM,
      straightEdge vertexNE vertexSE cornerMM cornerMP,
      straightEdge vertexSE vertexSW cornerMP cornerPP ]
  target := -1

example : squareClockwise.valid = true := by decide
example : squareClockwise.total = -1 := by decide
example : squareClockwise.increments = [0, 0, -1, 0] := by decide
example : squareClockwise.toWinding.valid = true := by decide

/-- Contorno que **não** envolve a origem: as caixas de valores ficam no
quadrado `[-40,-20] × [-10,10]`, inteiramente no semiplano esquerdo. Ele cruza
o eixo real negativo duas vezes, com sinais opostos, e o total é `0`. Regressão
mais forte que um contorno no semiplano direito, porque exercita os dois sinais
e o cancelamento. -/
def leftCornerMM : ComplexBox := { reLo := -41, reHi := -39, imLo := -11, imHi := -9 }
def leftCornerPM : ComplexBox := { reLo := -21, reHi := -19, imLo := -11, imHi := -9 }
def leftCornerPP : ComplexBox := { reLo := -21, reHi := -19, imLo := 9, imHi := 11 }
def leftCornerMP : ComplexBox := { reLo := -41, reHi := -39, imLo := 9, imHi := 11 }

def squareAvoidingOrigin : IntervalWindingCertificate where
  edges :=
    [ straightEdge vertexSW vertexSE leftCornerMM leftCornerPM,
      straightEdge vertexSE vertexNE leftCornerPM leftCornerPP,
      straightEdge vertexNE vertexNW leftCornerPP leftCornerMP,
      straightEdge vertexNW vertexSW leftCornerMP leftCornerMM ]
  target := 0

example : squareAvoidingOrigin.valid = true := by decide
example : squareAvoidingOrigin.total = 0 := by decide
example : squareAvoidingOrigin.increments = [0, -1, 0, 1] := by decide
example : squareAvoidingOrigin.toWinding.valid = true := by decide

/-- Regressão com racionais não inteiros escritos como `mkRat`, que é a forma
emitida pelo produtor externo. As caixas são as do quadrado anti-horário
divididas por `10`. -/
def scaledCornerPP : ComplexBox :=
  { reLo := mkRat 9 10, reHi := mkRat 11 10, imLo := mkRat 9 10, imHi := mkRat 11 10 }
def scaledCornerMP : ComplexBox :=
  { reLo := mkRat (-11) 10, reHi := mkRat (-9) 10, imLo := mkRat 9 10, imHi := mkRat 11 10 }
def scaledCornerMM : ComplexBox :=
  { reLo := mkRat (-11) 10, reHi := mkRat (-9) 10, imLo := mkRat (-11) 10, imHi := mkRat (-9) 10 }
def scaledCornerPM : ComplexBox :=
  { reLo := mkRat 9 10, reHi := mkRat 11 10, imLo := mkRat (-11) 10, imHi := mkRat (-9) 10 }

def scaledSquare : IntervalWindingCertificate where
  edges :=
    [ straightEdge vertexSW vertexSE scaledCornerPP scaledCornerMP,
      straightEdge vertexSE vertexNE scaledCornerMP scaledCornerMM,
      straightEdge vertexNE vertexNW scaledCornerMM scaledCornerPM,
      straightEdge vertexNW vertexSW scaledCornerPM scaledCornerPP ]
  target := 1

example : scaledSquare.valid = true := by decide
example : scaledSquare.total = 1 := by decide

/-! ### O quarto caso do semiplano esquerdo, exercitado

As caixas de valores ficam em torno do retângulo `[-30,-10] × [-10,10]`, inteiro
no semiplano esquerdo e com a origem **fora**. A primeira aresta liga os valores
`(-10,10)` e `(-30,10)` mas **mergulha** até `Im = -13/2` no meio do caminho: o
`edgeBox` não é estritamente `Im > 0` nem estritamente `Im < 0`, então os ramos
1 e 2 não se aplicam; e os dois nós têm `Im > 0`, então os casos de sinal oposto
também não. Antes da correção esta aresta devolvia `none` e forçava refinamento
do contorno; agora ela classifica como `0`, que é o valor demonstrado por
`edgeCrossing_eq_zero_of_sameSignIm`.

Este é o caso que aparece perto do cruzamento forçado do eixo real
(`docs/WINDING_DATA.md` §2.2), onde as arestas ficam inteiras em `Re < 0`.

As caixas usam `nodeBox` (`centro ± raio`) e o literal `-13/2`, isto é,
subtração, adição e divisão de `Rat` — tudo recomputado por `decide`. É a
verificação da correção da nota de cabeçalho. -/

def dipNodeA : ComplexBox := nodeBox (-10) 10 1
def dipNodeB : ComplexBox := nodeBox (-30) 10 1
def dipNodeC : ComplexBox := nodeBox (-30) (-10) 1
def dipNodeD : ComplexBox := nodeBox (-10) (-10) 1

/-- A aresta que mergulha: os dois nós com `Im > 0`, a caixa cruzando `Im = 0`,
tudo em `Re < 0`. -/
def dipEdge : WindingEdge where
  startPoint := vertexSW
  endPoint := vertexSE
  startBox := dipNodeA
  endBox := dipNodeB
  edgeBox := { reLo := -31, reHi := -9, imLo := -13 / 2, imHi := 11 }

-- A aritmética de `Rat` reduz no kernel: `nodeBox` calcula `centro ± raio`.
example : dipNodeA = { reLo := -11, reHi := -9, imLo := 9, imHi := 11 } := by decide

-- Nenhum dos dois primeiros ramos se aplica, e o terceiro não é de sinal oposto.
example : dipEdge.edgeBox.strictlyPosRe = false := by decide
example : dipEdge.edgeBox.strictlyPosIm = false := by decide
example : dipEdge.edgeBox.strictlyNegIm = false := by decide
example : dipEdge.edgeBox.strictlyNegRe = true := by decide
example : dipEdge.shapeOk = true := by decide
example : dipEdge.nondegenerate = true := by decide
example : dipEdge.edgeCrossing = some 0 := by decide

def dipContour : IntervalWindingCertificate where
  edges :=
    [ dipEdge,
      straightEdge vertexSE vertexNE dipNodeB dipNodeC,
      straightEdge vertexNE vertexNW dipNodeC dipNodeD,
      straightEdge vertexNW vertexSW dipNodeD dipNodeA ]
  target := 0

example : dipContour.increments = [0, 1, 0, -1] := by decide
example : dipContour.valid = true := by decide
example : dipContour.total = 0 := by decide
example : dipContour.toWinding.valid = true := by decide

/-- **Rejeição 1: aresta ambígua.** O mesmo quadrado anti-horário, mas com o
`edgeBox` do lado esquerdo grosseiro demais (ele contém a origem). Nenhum dos
três ramos se aplica, `edgeCrossing` devolve `none` e o certificado é
recusado: o contorno precisa ser refinado. -/
def ambiguousEdge : WindingEdge where
  startPoint := vertexSE
  endPoint := vertexNE
  startBox := cornerMP
  endBox := cornerMM
  edgeBox := { reLo := -11, reHi := 1, imLo := -11, imHi := 11 }

def ambiguousContour : IntervalWindingCertificate where
  edges :=
    [ straightEdge vertexSW vertexSE cornerPP cornerMP,
      ambiguousEdge,
      straightEdge vertexNE vertexNW cornerMM cornerPM,
      straightEdge vertexNW vertexSW cornerPM cornerPP ]
  target := 1

example : ambiguousEdge.shapeOk = true := by decide
example : ambiguousEdge.edgeCrossing = none := by decide
example : ambiguousContour.valid = false := by decide

/-- **Rejeição 2: encadeamento quebrado.** A aresta que ia das caixas
`(-10,10)` a `(-10,-10)` foi removida, de modo que nem os pontos nem as caixas
de arestas consecutivas casam. O contorno ainda "fecha" e a soma bate com o
alvo `0`; mesmo assim o certificado é recusado, porque a lista não descreve um
ciclo e a soma telescópica não significaria nada. -/
def brokenChainContour : IntervalWindingCertificate where
  edges :=
    [ straightEdge vertexSW vertexSE cornerPP cornerMP,
      straightEdge vertexNE vertexNW cornerMM cornerPM,
      straightEdge vertexNW vertexSW cornerPM cornerPP ]
  target := 0

example : brokenChainContour.total = 0 := by decide
example : (∀ edge ∈ brokenChainContour.edges, edge.classified = true) := by decide
example : closed brokenChainContour.edges = true := by decide
example : chained brokenChainContour.edges = false := by decide
example : brokenChainContour.valid = false := by decide

/-- **Rejeição 3: contorno aberto.** Tirando a última aresta do quadrado
anti-horário, o encadeamento continua correto mas o contorno não fecha. -/
def openContour : IntervalWindingCertificate where
  edges :=
    [ straightEdge vertexSW vertexSE cornerPP cornerMP,
      straightEdge vertexSE vertexNE cornerMP cornerMM,
      straightEdge vertexNE vertexNW cornerMM cornerPM ]
  target := 1

example : chained openContour.edges = true := by decide
example : closed openContour.edges = false := by decide
example : openContour.valid = false := by decide

/-- **Rejeição 4: soma diferente do alvo.** -/
def wrongTargetContour : IntervalWindingCertificate where
  edges := squareCounterClockwise.edges
  target := 2

example : wrongTargetContour.valid = false := by decide

/-! ### As duas rejeições que só a geometria pega

Estas duas são o motivo de os pontos existirem: antes de o certificado carregar
`startPoint`/`endPoint`, os dois contornos abaixo eram **aceitos**. -/

/-- **Rejeição 5: as caixas fecham, os pontos não.**

Exatamente o quadrado anti-horário, com um único campo diferente: a última
aresta termina em `vertexSE` em vez de voltar a `vertexSW`. Tudo o que o
certificado antigo sabia olhar continua perfeito — as caixas encadeiam, a caixa
final coincide com a inicial, as quatro arestas classificam, a soma é `1` e bate
com o alvo. No plano `z`, porém, a curva não fecha: ela termina num vértice
diferente daquele em que começou, e não há contorno nenhum. -/
def boxChainOnlyContour : IntervalWindingCertificate where
  edges :=
    [ straightEdge vertexSW vertexSE cornerPP cornerMP,
      straightEdge vertexSE vertexNE cornerMP cornerMM,
      straightEdge vertexNE vertexNW cornerMM cornerPM,
      straightEdge vertexNW vertexSE cornerPM cornerPP ]
  target := 1

-- Tudo o que o certificado antigo verificava continua valendo:
example : chained boxChainOnlyContour.edges = true := by decide
example : (∀ edge ∈ boxChainOnlyContour.edges, edge.shapeOk = true) := by decide
example : (∀ edge ∈ boxChainOnlyContour.edges, edge.classified = true) := by decide
example : (∀ edge ∈ boxChainOnlyContour.edges, edge.nondegenerate = true) := by decide
example : boxChainOnlyContour.total = 1 := by decide
-- As caixas fecham; os pontos não. É só isto que recusa o certificado.
example :
    ((boxChainOnlyContour.edges.getLast?).map WindingEdge.endBox
      = (boxChainOnlyContour.edges.head?).map WindingEdge.startBox) := by decide
example :
    ¬ ((boxChainOnlyContour.edges.getLast?).map WindingEdge.endPoint
      = (boxChainOnlyContour.edges.head?).map WindingEdge.startPoint) := by decide
example : closed boxChainOnlyContour.edges = false := by decide
example : boxChainOnlyContour.valid = false := by decide

/-- **Rejeição 6: aresta degenerada.**

Um ciclo de cinco arestas em que a segunda vai de `vertexSE` para `vertexSE`.
Ela encadeia com as vizinhas, classifica (a caixa está toda em `Im > 0`, ramo
2) e contribui `0` para a soma, de modo que encadeamento, fechamento,
classificação e total continuam todos corretos. Uma aresta de comprimento zero
não é uma aresta do contorno: ela ocupa uma posição na cadeia sem descrever
segmento nenhum, e é assim que um produtor pode fabricar encadeamento onde não
há geometria. -/
def degenerateContour : IntervalWindingCertificate where
  edges :=
    [ straightEdge vertexSW vertexSE cornerPP cornerMP,
      straightEdge vertexSE vertexSE cornerMP cornerMP,
      straightEdge vertexSE vertexNE cornerMP cornerMM,
      straightEdge vertexNE vertexNW cornerMM cornerPM,
      straightEdge vertexNW vertexSW cornerPM cornerPP ]
  target := 1

-- Tudo o mais passa:
example : chained degenerateContour.edges = true := by decide
example : closed degenerateContour.edges = true := by decide
example : (∀ edge ∈ degenerateContour.edges, edge.shapeOk = true) := by decide
example : (∀ edge ∈ degenerateContour.edges, edge.classified = true) := by decide
example : degenerateContour.increments = [0, 0, 1, 0, 0] := by decide
example : degenerateContour.total = 1 := by decide
-- E ainda assim é recusado, pela aresta de comprimento zero:
example : (straightEdge vertexSE vertexSE cornerMP cornerMP).nondegenerate = false := by
  decide
example : degenerateContour.valid = false := by decide

#print axioms ComplexBox.mem_of_contains
#print axioms ComplexBox.contains_hull_left
#print axioms WindingEdge.increment_admissible
#print axioms WindingEdge.edgeCrossing_eq_zero_of_sameSignIm
#print axioms WindingEdge.increment_eq_zero_of_sameSignIm
#print axioms WindingEdge.zero_not_mem_edgeBox
#print axioms IntervalWindingCertificate.valid_sound
#print axioms IntervalWindingCertificate.total_eq_target_of_valid
#print axioms IntervalWindingCertificate.geometryClosed
#print axioms IntervalWindingCertificate.nondegenerate_edges
#print axioms IntervalWindingCertificate.zero_not_mem_edgeBoxes
#print axioms IntervalWindingCertificate.zero_not_mem_nodeBoxes
#print axioms IntervalWindingCertificate.toWinding_valid
#print axioms IntervalWindingCertificate.toWinding_total_eq_target

end Choptuik.WindingInterval
