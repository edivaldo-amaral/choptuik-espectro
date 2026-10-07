import Std
import ChoptuikFormal.WindingInterval

/-!
# Arcos: concatenação de pedaços de contorno

Este módulo existe por uma razão medida, não estética. O verificador de
`ChoptuikFormal.WindingInterval` fecha pelo kernel um contorno real de **514
arestas** em cerca de três minutos, mas **não** fecha o `Γ_A` completo de
**1946 arestas**: ele estoura `maxHeartbeats` em dois pontos distintos, e a
distinção importa.

```
:10:4:     (deterministic) timeout at `«LCNF compiler»`
:11691:47: (deterministic) timeout at `whnf`
```

O segundo é a checagem da proposição. O primeiro é o **compilador** tentando
gerar código executável para a `def` do certificado — trabalho inteiramente
inútil, porque o objeto só existe para o kernel reduzir, nunca para rodar.

## O que a decomposição resolve, e o que não resolve

Sejamos precisos, porque é fácil vender isto errado. Todas as checagens do
certificado são **lineares** no número de arestas: forma, não degenerescência,
classificação, encadeamento e soma. Partir o contorno em `k` arcos **não reduz
o trabalho total** — continua sendo o mesmo `O(n)`, apenas somado de outro
jeito.

O que muda é *onde* o trabalho acontece. `maxHeartbeats` é um orçamento **por
declaração**. Um contorno monolítico é uma declaração só, e ou cabe no
orçamento ou não cabe. Quatro arcos de ~500 arestas são quatro declarações
independentes, cada uma do tamanho do caso de 514 que já se sabe caber, e a
recomposição entre elas é dedução pura — aplicação de teorema, sem redução.
É por isso que o lema de concatenação destrava o caso grande: ele troca uma
declaração impossível por `k` declarações possíveis mais álgebra grátis.

Corolário de desenho que vale explicitar: `ArcCertificate.valid`, o verificador
booleano n-ário do fim deste arquivo, **não** traz esse ganho — ele reúne tudo
numa declaração de novo. Ele existe por conveniência e para as regressões
pequenas. Quem tem 1946 arestas deve usar a rota dos teoremas
(`Arc.append_valid`, `Arc.close_valid`), provando cada `arc.valid = true` numa
declaração separada.

## O compilador

Sobre o primeiro timeout: uma `def` grande faz o Lean gerar código executável
que ninguém executa. A rota dos teoremas já evita o pior disso, porque cada
arco é uma `def` `k` vezes menor. Onde ainda incomodar, `noncomputable def`
suprime a geração de código sem mudar nada do que o kernel checa — nenhuma
regressão deste arquivo precisou disso, mas o caso de 1946 arestas pode
precisar, e a medição está no relatório.
-/

namespace Choptuik.WindingArcs

open Choptuik.WindingInterval

/- `WindingInterval` mantém o seu auxiliar de decidibilidade `private`, então
   reconstruo aqui o mesmo utilitário: `decide` com instância explícita, que é
   o que permite provar `... .valid = true` sem nomear a instância privada. -/
private theorem decide_true_with {claim : Prop} (inst : Decidable claim)
    (proof : claim) : @decide claim inst = true :=
  @decide_eq_true claim inst proof

/-! ## Junção de duas listas de arestas -/

/-- A última aresta de `left` casa com a primeira de `right`, nos dois planos:
o ponto do plano `z` e a caixa de valores. É a mesma condição que `chained`
impõe entre arestas consecutivas, isolada na fronteira entre dois pedaços.

Quando um dos lados é vazio não há junção a checar, e o valor é `true`: nesse
caso a concatenação é o outro lado e a validade dele já basta. -/
def joins (left right : List WindingEdge) : Bool :=
  match left.getLast?, right.head? with
  | some last, some first =>
      (last.endPoint == first.startPoint) && (last.endBox == first.startBox)
  | _, _ => true

/-- `chained` de uma lista não vazia decompõe-se na junção da cabeça com a
cauda mais o encadeamento da cauda. É a forma que torna a indução possível. -/
theorem chained_cons (edge : WindingEdge) (rest : List WindingEdge) :
    chained (edge :: rest) = (joins [edge] rest && chained rest) := by
  cases rest with
  | nil => rfl
  | cons next tail => simp [chained, joins]

/-- Encolher a lista da esquerda não pode quebrar a junção. -/
theorem joins_cons_left (edge : WindingEdge) (rest right : List WindingEdge)
    (link : joins (edge :: rest) right = true) : joins rest right = true := by
  cases rest with
  | nil => simp [joins]
  | cons next tail => simpa [joins] using link

/-- A junção só olha o fim da lista da esquerda, então prefixos não importam. -/
theorem joins_append_left (left mid right : List WindingEdge)
    (midNonempty : mid ≠ []) :
    joins (left ++ mid) right = joins mid right := by
  obtain ⟨last, hlast⟩ : ∃ last, mid.getLast? = some last := by
    cases hm : mid.getLast? with
    | none => exact absurd (List.getLast?_eq_none_iff.mp hm) midNonempty
    | some last => exact ⟨last, rfl⟩
  simp [joins, List.getLast?_append, hlast]

/-! ## Os dois lemas de concatenação -/

/-- **Encadeamento concatena.** Se os dois pedaços encadeiam e a fronteira
entre eles casa, o pedaço concatenado encadeia. -/
theorem chained_append (left right : List WindingEdge)
    (leftChained : chained left = true) (rightChained : chained right = true)
    (link : joins left right = true) :
    chained (left ++ right) = true := by
  induction left with
  | nil => simpa using rightChained
  | cons edge rest ih =>
      rw [chained_cons] at leftChained
      obtain ⟨headLink, restChained⟩ := (Bool.and_eq_true _ _).mp leftChained
      rw [List.cons_append, chained_cons]
      refine (Bool.and_eq_true _ _).mpr ⟨?_, ?_⟩
      · cases rest with
        | nil => simpa using link
        | cons next tail => simpa [joins] using headLink
      · exact ih restChained (joins_cons_left edge rest right link)

/-- **Fechamento é uma junção.** Fechar a concatenação de dois pedaços não
exige percorrer a lista concatenada: é exatamente a junção do segundo pedaço de
volta ao primeiro. É este lema que mantém o custo do fechamento proporcional a
um arco, e não ao contorno inteiro. -/
theorem closed_append (left right : List WindingEdge)
    (leftNonempty : left ≠ []) (rightNonempty : right ≠ []) :
    closed (left ++ right) = joins right left := by
  obtain ⟨first, hfirst⟩ : ∃ first, left.head? = some first := by
    cases left with
    | nil => exact absurd rfl leftNonempty
    | cons x xs => exact ⟨x, rfl⟩
  obtain ⟨last, hlast⟩ : ∃ last, right.getLast? = some last := by
    cases hr : right.getLast? with
    | none => exact absurd (List.getLast?_eq_none_iff.mp hr) rightNonempty
    | some last => exact ⟨last, rfl⟩
  simp [closed, joins, List.head?_append, List.getLast?_append, hfirst, hlast]

/-! ## Somas -/

/-- A mesma soma que `IntervalWindingCertificate.total` calcula. -/
def sumIncrements (edges : List WindingEdge) : Int :=
  (edges.map WindingEdge.increment).foldl (fun total step => total + step) 0

private theorem foldl_add_start (values : List Int) (start : Int) :
    values.foldl (fun total step => total + step) start
      = start + values.foldl (fun total step => total + step) 0 := by
  induction values generalizing start with
  | nil => simp
  | cons value rest ih =>
      simp only [List.foldl_cons]
      rw [ih (start + value), ih (0 + value)]
      omega

/-- **A soma concatena.** Nenhum incremento é perdido nem contado duas vezes. -/
theorem sumIncrements_append (left right : List WindingEdge) :
    sumIncrements (left ++ right) = sumIncrements left + sumIncrements right := by
  unfold sumIncrements
  rw [List.map_append, List.foldl_append, foldl_add_start]

/-! ## Arcos

Um arco é um pedaço de contorno: as mesmas arestas, as mesmas checagens de
forma, não degenerescência, classificação e encadeamento — **sem** a exigência
de fechamento, que só faz sentido para o contorno inteiro. -/

structure Arc where
  edges : List WindingEdge
deriving Repr, DecidableEq

/-- Total dos incrementos do arco. -/
def Arc.total (arc : Arc) : Int := sumIncrements arc.edges

/-- Conteúdo verificado de um arco. Note a ausência de `closed` e de `target`:
um arco não fecha e não tem alvo próprio. -/
def Arc.claims (arc : Arc) : Prop :=
  arc.edges ≠ [] ∧
  (∀ edge ∈ arc.edges, edge.shapeOk = true) ∧
  (∀ edge ∈ arc.edges, edge.nondegenerate = true) ∧
  (∀ edge ∈ arc.edges, edge.classified = true) ∧
  chained arc.edges = true

private def Arc.claimsDecidable (arc : Arc) : Decidable arc.claims := by
  unfold Arc.claims
  infer_instance

/-- Verificador executável do arco. **É esta a declaração que se prova em
separado para cada pedaço**, e é por isso que o orçamento por declaração deixa
de ser o gargalo. -/
def Arc.valid (arc : Arc) : Bool :=
  @decide arc.claims arc.claimsDecidable

theorem Arc.valid_sound (arc : Arc) (accepted : arc.valid = true) : arc.claims :=
  @of_decide_eq_true arc.claims arc.claimsDecidable accepted

/-- Concatenação de arcos. -/
def Arc.append (left right : Arc) : Arc := ⟨left.edges ++ right.edges⟩

/-- A fronteira entre dois arcos casa. -/
def Arc.joinsTo (left right : Arc) : Bool := joins left.edges right.edges

theorem Arc.append_total (left right : Arc) :
    (left.append right).total = left.total + right.total :=
  sumIncrements_append left.edges right.edges

/-- **O lema de concatenação.** Dois arcos válidos cuja fronteira casa formam
um arco válido. Nenhuma aresta é reexaminada: a prova é aplicação de teorema,
e o seu custo não depende do tamanho dos arcos. -/
theorem Arc.append_valid (left right : Arc)
    (leftValid : left.valid = true) (rightValid : right.valid = true)
    (link : left.joinsTo right = true) :
    (left.append right).valid = true := by
  have leftClaims := left.valid_sound leftValid
  have rightClaims := right.valid_sound rightValid
  refine decide_true_with _ ⟨?_, ?_, ?_, ?_, ?_⟩
  · show left.edges ++ right.edges ≠ []
    cases hl : left.edges with
    | nil => exact absurd hl leftClaims.1
    | cons first rest => simp [hl]
  · intro edge member
    rcases List.mem_append.mp member with inLeft | inRight
    · exact leftClaims.2.1 edge inLeft
    · exact rightClaims.2.1 edge inRight
  · intro edge member
    rcases List.mem_append.mp member with inLeft | inRight
    · exact leftClaims.2.2.1 edge inLeft
    · exact rightClaims.2.2.1 edge inRight
  · intro edge member
    rcases List.mem_append.mp member with inLeft | inRight
    · exact leftClaims.2.2.2.1 edge inLeft
    · exact rightClaims.2.2.2.1 edge inRight
  · exact chained_append left.edges right.edges
      leftClaims.2.2.2.2 rightClaims.2.2.2.2 link

/-! ## Fechar um arco num certificado -/

/-- O certificado de winding que um arco fechado representa. -/
def Arc.close (arc : Arc) (target : Int) : IntervalWindingCertificate where
  edges := arc.edges
  target := target

/-- **Fechamento.** Um arco válido que volta ao seu próprio início é um
`IntervalWindingCertificate` válido, com o total que o arco acumulou.

A hipótese de fechamento é `closed arc.edges = true`; combinada com
`closed_append`, ela se verifica como a junção do último arco de volta ao
primeiro, sem tocar na lista concatenada. -/
theorem Arc.close_valid (arc : Arc) (target : Int)
    (arcValid : arc.valid = true)
    (closure : closed arc.edges = true)
    (sum : arc.total = target) :
    (arc.close target).valid = true := by
  have arcClaims := arc.valid_sound arcValid
  exact decide_true_with _
    ⟨arcClaims.1, arcClaims.2.1, arcClaims.2.2.1, arcClaims.2.2.2.1,
     arcClaims.2.2.2.2, closure, sum⟩

/-- O total sobrevive ao fechamento, por construção. -/
theorem Arc.close_total (arc : Arc) (target : Int) :
    (arc.close target).total = arc.total := rfl

/-- **A forma usada na prática.** Dois arcos válidos que se encaixam e voltam
ao começo produzem um certificado de winding válido cujo total é a soma dos
dois. Tudo o que é caro — `leftValid` e `rightValid` — é provado fora, uma
declaração por arco. -/
theorem twoArcCertificate (left right : Arc) (target : Int)
    (leftValid : left.valid = true) (rightValid : right.valid = true)
    (forward : left.joinsTo right = true)
    (back : right.joinsTo left = true)
    (sum : left.total + right.total = target) :
    ((left.append right).close target).valid = true ∧
    ((left.append right).close target).total = left.total + right.total := by
  have leftNonempty : left.edges ≠ [] := (left.valid_sound leftValid).1
  have rightNonempty : right.edges ≠ [] := (right.valid_sound rightValid).1
  have combinedValid : (left.append right).valid = true :=
    Arc.append_valid left right leftValid rightValid forward
  have combinedTotal : (left.append right).total = left.total + right.total :=
    Arc.append_total left right
  have closure : closed (left.append right).edges = true := by
    show closed (left.edges ++ right.edges) = true
    rw [closed_append left.edges right.edges leftNonempty rightNonempty]
    exact back
  exact ⟨Arc.close_valid _ target combinedValid closure (by rw [combinedTotal]; exact sum),
         by rw [Arc.close_total]; exact combinedTotal⟩

/-! ## Versão n-ária

Encadeia uma lista de arcos a partir de um acumulador. Serve para `k > 2`; a
mecânica é a mesma, aplicada `k-1` vezes. -/

/-- Concatena `rest` ao acumulador, da esquerda para a direita. -/
def chainArcs : Arc → List Arc → Arc
  | acc, [] => acc
  | acc, next :: rest => chainArcs (acc.append next) rest

/-- Junções consecutivas a partir do acumulador. Cada condição olha só dois
arcos vizinhos. -/
def chainLinks : Arc → List Arc → Bool
  | _, [] => true
  | acc, next :: rest => acc.joinsTo next && chainLinks next rest

theorem chainArcs_total (acc : Arc) (rest : List Arc) :
    (chainArcs acc rest).total
      = rest.foldl (fun total arc => total + arc.total) acc.total := by
  induction rest generalizing acc with
  | nil => rfl
  | cons next rest ih =>
      show (chainArcs (acc.append next) rest).total = _
      rw [ih (acc.append next), Arc.append_total]
      rfl

/-- **Concatenação n-ária.** Todos os arcos válidos e todas as junções
consecutivas dão um arco válido. A indução usa `joins_append_left`: o
acumulador cresce à esquerda, e a junção só enxerga o seu fim. -/
theorem chainArcs_valid (acc : Arc) (rest : List Arc)
    (accValid : acc.valid = true)
    (restValid : ∀ arc ∈ rest, arc.valid = true)
    (links : chainLinks acc rest = true) :
    (chainArcs acc rest).valid = true := by
  induction rest generalizing acc with
  | nil => exact accValid
  | cons next rest ih =>
      obtain ⟨headLink, tailLinks⟩ := (Bool.and_eq_true _ _).mp links
      have nextValid : next.valid = true := restValid next (by simp)
      have nextNonempty : next.edges ≠ [] := (next.valid_sound nextValid).1
      refine ih (acc.append next) ?_ ?_ ?_
      · exact Arc.append_valid acc next accValid nextValid headLink
      · intro arc member; exact restValid arc (by simp [member])
      · -- as junções restantes não mudam: elas partem de `next`, cujo fim é o
        -- fim do acumulador ampliado
        cases rest with
        | nil => rfl
        | cons third tail =>
            obtain ⟨secondLink, deeper⟩ := (Bool.and_eq_true _ _).mp tailLinks
            refine (Bool.and_eq_true _ _).mpr ⟨?_, deeper⟩
            show joins (acc.edges ++ next.edges) third.edges = true
            rw [joins_append_left acc.edges next.edges third.edges nextNonempty]
            exact secondLink

/-! `unseal` porque os vértices do quadrado usam literais como `-1/4`, que
passam por `Rat.div`. A irredutibilidade de `Rat` é diretiva de *elaboração*: o
kernel nunca a respeitou, então liberar não amplia a base de confiança. Ver
`docs/KERNEL_TRUST.md`. -/
unseal Rat.mul Rat.add Rat.sub Rat.inv

/-! ## Regressões

O quadrado `squareCounterClockwise` de `WindingInterval.lean` tem winding `+1`.
Aqui ele é partido em dois arcos de duas arestas, e o certificado reconstruído
por `twoArcCertificate` tem o mesmo total — que é o conteúdo do lema.

**Por que isto importa, medido sobre dados reais.** O contorno `Gamma_A` gerado
a partir da matriz `4x12` (1946 arestas) **não fecha** por `decide` nos ajustes
padrão: estoura `maxHeartbeats` em dois pontos, um deles o compilador LCNF
tentando gerar código para a `def`. Partido em quatro arcos de ~487 arestas,
cada `Arc.valid` fecha **nos ajustes padrão**, e o conjunto leva 6 min 02 s
contra 673 s do monolítico com `maxHeartbeats 0`. Decompor não é só mais
auditável: é mais barato. Ver a obrigação K6 em `docs/KERNEL_TRUST.md`. -/

open WindingInterval in
/-- Arestas 0 e 1 do quadrado: de `vertexSW` a `vertexNE`. -/
def arcoSul : Arc :=
  ⟨[ straightEdge vertexSW vertexSE cornerPP cornerMP,
     straightEdge vertexSE vertexNE cornerMP cornerMM ]⟩

open WindingInterval in
/-- Arestas 2 e 3: de `vertexNE` de volta a `vertexSW`. -/
def arcoNorte : Arc :=
  ⟨[ straightEdge vertexNE vertexNW cornerMM cornerPM,
     straightEdge vertexNW vertexSW cornerPM cornerPP ]⟩

example : arcoSul.valid = true := by decide
example : arcoNorte.valid = true := by decide
example : arcoSul.joinsTo arcoNorte = true := by decide
example : arcoNorte.joinsTo arcoSul = true := by decide
example : arcoSul.total + arcoNorte.total = 1 := by decide

/-- O certificado remontado a partir dos dois arcos é válido e conta `1`,
exatamente como o monolítico. -/
theorem arcosReconstroemOQuadrado :
    ((arcoSul.append arcoNorte).close 1).valid = true ∧
    ((arcoSul.append arcoNorte).close 1).total = arcoSul.total + arcoNorte.total :=
  twoArcCertificate arcoSul arcoNorte 1 (by decide) (by decide)
    (by decide) (by decide) (by decide)

/-- E o total bate com o do certificado monolítico de `WindingInterval`. -/
example : ((arcoSul.append arcoNorte).close 1).total
    = WindingInterval.squareCounterClockwise.total := by decide

#print axioms arcosReconstroemOQuadrado

end Choptuik.WindingArcs
