import Std

/-!
# Composição exata do bound de cauda `‖C - C_N‖ ≤ ε_N` (obrigação C1)

Este módulo **não** demonstra nenhuma estimativa de EDP. Ele faz exatamente uma
coisa, e faz com aritmética racional exata verificada pelo kernel: dada a tripla
de constantes que as três hipóteses de EDP (batizadas `P1`, `P2`, `P3` em
`docs/TAIL_BOUND.md`) teriam de fornecer, ele **compõe** essas constantes na
estimativa de truncamento

```
ε_N(N_F, N_C)
  = M_J · A_0 · M_0 · ( c_F · q_F^(N_F+1) + c_C · q_C^(N_C+1) )
```

e decide, sem ponto flutuante, se `ε_N ≤ alvo`.

Significado de cada fator (derivação completa em `docs/TAIL_BOUND.md`):

* `M_0` (`resolventBound`, hipótese **P1**): `‖R_0‖ ≤ M_0` no espaço base,
  `R_0 = P(z_0)⁻¹`;
* `A_0` (`analyticBound`, hipótese **P2**): estimativa a priori analítica
  `‖u‖_{σ,ρ} ≤ A_0 ‖u‖` para `u` na imagem de `R_0`;
* `M_J` (`pencilBound`, hipótese **P3**): `‖J‖ ≤ M_J` na escala analítica;
* `q_F` (`decayFourier`): **majorante racional** de `e^{-σ}`, onde `σ > 0` é a
  largura da faixa de analiticidade em `τ_R`;
* `q_C` (`decayChebyshev`): **majorante racional** de `1/ρ`, onde `ρ > 1` é o
  parâmetro da elipse de Bernstein em `ξ`;
* `c_F`, `c_C` (`fourierTailFactor`, `chebyshevTailFactor`): constantes da
  realização de norma escolhida. Na realização de Wiener (ℓ¹ com pesos) valem
  `c_F = c_C = 1`; na realização de norma do supremo via Cauchy--Bernstein
  valem as somas geométricas `supFactorFourier`/`supFactorChebyshev` abaixo.

Fornecer `q_F` e `q_C` **racionais** é parte da obrigação `P2`: uma prova que
só diga "existe `σ > 0`" não alimenta este verificador. `e^{-σ}` e `1/ρ` são
irracionais em geral; a obrigação analítica tem de terminar com um majorante
racional explícito, e é esse majorante que entra aqui.

Nada neste arquivo é um teorema sobre o operador de Choptuik. O que o kernel
aceita é: *se* `P1`, `P2` e `P3` valerem com estas constantes, *então* a
desigualdade aritmética anunciada vale.
-/

namespace Choptuik.TailBound

/- No Lean 4.33 `Rat.mul`, `Rat.add`, `Rat.sub` e `Rat.inv` são `@[irreducible]`,
   de modo que a tática `decide` empaca ao avaliar aritmética racional — um dos
   dois bloqueios diagnosticados em `docs/KERNEL_TRUST.md` (o outro, `Array`,
   não afeta este módulo, que não usa `Array`). `unseal`
   apenas devolve a transparência normal a essas quatro definições: é uma
   anotação de redutibilidade, não acrescenta hipótese nenhuma e o kernel
   continua reconferindo cada redução. Portanto o bloqueio de `Rat` não é
   estrutural: é atributo de redutibilidade. Assim todas as regressões abaixo usam
   `decide` (kernel) e nenhuma usa `native_decide` (que confiaria no compilador
   fora do kernel, via `Lean.ofReduceBool`). -/
unseal Rat.mul Rat.add Rat.sub Rat.inv

/-- Constantes das três hipóteses de EDP mais os truncamentos.

`fourierTruncation = N_F` guarda os modos de Fourier `|n| ≤ N_F` em `τ_R`;
`chebyshevTruncation = N_C` guarda os graus de Chebyshev `m ≤ N_C` em `ξ`. -/
structure TailData where
  /-- `M_0` de **P1**. -/
  resolventBound : Rat
  /-- `A_0` de **P2**. -/
  analyticBound : Rat
  /-- `M_J` de **P3**. -/
  pencilBound : Rat
  /-- Majorante racional de `e^{-σ}`; parte de **P2**. -/
  decayFourier : Rat
  /-- Majorante racional de `1/ρ`; parte de **P2**. -/
  decayChebyshev : Rat
  /-- `c_F`: constante de soma da cauda de Fourier na norma escolhida. -/
  fourierTailFactor : Rat
  /-- `c_C`: constante de soma da cauda de Chebyshev na norma escolhida. -/
  chebyshevTailFactor : Rat
  /-- `N_F`. -/
  fourierTruncation : Nat
  /-- `N_C`. -/
  chebyshevTruncation : Nat
deriving Repr, DecidableEq

/-- Cauda de Fourier `c_F · q_F^(N_F+1)`.

O expoente é `N_F + 1`, e não `N_F`: o truncamento guarda `|n| ≤ N_F`, logo o
primeiro modo descartado é `|n| = N_F + 1`. -/
def TailData.fourierTail (data : TailData) : Rat :=
  data.fourierTailFactor * data.decayFourier ^ (data.fourierTruncation + 1)

/-- Cauda de Chebyshev `c_C · q_C^(N_C+1)`, pela mesma razão de indexação. -/
def TailData.chebyshevTail (data : TailData) : Rat :=
  data.chebyshevTailFactor * data.decayChebyshev ^ (data.chebyshevTruncation + 1)

/-- `ε_N = M_J · A_0 · M_0 · (cauda de Fourier + cauda de Chebyshev)`.

A soma das duas caudas vem de cobrir o complemento do retângulo de truncamento
pela união `{|n| > N_F} ∪ {m > N_C}`; a desigualdade é subaditiva, não uma
igualdade. -/
def TailData.epsilon (data : TailData) : Rat :=
  data.pencilBound * data.analyticBound * data.resolventBound *
    (data.fourierTail + data.chebyshevTail)

/-- Condições sem as quais a composição acima não significa nada: decaimentos
estritamente dentro de `(0,1)` e constantes não negativas. -/
def TailData.admissible (data : TailData) : Prop :=
  0 < data.decayFourier ∧ data.decayFourier < 1 ∧
  0 < data.decayChebyshev ∧ data.decayChebyshev < 1 ∧
  0 ≤ data.resolventBound ∧ 0 ≤ data.analyticBound ∧ 0 ≤ data.pencilBound ∧
  0 ≤ data.fourierTailFactor ∧ 0 ≤ data.chebyshevTailFactor

/-- Realização de Wiener (norma ℓ¹ com pesos nos coeficientes): a projeção de
truncamento tem norma `1` e a soma da cauda já está absorvida pelos pesos, de
modo que não sobram denominadores geométricos. -/
def wienerFactor : Rat := 1

/-- Realização de norma do supremo (Cauchy em `τ_R` + Bernstein em `ξ`):
`c_F = 4 / ((1-q_F)(1-q_C))`, vindo de `2` (Bernstein) `× 2` (dois sinais de
`n`) `× Σ_{n>N_F} q_F^{n-N_F-1} × Σ_{m≥0} q_C^m`. -/
def supFactorFourier (decayF decayC : Rat) : Rat :=
  4 / ((1 - decayF) * (1 - decayC))

/-- Companheira de `supFactorFourier`:
`c_C = 2(1+q_F) / ((1-q_F)(1-q_C))`, com `Σ_{n ∈ ℤ} q_F^{|n|} = (1+q_F)/(1-q_F)`. -/
def supFactorChebyshev (decayF decayC : Rat) : Rat :=
  2 * (1 + decayF) / ((1 - decayF) * (1 - decayC))

/-- Variante bilateral do truncamento de Galerkin, `C_N = Π C Π`:
`‖C - Π C Π‖ ≤ ‖(I-Π) C‖ + ‖Π‖ · ‖C (I-Π)‖`. O segundo termo **não** segue de
`P2`; ele exige a hipótese dual `P2*` (regularidade analítica do adjunto), cujas
constantes devem ser passadas em `adjointSide`.

`TailData.epsilon` sozinho é o truncamento unilateral `C_N = Π C`, que tem o
mesmo determinante de Fredholm finito (Lema 2.1 de `docs/TAIL_BOUND.md`) e só
consome `P2`.

**Atenção, isto mudou.** Este docstring dizia antes que o campo `tailError` do
`Tile` tinha de receber `twoSidedEpsilon`, e não `epsilon`. Aquilo era correto
enquanto `Tile.transferMajorant` usava a forma bilateral `max(1, ‖V‖/(1-η))`,
que só é majorante porque `I + λΠCΠ` age como identidade em `(ran Π)^⊥`.

O `Tile` foi reescrito desde então na forma **unilateral**
`maxDistance · tailError · (1 + maxDistance · inverseBound · operatorBound)`,
justamente para dispensar `P2*`. O tile passou portanto a corresponder a
`C_N = Π C`, que é o truncamento de `TailData.epsilon`. Hoje é **`epsilon`** que
casa com o tile, e é com ele que `ChoptuikFormal.TailTileBridge` ancora o
`tailError`. Usar `twoSidedEpsilon` continua **sólido** (a soma domina cada
parcela), apenas mais pessimista — e reintroduz a dependência de `P2*` que o
projeto eliminou por decisão. Ver `docs/TAIL_BOUND.md`, §2 e §6, e a linha `P2*`
de `docs/PROOF_OBLIGATIONS.md`. -/
def twoSidedEpsilon (rangeSide adjointSide : TailData) : Rat :=
  rangeSide.epsilon + adjointSide.epsilon

/-- Um pedido concreto: estas constantes produzem `ε_N` abaixo de `target`?

`target` é o que o restante do certificado exige, tipicamente
`target = 1 / (sup_Γ |z-z₀| · sup_Γ ‖(I+(z-z₀)C_N)⁻¹‖)`, para que o majorante
de transferência do tile fique `< 1`. -/
structure TailCertificate where
  data : TailData
  target : Rat
deriving Repr, DecidableEq

/-- Conteúdo matemático conferido para um pedido. -/
def TailCertificate.claims (certificate : TailCertificate) : Prop :=
  certificate.data.admissible ∧
  0 < certificate.target ∧
  certificate.data.epsilon ≤ certificate.target

private def TailCertificate.claimsDecidable
    (certificate : TailCertificate) : Decidable certificate.claims := by
  unfold TailCertificate.claims TailData.admissible
  infer_instance

/-- Verificador executável, exato em `Rat`. -/
def TailCertificate.valid (certificate : TailCertificate) : Bool :=
  @decide certificate.claims certificate.claimsDecidable

/-- O verificador não pode aceitar um pedido cujas desigualdades racionais sejam
falsas. -/
theorem TailCertificate.valid_sound
    (certificate : TailCertificate) (accepted : certificate.valid = true) :
    certificate.claims :=
  @of_decide_eq_true certificate.claims certificate.claimsDecidable accepted

/-- Consequência usada pelo resto do certificado: o `ε_N` composto respeita o
alvo. Continua condicional a `P1`, `P2` e `P3`, que este arquivo não prova. -/
theorem TailCertificate.epsilon_le_target_of_valid
    (certificate : TailCertificate) (accepted : certificate.valid = true) :
    certificate.data.epsilon ≤ certificate.target :=
  (certificate.valid_sound accepted).2.2

/-- Os decaimentos aceitos são de fato contrações: sem isto, `q^(N+1)` não é
cauda de nada. -/
theorem TailCertificate.decay_lt_one_of_valid
    (certificate : TailCertificate) (accepted : certificate.valid = true) :
    certificate.data.decayFourier < 1 ∧ certificate.data.decayChebyshev < 1 :=
  ⟨(certificate.valid_sound accepted).1.2.1,
   (certificate.valid_sound accepted).1.2.2.2.1⟩

/-! ## Regressões com os truncamentos reais do projeto

As constantes abaixo **não são certificadas**. Elas formam o conjunto de
hipóteses batizado `H-opt` em `docs/TAIL_BOUND.md`:

* `M_0 = 10`, `A_0 = 4`, `M_J = 1`;
* `q_F = q_C = 2/5`, escolhido para ser compatível com a razão de convergência
  observada no localizador em ponto flutuante (diagnóstico, jamais prova).

O alvo `1/150` corresponde a `sup_Γ |z-z₀| = 3` e `sup_Γ ‖(I+(z-z₀)C_N)⁻¹‖ = 50`
no lema de transferência de `docs/PROOF_OBLIGATIONS.md`. Também é hipótese.

O resultado é honesto e negativo para as malhas de hoje: com `H-opt`, as três
malhas efetivamente calculadas (4×12, 6×18, 8×24) **falham**, e a última falha
por pouco — fica acima do alvo por um fator maior que `157/100`. -/

/-- Constantes hipotéticas `H-opt`, realização de Wiener. -/
def optimisticGrid (fourier chebyshev : Nat) : TailData where
  resolventBound := 10
  analyticBound := 4
  pencilBound := 1
  decayFourier := 2 / 5
  decayChebyshev := 2 / 5
  fourierTailFactor := wienerFactor
  chebyshevTailFactor := wienerFactor
  fourierTruncation := fourier
  chebyshevTruncation := chebyshev

/-- Alvo hipotético imposto pelo lema de transferência de Neumann. -/
def neumannTarget : Rat := 1 / 150

def optimisticCertificate (fourier chebyshev : Nat) : TailCertificate where
  data := optimisticGrid fourier chebyshev
  target := neumannTarget

-- O valor exato composto para a malha 4×12.
example : (optimisticGrid 4 12).epsilon = 100065536 / 244140625 := by decide

-- As três malhas realmente calculadas no repositório não fecham sob `H-opt`.
example : (optimisticCertificate 4 12).valid = false := by decide
example : (optimisticCertificate 6 18).valid = false := by decide
example : (optimisticCertificate 8 24).valid = false := by decide

-- Quantificando a falha da melhor malha: `ε` excede o alvo por mais de 1,57×.
example : neumannTarget * (157 / 100) < (optimisticGrid 8 24).epsilon := by decide

-- Ainda sob `H-opt`, bastaria subir o truncamento de Fourier de 8 para 9.
example : (optimisticCertificate 9 24).valid = true := by decide
example : (optimisticCertificate 12 24).valid = true := by decide

/-- Mesmas constantes de EDP, mas na realização de norma do supremo, onde as
somas geométricas aparecem explicitamente: `c_F = 100/9` e `c_C = 70/9`. -/
def supGrid (fourier chebyshev : Nat) : TailData where
  resolventBound := 10
  analyticBound := 4
  pencilBound := 1
  decayFourier := 2 / 5
  decayChebyshev := 2 / 5
  fourierTailFactor := supFactorFourier (2 / 5) (2 / 5)
  chebyshevTailFactor := supFactorChebyshev (2 / 5) (2 / 5)
  fourierTruncation := fourier
  chebyshevTruncation := chebyshev

def supCertificate (fourier chebyshev : Nat) : TailCertificate where
  data := supGrid fourier chebyshev
  target := neumannTarget

example : supFactorFourier (2 / 5) (2 / 5) = 100 / 9 := by decide
example : supFactorChebyshev (2 / 5) (2 / 5) = 70 / 9 := by decide

-- A realização de norma do supremo é estritamente pior pelo fator `c_F = 100/9`:
-- 8×24 fica 17× acima do alvo (contra 1,57× na realização de Wiener) e 12×36
-- passa com folga de apenas 2,2×.
example : (supCertificate 8 24).valid = false := by decide
example : (supCertificate 12 36).valid = true := by decide
example : neumannTarget * 17 < (supGrid 8 24).epsilon := by decide

/-- Conjunto pessimista `H-pes`: mesma estrutura, decaimento de Fourier mais
lento (`q_F = 3/5`) e constante analítica maior (`A_0 = 40`). O requisito salta
de `N_F = 9` para algo como `N_F = 24`: uma mudança modesta na faixa de
analiticidade multiplica por vários o tamanho da matriz necessária. -/
def pessimisticGrid (fourier chebyshev : Nat) : TailData where
  resolventBound := 10
  analyticBound := 40
  pencilBound := 1
  decayFourier := 3 / 5
  decayChebyshev := 1 / 2
  fourierTailFactor := wienerFactor
  chebyshevTailFactor := wienerFactor
  fourierTruncation := fourier
  chebyshevTruncation := chebyshev

def pessimisticCertificate (fourier chebyshev : Nat) : TailCertificate where
  data := pessimisticGrid fourier chebyshev
  target := neumannTarget

example : (pessimisticCertificate 8 24).valid = false := by decide
example : (pessimisticCertificate 20 32).valid = false := by decide
example : (pessimisticCertificate 24 48).valid = true := by decide

/-! ## Regressões negativas: o verificador tem de recusar dados degenerados -/

/-- Decaimento `1` não é decaimento: a série não converge e o truncamento não
tem cauda somável. -/
def degenerateDecay : TailCertificate where
  data := { optimisticGrid 12 24 with decayFourier := 1 }
  target := neumannTarget

example : degenerateDecay.valid = false := by decide

/-- Alvo não positivo é recusado mesmo com `ε_N` minúsculo. -/
def degenerateTarget : TailCertificate where
  data := optimisticGrid 40 48
  target := 0

example : degenerateTarget.valid = false := by decide

/-- Constante analítica negativa é recusada (seria um "bound" sem significado
que baixaria `ε_N` artificialmente). -/
def degenerateConstant : TailCertificate where
  data := { optimisticGrid 12 24 with analyticBound := -4 }
  target := neumannTarget

example : degenerateConstant.valid = false := by decide

/-- A variante bilateral soma as duas caudas; com `H-opt` em 12×24 ela ainda
cabe no alvo, o que mostra que a diferença entre `Π C` e `Π C Π` não é o
gargalo — o gargalo é `P2`. -/
example : twoSidedEpsilon (optimisticGrid 12 24) (optimisticGrid 12 24) ≤
    neumannTarget := by decide

#print axioms TailCertificate.valid_sound
#print axioms TailCertificate.epsilon_le_target_of_valid
#print axioms TailCertificate.decay_lt_one_of_valid

end Choptuik.TailBound
