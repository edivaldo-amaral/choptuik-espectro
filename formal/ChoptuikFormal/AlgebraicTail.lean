import ChoptuikFormal.RationalCertificate

/-! Aritmetica da rota de cauda algebrica. Os bounds de operadores estao
demonstrados em papel em docs/ALGEBRAIC_TAIL.md, sob hipoteses explicitas. -/
namespace Choptuik.AlgebraicTail
open RationalCertificate
unseal Rat.mul Rat.add Rat.sub Rat.inv
set_option exponentiation.threshold 512

def perturbationMajorant : Rat :=
  (17/100) * (80/9) + 21 + 6 * ((1/2)^25 + (1/2)^277)

def freeBound (base : Rat) : Rat := 18 / (1/6 + base)
def resolventBound (base : Rat) : Rat :=
  freeBound base / (1 - 23 * freeBound base)
def outputTail (m n : Nat) : Rat :=
  18 * qmax (2 / (m : Rat)) (1 / (828 + (1/6) * ((n : Rat) + 1))) /
    (1 - 23 * freeBound 828)

theorem perturbation_lt_23 : perturbationMajorant < 23 := by decide
theorem baseNeumannMargin : 23 * freeBound 828 < (1/2 : Rat) := by decide
theorem baseResolventBound : resolventBound 828 = (108/2485 : Rat) := by decide
theorem rightBoundaryMargin : 23 * freeBound 414 = (2484/2485 : Rat) := by decide
theorem rightBoundaryAccepted : 23 * freeBound 414 < 1 := by decide
theorem outputTailPositive : 0 < outputTail 12 36 := by decide
theorem outputTailRefinement : outputTail 24 72 < outputTail 12 36 := by decide

/-- Cauda de ENTRADA para o precondicionador livre; requer m>=1,n>=2. -/
def inputTail (base : Rat) (m n : Nat) : Rat :=
  18 * qmax (2 / (m : Rat)) (1 / (base + (1/6) * ((n : Rat) - 1)))

/-- Primeiro fator do majorante, antes do custo da inversa finita. -/
def transferFloor (base distance : Rat) (m n : Nat) : Rat :=
  (23 + distance) * inputTail base m n

theorem actual4x12Floor : transferFloor 0 1 4 12 = (2592/11 : Rat) := by decide
theorem actual12x36Floor : transferFloor 0 1 12 36 = (2592/35 : Rat) := by decide
theorem allExistingMeshesFailFloor :
    ([ (4,12), (6,18), (8,24), (10,30), (12,36) ] : List (Nat × Nat)).all
      (fun mn => decide (1 < transferFloor 0 1 mn.1 mn.2)) = true := by decide

/-- Nem g=0 faria essas malhas passarem; esta fronteira e estrita. -/
theorem floorBoundary : transferFloor 0 1 864 2593 = 1 := by decide
theorem floorFirstIntegerPass : transferFloor 0 1 865 2594 < 1 := by decide

#print axioms perturbation_lt_23
#print axioms rightBoundaryAccepted
#print axioms allExistingMeshesFailFloor
end Choptuik.AlgebraicTail
