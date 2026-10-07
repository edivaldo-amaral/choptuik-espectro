import ChoptuikFormal.RationalCertificate

/-! Verificador escalar do majorante por blocos. A identificacao de a,b,c,d
com normas de blocos de um operador infinito e uma obrigacao externa.
Nao certifica essa identificacao nem a equivalencia com o quociente fisico. -/
namespace Choptuik.FrequencyBlocks
open RationalCertificate
unseal Rat.mul Rat.add Rat.sub Rat.inv

structure Bounds where
  a : Rat
  b : Rat
  c : Rat
  d : Rat
  weight : Rat
deriving Repr

def Bounds.inequalities (e : Bounds) : Prop :=
  0 ≤ e.a ∧ 0 ≤ e.b ∧ 0 ≤ e.c ∧ 0 ≤ e.d ∧ 0 < e.weight ∧
  e.a + e.weight * e.c < 1 ∧ e.d + e.b / e.weight < 1

instance (e : Bounds) : Decidable e.inequalities := by
  unfold Bounds.inequalities
  infer_instance

def Bounds.valid (e : Bounds) : Bool := decide e.inequalities

theorem Bounds.valid_sound (e : Bounds) (h : e.valid = true) : e.inequalities :=
  of_decide_eq_true h

def asymmetric : Bounds := ⟨1/10, 4, 1/100, 1/5, 10⟩
theorem asymmetricAccepted : asymmetric.valid = true := by decide
theorem unweightedRejected : ({ asymmetric with weight := 1 }).valid = false := by decide

/-- Bounds conservadores do experimento com inversas separadas. -/
def adjacentShell : Bounds := ⟨1/100000000, 441/100, 227/10, 1/100000000, 1⟩
theorem adjacentShellRejected : adjacentShell.valid = false := by decide
theorem adjacentCouplingDominates :
    (1-adjacentShell.a)*(1-adjacentShell.d) < adjacentShell.b*adjacentShell.c := by decide

/-- O fundo de referencia tem suporte finito, mas o fundo exato nao. -/
def backgroundError : Rat := (1/2)^25 + (1/2)^277
set_option exponentiation.threshold 512
theorem distantCouplingStillPositive : 0 < 6*backgroundError := by decide
theorem distantCouplingBound : 6*backgroundError < (1/5000000 : Rat) := by decide

#print axioms Bounds.valid_sound
#print axioms asymmetricAccepted
end Choptuik.FrequencyBlocks
