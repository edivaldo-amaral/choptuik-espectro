import ChoptuikFormal.RationalCertificate

/-! Contraexemplo exato: comprimir a inversa nao e inverter a compressao.
Ele impede identificar a resolvente finita com o bloco da resolvente total
sem uma estimativa adicional dos acoplamentos. -/
namespace Choptuik.ResolventCompression
open RationalCertificate
unseal Rat.mul Rat.add Rat.sub Rat.inv

def pencil : RatMatrix := ⟨2, 2, [2, 1, 1, 2]⟩
def resolvent : RatMatrix := ⟨2, 2, [2/3, -1/3, -1/3, 2/3]⟩
def compressedPencil : RatMatrix := ⟨1, 1, [2]⟩
def finiteResolvent : RatMatrix := ⟨1, 1, [1/2]⟩

theorem fullInverse : pencil.mul resolvent = RatMatrix.identity 2 := by decide
theorem finiteInverse :
    compressedPencil.mul finiteResolvent = RatMatrix.identity 1 := by decide
theorem compressionMismatch :
    resolvent.entry 0 0 ≠ finiteResolvent.entry 0 0 := by decide
theorem missingError :
    resolvent.entry 0 0 - finiteResolvent.entry 0 0 = (1/6 : Rat) := by decide

#print axioms compressionMismatch
#print axioms missingError
end Choptuik.ResolventCompression
