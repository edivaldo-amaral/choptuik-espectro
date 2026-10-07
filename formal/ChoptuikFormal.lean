import ChoptuikFormal.RationalCertificate
import ChoptuikFormal.WindingInterval
import ChoptuikFormal.TailBound
import ChoptuikFormal.DeterminantBound
import ChoptuikFormal.TailTileBridge
import ChoptuikFormal.WindingArcs
import ChoptuikFormal.RealTile
import ChoptuikFormal.TileBlocks
import ChoptuikFormal.ResolventCompression
import ChoptuikFormal.AlgebraicTail
import ChoptuikFormal.FrequencyBlocks

/-!
Um ledger formal da redução da contagem espectral.

Isto não postula a parte analítica que ainda falta. Em vez disso, expõe cada
obrigação como hipótese de um certificado. Assim, o kernel só permite concluir
"um modo" quando existência, equivalência física, Fredholm/cauda, winding e o
setor oposto tiverem sido fornecidos.
-/

namespace Choptuik

/-- Os dois autoespaços da involução de meio eco. -/
inductive HalfEchoSector where
  | same
  | opposite
deriving DecidableEq, Repr

/-- Multiplicidade algébrica física no semiplano instável, por setor. -/
structure SectorCount where
  same : Nat
  opposite : Nat
deriving DecidableEq, Repr

def SectorCount.total (count : SectorCount) : Nat :=
  count.same + count.opposite

/-- Obrigações matemáticas que um certificado espectral completo deve ligar
ao sistema Einstein--campo escalar, e não apenas a uma matriz truncada. -/
structure SpectralObligations (physicalMultiplicity : Nat) where
  backgroundExists : Prop
  backgroundCertified : backgroundExists
  gaugeFixedEquivalentToPhysicalQuotient : Prop
  gaugeCertified : gaugeFixedEquivalentToPhysicalQuotient
  analyticFredholmPencil : Prop
  fredholmCertified : analyticFredholmPencil
  quantitativeTailBound : Prop
  tailCertified : quantitativeTailBound
  noSpectrumOnCountingBoundary : Prop
  boundaryCertified : noSpectrumOnCountingBoundary
  sameSectorWinding : Nat
  sameWindingCertified : sameSectorWinding = 1
  oppositeSectorWinding : Nat
  oppositeWindingCertified : oppositeSectorWinding = 0
  countTransfer : physicalMultiplicity = sameSectorWinding + oppositeSectorWinding

/-- A conclusão final é puramente dedutiva depois que todas as obrigações
analíticas e computacionais do certificado forem preenchidas. -/
theorem exactlyOnePhysicalUnstableMode
    {physicalMultiplicity : Nat}
    (certificate : SpectralObligations physicalMultiplicity) :
    physicalMultiplicity = 1 := by
  calc
    physicalMultiplicity =
        certificate.sameSectorWinding + certificate.oppositeSectorWinding :=
      certificate.countTransfer
    _ = 1 + 0 := by
      rw [certificate.sameWindingCertified, certificate.oppositeWindingCertified]
    _ = 1 := rfl

/-! ## Versão certificada da contabilidade

Em `SpectralObligations` os dois windings são números naturais *asseverados*:
nada no tipo exige que venham de algum lugar. A versão abaixo substitui os dois
campos por certificados intervalares aceitos, de modo que o winding passa a ser
**calculado** a partir de caixas racionais verificadas pelo kernel
(`ChoptuikFormal.WindingInterval`), e não fornecido por um oráculo externo.

As obrigações genuinamente analíticas continuam sendo hipóteses — é isso que
elas são. O que deixou de ser hipótese é (i) a contagem discreta e (ii) o
controle de cauda somado à invertibilidade sobre o contorno, que agora entram
como um `JustifiedBoundary` aceito: cada tile satisfaz as desigualdades de
Neumann **e** tem o seu `tailError` ancorado no `epsilon_N` calculado das
constantes P1/P2/P3. -/

open WindingInterval in
/-- As mesmas obrigações, com a contagem por winding substituída por dois
certificados intervalares aceitos. -/
structure CertifiedSpectralObligations (physicalMultiplicity : Nat) where
  backgroundExists : Prop
  backgroundCertified : backgroundExists
  gaugeFixedEquivalentToPhysicalQuotient : Prop
  gaugeCertified : gaugeFixedEquivalentToPhysicalQuotient
  analyticFredholmPencil : Prop
  fredholmCertified : analyticFredholmPencil
  tailAndInvertibility : TailTileBridge.JustifiedBoundary
  tailAndInvertibilityAccepted : tailAndInvertibility.valid = true
  noSpectrumOnCountingBoundary : Prop
  boundaryCertified : noSpectrumOnCountingBoundary
  sameSectorCertificate : IntervalWindingCertificate
  sameSectorAccepted : sameSectorCertificate.valid = true
  sameSectorTarget : sameSectorCertificate.target = 1
  oppositeSectorCertificate : IntervalWindingCertificate
  oppositeSectorAccepted : oppositeSectorCertificate.valid = true
  oppositeSectorTarget : oppositeSectorCertificate.target = 0
  countTransfer : (physicalMultiplicity : Int) =
    sameSectorCertificate.total + oppositeSectorCertificate.total

/-- O winding de cada setor é lido do certificado intervalar, não postulado. -/
theorem sameSectorWindingIsOne
    {physicalMultiplicity : Nat}
    (obligations : CertifiedSpectralObligations physicalMultiplicity) :
    obligations.sameSectorCertificate.total = 1 := by
  rw [obligations.sameSectorCertificate.total_eq_target_of_valid
        obligations.sameSectorAccepted, obligations.sameSectorTarget]

theorem oppositeSectorWindingIsZero
    {physicalMultiplicity : Nat}
    (obligations : CertifiedSpectralObligations physicalMultiplicity) :
    obligations.oppositeSectorCertificate.total = 0 := by
  rw [obligations.oppositeSectorCertificate.total_eq_target_of_valid
        obligations.oppositeSectorAccepted, obligations.oppositeSectorTarget]

/-- O controle de cauda deixa de ser postulado: cada tile do contorno satisfaz
as desigualdades de Neumann e tem `tailError` justificado pelas constantes. -/
theorem tailIsJustifiedOnEveryTile
    {physicalMultiplicity : Nat}
    (obligations : CertifiedSpectralObligations physicalMultiplicity) :
    ∀ tile ∈ obligations.tailAndInvertibility.boundary.tiles,
      tile.inequalities ∧
        obligations.tailAndInvertibility.tailData.epsilon ≤ tile.tailError :=
  obligations.tailAndInvertibility.every_tile_justified
    obligations.tailAndInvertibilityAccepted

/-- Conclusão final a partir de certificados verificados pelo kernel. -/
theorem exactlyOnePhysicalUnstableModeCertified
    {physicalMultiplicity : Nat}
    (obligations : CertifiedSpectralObligations physicalMultiplicity) :
    physicalMultiplicity = 1 := by
  have transfer := obligations.countTransfer
  have same := sameSectorWindingIsOne obligations
  have opposite := oppositeSectorWindingIsZero obligations
  rw [same, opposite] at transfer
  omega

/-- Contabilidade separada impede esquecer o setor de simetria oposta. -/
theorem sectorCountIsOne
    (count : SectorCount)
    (sameCount : count.same = 1)
    (oppositeCount : count.opposite = 0) :
    count.total = 1 := by
  unfold SectorCount.total
  rw [sameCount, oppositeCount]

#print axioms exactlyOnePhysicalUnstableMode
#print axioms exactlyOnePhysicalUnstableModeCertified
#print axioms tailIsJustifiedOnEveryTile
#print axioms sectorCountIsOne

end Choptuik
