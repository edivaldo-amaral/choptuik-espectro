# T2: enunciado final (02/10/2026)

## Teorema

Seja ω* a solução crítica autossemelhante discreta (DSS) de Reiterer–Trubowitz (RT) para o colapso
**esfericamente simétrico** de um campo escalar sem massa, com expoente de eco μ*
(|μ* − μ_RefA| <= 3,9·10⁻¹¹). Seja L(s) o operador de RT linearizado em torno de ω*, no problema de
período completo, que é a soma direta dos setores A e B. Vale a modelagem física da seção seguinte
(simetria esférica, H-reg e identificação do fundo). H-gauge, H-rec e H-cone agora são lemas provados
(G, R e C).

1. **Contagem.** Em cada classe de Floquet do problema de período completo (faixa semiaberta de largura
   **1/2** em Im s), as raízes de L em Re s > 0, contadas com multiplicidade algébrica, são exatamente
   quatro, todas reais:

   | raiz | setor | multiplicidade | natureza |
   |---|---|---|---|
   | μ* ≈ 0,16831 | A | 1 | gauge: autoespaço span(g*), g* = (Z* + 1)ω*, Dg* = 0 |
   | s_K ≈ 0,40102 | A | 1 | viola as constraints |
   | ≈ 0,73318 | A | 1 | satisfaz as constraints em todo o espaço de raízes |
   | s_B ≈ 0,04142 | B | 1 | viola as constraints |

   Nas coordenadas de L_A, que é o que os certificados usam (Lema 5.1 de `SECTOR_B.md`), a raiz do setor
   B aparece em s_B + i/2 ≈ 0,0414 + 0,5i, e a faixa é Q̃ = {−1/4 < Im s <= 3/4}, de largura 1.

2. **Modos físicos.** Módulo gauge, e só com perturbações que satisfazem as constraints, o número de
   **modos instáveis físicos é 1**, a raiz real s ≈ 0,7332. Aqui "gauge" é o gauge linearizado padrão:
   L_X com X campo C¹ no domínio de RT.

   Se os gauges forem restritos aos que fixam o vértice e o cone, a contagem é 2. O modo extra, em μ*, é
   tangente às soluções exatas obtidas transladando o instante T* da singularidade, todas isométricas ao
   fundo (Lema C, `GAUGE_GLOBAL.md` §6).

3. **Eixo imaginário: só a fase, entre as raízes de L.** Em Re s = 0, a única raiz de L é s ≡ 0,
   algebricamente simples, com autovetor ∂_τω*: a translação de fase τ ↦ τ + c, de gauge, n = 1 na
   classificação de `GAUGE_DOMAIN_CLASSIFICATION.md`.

   Isto **não** é uma classificação dos modos físicos neutros. Em s = 0 há duas perturbações físicas que
   não são de gauge e que L não vê, porque δω = 0:
   - a escala global, δg = 2ḡ;
   - a translação do escalar, δφ = const.

   Ambas são simetrias exatas das equações. Os Lemas G e R valem só para Re s > 0. (Revisão leve,
   `.codex-runs/2026-10-03-revisao-CRH/`, L1.)

   *Prova (certificados existentes, `C1_REAVALIACAO.md` §13):*
   - Os Rouchés das faixas A e B cobrem a aresta Re s = 0 para Im s ∈ [−1/4, 3/4], com discos em que o
     operador deflacionado L̃ = L + E é invertível, pelas cotas t_p certificadas e θ < 1. Isso é um período inteiro de L_A, e Σ_A é invariante
     por s ↦ s + i.
   - Pela fórmula de multiplicidade da deflação de posto finito, mult_L̃(s) = mult_L(s) + ord_s(q/(s(s − μ))),
     e q não se anula em Re s >= 0.
   - Logo, para s ≠ 0 no eixo, mult_L(s) = mult_L̃(s) = 0. Em s = 0, mult_L(0) = mult_L̃(0) + 1 = 1, porque
     s = 0 está no interior de dois discos certificados, com folgas 0,018 e 0,005.

**Natureza do enunciado.** É um teorema **espectral**: modos de Floquet analíticos e, pelo Lema R e
pela volta, modos físicos de Floquet das equações linearizadas. Não afirma nada sobre o semigrupo de evolução nem sobre
estabilidade linear no sentido de evolução.

**Localizações.**
- μ*, s_K e s_B são certificadas por NK: raios 6·10⁻⁵ (μ*, S4), 3,8·10⁻⁴ (s_K, S3b) e 3,7·10⁻⁵ (s_B).
- 0,7332: NK de L (03/10, `build/nk733/C3.json`): **λ* ∈ [0,73316950; 0,73318982]**
  (|λ* − 0,7331796596606924| <= 1,02·10⁻⁵), algebricamente simples, θ = 0,806 e raio 2,2·10⁻⁵.
  Como N_A = 3 e μ* e s_K estão certificadas em outros discos, esta é a terceira raiz da faixa A.

**Realidade das raízes.** Pelo L-real, as raízes no toro são as de RT, onde a matriz é real. Lá vale
Σ_A = conj(Σ_A), e o setor B tem a simetria s ↦ conj(s) + i (`SECTOR_B.md` §7).
- Faixa A: cada raiz é simples e única na sua faixa, e a faixa é invariante pela conjugação. Logo a
  raiz é fixa pela conjugação, isto é, real.
- Faixa B: a faixa é invariante por s ↦ conj(s) + i, cujo conjunto fixo é Im s = 1/2. Logo
  Im = 1/2 em L_A, e s_B é real.

**Conferência externa** (não faz parte da prova). Com o período do eco de Choptuik, Δ ≈ 3,4453, o
tempo logarítmico é τ_log = (Δ/4π)τ. A raiz física dá λ = 4πλ*/Δ ∈ [2,67416; 2,67424], ou seja γ = 1/λ ≈ 0,3739,
o expoente crítico conhecido do campo escalar. A raiz de gauge dá 4πμ*/Δ ≈ 0,614: o expoente de um
modo de gauge depende das coordenadas, e não é o "1" do gauge de fatiamento usual.

## As hipóteses

### Antes hipóteses, agora provadas

- ~~(H-rec) Reconstrução~~: **as duas direções estão provadas.**
  - **Ida (Lema R, `HREC_IDA.md`, 03/10, revisado):** todo modo físico de Floquet com Re s > 0, na classe
    H-reg, entra no ansatz de RT por um gauge analítico, com h ∈ D(Y+), L(s)h = 0 e C(s)h = 0.
    - A ressonância s ≡ μ fecha pela simplicidade certificada de μ*.
    - O espaço espectral vem da recuperação de raio a partir de qualquer peso > 1, para todo s.
    - As cadeias físicas têm comprimento 1 módulo gauge.
  - **Volta (revisada em 03/10):** todo elemento de ker L(s) ∩ ker C(s) em D(Y+), com Re s > 0, vem de um
    modo físico de Floquet real-analítico em M̂ (`GEOMETRIC_FLOQUET_RECONSTRUCTION.md`, atualização de
    03/10, com o Lema V1: L e C juntos dão os dez δΩ nulos).
- ~~(H-gauge) Completude do gauge residual~~: **provada em 03/10** (Lema G, `GAUGE_GLOBAL.md`).
  - Para todo campo de gauge C¹ em M × S² que preserva o ansatz, os modos de gauge de Floquet com
    Re s > 0 formam ℂ·g*, em s ≡ μ.
  - A prova é global em todo o domínio de RT e não usa analiticidade dos germes.
  - Usa um fato certificado, F1: ω4*(·, 0) ≢ 0. O Lema G′ admite gauges só contínuos no cone, com δω
    analítico, usando F1′.
  - Revisado em 03/10 (`.codex-runs/2026-10-03-revisao-G/`): não quebrou.
  - Definição de gauge usada em T2: X de classe C¹ em M̂, ou contínuo com δω_X analítico, com μ e as
    coordenadas (τ, x) fixos.

### Modelagem física

- **Simetria esférica:** o resultado é sobre perturbações esfericamente simétricas.
- ~~(H-cone) Cone móvel~~: **não é mais hipótese** (Lema C, `GAUGE_GLOBAL.md` §6, 03/10).
  - Com o gauge linearizado padrão, g* é de gauge e a contagem é 1.
  - Com gauges que fixam o cone, a contagem é 2, e o modo extra é a translação de T*, uma simetria exata.
  - O teorema enuncia as duas.
- **(H-reg′) Regularidade:** a classe física de perturbações é **C^∞ em τ e real-analítica em ξ**, com
  vizinhança complexa uniforme em τ, de raio arbitrariamente pequeno (`HREC_IDA.md`, R5 com peso de
  Fourier 1, 03/10).
  - Pelo L-real e pela R5, a contagem não depende do raio.
  - Fica fora a regularidade só C^∞ em ξ, através do cone e no eixo: é um problema de regularidade
    aberto.
- **Identificação do fundo:** o fundo de RT, com μ e as coordenadas (τ, ξ) fixos, é a solução crítica de
  Choptuik.

## A cadeia da prova

| passo | afirmação | onde | natureza | revisão |
|---|---|---|---|---|
| 0 | existência de ω* na bola de RT | artigo de RT | prova assistida publicada | — |
| 1 | L(s) = L(0) + s·I; C(s) = C0 + sC1 | `T2_HIPOTESES.md` §3.1 (L-afim) | papel | — |
| 2 | redução B → A: Σ_B = Σ_A − i/2, constraints transportadas | `SECTOR_B.md` §2–5; L-iso, L-constr | papel | — |
| 3 | família de Fredholm analítica de índice 0, raízes isoladas | `RADIUS_RECOVERY.md`, F1 | papel, auditado 16/09 | parcial |
| 4 | contagens do toro = contagens na realização de RT | L-real | certificado racional + papel | — |
| 5 | nada em Re s >= 2,93 | F3 | Arb | — |
| 6 | N(L; faixa A, Re s > 0) = 3 (zeros de H̃ = L̃Q0(s) = zeros de L̃, porque Q0(s) é invertível e analítica em Re s > −μ; Brauer: N(L) = N(L̃) + 1) | `C1_REAVALIACAO.md` §6, §9–10 | certificado | revisado 01/10 |
| 7 | N(L; faixa B, Re s > 0) = 1 | §11 | certificado | revisado 02/10 |
| 8 | μ* simples, autoespaço span(g*), Dg* = 0 | `S4_GAUGE.md` | certificado + lema de Corr | revisado 26/09 |
| 9 | s_K simples, autovetor viola as constraints | `S3B_ROTA.md` | certificado | revisado 25/09 |
| 10 | raiz B simples, autovetor viola as constraints | §11 | certificado | revisado 02/10 |
| 11 | K injetivo em toda a faixa A com Re s >= 0, exceto em s_K: S3a (0 <= Re s <= 1,765 fora do disco de S3b), Rouché de K no disco (só s_K) e R♯ <= 1,765 (Re s >= 1,765). Logo vale na raiz ~0,7332, qualquer que seja sua posição exata. Os certificados de K estão no toro♯(129/128, 9/8), e C(s)h cai nesse espaço pela perda de raio de `SHARP_DOMAIN.md` de (65/64, 5/4) para (129/128, 9/8); ver o fechamento de 03/10 ali | `S3_CONDICIONAMENTO_SHARP.md` §9–10, `S3B_ROTA.md` | certificado | revisado 24 e 25/09 |
| 11a | 0,7332 simples, \|λ* − 0,73318\| <= 1,02·10⁻⁵ | `build/nk733/` (C3.json, MANIFEST.sha256) | certificado (NK bordejado) | mesma cadeia da faixa B, revisada 02/10; esta instância sem revisão própria |
| 11b | Lema G (e G′): os modos de gauge de Floquet com Re s > 0 são ℂ·e^{μτ}g* (s ≡ μ), para gauges C¹, ou contínuos com δω analítico | `GAUGE_GLOBAL.md`, `gauge_global_centro.py` (F1, F1′) | papel + racional | revisado 03/10 |
| 11c | Lema R (H-rec, ida): modos físicos de Floquet ↔ (ker L ∩ ker C)/gauge, com cadeias | `HREC_IDA.md`, `radius_recovery_geral.py` | papel + racional | revisado 03/10 |
| 11d | Lema C: e^{μτ}g* é tangente às translações de T*, soluções exatas isométricas ao fundo | `GAUGE_GLOBAL.md` §6 | papel | — |
| 12 | KC = ML; quociente dim E − rank(D\\|E) − [gauge] | `SHARP_LINEARIZED_IDENTITY.md`, `CONSTRAINT_ROOT_QUOTIENT.md` | papel | auditado 16/09 |

**Conta:**
- μ*: E = span(g*) é todo de gauge (Lema G; X0 é um campo analítico em M̂), logo
  dim E − rank(D|E) − 1 = 1 − 0 − 1 = 0. Com gauges que fixam o cone, seria 1 (Lema C).
- s_K e a raiz B: rank(D|E) = dim E = 1, então 0.
- 0,7332: rank(D|E) = 0, porque K é injetivo ali (passo 11, também ao longo de cadeias), e o autovetor
  não é de gauge (Lema G, s ≢ μ). Então 1.
- Total: **1**. Pelo Lema R (ida, com cadeias) e pela volta, é o número de modos instáveis físicos de
  Floquet, módulo gauge.

## O que não é coberto

- **A modelagem acima** (simetria esférica, H-reg, identificação do fundo), que não se reduz ao
  operador.
- **A volta da reconstrução** foi revisada em 03/10. A revisão apontou que C(s)h cai num espaço de raio
  menor; isso se fecha porque os certificados de K estão exatamente nesse espaço (`SHARP_DOMAIN.md`,
  fechamento).
- **Os passos 1–3 e 12 são provas em papel**, feitas em parte no interlúdio de 14–16/09 e auditadas em
  16/09, sem formalização em Lean.
- **Revisão independente de T2** (02–03/10, `.codex-runs/2026-10-02-revisao-T2/`): não quebrou o resultado. O L-real, a redução B → A, o quociente e KC = ML conferem. As lacunas de enunciado e de citação (1–11) estão corrigidas neste documento e em `T2_HIPOTESES.md`.
- **A localização de 0,7332 (passo 11a) não passou por revisão própria.** Usa a mesma cadeia NK da raiz B,
  que foi revisada; para T2 basta a contagem, porque K é injetivo em toda a faixa A menos s_K (passo 11).
