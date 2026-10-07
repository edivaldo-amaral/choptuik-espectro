# Perturbações não esféricas do espaço de Choptuik (Fases 0 e 1, 03/10/2026)

**Objetivo.** Estender T2 além da simetria esférica. O alvo é um teorema T3: nenhum modo instável com
ℓ >= 1, módulo gauge. Junto com T2, isso daria "exatamente 1 modo instável" para perturbações lineares
gerais, no sentido espectral.

**Por quê.** RT (artigo publicado como "Choptuik's critical spacetime exists", Comm. Math. Phys. 2019)
dizem que os resultados numéricos sobre perturbações "parecem inconclusivos". Citam:
- Martín-García & Gundlach (MG), Phys. Rev. D 59, 064031 (1999), arXiv gr-qc/9809059: análise linear
  numérica, "all nonspherical perturbations decay";
- Choptuik, Hirschmann, Liebling & Pretorius (2003): colapso axissimétrico não linear;
- Gundlach & Martín-García, Living Rev. Rel. 10 (2007), §3.7 e §5.3.

**Dados de MG** (Tabela 1; λ = κ + iω, período Δ ≈ 3,445; convertemos com s = λK/(2π), K = Δ/2):

| setor | κΔ (MG) | ωΔ/2π (MG) | Re s |
|---|---|---|---|
| ℓ = 0 par | +9,21 | 0 | +0,733 (a nossa raiz física) |
| ℓ = 1 par | −0,30 | 0 | −0,024 |
| **ℓ = 2 par** | **−0,07** | 0,3 | **−0,0056** |
| ℓ = 2 ímpar | −2,30 | 1,9 | −0,183 |
| ℓ = 3 par | −1,65 | 1,6 | −0,131 |
| ℓ = 3 ímpar | −3,28 | 1,4 | −0,260 |

**O ponto crítico é o ℓ = 2 par,** a 0,006 do eixo. É por isso que a questão é "inconclusiva": as
diferenças finitas de MG, com até 1600 pontos, não excluem com folga um sinal trocado.

## Fase 0: o setor axial no fundo de RT

Equação de Gerlach–Sengupta (MG, eq. 57, sem fonte para o campo escalar):
[r⁻²(r⁴Π)^{|A}]_{|A} = (ℓ − 1)(ℓ + 2)Π.

No fundo de RT valem:
- g^{AB} = e^{2ζ}·½(D⁺⊗D⁻ + D⁻⊗D⁺) e √|g| = e^{−2ζ}/μ;
- a divergência de c·D^σ é (D^σ + σμ)c;
- D^σ log r = μ/ξ − ω2^σ, que **não depende de ζ**.

Com G := r²Π, a equação fica

    (D⁺ + μ)(D⁻G + 2v⁻G) + (D⁻ − μ)(D⁺G + 2v⁺G) = 2(ℓ − 1)(ℓ + 2)(W/ξ)²G,

com v^σ = μ/ξ − ω2^σ e W = μ + ξ(ω1 − ω2⁻ − ω2⁺)/2. **Só entram ω e μ.** No centro, a análise indicial dá
G ~ ξ^ℓ (o outro expoente é −ℓ − 1). Minkowski é ω = 0, W = μ.

## Fase 1: numérica (não rigorosa), `scripts/ns_axial.py`

**Discretização.** Floquet G = e^{sτ}g, problema quadrático K0 + sK1 − 2s² = 0. Colocação em Fourier (τ,
período 4π, N_τ ímpar) e Chebyshev–Gauss em ξ ∈ [−ξ*, ξ*].

**Três cuidados, sem os quais aparecem autovalores espúrios com Re s > 0:**
1. Regularidade no centro imposta pela base G = ξ^ℓ·Σ a_k T_{2k}(ξ/ξ*).
2. **O cone dentro do domínio**: ξ* = 1,02 < 41/40. Com o domínio terminando no cone, a colocação aproxima
   as soluções singulares (1 − ξ)^α, um espectro contínuo de espaços não analíticos, cujo limiar o fundo
   desloca para Re s > 0.
3. Filtro de suavidade dos autovetores, pela fração dos coeficientes de Chebyshev na metade superior.

**Teste no plano.** s = −2μ, −3μ, −4μ, −5μ, exatos, sem espúrios.

**Fundo RefA:**

| ℓ (ímpar) | resolução | s menos amortecido | κΔ | MG |
|---|---|---|---|---|
| 2 | 13×32 / 17×40 | −0,18469 ± 0,0410i / −0,18459 ± 0,0414i | −2,320 | −2,30 |
| 3 | 17×40 | −0,26194 ± 0,1824i | −3,292 | −3,28 |
| 4 | 17×40 | −0,34144 ± 0,1360i | −4,291 | — |

ωΔ/2π = 2·Im s, módulo 1, bate com o de MG módulo 1.

**Conclusão parcial.** A formulação axial no fundo de RT reproduz MG a ~1%. O setor axial não tem nada
perto do eixo.

## Fase 0 e 1 do setor polar (par), ℓ >= 2

**Derivação** (`scripts/ns_polar_deriva.py`, sympy). Partimos das equações covariantes de Gerlach–Sengupta
pares, que estão no fonte TeX de MG: (inveqAB), (inveqAeven), (inveqwave) e (inveq3), com T² = 0, logo
k_AB sem traço.
- **Fundo de RT:** g_AB = e^{−2ζ}[[−(1 − ξ²), −ξ/μ], [−ξ/μ, 1/μ²]] e r = e^{−ζ}ξ/W. A normalização de RT
  entra trocando φ → φ/(2√π).
- **Incógnitas:** k_AB = e^{−2ζ}(A θ⁺θ⁺ + B θ⁻θ⁻), com θ^σ as formas nulas, B = (−1)^ℓ·PA, mais k e Φ.
- **Conferências:**
  - as seis equações são homogêneas em e^ζ, logo só ∂ζ entra, via ω3;
  - a reconstrução a partir dos coeficientes é exata;
  - a simetria P (EAm(ξ) = −EAp(−ξ), e a paridade da onda, do traço e dos vínculos) vale a 10⁻¹⁵.

  Um bug de substituição (∂²ζ zerado por `subs` sequencial) foi pego e corrigido com `xreplace`.
- **Fundo numérico:** consistência ∂_τ∂_ξ pelos dois caminhos a 10⁻⁹.

**Numérica** (`scripts/ns_polar.py`, `scripts/ns_polar_alvo.py`):
- subsistema de evolução: onda de Φ, traço e EAm, como em MG;
- colocação com paridade no centro e cone interior (ξ* = 1,02);
- *shift-invert* perto do alvo;
- **filtro pelos vínculos** (cpp, cmm), indispensável. O subsistema tem modos espúrios que violam os
  vínculos (resíduo ~1), inclusive com Re s > 0 e um em +0,0020 + 0,1475i, colado no modo físico. Eles
  mudam com a resolução.
- No plano: s = −2μ, −3μ, … exatos.

**Fundo RefA** (`build/naoesf/l2_alvo*.txt`, `l3_alvo.txt`; laboratorio2):

| ℓ (par) | resolução | s | κΔ | ωΔ/2π | vínculo | MG |
|---|---|---|---|---|---|---|
| 2 | 25×40 | −0,0047249 ± 0,146891i | −0,0594 | 0,294 | 0,84 | |
| 2 | 29×48 | −0,0047172 ± 0,147424i | −0,0593 | 0,295 | 0,53 / 0,24 | |
| 2 | 33×56 | −0,0047800 ± 0,147342i | −0,0601 | 0,295 | 0,18 / 0,055 | |
| 2 | 37×64 | **−0,0047637 ± 0,147345i** | **−0,0599** | **0,2947** | 0,18 / 0,024 | −0,07; 0,3 |
| 3 | 33×56 | −0,131735 ± 0,30747i | −1,655 | 0,615 | 0,10 / 0,08 | −1,65; 1,6 |
| 3 | 33×56 | −0,26742 ± 0,31441i | −3,361 | — | 0,016 | — |

Os dois números de "vínculo" são o modo e o seu alias em Im s + 1/2.

**Conclusão da Fase 1.** A formulação no fundo de RT reproduz MG nos dois setores, axial e polar, para
ℓ = 2 e 3. **O modo polar ℓ = 2 é amortecido, mas fica a só 0,0048 do eixo.** Com espectral, a precisão é
bem maior que a de MG: o valor estabiliza em ~2·10⁻⁵.

**Consequências para a fase rigorosa:**
1. **Localização.** Um NK no modo ℓ = 2 com raio ≪ 0,0048 (os NK esféricos deram 10⁻⁵ a 4·10⁻⁴), mais
   deflação de Brauer, mais Rouché/invertibilidade em Re s >= 0.
2. **Vínculos.** A estrutura de vínculos tem de entrar no certificado, como KC = ML no caso esférico,
   porque o subsistema de evolução tem raízes que violam os vínculos.
3. **ℓ grande.** É preciso uma cota uniforme.

O tamanho comparável é o do projeto esférico inteiro.

## Próximos passos

1. ~~Setor polar~~ e ~~ℓ = 2 com precisão~~: **feitos** (acima). Re s(ℓ = 2) = −0,00476.
2. **Formulação bem condicionada** para a fase rigorosa: Galerkin nos coeficientes de RT com o inverso livre
   (como L), e os vínculos.
3. **ℓ = 1:** gauge (translação do centro) e o modo físico.
4. **Cota uniforme para ℓ grande.**
