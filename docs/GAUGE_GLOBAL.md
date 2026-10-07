# Completude global do gauge residual: Lema G (03/10/2026)

Substitui a hipótese H-gauge de `T2_ENUNCIADO.md` por um lema provado. `GAUGE_RT_ACTION.md` provava a
completude **local**, para germes analíticos. Aqui a prova é **global** em todo o domínio de RT, vale para
campos de gauge só de classe C¹ e usa um único fato certificado (F1).

**Revisado em 03/10** (`.codex-runs/2026-10-03-revisao-G/`): o revisor não quebrou o lema e refez F1 com
código próprio, chegando às mesmas frações exatas. As ressalvas R1–R7 estão incorporadas abaixo, e o
reforço que ele propôs virou o Lema G′ (§2.1).

## 1. Objetos

- **Domínio.** M = {(u_-, u_+): u_+ < 0, u_+ < u_- < −u_+/81}, o domínio do Teorema 1 de RT. Em (τ, ξ) é
  0 < ξ < 41/40, τ ∈ ℝ. A variedade física é M × S², com o eixo u_+ = u_- (ξ = 0) como pontos regulares.
  - **M contém o cone u_- = 0 (ξ = 1).**
  - Θ(u_-, u_+) = e^{−2πμ}(u_-, u_+), isto é, τ ↦ τ + 2π.
- **Fundo.** (ḡ, φ*) é real-analítico em M (ḡ é a métrica; g* fica reservado para o perfil de gauge
  (Z + 1)ω*, como em S4 e T2), com Θ*ḡ = e^{−2K}ḡ e φ*∘Θ = −φ*. Na forma de RT,
  g = e^{−2ζ}(−2(du_-du_+ + du_+du_-)/(μ²(u_+ + u_-)²) + ξ²(μ + Qξ²)⁻² g_{S²}).
- **Variáveis de RT.** Ω: (g, φ) na forma do ansatz ↦ ω = (ω1, ω2^±, ω3^±, ω4^±), definidas pelas
  equações (ekjehee) de RT. Valem:
  - Ω(e^{2c}g, φ) = Ω(g, φ) para c constante (ζ ↦ ζ − c);
  - a translação Θ preserva o ansatz;
  - ω4^σ = D^σφ, com D^σ = σ∂_τ + μ(1 + σξ)∂_ξ **fixos** (dependem só de μ).

  Na codificação de RT, ω_i^- = ω_i e ω_i^+ = −Pω_i (i = 2, 3, 4), com (Pf)(τ, ξ) = f(τ, −ξ).
- **Gauge.**
  - Uma perturbação linear (δg, δφ) é de gauge se δg = L_X ḡ e δφ = L_X φ* = X(φ*), para um campo
    vetorial X de classe C¹ em M × S², incluindo o eixo, possivelmente complexo.
  - X é **admissível** se (δg, δφ) está na forma linearizada do ansatz **com μ = μ\* e as coordenadas
    (τ, x) fixas**. Nesse caso, δω_X := DΩ(ḡ, φ*)[L_X ḡ, L_X φ*].
  - μ fixo é essencial: com δμ ≠ 0, as direções nulas ker(σdξ − μ(1 + σξ)dτ) mudam, e δg_{±±} deixa de
    se anular nos u_± antigos (revisão, R1). Em T2, μ e as coordenadas fixos fazem parte da identificação
    do fundo.
  - δω_X é um **modo de Floquet** com expoente s se δω_X∘Θ² = e^{4πs}δω_X, que é a condição para modos
    e^{sτ}h com h de período 4π, cobrindo os setores A e B.

## 2. Enunciado

**Lema G.** Seja X admissível e δω_X ≠ 0 um modo de Floquet com Re s > 0. Então e^{4π(s−μ)} = 1 e
δω_X = c·G*, com c ≠ 0 e G* := δω_{X0} = e^{μτ}g*. Aqui g* = (Z + 1)ω* é o **perfil** (notação de S4 e
T2), e X0 = ∂_{u_-} + ∂_{u_+} é a translação conjunta dos nulos (`GAUGE_RT_ACTION.md`).

Ou seja, **o espaço dos modos de gauge de Floquet com Re s > 0 é ℂ·e^{μτ}g\***, um único modo, que tem
expoente μ e pertence ao setor A.

### 2.1 Lema G′: gauge só contínuo no cone

**Lema G′** (proposto pela revisão). Valem as mesmas conclusões se X for apenas **contínuo** em M × S², de
classe C¹ fora do cone, e δω_X for analítico na faixa de RT, que é a classe de T2. Aqui usa-se o fato
F1′: ω4^+(·, 1) = −ω4(·, −1) ≢ 0.

*Prova.*
- **Passo 2.** Continua valendo. {u_+ = w} ∩ M atravessa o cone, e X^{u_+} é constante de cada lado e
  contínuo, logo constante na curva toda. Do mesmo modo, f_- é contínua em 0.
- **Passos 3–6.** Não usam diferenciabilidade no cone. A injetividade usa só a continuidade de f.
- **Passo 7, Re s >= μ.** Só usa continuidade.
- **Passo 7, 0 < Re s < μ.**
  1. δφ = X(φ*) é contínuo, e D^σδφ = δω4^σ é analítico, logo δφ é analítico.
  2. No eixo, −(e^{μτ}/μ)f(−e^{−μτ})ω4(τ, 0) é analítico. ω4(·, 0) tem zeros isolados e f é contínua, logo
     as singularidades são removíveis e f é analítica em (−∞, 0).
  3. Perto de um ponto do cone (τ0, 1) com ω4^+(τ0, 1) ≠ 0, que existe por F1′,
     f(u_-) = [2μe^{−μτ}δφ + f(u_+)ω4^-]/ω4^+ é analítica. Logo f é analítica em u = 0.
  4. f(Λu) = ρf(u) dá a_n(Λ^n − ρ) = 0 para os coeficientes de Taylor. Então f = a_n u^n perto de 0 e,
     pela homogeneidade, em todo ℝ, com s ≡ μ(1 − n). Re s > 0 força n = 0.

  ∎

Assim a classe de gauge exigida pela H-rec (ida) pode ser tanto "C¹ em M̂" (Lema G) quanto "contínuo, com
δω_X na classe analítica" (Lema G′). Os gauges que `HREC_IDA.md` constrói são analíticos, e estão nas
duas.

**Consequências para T2:**
- Nas raízes s_K, s_B e ≈ 0,7332, nenhum vetor não nulo do espaço de raízes é de gauge, porque s ≢ μ.
- Em μ*, o espaço de raízes, span(g*) por S4, é todo de gauge.
- Como as quatro raízes são simples, não há cadeias de Jordan a examinar.

## 3. Prova

**Passo 0 (simetria esférica).** Se δ = L_X(ḡ, φ*) é SO(3)-invariante, então δ = L_{X̄}(ḡ, φ*), com X̄ a
média de X sobre SO(3) pela medida de Haar, porque ḡ e φ* são invariantes e L comuta com a média.
- X̄ é C¹ e SO(3)-invariante.
- Os campos invariantes em S² são nulos, logo X̄ = X^τ∂_τ + X^ξ∂_ξ.
- No eixo, o único vetor invariante é a direção τ, logo X̄ é tangente ao eixo.

Daqui em diante X = X̄.

**Passo 1 (ansatz ⇒ equações de transporte).** Em coordenadas nulas, g_{++} = g_{--} = 0 e
g_{+-} = −Ω² ≠ 0. Então (L_X g)_{++} = 2g_{-+}∂_{u_+}X^{u_-} e (L_X g)_{--} = 2g_{+-}∂_{u_-}X^{u_+}.
- A forma do ansatz exige δg_{++} = δg_{--} = 0, porque a parte 2D tem de continuar conforme a
  du_-du_+.
- Logo ∂_{u_+}X^{u_-} = 0 e ∂_{u_-}X^{u_+} = 0 em M.
- As componentes angulares e mistas não restringem nada para X esfericamente simétrico.
- A forma do ansatz exige mais do que δg_{±±} = 0, por exemplo a regularidade de δQ no eixo. O lema só
  usa essa condição, que é necessária (revisão, R5).

**Passo 2 (globalidade).** As curvas de nível são conexas em M:
- {u_- = v} ∩ M é o intervalo u_+ ∈ (−∞, min(v, 0, −81v)) da reta u_- = v;
- {u_+ = w} ∩ M é o intervalo u_- ∈ (w, −w/81).

Logo X^{u_-} = f_-(u_-) e X^{u_+} = f_+(u_+) em todo M, com f_- ∈ C¹(ℝ), porque u_- percorre ℝ em M,
inclusive v = 0 no cone, e f_+ ∈ C¹((−∞, 0)). A regularidade C¹ vem de restringir X a uma reta
τ = τ0, que é transversal e onde u_- é monótona em ξ.

**Passo 3 (eixo).** X é tangente ao eixo u_+ = u_-, logo X(u_+ − u_-) = 0 ali, e f_+(u) = f_-(u) para todo
u < 0, que é a imagem do eixo. Seja f := f_-; então f_+ = f|_{(−∞,0)}. Em (τ, ξ):

    δτ = e^{μτ}[f(u_+) + f(u_-)]/(2μ),   δξ = e^{μτ}[(1 + ξ)f(u_-) − (1 − ξ)f(u_+)]/2.

**Passo 4 (equivariância), em versão linear, sem fluxos** (revisão, R3).
1. Θ preserva o ansatz: em (τ, x), Θ é τ ↦ τ + 2π, e o frame tem coeficientes independentes de τ. Logo
   Ω(Θ*(g, φ)) = Ω(g, φ)∘Θ. Derivando,
   DΩ(ḡ, φ*)[δ]∘Θ² = DΩ(e^{−4K}ḡ, φ*)[(Θ²)*δ].
2. Pela invariância conforme, DΩ(e^{2c}g, φ)[e^{2c}δg, δφ] = DΩ(g, φ)[δg, δφ].
3. (Θ²)*(L_X ḡ) = e^{−4K}L_{(Θ²)*X}ḡ, e (Θ²)*(X(φ*)) = ((Θ²)*X)(φ*), porque φ*∘Θ² = φ*.
4. Juntando, δω_X∘Θ² = δω_{(Θ²)*X}, e (Θ²)*X é C¹ e admissível.

Como u∘Θ² = Λu, com Λ = e^{−4πμ}, o campo (Θ²)*X tem X^{u} = f̃(u) = Λ⁻¹f(Λu). Com Θ² não aparece o
sinal de φ.

**Passo 5 (injetividade: δω_X = 0 ⇒ X = 0).** Usa só a componente ω4.
- δω4^σ = D^σ(X(φ*)) = 0 para σ = ± implica X(φ*) = c2, constante em M, que é conexo.
- Como D^+(u_+) = 0, D^+(u_-) = 2μe^{−μτ}, D^-(u_-) = 0 e D^-(u_+) = −2μe^{−μτ}, vale

      X(φ*) = (e^{μτ}/2μ)·[f(u_-)ω4^+ − f(u_+)ω4^-].

- **Eixo.** Em ξ = 0, u_± = −e^{−μτ}, ω4^-(τ, 0) = ω4(τ, 0) e ω4^+(τ, 0) = −ω4(τ, 0). Logo

      −(e^{μτ}/μ)·f(−e^{−μτ})·ω4(τ, 0) = c2   para todo τ.

  - ω4 é antiperiódica: ω4(τ + 2π, 0) = −ω4(τ, 0). Pelo valor intermediário, ω4(·, 0) se anula em algum
    τ, logo c2 = 0.
  - Por **F1**, ω4(·, 0) ≢ 0. Como é real-analítica em τ, seus zeros são isolados. Então f = 0 num
    subconjunto denso de (−∞, 0) e, por continuidade, em todo (−∞, 0).
- **Exterior do cone.** Em M ∩ {ξ > 1}, u_+ < 0, logo f(u_+) = 0 e f(u_-)ω4^+ = 0.
  - Se f(v) ≠ 0 para algum v > 0, então f ≠ 0 num intervalo I ∋ v, e ω4^+ se anula no aberto não
    vazio {ξ > 1, u_- ∈ I} ∩ M.
  - Por analiticidade, ω4^+ ≡ 0 em M. No eixo isso dá ω4(·, 0) ≡ 0, o que contradiz F1.
  - Logo f = 0 em (0, ∞) e, por continuidade, f(0) = 0.

Para X complexo, aplica-se o mesmo às partes real e imaginária, porque δω_{Re X} e δω_{Im X} são reais.

**Passo 6 (Floquet ⇒ homogeneidade).** Se δω_X∘Θ² = e^{4πs}δω_X, os Passos 4 e 5 dão
(Θ²)*X = e^{4πs}X, isto é,

    f(Λu) = ρ·f(u)  para todo u ∈ ℝ,   ρ = e^{4π(s−μ)},   |ρ| = e^{4π(Re s − μ)}.

**Passo 7 (regularidade no cone ⇒ f constante).** f é C¹ em u = 0, porque o cone está dentro de M.
- **Re s > μ, ou |ρ| > 1.** f(Λ^j u) = ρ^j f(u). O lado esquerdo tende a f(0), e |ρ^j| → ∞, logo f(u) = 0.
- **Re s = μ, ou |ρ| = 1.**
  - Se ρ ≠ 1, a sequência ρ^j não converge, porque |ρ^{j+1} − ρ^j| = |ρ − 1|. Logo f(u) = 0.
  - Se ρ = 1, f(u) = f(Λ^j u) → f(0), e f é constante.
- **0 < Re s < μ, ou Λ < |ρ| < 1.** f(0) = lim ρ^j f(u) = 0. Além disso,
  f'(0) = lim f(Λ^j u)/(Λ^j u) = lim (ρ/Λ)^j f(u)/u, com |ρ/Λ| > 1. Logo f(u) = 0.

Assim, para Re s > 0, ou f ≡ 0, o que daria δω_X = 0 e foi excluído, ou ρ = 1 e f ≡ c.

**Passo 8.** f ≡ c dá X = cX0 e δω_X = c·G*. O vetor G* ≠ 0 pelo Passo 5. A ação de X0 nos campos de RT
foi derivada e testada em `GAUGE_RT_ACTION.md` (`test_gauge_action.py`). ∎

## 4. O fato certificado F1

**F1: ω4*(·, 0) ≢ 0.** Na base de RT, ω(τ, ξ) = Σ v_mn e^{imτ/2 + inθ}, com θ = arccos ξ. O coeficiente de
Fourier m = 1 de ω4(·, 0) é c_1 = Σ_{n∈ℤ} v_{1n}i^n.
- Na bola de RT, ω* = RefA + RefB + Corr com ‖RefB‖ <= 2⁻²⁵ e ‖Corr‖ <= 2⁻²⁷⁷. A norma tem pesos >= 1 e
  ‖z‖_C >= |z|, logo |c_1(ω*) − c_1(RefA)| <= ε = 2⁻²⁵ + 2⁻²⁷⁷.
- Em racionais exatos (`scripts/gauge_global_centro.py`, `build/t2/gauge_global_centro.json`):
  c_1(RefA4) = 0,072916 + 0,493046i, logo **|c_1(ω4*(·, 0))| >= 0,493046 − 3·10⁻⁸ > 0**.
- O script confere, com um assert, que o campo de índice 3 de `RefA.dat` só tem m ímpar, isto é, que é ω4.
- **F1′** (para o Lema G′): em ξ = −1 (θ = π), c_1 = Σ_n v_{1n}(−1)^n = −0,079473 + 0,097753i, logo
  |c_1(ω4*(·, −1))| >= 0,097753 − 3·10⁻⁸ > 0.
- **Revisão:** recalculado com código próprio, com as mesmas frações exatas. Pela rota via `RefAplusB.dat`
  direto, |c_1(ω*)| >= 0,4930463216. ‖RefB‖ = 2,17·10⁻⁸ <= 2⁻²⁵.

## 5. Observações

1. **O papel do cone.** O lema usa que o cone u_- = 0 está dentro do domínio e que o gauge é C¹ através
   dele.
   - Se o gauge só precisasse ser regular no interior {ξ < 1}, f(u) = (−u)^{1−s/μ} daria um "modo de
     gauge" para todo s.
   - Heuristicamente, é a mesma regularidade através do cone que seleciona o espectro. Um campo singular
     no cone não é um difeomorfismo de M; o Lema G′ mostra que basta a continuidade, com δω regular.
2. **Cone fixo e cone móvel: ver o Lema C (§6).** Fixar o cone significa f(0) = 0. Isso exclui f
   constante, e então não há nenhum modo de gauge com Re s > 0. O §6 mostra que e^{μτ}g* é tangente a
   soluções exatas isométricas ao fundo, de modo que as duas contagens são teoremas.
3. **Simetrias que não são difeomorfismos:** a escala constante ζ ↦ ζ + c e a translação φ ↦ φ + c.
   - Elas estão no núcleo de DΩ: δω = 0 implica δQ = 0, δζ constante e δφ constante, pelas equações
     (dhjhfkfhfe1–4) de RT.
   - Não geram modos com Re s > 0.
4. **Regularidade usada.**
   - X é C¹ em M × S², o que é a definição de difeomorfismo.
   - Para Re s >= μ, que inclui a raiz física 0,7332, basta X contínuo no cone. Para todo Re s > 0 basta
     X contínuo com δω_X analítico (Lema G′).
   - Nenhuma analiticidade dos germes é usada: isso fecha "a exclusão de germes não analíticos".
5. **O que o Lema G não cobre: H-rec (ida).** O lema prova a unicidade do gauge até ℂg*. Não prova que toda
   perturbação física possa ser posta no ansatz, que é a existência da fixação de gauge.

## 6. Lema C: o modo de μ* é uma simetria exata (03/10/2026)

Seja ψ_ε(u_-, u_+) = (u_- + ε, u_+ + ε), a translação conjunta das coordenadas nulas. Ela é gerada por
X0 = ∂_{u_-} + ∂_{u_+} = e^{μτ}(μ⁻¹∂_τ + ξ∂_ξ).

**Lema C.**
1. Para cada ε, (ḡ_ε, φ_ε) := ψ_ε*(ḡ, φ*) é uma solução **exata** das equações de Einstein–campo escalar
   em ψ_ε⁻¹(M), **isométrica** a (ḡ, φ*) via ψ_ε.
2. Essa solução está na forma do ansatz de RT, com o mesmo μ e as mesmas coordenadas (τ, ξ). Ela é
   autossemelhante discreta em torno do vértice deslocado (−ε, −ε), isto é, é a mesma solução de
   Choptuik com o instante T* da singularidade transladado.
3. As variáveis de RT ω_ε := Ω(ḡ_ε, φ_ε) resolvem todas as equações Ω_i^σ = 0, e
   ω_ε = ω* + ε·e^{μτ}g* + O(ε²), uniformemente em compactos de M.

Logo e^{μτ}g* é o vetor tangente a uma curva de soluções exatas, todas isométricas ao fundo.

*Prova.*
- **(1)** As equações são invariantes por difeomorfismos, e ψ_ε é um difeomorfismo de ℝ² que se estende
  trivialmente a ℝ² × S².
- **(2)** ψ_ε preserva du_± e o eixo, porque u_+ = u_- se e só se u_+ + ε = u_- + ε. Logo a parte 2D de
  ḡ_ε continua conforme a du_-du_+. Com a = −2e^{−2ζ}/(μ²(u_+ + u_-)²), defina ζ_ε por
  e^{−2ζ_ε} = e^{−2ζ∘ψ_ε}(u_+ + u_-)²/(u_+ + u_- + 2ε)². Defina W_ε por R∘ψ_ε = e^{−ζ_ε}ξ/W_ε e ponha
  φ_ε = φ*∘ψ_ε. W_ε > 0 para ε pequeno em compactos. A regularidade no eixo vem de a métrica ser
  regular, como em R3 de `HREC_IDA.md`.
- **(3)** A identidade δω_{X0} = e^{μτ}(Z + 1)ω* foi conferida nas 7 funções em
  `.codex-runs/2026-10-03-revisao-CRH/c1_lema_c.py`; `test_gauge_action.py` testa a covariância. Pelo
  "Claim" de RT, as equações implicam Ω_i^σ = 0. A derivada em ε = 0 é, por definição,
  DΩ[L_{X0}(ḡ, φ*)] = δω_{X0} = e^{μτ}g* (§2, Passo 8).

∎

**As duas contagens.** Com o Lema G, sobram duas definições razoáveis de "modo físico módulo gauge", e
as duas dão teoremas.
- **Gauge linearizado padrão:** L_X com X campo C¹ em M̂, como no Lema G. X0 é um campo analítico em M̂,
  logo e^{μτ}g* é de gauge, e a contagem é **1**.
- **Gauges restritos aos que fixam o vértice e o cone,** com f(0) = 0: e^{μτ}g* deixa de ser de gauge, e
  a contagem é **2**. Para perturbações gerais, fora do ansatz, a definição consistente é a da revisão
  leve (L4): δg_{++} = 0 no cone, módulo campos X com X^{u_-} = 0 no cone. Mas o modo extra é, pelo Lema C, tangente às translações de T*, que são a mesma
  solução crítica deslocada no tempo.

A hipótese H-cone deixa de ser necessária. A convenção da literatura de fenômenos críticos, em que a
mudança de T* não é instabilidade, coincide com a definição padrão de gauge linearizado.

