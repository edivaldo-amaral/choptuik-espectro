# Revisão independente do manuscrito de T2 — relatório

Revisor: agente Claude (modelo `claude-fable-5-1`), contexto limpo, 08/10/2026.
Objeto: `paper/main.tex`, `paper/macros.tex`, `paper/refs.bib`, `paper/sec/*.tex` (commit 82fc9fc).
Regras seguidas: TAREFA.md desta pasta. Escrita apenas nesta pasta.

## R0. Manifesto do snapshot (início)

`r0_manifesto_inicio.txt`: `sha256sum -c` sobre as 1159 entradas de
`.snapshots/2026-10-08-revisao-artigo/MANIFEST.sha256` — nenhuma linha diferente de `OK` (15:31 de 08/10).


## R1. Auditoria dos números

Artefatos: `r1_numeros.py` (lê só JSON/TXT de `build/` e os dados de RT; parser próprio de `RefA.dat` e
`RefAplusB.dat`; λ e 1/λ em Arb com o K de 80 dígitos de RT), saída `r1_saida.txt` (125 linhas, tabela
"impresso | exato | tipo | seguro?"). Cobertura dos tiles refeita em `r6_cobertura_disco0125.txt` e
`r6_cobertura_disco005.txt` (ver R6).

Resultado: 115 linhas conferem e estão do lado seguro, inclusive **todos** os enclausuramentos das raízes, a tabela
de NK inteira (θ, Z2, Y, r, erro de λ, para as quatro raízes; condições do teorema de NK refeitas em racionais), os
discos (96 + 78, 175 pontos em 22 arquivos, θ máx. 0,998414 e 0,999074, raios 0,00285–0,4995), as janelas
(0,148785 e 0,151664), ε_B, t ≤ 559,151, t·ε_B ≤ 0,03337, testemunhos (0,169758 e 0,354694; c_ref e perdas), a
identificação de μ (dist 1,8155e-5 < r 6,0636e-5; raio de unicidade 3,26e-4), os tiles (154; θ ≤ 0,9961189), o
Rouché de K (raio 0,05, 16 tiles, θ ≤ 0,81256, autovalor 0,40102473111), F1 (recalculado do zero a partir de
`RefA.dat`, idêntico ao JSON), R5 (caixa 2^17×2^20 para S = 200, κ2' = 65/64) e as constantes 8385 e 1440 (exatas).

**λ e 1/λ (Arb, K de RT ± 1e-80, s* = λ0 ± erro exato do C3.json):** λ ∈ [2,6740406, 2,6741147],
1/λ ∈ [0,37395553, 0,37396590]. Os intervalos impressos [2,674040, 2,674115] e [0,373955, 0,373966] os contêm.
Correto.

**Casos em que o número impresso NÃO está do lado seguro ou não confere com o arquivo (10 linhas de `r1_saida.txt`, agrupadas em 8 itens):**

| # | impresso (local) | arquivo | comentário |
|---|---|---|---|
| N1 | ‖I − U*U‖_F ≤ 1,30·10⁻⁶ (counting.tex:251) | `contagem.json` eps_G = 1,30450·10⁻⁶ | arredondado para baixo; deveria ser 1,31·10⁻⁶. Sem efeito em t·ε_B. |
| N2 | \|e_ij\| ≤ 5,5·2,51·10⁻⁸ (counting.tex:128) | `schur.json` norma_phi = [5,50128, 0,78886] | ‖φ_0‖ = 5,50128 > 5,5. Sem efeito (1−10⁻⁶ continua). |
| N3 | "0,2665 in Y+" (setting.tex:309; app_paper_lemmas.tex:72) | `radius_transfer_s3_1.json` strong_defect = 546021/2048500 = 0,266547 | 0,266547 > 0,2665; deveria ser 0,2666. Sem efeito (< 1). |
| N4 | "max Re W(−J) ≤ −0,0232" e "Re s > 2,926" (counting.tex:14,59,74) | `free-part-toro.json` R_livre = −0,0231960 (μ ≥ 1/6) | Com o arquivo, a exclusão é Re s > 2,926002. Com o μ verdadeiro (−μλ0 = −0,023424) vale 2,925774 ≤ 2,926: o impresso é verdadeiro, mas não é o que o certificado grava. |
| N5 | "minimum slack 0,0771 for L and 0,0535 for K" (counting.tex:68) | folga_lambda máx. = 0,077197 (L), 0,053578 (K) | a folga é o **máximo** somado à cota dos centros, não o mínimo; e está arredondada para baixo. Descritivo. |
| N6 | grade 2048×4096 para B e B♯ (counting.tex:64; computation.tex:48) | `symbol-K-1024x2048.json`: nphi 1024, ntheta 2048 | a grade de K é 1024×2048. |
| N7 | Im c1(ω_A,4(·,−1)) = "0,097753…" e "\|c1\| ≥ 0,097753 − 3·10⁻⁸" (physical_modes.tex:94,100) | 0,0977527394… | o dígito impresso está arredondado para cima. A desigualdade do Lema 4.4 continua verdadeira porque \|c1\| = 0,1260, mas o "…" e a cota escritos estão errados (deveria ser 0,097752). Idem estilo: "0,072916…" e "−0,079473…" são arredondamentos, não truncamentos (0,0729159…, −0,0794728…). |
| N8 | (3,9642 + 3,1)·0,0525 < 0,371 (app_paper_lemmas.tex:68) | `l_real_toro.json` nB = 3,964101 | o certificado gravado usou ‖B‖ = 3,9641 + 10⁻⁶, **abaixo** da cota certificada 3,9641333 do símbolo. Com a cota certificada o defeito é 0,370838 < 0,371: a conclusão do artigo vale, mas o JSON não a certifica com o ‖B‖ certificado. |

Observações de fidelidade nos números (não são de arredondamento):
- `β+ ≤ 5191/500` é chamado de "published bound of [RT]" (app_paper_lemmas.tex:71; setting.tex:310 "from the
  published bounds"; physical_modes.tex:230). Não é publicado por RT: é calculado pelo projeto
  (`scripts/independent_component_beta.py`, 10,38108 na norma com pesos η) a partir dos dados e da bola de RT. Além
  disso, esse cálculo **não** aparece no inventário de certificados (`build/reproducao/inventario.json`, passo 4)
  nem no verificador de nível 1. Rodei o script (só leitura, sem `--out`): β_η = 10,3810817 ≤ 5191/500 = 10,382
  (`r1_beta_saida.txt`), e ele confere exatamente com `build/spectrum/component-exterior.json`. O número vale; a
  atribuição a RT e a rastreabilidade não.
- "89 sample points in the A band" (counting.tex:229): pela função `amostras` de `rouche_L_janelas.py`, são 82
  pontos distintos (101 avaliações por região); 89 é o tamanho do cache, que inclui os 10 centros. Cosmético.
- O `Z2` impresso na tabela de NK é o `Z2` do script, que majora ‖A D²F‖ (constante de Lipschitz de A·DF, ver
  docstring e linha `Z2 = 2*v0*max(...)` de `nk_L3_C.py`), enquanto o texto (roots.tex:35) define Z2 ≥ ½‖A D²F‖.
  É conservador (fator 2), não é erro; mas o r impresso é o calculado com esse Z2 maior.

## R2. Teorema contra prova

Numeração conferida em `paper/main.aux` (Lema 2.6 = Lreal, Prop. 2.5 = Fredholm, Lema 4.4 = F1, Lema 4.5 = gauge,
Lema 4.9 = R, Thm 4.11, Prop 4.12, Lema 5.3 = F3, Lema 5.4 = Brauer, …, Prop 6.7 = realidade).

### Teorema Principal (main_theorem.tex:26–59)

| afirmação | onde é provada | cobre exatamente? |
|---|---|---|
| (1) mod i/2·Z, exatamente 4 raízes em Re s > 0, com multiplicidade algébrica, L em Y+ ou no toro | Cor 5.2 ← Thm 5.1 (i)–(iii) + Lema 2.3; Lema 2.6 para Y+ | **Sim.** Raízes em Y+ ⊂ raízes no toro (inclusão Y+ ⊂ T com norma ≤ 1 leva o domínio máximo no domínio máximo); no toro todas as raízes com Re s > 0 têm representante na faixa Q̃ com \|s\| ≤ \|3 + 0,75i\| = 3,092 < 3,1, onde o Lema 2.6 iguala núcleos e cadeias; U = e^{iτ/2} é isomorfismo de Y+ (App. A.1), logo vale para todos os representantes. |
| cada raiz algebricamente simples, com representante real | Prop 6.2 + Lema 6.1 (s_K, s*, s_B+i/2), Prop 6.3 (μ), Prop 6.7 | **Sim** (ver R6). |
| enclausuramentos da tabela | C3.json (R1) | **Sim**, todos do lado seguro (R1). |
| autovetor em μ é g* = (1+Z)ω̄, puro gauge | Prop 6.3 (+ App. A.5, Lema 4.7) | **Sim.** |
| autovetores em s_B, s_K violam as constraints | Prop 6.4 (+ Lema 2.3 para s_B) | **Sim**; em s_B a passagem L_A(s_B+i/2) → L_B(s_B) é U, e U^{-1}C(s)U = C(s+i/2) leva C ≠ 0 em C ≠ 0. |
| autovetor em s* satisfaz as constraints | Cor 6.6 ← Prop 6.5 + Lema 2.4 | **Sim**, mas o texto da Cor 6.6 aplica o Lema 2.4 (enunciado em Y(a,b)) ao autovetor do **toro**; falta citar o Lema 2.6 (\|s*\| ≤ 3,1 ⇒ h ∈ Y+). Ver achado L-2. |
| no eixo imaginário só s ≡ 0, simples, autovetor ∂τω̄ | Thm 5.1(iii) + observação após (17) + Cor 5.2 | **Sim**: L̃_A invertível no segmento {Re s = 0, −1/4 ≤ Im s ≤ 3/4} ⊂ ∂Ω_A ∪ ∂Ω_B e q sem zeros em Re s ≥ 0 dão mult_L(0) = 1 e 0 no resto; mult_L(0) = mult_{L_A}(0) + mult_{L_A}(i/2) = 1 + 0. |
| (2) modos físicos com Re s > 0 mod gauge: dimensão 1, expoente s*, sem modos generalizados | Thm 4.11 + Prop 4.12 + §6.6 | **Sim**, com as hipóteses de Def 4.1/4.2. Em s* a raiz de L (espaço completo X) é simples porque s*+i/2 não é raiz de L_A. |
| cone fixo: quociente de dimensão 2; classe extra com expoente μ, tangente à família de soluções transladadas, isométricas | Prop 4.12 (parte do cone) + Lema 4.7 | **Sim** (ver R4.3). O Lema 4.7 vale "em compactos de M̂" — o enunciado do teorema não diz isso, mas é o sentido natural. |

Hipóteses: o teorema remete à Def 4.1 (C^∞ em τ, analítico em ξ numa vizinhança complexa E de [−1,1]). M̂ vai até
ξ < 41/40, e a Def 4.1 só pede analiticidade perto de [−1,1]; ver a lacuna L-1 em R4. O resto (Def 4.2 (a)/(b),
Def 4.3) é usado exatamente como enunciado.

### Versão informal (introduction.tex:48–60)
Coincide com o Teorema Principal, com duas diferenças de redação: "all real and simple" (o teorema diz "has a
real representative", que é o correto para s_B, cujo representante em L_A é s_B + i/2) e "Choptuik's spacetime" (o
teorema é sobre a solução de RT; o parágrafo anterior diz isso). Coberta pelas mesmas provas.

### Resumo (main.tex:29–44)
Todas as afirmações têm prova no texto: "exactly one unstable mode, modulo gauge" (§6.6 parte 2), "real and simple"
(Props 6.2, 6.7), os intervalos de λ e 1/λ (Obs. 3.1; refeitos em Arb, R1), "four roots with positive real part"
(Cor 5.2). "a rank-two deflation of the gauge roots": as raízes deflacionadas são 0 (translação de fase) e μ; chamar
0 de "gauge root" é coerente com §5.2 ("the two gauge roots 0 and μ") e com a Obs. 3.3, mas o Lema 4.5 só trata
Re s > 0. "The passage from the geometric problem to the operator rests on a global completeness result for the
gauge" — Lema 4.5 (global em M̂, usa o cone dentro do domínio). Nada no resumo excede o que se prova.

## R3. §2 contra RT (fonte `choptuik.tex` de arXiv:1203.3766v1, extraído no scratchpad)

Seções de RT v1 (`\def\thesubsection{\arabic{subsection}}`, linhas do TeX): §1 Introduction (l.119), §2 Einstein
equations coupled to scalar field (l.600), §3 Fourier–Chebyshev series (l.899), §4 Approximate solution yields true
solution (l.1396), §5 Parameter values and computer assisted results (l.2140), §6 Data structures (l.2510).

**Fórmulas** (conferidas linha a linha):
- D^σ = σ∂τ + μ(1+σξ)∂ξ — RT eq. (khdkjhdk…), l.≈840. ✔
- ω (eq. (3) do artigo) = RT (ekjehee), l.≈655, idêntico; ω_i^- = ω_i, ω_i^+ = −Pω_i (2Dpre). ✔
- Forma Ω_i^σ = ξ(D^σ+σμ)(·) + μ(linear) + ξ(quadrático): RT (eezuuirzr). ✔
- F = SΩ^+ = O_μω + μSΓ_1ω + SΞΓ_2(ω,ω): RT define Ω^+ = ΞH_μ v + μΓ_1 v + ΞΓ_2(v,v) e mostra
  bold-J_μ = SΞH_μ; bold-J_μ = diag(K_μ, J_μ, J_μ, J_μ) = ∂τ + μ diag(ξ∂ξ+1, (1+ξ)∂ξ+1, …), pois
  K̸ = ½(J̸ + PJ̸P) = ξ∂ξ + 1 (conta: P(1+ξ)∂ξP = −(1−ξ)∂ξ). Logo D_O do artigo ✔.
- H_μ: linhas (J,0,0,0), (0,J,0,0), (0,−JP,0,0), (0,0,J,0), (0,0,0,J), com J_μ = TauDerivative + μ(…): ∂τ só entra
  como H_1∂τ com H_1 constante de identidades e reflexões ✔; Γ_1, Γ_2, S, S♯, Ξ não têm derivadas ✔ (definições em
  §3). Portanto "∂τ só aparece em O_μ" ✔.
- S = R_ξ∘(½(1+P), comps 2, 4, 5); S♯ = (½(1−P) na comp. 1, −P na comp. 2♯) ✔ (App. A.4: c = (Oq_1, −Pq_2♯) ✔).
- Normalização Re(ω̄_4)_{10} = 0 (espaço "Gauged") ✔; pesos (65/64)^{|m|}(5/4)^{|n|} ✔; norma de RT usa
  ‖z‖_C = |Re z| + |Im z| ≥ |z| e soma das componentes sem pesos η (o artigo usa |·| e pesos η: normas equivalentes,
  e a desigualdade usada em F1 continua válida porque ‖·‖_C ≥ |·|).
- Identidade ♯: RT (6383hi3h39z), §3: O♯_μω♯ + μR_ξΓ♯_1ω♯ + Γ♯_2(ω,ω♯) = 0 com Γ♯_1w = (0, w_1), sob
  SΩ^+(μ,ω) = 0 ✔ = eq. (9). A forma off-shell (App. A.4) é a identidade geral de RT §2.3 (l.≈790), com
  Ω_i = Ω_i^- ✔ (o artigo usa e = −Pq = Ω^-). Conferi que T_s do App. A.4 reproduz os termos de RT e que
  R_ξT_sIc dá a parte principal diag(∂τ+s+μ(ξ∂ξ+1), J_s) ✔ (conta: ½(1−P)J_s c_1 = (∂τ+s+μ(ξ∂ξ+1))c_1 para c_1 ímpar).
- Cotas (7): |μ − (μ_A+μ_B)| ≤ 2^{-277} e ‖ω̄ − (ω_A+ω_B)‖ ≤ 2^{-277}: RT (kdkkkskkskksksks), **§5**, l.≈2353.
  ‖ω_B‖ ≤ 2^{-25}: hipótese (dkjfhk2) de **§4** com 𝓛_2 = 2^{-25} fixado e verificado em **§5**. ✔ (valores).
  Também |μ_B| ≤ 43·2^{-40} (§5), donde |μ − μ_A| ≤ 3,9·10⁻¹¹ (conferido: 3,881·10⁻¹¹, `r1_saida.txt`).

**Citações de seção:**

| local | citação | confere? |
|---|---|---|
| introduction.tex:35 | [§1] citação literal | ✔ (§1; texto idêntico até "inconclusive", ver R8) |
| setting.tex:17 | [Main result] | ✔ (§1) |
| setting.tex:71 | [§2] ω e codificação ± | ✔ |
| setting.tex:113 | [§3] fórmulas explícitas | ✔ |
| setting.tex:119 | [§3–4] construção por contração, pesos, perto dos dados; e logo após "Their final bounds are (7)" | **parcial**: a contração é §4 e o cenário de séries §3, mas os pesos κ, os dados e as cotas finais 2^{-277}, 2^{-25} estão em **§5** (e na visão geral de §1). Deveria ser "§§4–5". |
| setting.tex:159 | [§§2–3] H_μ | ✔ (§3) |
| setting.tex:236 | [§3] sistema ♯ | ✔ |
| app_formulas.tex:11 | [§2] Ω_1^+, Ω_2♯^+ | ✔ |
| app_paper_lemmas.tex:28 | [§3] inversas explícitas de J | ✔ (eqs. dkhddjjjjjff, §3; as cotas estão em §4) |
| app_paper_lemmas.tex:45 | [§4] R_ξ limitado para b > 1 | ✔ (§4, estimativas de U, V; ‖R_ξ‖ ≤ 2κ_2/(κ_2²−1) = D(κ_2) do App. A.3) |
| app_paper_lemmas.tex:91 | [§2] identidade off-shell | ✔ (§2.3) |
| physical_modes.tex:304 | [§2] fórmulas de curvatura | ✔ (dfkhdjhsdshkfd) |
| app_paper_lemmas.tex:71; setting.tex:310; physical_modes.tex:230 | "published bound(s) of [RT]" para β+ | **não**: 5191/500 não é publicado por RT (ver R1) |

**Lema 2.2 refeito.** F(μ,ω) = O_μω + μSΓ_1ω + SΞΓ_2(ω,ω), com O_μ = ∂τ + μD_O e Γ_1, Γ_2, S, Ξ algébricos em τ
(reflexões, multiplicação por ξ, R_ξ e produtos pontuais). Para h 4π-periódico,
D_ωF(μ,ω̄)[e^{sτ}h] = O_μ(e^{sτ}h) + μSΓ_1(e^{sτ}h) + 2SΞΓ_2(ω̄, e^{sτ}h)
= e^{sτ}[(∂τ + s)h + μD_O h + μSΓ_1h + 2SΞΓ_2(ω̄,h)], porque multiplicar por e^{sτ} comuta com tudo o que age só em
ξ e com produtos pontuais, e e^{-sτ}∂τe^{sτ} = ∂τ + s. Logo L(s) = J + B + s com J = O_μ, B = μSΓ_1 + 2SΞΓ_2(ω̄,·).
Para C: F♯ = S♯Ω^+ com Ω^+ = ΞH_μω + μΓ_1ω + ΞΓ_2(ω,ω) e H_μ = H_1∂τ + (termos em ξ), donde
C(s) = C_0 + sC_1, C_1 = S♯ΞH_1. ✔

**C_1 = (0, −ξh_2).** H_1h = (h_1, h_2, −Ph_2, h_3, h_4); ΞH_1h = (ξh_1, ξh_2, −ξPh_2, ξh_3, ξh_4);
S♯ toma (½(1−P)(ξh_1), −P(−ξPh_2)). Como h_1 é ímpar, ξh_1 é par e ½(1−P)(ξh_1) = 0; e P(ξPh_2) = (Pξ)(PPh_2) = −ξh_2.
Logo C_1h = (0, −ξh_2) ✔. Pela fórmula de RT, Ω_2♯^+ = ξ(D^++μ)ω_2^+ + …, o termo em s de δΩ_2♯^+ é ξh_2^+ = −ξPh_2, e
−P(−ξPh_2) = −ξh_2, o mesmo ✔. Também conferi δΩ_1^+ e δΩ_2♯^+ do App. B contra (eezuuirzr) com σ = + (o
polinômio 𝒬_1 bate termo a termo) e M_1 = C_1 (App. B) com a relação da §4.5.

## R4. §4, partes reescritas

### R4.1 Lema 4.5 (completude do gauge) — prova refeita; correta, com uma hipótese implícita

Passo 0 (média). Para R ∈ SO(3): R*(L_Xḡ) = L_{R*X}(R*ḡ) = L_{R*X}ḡ. Se L_X(ḡ,φ̄) é invariante, L_{R*X}(ḡ,φ̄) = L_X(ḡ,φ̄)
para todo R e, integrando (X é C¹, derivação e média comutam), L_{X̄}(ḡ,φ̄) = L_X(ḡ,φ̄). Um campo SO(3)-invariante em
R × R³ tem parte espacial b(τ,ξ)x (o estabilizador SO(2) de x ≠ 0 só fixa múltiplos de x), logo X̄ = X^τ∂τ + X^ξ∂ξ e,
por continuidade em x = 0, X^ξ → 0 no eixo. **Hipótese implícita:** a invariância de L_X(ḡ,φ̄) não está no enunciado.
Ela segue de "em forma RT" só se "forma RT" pressupõe simetria esférica (a definição na §4.1 é dada para
perturbações, que são esfericamente simétricas pela Def 4.1). Nos usos ((R6), Lema 4.9(3), Prop 4.12) L_X = δ' − δ é
esfericamente simétrico, então a prova vale onde é usada. Correção: acrescentar "L_X(ḡ,φ̄) esfericamente simétrico"
ao enunciado (achado I-8).

Passo 1. (L_Xḡ)_{++} = X^c∂_cḡ_{++} + 2ḡ_{c+}∂_+X^c = 2ḡ_{+−}∂_{u+}X^{u−} (ḡ_{++} ≡ 0, ḡ_{A+} = 0) ✔. Os conjuntos de nível
{u_- = v} ∩ M = {u_+ < min(0, v, −81v)} e {u_+ = w} ∩ M = {w < u_- < −w/81} são intervalos ✔; f_- ∈ C¹(R) porque
u_- assume todos os valores reais em M, inclusive 0 ✔.
Passo 2. X̄^ξ = 0 no eixo ⇒ X̄(u_+ − u_-) = X̄(−2ξe^{−μτ}) = 0 em ξ = 0 ⇒ f_+ = f_- em (−∞,0) ✔.
Passo 3. ψ = Θ², u∘ψ = Λu: du∘dψ = Λdu, logo (ψ*X)^u = Λ^{-1}X^u∘ψ = Λ^{-1}f(Λu) ✔; δω_X∘Θ² = δω_{(Θ²)*X} usa
(Θ²)*(L_Xḡ) = e^{−4K}L_{(Θ²)*X}ḡ e a invariância de ω por ζ → ζ + const (dita explicitamente em RT §2.1) ✔.
Passo 4. D^+ = 2μe^{−μτ}∂_{u−}, D^- = −2μe^{−μτ}∂_{u+} (conferido: D^+u_- = 2μe^{−μτ}, D^+u_+ = 0, D^-u_+ = −2μe^{−μτ},
D^-u_- = 0), donde X(φ̄) = (e^{μτ}/2μ)(f(u_-)ω̄_4^+ − f(u_+)ω̄_4^-) ✔; no eixo ω̄_4^±(τ,0) = ∓ω̄_4(τ,0) ✔ e
X(φ̄) = −μ^{-1}e^{μτ}f(−e^{−μτ})ω̄_4(τ,0) = c_2 ✔. c_2 = 0 (função real antiperiódica contínua tem zero) ✔; F1 + zeros
isolados ⇒ f = 0 em (−∞,0] ✔; fora do cone o conjunto {u_- ∈ I, u_+ < −81·sup I} é aberto não vazio e
analiticidade de ω̄_4^+ = −Pω̄_4 no cilindro conexo leva a ω̄_4^+ ≡ 0, contra F1 ✔. D^± são linearmente independentes
em todo ponto (det = 2μ), então D^±(Xφ̄) = 0 dá Xφ̄ constante ✔.
Passo 5. Linearidade de X ↦ δω_X + Passo 4 aplicado a (Θ²)*X̄ − e^{4πs}X̄ (do mesmo tipo) ⇒ f(Λu) = ρf(u),
ρ = Λe^{4πs} = e^{4π(s−μ)} ✔.
Passo 6. Re s > 0 ⇔ \|ρ\| > Λ; os três casos (\|ρ\| > 1 ⇔ Re s > μ; = 1; Λ < \|ρ\| < 1) esgotam Re s > 0 ✔. No
terceiro, f(0) = 0 e f'(0) = lim(ρ/Λ)^j f(u)/u com \|ρ/Λ\| > 1 força f(u) = 0 ✔. Conclusão ρ = 1, f ≡ c, X̄ = cX_0 ✔, e
X − cX_0 é Killing com (X − cX_0)(φ̄) = 0 ✔. ("Equivalently" é exagero: "a média de X é cX_0" implica as outras duas
afirmações, mas a recíproca só vale para a média; cosmético.)
Variante (b). Passos 1–5 só usam continuidade de f em 0 ✔, e o Passo 6 para Re s ≥ μ também ✔. Para 0 < Re s < μ:
δφ = X(φ̄) contínuo com D^σδφ analíticos ⇒ a 1-forma é fechada e o primitivo analítico difere de δφ por constantes
nos dois lados do cone, iguais por continuidade ⇒ δφ analítico ✔. No eixo, f(−e^{−μτ}) = −μe^{−μτ}δφ(τ,0)/ω̄_4(τ,0) é
quociente de analíticas, contínuo ⇒ singularidades removíveis ⇒ f analítica em (−∞,0) ✔. Perto de (τ0,1) com
ω̄_4^+(τ0,1) = −ω̄_4(τ0,−1) ≠ 0 (F1'), f(u_-) é analítica e u_- é submersão ⇒ f analítica em 0 ✔. Taylor:
a_k(Λ^k − ρ) = 0; ρ = Λ^k ⇔ s ≡ μ(1−k) (mod i/2·Z), Re s = μ(1−k) > 0 ⇒ k = 0 ⇒ Re s = μ, impossível em 0 < Re s < μ;
logo f ≡ 0 perto de 0 e, por f(u) = ρ^{-j}f(Λ^ju), em toda parte ✔. Observação 4.6: f = (−u)^{1−s/μ} satisfaz
f(Λu) = Λ^{1−s/μ}f(u) = e^{4π(s−μ)}f(u) ✔, e não é admissível em (b) porque δφ não seria analítico através do cone ✔.

### R4.2 Passo R7 do Lema 4.9 — fórmula e restrições de jato: corretas

Derivação: se δω' = e^{sτ}Σ_{j<k}(τ^j/j!)h_{k−1−j}, então, como ∂τ entra no operador selecionado só com
coeficiente 1 e nas constraints só com C_1 (Lema 2.2), o operador linearizado manda e^{sτ}τ^jh em
e^{sτ}(τ^jL(s)h + jτ^{j−1}h). Somando e reindexando, o coeficiente de τ^j/j! é L(s)h_{k−1−j} + h_{k−2−j}; com
i = k−1−j: L(s)h_i = −h_{i−1} (h_{−1} = 0) ✔. Nas constraints: C(s)h_i + C_1h_{i−1} = 0 ✔ = (16). Para k = 2 dá
e^{sτ}(h_1 + τh_0) ✔. Com D = C_0 − C_1L(0): Dh_i = C_0h_i − C_1(−s_0h_i − h_{i−1}) = C(s_0)h_i + C_1h_{i−1} ✔.

Teste numérico (`r4_jato.py`, saída `r4_jato_saida.txt`): família afim de matrizes com L(s) = A + s (A com bloco de
Jordan de tamanho 3 em −s_0 e mais três autovalores), K(s) = K_0 + s, C(s) = C_0 + sC_1 e M(s) = M_0 + sM_1 construídos
para que K(s)C(s) = M(s)L(s) (resíduo ≤ 2,8·10⁻¹⁴). Para cadeias de comprimento 1, 2, 3 e quatro escolhas de D:
u(τ) da fórmula do artigo resolve u' + Au = 0 (≤ 1,1·10⁻¹⁵); a constraint c(τ) = C_1u' + C_0u coincide com
e^{s_0τ}Σ(τ^j/j!)[C(s_0)h_{k−1−j} + C_1h_{k−2−j}] (≤ 4,9·10⁻¹⁵); Dh_i = jato_i (≤ 2,9·10⁻¹⁵); e a dimensão física
medida por amostragem de c(τ) é igual a dim E − posto(D|E) nos quatro casos (1, 2, 3, 0). Nos dois exemplos 2×2, a
contagem "ingênua" (C(s_0)h_j = 0 para cada j) dá 2 contra 1 verdadeiro e 1 contra 2 verdadeiros: a observação após
o Thm 4.11 está certa nos dois sentidos.

### R4.3 Proposição 4.12 (multiplicidade física), inclusive o cone fixo — correta

Fórmula: o espaço de modos generalizados em forma RT com expoente s_0 é {e^{−Aτ}v : v ∈ E} (u(τ) finito em τ porque
A + s_0 é nilpotente em E); por K(0)D = DA (de KC = ML, comparando potências de s: M_1 = C_1, M_0 = C_0 + K(0)C_1 − C_1L(0),
K(0)C_0 = M_0L(0) ⇒ K(0)D = DA ✔), D(E) é K(0)-invariante e de dimensão finita e De^{−Aτ}v = e^{−K(0)τ}Dv, que é ≡ 0 sse
Dv = 0. Logo modos físicos ≅ ker(D\|E) e dim = dim E − posto(D\|E) ✔; o quociente por G_{s_0} ∩ E exige
G_{s_0} ∩ E ⊂ ker(D\|E), que vale porque C(μ)g* = 0 (App. A.5) — o texto não diz isso explicitamente (cosmético).

Cone fixo. (i) A fixação fica tangente ao cone? Sim: δg_{++} = 0 em {u_- = 0} ⇒ f_-(·,1) = 0; em ξ = 1 a equação de
y_- é (∂τ + s − μ)y_-(·,1) = 0 ⇒ y_-(·,1) = 0 fora de ressonância; na ressonância o múltiplo livre de e^{im_0τ/2} é posto
em zero ⇒ X^{u−} = e^{(s−μ)τ}y_- = 0 no cone ✔. (ii) A obstrução some? Sim: κ = (f_-(·,1))_{m_0} = 0 ✔. (iii) Vale com
G_s = 0? Sim: duas fixações que preservam o cone diferem por X tangente com L_X em forma RT; pelo Lema 4.5, X̄ = cX_0 e
X̄(u_-) = c no cone (a média de um campo tangente ao cone é tangente) ⇒ c = 0 ⇒ L_X = 0; logo a fixação é única. Se
h = cg* então δ' = cL_{X_0} e δ = cL_{X_0} − L_X, que só é puro gauge (no sentido do cone) se c = 0, porque L_{X_0} = L_Y com Y
tangente daria δω_{X_0−Ȳ} = 0 ⇒ Ȳ = X_0 (Passo 4), contradição. Sobrejetividade: os modos do Lema 4.10 estão em
forma RT, logo δg_{++} ≡ 0 ✔. Também conferi que somar L_X com X tangente preserva δg_{++} = 0 no cone
((L_Xḡ)_{++} = 2ḡ_{+−}∂_{u+}X^{u−} e ∂_{u+} é tangente ao cone) ✔.

### R4.4 Teorema 4.11 (correspondência) — correto, com uma lacuna de domínio

- Bem definido no quociente: duas fixações diferem por campo admissível cujo δω ∈ G_s (Lema 4.5) ✔; δ ~ δ + L_Y
  (Y de classe (a) ou (b)) dá X + Y de classe (b) com L_{X+Y} em forma RT e RT-variáveis analíticas ✔.
- Linear ✔ (fórmula (15) linear em f; na ressonância, escolha canônica do múltiplo de e^{im_0τ/2}).
- Injetivo: h ∈ G_s ⇒ δ puro gauge, por (R6) ✔ (δω' = 0 ⇒ δQ = 0, δζ e δφ constantes de tipo Floquet com Re s > 0 ⇒ 0).
- Sobrejetivo: Lema 4.10 ✔ (conferi a, b: (∂τ+s)z = ½[(1−ξ)h^+ − (1+ξ)h^-], ∂ξz = (h^+ + h^-)/(2μ); paridades: a par, b
  ímpar ⇒ z par ✔). G_s ⊂ ker L(s) ∩ ker C(s) por (8) ✔.
- Transporte (15): verifiquei T_c(s) y = f pela identidade (∂τ + μ(ξ−c)∂ξ)(g∘Φ_t) = −d/dt(g∘Φ_t) com
  Φ_t = (τ−t, c + e^{−μt}(ξ−c)); o integrando é O(e^{−Re s·t}) porque f − f(·,c) se anula em ξ = c ✔; unicidade pela
  expansão em (ξ−c)^n com (∂τ + s + μ(n−1))y_n = 0 ✔ (n = 0 exige s ∉ μ − (i/2)Z).
- (R2): o termo secular κτe^{im_0τ/2} resolve a parte ressonante ✔; L_{τX_0} = τL_{X_0} + 2dτ⊙X_0^♭ ✔;
  δω[A + τB] = (δω[A] + termos sem derivada de B) + τδω[B], com δω[B] = κe^{μτ}g* = e^{s_0τ}κe^{im_0τ/2}g* ✔;
  e^{−im_0τ/2}L(s_0)e^{im_0τ/2} = L(μ) ✔ ⇒ cadeia de comprimento 2 de L_A em μ, contra a Prop 6.3 ✔.
- (R3): refiz δζ = −½μ²e^{2ζ}(γξ² + d) (coeficiente de dξ² de ḡ_2D = e^{−2ζ}(−(1−ξ²)dτ² − (2ξ/μ)dτdξ + dξ²/μ²)) e
  δQ = −½We^{2ζ}(d(2μQ + ξ²Q²) − μ²γ) (de δW = −W(δζ + δb/2b), b = e^{−2ζ}ξ²/W², W² − μ² = ξ²(2μQ + ξ²Q²)) ✔.

**Lacuna L-1 (reparável).** A Def 4.1(iii) só pede analiticidade numa vizinhança complexa E de [−1,1], mas os modos
vivem em M̂ = R × B_{41/40}. As passagens (R3)–(R6) identificam δω' em **todo** M̂ com a série de h (coeficientes de
Chebyshev, que só "veem" [−1,1]). Para ξ ∈ (1+ε, 41/40) isso precisa de um argumento que o texto não dá: além do cone
as duas direções características (dξ/dτ = μ(1+ξ) e μ(ξ−1)) são positivas, as superfícies ξ = const são espaciais, e a
solução do sistema RT linearizado em ξ > 1+ε é determinada pelos dados em ξ = 1+ε (domínio de dependência). Com isso
(R6) ("δω' = 0 ⇒ δ' = 0 em M̂") e a injetividade valem. Correção: uma frase com esse argumento, ou pedir E ⊃ [−41/40, 41/40].
A construção de X em (R1) não tem esse problema: as características de T_{±1} ficam em (−1, 41/40).

## R5. §5

**Rouché de Gohberg–Sigal na forma de homotopia (counting.tex:84–88).** A forma de GS 1971 é: se A é holomorfa,
Fredholm (finitamente meromorfa) em Ω̄, invertível em ∂Ω, e ‖A^{-1}B‖ < 1 em ∂Ω, então A e A+B têm o mesmo número
de zeros. A forma de homotopia segue por compacidade de ∂Ω × [0,1]: se F_t = R + t(H−R) é invertível em ∂Ω para todo t,
(t,s) ↦ ‖F_t(s)^{-1}‖ é contínua e limitada, e passos t → t' com \|t'−t\|‖F_t^{-1}(H−R)‖ < 1 aplicam GS sucessivamente.
**Mas** o passo a passo usa GS em cada F_t, e GS pede F_t holomorfa e Fredholm em Ω̄; o enunciado só supõe isso para H e R. Não encontrei contraexemplo (para somas diretas de blocos finitos o princípio do mínimo do módulo aplicado ao determinante impede a perda de Fredholmness no interior), mas também não é uma consequência imediata, e o enunciado clássico supõe a família Fredholm. Na aplicação não há problema: H − R = (blocos fora da diagonal) + P_far(B+E)Q_0P_far é compacto, então F_t = R + t·compacto é Fredholm de índice 0 para todo t. Correção: acrescentar "com H − R compacto" (ou "R + t(H−R) Fredholm de índice 0 em Ω̄") ao enunciado (achado I-6). O critério suficiente ‖R^{-1}(H−R)‖ < 1 está certo.

**Lema 5.4 (Brauer) — correto.** L(s)^{-1}x_0 = x_0/s e L(s)^{-1}x_1 = x_1/(s−μ) ⇒ L^{-1}(L̃−L) tem posto 2 com imagem em
span{x_0,x_1}; det(I + MD^{-1}) = det(D + M)/det D com D = diag(s, s−μ), M_ij = δ_jφ_i^*x_j ⇒ q(s)/(s(s−μ)) ✔. Com
φ_i^*x_j = δ_ij + e_ij e δ_0 = 1: q = (s + 1 + e_00)(s − μ + δ_1(1+e_11)) − δ_1e_01e_10 ✔ (o texto escreve δ_0δ_1, igual).
Multiplicidade: fatoração L̃ = L(I + L^{-1}(L̃−L)), multiplicatividade do resíduo logarítmico de GS para funções
finitamente meromorfas e ord det para perturbação de posto finito da identidade ✔. Teste em matrizes (`r5_brauer.py`,
`r5_brauer_saida.txt`): identidade do determinante ≤ 1,9·10⁻¹² e raízes(L̃) = raízes(L) − {0, μ} + zeros(q) como
multiconjuntos, inclusive com 0 de multiplicidade 2 (bloco de Jordan). Os fatores de q: \|s + 1 + e_00\| ≥ 1 − \|e_00\| e
\|s − μ + δ_1(1+e_11)\| ≥ 1 + μ_A − μ − δ_1\|e_11\| em Re s ≥ 0, ambos ≥ 1 − 2·10⁻⁷ com ‖φ_0‖ = 5,50128 (não 5,5; R1 N2)
e ‖x_1 − x_1^A‖ ≤ 2,5056·10⁻⁸ ✔. **Nenhum script** do verificador confere q, apesar de §7.1(a) dizer que "the
deflation polynomial q" está entre as desigualdades finais feitas em racionais (achado E-2). A conta é de mão e está
certa.

**Lema 5.5 (bloco finito) — correto como enunciado, mas não é o que o certificado usa.** J + s preserva m e é
triangular superior em n (RT: J_μ = V_{1,−1}(U_{0,·} + U_{1,·}), e U[k,f] tem (U v)_n = f v_{n+k}); logo preserva o
subespaço {n < 128}, e Q_0(s)P_Z = P_ZQ_0P_Z = (P_Z(J+s)P_Z)^{-1} ✔; dimensão 14016 = 31·64 (h_1: m par, n ímpar) +
2·31·128 (h_2, h_3) + 32·128 (h_4, m ímpar) ✔. Porém o certificado não usa P_ZL̃P_Z do operador verdadeiro: usa a matriz
de ponto flutuante A = W L̃^A_Z(0) W^{-1} montada com ω_A e x_j^A, "tomada exata como dado" (docstrings de
`rouche_L_contagem.py` e `rouche_L_rig.py`), e a diferença para o operador verdadeiro vai para E. A lógica fecha (a
referência pode ser qualquer família analítica cujos zeros se contem), mas o texto da §5.3 e o Lema 5.5 dizem que a
referência é a truncação verdadeira e que os zeros são autovalores de L̃_Z(0). Correção: definir R com o bloco
(𝖠 + s)(P_Z(J+s)P_Z)^{-1}, 𝖠 a matriz dos dados, e dizer que a diferença entra nas entradas (achado I-5). Mesmo item:
na §5.3 "t_p ≥ ‖Ĥ_ZZ(p)^{-1}‖" está errado — o código certifica t_p ≥ ‖(𝖠 + p)^{-1}‖ (`rouche_L_rig.py`, docstring:
"t_p >= ||(A + p)^-1||, por Loewner"), que é o que a §5.5 usa em t(s) ≤ t_p/(1 − rt_p) (1-Lipschitz de σ_min(𝖠+s)).

**Lema 5.6 (Perron) — correto.** Bloco (i,j) de (I − VR) − tVE: ‖·‖ ≤ ρ_i + t n_ii ≤ ρ_i + n_ii (i = j), ≤ t n_ij ≤ n_ij
(i ≠ j): é aqui que entra a monotonicidade em t ∈ [0,1]. Com w > 0 e Nw ≤ (θ+ε)w (vetor de Perron de N + ε𝟙, positivo
mesmo se N for redutível), na norma max_i ‖x_i‖/w_i a norma é ≤ θ + ε < 1 ⇒ V(R+tE) invertível ⇒ R+tE injetivo; Fredholm
de índice 0 ⇒ invertível ✔. Nos certificados o vetor v é dado e Nv ≤ θv é verificado em racionais, o que dispensa o
argumento de Perron–Frobenius ✔ (conferido para os 16 tiles do Rouché de K, R6).

**Lema 5.7 (janelas) — correto.** F_k(s) = I − V_J(c_k)Ĥ_JJ(s) é analítica em Re s > −μ; ‖F_k‖ é subarmônica ✔; o bordo
de cada retângulo é coberto por discos (q, r) com q no ponto médio de cada subsegmento e r = meio comprimento·(1+10⁻⁹)
(`amostras` em `rouche_L_janelas.py`) ✔; os 10 retângulos de cada faixa cobrem [0,3]×[−1/4,1/4] e [0,3]×[1/4,3/4]
(conferido nas listas `REGIOES`, `REGIOES_B`) ✔; Ĥ_JJ é matriz finita, então ‖F_k‖ < 1 ⇒ invertível ✔. Na faixa B o
centro do retângulo [0,8;1,2]×[1/4,3/4] é c = 1 (fora do retângulo); o princípio do máximo não exige c dentro ✔.

**Lema 5.8 (contagem finita) — correto.** ‖𝖠 − 𝖡‖ ≤ ‖𝖠U − UT‖‖U^{-1}‖ e σ_min(U)² ≥ 1 − ‖I − U*U‖_2 ≥ 1 − ‖I − U*U‖_F ✔;
t·ε_B < 1 ⇒ 𝖠 + s + t'(𝖡−𝖠) invertível em ∂Ω ⇒ det não se anula ao longo da homotopia e o índice de giro é constante ✔.

**O Teorema 5.1 (ii)–(iii) segue?** Sim, dos lemas e certificados descritos, com as correções de redação acima:
N(L̃_A,Ω) = N(R,Ω) (Perron nos 174 discos + cobertura racional), N(R,Ω) = N(Ĥ_ZZ,Ω) (janelas sem zeros em Ω̄), N(Ĥ_ZZ,Ω_A) = 2 e
N(Ĥ_ZZ,Ω_B) = 1 (Lema 5.8: t·ε_B ≤ 0,0334), e (17) dá 3 e 1. (iii) segue da invertibilidade de L̃_A em ∂Ω e de q sem zeros.

## R6. §6

**Lema 6.1 (bordejado ⇔ simples) — correto.** DF(y,λ)[η,ν] = (L(λ)Q_0η + νQ_0y, ℓ(Q_0η)). (⇒) se DF injetivo: w ∈ ker L,
w' = w − ℓ(w)h_* ∈ ker com ℓ(w') = 0 ⇒ DF[Q_0^{-1}w', 0] = 0 ⇒ w' = 0; se L k = −h_*, k' = k − ℓ(k)h_* dá
DF[Q_0^{-1}k', 1] = 0, contradição. (⇐) como no texto. DF = [[L(λ)Q_0, Q_0y],[ℓQ_0, 0]] é perturbação de posto finito de
diag(L(λ)Q_0, 0), que tem índice 0 ✔.

**Teorema de NK (roots.tex:31–40) — correto**, com duas observações: (a) a conclusão "F tem zero único" a partir de
"T tem ponto fixo único" exige 𝒜 injetiva; aqui 𝒜 = diag(inversa exata da matriz bordejada, V_J, I) é invertível, mas o
texto não diz (cosmético); (b) com Z_2 ≥ ½‖𝒜D²F‖: ‖T(x) − x_0‖ ≤ Y + θr + Z_2r² ≤ r no menor raiz r de
Z_2r² − (1−θ)r + Y, e ‖DT(x)‖ ≤ θ + 2Z_2r = 1 − √((1−θ)² − 4Z_2Y) < 1 ✔; unicidade em B(x_0,ρ) com θ + 2Z_2ρ < 1 ✔. O
script (`nk_L3_C.py`) usa um Z_2 que majora ‖𝒜D²F‖ inteiro (fator 2 a mais) e verifica Y + θr + Z_2r² ≤ r e
θ + 2Z_2r < 1 em racionais; conservador. As quatro condições foram refeitas em `r1_numeros.py` ✔.

**Prop 6.3 (μ identificado) — correta.** ‖x_g − x_0‖ ≤ 1,8155·10⁻⁵ < r = 6,0636·10⁻⁵ (E.json), e o raio de unicidade
(1−θ)/(2Z_2) = 3,26·10⁻⁴ dá folga ainda maior. Redação: "‖x_g − x_0‖ ≤ 1,82·10⁻⁵ < r, with r ≤ 6,07·10⁻⁵" usa uma cota
**superior** de r para concluir dist < r, o que não segue dos números impressos; trocar por "r ≥ 6,06·10⁻⁵" ou citar o
raio de unicidade (cosmético C-4).

**Prop 6.5 — a região certificada cobre o que a Cor 6.6 usa? Sim, mas não pelo que o verificador confere.**
- Os tiles: os JSON `build/cobertura_lab/grande-*.json` gravam `disco = [0.401, 0.05]` (o planejador excluiu o disco de
  raio 0,05). Mas `scripts/check_cover.py` tem `--disco 0.401,0.125` por padrão, e o verificador de nível 1 o chama sem
  `--disco` (`verifica_T2.py`, `check_s3a`). Também a revisão independente de 24/09 (`2026-09-24-revisao-agente`) declara
  ter conferido a injetividade "fora de |s − 0,401| < 1/8". Ou seja, a coroa 0,05 ≤ \|s − 0,401\| < 1/8, que o Rouché de K
  (círculo de raio 0,05) **não** cobre, só estava coberta por construção do planejador (`tile_cover.py`, cujo padrão também é 1/8, foi rodado com `--disco 0.401,0.05`); nenhuma checagem a conferiu.
- Refiz a cobertura com o raio do artigo: `check_cover.py build/cobertura_lab build/cobertura_local --disco 0.401,0.05` →
  "154 discos certificados; células resolvidas: 2871; sem cobertura: 0; COBERTURA COMPLETA" (`r6_cobertura_disco005.txt`;
  com o padrão 1/8: 1770 células, `r6_cobertura_disco0125.txt`). Há tiles com centro a 0,038 de 0,401 e raio 0,0147. Logo
  a Prop 6.5 **vale como enunciada**; o defeito é do verificador e da descrição do que foi checado (achado E-1).
- O Rouché de K: 16 tiles de raio r_Z = 0,0098017238 centrados no círculo \|s − 0,401\| = 0,05; o pior ponto do círculo
  fica a 2·0,05·sin(π/32) = 0,0098017140 do centro mais próximo — cobertura com folga de 9,8·10⁻⁹ (`r5_brauer.py`) ✔.
  Recalculei em racionais o Perron das 16 matrizes `N_rouche` com vetor de Perron próprio: máx. 0,8125577 ≤ 0,8126 ✔.
  O nível 1 (`check_rouche_K`) só lê as marcas `certificado`/`ok`, não refaz essa conta (achado E-2).
- Semiplanos: tiles e Re s > 1,765 cobrem 0 ≤ Im s ≤ 1/4; o círculo é coberto inteiro. A Cor 6.6 usa s* (real, a 0,332
  de 0,401: fora do disco de 0,05, dentro dos tiles) e s_K (real, \|s_K − 0,401\| ≤ 1,72·10⁻⁴: dentro do disco) ✔. O
  enunciado "exceto em um ponto s_K'" é frouxo se s_K' tivesse Im < 0, mas a Cor 6.6 resolve (s_K é raiz real de K no
  disco, e o disco tem uma só) ✔.
- Domínio: Cor 6.6 diz que h está no domínio de L no **toro** e invoca o Lema 2.4, que é enunciado em Y(a,b). O caminho
  correto é: \|s*\| ≤ 3,1 ⇒ h ∈ Y+ (Lema 2.6) ⇒ C(s*)h ∈ domínio de K em Y♯(129/128, 9/8) ⊂ T♯ (Lema 2.4 + App. A.4) (achado L-2).

**Prop 6.7 (realidade) — correta.** 𝒞: h_mn ↦ conj(h_{−m,−n}) (= conj(h_{−m,n}) do App. A.1, pois h_{m,−n} = h_mn) é
isometria em Y(a,b) (pesos pares em m e n) e **não** é sequer limitada no toro (‖𝒞h‖² = Σa^{−2m}\|h_m\|²), por isso o
argumento passa por Y+ ✔. 𝒞L(s)𝒞 = L(s̄) porque ω̄ é real e (𝒞∂τ𝒞h)_m = (im/2)h_m ✔; 𝒞 preserva setores (paridade de m)
✔. Lema 2.6 transfere raízes e multiplicidades para \|s\| ≤ 3,1, e s̄_K, s̄*, s̄_B + i estão nesse disco ✔. Em Ω_A, os três
discos NK são disjuntos e contêm uma raiz simples cada; como N = 3, não há outra raiz em Ω_A, logo s̄_K = s_K ✔. Em Ω_B
o mapa s ↦ s̄ + i preserva Ω_B (Im ↦ 1 − Im) e a única raiz é fixa ⇒ Im = 1/2 ✔.

**Montagem (§6.6) — correta.** Parte (1): Cor 5.2 + Props 6.2, 6.3, 6.7; parte (2): Prop 4.12 forma simples em cada raiz
(espaço de raízes de L em X de dimensão 1 em cada uma, porque s + i/2 não é raiz de L_A para s ∈ {μ, s_K, s*} e s_B não é
raiz de L_A); Thm 4.11. O intervalo de s* "λ_0 ± 1,0163·10⁻⁵, rounded outward to eight digits" ✔ (R1).

## R7. §7 contra os arquivos

Artefatos: `r7_digests.py` / `r7_digests_saida.txt` (cópia do `digest_certificados.py` que escreve só aqui, mais
sha256 recalculado de cada arquivo do inventário), `r1_saida.txt`, leitura de `build/reproducao/*`.

| afirmação (local) | arquivo | resultado |
|---|---|---|
| 175 pontos em 22 arquivos; 174 discos de raio positivo e 1 de raio 0 (computation.tex:72) | `discos_v2/*`, `faixaB/discos/*` | ✔ (96 + 79 pontos; o de raio 0 é `faixaB/discos/disco_0_75.json`, p = 0,291688 + 0,75i) |
| pior Perron 0,999074, diferenças ≤ 7,8·10⁻¹⁶ (computation.tex:82) | `build/reproducao/discos/*` | ✔ (0,99907378; 7,77·10⁻¹⁶; 175 pontos rechecados, todos ok) |
| nível 2: matriz idêntica bit a bit | `nivel2_A.json` | ✔ (sha256 igual, dmax = 0) |
| nível 2: resíduo de Schur idêntico | `nivel2_contagem.json` vs `rouche_L_rig/contagem.json` | ✔ (JSON idênticos) |
| nível 2: cauda do centro 3,0 idêntica em todo campo, 5,3 h | `nivel2_T3_3_0.json` vs `cauda/T3_3_0.json`; log | ✔ (23 campos, nenhum diferente; 19135 s) |
| nível 1: "all checks pass" | `build/reproducao/nivel1_*.json` | ✔ todos `ok: true` |
| nível 1 refaz "in exact rational arithmetic the final inequalities of every certificate" | `verifica_T2.py` | **exagero**: `simbolo_K` e `rouche_K` só leem o JSON gravado (o próprio código diz "recalculo … nivel 2"); `s3a` confere a cobertura com o disco de 1/8 (R6); não há verificação de q (R5) nem da parte livre −0,0232/0,325 (`free-part-toro.json`) |
| "symbol bounds, recombined from the stored bands" | idem | só para L (4 faixas); K é um arquivo único |
| Tabela 1: "symbols of B and B♯ (Arb, 2048×4096 grid)" | `symbol-K-1024x2048.json` | **errado para B♯** (1024×2048) |
| Tabela 1: 13 centros de fatores e de caudas | `rouche_L_rig/F`, `cauda` (10) + `faixaB/cauda` (3) | ✔ |
| Tabela 1: 154 tiles; F1 em segundos | `cobertura_lab/*`; nível 1 F1 31 s | ✔ |
| Apêndice C: digests | recalculados | ✔ os 11 digests e contagens de arquivos conferem com a tabela e com `digests.json`; 0 arquivos com sha diferente do inventário; ausentes no notebook só os pesados (`A.npy` e 3 `F_*.npy` da faixa B), como esperado |
| Apêndice C: "the list of files on which the corresponding certificate depends" | `inventario.json` | **incompleto**: (i) a linha da Prop 6.5 (13 arquivos) não contém os arquivos do Rouché de K (`build/s3b/rouche_K_circ.json`, `rouche_K_contagem_20x80.json`), que estão na linha de s_K; (ii) `build/spectrum/free-part-toro.json` (parte livre do Lema 5.3/Thm 5.1(i)), `build/rouche_L_rig/z/lab2/schur.json` (normas de φ e x usadas por todos os discos e janelas) e a fonte de β+ = 5191/500 não estão em linha nenhuma; (iii) `build/tiles_rig/tile-0_2-r0_04.json` está na linha da Prop 6.5 mas não é lido por `check_cover.py` |
| "a script downloads [RT data] and checks their sha256" | passo 0 do inventário (3 arquivos) | ✔ (sha de `RefA.dat`, `RefAplusB.dat` e do tarball conferem no disco) |
| `\todo{Repository and Zenodo record become public with the preprint.}` (computation.tex:102) | — | aparece em vermelho no PDF; remover antes da submissão; "Draft of \today" idem |

Não verifiquei (sem artefato no repositório): "five commodity computers, 4–6 cores, 7–31 GB", "a few days", "about two
hours on four cores" (a soma dos `tempo_s` gravados é compatível).

## R8. Literatura

- **Citação de RT na §1.2** (introduction.tex:32–35): idêntica ao TeX de RT v1, §1 (l.≈274), até "inconclusive"; o
  original continua com ", see [MG1] and [CHLP], and Sections 3.7 and 5.3 of [GM]". ✔ Também ✔: "We have taken the
  freedom…not formally stated local uniqueness" (Obs. 3.2(d)); "consistent with 3.445452402(3) for 2K in [MG2]" (§2.1).
- **Gundlach (1997), PDF desta pasta:** "It is real, with λ1 = −2.674 ± 0.009, corresponding to a critical exponent
  γ = 0.374 ± 0.001" e "γ = −1/λ1 = 0.374 ± 0.001" (Sec. e App. F); os modos crescentes estão "in the left half plane of
  λi, corresponding to modes that grow towards the singularity τ = −∞". Logo a convenção de sinal da Obs. 3.1 ✔ e
  "our interval lies within his error bar" ✔ ([2,674040, 2,674115] ⊂ [2,665, 2,683]). A Tabela II dele dá
  λ1 = −2,67410 com 1601 pontos, dentro do nosso intervalo.
- **Choptuik (1993):** γ ≈ 0,37 no resumo (PRL 70, 9) ✔; Δ ≈ 3,44 ✔.
- **Martín-García e Gundlach (1999)** (PDF de arXiv:gr-qc/9809059, Tabelas I e II): par ℓ = 2: κΔ = −0,07 (ωΔ/2π = 0,3)
  nas resoluções altas; ímpar ℓ = 2: κΔ = −2,30 ✔. A frase do texto ("least damped mode … κΔ ≈ −0,07") ✔. Conferência
  interna da §8.2: κΔ = 4π Re s = 4π·(−0,004764) = −0,0599 ✔ e ωΔ/2π = 2 Im s = 0,295 ≈ 0,3 ✔ (cálculo não rigoroso,
  apresentado como tal).
- **Reiterer (2015), arXiv:1510.05310:** resumo: "We identify a class of time-periodic linear symmetric hyperbolic
  equations that are finite codimension stable, because an associated operator has compact resolvent, sufficiently far to
  the right in the complex plane. This paper is an attempt to capture abstractly the observation in numerical general
  relativity that some discretely self-similar spacetimes, such as Choptuik's critical spacetime, are finite codimension
  stable." A §8.1 diz "with compact resolvents in the right half plane": o original diz "sufficiently far to the right";
  ajustar (imprecisão leve, I-12). O resto da frase ✔.
- **Metadados** (`r8_crossref.txt`, Crossref por DOI): os 18 DOIs da bibliografia conferem em título, revista, volume,
  páginas e ano (RT 2019: CMP 368(1), 143–186 ✔; Gundlach PRD 55, 695–713 ✔; GS 1971 Math. USSR Sb. 13, 603–625 ✔; Rump
  2006 BIT 46, 433–452 ✔; …). Detalhes menores: Lanford: Crossref dá "Lanford III O." ✔; a bibliografia de RT dá
  "427–434" para Lanford, Crossref e o artigo dão 427–435 ✔ (o artigo está certo).

## R9. §7.5 e §7.6 contra o repositório (feito depois de R1–R8; li só título e escopo dos relatórios)

Pastas `.codex-runs/*revisao*` (exceto esta), com a primeira(s) linha(s) do `RELATORIO.md`:

| pasta | título / escopo | corresponde a (§7.5) |
|---|---|---|
| 2026-09-24-revisao-s3a | "Revisão independente do certificado S3a"; "Não executei a cobertura inteira" | rodada 1 dos tiles; pela tarefa da rodada 2, "descrevendo o código em vez de verificá-lo". `paper/AI_USE.md` diz que a checagem de 24/09 com Codex `gpt-5.6-luna` "turned out to be superficial and is not counted" — compatível com esta pasta |
| 2026-09-24-revisao-s3a-r2 | **sem RELATORIO.md** (só TAREFA.md: "rodada 2") | — |
| 2026-09-24-revisao-agente | "Revisão independente de S3a (tiles de Perron)"; escopo: injetividade "fora de \|s − 0,401\| < 1/8" | tiles (Prop 6.5) |
| 2026-09-25-revisao-s3b | S3b | certificado em s_K (com Rouché de K) |
| 2026-09-25-revisao-s4 | S4 | certificado em μ |
| 2026-10-01-revisao-c1 | Rouché de L, faixa A | bandas |
| 2026-10-02-revisao-B | faixa B (C4) | bandas, inclusive s_B |
| 2026-10-02-revisao-T2 | T2; "agente Claude (Fable 5.1)" | lemas estruturais e cadeia |
| 2026-10-03-revisao-CRH | "Lema C, eixo imaginário e H-reg′"; "Claude Fable 5.1" | "modelling assumptions" (o commit 529a196 chama de "modelagem reduzida (Lema C, eixo imaginario, H-reg')") |
| 2026-10-03-revisao-G | Lema G | gauge |
| 2026-10-03-revisao-R | Lema R (ida) | reconstrução direta |
| 2026-10-03-revisao-V | a volta | recíproca |

Achados:
- **"There were ten such checks … the tiles of Prop 6.5, in two rounds"** (computation.tex:117–119). Há 10 relatórios
  contáveis de agentes Claude (agente, s3b, s4, c1, B, T2, CRH, G, R, V) — dez, como diz o texto — mas só **um** é sobre os
  tiles. A "segunda rodada" só pode ser a de `revisao-s3a` (Codex, que o `AI_USE.md` diz não contar) ou a de
  `revisao-s3a-r2` (sem relatório). Se as duas rodadas contam, as checagens são 11 e uma não é de `claude-fable-5-1`
  (contra §7.6); se não contam, "in two rounds" é falso. Corrigir (achado E-4).
- A checagem dos tiles conferiu a região fora do disco de raio **1/8**, não 0,05 (ver R6/E-1); o texto da §7.5 não diz
  isso.
- "it was not allowed to read the author's derivations or conclusions, only the precisely stated claims and the code and
  data" (computation.tex:111–112): as revisões G, R e V tinham como objeto exatamente documentos de derivação do autor
  (`docs/GAUGE_GLOBAL.md`, `docs/HREC_IDA.md`, `docs/GEOMETRIC_FLOQUET_RECONSTRUCTION.md`, conforme o título/escopo), e a
  de tiles leu "S3 §10.4–10.13". Melhor: "only the statements and derivations under review, the code and the data, not
  the author's conclusions, notes or earlier audits" (imprecisão leve, dentro de I-10).
- Exemplos citados ("an early version of the tiles used the dyadic approximation of μ"; "the deflation error outside the
  finite box had been omitted"): há rastro no código — `check_cover.py` ("mu verdadeiro de RT … (revisao independente,
  24/09)") e `rouche_L_disco.py`/`rouche_L_janelas.py` ("deflação verdadeira fora de Z (revisão de 01/10, L1)") ✔.
- **Identificadores de modelo.** `paper/AI_USE.md` e §7.6 coincidem: Codex `gpt-5.6-luna`, `gpt-5.6-sol`,
  `gpt-6-astra`, `gpt-reserve`; Claude `claude-opus-5` (11–23/09), `claude-opus-5-5` (23/09–08/10); checagens
  `claude-fable-5-1`; ChatGPT "Sol 6.1". Trailers do `git log`: 51 commits "Claude Opus 5" (17–23/09) e 159 "Claude Opus
  5.5" (23/09–08/10) — consistentes com as datas do AI_USE.md (o histórico começa em 17/09, então 11–16/09 não aparece no
  git, o que não é contradição). Nenhum trailer de Fable (as revisões não fizeram commits) nem de Codex (esperado). Os
  relatórios T2 e CRH se identificam como "Fable 5.1"; c1 e B como "outro modelo"; os demais não dizem o modelo nas
  primeiras linhas, então "all checks by claude-fable-5-1" não é verificável só por elas.

## R10. Síntese

Gravidades como na TAREFA. "Sem efeito no teorema" quer dizer que o Teorema Principal continua provado depois da
correção sugerida, sem conta nova (salvo onde indicado).

### Erros (afirmação falsa ou não sustentada) — 4, nenhum atinge o Teorema Principal

| id | local | achado | evidência | correção |
|---|---|---|---|---|
| E-1 | computation.tex:74–78 (§7.3), 13 e 115 (§7.1, §7.5); app_certificates_table.tex:13 | O verificador de nível 1 e a revisão independente dos tiles conferiram a cobertura da região **fora do disco de raio 1/8**, não de 0,05 como diz a Prop 6.5; a coroa 0,05 ≤ \|s − 0,401\| < 1/8, que o Rouché de K não cobre, só estava coberta pelo planejador que gerou os tiles (`tile_cover.py` rodado com `--disco 0.401,0.05`, como gravado nos JSON), não por nenhuma das checagens. A Prop 6.5 **vale** (conferi agora). | `check_cover.py` (padrão `--disco 0.401,0.125`), `verifica_T2.py::check_s3a`, título/escopo de `2026-09-24-revisao-agente`; `r6_cobertura_disco005.txt` (2871 células, COMPLETA) vs `r6_cobertura_disco0125.txt` (1770) | passar `--disco 0.401,0.05` em `check_s3a` (ou mudar o padrão), rodar de novo, e dizer na §7.5 o que a revisão de tiles cobriu |
| E-2 | computation.tex:10–13, 69–80 | "Exact rational arithmetic: … the deflation polynomial q" e "Level 1 re-checks, in exact rational arithmetic, the final inequalities of every certificate": não há código que verifique q; `simbolo_K` e `rouche_K` só leem o JSON gravado; a parte livre (−0,0232 e 0,325) não é conferida | `verifica_T2.py` (`check_simbolo_K`, `check_rouche_K`); busca sem resultado por verificação de q; R5, R7 | escrever o que cada check faz; ou acrescentar os checks (eu refiz o Perron do Rouché de K em racionais: 0,8125577; e q à mão) |
| E-3 | app_paper_lemmas.tex:71; setting.tex:310; physical_modes.tex:230–231 | β+ ≤ 5191/500 chamado de "published bound of [RT]": não é de RT; é calculado pelo projeto e não está no inventário nem no nível 1 | `scripts/independent_component_beta.py`; `r1_beta_saida.txt` (β_η = 10,3810817 ≤ 10,382, confere) | "a bound computed from the data and the ball of [RT] (script …)", e pôr o script/saída no inventário e no nível 1 |
| E-4 | computation.tex:117–119; paper/AI_USE.md | "ten such checks … the tiles of Prop 6.5, in two rounds": só há um relatório contável sobre os tiles; a outra "rodada" é a de Codex que o AI_USE diz não contar, ou a pasta `revisao-s3a-r2` sem relatório | R9 | "the tiles (one check; an earlier attempt with another tool was superficial and is not counted)" |

### Lacunas (prova incompleta, reparável) — 2

| id | local | achado | evidência | correção |
|---|---|---|---|---|
| L-1 | physical_modes.tex:35–39 (Def 4.1(iii)); 296–317 ((R3)–(R6)) | Analiticidade só numa vizinhança E de [−1,1], mas M̂ vai até ξ < 41/40; (R3)–(R6) identificam δω' em todo M̂ com a série de h sem dizer por quê na faixa 1+ε < ξ < 41/40 | R4.4 | uma frase: além do cone ξ = const é espacial (as duas características têm dξ/dτ > 0) e a solução em ξ > 1+ε é determinada pelos dados em ξ = 1+ε; ou pedir E ⊃ [−41/40, 41/40] |
| L-2 | roots.tex:152–154 (Cor 6.6) | aplica o Lema 2.4 (enunciado para Y(a,b)) ao autovetor do **toro** | R6 | inserir "by Lemma 2.6, h ∈ Y+ since \|s*\| ≤ 3,1" (idem para h_* em s_K) |

### Imprecisões (enunciado ou número a ajustar) — 13

| id | local | achado | evidência | correção |
|---|---|---|---|---|
| I-1 | counting.tex:251; :128; setting.tex:309 e app_paper_lemmas.tex:72; physical_modes.tex:94,100 | quatro números impressos do lado **inseguro** (literalmente falsos, sem efeito): ‖I−U*U‖_F ≤ 1,30·10⁻⁶ (é 1,3045·10⁻⁶); ‖φ_i‖ ≤ 5,5 (é 5,5013); "≤ 0,2665" (é 0,266547); Im c_1(·,−1) "0,097753…" (é 0,0977527…) | `r1_saida.txt` | 1,31·10⁻⁶; 5,51; 0,2666; 0,097752 |
| I-2 | counting.tex:14,59,74 | −0,0232 e 2,926 não saem do certificado gravado (−0,023196 com μ ≥ 1/6 ⇒ 2,926002); valem com o μ verdadeiro (−0,023424 ⇒ 2,92577) | `free-part-toro.json`; R1 N4 | dizer "with μ = μ_RT" ou imprimir 2,927 |
| I-3 | app_paper_lemmas.tex:68; `build/t2/l_real_toro.json` | o certificado L-real usou ‖B‖ = 3,964101, abaixo da cota certificada 3,9641333; com a cota certificada dá 0,370838 < 0,371 | R1 N8 | refazer `l_real_toro.py --nB 3.9642` (barato) |
| I-4 | counting.tex:64,68; computation.tex:48 | grade de B♯ é 1024×2048, não 2048×4096; "minimum slack 0,0771/0,0535" é a folga **máxima** somada (0,077197/0,053578), arredondada para baixo | `symbol-K-1024x2048.json`, `symbol-L-faixa*.json` | corrigir a grade; "the Lipschitz slack added is at most 0,0772 for L and 0,0536 for K" |
| I-5 | counting.tex:146–164, 187–189 (§5.3, Lema 5.5) | a referência certificada usa a matriz de ponto flutuante 𝖠 (dados ω_A, x_j^A), "tomada exata", e não P_ZL̃P_Z verdadeiro; e t_p majora ‖(𝖠+p)^{-1}‖, não ‖Ĥ_ZZ(p)^{-1}‖ | docstrings de `rouche_L_contagem.py`, `rouche_L_rig.py`; R5 | definir o bloco Z da referência com 𝖠 e dizer que a diferença entra em E; trocar a definição de t_p |
| I-6 | counting.tex:84–88 | Rouché de GS em homotopia: falta a hipótese de que R + t(H−R) é Fredholm de índice 0 (verdadeira aqui, H − R compacto) | R5 | acrescentar "H − R compact" |
| I-7 | setting.tex:113–119 | "[RT, §3–4]" para pesos, dados e cotas finais (7): estes estão em RT §5 (e §4 para 𝓛_2) | R3 | "[RT, §§4–5]" |
| I-8 | physical_modes.tex:104–112 (Lema 4.5) | o Passo 0 usa que L_X(ḡ,φ̄) é SO(3)-invariante, o que o enunciado não diz | R4.1 | "such that L_X(ḡ,φ̄) is spherically symmetric and in RT form" |
| I-9 | setting.tex:278–289 (Prop 2.5) vs Lema 4.8 e App. A.3 | a Prop 2.5 é enunciada só para a ∈ {65/64, 129/128}, b ∈ {5/4, 9/8}, mas o Lema 4.8 usa J^{-1} e Fredholm em Y(κ1',κ2') com κ1' = 1 e κ2' qualquer em (1, 9/8] | R4 | enunciar para 1 ≤ a ≤ 65/64, 1 < b ≤ 5/4 (a prova do App. A.2 já cobre) |
| I-10 | computation.tex:110–113 | "not allowed to read the author's derivations": G, R e V revisaram justamente documentos de derivação, e a de tiles leu "S3 §10.4–10.13" | R9 | "only the statements and derivations under review, code and data; not the author's conclusions or earlier audits" |
| I-11 | app_certificates.tex:43–49; app_certificates_table.tex:13 | o inventário por alegação está incompleto: os arquivos do Rouché de K estão na linha de s_K, não na da Prop 6.5; faltam `free-part-toro.json`, `rouche_L_rig/z/lab2/schur.json` e a fonte de β+; `tile-0_2-r0_04.json` está listado mas não é usado | R7 | completar `inventario.json` e regerar a tabela |
| I-12 | discussion.tex:17–19 | Reiterer (2015): "compact resolvents in the right half plane"; o resumo diz "sufficiently far to the right in the complex plane" | R8 | "far enough to the right" |
| I-13 | introduction.tex:54–55 | versão informal: "four roots … all real"; o teorema diz "has a real representative" (s_B aparece como s_B + i/2) | R2 | "each with a real representative" |

### Cosméticos — 9

| id | local | achado |
|---|---|---|
| C-1 | computation.tex:102; main.tex:25 | `\todo{…}` aparece em vermelho no PDF; "Draft of \today" |
| C-2 | counting.tex:229 | "89 sample points": são 82 pontos de amostra distintos (101 avaliações); 89 inclui os 10 centros (`amostras` de `rouche_L_janelas.py`) |
| C-3 | roots.tex:35 vs tabela | o Z2 impresso majora ‖𝒜D²F‖, não ½‖𝒜D²F‖ (conservador por fator 2); dizer |
| C-4 | roots.tex:94 | "‖x_g − x_0‖ ≤ 1,82·10⁻⁵ < r, with r ≤ 6,07·10⁻⁵" conclui dist < r de uma cota superior de r; usar r ≥ 6,06·10⁻⁵ ou o raio de unicidade 3,26·10⁻⁴ |
| C-5 | physical_modes.tex:100 | "0,072916…", "−0,079473…" são arredondamentos (0,0729159…, −0,0794728…), não truncamentos |
| C-6 | setting.tex:216–218 | Q̃ = {−1/4 < Im ≤ 3/4} e "lower half {\|Im s\| ≤ 1/4}" (inclui Im = −1/4, fora de Q̃) |
| C-7 | physical_modes.tex:109 | "Equivalently, the SO(3)-average of X is cX_0": é consequência, não equivalência |
| C-8 | roots.tex:31–40 | o teorema de NK precisa de 𝒜 injetiva (aqui é invertível por construção); dizer |
| C-9 | physical_modes.tex:419–421 | o quociente ker(D\|E)/(G ∩ E) pressupõe G ∩ E ⊂ ker(D\|E) (vale porque C(μ)g* = 0); dizer |

### Os cinco achados mais importantes
1. **E-1** — a cobertura da coroa 0,05 ≤ \|s − 0,401\| < 1/8 na Prop 6.5 só era garantida pelo planejador dos tiles;
   o nível 1 e a revisão independente conferiram com 1/8. Refeita aqui com sucesso (`r6_cobertura_disco005.txt`).
2. **E-2/I-11** — a §7 superestima o que o verificador de nível 1 refaz (q, símbolo de K, Rouché de K, parte livre), e o
   inventário por alegação está incompleto.
3. **L-1** — a Def 4.1 só dá analiticidade perto de [−1,1], mas a correspondência é afirmada em M̂ (até 41/40); falta o
   argumento de domínio de dependência além do cone.
4. **I-5** — a §5.3/Lema 5.5 descrevem uma referência (truncação verdadeira, t_p = ‖Ĥ_ZZ^{-1}‖) diferente da certificada
   (matriz de dados 𝖠, t_p = ‖(𝖠+p)^{-1}‖); a lógica fecha, o texto não descreve o que foi feito.
5. **I-1/I-2/I-3/E-3** — números auxiliares impressos do lado inseguro ou não sustentados pelo arquivo correspondente
   (1,30·10⁻⁶; 5,5; 0,2665; 0,097753; −0,0232/2,926; ‖B‖ no L-real; β+ "publicado").

### Veredicto

**Sim: o artigo, tal como está, sustenta o Teorema Principal**, no sentido de que não encontrei afirmação falsa ou
lacuna que comprometa nenhuma das partes (1) e (2). Todos os enclausuramentos impressos das raízes, a tabela de NK, os
θ dos discos, janelas, ε_B·t, testemunhos, a identificação de μ e os intervalos de λ e 1/λ (refeitos em Arb com o K de
RT) conferem com os certificados e estão do lado seguro. As provas reescritas (Lema 4.5 com a variante (b), passo R7,
Prop 4.12 inclusive o cone fixo, Thm 4.11, Lemas 5.4–5.8, Lema 6.1, Props 6.3, 6.5, 6.7) estão corretas; refiz as
derivações e testei R7/Prop 4.12 e o Lema 5.4 numericamente. O único ponto em que a sustentação dependia de algo não
verificado — a cobertura da coroa na Prop 6.5 — foi conferido nesta revisão e vale. Antes da submissão é preciso
corrigir E-1 a E-4 (afirmações sobre verificação e atribuição), fechar L-1 e L-2 (uma frase cada) e ajustar os números
de I-1 a I-3.

## Manifesto do snapshot (fim) e artefatos

`r0_manifesto_fim.txt` (16:10 de 08/10): `sha256sum -c` sobre as 1159 entradas — nenhuma linha diferente de `OK`;
`git status --untracked-files=no` vazio (nenhum arquivo rastreado modificado). Não usei o laboratorio2 nem outra máquina.
Scripts do repositório foram só lidos ou executados sem gravar (`check_cover.py`, `independent_component_beta.py` sem
`--out`); o `digest_certificados.py` foi copiado para `r7_digests.py`, que grava só aqui. RT foi extraído no
scratchpad da sessão. Não li nenhum arquivo proibido; em `.codex-runs/` li só título e escopo dos relatórios (R9) e as
cinco primeiras linhas de `2026-09-24-revisao-s3a-r2/TAREFA.md` (essa pasta não tem relatório).

Artefatos nesta pasta: `r0_manifesto_inicio.txt`, `r0_manifesto_fim.txt`; `r1_numeros.py` → `r1_saida.txt`;
`r1_beta_saida.txt`; `r4_jato.py` → `r4_jato_saida.txt`; `r5_brauer.py` → `r5_brauer_saida.txt`;
`r6_cobertura_disco0125.txt`, `r6_cobertura_disco005.txt`; `r7_digests.py` → `r7_digests_saida.txt`; `r8_crossref.txt`.
