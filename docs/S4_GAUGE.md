# S4: a raiz de gauge em s = μ (25/09/2026)

## O que S4 precisa

A contagem física (T2) é: raízes de L no setor A com multiplicidade (winding, C3: **3**),
menos as que violam as constraints (S3b: exatamente uma, `s_K ≈ 0,401`, simples), menos o
gauge. A classificação do gauge residual (`GAUGE_DOMAIN_CLASSIFICATION.md`,
`GAUGE_RT_ACTION.md`) dá um único gauge instável: o modo n = 0, **s = μ**, translação conjunta
das coordenadas nulas (desloca o instante T* da singularidade; move o cone), com gerador exato

```
g = (Z + 1) ω*,   Z = μ⁻¹ ∂τ + ξ ∂ξ,   L_*(μ*) g = 0,   D g = 0 (constraints),
```

e, para qualquer fundo w, `L_w(μ)(Z + 1)w = (Z + 2) F_w` (F_w = defeito de w). Pela álgebra de
`CONSTRAINT_ROOT_QUOTIENT.md`, na raiz μ a dimensão física é `dim E − rank(D|E) − 1`: se
**dim E = 1** (μ algebricamente simples), E = span(g), D g = 0 e a contribuição física é 0.

**Parte computacional (FECHADA em 25/09):** μ* é raiz algebricamente simples de L (seção abaixo).
**Parte conceitual (redigida abaixo):** a realização com cone móvel, em que o modo n = 0 é
quocientado; fixação global e reconstrução como hipóteses explícitas (H-rec, H-gauge, H-cone).

## Parte computacional: NK centrado no vetor de gauge

1. **Vetor de gauge na base de L** (`gauge_vector_L.py`): `g_A = (Z + 1) ω_RefA`, com ∂τ → i m/2
   e ξ∂ξ nos coeficientes simétricos. Conferências: |cos(g_A, autovetor da truncagem perto de
   μ)| = 0,99999992 (12×48); **‖L(μ) g_A‖/‖g_A‖ = 4,4·10⁻¹⁰** sem pesos, na saída completa (o
   defeito de RefA; a auditoria de 16/09 mediu 1,1·10⁻¹⁰), 4,3·10⁻⁷ nos pesos de L. A cauda de
   g_A fora de Z = 32×128 é **2,6·10⁻⁷** relativa.
2. **RefB** (`gauge_refB.py`, diferença EXATA RefAplusB − RefA em inteiros, 450×1350 modos):
   `‖RefB‖_L = 2,35·10⁻¹⁰`, **`‖(Z + 1) RefB‖_L = 2,45·10⁻⁸`**; μ_RefA + μ_RefB bate com o μ de
   80 dígitos (|μ_RefB| = 3,881·10⁻¹¹).
3. **NK de L em λ0 = μ_RefA**, mesma arquitetura de S3b (Z' 32×128 bordejado com
   `h0 = P_Z g_A` normalizado — `nk_L2_Z.py --h0` —, cauda com bloco-Jacobi em 26 janelas,
   montagem de 3 níveis): autopar único x* na bola e simplicidade algébrica.
   - **Etapa A rigorosa** (`build/s4/A.txt`, `rig/Z.json`; laboratório, 1123 s):
     `‖M̂⁻¹‖ <= 52,14972` (Loewner), `a <= 0,40461`, `ρ_Z <= 1,06·10⁻⁵`, `Y_Z <= 3,00·10⁻⁶`,
     `ε_far <= 0,0154`, `‖Q_ZZ‖ <= 3,3950`; erro do fundo `rt = 4,7·10⁻⁸` (RefB), `ε_μ = 2,8·10⁻⁸`.
     Normas restritas de Q0 (`rig/Q.json`): `q_ZT = 0,18412`, `q_TT = 0,32193`, extra 0,17892.
   - **Estimativa da cauda** (`build/s4/nk3_mu.json`, 7603 s): **θ0 = 0,860** (pesos 1), raio NK
     ~5,7·10⁻³. A cauda é quase a mesma de 0,401 (janela ±32..47: ‖V‖ = 1,741 contra 1,692,
     acoplamento a Z' 0,36–0,38 contra 0,36); o que muda em μ é o nível Z' (‖M̂⁻¹‖ 52 contra 21).
   - **Etapa B rigorosa** (Loewner por janela): laboratório 11 janelas em ordem direta, laboratorio2
     15 em ordem reversa (checkpoints complementares `rig/T3.lab.json`, `rig/T3.lab2.json`); junção e
     far no laboratorio2 (`B_final.txt`, 1219 s) → `rig/T3.json`. Janelas: ‖V‖ de 1,045 a 1,741,
     acoplamento a Z' até 0,394 (±32..47), ρ até 9,6·10⁻³ (central).
   - **Etapa C** (`C.txt`, `rig/C3.json`): Perron de 3 níveis em racionais, **θ = 0,880172**,
     Z2 = 183,7, Y = 6,59·10⁻⁶, discriminante (1 − θ)² − 4 Z2 Y = 0,0095, **r = 6,06·10⁻⁵**:
     zero único (y*, λ*), |λ* − μ_A| <= 2,8·10⁻⁵, DF(x*) invertível ⇒ λ* algebricamente simples.
4. **Identificação** (`nk_L3_E.py`, `E.txt`, `rig/E.json`): o autopar de gauge exato x_g = (Q0⁻¹ g_n, μ*),
   g_n = g*/⟨h0, g*⟩, está na bola, logo x_g = x* e μ* é simples. Pela equação de autovalor,
   `y_g − y0 = −B_*(g_n − h0) − L_*(λ0) h0 + (λ0 − μ*) g_n − δμ J1 (g_n − h0)`, então
   `‖y_g − y0‖ <= (‖B‖ ‖g_n − h0‖ + ‖L_*(λ0) h0‖ + |λ0 − μ*| ‖g_n‖)/(1 − |δμ| ‖J1 Q0‖)`, sem
   derivadas do erro de RT. `g* − g_A = (Z + 1)(RefB + Corr) + (μ*⁻¹ − μ_A⁻¹) ∂τ RefA`: RefB
   calculado; Corr ω (~2⁻²⁷⁷ na bola de RT) entra por ‖(Z + 1) O_μ⁻¹‖ e pela equação de ponto
   fixo de RT (`O Corr = termos em SSPACE`). Estimativa: ‖y_g − y0‖ ~ 10⁻⁶ contra r ~ 10⁻⁴.
   - **Lema de Corr (revisão de S4, L2 e L3; `corr_gauge_bound.py`, `build/s4/corr_gauge.json`).**
     Pela definição de C em RT (fdrrihifuhr78), com SΩ⁺(μ*, ω*) = 0 e O_μ = ∂τ + μ D_O afim em μ:
     O_{μ*} Corr = G := −SΩ⁺(Ref) − Corr_μ(D_O + SΓ1)Ref − μ* SΓ1 Corr − 2 SΞΓ2(Ref, Corr) −
     SΞΓ2(Corr, Corr). Com as constantes publicadas por RT (𝓛3 = 2⁻²⁹⁴, 𝓛1 = 10,5, 𝓛2 = 2⁻²⁵,
     𝓚1 = 80/9, bola R = 2⁻²⁷⁷ com κ* = 1) e ‖RefA‖_SSPACE = 14,425 (exata dos dados),
     ‖G‖ <= 5,5·10⁻⁸⁰. Na componente B, D_OB = ξ∂ξ + 1, logo μ*(Z* + 1) = O_B e
     (Z* + 1)Corr_1 = G_1/μ*. Nas componentes A, D_OA = (1 + ξ)∂ξ + 1. Aí ∂τ Corr = ∂τ O⁻¹ G, com
     ‖∂τ O⁻¹‖ <= 18 pela forma fatorada de O⁻¹ e pelas cotas de banda de RT. Como (1 + ξ)∂ξ =
     (1 + S)(1 − S)⁻¹ n, vale ‖n Corr‖ <= 9 ‖(1 + ξ)∂ξ Corr‖, e ξ∂ξ = D_OB − 1 fecha a conta.
     Resultado: **‖(Z* + 1) Corr‖_L <= 1,1·10⁻⁷⁵** (‖·‖_L <= √2 ‖·‖_SSPACE), contra os 10⁻⁶⁰ usados.
     **Domínio:** g* = g_A + e com ‖e‖_L finito, logo g* ∈ Y. L_*(μ*) g* = 0 vale coeficiente a
     coeficiente, então (J + λ0) g_n ∈ Y. J + λ0 é injetivo no domínio máximo (as soluções
     homogêneas (1 + ξ)^{−(1 + (λ0 + im/2)/μ)} e ξ^{−(…)} são singulares dentro da elipse), logo
     y_g = Q0⁻¹ g_n ∈ Y. Isso dispensa `RADIUS_RECOVERY.md`.
   - **Critério** (`nk_L3_E.py`): F é quadrática, então ‖DT(x)‖ <= θ + 2 Z2 ‖x − x0‖_v para
     todo x; com ρ = max(‖x_g − x0‖_v, r), se θ + 2 Z2 ρ < 1, T é contração em B(x0, ρ) e o único
     zero lá é x* (que está em B(x0, r)); logo x_g = x*. A escala de unicidade é (1 − θ)/(2 Z2)
     (S3b: 1,9·10⁻³), bem maior que o r mínimo.
   - **Prévia** (sem C3; `CHOPTUIK_EPS_FUNDO_L=4.7e-08`): ‖P_Z g_A‖ = 2,60435, cauda 6,9·10⁻⁷;
     ‖g* − g_A‖ <= 2,5·10⁻⁸; ‖g_n − h0‖ <= 2,8·10⁻⁷; ‖L_A(λ0) h0‖ = 3,71·10⁻⁶ (arredondamento
     calculado, 6,7·10⁻¹¹); com RT e μ, 3,76·10⁻⁶; **‖y_g − y0‖ <= 4,9·10⁻⁶**.
   - **Resultado:** com v = (0,464; 1; 0,269), **‖x_g − x0‖_v <= 1,82·10⁻⁵ <= r = 6,06·10⁻⁵**
     (até o raio mínimo basta); θ + 2 Z2 ρ = 0,902 < 1, unicidade até 3,26·10⁻⁴. **IDENTIFICADO:
     x_g = x*, então μ* é raiz algebricamente simples de L, com autoespaço span(g*).** Como D g* = 0,
     pela álgebra de `CONSTRAINT_ROOT_QUOTIENT.md` a contribuição física de μ é 0.

## Parte conceitual: a realização com cone móvel

**O que o modo n = 0 é.** As coordenadas de RT são centradas no evento de acumulação
(u₊ = u₋ = 0), e u_σ = −(1 + σξ)e^{−μτ}. A translação conjunta u_± → u_± + ε desloca esse
evento ao longo do centro, ou seja, muda o instante T* da singularidade, e com ele o cone de luz
passado. Na linearização em torno do fundo, essa translação aparece como o campo
X = e^{μτ}(μ⁻¹∂τ + ξ∂ξ), e seu efeito nos campos de RT é exatamente g = (Z + 1)ω*, em s = μ
(`GAUGE_RT_ACTION.md`). É o "modo de gauge" da literatura de fenômenos críticos: uma solução crítica com outro T* é a mesma
solução. O seu expoente depende das coordenadas: aqui s = μ*, que no tempo logarítmico τ_log = (Δ/4π)τ
dá 0,614, e não o "1" do gauge de fatiamento usual (revisão de T2, 03/10).

**Qual realização a contagem usa.** O espaço espectral de RT (campos em (τ, ξ) fixos, nos pesos
(65/64, 5/4)) CONTÉM g, no domínio forte (`RADIUS_RECOVERY.md`). Logo a realização espectral é
naturalmente a de **cone móvel**: nela a mudança de T* é uma perturbação admissível, e deve ser
quocientada. Na realização de cone fixo (f(0) = 0), g é excluído e não há gauge positivo nesta
classe (`GAUGE_DOMAIN_CLASSIFICATION.md`). As duas não podem ser misturadas numa mesma contagem.
A contagem de S1–S4 é feita inteira na realização espectral, isto é, de cone móvel.

**A faixa fundamental** (revisão de S4, L1). Σ_A é invariante por s ↦ s + i (`SECTOR_B.md` §4;
no código, J_m(s + i) = J_{m+2}(s) e os blocos de B só dependem de d): μ + ij, j ∈ ℤ, são todas
raízes de gauge, e "raízes em Re s > 0" sem faixa é um número infinito. Pelo Lema 5.1 de
`SECTOR_B.md`, contar L_A em Q̃ = {−1/4 <= Im s < 3/4} equivale a contar os dois setores do
problema de período completo. Q̃ = faixa de Γ_A (|Im s| <= 1/4; C3) ∪ faixa B (1/4 < Im s < 3/4;
C4, onde fica a raiz 0,0414 + 0,5i que viola as constraints). Os aliases de μ têm Im inteiro:
em Q̃ só μ aparece, e na faixa B não há gauge.

**A conta, nessa realização.** Com N_A = número de raízes de L_A em Re s > 0, |Im s| <= 1/4, com
multiplicidade algébrica, na realização RT com sinal (winding; C1–C4, T1, F2–F3):

```
N_físico(faixa A) = N_A − 1 (s_K, viola as constraints: S3b) − 1 (μ, gauge n = 0: S4) ,
```

e a faixa B entra à parte (C4: multiplicidade física 0). Isso vale desde que:

- (i) s_K e μ sejam algebricamente simples (S3b; S4, parte computacional);
- (ii) em μ o autoespaço seja span(g) com D g = 0 (S4: identificação). A simplicidade e D g = 0
  são necessárias juntas: sem simplicidade, uma cadeia sobre g poderia contribuir;
- (iii) nenhuma outra raiz na faixa seja gauge. A classificação de Floquet do gauge residual
  analítico dá expoentes Re s = μ(1 − n), n >= 0; o único positivo é n = 0. A fase (n = 1) está
  em s = 0, fora de Re s > 0, e os demais têm Re s < 0. A completude vale para germes analíticos
  f(u). A versão global, para gauges C¹, é o Lema G de `GAUGE_GLOBAL.md` (03/10);
- (iv) nas demais raízes s0 da faixa, rank(D|E_{s0}) = 0, isto é, o espaço de raízes inteiro
  satisfaz as constraints. Onde K(s0) é injetivo, isso segue de K C = M L, também ao longo de
  cadeias: de L h0 = 0 e L h1 + L' h0 = 0 vêm C h0 = 0 e K(C h1 + C' h0) = 0 (derivando K C = M L),
  logo C h1 + C' h0 = 0. K é injetivo na região de S3a e no disco de S3b fora de s_K; em s_K,
  rank(D|E) = dim E = 1 (S3b: simples, testemunho > 0). S3a cobre 0 <= Re s <= 1,765,
  |Im s| <= 1/4; o resto da faixa A (1,765 < Re s < R, com R <= 2,93 de F3) e a faixa B pertencem
  às obrigações abertas F3 e C4;
- (v) N_A seja contado com multiplicidade algébrica no mesmo espaço, setor e realização de S4
  (espaço com sinal, pesos de RT, cone móvel). É o que C3 e T1 devem entregar.

**Hipóteses que ficam explícitas** (físicas, não computacionais; enunciadas no teorema final):

- **(H-rec) Reconstrução:** toda perturbação linear física (solução das equações de
  Einstein–campo escalar linearizadas, esfericamente simétrica, analítica nos domínios de RT)
  pode ser posta no ansatz de RT por um gauge admissível; e, reciprocamente, todo elemento de
  ker L(s) que satisfaz as constraints vem de uma perturbação física.
  `GEOMETRIC_FLOQUET_RECONSTRUCTION.md` faz a recíproca dentro do ansatz para Re s > 0.
- ~~(H-gauge) Completude do gauge residual~~: **provada em 03/10** (Lema G, `GAUGE_GLOBAL.md`):
  globalmente em M, para gauges C¹, sem analiticidade dos germes.
- **(H-cone) Cone móvel: virou o Lema C em 03/10** (`GAUGE_GLOBAL.md` §6). Texto original: a noção física de instabilidade quocienta a mudança de T*. É a
  convenção da literatura de fenômenos críticos; a de cone fixo daria uma contagem diferente
  e não é a usada.
- **(H-reg) Regularidade** (revisão de S4, L6): a classe física de perturbações é a analítica
  nos raios de RT, ou, equivalentemente para a contagem, o espectro em Re s > 0 não depende do
  espaço. `RADIUS_RECOVERY.md` cobre só a troca entre dois raios analíticos, e noutra realização.

Com (H-rec), (H-gauge), (H-cone) e (H-reg), com os certificados e com (iv) e (v) fechados
pelas obrigações abertas, a contagem física na faixa A é N_A − 2.

## O erro do fundo com RefB explícito (25/09, noite)

Em μ, o nível bordejado é ~2,5× pior condicionado que em 0,401 (`‖M̂⁻¹‖ <= 52,15` contra 20,67),
porque a raiz de fase s = 0 está a 0,168. Com o erro de RT tratado como 40 ε_ω = 1,19·10⁻⁶, o piso
do resíduo seria ~6·10⁻⁵ e o NK não fecharia. Mas RefB é conhecido: `fundo_refB.py` calcula, da
diferença EXATA RefAplusB − RefA, a cota de Young **‖B_L(RefB)‖ <= 4,67·10⁻⁸** (todas as 899
frequências |d| <= 449). O fundo verdadeiro é RefA + RefB + Corr, com ‖Corr‖ ~ 2⁻²⁷⁷ (termo 1e-60,
generoso). As etapas de S4 usam `CHOPTUIK_EPS_FUNDO_L=4.7e-08` no lugar de 40 ε_ω (registrado nas
saídas via `dB_L['rt']`). O μ verdadeiro continua em `eps_mu_L`.
