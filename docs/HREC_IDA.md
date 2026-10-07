# H-rec (ida): todo modo físico de Floquet entra no ansatz de RT (Lema R, 03/10/2026)

**Estado: provado e revisado.**
- A revisão independente de 03/10 (`.codex-runs/2026-10-03-revisao-R/`) não quebrou o lema e refez com
  artefatos próprios as peças de maior risco: R1, R2, R4 e R5.
- As lacunas L1–L9 eram de escopo e redação e estão incorporadas abaixo:
  - L1: todo s, não só |s| <= 3,1;
  - L9: cadeias, na peça nova R7;
  - L4: R2 completo;
  - L6: R3 completo.

## 1. Definições e enunciado

**Notação.**
- ḡ é a métrica do fundo, e g* = (Z + 1)ω* é o perfil de gauge (S4).
- M̂ = M × S² mais o eixo, com M o domínio de RT, que contém o cone.
- Θ²: τ ↦ τ + 4π. Pelo Teorema 1 de RT, (Θ²)*ḡ = e^{−4K}ḡ e φ*∘Θ² = φ*.

**Modo físico de Floquet (definição).** É uma perturbação linear δ = (δg, δφ), esfericamente simétrica,
tal que:
- δ resolve as equações de Einstein–campo escalar linearizadas em torno de (ḡ, φ*) em M̂;
- δ está na classe **H-reg′** (03/10). Escrevendo δ pelos coeficientes de R3 (α, β, γ, d e δφ, em
  (τ, x) cartesianas, com o fator de Floquet removido), esses perfis estão em
  C^∞(ℝ/4πℤ; A(E)), para alguma vizinhança complexa E de [−1, 1] no plano de ξ, com o domínio de RT
  indo até ξ = 41/40. Aqui A(E) são as funções holomorfas limitadas em E. Isto é: **C^∞ em τ e
  real-analítica em ξ**, com vizinhança complexa uniforme em τ. A versão anterior pedia analiticidade
  real conjunta em M̂, o que é mais forte.
- δ satisfaz (Θ²)*δg = e^{4πs}e^{−4K}δg e (Θ²)*δφ = e^{4πs}δφ.

O fator e^{−4K} é o do fundo: isso equivale a δζ, δW/W e δφ serem e^{sτ} vezes funções 4π-periódicas.
Um modo **generalizado** de ordem k (uma cadeia física) é uma solução da forma e^{sτ}Σ_{j<k}τ^j·(perfis
4π-periódicos analíticos), no mesmo sentido.

**Gauge.** δ ↦ δ + L_X(ḡ, φ*), com μ e as coordenadas (τ, x) fixos, e X:
- de classe C¹ em M̂ (Lema G de `GAUGE_GLOBAL.md`); ou
- contínuo em M̂, C¹ fora do cone e com δω_X na classe analítica (Lema G′).

Os gauges construídos abaixo são analíticos (às vezes seculares em τ, em R2 e R7), logo estão nas duas
classes.

**Lema R (H-rec, ida).** Seja δ um modo físico de Floquet com Re s > 0.
1. Existe X real-analítico em M̂ e Floquet tal que δ' = δ + L_X(ḡ, φ*) está na forma do ansatz de RT.
2. As variáveis de RT de δ' são δω = e^{sτ}h, com h ∈ D(Y+), L(s)h = 0 e C(s)h = 0. Em particular, pela
   inclusão Y+ ⊂ toro, h é um vetor do núcleo na realização do toro, onde os certificados contam.
3. A classe de δ módulo gauge determina h módulo a direção de gauge. Essa direção é 0 se s ≢ μ, e é
   ℂ·e^{(μ−s)τ}g* se s ≡ μ mod (i/2)ℤ, isto é, o modo e^{μτ}g* escrito com o rótulo s. Além disso, h está
   nessa direção se e só se δ é de gauge.
4. (R7) Uma cadeia física de comprimento k dá uma cadeia de L de comprimento k módulo gauge.

**Consequência.** Com a volta (`GEOMETRIC_FLOQUET_RECONSTRUCTION.md`), a correspondência

    {modos físicos com expoente s}/gauge  ↔  (ker L(s) ∩ ker C(s)) / (direção de gauge)

é bijetiva para todo Re s > 0. Para cadeias, a correspondência é com as cadeias de L que satisfazem as
constraints de jato (C(s0)h0 = 0, C(s0)h_j + C1h_{j−1} = 0), e preserva o comprimento
(`GEOMETRIC_FLOQUET_RECONSTRUCTION.md`, atualização de 03/10). A contagem física de T2 é, portanto,
a contagem de modos físicos de Floquet, sob H-reg. H-cone virou o Lema C (`GAUGE_GLOBAL.md` §6).

## 2. Prova

### R1. Transporte nulo, fora da ressonância (`GLOBAL_NULL_GAUGE_TRANSPORT.md`)

- **As equações.** Para remover δg_{±±}, resolve-se ∂_{u_+}X^{u_-} = −δg_{++}/(2g_{+-}) e
  ∂_{u_-}X^{u_+} = −δg_{--}/(2g_{+-}). Com X^± = e^{(s−μ)τ}y_± e
  f_- := −μe^{−sτ}δg_{++}/g_{+-}, f_+ := −μe^{−sτ}δg_{--}/g_{+-}, isso vira T_1(s)y_- = f_- e
  T_{−1}(s)y_+ = f_+, com T_c(s) = ∂_τ + μ(ξ − c)∂_ξ + s − μ.
- **Periodicidade.** As duas componentes escalam do mesmo jeito sob Θ²:
  δg_{++}∘Θ² = e^{8πμ}e^{4πs}e^{−4K}δg_{++} e g_{+-}∘Θ² = e^{8πμ}e^{−4K}g_{+-}. Logo f_∓ são
  4π-periódicas. (Revisão, R0.)
- **Solução.** Fora de s ∈ μ − (i/2)ℤ, com Re s > 0, a solução periódica analítica é única. É dada por
  (*): traço no cone mais integral característica, com caminho dentro de um domínio convexo. Converge
  mesmo para 0 < Re s < μ, porque o integrando subtraído se anula em c.
- **Centro.** Pela regularidade cartesiana de δg, f_- e Pf_+ colam analiticamente em ξ = 0. A reflexão dá
  y_+ = Py_-, com A par e B ímpar, isto é, X regular no eixo.
- **Revisão:** derivação de T_{±1} em sympy (razão exata 1/(2μ)); teste numérico próprio com fonte
  analítica não polinomial, com resíduo de 10⁻¹¹ a 10⁻⁸, inclusive na ressonância
  (`r1_transporte_numerico.py`); `test_global_null_gauge` 7/7.

### R2. A ressonância s ≡ μ: fechada pela simplicidade de μ* (S4)

Seja s0 = μ − im0/2 e κ := (f_-(·, 1))_{m0}.
1. **Necessidade e solução secular.**
   - Em ξ = 1, a equação é (∂_τ − im0/2)y(·, 1) = f_-(·, 1). Para y periódico é preciso κ = 0.
   - Se κ ≠ 0, como T_1(s0)(τe^{im0τ/2}) = e^{im0τ/2}, vale y_- = y_reg + κτe^{im0τ/2}, com y_reg dado por
     (*) no complemento.
   - Como f_+(·, −1) = f_-(·, 1), a mesma constante vale para y_+.
2. **Parte secular de X.** X = X_reg + κτX0, com X0 = ∂_{u_-} + ∂_{u_+} = e^{μτ}(μ⁻¹∂_τ + ξ∂_ξ), regular
   no eixo. A mesma constante nas duas componentes é necessária: κ(∂_{u_-} + c∂_{u_+}) só é regular no eixo
   se c = 1.
3. **δ' = A + τB.**
   - L_{τX0}ḡ = τL_{X0}ḡ + 2dτ⊙X0^♭. Ponha B := κL_{X0}(ḡ, φ*) e A := δ + L_{X_reg}(ḡ, φ*) + 2κdτ⊙X0^♭.
   - B está no ansatz, porque (L_{X0}ḡ)_{±±} = 0. Logo A_{±±} = δ'_{±±} − τB_{±±} = 0, e A também está.
   - A é Floquet, porque (Θ²)*X0 = e^{4πμ}X0 = e^{4πs0}X0.
4. **Forma de δω'.** A extração de R3 é pontual e linear, e δω é de primeira ordem, com D^σ(τ) = σ. Logo
   δω[A + τB] = δω[A] + τδω[B] + σ·(variáveis de B). Com δω[L_{X0}(ḡ, φ*)] = e^{μτ}(Z + 1)ω*,
   δω' = e^{s0τ}(h1 + τh0), com h0 = κe^{im0τ/2}g* ≠ 0 e h1 4π-periódico analítico.
5. **Cadeia.**
   - δ' resolve as equações, porque o gauge preserva soluções e R4 é local.
   - Pelo L-afim, S δΩ^+[τψ] = τ S δΩ^+[ψ] + ψ.
   - Logo 0 = e^{s0τ}[L(s0)h1 + h0 + τL(s0)h0]; tomando τ ↦ τ + 4π, L(s0)h0 = 0 e L(s0)h1 = −h0.
6. **Reetiquetamento e setor.**
   - L(s)(e^{iατ}h) = e^{iατ}L(s + iα)h. Com α = −m0/2: L(μ)h1' = −κg*, com h1' = e^{−im0τ/2}h1
     4π-periódico. É uma cadeia em s = μ exato.
   - T = σ·(translação de 2π), com σ = diag(1, 1, 1, −1), comuta com L, e g* ∈ A. Projetando com
     P_A = (1 + T)/2: L_A(μ)P_A h1' = −κg*.
7. **Contradição.**
   - Por R5 (cadeias, em s = μ), P_A h1' ∈ D(Y+) ⊂ domínio do toro, pela inclusão ‖x‖_toro <= ‖x‖_+/243,
     com o mesmo Q = J_μ⁻¹.
   - S4 certificou no toro ℓ² com sinal, setor A, que o NK bordejado F(y, λ) = (L(λ)Q0y, ⟨h0, Q0y⟩ − 1)
     tem zero único (g_n, μ*) com DF(x*) invertível. O raio é 6,06·10⁻⁵, com unicidade até 3,26·10⁻⁴, e
     μ* = μ, porque L(μ)g* = 0.
   - Uma cadeia L(μ)h1 = −g_n dá o vetor (Q0⁻¹h1 − ⟨h0, h1⟩y*, 1) no núcleo de DF(x*). Isso contradiz a
     invertibilidade.
   - **Logo κ = 0.** Na ressonância, a fixação de gauge é Floquet e única módulo κ'X0, e R1 vale.
- **Revisão:** sympy com funções genéricas (`r2_gauge_x0_e_cadeia.py`). Confere a forma de δω',
  δω[L_{X0}] = e^{μτ}(Z + 1)ω nas 7 funções, a ausência de termos em τ que escapem e (L_{X0}ḡ)_{±±} = 0.

### R3. Do ansatz às variáveis de RT (prova completa da revisão)

- **Fundo.** a = −2e^{−2ζ}/(μ²(u_+ + u_-)²) e b = e^{−2ζ}ξ²/W², com W = μ + ξ²Q. A parte 2D é
  ḡ_2D = e^{−2ζ}[−(1 − ξ²)dτ² − (2ξ/μ)dτ⊙dξ + dξ²/μ²].
- **Perturbação no ansatz.** δg = (δa/a)ḡ_2D + δb·g_{S²}. Daí δζ = −δa/(2a), δW = −W(δζ + δb/(2b)),
  δQ = δW/ξ² e δφ' = δφ + X(φ*).
- **Eixo.** Uma 2-forma simétrica O(3)-invariante analítica é α dτ² + 2βdτ(x·dx) + γ(x·dx)² + d|dx|², com
  α, β, γ, d analíticas em (τ, ξ²). Então:
  - δb = dξ²;
  - o coeficiente de dξ² dá δζ = −(μ²e^{2ζ}/2)(γξ² + d);
  - δζ + δb/(2b) = ξ²(e^{2ζ}/2)[d(2μQ + ξ²Q²) − μ²γ].

  Logo δQ = −W(e^{2ζ}/2)[d(2μQ + ξ²Q²) − μ²γ], analítica e par. As componentes dτ² e dτdξ não dão
  condição extra, porque δg_{±±} = 0 já força a parte 2D a ser conforme.
- **Codificação.** Para F par, P(D^-F) = −D^+F. Logo a codificação δω_i^+ = −Pδω_i^- é automática.

δω é Floquet com expoente s e analítico numa vizinhança do cilindro real |ξ| <= 1 (até |ξ| < 41/40).

### R4. As equações (`GEOMETRIC_FLOQUET_RECONSTRUCTION.md`, parte recíproca)

δ' está no ansatz e resolve as equações linearizadas. As definições dão, linearizadamente,
Ω1^σ = Ω2^σ + Ω2♯^σ e a igualdade entre os sinais para Ω2, Ω3 e Ω4. As fórmulas de curvatura
(dfkhdjhsdshkfd) de RT forçam os cinco resíduos restantes a zero fora do centro, e a analiticidade
estende isso ao centro. Logo δΩ^σ = 0 para todos: L(s)h = 0 e C(s)h = 0. O argumento é local e vale
também para perturbações seculares (R2, R7).

**Revisão.**
- RT enunciam as fórmulas de curvatura sem derivá-las. A revisão calculou o Ricci 4D direto da métrica,
  com (ζ, Q, φ) genéricos, em 3 jatos racionais, e as fórmulas e as 5 identidades de definição batem
  exatamente.
- O controle negativo detecta as 4 alterações introduzidas (`r4_curvatura_rt.py`, `r4_controle_negativo.py`,
  que precisam de sympy).
- `test_geometric_reconstruction` 5/5; ele inverte a matriz das cinco combinações transcrita.

### R5. Espaço espectral: recuperação de raio de qualquer peso > 1, para qualquer s

- **Pesos.** h está em C^∞(ℝ/4πℤ; A(E)).
  - Os coeficientes de Fourier h_m decaem mais rápido que qualquer potência de |m|, uniformemente em E.
  - Os coeficientes de Chebyshev decaem como ρ^{−n}, com ρ o parâmetro de uma elipse de Bernstein
    contida em E.

  Logo h ∈ Y_{(1, κ2')} para todo 1 < κ2' < ρ: peso **1** em Fourier, sem crescimento exponencial. Além
  disso, J_μh = −(B + s)h ∈ Y_{(1, κ2')}. O mesmo vale para os perfis de uma cadeia.
- **κ1' = 1 é permitido.**
  - Q = J_μ⁻¹ preserva m, então ‖QT‖ é o supremo sobre os modos, e free_tail não depende do peso de
    Fourier.
  - A norma de B com peso de Fourier 1 é menor ou igual à norma com peso κ1 > 1, porque a convolução pelos
    coeficientes do fundo tem norma <= Σ_k|b_k| <= Σ_k κ1^{|k|}|b_k|.

  Logo as cotas abaixo valem com κ1' = 1, sem mudança.
- **Lema R5.** Para todo κ' = (κ1', κ2') com 1 <= κ1' <= 65/64 e 1 < κ2' <= 9/8, e **todo S > 0**, ker L(s) e as cadeias de Jordan em D(Y_{κ'}) estão
  em D(Y+), para |s| <= S.
- **Prova.** É o argumento de `RADIUS_RECOVERY.md` (existência forte e unicidade fraca, com cadeias por
  indução), com Q = J_μ⁻¹ e H(s) = I + (B + s)Q. Basta uma caixa G com contração nos dois espaços.
  - **Fraco.** (β(κ2') + S)·free_tail(M, N, κ2'), com β(κ2') = β+·D(κ2')/D(5/4).
    - D(k) = 2k/(k² − 1) é 𝒦1♯ de RT, nas estimativas auxiliares do sistema ♯.
    - A reponderação vale porque a parte linear de β+ é proporcional a 40/9.
    - free_tail não depende de κ1'.
  - **Cota explícita.** free_tail <= max((3/2)(1 + 2μ_max)/(b(1 − r)), (3/2)(1 + 4/(1 − r))/(μ_min(N + 1))),
    com r = 1/κ2' e b = M/2. Ela tende a 0, logo para todo κ2' > 1 e todo S existe G.
  - **Forte.** (β+ + S)·structured_input_tail(M, N) decresce com G.
- **Artefato.** `scripts/radius_recovery_geral.py` (`build/t2/radius_recovery_geral.json`) confere
  free_tail <= cota e dá caixas explícitas para S ∈ {3,1; 10; 50; 200} e κ2' até 1 + 2⁻¹². Exemplos:
  S = 200, κ2' = 1 + 2⁻⁶ → G = 2¹⁷×2²⁰.
- **Revisão:** free_tail recalculado por um caminho independente em 4 pesos e 3 caixas; free_tail <= cota
  provado em papel e conferido em 4000 amostras.
- **Uso de F3.** Para um modo físico com Re s >= 2,93, h ∈ D(Y+) ⊂ toro está no núcleo na realização em
  que F3 exclui raízes. Logo esse modo não existe.

### R6. Unicidade módulo gauge

- **δω' = 0.** As equações (dhjhfkfhfe1–4) de RT dão (δζ, δQ, δφ) = (c1, 0, c2), constantes. Um perfil
  Floquet com Re s > 0 não é constante, logo δ' = 0 e δ é de gauge.
- **δω' na direção de gauge.** Se δω' = c·e^{(μ−s)τ}g*, isto é, c vezes o modo e^{μτ}g*, então
  δ + L_{X − cX0}(ḡ, φ*) tem δω = 0 e cai no caso anterior.
- **Dois gauges.** Dois gauges que levam δ ao ansatz diferem por um campo admissível da classe declarada.
  Pelo Lema G ou G′, os δω diferem por um elemento da direção de gauge.

### R7. Cadeias físicas (lacuna L9 da revisão)

Seja δ um modo generalizado de ordem k, com perfis δ_j. O transporte de jatos de
`GLOBAL_NULL_GAUGE_TRANSPORT.md` resolve T_c(s)y_j = f_j − y_{j−1}, isto é, gauges polinomiais em τ.
- **Fora da ressonância:** cada ordem é solúvel.
- **Na ressonância:** a compatibilidade de cada ordem é forçada como em R2. Um termo secular a mais
  alongaria a cadeia de L em μ, e μ* é simples.

O δ' resultante está no ansatz e, por R3–R5, δω' = e^{sτ}Σ_j τ^j h_j, com (h_0, ..., h_{k−1}) uma cadeia
de L(s) em D(Y+).
- Se a cadeia de L tem comprimento menor, por exemplo h_0 = 0, então, por R6 aplicado à ordem mais baixa,
  δ_0 é de gauge, e a cadeia física encurta módulo gauge.
- Como as quatro raízes de L em Re s > 0 são simples, **toda cadeia física tem comprimento 1 módulo
  gauge**. A multiplicidade física é a multiplicidade de L módulo gauge.

## 3. O que fica

- **H-reg′ (modelagem):** a classe física é C^∞ em τ e real-analítica em ξ, com vizinhança complexa
  uniforme.
  - A analiticidade em τ não é mais exigida, porque R5 vale com peso de Fourier 1.
  - As outras peças preservam C^∞ em τ:
    - R1, porque o inverso temporal é o multiplicador 1/(s − μ + im/2) e a integral característica
      comuta com ∂_τ;
    - R2 e R7, que têm a mesma estrutura;
    - R3 e R4, que são locais;
    - R6.
  - **Fica fora a regularidade só C^∞ em ξ**, através do cone e no eixo. Provar que um modo C^∞ com
    Re s > 0 é analítico em ξ é um problema de regularidade aberto. O sistema é hiperbólico, periódico em
    τ e acoplado. A cota livre radial piora como 1/(κ2' − 1), então a máquina de R5 não alcança κ2' = 1.
  - A afirmação "Floquet módulo gauge implica Floquet num gauge analítico" não está provada; exigiria
    resolver (Θ²)*Z − e^{4πs}Z = Y. Por isso a definição pede Floquet exato.
- **H-cone:** não é mais hipótese. Pelo Lema C (`GAUGE_GLOBAL.md` §6), a direção de gauge é tangente às
  translações de T*; as duas contagens, 1 e 2, são teoremas.
- **Natureza espectral:** o lema trata de modos de Floquet. Não trata da completude dos modos para o
  problema de evolução.
- **A volta** (`GEOMETRIC_FLOQUET_RECONSTRUCTION.md`) foi revisada em 03/10 (`.codex-runs/2026-10-03-revisao-V/`).
  Os certificados citados (S4, F3, L-real) têm revisões próprias.
