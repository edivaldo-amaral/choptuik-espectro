# Revisão independente leve: Lema C, eixo imaginário e H-reg′ (03/10/2026)

Revisor: sub-agente de contexto limpo (Claude Fable 5.1). Li só os arquivos permitidos em `TAREFA.md`. De
`docs/C1_REAVALIACAO.md` li só as §6, §9, §11 (com a 11.1) e §13. Usei só o notebook, sem ssh, com
`.venv/bin/python`, `OPENBLAS_NUM_THREADS=2` e sympy 1.14 instalado com `pip --target` no scratchpad. O
artigo de RT foi extraído no scratchpad. Não alterei nenhum arquivo existente.

## Manifesto, no início

- `.snapshots/2026-10-03-revisao-CRH/MANIFEST.sha256` (sha256 e44d2816…a04b): **240/240 OK**, conferido
  com `LC_ALL=C sha256sum -c`.
- `COMMIT.txt` = 529a19696d3b790d1d97819acdd415b33a1edb02 = HEAD no início.

## Veredictos

| item | veredicto | artefato |
|---|---|---|
| **C1** Lema C | **OK com ressalva** | `c1_lema_c.py` → `c1_saida.txt`: 20/20 conferências simbólicas com sympy. (a) ψ_ε em (τ, ξ) preserva du_± e o eixo. (b) Fator (u_+ + u_-)². (c) ψ_ε*ḡ está na forma de RT com ζ_ε, Q_ε e W_ε = W∘ψ_ε. (d) d/dε ω_ε = e^{μτ}(Z + 1)ω nas sete funções, além de δζ e δQ. (e) X0. Também rodei `test_gauge_action` (3/3 OK). |
| **C2** eixo imaginário | **OK com ressalva** para a afirmação sobre L. **LACUNA (média)** na leitura física ("não há modos neutros físicos") | `c2_eixo.py` → `c2_saida.txt`: cobertura racional do segmento fechado, discos que contêm 0, i/4, i/2, ±… e cota de q. `c2_formula.py` → `c2_formula_saida.txt`: a fórmula de deflação em dimensão finita, com casos de Jordan e de raiz dupla. `c2_neutros.py` → `c2_neutros_saida.txt`: os modos de escala e de translação de φ. |
| **C3** H-reg′ (κ1′ = 1) | **OK** (ressalva trivial) | `c3_hreg.py` → `c3_saida.txt`: assinatura de free_tail; majorantes de B recalculados com peso de Fourier 1 contra o original; convolução com pesos; perfil C^∞ não analítico em τ. As derivações estão abaixo. |

## C1: Lema C

**Derivação** (conferida em `c1_lema_c.py`).
- u_σ = −(1 + σξ)e^{−μτ}. Em (τ, ξ), ψ_ε é
  - T = −μ⁻¹log(e^{−μτ} − ε);
  - Ξ = ξe^{−μτ}/(e^{−μτ} − ε).

  Valem u_σ(T, Ξ) = u_σ + ε, e Ξ = 0 se e só se ξ = 0. Logo ψ_ε preserva du_± e o eixo. **[OK]**
- A parte radial de RT é e^{−2ζ}(−dτ² + (dξ − μξdτ)²/μ²) = −2e^{−2ζ}(du_-du_+ + du_+du_-)/(μ²(u_+ + u_-)²).
  **O fator (u_+ + u_-)² confere.** **[OK]**
- Com S = u_+ + u_-:
  - e^{−2ζ_ε} = e^{−2ζ∘ψ_ε}S²/(S + 2ε)², exatamente a fórmula do §6;
  - a parte radial de ψ_ε*ḡ é a forma de RT com ζ_ε;
  - **W_ε = W∘ψ_ε exatamente.** O fator S/(S + 2ε) de e^{−ζ_ε} cancela o de Ξ/ξ;
  - Q_ε = Q∘ψ_ε·(Ξ/ξ)², com Ξ/ξ = e^{−μτ}/(e^{−μτ} − ε), que não depende de ξ.

  Logo Q_ε é par e analítica no eixo sempre que Q for, e W_ε > 0 vale em todo o domínio transladado, não
  só para ε pequeno. A regularidade no eixo é imediata e não precisa de R3. **[OK]**
- **d/dε ω_ε|_{ε=0} = e^{μτ}(Z + 1)ω nas sete funções** ω1, ω2^±, ω3^± e ω4^±, pelas definições
  (ekjehee) de RT.
  - A identidade depende só do 2-jato de (ζ, Q, φ) no ponto, e foi conferida com polinômios de grau 3 de
    coeficientes simbólicos, cujo 2-jato num ponto é arbitrário.
  - Na primeira tentativa, com funções genéricas, o sympy deixou objetos `Subs` com ε dentro. Era um
    artefato da ferramenta, não da identidade.
  - Também conferem δζ = e^{μτ}(Zζ − 1), δQ = e^{μτ}(ZQ + 2Q) e X0 = e^{μτ}(μ⁻¹∂_τ + ξ∂_ξ). **[OK]**
- **X0 é C¹ (de fato analítico) em M̂.** Em (τ, x) cartesianas, X0 = e^{μτ}(μ⁻¹∂_τ + x·∂_x), analítico em
  ℝ × ℝ³, inclusive no eixo e no cone. X0 é admissível, porque ψ_ε fica no ansatz, e é Floquet:
  (Θ²)*X0 = e^{4πμ}X0. **[OK]**
- **"Contagem 1 com o gauge padrão" segue**, dados o Lema G e S4 (E = span(g*), simples). e^{μτ}g* = δω_{X0}
  é de gauge, e dim E − rank(D|E) − 1 = 0 em μ*. **[OK]**
- **"Contagem 2 com cone fixo" segue no ansatz.** Com gauges admissíveis com f(0) = 0, o Lema G não deixa
  nenhum modo de gauge em Re s > 0, e g* (com C(μ)g* = 0) conta. Ver a ressalva R-C1b.

**Ressalvas de C1.**
- **R-C1a (baixa, rastreabilidade).**
  - O Passo 8 e o §6(3) de `GAUGE_GLOBAL.md` citam `test_gauge_action.py` para a ação de X0 nas variáveis
    de RT.
  - O teste confere outra identidade: a covariância **fora da solução das equações**,
    DΩ_eq(w)[(Z + 1)w] = (Z + 1)Ω_eq(w), via `linearized_omega` e `omega` de `sharp_propagation`.
  - Não confere δω_{X0} = DΩ_var[L_{X0}(ḡ, φ*)] = e^{μτ}(Z + 1)ω, que é a mudança de variáveis.
  - Esta última está agora conferida por `c1_lema_c.py`.
- **R-C1b (baixa, definição).** "Gauges que fixam o cone" só está definido para X admissível, via f(0) = 0.
  - No nível de perturbações gerais, a definição ingênua não funciona: "módulo X tangente ao cone" dá um
    quociente de dimensão infinita em cada expoente, porque X^{u_-}|_cone = e^{(s−μ)τ}Y^ξ(τ, 1) é livre.
  - A definição consistente pede perturbações que mantêm u_- = 0 nulo (δg_{++}|_cone = 0) e gauges com
    X^{u_-}|_cone = 0.
  - Com ela, o transporte de R1 no cone, (∂_τ + s − μ)y_-(·, 1) = 0, dá y_-(·, 1) = 0 para s ≢ μ. Só em
    s ≡ μ sobra κ′X0, e a contagem é 2.
  - Convém escrever essa definição no enunciado.
- **R-C1c (trivial).** Para ε < 0, parte de ψ_ε⁻¹(M) tem u_+ + u_- >= 0 e sai da carta (τ, ξ). O lema
  vale em compactos de M, como está escrito.

## C2: eixo imaginário

**Cobertura** (`c2_eixo.py`, racionais exatos a partir dos floats gravados, que são os valores usados no
Perron racional).
- 175 entradas lidas: 96 da faixa A (`discos_v2`) e 79 da faixa B.
  - Todas têm `ok` e θ racional < 1.
  - Uma tem r = 0 (0,291688 + 0,75i), a L1 já conhecida, e foi descartada.
  - Ficam **174 discos válidos**.
- A meia-corda de cada disco na reta Re s = 0 foi cotada por baixo em racionais (h_lo² <= r² − a²).
  - **A união cobre o segmento fechado Re s = 0, Im s ∈ [−1/4, 3/4]**, com 76 discos que tocam o eixo.
  - Uma cadeia gulosa de 69 cordas tem sobreposição mínima de 1,06·10⁻⁵ e vai de −0,2656 a 0,7517.
  - **[OK]**
- Cantos e junção:
  - −i/4: folga 0,0156;
  - i/4: 4 discos, folga 0,0156;
  - i/2: folga 0,0004, num disco de r = 0,0028 junto da raiz B (0,0414 + 0,5i);
  - 3i/4: folga 0,0097.
- Validade de um disco no eixo. A docstring de `rouche_L_disco.py` diz que o Perron vale para todo s com
  |s − p| <= r e Re s >= 0, onde valem as formas fechadas de Q0. Os pontos do eixo têm Re s = 0. **[OK]**
- **Não há disco centrado em s = 0.** `discos_v2/disco_0.json` é o arquivo do **centro de cauda** c = 0.
  - s = 0 está no interior de dois discos:
    - p = 0,004714i, r = 0,022676, folga 0,01796, θ = 0,99216;
    - p = −0,017065i, r = 0,022210, folga 0,00515, θ = 0,99737.
  - Isso basta para L̃(0) invertível. Ver L2.
- Conferência de sanidade na referência finita deflacionada (`diagT_lab2.npy`):
  - nenhum autovalor no eixo com Im ∈ [−1/4, 3/4];
  - os transladados não deflacionados ±i, ±2i, … aparecem, como esperado;
  - em −1/4 < Im <= 3/4 e Re > −0,3 ficam só 0,0414 + 0,5i, 0,401, 0,733 e raízes com Re < 0.

**A fórmula mult_L̃(s) = mult_L(s) + ord_s(q/(s(s − μ))).**
- Vale pelo teorema de Gohberg–Sigal (resíduo logarítmico, aditividade M(AB) = M(A) + M(B) para funções
  operador finitamente meromorfas de Fredholm), aplicado a L̃ = L(I + L⁻¹E).
- Para I + F com F de posto finito, M(I + F) = ord det(I + F).
- **Hipóteses**, todas satisfeitas aqui:
  1. L(s) = L(0) + sI é uma família analítica de Fredholm de índice 0 numa vizinhança conexa do segmento.
     Pelo `RADIUS_RECOVERY.md`, H(s) = I + (B + s)Q com Q compacto, para todo s. Q0(s) é invertível e
     analítica em Re s > −μ, então as multiplicidades de L e de H̃ = L̃Q0 coincidem.
  2. L(s0) é invertível em algum ponto. Vale por F3 e Neumann em s real grande, ou num ponto do contorno
     onde L̃ é invertível e det ≠ 0, ∞. Então L⁻¹ é finitamente meromorfa.
  3. E = Σ δ_j x_j φ_j* tem posto 2, é constante e limitada.
  4. Os x_j são **autovetores exatos** no domínio, com L(s)x_j = (s − λ_j)x_j. Só assim
     det(I_2 + [φ_i*L⁻¹δ_j x_j]) = q(s)/(s(s − μ)) exatamente. Com φ_i*x_j = δ_ij + e_ij, sai
     q = (s + δ0(1 + e00))(s − μ + δ1(1 + e11)) − δ0δ1e01e10, que confere com o §9.
- **Os vetores são os exatos?** Sim, no operador certificado.
  - O código monta a deflação com x_j^A (RefA) e soma a diferença δ_k(x_k − P_Z x_k^A)φ_k* como
    perturbação: dEZ nas linhas de Z e dET nas de T e do far, em `rouche_L_disco.py` (L1 da revisão de 01/10).
  - Os erros usados são ‖∂_τRefB‖ = 1,4563·10⁻⁹ (mais 10⁻⁶⁰ para Corr), e = 2,5056·10⁻⁸ (`build/s4/rig/E.json`)
    e a cauda `cauda_x` (`schur.json`).
  - Esses números são herdados de S4 e da revisão de S4; não os refiz.
  - `c2_formula.py` mostra, em dimensão finita, por que a exatidão importa: com vetor de fase inexato
    (|η| = 10⁻⁶), det(L̃)/det(L)·s(s − μ) − q deixa de se anular perto de 0 (3,5·10⁻⁷).
- **Teste da lógica** (`c2_formula.py`):
  - a identidade det L̃/det L = q/(s(s − μ)) confere com erro de 9·10⁻¹⁵;
  - as raízes de L̃ são as de L, menos 0 e μ, mais as duas de q (≈ −1);
  - com 0 simples, σ_min(L̃(0)) = 0,18;
  - com 0 duplo (bloco de Jordan ou semissimples), L̃(0) é singular (σ_min ~ 10⁻¹⁵).

  Logo "L̃(0) invertível" certifica mult_L(0) = 1. **[OK]**
- **q não se anula em Re s >= 0.**
  - |e_ij| <= max‖φ_i‖·max‖P_Z(x_j − x_j^A)‖ <= 5,5013·2,5056·10⁻⁸·1,001 = 1,38·10⁻⁷.
  - Os fatores satisfazem |s + 1 + e00| >= 1 − 1,38·10⁻⁷ e |s + 1 + (μ_A − μ*) + δ1e11| >= 1 − 3,9·10⁻¹¹ − 1,17·1,38·10⁻⁷.
  - Então |q| >= 0,9999997.
  - q(0) ≈ 1 ≠ 0, logo o polo de det(I + L⁻¹E) em 0 é simples. **[OK]**
- **Período e setores.**
  - L(s)(e^{iατ}h) = e^{iατ}L(s + iα)h, porque L(s) = L(0) + s vem de ∂_τ → ∂_τ + s.
  - Com α = 1 (m → m + 2) o setor é preservado (T = σ·translação de 2π comuta). A multiplicação por e^{iτ}
    é um isomorfismo do espaço com pesos κ1^{|m|} (fator <= κ1²) e do domínio. Logo Σ_A e as
    multiplicidades são i-periódicos.
  - Com Σ_B = Σ_A − i/2: Σ_full ∩ iℝ = iℤ ∪ i(ℤ + 1/2) = (i/2)ℤ ≡ 0 mod i/2.
  - Em 0: mult_full(0) = mult_A(0) + mult_B(0) = 1 + mult_A(i/2) = 1 + 0. i/2 está coberto, folga 4·10⁻⁴.
  - Com a convenção oposta (Σ_A + i/2), o resultado é o mesmo, por i-periodicidade. **[OK]**
- **Furos procurados.**
  - Cantos: cobertos.
  - s = 0 sobre o contorno de Ω_A: para o §13 não importa, porque só se usa invertibilidade pontual. Para
    a contagem do item 1, o recuo por meio-círculo do §6 está correto: L̃ é invertível num disco em torno
    de 0, e L só tem a raiz 0 ali.
  - Convenção dos setores: indiferente no eixo.

**Conclusão sobre L.** A única raiz de L (período completo) em Re s = 0 é s ≡ 0 mod i/2, algebricamente
simples, com autovetor ∂_τω*, que é de gauge (n = 1). **Confere.**

**A leitura física não confere: LACUNA L1 abaixo.**

## C3: H-reg′

- **free_tail não depende do peso de Fourier.**
  - A assinatura é `free_tail(m, n, radial, mu_min, mu_max)`, sem κ1. A fórmula só usa M, N, κ2 e μ.
  - Motivo: Q = J_μ⁻¹ e T = I − G são bloco-diagonais em m, porque o principal tem coeficientes que não
    dependem de τ.
  - Na norma ℓ^p em m de normas radiais, com pesos escalares κ1^{|m|}, a norma de um operador bloco-diagonal
    é sup_m ‖(QT)_m‖_rad, que não depende de κ1. **[OK]**
- **A norma de B com peso 1 é <= a com peso κ1 > 1.**
  - B não contém ∂_τ (está em J). É multiplicação pelo fundo, em convolução em m, mais operadores radiais.
  - Young com peso submultiplicativo: ‖b * x‖_κ <= Σ_k κ^{|k|}|b_k|·‖x‖_κ, e a cota é monótona em κ.
    Numericamente (c3 §3), as normas em ℓ² com peso são 3,754 (κ = 1), 3,822 (65/64) e 4,853 (1,2), todas
    abaixo das cotas.
  - Recalculei os majorantes **do código** com `fourier = 1`:
    - `component_majorant` (65/64 → 1): todas as entradas diminuem, razão máxima 0,9904, cota de coluna
      19,688 <= 20,273;
    - `sharp_majorant` (129/128 → 1): razão máxima 0,9966, cota 4,552 <= 4,622.
  - Todos os índices m de RefA são >= 0 (m ∈ [0, 40]). Se houvesse m < 0, fourier^m < 1 quebraria a
    monotonia. **[OK]**
- **Um perfil C^∞(ℝ/4πℤ; A(E)) está em Y_{(1, κ2′)}, junto com J_μh.**
  - Em m: ‖h_m‖_{A(E)} <= ‖∂_τ^k h‖_∞(2/|m|)^k para todo k, decaimento superpolinomial uniforme em E.
  - Em n: f(cos θ), com f analítica numa elipse de Bernstein E_ρ ⊂ E, tem coeficientes <= 2‖f‖ρ^{−|n|}.
  - Daí Σ_{m,n} κ2′^{|n|}|h_mn|·(pesos polinomiais) < ∞ para κ2′ < ρ.
  - J_μ multiplica por im/2 + μ(n + 1) (mais a parte triangular em n), que é polinomial. Logo
    J_μh ∈ Y_{(1, κ2′)}, ou, de modo equivalente, ∂_τh e ξ∂_ξh estão em C^∞(A(E′)) com E′ ⊂ E, por Cauchy.
  - Exemplo (c3 §4): h = exp(−1/sin²(τ/4))/(2 − ξ).
    - Em m, |g_m| ~ exp(−1,7 m^{2/3}) (Gevrey 3/2). Isso é superpolinomial e **não geométrico**: h está em
      Y_{(1, κ2′)} e **não** em Y_{(κ1, ·)} para nenhum κ1 > 1.
    - Em n, o decaimento é exatamente ρ^{−n}, com ρ = 2 + √3.

    É exatamente o caso que κ1′ = 1 passa a admitir. **[OK]**
- **R5 com κ1′ = 1.** A inclusão Y+ ⊂ Y_{(1, κ2′)} é contínua, e G é finito. A cota fraca
  (β(κ2′) + S)·free_tail não usa κ1′, e β(κ2′) só reescala a divisão por ξ (D(κ2′)/D(5/4)). Existência
  forte e unicidade fraca valem sem mudança. **[OK]**
- **As peças preservam C^∞ em τ.**
  - **R1.** T_c(s) = ∂_τ + μ(ξ − c)∂_ξ + s − μ tem coeficientes que não dependem de τ. Por modo:
    - a_m = im/2 + s − μ, e y_m(c) = f_m(c)/a_m (o traço no cone);
    - z_m(ξ) = μ⁻¹∫_0^1 t^{a_m/μ − 1}g_m(c + t(ξ − c))dt, com g_m = f_m − f_m(c).

    Daí |z_m| <= |ξ − c|·sup|g_m′|/Re s e |y_m(c)| <= |f_m(c)|/dist(s, μ − (i/2)ℤ). Então
    ‖y_m‖_{A(E′)} <= C(s)‖f_m‖_{A(E)}, **uniformemente em m**, com E convexo (estrelado em ±1). O
    decaimento rápido em m, isto é, C^∞ em τ, se preserva. A conta usa Re s > 0: com Re s = 0 a integral
    diverge, e há ressonância em a_m + μ = 0, isto é, s = −im/2. Isso é relevante para L1.
  - **R2 e R7.** Mesma estrutura, com termo secular polinomial em τ.
  - **R3.** Algébrica e local. A divisão por ξ² de funções pares que se anulam a ordem 2 preserva
    A(E) uniformemente em τ.
  - **R4.** Local. Estender δΩ = 0 de ξ ≠ 0 ao centro só precisa de continuidade.
  - **R6.** Usa as equações (dhjhfkfhfe1–4), sem regularidade extra.
  - Os gauges resultantes são C^∞ em τ e analíticos em ξ, logo C¹, e estão na classe do Lema G. **[OK]**
- **Região 1 < ξ < 41/40 fora de E.** Se E não chega a 41/40, o perfil ali é a solução do problema de
  Cauchy a partir de ξ = 1 + δ: g^{−1}(dξ, dξ) ∝ μ²(1 − ξ²) < 0, então ξ é temporal além do cone. h = 0
  em E força h = 0 até 41/40. A injetividade não tem furo.
- **Ressalva trivial.** A docstring de `scripts/radius_recovery_geral.py` ainda diz "1 < k1′ <= 65/64". O
  código não usa κ1′, então nada muda.

## Lacunas

| # | gravidade | onde | lacuna | correção sugerida |
|---|---|---|---|---|
| **L1** | **média** (enunciado; não afeta os itens 1–2 nem a contagem instável) | `C1_REAVALIACAO.md` §13 ("Não há modos neutros físicos"); título do item 3 de `T2_ENUNCIADO.md` ("sem modos neutros além da fase") | Ver os quatro pontos abaixo da tabela. | Restringir o item 3 a L ("a única raiz de L em Re s = 0 é a fase, simples"). Se quiser a leitura física: "módulo gauge **e módulo as simetrias de escala e de translação de φ**, invisíveis a L", condicionada a estender o Lema R a Re s = 0, com a ressonância de fase tratada à parte. |
| L2 | baixa (rastreabilidade) | §13 e item 3 de T2: "há um disco centrado em s = 0 (`disco_0.json`)" | Não há disco com p = 0. `disco_0.json` é o arquivo do centro de cauda c = 0. | Citar os dois discos que contêm 0: p = 0,004714i, r = 0,022676, folga 0,018; e p = −0,017065i, r = 0,022210, folga 0,005. A conclusão não muda. |
| L3 | baixa (redação) | §13: "cada ponto está num disco em que t·ε_B < 1, logo H̃(s) é invertível" | Atribuição errada do mecanismo. No disco, H̃ é invertível por t_p (Loewner: A(p) invertível, t_s = t_p/(1 − rt_p)) e θ < 1 (Perron racional). t·ε_B < 1 é a homotopia A → B da contagem finita. | Trocar por "t_p certificado e θ < 1". |
| L4 | baixa (definição) | `GAUGE_GLOBAL.md` §6, item 2 de T2 | "Gauges que fixam o cone" só está definido para X admissível (f(0) = 0). Para perturbações gerais, "módulo X tangente ao cone" dá um quociente de dimensão infinita. | Escrever a definição: perturbações com δg_{++}\|_cone = 0, módulo X com X^{u_-}\|_cone = 0. Com ela, R1 dá a contagem 2. |
| L5 | baixa (rastreabilidade) | `GAUGE_GLOBAL.md` Passo 8 e §6(3) | `test_gauge_action.py` não testa δω_{X0} = e^{μτ}(Z + 1)ω (mudança de variáveis). Testa a covariância fora da solução das equações de campo. | Citar `c1_lema_c.py` desta revisão (ou o r2 da revisão de R). |
| L6 | trivial | `scripts/radius_recovery_geral.py`, docstring | Ainda diz "1 < k1′". | Atualizar para 1 <= k1′. |

**Os pontos de L1:**
1. Com a própria definição de gauge do projeto (só difeomorfismos, §5, obs. 3, de `GAUGE_GLOBAL.md`), há
   **dois modos físicos neutros, não de gauge**, em s ≡ 0:
   - a **escala** δ = (2ḡ, 0), com δζ = −1, no setor A;
   - a **translação do escalar** δ = (0, 1), no setor B do período completo.
2. Os dois resolvem as equações linearizadas: os Christoffel e o Ricci não mudam por escala constante, e
   □_{e^{2c}g} = e^{−2c}□_g. São Floquet com s = 0, estão em H-reg′ e são invisíveis a L, porque δω = 0
   (`c2_neutros.py`).
3. Não são de gauge. Estão na forma do ansatz, logo um X com L_X = δ seria admissível, com δω_X = 0. O
   Passo 5 do Lema G, que não usa Floquet nem Re s, dá X = 0. Contradição.
4. Além disso, o Lema R (ida) e o Lema G só valem para Re s > 0:
   - a integral característica de R1 diverge em Re s = 0, com ressonância de fase em s = −im/2;
   - R6 usa Re s > 0 justamente para excluir os perfis constantes, que são esses dois modos.

**Herdado, não verificado aqui:**
- a montagem Perron de cada disco (θ < 1); conferi só os θ racionais gravados e a cobertura;
- t_p, as caudas e o certF;
- e = 2,5·10⁻⁸ (S4, etapa E) e ‖∂_τRefB‖;
- o caráter Fredholm na realização do toro (`RADIUS_RECOVERY.md`, "auditado, revisão parcial" em T2).

## Artefatos (sha256)

```
90799f3699085358648f06a1e54f6238667ae9c40a36ac746bcd9183f0867f81  c1_lema_c.py
52d6adbb7c5163ae25cad745f27fdc6367a528f7964f5a75b32ca49502d1048d  c1_saida.txt
82b08380c2afbc10b5af0729a4851dbd54148c65fce8babe31170411a0961bfe  c2_eixo.py
f9dd7b90668424ff711470297804c5d1c6f341dcf8ef73c90f014f611be8fb36  c2_saida.txt
7c2dbfc55cb5892062ee813dace03e20e951727436b983d35f55e11a1def9404  c2_formula.py
a96cbe698bc9845e0670b4296346a6121ef54c1408e16e5ef8b30fd75c581b68  c2_formula_saida.txt
ee2550ac1bfaba343126adfeb8d43300a3131e6d03d30163157ce31304ba1ce0  c2_neutros.py
738965bb69a6d76c646c2be4dd6a396bc2b2810b5a412c9b29e9d76886255b59  c2_neutros_saida.txt
346c8bdff8d06ae18c4d510d8f2a0642f7c5b46e2c556119907a440301697aaa  c3_hreg.py
efa3c94cb6e99401639ec0f0e3a31d3b750030b2be13abec2beff4a1b8dff54c  c3_saida.txt
```

**Como rodar.**
- `c2_eixo.py` e `c2_formula.py` usam só o `.venv`.
- `c1_lema_c.py` e `c2_neutros.py` precisam de sympy: `PYTHONPATH=<scratch>/pylib`.
- `c3_hreg.py` roda da raiz do projeto, porque importa `scripts/`.

## Manifesto, no fim

- `MANIFEST.sha256`: **240/240 OK** (`LC_ALL=C sha256sum -c`, rc = 0).
- **O HEAD mudou durante a revisão, sem ação minha:** de 529a196 para bd5dff0 ("Não esférico, Fases 0–1").
  - O commit só acrescenta `docs/NAO_ESFERICO.md` e `scripts/ns_axial.py`, nenhum deles revisado aqui.
  - Também apareceram `scripts/ns_polar*.py` não rastreados, de outra sessão.
  - Nenhum arquivo do manifesto mudou.
- Nenhum .pyc dos módulos que importei foi regravado: as execuções usaram `PYTHONDONTWRITEBYTECODE=1`, e
  os .pyc de `test_gauge_action` e `sharp_propagation` são de 17/09. Os únicos .pyc novos em
  `scripts/__pycache__` são `ns_axial` e `ns_polar`, da outra sessão.
