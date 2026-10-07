# Revisão independente do Lema R (H-rec, ida) — relatório

- **Revisor:** agente independente (Claude), 03/10/2026, das 11:31 às ~12:00 (−04).
- **Objeto:** `docs/HREC_IDA.md` no commit 8687cc39fc396124ec460491f271aefbad3e036b.
- **Objetivo:** tentar quebrar o lema.
- **Resultado:** **não quebrou.** Não achei erro matemático em nenhuma das seis peças. As peças de maior
  risco se sustentam:
  - R2, a ressonância resolvida pela simplicidade de μ*;
  - R5, a recuperação de raio a partir de qualquer peso;
  - R1 e R4, as provas herdadas.

  Há lacunas de escopo e de redação, todas consertáveis. A mais relevante (L1) é a diferença entre
  |s| <= 3,1 no lema e "Re s > 0" na bijeção e em T2.

## Manifesto

- **Início:** `LC_ALL=C sha256sum -c .snapshots/2026-10-03-revisao-R/MANIFEST.sha256` deu 241/241 OK e código
  de saída 0 (`manifest_inicio.txt`). HEAD = 8687cc3, igual a `COMMIT.txt`.
- **Fim:** ver a última seção.

## Tabela de veredictos

| peça | veredicto | artefato | o que foi feito |
|---|---|---|---|
| R0 enunciado e definições | **OK com ressalva** | derivação §R0 | Verifiquei o fator e^{−4K}, a periodicidade de f_∓, a classe de gauge e o item 3. Ressalvas L1, L2, L3, L5 e L9 |
| R1 transporte | **OK** | `r1_transporte_numerico.py` (+ `.log`); derivação §R1; `test_global_null_gauge` (7/7 OK) | Rederivei T_{±1}, o sinal e o fator μ (sympy). Provei que (*) resolve a equação, converge para 0 < Re s < μ e é analítica. Conferi a unicidade por Taylor, a reflexão, a paridade A/B e a colagem de f_- com Pf_+. Teste próprio com fonte não polinomial, inclusive a ressonância m0 = 0 e m0 = 1 |
| R2 ressonância | **OK com ressalva** (redação: L4, L5) | `r2_gauge_x0_e_cadeia.py` (+ `.log`); derivação §R2 | Derivei cada passo. Conferi em sympy, com funções genéricas, quatro pontos: X0, B no ansatz, δω[L_{X0}] = e^{μτ}(Z+1)ω para as 7 funções, e a forma afim em τ da extração. Conferi a cadeia, o reetiquetamento m0 ≠ 0, a projeção no setor A, o espaço de S4 (toro ℓ² com sinal) e a inclusão D(Y+) ⊂ toro |
| R3 extração e eixo | **OK** (o esboço vira prova; L6) | derivação completa §R3 | Prova completa de δζ + δb/(2b) = ξ²·(par analítica) e de δQ analítica e par. Também mostrei que a codificação de RT (ω^+ = −Pω^-, ω1 ímpar) vale automaticamente |
| R4 equações | **OK** (L7) | `r4_curvatura_rt.py` (+ `.log`), `r4_controle_negativo.py` (+ `.log`); `test_geometric_reconstruction` (5/5 OK) | Conferi de forma independente as fórmulas de curvatura de RT (dfkhdjhsdshkfd) e as 5 identidades de definição. Usei Ricci 4D calculado direto da métrica, com (ζ, Q, φ) genéricos, de forma exata em 3 jatos racionais. O controle negativo detecta 4/4 alterações. Inverti a matriz 5×5 à mão. R4 é local e vale para perturbações seculares |
| R5 recuperação de raio | **OK** (L1 de escopo) | `r5_free_tail_independente.py`, `r5_cota_e_extensao.py`, `r5_saida/` (+ logs); derivação §R5 | `radius_recovery_geral.py` reproduz o JSON byte a byte. Recalculei a cauda livre de forma independente (inversa triangular na base T_n, sem o formalismo de bandas) em 4 pesos e 3 caixas: free_tail sempre majora. Provei free_tail <= cota em papel e em 4000 amostras. Achei D(k) em RT (𝒦1♯). Conferi a reponderação de β (a parte linear é ∝ 40/9). Mostrei que a extensão a |s| <= S para S qualquer funciona |
| R6 unicidade módulo gauge | **OK com ressalva** (L3, L5) | derivação §R6 | Conferi os três itens. O Lema G se aplica a X1 − X2 + Y sem exigir que o campo seja Floquet, só δω |

Nenhuma peça recebeu LACUNA grave ou ERRO.

---

## R0. Enunciado e definições

**O fator e^{−4K}.** É o natural. RT (Main result) dá Θ*ḡ = e^{−2K}ḡ e ζ∘Θ = ζ + K. Logo
(Θ²)*ḡ = e^{−4K}ḡ, e um δg "relativo ao fundo" com expoente s satisfaz (Θ²)*δg = e^{4πs}e^{−4K}δg. Conferi
que isso equivale a δζ, δW/W e δφ serem e^{sτ}·(4π-periódicas):
- δζ = −δa/(2a);
- a(Θ²p) = e^{8πμ}e^{−4K}a(p), porque u∘Θ² = e^{−4πμ}u e a componente g_{+-} escala com e^{8πμ};
- δa escala com o mesmo fator vezes e^{4πs}.

Para φ não há fator, porque φ*∘Θ² = φ*. Coerente com a definição de Floquet de `GAUGE_GLOBAL.md`
(δω∘Θ² = e^{4πs}δω).

**Periodicidade de f_∓.** Com f_- = −μe^{−sτ}δg_{++}/g_{+-}:
- δg_{++}(Θ²p) = e^{8πμ}e^{4πs}e^{−4K}δg_{++}(p);
- g_{+-}(Θ²p) = e^{8πμ}e^{−4K}g_{+-}(p).

Logo f_- é 4π-periódica. HREC não define f_∓ (remete a `GLOBAL_NULL_GAUGE_TRANSPORT.md`): ver L8.

**Hipóteses escondidas ou excessivas na definição de "modo físico":**
- **"Equivale a dizer 'Floquet num gauge analítico'"** é uma afirmação sem prova (L2). Ela não é usada no
  lema, que assume Floquet exato.
  - Para valer, seria preciso resolver a equação cohomológica (Θ²)*Z − e^{4πs}Z = Y para o gauge Y que
    mede a falha de Floquet, com Z analítico.
  - Isso não é trivial nas ressonâncias.
- **Analiticidade em M̂ inteiro, até ξ = 41/40.** É mais do que a prova usa. R1 e R5 só precisam de
  analiticidade numa vizinhança complexa de 0 <= ξ <= 1, eixo e cone incluídos, e um pouco além de ξ = 1
  para o transporte com c = 1. R5 recupera o resto. Isso não é furo: a definição poderia ser
  enfraquecida. A analiticidade **através do cone** é usada de fato (transporte com centro c = 1, e R5) e
  é o conteúdo de H-reg.
- **Só autovetores de Floquet.** A definição exclui soluções τ-polinomiais × Floquet, isto é, cadeias
  físicas. T2 conta "com multiplicidade algébrica" (L9).
- O resto (simetria esférica, μ e coordenadas fixos, fundo de RT) está declarado em T2.

**Classe de gauge.** "C¹ em M̂ (Lema G)" coincide com `GAUGE_GLOBAL.md`. A versão G′ em HREC diz
"contínuo em M̂ com δω_X na classe analítica", mas o Lema G′ exige também **C¹ fora do cone**: os Passos
1–2 diferenciam X (L3). Os gauges construídos em R1 e R2 são analíticos e estão nas duas classes.

**O item 3 e a bijeção.**
- O item 3 segue de R6 e do Lema G. Ressalva de notação (L5): em s = μ − im0/2 com m0 ≠ 0, a direção de
  gauge no perfil é ℂ·e^{(μ−s)τ}g* = ℂ·e^{im0τ/2}g*, não ℂg*.
- A bijeção do §1 usa, além de R1–R6:
  - a volta (`GEOMETRIC_FLOQUET_RECONSTRUCTION.md`, não revisada aqui);
  - a passagem de ker L ∩ ker C para δΩ = 0 completo, que é o argumento ♯ de RT linearizado.
- A bijeção é afirmada "para Re s > 0", mas o lema só cobre |s| <= 3,1 (L1).

## R1. Transporte

**Derivação de T_{±1}** (sympy, em `r1_transporte_numerico.py`, e à mão).
- Com u_σ = −(1 + σξ)e^{−μτ}, a inversa da jacobiana dá ∂_{u_+} = (e^{μτ}/2μ)[∂_τ + μ(ξ − 1)∂_ξ].
- (L_X ḡ)_{++} = 2g_{+-}∂_+X^-, porque ḡ_{++} ≡ 0. Logo δg'_{++} = 0 ⇔ ∂_+X^- = −δg_{++}/(2g_{+-}).
- Com X^- = e^{(s−μ)τ}y_-, vem ∂_+X^- = (e^{sτ}/2μ)·T_1(s)y_-, conferido simbolicamente: a razão é
  exatamente 1/(2μ).
- Portanto T_1(s)y_- = −μe^{−sτ}δg_{++}/g_{+-} = f_-. O sinal e o fator μ conferem com a definição de f_-
  em `GLOBAL_NULL_GAUGE_TRANSPORT.md`. O caso σ = − é análogo.

**A fórmula (*) resolve a equação.** Seja Φ_t(τ, ξ) = (τ − t, c + e^{−μt}(ξ − c)). Vale
(∂_τ + μ(ξ − c)∂_ξ)(g∘Φ_t) = −(d/dt)(g∘Φ_t). Daí:
- T_c(s)∫_0^∞ e^{−(s−μ)t}g∘Φ_t dt = −∫ (d/dt)[e^{−(s−μ)t}g∘Φ_t] dt = g − lim_{t→∞}(…) = g.
- O limite é nulo porque |g∘Φ_t| <= D_c·L·e^{−μt}: g se anula em ξ = c, e L = ‖∂_ξ f‖.
- O integrando é <= D_c·L·e^{−Re s·t}, logo converge para todo Re s > 0, inclusive 0 < Re s < μ.
- Sem subtrair o traço, o integrando cresceria como e^{(μ − Re s)t}.
- O termo R_c f_c resolve o traço. Somando, T_c y = f_c + g = f.

**Analiticidade.** Para c ∈ Ω convexo, o caminho c + e^{−μt}(ξ − c) é combinação convexa de ξ e c, logo
fica em Ω. Com convergência uniforme, a solução é analítica em (τ, ξ) numa faixa × Ω.

**Unicidade fora da ressonância.**
- O coeficiente de Taylor y_n em ξ = c satisfaz (∂_τ + s + μ(n − 1))y_n = 0.
- Uma solução periódica não nula exige s + μ(n − 1) + im/2 = 0.
- Para n >= 1, a parte real é >= Re s > 0. Para n = 0, isso é a ressonância s ∈ μ − (i/2)ℤ.
- A identidade analítica estende a unicidade a Ω. A unicidade é na classe analítica; a unicidade entre
  gauges C¹ é o papel de R6 e do Lema G.

**Colagem e reflexão.** Para δg analítico na forma cartesiana α dτ² + 2βdτ(x·dx) + γ(x·dx)² + d|dx|², com
α, β, γ, d funções de (τ, ξ²), e V = a∂_τ + b∂_ξ:

    δg(V, V) = αa² + 2βξab + (γξ² + d)b².

Com ∂_± (b = e^{μτ}(ξ ∓ 1)/2), a expressão de δg_{++}(τ, ξ) é um polinômio em ξ com coeficientes pares.
Ele é idêntico a δg_{--}(τ, −ξ). Como g_{+-} = −e^{−2ζ}e^{2μτ}/(2μ²) é par:
- f_- se estende analiticamente a |ξ| < 41/40;
- f_-(τ, −ξ) = f_+(τ, ξ).

Essa é a colagem pedida, e é exatamente a regularidade cartesiana de δg. Além disso:
- PT_1P = T_{−1}, e pela unicidade y_+ = Py_-;
- A = (y_+ + y_-)/(2μ) é par e B = ((1 + ξ)y_- − (1 − ξ)y_+)/2 é ímpar, logo
  X = e^{sτ}(A∂_τ + (B/ξ)x·∂_x) é analítico no eixo.

**Testes.**
- `test_global_null_gauge`: 7/7 OK. São polinômios racionais: inversão modo a modo, ressonância e alias,
  obstrução de cadeia, reflexão e cone fixo.
- **Teste próprio** (`r1_transporte_numerico.py`):
  - **Fonte.** f_- = e^{0,7ξ}(1 + 0,3cos(τ/2)) + sin τ/(3 − ξ) + 0,2i·e^{iτ/2}/(2,5 + ξ), analítica e não
    polinomial.
  - **Resolução.** y vem de (*) por quadratura. O traço do cone entra por FFT e a diferença
    f(·, ξ) − f(·, c) é calculada sem cancelamento.
  - **Pontos.** ξ ∈ {−0,9; 0; 0,5; 1; 1,02}, este último além do cone.
  - **Expoentes:** s = 0,09 + 0,3i e s = 0,15 − 1,1i, ambos com 0 < Re s < μ, e s = 0,7332.
  - **Resíduo** |T_1 y − f|, por diferenças finitas de 4ª ordem:
    - <= 1,5·10⁻⁸ no caso de decaimento lento (Re s = 0,09), limitado pela quadratura;
    - <= 3·10⁻¹¹ nos outros.
  - **Periodicidade:** |y(τ + 4π) − y| <= 3·10⁻¹⁵.
  - **Reflexão e paridades:** exatas em ponto flutuante.
  - **Ressonância.** Em m0 = 0 (κ = 2,01) e m0 = 1 (κ = 0,302 + 0,057i), a solução no complemento resolve
    T_1(s0)y_reg = f − κe^{im0τ/2} com resíduo 4·10⁻¹², e T_1(s0)(τe^{im0τ/2}) = e^{im0τ/2}.
  - **Saída:** "R1 NUMERICO: CONFERE".

## R2. A ressonância s ≡ μ

Seja s0 = μ − im0/2 (Re s0 = μ > 0) e κ := (f_-(·, 1))_{m0}, com f(τ) = Σ f_m e^{imτ/2}.

1. **Necessidade e identidade secular.**
   - Avaliando T_1(s0)y = f_- em ξ = 1: (∂_τ − im0/2)y(·, 1) = f_-(·, 1), e o coeficiente m0 do lado
     esquerdo é 0. Logo κ = 0 é necessário para y periódico.
   - T_1(s0)(τe^{im0τ/2}) = e^{im0τ/2} + (im0/2)τe^{im0τ/2} − (im0/2)τe^{im0τ/2} = e^{im0τ/2}: o termo
     radial some, porque a função não depende de ξ.
   - Com y_- = y_reg + κτe^{im0τ/2}, o traço de f_- − κe^{im0τ/2} tem coeficiente m0 nulo, e y_reg existe
     por (*) no complemento, conforme o teste de R1.
   - Para y_+: f_+(·, −1) = f_-(·, 1), logo a mesma constante κ, e y_+ = Py_reg + κτe^{im0τ/2}.
2. **Parte secular de X.**
   - X^± = e^{(s0−μ)τ}y_± = e^{−im0τ/2}y_{reg,±} + κτ.
   - Logo X = X_reg + κτ(∂_{u_-} + ∂_{u_+}) = X_reg + κτX0, com a **mesma** constante nas duas componentes.
   - Isso é necessário: para κ(∂_{u_-} + c∂_{u_+}), B ∝ (1 − c) + (1 + c)ξ, que só é ímpar (regular no
     eixo) se c = 1.
   - Sympy (1): X0 = e^{μτ}(μ⁻¹∂_τ + ξ∂_ξ), regular no eixo, porque ξ∂_ξ = x·∂_x.
3. **δ' = A + τB.**
   - L_{τX0}ḡ = τL_{X0}ḡ + 2dτ⊙X0^♭, e L_{τX0}φ* = τX0(φ*).
   - Ponha B := κL_{X0}(ḡ, φ*) e A := δ + L_{X_reg}(ḡ, φ*) + 2κdτ⊙X0^♭.
   - **B está no ansatz:** sympy (2) dá (L_{X0}ḡ)_{±±} = 0. Isso era esperado, porque X0 tem componentes
     constantes em u e ḡ_{±±} ≡ 0.
   - **A está no ansatz:** A_{±±} = δ'_{±±} − τB_{±±} = 0.
   - Sympy (5) confere que a parte ++ do termo extra é 2g_{+-}∂_+τ ≠ 0. É ela que cancela a fonte
     ressonante; o termo foi absorvido em A e não escapa.
   - **A é Floquet:**
     - (Θ²)^*X0 = e^{4πμ}X0, e e^{4πμ} = e^{4πs0}, porque m0 é inteiro;
     - dτ é invariante;
     - δ e X_reg são Floquet.
   - **A é analítico no eixo:** X0 e τ são regulares.
4. **Forma de δω'.**
   - A extração (δζ, δQ, δφ) é pontual e linear sobre funções, logo ela aplicada a A + τB dá V_A + τV_B.
   - δω é de 1ª ordem, com D^σ(τ) = σ. Sympy (4) confere, para as 7 funções e com A, B genéricos:
     δω[A + τB] = δω[A] + τδω[B] + σ·(variáveis de B). Não há outros termos em τ.
   - Sympy (3) confere, fora das equações, para as 7 funções: δω[L_{X0}(ḡ, φ*)] = e^{μτ}(Z + 1)ω.
   - Logo δω' = e^{s0τ}(h1 + τh0), com h0 = κe^{(μ−s0)τ}g*_perfil = κe^{im0τ/2}g*_perfil ≠ 0. O perfil h1
     é 4π-periódico e analítico.
5. **Cadeia.**
   - δ' resolve as equações linearizadas: δ resolve, e L_X(ḡ, φ*) resolve para todo X liso, por
     invariância por difeomorfismos.
   - R4 é local e não usa Floquet: todos os δΩ se anulam.
   - Em Ω_i^+, a única derivada temporal é ξ(D^+ + μ)ω_i^-, com coeficiente ξ em ∂_τ; os termos quadráticos
     e os termos μ·(linear) não têm derivada.
   - Após S = Ξ^{−1,reg}∘seleção, temos Ξ^{−1,reg}Ξ = 1. Na componente 1, (1 + P)/2 preserva ξh1, porque h1
     é ímpar. Logo S δΩ^+[τψ] = τ S δΩ^+[ψ] + ψ, que é o L-afim.
   - Portanto 0 = e^{s0τ}[L(s0)h1 + h0 + τL(s0)h0]. Tomando τ ↦ τ + 4π, vem L(s0)h0 = 0 e L(s0)h1 = −h0.
6. **Reetiquetamento e setor (m0 ≠ 0).**
   - L(s)(e^{iατ}h) = e^{iατ}L(s + iα)h, porque ∂_τ ocorre só com coeficiente 1.
   - Com α = −m0/2: L(μ)h1' = −κg*, com h1' = e^{−im0τ/2}h1, que ainda é 4π-periódico. É uma cadeia em
     s = μ exato, o μ de RT.
   - T h = σh(· + 2π), com σ = diag(1, 1, 1, −1), comuta com L, porque ω* ∈ A e as equações são
     equivariantes por ω4 ↦ −ω4. Além disso, g* = (Z + 1)ω* ∈ A.
   - Projetando com P_A = (1 + T)/2: L_A(μ)P_A h1' = −κg*. O argumento vale para m0 par ou ímpar.
7. **Espaço e contradição.**
   - **Em que espaço vive a cadeia.** h1' é analítico periódico. Por R5, aplicado à cadeia em s = μ
     (|μ| <= 3,1), h1' ∈ D(Y+), com o lado direito κg* ∈ D(Y+) (RADIUS_RECOVERY).
   - **Em que realização S4 certificou.** No toro ℓ² com sinal, setor A
     (`signed_operator_L.py`: "Norma do toro de raios (65/64, 5/4): peso κ1^{2m}ω2(n)"; T2_HIPOTESES, L-real).
   - **Inclusão.** Y+ ⊂ toro, com ‖x‖_toro <= ‖x‖_+/243, e Q = J_μ⁻¹ é o mesmo operador nos dois espaços.
     Logo P_A h1' está no domínio do toro.
   - **Em qual disco.** S4 identifica o zero do NK bordejado com o par exato (g_n, μ*) e dá DF(x*)
     invertível. O raio é 6,06·10⁻⁵ e a unicidade vale até 3,26·10⁻⁴.
   - **Contradição.** Com F(y, λ) = (L(λ)Q0y, ⟨h0, Q0y⟩ − 1) e uma cadeia L(μ)h1 = −g_n, o vetor
     (Q0⁻¹h1 − ⟨h0, h1⟩y*, 1) está no núcleo de DF(x*). Isso contradiz a invertibilidade. Logo κ = 0.

**Furos procurados e não encontrados:**
- termos em τ que escapem: não há, por (4) e (5);
- R4 em perturbações seculares: é local, vale;
- a realização de h1: D(Y+) ⊂ toro;
- m0 ≠ 0: resolvido por reetiquetamento;
- regularidade no eixo do gauge secular: X0 é regular;
- a constante κ igual nos dois lados: vale.

**O que falta é só redação (L4):** os passos 6 e 7 não estão escritos em HREC (estão como "a conferir"),
e o espaço de S4 (toro) não é citado.

## R3. Do ansatz às variáveis de RT: prova completa

**Fundo.** Ponha a = −2e^{−2ζ}/(μ²(u_+ + u_-)²) = −e^{−2ζ}e^{2μτ}/(2μ²) e b = e^{−2ζ}ξ²/W², com W = μ + ξ²Q.
Com du_σ = e^{−μτ}(μ(1 + σξ)dτ − σdξ), a parte 2D do fundo é

    ḡ_2D = a(du_-du_+ + du_+du_-) = e^{−2ζ}[−(1 − ξ²)dτ² − (2ξ/μ)dτ⊙dξ + dξ²/μ²],

conferida no script R4 pela inversa do frame de RT. No eixo vale e^{−2ζ}(−dτ² + dξ²/μ²).

**Perturbação no ansatz.** Uma perturbação esfericamente simétrica no ansatz tem
δg = (δa/a)ḡ_2D + δb·g_{S²} e δφ. Daí:
- δζ = −δa/(2a);
- de log b = −2ζ + 2log ξ − 2log W vem δW = −W(δζ + δb/(2b));
- δQ = δW/ξ².

**Regularidade no eixo.** Seja δ' analítico em coordenadas cartesianas. Uma 2-forma simétrica O(3)-invariante
e analítica se escreve δg = α dτ² + 2βdτ(x·dx) + γ(x·dx)² + d|dx|², com α, β, γ, d analíticas em (τ, ξ²).
Isso é o resultado padrão de invariantes: cada termo homogêneo de Taylor é um polinômio O(3)-equivariante
gerado por δ_ij e x_ix_j sobre os invariantes |x|²; a forma ε_{ijk}x^k dx^i⊙dx^j é nula por simetria.
Usando dξ = x·dx/ξ e |dx|² = dξ² + ξ²g_{S²}:
- δb = dξ²;
- o coeficiente de dξ² da parte 2D é γξ² + d.

Comparando com (δa/a)ḡ_2D, o coeficiente de dξ² dá

    −2δζ e^{−2ζ}/μ² = γξ² + d   ⇒   δζ = −(μ²e^{2ζ}/2)(γξ² + d),

que é par e analítica. Com δb/(2b) = dW²e^{2ζ}/2 e W² − μ² = ξ²(2μQ + ξ²Q²):

    δζ + δb/(2b) = (e^{2ζ}/2)[d(W² − μ²) − μ²γξ²] = ξ²·(e^{2ζ}/2)[d(2μQ + ξ²Q²) − μ²γ].

Logo δQ = −W·(e^{2ζ}/2)[d(2μQ + ξ²Q²) − μ²γ], analítica e par. As componentes dτ² e dτdξ não dão condição
extra, porque δg_{±±} = 0 já força a parte 2D a ser conforme a ḡ_2D. ∎

**Codificação de RT.** Para F par, P(D^-F) = −D^+F. Logo δω_i^+ = −Pδω_i^- (i = 2, 3, 4) e
P(δω1) = −δω1 valem automaticamente para potenciais pares. A "colagem" no eixo da codificação de RT é,
portanto, consequência da analiticidade cartesiana. Os perfis vivem num cilindro analítico em
|ξ| < 41/40.

## R4. As equações

- **A fórmula de RT, conferida de forma independente** (`r4_curvatura_rt.py`).
  - Montei a métrica 4D a partir do frame de RT, com ζ, Q e φ funções simbólicas genéricas, fora das
    equações, e calculei o Ricci pelos símbolos de Christoffel.
  - Comparei as componentes de frame (θθ, rr, 00, 0r) de Ric − 2dφ⊗dφ e e^{−2ζ}□φ com o lado direito de
    (dfkhdjhsdshkfd), com os Ω definidos por (eezuuirzr) e os ω por (ekjehee).
  - As **cinco** diferenças e as **cinco identidades de definição** são **exatamente 0** em 3 jatos
    racionais aleatórios. As identidades são Ω1^σ − Ω2^σ − Ω2♯^σ (σ = ±) e Ω^-_i − Ω^+_i (i = 2, 3, 4).
  - Os lados esquerdos não são nulos, por exemplo −3,69 e 9,64, logo o teste não é vácuo.
  - **Controle negativo** (`r4_controle_negativo.py`): alterar um coeficiente em cada fórmula dá diferenças
    de 2,18; 4,36; 2,59; −1,59. As quatro alterações foram detectadas.
- **Inversão.** Substituindo as identidades, com Ω2^± = a, Ω3^± = c, Ω4^± = d, Ω2♯^σ = b_σ e
  Ω1^σ = a + b_σ, as cinco combinações geométricas são:
  - δ_ij: 4a;
  - x^ix^j: 2(b_- + b_+) + 4c;
  - 00: (b_- + b_+) − 2a − 2c;
  - 0i: b_+ − b_-;
  - □: 2d.

  Todas nulas implicam a = d = 0, b_+ = b_- e c = b = 0. A matriz coincide com a de
  `test_geometric_reconstruction.py`.
- **Linearização.** O lado direito é linear nos Ω, com coeficientes que dependem só de x, e o fundo tem
  Ω = 0 e Ric − 2dφdφ = 0. Logo as equações linearizadas (para ξ ≠ 0) equivalem às cinco combinações de δΩ
  nulas. As identidades valem fora das equações, logo valem linearizadas. Assim todos os δΩ se anulam para
  ξ ≠ 0, e no eixo por continuidade.
- **O que `test_geometric_reconstruction.py` testa** (5/5 OK), com polinômios racionais:
  1. a identidade de W, D^σW − W(ω2^σ − ω3^σ) = (Ω1^σ − Ω2^σ − Ω2♯^σ)/2, e sua linearização, fora das
     equações;
  2. as diferenças de resíduos −2μξ(T_s b − ∂_ξ a) = −PδΩ − δΩ;
  3. a reconstrução de (q, z, p) a partir de h;
  4. a reconstrução de cadeias;
  5. a invertibilidade da matriz 5×5 transcrita.

  Ele **não** testa a fórmula de curvatura (dfkhdjhsdshkfd) em si; o meu script cobre isso (L7).
- **Perturbações seculares.** R4 vale para elas. O argumento é pontual (identidades e curvatura em cada
  ponto) e só exige δ' de classe C², no ansatz. Floquet não entra.

## R5. Recuperação de raio

- **`radius_recovery_geral.py`.** Rodei a partir da pasta da revisão (`r5_saida/`) para não sobrescrever
  `build/t2`. O JSON é idêntico ao do repositório (`diff` vazio). Caixas:
  - 9/8 → 1024×16384, com fraco 0,603 e forte 0,197;
  - 1 + 2⁻⁸ → 2²⁰×2²³, com fraco 0,443;
  - 1 + 2⁻¹² → 2²⁸×2³¹, com fraco 0,438.
- **free_tail independente** (`r5_free_tail_independente.py`).
  - Os operadores J = im/2 + μ((1 + ξ)∂_ξ + 1) e K = im/2 + μ(ξ∂_ξ + 1) foram montados na base T_n com
    `numpy.polynomial.chebyshev`, e a inversa triangular com `solve_triangular`. A diagonal confere,
    μ(n + 1) + im/2, e a inversa tem resíduo < 10⁻⁹.
  - Norma: |c0| + Σ k^n|c_n|, a norma simétrica de RT.
  - Calculei o supremo das somas de coluna na cauda amostrada: |m| <= 40 e alguns até 400, n <= 700. Duas
    versões:
    - blocos reais |Re q| + |Im q|, que é a norma exata nos graus de liberdade reais;
    - (3/2)|q|, que é o que free_tail majora.
  - Para os pesos 9/8, 17/16, 33/32 e 5/4 e as caixas 8×24, 16×64 e 32×128, **free_tail majora em todos os
    casos**. Exemplos:
    - k2 = 9/8, 32×128: amostra 0,929 (real) e 1,105 (3/2|q|), contra free_tail 1,186;
    - k2 = 33/32, 32×128: 2,88 e 3,40, contra 4,53.

    É razoavelmente justo e não é vácuo. Em 5/4 a amostra também fica abaixo de structured_input_tail
    (0,628 em 32×128).
- **Prova em papel de free_tail para peso genérico r = 1/k2.** Uso a estrutura de coluna de RT:
  J_μ⁻¹ = Δ⁻¹[1, f]Λ⁻¹(1 + Δ_{1,−1}), com f_j = (−im/2 + μj)/(im/2 + μ(j + 1)) e |f_j| <= 1, e K com |f| = 1.
  - A coluna n tem diagonal 1/|D_n|, primeiro termo descendente 2μn/(|D_n||D_{n−k}|), e os seguintes
    multiplicados por produtos de |f| <= 1.
  - A razão de pesos é r^{jk}, logo a soma é <= 1/|D_n| + (2μn/(|D_n||D_{n−k}|))·r^k/(1 − r^k).
  - **Cauda de Fourier, k = 1:** |D_n||D_{n−1}| >= x² + b² >= 2bx, com x = μn, o que dá 1/(b(1 − r)).
  - **Cauda de Fourier, k = 2:** com x = μ(n − 1), dá 1/(b(1 − r²)) + 2μ_max r²/(b²(1 − r²)).
  - **Cauda radial:** |D_j| >= μ(j + 1), o que dá (1/(μ(n + 1)))·max(1 + 2r/(1 − r), 1 + (2n/(n − 1))r²/(1 − r²)),
    decrescente em n.
  - Multiplicando por 3/2, obtém-se exatamente `free_tail`. **Nada depende de r = 8/9**, e Q preserva m.
- **Cota explícita.** Para b >= 1 e N >= 2:
  - 1/(b(1 − r²)) <= 1/(b(1 − r)) e 2μ_max r²/(b²(1 − r²)) <= 2μ_max/(b(1 − r));
  - 2r/(1 − r) <= 4/(1 − r) e (2N/(N − 1))r²/(1 − r²) <= 4/(1 − r).

  Logo free_tail <= cota. Conferido em 4000 amostras exatas (`r5_cota_e_extensao.py`, 0 violações). A cota
  tende a 0 quando M, N → ∞, para r < 1 fixo.
- **D(k) no artigo.** Ξ^{−1,reg} = 2Δ_{1,1}Δ⁻¹_{2,1} ("Basic linear operators", com a propriedade eq:sid4).
  - A cota ‖Ξ^{−1,reg}Γ1♯‖ <= 𝒦1♯ := 2κ2⁻¹(1 − κ2⁻²)⁻¹ = 2κ2/(κ2² − 1) = D(κ2) está nas "Auxiliary
    estimates for (U♯)⁻¹★♯".
  - Também ‖SΓ1‖ <= 𝒦1 = 2D(κ2) (eq. fdhkjhjdh4z894z89rhiu). Na tabela de parâmetros, 𝒦1♯ = 40/9 e
    𝒦1 = 80/9.
  - A cota vem do lema geral de bandas, ‖Δ⁻¹[k, f]‖ <= C_f(1 − κ2^{−k})⁻¹, válido para todo κ2 > 1.
- **Reponderação de β.** β+ = 5191/500 é max_j Σ_i η_i(A_nl + A_lin)_ij/η_j mais o erro do fundo
  (`COMPONENT_EXTERIOR.md`).
  - A_lin tem entradas múltiplas de (40/9)μ_max = D(5/4)μ_max (`component_exterior.gamma1_majorant`), logo
    escala exatamente com D(k')/D(5/4).
  - A_nl vem de normas de campos do fundo, que não crescem com pesos menores (álgebra ℓ¹ ponderada), e o
    mesmo vale para o erro.
  - Como D(k')/D(5/4) >= 1, β(k') <= β+·D(k')/D(5/4).
- **Existência forte e unicidade fraca.** Com x = J_μh ∈ Y_{κ'}, Gx é finito e QG ⊂ G, porque J⁻¹ é
  triangular em n. Logo T H G Gx ∈ Y+.
  - Neumann forte dá y ∈ Y+ ⊂ Y_{κ'}.
  - A unicidade fraca dá Tx = y, logo x ∈ Y+.
  - Para cadeias, o lado direito −Th_{j−1} − THGx_j é forte por indução.
  - Conferido; é o argumento revisado de `RADIUS_RECOVERY.md`.
- **Perfil analítico em algum Y_{κ'}.** Analiticidade numa vizinhança complexa de [0, 4π] × [−1, 1] dá
  decaimento geométrico de Fourier–Chebyshev, porque a elipse de Bernstein E_ρ cabe na vizinhança para
  ρ > 1 pequeno.
  - J_μh = (∂_τ + μD_O)h é analítico num domínio menor (Cauchy).
  - J_μ é injetivo nos analíticos: as soluções homogêneas (1 + ξ)^{−1−im/(2μ)} e ξ^{−1−im/(2μ)} são
    singulares em [−1, 1].
  - Logo h = Q(J_μh) ∈ D(Y_{κ'}).
- **Extensão a |s| arbitrário.** (β + S)·cota → 0 e (β+ + S)·structured_input_tail → 0. Há caixa para
  S = 10, 50 e 200 (`r5_cota_e_extensao.log`), o que conserta L1.

## R6. Unicidade módulo gauge

1. **δω' = 0.** As linearizações de (dhjhfkfhfe1–4) dão:
   - δQ = (δω1 − δω2^- − δω2^+)/(2ξ) = 0;
   - D^±δζ = 0, logo δζ = c1, porque D^± gera o espaço tangente e M é conexo;
   - D^±δφ = 0, logo δφ = c2.

   Floquet com Re s > 0 dá c(1 − e^{4πs}) = 0, logo c = 0. Então δ' = (−2δζḡ − …, δφ) = 0, e
   δ = −L_X(ḡ, φ*) é de gauge, com X analítico, logo C¹.
2. **δω' = c·e^{μτ}g*.** δ + L_{X − cX0} tem δω = 0, por linearidade e por δω[L_{X0}] = e^{μτ}g*
   (sympy (3)), e cai no caso 1.
3. **Dois gauges.** Sejam X1 e X2 dois gauges que levam representantes da mesma classe [δ] ao ansatz. Os
   representantes diferem por L_Y, com Y na classe de gauge.
   - Z = X1 − X2 + Y é admissível: só se usa δg_{±±} = 0, como no Passo 1 do Lema G.
   - Z é C¹ (ou de classe G′).
   - δω_Z = e^{sτ}(h1 − h2) é Floquet.
   - O Lema G não exige Z Floquet, só δω_Z. Logo δω_Z ∈ ℂ·e^{μτ}g*, e é nulo se s ≢ μ.
   - "h ∈ ℂ(e^{(μ−s)τ})g* ⇔ δ de gauge": (⇐) pelo Lema G aplicado a X + Y; (⇒) pelo item 2.

   Para G′ vale a ressalva L3.

---

## Lacunas

| id | gravidade | descrição | conserto proposto |
|---|---|---|---|
| L1 | **moderada** (escopo) | O lema exige \|s\| <= 3,1, mas o §1 afirma a bijeção "para Re s > 0", e T2 conclui "1 modo físico" em todo Re s > 0. Com o reetiquetamento de Floquet (Im s em (−1/4, 3/4]), \|s\| <= 3,1 cobre só Re s <~ 3,0. F3 exclui raízes de L com Re s >= 2,93, mas, para aplicá-lo a um modo **físico**, é preciso pôr o perfil no domínio do toro, isto é, R5 naquele s | Enunciar R5 para \|s\| <= S qualquer: a prova não muda, e a caixa existe para todo S, como mostra `r5_cota_e_extensao.log` com S = 10, 50 e 200. Depois usar Y+ ⊂ toro e F3. Ou restringir a bijeção e T2 a \|s\| <= 3,1 |
| L2 | menor | "Equivale a dizer 'Floquet num gauge analítico'" não está provado: exige resolver uma equação cohomológica para o gauge, inclusive na ressonância | Retirar o "equivale", ou provar |
| L3 | menor | A classe G′ em HREC ("contínuo em M̂ com δω_X analítico") omite "C¹ fora do cone", que o Lema G′ usa (Passos 1–2) | Copiar a classe exata de `GAUGE_GLOBAL.md` §2.1 em HREC e T2 |
| L4 | menor (redação de R2) | Faltam por escrito o reetiquetamento h ↦ e^{−im0τ/2}h, a projeção P_A (g* ∈ A, [L, T] = 0), a realização de S4 (toro ℓ² com sinal, setor A) com a inclusão D(Y+) ⊂ toro, e o argumento do sistema bordejado (DF invertível ⇒ sem cadeia) | Incorporar os passos 6–7 de §R2 deste relatório |
| L5 | menor (notação) | "h módulo ℂg*" e "(ℂg* se s ≡ μ)": para s = μ − im0/2 com m0 ≠ 0, a direção de gauge no perfil é ℂe^{im0τ/2}g* | Escrever ℂ·e^{(μ−s)τ}g* |
| L6 | menor | R3 estava em esboço | Usar a prova completa de §R3 |
| L7 | menor | R4 depende da fórmula (dfkhdjhsdshkfd), cuja derivação RT omite ("not included"). O teste do repositório só inverte a matriz transcrita | Incorporar `r4_curvatura_rt.py` e o controle negativo ao repositório: eles conferem a fórmula e as identidades fora das equações |
| L8 | menor | HREC não define f_∓ e não deriva sua periodicidade a partir do fator e^{−4K} | Copiar a definição f_- = −μe^{−sτ}δg_{++}/g_{+-} e a conta de §R0 |
| L9 | menor a moderada (escopo de T2) | O Lema R trata só de autovetores de Floquet. T2 diz "contados com multiplicidade algébrica", o que, do lado físico, pede soluções τ-polinomiais × Floquet (cadeias físicas) | Estender a ida a jatos: transporte de jatos de `GLOBAL_NULL_GAUGE_TRANSPORT.md` e termos seculares como em R2. Como as quatro raízes são simples, uma cadeia física daria uma cadeia de L, o que é contradição. Ou enunciar T2 sem "algébrica" do lado físico |

Nenhuma lacuna invalida o lema na sua forma enunciada (|s| <= 3,1, autovetores de Floquet).

## O que não verifiquei

- **Os certificados numéricos citados:**
  - o NK de S4 (θ, r, identificação);
  - L-real (‖J⁻¹T‖_toro, a inclusão 1/243);
  - L-afim;
  - F3;
  - o valor β+ = 5191/500: só conferi sua composição;
  - a derivação de structured_input_tail: só a comparei com a amostra em k2 = 5/4.
- **A volta** (`GEOMETRIC_FLOQUET_RECONSTRUCTION.md`, parte direta) e a passagem ker L ∩ ker C ⇒ δΩ = 0
  completo (argumento ♯ de RT, linearizado). A bijeção do §1 depende delas.
- **O teorema de representação** de 2-formas O(3)-invariantes analíticas (α, β, γ, d analíticas em ξ²):
  citado como padrão.
- **A conferência de R4** é exata, mas em 3 jatos aleatórios, não uma simplificação simbólica completa.
  É uma identidade racional nos jatos; a evidência é do tipo Schwartz–Zippel.
- **O teste de R1** é numérico (quadratura), não rigoroso.
- **Efeito colateral:** rodar os testes do repositório pode ter gerado `.pyc` em `scripts/__pycache__/`, que
  está no `.gitignore`. Nenhum arquivo rastreado foi alterado (`git status` limpo em docs, scripts e
  build/t2).

## Artefatos (nesta pasta)

- `r1_transporte_numerico.py` / `.log`: R1 com fonte não polinomial, ressonância e derivação sympy de T_1.
- `r2_gauge_x0_e_cadeia.py` / `.log`: R2 (X0, B no ansatz, δω[L_{X0}] = e^{μτ}(Z + 1)ω, estrutura secular,
  termo ++).
- `r4_curvatura_rt.py` / `.log`: fórmulas de curvatura de RT e identidades, com Ricci 4D direto.
- `r4_controle_negativo.py` / `.log`: controle negativo.
- `r5_free_tail_independente.py` / `.log`: cauda livre recalculada.
- `r5_cota_e_extensao.py` / `.log`: free_tail <= cota (4000 amostras) e extensão a S arbitrário.
- `r5_radius_recovery_geral.log`, `r5_saida/build/t2/radius_recovery_geral.json`: rodada do script do autor,
  idêntica ao repositório.
- `manifest_inicio.txt`, `manifest_fim.txt`.

sympy foi instalado no scratchpad da sessão (`pip install --target`). O artigo de RT foi extraído no
scratchpad, não nesta pasta.

## Manifesto (fim)

`LC_ALL=C sha256sum -c .snapshots/2026-10-03-revisao-R/MANIFEST.sha256` às 11:55 (−04) deu 241/241 OK, com
código de saída 0 e nenhuma falha (`manifest_fim.txt`). HEAD continua 8687cc3. `git status` não mostra
alterações em `docs/`, `scripts/`, `build/t2/` nem `.cache/`.
