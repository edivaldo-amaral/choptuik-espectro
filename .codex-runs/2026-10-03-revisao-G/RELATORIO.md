# Revisão independente do Lema G (`docs/GAUGE_GLOBAL.md`), 03/10/2026

Revisor: agente independente (Claude), sem acesso aos documentos proibidos pela TAREFA. Objetivo:
tentar quebrar o Lema G. **Resultado: não quebrou.** A prova está correta. Há uma ressalva de
compatibilidade de hipóteses com T2 (R2, média) e ressalvas de redação (R1, R3–R7), nenhuma das quais
invalida o lema ou o uso que T2 faz dele.

## 0. Manifesto e ambiente

**Início (11:07).**
- `COMMIT.txt` = `b31da1cb6adea86e543d0c03a60c9a468cb65361`, igual a `git rev-parse HEAD`.
- `sha256sum -c MANIFEST.sha256`: **238/238 SUCESSO, 0 falhas** (`manifesto_inicio.txt`).
- O artigo de RT foi extraído do tarball do cache para `rt/`. `rt/anc/sourcecode/RefA.dat` e
  `rt/anc/RefAplusB.dat` são idênticos byte a byte aos do cache (`cmp`).

**Fim:** ver §7.

**Ambiente.**
- O `.venv` **não tem sympy nem mpmath**, ao contrário do que diz a TAREFA. Só tem numpy, scipy e
  python-flint. Instalei sympy 1.14 e mpmath 1.3 com `pip install --target` no scratchpad da sessão,
  fora do projeto, e rodei com `PYTHONPATH` apontando para lá e `OPENBLAS_NUM_THREADS=2`.
- Nenhum arquivo fora desta pasta foi criado ou alterado por mim. Os `.pyc` de `scripts/__pycache__`
  mantêm as datas de 17/09 e 01/10.
- Durante a revisão apareceram arquivos novos criados por outra pessoa: `scripts/radius_recovery_geral.py`,
  `docs/HREC_IDA.md` e `build/t2/radius_recovery_geral.json`. Não estão no manifesto e não os li.
- **Não rodei `scripts/gauge_global_centro.py`**, porque ele grava em `build/t2/`. Comparei o JSON
  existente com o meu cálculo.

## 1. Tabela de veredictos

| alegação | veredicto | artefato | o que foi feito |
|---|---|---|---|
| G1 definições | **OK com ressalva** (R1, R2) | §2.1 (derivação); `g_simbolico.py` [A]; `f1_proprio.json` (F1′) | Conferi "gauge", "admissível" e "Floquet" contra RT e T2. Listei as hipóteses escondidas: μ fixo não está dito no lema (R1), e a classe de gauge de H-rec tem de coincidir com a do lema (R2). Proponho um reforço que dispensa o C¹ no cone. |
| G2 Passo 0 | **OK** | §2.2 (derivação) | Média de Haar, classe C¹, ausência de campo invariante em S², tangência ao eixo. |
| G3 Passos 1–2 | **OK com ressalva** (R5) | `g_simbolico.py` [A], [B], [F], [H] | A métrica de RT foi derivada do frame e (L_X g)_{±±} conferidas com X genérico. "Forma do ansatz" exige mais que δg_{±±} = 0, mas essa condição é necessária, e é só ela que se usa. Conferi a conexidade das curvas de nível e o C¹ de f_±. |
| G4 Passo 3 | **OK** | `g_simbolico.py` [C], [H]; §2.4 | Imagem de u_+ = (−∞, 0), cobertura do eixo, δτ e δξ simbólicos, tangência por continuidade. |
| G5 Passo 4 | **OK com ressalva** (R3) | `g_simbolico.py` [D]; `g5_equivariancia_linear.py`; `g5_perfis_nulos.py` | O argumento por fluxos é informal. A versão linear vale, e eu a escrevi e conferi. f̃(u) = Λ⁻¹f(Λu) confere, e f = u^n dá μ(1 − n) para n = 0..6 e três valores de μ. |
| G6 Passo 5 | **OK** | `g_simbolico.py` [E]; `f1_proprio.json`; §2.6 | Fórmula de X(φ*), sinais de D^±(u_∓), eixo, antiperiodicidade, zeros isolados, aberto exterior não vazio. Nenhum contraexemplo: a única saída seria ω4*(·, 0) ≡ 0. |
| G7 F1 | **OK** (R4 cosmética) | `f1_proprio.py`, `f1_proprio.json`, `f1_refB_saida.txt` | c_1 recalculado com código próprio em racionais exatos, idêntico ao do autor. Base, simetria, índice 3 = ω4, bola e pesos conferidos no artigo. ‖RefB‖ <= 2⁻²⁵ recalculado exatamente. Via RefAplusB, a cota não depende de 𝓛2. |
| G8 Passos 6–7 | **OK** | §2.8 (derivação) | Os três regimes de \|ρ\| são completos para Re s > 0. ρ = 1 ⇔ s ≡ μ mod (i/2)ℤ. Cobre os setores A e B. O caso Re s = μ com ρ ≠ 1 também. Exemplos de nitidez. |
| G9 Passo 8 e T2 | **OK com ressalva** (R2, R6, R7) | `g_simbolico.py` [F], [G]; `test_gauge_action_saida.txt` | δω_{X0} = e^{μτ}(Z + 1)ω confere nas 7 funções, com ζ, Q e φ genéricos, e X0 é admissível. `test_gauge_action.py`: 3/3 OK. O uso em T2 (11b e "Conta") é exatamente o que o lema prova. Obs. 1 e 2 corretas; a última frase da Obs. 1 é heurística. |

## 2. Detalhamento por alegação

### 2.1 G1 (definições)

- **Gauge.** δg = L_X g*, δφ = X(φ*), com X C¹ na variedade 4D ℝ × B_{41/40}(0) ⊃ M × S² ∪ eixo,
  possivelmente complexo. É a noção correta de difeomorfismo linearizado.
  - Observação: X não precisa ser tangente a ∂M; X0 não é. No nível linear isso é irrelevante, porque
    a perturbação é um tensor no domínio fixo. Mas a noção é só infinitesimal: o fluxo de X não
    preserva M (ver R3).
- **Admissível.** Preserva a forma linearizada do ansatz de RT, com **μ = μ\* fixo e coordenadas
  (τ, x) fixas**.
  - A hipótese μ fixo é **essencial**. Se δμ ≠ 0, as direções nulas ker(σdξ − μ(1 + σξ)dτ) mudam, e
    nos u_± antigos aparece δg_{++} ≠ 0. O Passo 1 cairia.
  - No lema, μ fixo está implícito ("D^σ fixos, dependem só de μ"). T2 o declara ("Identificação do
    fundo: ... μ e as coordenadas (τ, ξ) fixos"), e L(s) = L(0) + sI age só em h. → R1.
- **Modo de Floquet.** A definição δω∘Θ² = e^{4πs}δω equivale a δω = e^{sτ}h com h de período 4π.
  - Cobre os dois setores, porque B também é 4π-periódico.
  - O expoente é definido mod (i/2)ℤ.
  - Conferi a contagem de g* nas faixas de T2. A cópia de g* com perfil de paridade B,
    e^{−iτ/2}(Z + 1)ω*, tem expoente μ + i/2 em B. Nas coordenadas de L_A ela fica em μ + i, fora de
    Q̃ = (−1/4, 3/4]. Logo cada classe de largura 1/2 contém **um** representante de g*, coerente com
    S4 ("na faixa B não há gauge").
- **Hipóteses escondidas, item a item.**
  - *μ fixo:* R1.
  - *O domínio é o M de RT e contém o cone:* sim. Em (τ, ξ), M = ℝ × (0, 41/40) (`g_simbolico.py`
    [H]), e o cone é ξ = 1.
  - *É legítimo exigir X C¹ através do cone?* Sim, porque é a definição de difeomorfismo da variedade
    em que o cone é interior. Além disso, isso **não é essencial para T2** (reforço abaixo).
  - *X complexo:* tratado. Como δω_X = δω_{Re X} + iδω_{Im X}, com as duas partes reais, a injetividade
    se aplica separadamente.
  - *δω no espaço de raízes sem X Floquet:* tratado pelo Passo 6. A injetividade força
    (Θ²)*X = e^{4πs}X, sem supor X Floquet.
    - Vetores generalizados (e^{sτ}(h1 + τh0)) **não** são cobertos pelo enunciado.
    - T2 não precisa deles, porque as quatro raízes são simples (S4, S3b, NK).
  - *Compatibilidade com H-rec:* T2 chama de gauge um elemento de E com e^{sτ}h = δω_X. Para que os
    modos de puro gauge físicos caiam nessa classe, o gauge de fixação de H-rec (ida) tem de ser C¹
    como no lema. Isso não está dito. → R2.
- **Reforço proposto (não é lacuna): o C¹ no cone pode ser trocado por "δω_X na classe analítica de
  T2".** Suponha X **contínuo** em M, C¹ fora do cone, e δω_X analítico na faixa de RT.
  - A continuidade no cone é necessária: o Passo 2 usa que {u_+ = w} ∩ M atravessa o cone, e sem
    continuidade X^+ poderia saltar ali.
  - Com ela, os Passos 0–6 valem sem diferenciabilidade no cone. A injetividade não a usa, porque o
    exterior usa só a continuidade de f.
  - Para Re s >= μ, o Passo 7 já só usa continuidade (Obs. 4). Para 0 < Re s < μ:
  1. δφ = X(φ*) é analítico, por ser o potencial de δω4^σ.
  2. No eixo, −(e^{μτ}/μ)f(−e^{−μτ})ω4(τ, 0) é analítico. Como f é contínua e os zeros de ω4(·, 0)
     são isolados, f é analítica em (−∞, 0).
  3. Perto de um ponto do cone (τ0, 1) com ω4^+(τ0, 1) ≠ 0,
     f(u_-) = [2μe^{−μτ}X(φ*) + f(u_+)ω4^-]/ω4^+ é analítica. Logo f é analítica em u = 0.
  4. Com f(Λu) = ρf(u), os coeficientes de Taylor satisfazem a_n(Λ^n − ρ) = 0. Então f = a_n u^n e
     s ≡ μ(1 − n). Re s > 0 força n = 0.

  O fato necessário é **F1′: ω4^+(·, 1) = −ω4(·, −1) ≢ 0**. Calculei o coeficiente m = 1 em ξ = −1
  (θ = π): c_1 = Σ_n v_{1n}(−1)^n = −0,0794729 + 0,0977527i para RefA (`f1_proprio.json`). Pela
  mesma bola, |c_1(ω*)| >= 0,0977527 − 3·10⁻⁸ > 0.

### 2.2 G2 (Passo 0)

- Seja δ = L_X(g*, φ*) admissível. A forma do ansatz é SO(3)-invariante, logo R*δ = δ para todo
  R ∈ SO(3).
- Como R*g* = g* e R*φ* = φ*, vale R*L_X(g*, φ*) = L_{R*X}(g*, φ*). Integrando em Haar,
  L_{X̄}(g*, φ*) = δ, com X̄ = ∫R*X dR. A derivada passa para dentro da integral: integrando C¹ em
  grupo compacto agindo suavemente, convergência dominada nos compactos. Logo X̄ é C¹ e invariante.
- No ponto (τ, ξ, ω) ∈ M × S², o estabilizador SO(2) age na parte tangente a S² por rotação, que só
  fixa o 0. Logo X̄ não tem componente angular, e X^τ, X^ξ só dependem de (τ, ξ).
- No eixo x = 0, SO(3) age em T_p = ℝ_τ ⊕ ℝ³ pela representação padrão em ℝ³. O único vetor
  invariante é a direção τ. Em cartesianas, X^i = X^ξ x^i/ξ com |X^x| = |X^ξ|, e a continuidade em
  x = 0 dá X^ξ(τ, ξ) → 0 quando ξ → 0. Essa é a tangência usada no Passo 3.
- Complexo: igual. Comutação com Θ: Θ comuta com SO(3), logo a média comuta com (Θ²)*.

### 2.3 G3 (Passos 1–2)

- **Métrica.** A partir do frame (ekhfkjfhkfhfkhf) de RT, calculei g^{-1} em (τ, ξ, θ, φ),
  inverti e passei a (u_-, u_+). Obtive g_{--} = g_{++} = 0, g_{-+} = −2e^{−2ζ}/(μ²(u_+ + u_-)²) e
  g_{θθ} = e^{−2ζ}ξ²/W², com W = μ + ξ²Q (`g_simbolico.py` [A], [A1]).
  - Isso confere a fórmula do "Main result".
  - Conferi também, em cartesianas, Σ e_i⊗e_i = μ² r̂r̂ + W²(I − r̂r̂).
  - D^σ age em funções radiais como σ∂_τ + μ(1 + σξ)∂_ξ.
- **Derivada de Lie** ([B]). Para uma métrica dupla-nula genérica G(du_-du_+ + du_+du_-) + R²g_{S²} e
  X = X^-(u)∂_- + X^+(u)∂_+:
  - (L_X g)_{++} = 2g_{-+}∂_+X^- e (L_X g)_{--} = 2g_{+-}∂_-X^+;
  - as componentes mistas são nulas, e a parte angular é proporcional a g_{S²};
  - com componentes angulares genéricas, (L_X g)_{++} não muda.
- **"Forma do ansatz" exige mais que δg_{±±} = 0** ([F]).
  - O ansatz linearizado é (δζ, δQ, δφ). δζ vem de δg_{-+} = −2δζ·g_{-+}, e δW de
    δ(r²) = r²(−2δζ − 2δW/W). Isso sempre se resolve para ξ > 0.
  - O que sobra é a **regularidade de δQ = δW/ξ² no eixo** (e paridade/analiticidade na realização).
    A positividade W > 0 é aberta e irrelevante no nível linear.
  - Para X_f com f_+ = f_- = f, obtive:
    - δζ = X(ζ) + (f(u_-) + f(u_+))/(u_+ + u_-) − (f′(u_-) + f′(u_+))/2;
    - δW = X(W) + W·E_f, com E_f = [(u_- − u_+)(f′(u_-) + f′(u_+)) − 2(f(u_-) − f(u_+))]/(2(u_- − u_+)).
  - Expandindo em d = u_- − u_+ = 2ξe^{−μτ}, E_f = f‴(c)d²/12 + o(d²). Logo δQ é regular se f ∈ C³
    perto do eixo, mas **com f só C¹ pode não ser**. Então a classe admissível é possivelmente menor
    que "C¹ com δg_{±±} = 0".
  - Isso **não afeta o lema**, porque o lema usa só a condição necessária δg_{±±} = 0, que conferi ser
    necessária: g_{±±} ≡ 0 para todo (ζ, Q) com μ fixo. Exemplos: f = 1, u, u² dão E_f = 0, e
    f = u³, u⁴ dão E_f/ξ² = 2e^{−2μτ} e −8e^{−3μτ}, regulares. → R5 (redação).
- **Curvas de nível** ([H] e conta direta). Em (τ, ξ), u_+ < 0 ⇔ ξ > −1, u_+ < u_- ⇔ ξ > 0 e
  u_- < −u_+/81 ⇔ ξ < 41/40.
  - {u_- = v} ∩ M = {u_+ < 0, u_+ < v, u_+ < −81v} = (−∞, min(0, v, −81v)): o mínimo é v para
    v < 0, 0 para v = 0 e −81v para v > 0. É um intervalo.
  - {u_+ = w} ∩ M = (w, −w/81) para w < 0, e vazio para w >= 0.
  - As duas famílias são conexas, logo ∂_+X^- = 0 ⇒ X^- = f_-(u_-) em todo M.
  - u_- percorre ℝ: v > 0 é atingido com u_+ < −81v.
- **C¹.** Na reta τ = τ0, u_- = −(1 − ξ)e^{−μτ0} é um difeomorfismo afim de ξ ∈ (0, 41/40) sobre
  (−e^{−μτ0}, e^{−μτ0}/40), e X^{u_-} = du_-(X) é C¹ em ξ. A união sobre τ0 cobre ℝ, logo
  f_- ∈ C¹(ℝ). f_+ ∈ C¹((−∞, 0)) da mesma forma.

### 2.4 G4 (Passo 3)

- u_+ = −(1 + ξ)e^{−μτ}, com ξ ∈ (0, 41/40) e τ ∈ ℝ, cobre (−∞, 0). O eixo é u_± = −e^{−μτ}, que
  cobre (−∞, 0).
- Em M, X(u_+ − u_-) = f_+(u_+) − f_-(u_-). Em (τ, ξ),
  X(u_+ − u_-) = −2e^{−μτ}(X^ξ − μξX^τ), que tende a 0 no eixo (§2.2). Tomando o limite
  (u_-, u_+) → (u, u) dentro de M, com f_- e f_+ contínuas, vale f_+(u) = f_-(u) para todo u < 0.
- δτ = e^{μτ}[f(u_+) + f(u_-)]/(2μ) e δξ = e^{μτ}[(1 + ξ)f(u_-) − (1 − ξ)f(u_+)]/2 conferem
  simbolicamente com f_± genéricas, e em ξ = 0 vale δξ = e^{μτ}[f_- − f_+]/2 ([C]).

### 2.5 G5 (Passo 4)

- **O argumento do documento, por fluxos, é informal.**
  - Ω(ψ_ε*(g*, φ*)) não está definido para ε ≠ 0, porque ψ_ε*(g*) só está na forma do ansatz em
    primeira ordem.
  - Além disso, o fluxo de um X C¹ não preserva M.
  - A conclusão vale pela **versão linear**, que dispensa fluxos:
    1. Θ preserva o ansatz: em (τ, x), Θ é τ ↦ τ + 2π, e o frame tem coeficientes independentes de
       τ. Logo Ω(Θ*(g, φ)) = Ω(g, φ)∘Θ para (g, φ) no ansatz. Derivando ao longo de uma curva no
       ansatz: DΩ(g*, φ*)[δ]∘Θ² = DΩ((Θ²)*(g*, φ*))[(Θ²)*δ] = DΩ(e^{−4K}g*, φ*)[(Θ²)*δ].
    2. Invariância conforme: Ω(e^{2c}g, φ) = Ω(g, φ) dá
       DΩ(e^{2c}g, φ)[e^{2c}δg, δφ] = DΩ(g, φ)[δg, δφ]. Conferi a invariância das 7 funções por
       ζ ↦ ζ + c ([G]).
    3. (Θ²)*(L_X g*) = L_{(Θ²)*X}((Θ²)*g*) = e^{−4K}L_{(Θ²)*X}g*, e
       (Θ²)*(X(φ*)) = ((Θ²)*X)(φ*), porque φ*∘Θ² = φ*.
    4. Juntando: δω_X∘Θ² = δω_{(Θ²)*X}. (Θ²)*X é C¹ e admissível, porque Θ² preserva o eixo e o
       ansatz.
  - → R3.
- **Conferência direta** (`g5_equivariancia_linear.py`). Com fundo genérico satisfazendo
  ζ(τ + 4π) = ζ + 2K e W, φ 4π-periódicos, as variações de X_f em τ + 4π são as de X_{f̃} em τ, para
  δζ, δW e δφ: **3/3 OK**.
- **Direção do pullback.** Em (τ, ξ), dΘ² = id, logo ((Θ²)*X)(τ, ξ) = X(τ + 4π, ξ). Conferi
  simbolicamente que isso é X_{f̃}, com f̃(u) = Λ⁻¹f(Λu) e Λ = e^{−4πμ}, nas duas componentes
  ([D]). Em u, F(u) = Λu, dF = Λ e (F*X)^u(u) = Λ⁻¹X^u(Λu).
- **Perfis f = u^n** (`g5_perfis_nulos.py`, que importa `gauge_null_profiles.null_generator`). Para
  μ ∈ {1/6, 17/100, μ_RefA} e n = 0..6:
  - o pushforward calculado direto da troca de coordenadas coincide com o `null_generator`;
  - s_n = μ(1 − n);
  - o multiplicador de (Θ²)* é e^{4πs_n}.

  **21/21 OK.** No meu [D], A_n e B_n de `GAUGE_DOMAIN_CLASSIFICATION.md` conferem para n = 0..5.

### 2.6 G6 (Passo 5)

- D^+(u_+) = 0, D^-(u_-) = 0, D^+(u_-) = 2μe^{−μτ} e D^-(u_+) = −2μe^{−μτ}. Com φ genérica,
  X(φ) = (e^{μτ}/2μ)[f(u_-)D^+φ − f(u_+)D^-φ] ([E]).
- **D^σ fixo.** No ansatz com μ fixo, ω4^σ = D^σφ com D^σ = σ∂_τ + μ(1 + σξ)∂_ξ independente de
  (ζ, Q) ([A1]). Logo δω4^σ = D^σ(δφ) = D^σ(X(φ*)). Como D^+ e D^- são independentes em cada ponto
  (determinante 2μ ≠ 0) e M = ℝ × (0, 41/40) é conexo, δω4 = 0 ⇒ X(φ*) = c2.
- **Eixo.**
  - Pela codificação de RT (§"From 2D to 4D"), ω4^- = ω4 e ω4^+ = −Pω4. Em ξ = 0, ω4^+ = −ω4.
  - Conferência independente: Pφ = φ (RT) dá ∂_ξφ(τ, 0) = 0, logo ω4^∓(τ, 0) = ∓∂_τφ(τ, 0), coerente.
  - X(φ*) é contínuo até o eixo, logo −(e^{μτ}/μ)f(−e^{−μτ})ω4(τ, 0) = c2.
- **Antiperiodicidade.**
  - ω4 é real e antiperiódica (2Dpre (a) de RT). Pelo valor intermediário, tem um zero τ0, logo
    c2 = 0.
  - Ilustração, não prova: RefA tem exatamente 2 trocas de sinal de ω4(·, 0) em [0, 4π) (amostragem
    de 2048 pontos), antiperiodicidade até 5·10⁻¹⁶ e parte imaginária nula (`f1_proprio.json`).
- **Analiticidade.**
  - ω4 é analítica em |Im τ| < 2 log κ1 × elipse (RT, "Properties of TF"), logo os zeros de
    ω4(·, 0) ≢ 0 são isolados.
  - Então f = 0 num denso de (−∞, 0) e, por continuidade, em todo (−∞, 0).
- **Exterior.**
  - Para v > 0, os pontos (v, u_+) com u_+ < −81v estão em M e têm ξ > 1 (u_- > 0 ⇔ ξ > 1). Logo
    {ξ > 1, u_- ∈ I} ∩ M é aberto e não vazio para qualquer intervalo I ⊂ (0, ∞).
  - Nesse aberto, f(u_+) = 0 e c2 = 0, logo f(u_-)ω4^+ = 0.
  - ω4^+ = −Pω4 é analítica em ℝ × (−41/40, 41/40), que é conexo, e anulá-la num aberto a anula em
    tudo. Em ξ = 0, isso dá ω4 ≡ 0, contradizendo F1.
- **Tentativa de contraexemplo.**
  - δω_X = 0 ⇒ δω4 = 0 ⇔ X(φ*) constante. Pelo acima, isso força f = 0 se ω4*(·, 0) ≢ 0.
  - Famílias candidatas:
    - f = u (fase, X = −μ⁻¹∂_τ): δφ = −μ⁻¹∂_τφ*, não constante no eixo, porque
      ∂_τφ*(τ, 0) = −ω4(τ, 0) ≢ 0;
    - f = 1: δω = g* ≠ 0, cujo coeficiente m = 1 no eixo é (1 + i/(2μ))c_1 ≠ 0;
    - f suportada em u > 0: o exterior obriga ω4^+ ≡ 0;
    - f suportada em u < 0: o eixo obriga f = 0.
  - Um contraexemplo exigiria um fundo com φ constante ao longo do eixo (ω4(·, 0) ≡ 0), como um
    fundo CSS com f = u. F1 exclui isso. **Nenhum contraexemplo.**

### 2.7 G7 (F1)

- **Código próprio** (`f1_proprio.py`, sem importar nada do projeto). Lê o formato documentado por RT
  (Data structures; conferi `Field_to_file` em `rt/anc/sourcecode/Field.c`: as linhas u são m e as
  colunas v são n), em inteiros e `Fraction`.
  - μ_RefA = 722873400·2⁻³², como no artigo.
  - Paridades conferidas: campo 0 só com (m par, n ímpar); campos 1 e 2 só com m par; **campo 3 só com
    m ímpar**; v_{0n} reais.
  - Re v_{1,0}(ω4) = 0 (espaço Gauged de RT).
  - **c_1(RefA4(·, 0)) = 43042006046198055061/2⁶⁹ + i·291043177964010739899/2⁶⁹ ≈ 0,0729160 + 0,4930463i**,
    **idêntico**, fração por fração, ao `build/t2/gauge_global_centro.json` do autor.
- **Convenções no artigo** (§"Fourier-Chebyshev series", (eq:operators)).
  - TF(v) = Σ v_mn e^{imτ/2 + inθ}, com θ = arccos ξ. ξ = 0 ⇔ θ = π/2 ⇔ e^{inθ} = iⁿ.
  - O espaço V impõe v_mn = conj(v_{−m,−n}) e v_mn = v_{m,−n}.
  - A fórmula c_m = v_{m0} + 2Σ_{L>=1}(−1)^L v_{m,2L} coincide com a identidade (eq;sidaux3) de RT
    para TF(v)|_{ξ→0}, que é (B_{0,δ0n} − B_{2,δ0n})BI_{2,1}v. Isso confere a convenção independentemente
    de θ.
  - O coeficiente de Fourier de e^{iτ/2} obtido por FFT numa malha de 64 pontos reproduz c_1.
- **Índice 3 = ω4.**
  - SSPACE = (TP ∩ XiOdd) × TP × TP × TAP, com ω4 na 4ª componente.
  - O código de RT aplica `Field_set_VGauged((r->comp)+3)`.
  - Os dados têm só m ímpar no campo 3.
- **Bola.**
  - (kdkkkskkskksksks): ‖ω − (RefA + RefB)‖_SSPACE <= 2⁻²⁷⁷, isto é, a bola **é** para
    ω* − (RefA + RefB).
  - (dkjfhk2) com 𝓛2 = 2⁻²⁵: ‖RefB‖_SSPACE <= 2⁻²⁵, **na mesma norma**.
  - A norma é Σ_{m,n∈ℤ} κ1^{|m|}κ2^{|n|}‖v_mn‖_C, com **pesos >= 1** e ‖z‖_C = |Re z| + |Im z| >= |z|.
  - Logo |c_1(ω*) − c_1(RefA)| <= Σ_n |δv_{1n}| <= ‖δω4‖_V <= 2⁻²⁵ + 2⁻²⁷⁷ ≈ 2,98·10⁻⁸ < 3·10⁻⁸.
- **Conferência de 𝓛2 (exata).** Calculei a norma SSPACE de RefB = RefAplusB − RefA, nas 4
  componentes, 450 × 1350 modos, em inteiros exatos: **‖RefB‖ = 2,1686·10⁻⁸ = 0,7277·2⁻²⁵ <= 2⁻²⁵**
  (`f1_refB_saida.txt`).
- **Rota que dispensa 𝓛2.** c_1(RefAplusB) foi calculado exatamente, e
  |c_1(RefAplusB) − c_1(RefA)| = 1,9·10⁻²¹. Então **|c_1(ω*)| >= |Im c_1(RefAplusB)| − 2⁻²⁷⁷ =
  0,4930463216…** Pela rota do autor: 0,4930462918 > 0.
- **Erro de documentação.** §4 de `GAUGE_GLOBAL.md` diz que "o script confere que o campo de índice 3
  de RefA.dat só tem m ímpar". `gauge_global_centro.py` **não tem** essa conferência: só faz
  `assert 0 <= i < num_m`. O fato é verdadeiro (meu script confere). → R4.

### 2.8 G8 (Passos 6–7)

- **Passo 6.**
  - δω_{(Θ²)*X} = δω_X∘Θ² = e^{4πs}δω_X = δω_{e^{4πs}X}, logo δω de (Θ²)*X − e^{4πs}X é 0.
  - Esse campo é C¹ e admissível, da forma f̂(u_±)∂_± com f̂(u) = Λ⁻¹f(Λu) − e^{4πs}f(u).
  - Pelo Passo 5, f̂ ≡ 0, isto é, f(Λu) = Λe^{4πs}f(u) = ρf(u), com ρ = e^{4π(s−μ)}, para todo u ∈ ℝ.
- **Regimes.** Com Λ = e^{−4πμ} ∈ (0, 1) e |ρ| = e^{4π(Re s − μ)}:
  - Re s > 0 ⇔ |ρ| > Λ;
  - Re s > μ ⇔ |ρ| > 1;
  - Re s = μ ⇔ |ρ| = 1;
  - 0 < Re s < μ ⇔ Λ < |ρ| < 1.

  Os três casos cobrem Re s > 0 sem buraco.
  - **|ρ| > 1:** ρ^j f(u) = f(Λ^j u) → f(0), que é limitado. Logo f(u) = 0. Só usa continuidade em 0.
  - **|ρ| = 1, ρ ≠ 1:** se f(u) ≠ 0, então ρ^j = f(Λ^j u)/f(u) converge. Mas
    |ρ^{j+1} − ρ^j| = |ρ − 1| > 0, então a sequência não é de Cauchy. Contradição, logo f(u) = 0. Só
    usa continuidade.
  - **ρ = 1:** f(u) = f(Λ^j u) → f(0), logo f é constante. Só usa continuidade.
  - **Λ < |ρ| < 1:**
    - f(0) = lim ρ^j f(u) = 0, ou diretamente f(0)(1 − ρ) = 0.
    - f é diferenciável em 0, e para u ≠ 0, f(Λ^j u)/(Λ^j u) = (ρ/Λ)^j f(u)/u → f′(0), finito.
    - Como |ρ/Λ| > 1, f(u) = 0. Usa diferenciabilidade em 0.
- **ρ = 1 ⇔ e^{4π(s−μ)} = 1 ⇔ 4π(s − μ) ∈ 2πiℤ ⇔ s ≡ μ mod (i/2)ℤ.** Para as raízes de T2,
  Re s ∈ {0,0414; 0,401; 0,733} ≠ μ, logo |ρ| ≠ 1 e ρ ≠ 1 sem depender de Im s.
- **Setores.** O argumento só usa Θ² e não vê o setor. O único modo, c·g*, tem perfil de paridade A
  (Z preserva paridades em τ e ξ). Seus aliases com perfil B ficam fora das faixas de T2 (§2.1).
- **Nitidez.** f(u) = |u|^a, com a = 1 − s/μ, satisfaz f(Λu) = ρf(u).
  - Para 0 < Re s < μ, temos 0 < Re a < 1: f é contínua, mas não diferenciável em 0. O C¹ é
    necessário nesse regime.
  - Para Re s >= μ, Re a <= 0 e a ≠ 0: f não é contínua em 0. A continuidade basta.

  Isso confirma a Obs. 4.

### 2.9 G9 (Passo 8 e consequências)

- **f ≡ 1.**
  - X0 = ∂_- + ∂_+ = e^{μτ}(μ⁻¹∂_τ + ξ∂_ξ) ([C]).
  - Derivei δζ, δW e δQ de L_{X0}g ([F]): δg_{±±} = 0, δζ = e^{μτ}(Zζ − 1), δW = X0(W) e
    δQ = e^{μτ}(ZQ + 2Q), que é regular e par. Logo X0 é admissível.
  - Com δφ = e^{μτ}Zφ, as 7 funções de (ekjehee) satisfazem **δω = e^{μτ}(Z + 1)ω**, com ζ, Q e φ
    genéricos ([G], 7/7). Isso confere `GAUGE_RT_ACTION.md`, inclusive [D^σ, X0] = e^{μτ}D^σ ([E]).
  - `test_gauge_action.py`: **3/3 OK**.
  - g* ≠ 0 pela injetividade, porque X0 ≠ 0.
- **Uso em T2.**
  - Passo 11b: "os modos de gauge de Floquet com Re s > 0 são ℂ·g* (s ≡ μ)" é exatamente o enunciado.
  - "Conta":
    - **0,7332:** E é simples. Se e^{sτ}h fosse δω_X, o lema daria e^{4π(s−μ)} = 1, impossível com
      Re s ≠ μ. Exatamente o que o lema prova.
    - **μ\*:** "E = span(g*) é todo de gauge (Lema G)". A parte "g* é gauge" vem da admissibilidade de
      X0, conferida acima. A parte "E = span(g*)" é S4, não revisto aqui. O lema não é necessário
      para essa linha, mas é coerente com ela.
    - **s_K e s_B:** violam as constraints, e o lema é irrelevante para elas.
  - A ressalva de compatibilidade de classes (R2) vale para o uso.
  - Notação: no lema, g* inclui e^{μτ}; em T2 e S4, g* = (Z* + 1)ω* é o perfil. → R6.
- **Observação 1.** Correta:
  - f(u) = (−u)^{1−s/μ} em u < 0 é suave perto do eixo, logo admissível no interior;
  - é Floquet com expoente s;
  - não é nula, pelo argumento do eixo.

  A frase "é a mesma regularidade através do cone que torna o espectro discreto" é heurística, não
  provada e não usada. → R7.
- **Observação 2.** Correta. Cone fixo ⇔ X(u_-)|_{u_-=0} = f(0) = 0, e f constante com f(0) = 0
  é 0. Então não há gauge com Re s > 0. É coerente com T2 ("a contagem passa a 2").
- **Observações 3–5.** Também corretas:
  - Obs. 3: δω = 0 ⇒ δQ = 0, δζ e δφ constantes, por (dhjhfkfhfe1–4);
  - Obs. 4: §2.8;
  - Obs. 5: o lema não prova a ida de H-rec.

## 3. Lacunas e ressalvas (gravidade e conserto)

Nenhuma é ERRO, e nenhuma invalida o Lema G.

| # | gravidade | descrição | conserto proposto |
|---|---|---|---|
| R2 | **média** (para a cadeia de T2, não para o lema) | O lema só vale para gauges C¹ em M × S², incluindo eixo e cone. H-rec (ida) fala em "gauge admissível" sem fixar a regularidade. Se o gauge de fixação de H-rec for menos regular, um puro gauge físico, depois de fixado, pode ser δω_Y com Y fora da classe do lema. | Declarar em H-rec (`T2_ENUNCIADO.md`, `T2_HIPOTESES.md`) que o gauge de fixação é C¹ no sentido do Lema G. Ou adotar o reforço do §2.1: X contínuo em M, C¹ fora do cone e δω_X na classe analítica de T2, com F1′ (ω4(·, −1) ≢ 0, \|c_1\| >= 0,0977). |
| R1 | baixa | "μ fixo e coordenadas (τ, x) fixas" não está na definição de "admissível". É essencial: com δμ ≠ 0, δg_{±±} ≠ 0 nos u_± antigos. | Acrescentar a §1 do lema: "na forma linearizada do ansatz com μ = μ* e coordenadas fixas". |
| R3 | baixa | Passo 4 por fluxos: Ω(ψ_ε*(g*, φ*)) não está definido para ε ≠ 0, e o fluxo não preserva M. | Usar a versão linear do §2.5 (naturalidade de DΩ sob Θ + invariância conforme + (Θ²)*L_X = e^{−4K}L_{(Θ²)*X}). |
| R5 | cosmética | O Passo 1 diz que a forma do ansatz "exige" δg_{±±} = 0. A forma também exige regularidade de δQ no eixo, o que com f só C¹ pode falhar (E_f ~ f‴d²/12). | Dizer que δg_{±±} = 0 é condição **necessária** e que só ela é usada. |
| R4 | cosmética | §4 afirma que `gauge_global_centro.py` confere "m ímpar no campo 3"; não confere. | Adicionar o `assert` ao script, ou corrigir o texto. O fato é verdadeiro (`f1_proprio.py`). |
| R6 | cosmética | g* = e^{μτ}(Z + 1)ω* no lema; g* = (Z* + 1)ω* em T2 e S4. | Uniformizar: "δω_{X0} = e^{μτ}g*". |
| R7 | cosmética | Obs. 1, "torna o espectro discreto", é heurística. | Marcar como comentário não provado. |

## 4. O que não foi verificado

- **O teorema de RT** (existência e bola 2⁻²⁷⁷) é tomado como publicado. Dos dados, só conferi 𝓛2
  (‖RefB‖ <= 2⁻²⁵, exato), μ_RefA e c_1 de RefAplusB. A contração e os arquivos INV/HPC não.
- **S4** (E = span(g*) em μ*, simplicidade) e as localizações e simplicidades das outras raízes (NK,
  Rouché): fora do escopo.
- **H-rec (ida), H-cone e H-reg** são hipóteses de T2. Não li `docs/HREC_IDA.md`, que apareceu durante
  a revisão e não está na lista permitida.
- **A analiticidade de ω\*** na faixa |Im τ| < 2 log κ1 × elipse e a reconstrução 2D→4D (paridades
  Pζ = ζ, Pφ = φ, PQ = Q) são tomadas de RT. Li as passagens, mas não refiz as provas.
- **O reforço do §2.1** é um esboço com F1′ calculado. Não é parte do lema e não foi revisado por
  terceiros.
- **`gauge_global_centro.py` não foi executado**, porque grava em `build/t2/`. Comparei o JSON
  existente.

## 5. Artefatos (nesta pasta)

| arquivo | conteúdo |
|---|---|
| `f1_proprio.py`, `f1_proprio.json`, `f1_proprio_saida.txt`, `f1_proprio_sem_refB.json` | F1 com código próprio: c_1 exato, paridades, Gauged, FFT, c_1 em ξ = −1 (F1′), c_1(RefAplusB), cotas inferiores |
| `f1_refB_saida.txt`, `f1_refB_log.txt` | Norma SSPACE exata de RefB = 2,1686·10⁻⁸ <= 2⁻²⁵ |
| `g_simbolico.py`, `g_simbolico_saida.txt` | **67/67** conferências simbólicas [A]–[H] |
| `g5_equivariancia_linear.py`, `g5_equivariancia_linear_saida.txt` | Equivariância linear: 3/3 |
| `g5_perfis_nulos.py`, `g5_perfis_nulos_saida.txt` | `null_generator` vs pushforward direto: 21/21 |
| `test_gauge_action_saida.txt` | `scripts/test_gauge_action.py`: 3/3 OK |
| `manifesto_inicio.txt`, `manifesto_fim.txt` | Conferência do manifesto |
| `rt/` | Artigo de RT extraído |

## 6. Conclusão

O Lema G está correto. A prova global (Passos 0–8) fecha:
- com X C¹ (contínuo no cone basta para Re s >= μ);
- sem analiticidade dos germes;
- com o único fato numérico F1, que reproduzi com código próprio em racionais exatos, com margem de
  0,49 contra um erro de 3·10⁻⁸.

O uso em T2 (11b e "Conta") é exatamente o que o lema prova. A ressalva que importa para a cadeia de
T2 é R2: alinhar a classe de regularidade do gauge de H-rec com a do lema, ou usar o reforço do §2.1.

## 7. Manifesto (fim)

- 11:27: `git rev-parse HEAD` = `b31da1cb6adea86e543d0c03a60c9a468cb65361`, igual a `COMMIT.txt`.
- `sha256sum -c MANIFEST.sha256`: **238/238 SUCESSO, 0 falhas** (`manifesto_fim.txt`).
- Nenhum arquivo rastreado em `docs/`, `scripts/` ou `build/t2` foi modificado (`git status`). Os três
  arquivos novos citados no §0 não foram criados por esta revisão.
- Tempo total: cerca de 20 min de relógio (11:07–11:27).
