# Revisão independente da volta da reconstrução (03/10/2026)

Objeto: `docs/GEOMETRIC_FLOQUET_RECONSTRUCTION.md` (a "volta").
Revisor: agente independente (Claude), com artefatos nesta pasta. Só o notebook, sem ssh,
`OPENBLAS_NUM_THREADS=2`. O sympy foi instalado no scratchpad da sessão (`PYTHONPATH`), não no `.venv`.
O artigo de RT foi extraído no scratchpad. Nenhum arquivo existente foi alterado.

**Resultado em uma linha.** A volta **não quebrou**. As identidades centrais (compatibilidade, log W,
fórmulas de curvatura de RT) foram provadas aqui com funções **genéricas**, fora das soluções. O L e o C
dos certificados coincidem ponto a ponto com S δΩ⁺ e S♯ δΩ⁺ transcritos do TeX. As lacunas são de
**ligação e de enunciado**, não de matemática:
- o documento não deriva V1 (de L h = 0 e C h = 0 para δΩ = 0) nem diz em que espaço h vive;
- o status do documento está desatualizado;
- as cadeias pedem as constraints de jato.

Há uma observação fora do escopo que condiciona V8: de onde vem C(s)h = 0 em 0,7332, com perda de raio.

## Tabela de veredictos

| | alegação | veredicto | artefato | o que foi feito |
|---|---|---|---|---|
| V0 | enunciado | **OK com ressalva** (L2, L3) | §V0 abaixo | enunciado preciso reconstruído; hipóteses escondidas listadas; o documento não fixa o espaço de h nem a classe da perturbação |
| V1 | seleção + constraints ⇒ os dez δΩ nulos | **OK com ressalva** (L1: ausente do documento) | `v123_identidades.py` (I5, I6, I7), `v1_L_C_pontual.py` | prova escrita; anulamento em ξ = 0 genérico; L e C do código = S δΩ⁺, S♯ δΩ⁺ do TeX a 10⁻¹⁵, com controles negativos |
| V2 | potenciais e compatibilidade | **OK** | `v123_identidades.py` (I3, I4) | identidades (3,3), (4,4), (2,2) provadas com fundo e perturbação genéricos (até h⁺, h⁻ independentes); controles com pares errados falham; T_s⁻¹ conferido |
| V3 | q, z, p e log W | **OK** | `v123_identidades.py` (I1, I2), `v3_reconstrucao_numerica.py` | identidade algébrica de log W genérica, σ = ±; linearizada e dividida por W; recuperação de h1..h4 derivada; ilustração numérica discriminante |
| V4 | regularidade | **OK** | §V4, `v4_W_refA.py` | paridades derivadas; centro removível; cone sem integral radial; elipse 41/40 conferida; W* >= 0,1683 na malha (não rigoroso) |
| V5 | métrica e Floquet | **OK** | `v7_cadeias.py` (V5) | δg⁻¹ do documento = derivada de g⁻¹; (Θ²)*δg = e^{4πs}e^{−4K}δg conferido simbolicamente |
| V6 | equações geométricas | **OK** | `v6_ricci_rt.py`, `testes_unittest.txt` | as fórmulas (dfkhdjhsdshkfd) de RT conferidas com o Ricci 4D **direto da métrica**, com ζ, Q, φ genéricos (não jatos); controle negativo; linearização e centro escritos |
| V7 | unicidade e cadeias | **OK com ressalva** (L4) | `v7_cadeias.py` | núcleo trivial; constantes; recorrência de jatos e compatibilidade de ordem 1 genéricas; termo secular = derivada em s |
| V8 | consequência para T2 | **OK com ressalva** (O1, L4) | §V8 | dá perturbação física, não nula, não de gauge (Lema G, s ≢ μ), na classe H-reg; correspondência bijetiva para autovetores; C(s)h = 0 vem de fora da volta |

## Conferência do manifesto

- **Início (12:23):** `sha256sum -c .snapshots/2026-10-03-revisao-V/MANIFEST.sha256` → 241/241 OK
  (`manifest_inicio.txt`). `COMMIT.txt` = `8597335479861129602c3fc641cac79f15425632` = HEAD.
- **Fim (12:47):** 241/241 OK (`manifest_fim.txt`). HEAD inalterado (`8597335`); nenhum arquivo rastreado modificado
  (`git status --porcelain` sem entradas além de não rastreados).

---

## V0. O enunciado

O documento não enuncia a volta como teorema. Ele fixa "o fundo exato RT, μ fixo, W > 0 no domínio real
e as normalizações/paridades RT" e **toma como hipótese** que "todos os resíduos δΩ^σ" são nulos. Não
menciona L(s), C(s), D(Y+), o toro nem a definição de modo físico de `HREC_IDA.md`. O enunciado abaixo
é o que os argumentos do documento, completados nesta revisão, de fato provam.

**Volta (forma precisa).** Seja s ∈ ℂ com Re s > 0. Basta, na verdade, s ∉ (i/2)ℤ; Re s > 0 é usado
também para excluir as constantes. Seja h = (h1, h2, h3, h4) tal que:
- h é 4π-periódico em τ, com h1 ímpar em ξ;
- h está em D(Y+), complexificado, setor A, B ou soma;
- L(s)h = 0 e C(s)h = 0, onde L(s)h = S·δΩ⁺[h] e C(s)h = S♯·δΩ⁺[h], com
  δΩ⁺[h] := e^{−sτ} (d/dε) Ω⁺(μ, ω* + εe^{sτ}h)|_{ε=0} no fundo exato de RT.

Então existem q, z, p únicos, 4π-periódicos, pares em ξ e real-analíticos em ℝ × (−41/40, 41/40), tais
que a perturbação no ansatz (δQ, δζ, δφ) = e^{sτ}(q, z, p) tem estas propriedades:
1. Suas variáveis de RT, pelas definições (ekjehee) linearizadas, são δω = e^{sτ}h.
2. A perturbação métrica δg, dada pelo frame (δe0 = 0, δe_i como no documento), e δφ são
   real-analíticas em M̂ = ℝ × B_{41/40}(0), em coordenadas cartesianas, o que inclui o eixo e o cone
   ξ = 1.
3. (δg, δφ) resolve as equações de Einstein–campo escalar linearizadas em torno de (ḡ, φ*), em M̂.
4. (Θ²)*δg = e^{4πs}e^{−4K}δg e (Θ²)*δφ = e^{4πs}δφ.

Em suma, δ é um "modo físico de Floquet" no sentido de `HREC_IDA.md` §1. A aplicação h ↦ δ é linear e
injetiva. Para cadeias, vale o mesmo com e^{sτ}Σ(τʲ/j!)·(perfis), **desde que** os resíduos de **jato**
se anulem (ver L4).

**Onde h vive.**
- A volta só usa que TF(h) é analítico numa vizinhança complexa de ℝ × (−41/40, 41/40), periódico em τ.
- D(Y+) (ℓ¹, raios 65/64 e 5/4) dá isso: a série de Chebyshev converge absolutamente na elipse fechada
  E_{5/4}, de semieixo (5/4 + 4/5)/2 = 41/40, e a de Fourier converge na faixa |Im τ| < 2 log(65/64).
- O toro com sinal (ℓ², mesmos raios) também dá isso. Pela desigualdade de Cauchy–Schwarz, ℓ² com peso
  ρⁿ está contido em ℓ¹ com peso ρ'ⁿ para todo ρ' < ρ, e portanto TF(h) é analítico na elipse **aberta**
  E_{5/4}.
- O lema L-real coloca de todo modo os vetores do núcleo no toro em D(Y+), para |s| <= 3,1, o que
  inclui 0,7332.

**Hipóteses escondidas** (todas satisfeitas, mas não enunciadas no documento):
- (a) **Fundo exato:** Ω(ω*) = 0 (E1). Ele é usado para linearizar a identidade de log W sem o termo de
  resíduo e para matar E*(δe, e) em V6.
- (b) **W* > 0 no cilindro real |ξ| < 41/40:** é a condição (d) de 2Dpre, que RT conferem com RefA.
- (c) **As fórmulas de curvatura (dfkhdjhsdshkfd):** RT as enunciam sem prova. Foram conferidas aqui com
  funções genéricas (V6).
- (d) **μ e as coordenadas (τ, x) fixos:** sem δμ.
- (e) **Para cadeias**, as constraints de jato, e não C(s0)h_j = 0 (L4).
- (f) **"Físico"** quer dizer: na classe H-reg, real-analítica em M̂, e Floquet exato com o fator
  e^{−4K}.

**O sentido de Floquet.**
- ζ*∘Θ = ζ* + K, logo ζ*∘Θ² = ζ* + 2K. Q* é periódico, e o frame tem coeficientes independentes de τ.
- Pela fórmula do documento, δg⁻¹ = e^{sτ}e^{2ζ*}[...]. Então δg = −ḡ δg⁻¹ ḡ = e^{sτ}e^{−2ζ*}(ĝ[...]ĝ),
  com ĝ = e^{2ζ*}ḡ invariante.
- Daí (Θ²)*δg = e^{4πs}e^{−4K}δg e (Θ²)*δφ = e^{4πs}δφ, que é exatamente a definição de `HREC_IDA.md`.
  Conferido simbolicamente em `v7_cadeias.py` (V5).

## V1. Seleção + constraints ⇒ os dez δΩ^σ nulos

**O documento não prova isto.** Ele assume "todos os resíduos δΩ^σ nulos". A ligação está na Remark de
RT após 2DprobSeries ("SΩ⁺ = 0 e S♯Ω⁺ = 0 ⟺ Ω⁺ = 0") e em `SHARP_LINEARIZED_IDENTITY.md`
("e = Ic + RF"), que o documento não cita. Prova, escrita para a linearização:

1. **Definições.** S = Ξ⁻¹∘(½(1+P) na comp. 1, comps. 2, 4, 5) e S♯ = (½(1−P) na comp. 1, −P na
   comp. 3), com Ξ⁻¹f = (f − f|_{ξ=0})/ξ (RT, eq. sid4). As componentes de Ω⁺ são
   (Ω1, Ω2, Ω2♯, Ω3, Ω4)⁺.
   - L(s)h = S δΩ⁺[h] e C(s)h = S♯ δΩ⁺[h]. Conferido **ponto a ponto** contra o código dos
     certificados (`signed_operator_L`, `signed_constraints`), com fundo RefA, h aleatório do setor A
     e s = 0,7 + 0,3i: diferença relativa <= 10⁻¹⁵ nas 4 + 2 componentes.
   - Os controles negativos (s + 0,01; um coeficiente de Ω3 alterado) dão discrepâncias de
     1,5·10⁻¹ e 1,3 (`v1_L_C_pontual_saida.txt`).
   - A parte livre de C é exatamente μ(ξ∂ξ + 2)h1 e −ξ(∂τ + s + μ)h2 + μ(ξ − ξ²)∂ξh2 + 2μ(h2 − h3),
     como em `signed_constraints.py` (I7, simbólico).
2. **Anulamento em ξ = 0.** As componentes selecionadas ½(1+P)δΩ1⁺, δΩ2⁺, δΩ3⁺ e δΩ4⁺ se anulam em
   ξ = 0 para todo h com h1 ímpar e todo fundo codificado. Isso foi conferido com funções genéricas
   (I6).
   - O motivo: os termos principal e quadrático têm fator ξ, e as componentes 1, 2, 4, 5 de Γ1 são
     ímpares.
   - Controle: a parte ímpar de δΩ1⁺ e δΩ2♯⁺ **não** se anulam em ξ = 0.
   - Logo Ξ⁻¹f = 0 ⇒ f(τ, ξ) = f(τ, 0) = 0. **Nada se perde na divisão regularizada.**
3. **De L e C a δΩ⁺ = 0.**
   - L(s)h = 0 ⇒ ½(1+P)δΩ1⁺ = δΩ2⁺ = δΩ3⁺ = δΩ4⁺ = 0.
   - C(s)h = 0 ⇒ ½(1−P)δΩ1⁺ = 0 e PδΩ2♯⁺ = 0, isto é, δΩ2♯⁺ = 0.
   - Somando as partes par e ímpar, δΩ1⁺ = 0. Logo os cinco δΩ⁺ são nulos.
4. **Os sinais −.** Com a codificação ω_i⁺ = −Pω_i e h_i⁺ = −Ph_i (i = 2, 3, 4), e ω1, h1 ímpares,
   vale PΩ_i^σ = −Ω_i^{−σ} e PδΩ_i^σ = −δΩ_i^{−σ}. Isso foi conferido com funções genéricas, não
   linear e linearizado, para os 5 × 2 casos (I5); o fator e^{sτ} comuta com P. Logo δΩ⁻ = −PδΩ⁺ = 0.
   **Os dez δΩ^σ são nulos.**
5. **Projeções de paridade.** O setor A ou B é preservado por tudo (as operações são diagonais em m ou
   não mexem em τ). A soma A ⊕ B é o espaço 4π-periódico inteiro. Nada se perde.

**Observação estrutural**, confirmada numericamente:
- δΩ2⁺, δΩ3⁺ e δΩ4⁺ são selecionados por L. Logo as compatibilidades (2,2), (3,3) e (4,4) de V2 já
  seguem de L h = 0, sem C.
- C entra **só** pela identidade de log W, via a parte ímpar de δΩ1 e via δΩ2♯, isto é, na
  recuperação de h2.
- Em `v3_reconstrucao_numerica.py`, os autovetores de Galerkin dão o seguinte:

  | raiz | D_s^± z − h3 | e2 = |D_s⁻(z + ξ²q/W*) − h2| (rel.) | |C(s)h|/|h| |
  |---|---|---|---|
  | s_K (viola C) | ~10⁻¹² | ≈ 1,17, estável nas três truncagens | 0,285 |
  | μ* | ~10⁻¹² | 3,4·10⁻² → 6,1·10⁻³ → 1,2·10⁻³ (8×24 → 10×30 → 12×36) | acompanha e2 |
  | 0,7332 | ~10⁻¹² | 3,6·10⁻² → 5,7·10⁻³ → 7,8·10⁻⁴ (8×24 → 10×30 → 12×36) | acompanha e2 |

  É uma ilustração, não prova. O limiar de 10⁻² do script é arbitrário: em 8×24 ele imprime "PADRAO
  INESPERADO" só porque e2 = 3,5·10⁻²; a tendência é o que importa.

## V2. Potenciais e compatibilidade

- **Equivalência.** D_s^σ z = σT_s z + μ(1 + σξ)∂ξz. A soma das duas equações dá 2μ∂ξz = h⁺ + h⁻,
  isto é, ∂ξz = b. Substituindo, T_s z = h⁺ − (1 + ξ)(h⁺ + h⁻)/2 = a. A recíproca é imediata. Conferido
  em I4, nos dois sentidos.
- **T_s⁻¹.**
  - No modo e^{imτ/2}, T_s⁻¹ multiplica por 1/(s + im/2). A integral ∫₀^∞ e^{−st}a(τ − t)dt dá
    exatamente isso, por integração por partes, com Re s > 0.
  - |1/(s + im/2)| <= 1/Re s. Com ‖z‖_C = |Re z| + |Im z| <= √2|z|, a norma nos DOFs reais é
    <= √2/Re s <= 2/Re s. A cota do documento está certa, embora folgada.
  - T_s⁻¹ é diagonal em m, logo preserva o setor e a paridade temporal. Não age em ξ, logo preserva a
    paridade radial. É holomorfa em s para Re s > 0 (de fato, para s ∉ −(i/2)ℤ).
  - **Ausência de obstrução de período.** A primitiva F de uma 1-forma fechada e^{sτ}×(periódica), na
    cobertura, satisfaz F∘Θ² = e^{4πs}F + c. Então F + c/(e^{4πs} − 1) é Floquet exata, e é única,
    porque e^{4πs} ≠ 1. Isso equivale à fórmula de T_s⁻¹ e mostra por que s = 0 é excluído.
  - Se T_s b = ∂ξa, então ∂ξ(T_s⁻¹a) = T_s⁻¹∂ξa = T_s⁻¹T_s b = b, porque T_s⁻¹ comuta com ∂ξ e T_s é
    injetivo nas periódicas.
- **Identidades de compatibilidade.** −2μξ(T_s b_i − ∂ξa_i) = δΩ_j⁻ − δΩ_j⁺ para
  (i, j) = (3,3), (4,4), (2,2).
  - **Provadas com sympy** (I3) com fundo ω **genérico** (sem Ω(ω) = 0) e perturbação genérica, até
    com h_i⁺ e h_i⁻ independentes. São, portanto, identidades **fora das soluções**, como o documento
    afirma.
  - Derivação em papel: os termos μ e quadráticos de Ω2, Ω3 e Ω4 não dependem de σ e cancelam na
    diferença. Sobra ξ[(D_s⁻ − μ)h⁺ − (D_s⁺ + μ)h⁻] = −2μξ(T_s b − ∂ξa), usando
    μ(1−ξ)∂ξh⁺ − μ(1+ξ)∂ξh⁻ = 2μ∂ξa + μ(h⁺ + h⁻). É a mesma conta do fecho de α em RT ("From 2D to
    4D").
  - Controles: os pares trocados (4,3), (2,4) e (3,2) falham, δΩ1⁻ − δΩ1⁺ não se anula, e um Ω3
    alterado quebra a identidade (3,3).
- **Aviso sobre os testes do autor.** `test_closed_form_residuals_without_on_shell_assumption` confere
  as três identidades numa **única instância** polinomial (com as tabelas C de RT), o que não é prova.
  A frase "verificadas diretamente com polinômios racionais" deve ser lida assim. A prova genérica está
  em I3.

## V3. Reconstrução e log W

- **A identidade algébrica.** Com W := μ + ξ(ω1 − ω2⁻ − ω2⁺)/2, isto é, Q dado por (dhjhfkfhfe1), vale
  D^σW − W(ω2^σ − ω3^σ) = (Ω1^σ − Ω2^σ − Ω2♯^σ)/2.
  - **Provada com funções genéricas, para σ = + e σ = −, com ω⁺ e ω⁻ independentes** (I1).
  - O teste do autor só confere σ = + numa instância.
- **A linearização.** Com δW = ξ(h1 − h2⁻ − h2⁺)/2 = ξ²q:
  - versão sem dividir (I2): D_s^σδW − δW(ω2^σ − ω3^σ) − W(h2^σ − h3^σ) = (δΩ1 − δΩ2 − δΩ2♯)^σ/2;
  - versão dividida por W (I2b): vale com o termo −(D^σW − W(ω2^σ − ω3^σ))δW/W² que o documento
    menciona. No fundo exato esse termo é zero, e D_s^σ(ξ²q/W*) − (h2^σ − h3^σ) =
    (δΩ1 − δΩ2 − δΩ2♯)^σ/(2W*). Confere com o documento.
- **Recuperação de h1, h2, h3, h4**, dados os dez δΩ nulos (V1):
  1. h3^σ = D_s^σz, por V2 com (3,3); h4^σ = D_s^σp, com (4,4).
  2. h2^σ = h3^σ + D_s^σ(ξ²q/W*) = D_s^σ(z + ξ²q/W*), pela identidade de log W linearizada com
     δΩ1 − δΩ2 − δΩ2♯ = 0. **Aqui, e só aqui, C é necessário.**
  3. h1 = 2ξq + h2⁻ + h2⁺, pela definição de q.

  E isso é exatamente δω da perturbação:
  - ω1 = 2ξQ + ω2⁻ + ω2⁺ (os ±σμ de (ekjehee) cancelam);
  - ω2^σ = D^σ(ζ + log W) − σμ, com δ log W = ξ²δQ/W*;
  - ω3^σ = D^σζ − σμ e ω4^σ = D^σφ;
  - D^σ não depende de (Q, ζ), porque x^k e_k/ξ = μ∂ξ.

  **As fórmulas de q, z, p estão corretas.**
- **Os testes do autor.**
  - `test_reconstruct_definitions_at_constant_W` reconstrói **só com W constante** (fundo Q* = 0) e
    T_s⁻¹ polinomial (série de Neumann finita).
  - O caso com W* não constante só é coberto pela identidade de log W (teste 1).
  - A ilustração numérica de V1 cobre o caso real, com W* de RefA.

## V4. Regularidade

- **Paridade.** h1 é ímpar, e h2 − Ph2 também, logo o numerador de q é ímpar e q é par. a_i é par: com
  h = e + o, a = −e − ξo. b_i = o/μ é ímpar. Como T_s⁻¹ não age em ξ, z e p são pares. Funções pares
  real-analíticas de ξ são funções analíticas de ξ² = |x|², logo analíticas em coordenadas cartesianas.
- **Centro removível.** O numerador é ímpar e analítico, logo se anula em ξ = 0, e Ξ⁻¹ (regularizado)
  coincide com a divisão. q ∈ Y+, porque Ξ⁻¹ é limitado (40/9 nos pesos de RT).
- **Cone.**
  - z e p vêm de T_s⁻¹, que é multiplicação diagonal em m e não envolve ∂ξ nem integração em ξ.
  - q é algébrico, e ξ²q/W* precisa só de W* ≠ 0 real.
  - Portanto **não há integral radial atravessando ξ = 1**. A analiticidade em ξ = 1 é herdada de h.
    Confere.
- **W* > 0.**
  - W* = (2μ + ξ(ω1 − ω2 + Pω2))/2 > 0 em |ξ| < 41/40 é a condição (d) de 2Dpre, que RT conferem com
    RefA e erro <= 2⁻²⁴.
  - Conferência numérica própria, não rigorosa (`v4_W_refA.py`, malha 800 × 821): min W* = 0,168307 = μ,
    em ξ = 0; min W* = 0,4285 em |ξ − 1| < 0,01. A margem é enorme diante de (41/80)·3·2⁻²⁴ ≈ 9·10⁻⁸.
  - As paridades do fundo conferem: ω1 só tem n ímpar, ω1..3 só m par e ω4 só m ímpar.
  - O documento evita, com razão, afirmar 1/W* na elipse complexa: basta analiticidade perto dos
    compactos reais.
- **A elipse.** O semieixo é (5/4 + 4/5)/2 = 41/40. Confere (RT, linha 320: "|ξ| < 41/40", que dá
  u_- < −u_+/81).
  - Em (τ, ξ), M = {u_+ < 0, u_+ < u_- < −u_+/81} é {0 < ξ < 41/40}: de u_+ < u_- sai ξ > 0, e de
    u_- < −u_+/81 sai 80ξ < 82.
  - **Sim:** h ∈ Y+ é real-analítico em todo ℝ × [0, 41/40) e, pela paridade, em M̂ inteiro, com o eixo.
  - A vizinhança complexa afina quando ξ → 41/40 (a elipse tem espessura zero no vértice). Isso basta,
    porque M é aberto.
  - Em ξ = 41/40 há só continuidade, e esse ponto não está em M.

## V5. A métrica

- **O frame.**
  - Do frame de RT: e0 = ∂τ + μξ∂ξ, x^k e_k/ξ = μ∂ξ (os termos em Q cancelam), e t·e = W t·∂ para t
    tangente.
  - Logo g = e^{−2ζ}[−(1 − ξ²)dτ² − (2ξ/μ)dτdξ + dξ²/μ²] + e^{−2ζ}ξ²W⁻²g_{S²}.
  - Ortonormalidade a menos de e^{−2ζ}: conferida em `v6_ricci_rt.py`.
- **δg⁻¹.**
  - δe0 = 0, porque e0 não depende de (Q, ζ).
  - δe_i = q Σ_k x^k(x^k∂_i − x^i∂_k), de ∂e_i/∂Q.
  - A fórmula do documento para δg⁻¹ é exatamente (d/dε)g⁻¹(Q + εe^{sτ}q, ζ + εe^{sτ}z), conferida
    simbolicamente em cartesianas com Q, ζ, q, z genéricos (`v7_cadeias.py`, V5).
  - δg = −ḡδg⁻¹ḡ e δφ = e^{sτ}p.
- **Floquet.** Com ζ* = Z_p + Kτ/(2π), Z_p 4π-periódico, vale (Θ²)*δg = e^{4πs}e^{−4K}δg. Conferido
  simbolicamente no ponto (x1, 0, 0), x1 genérico, o que basta pela simetria esférica. Ver também V0.

## V6. As equações geométricas

- **As fórmulas de RT.** RT enunciam (dfkhdjhsdshkfd) sem derivação. `v6_ricci_rt.py` calcula o Ricci 4D
  **direto da métrica**, com Christoffel, em (τ, ξ, θ, ϕ), com ζ, Q, φ funções **genéricas**. Confere,
  com resíduo exatamente 0:
  - Ric(e_r, e_r) − 2(e_rφ)² = (A + B)/(4ξ) e Ric(e_θ, e_θ) = A/(4ξ), que são as duas partes de
    δ_ij A/(4ξ) + x^ix^j B/(4ξ³);
  - Ric(e0, e0) − 2(e0φ)²;
  - Ric(e0, e_r) − 2e0φ e_rφ, que é a parte x^i/(2ξ²);
  - Ric(e0, e_θ) = Ric(e_r, e_θ) = 0;
  - e^{−2ζ}□φ = (Ω4⁻ + Ω4⁺)/(2ξ);
  - as cinco identidades de definição do Claim de RT, Ω1^σ − Ω2^σ − Ω2♯^σ = 0 e Ω_j⁻ = Ω_j⁺ (j = 2, 3,
    4), para ω dados por (ekjehee).

  Controle negativo: mudar um coeficiente de Ω1 dá resíduo não nulo. Isto é mais forte que a conferência
  por 3 jatos da revisão R4 da ida. Os Ω vêm da minha transcrição do TeX, independente das tabelas C.
- **A linearização** (escrita aqui).
  - As identidades valem para todo (Q, ζ, φ) com W ≠ 0. São identidades polinomiais nos jatos, logo
    valem também para direções complexas.
  - Derive ao longo de (Q* + εe^{sτ}q, ζ* + εe^{sτ}z, φ* + εe^{sτ}p). Com E := Ric − 2dφ⊗dφ:
    δ[E(e_a, e_b)] = (δE)(e_a, e_b) + E*(δe_a, e_b) + E*(e_a, δe_b) = F_ab(δΩ).
  - No fundo exato, E* = 0. Por V3, δω = e^{sτ}h, logo δΩ = e^{sτ}·(δΩ de h) = 0.
  - Então (δE)(e_a, e_b) = 0 para ξ ≠ 0, no frame, que é base porque W* > 0. As demais componentes saem
    pela simetria esférica.
  - Para o campo escalar: δ(e^{−2ζ}□φ) = e^{−2ζ*}δ(□φ) − 2δζ e^{−2ζ*}□φ* = e^{−2ζ*}δ(□φ).
- **O centro.** δg, δφ e o fundo são analíticos em (τ, x), logo δE e δ□φ são contínuos. Como se anulam
  para ξ ≠ 0, anulam-se no eixo.
- **`test_geometric_reconstruction`, 5/5** (`testes_unittest.txt`). O que cada teste confere:
  1. a identidade de log W e sua linearização, numa instância polinomial, só σ = +;
  2. as três compatibilidades, numa instância;
  3. a reconstrução com **W constante** e T_s⁻¹ polinomial;
  4. a recorrência de jatos, com W constante;
  5. a invertibilidade da matriz 5 × 5 das combinações geométricas, que serve à **ida** (R4), não à
     volta.

  Conferi as cinco linhas dessa matriz a partir de (dfkhdjhsdshkfd), com Ω2♯^σ = Ω1^σ − Ω2 e sinais
  iguais: [0,0,4,0,0], [2,2,−4,4,0], [1,1,−4,−2,0], [−1,1,0,0,0], [0,0,0,0,2]. A matriz é invertível,
  por eliminação: r1 dá Ω2, r5 dá Ω4, r2 − 2r3 dá Ω3, e r3, r4 dão Ω1^±. Os testes não conferem V1 nem
  as fórmulas de curvatura. `test_root_constraints`, 3/3.

## V7. Unicidade e cadeias

- **Núcleo trivial.** A aplicação h ↦ (q, z, p) é linear, e h = 0 dá q = 0 e z = p = T_s⁻¹0 = 0.
  Também é injetiva no sentido útil: se δω(q, z, p) = 0, então:
  - T_s z = 0 e z_ξ = 0, logo z = 0, porque z é periódico e e^{4πs} ≠ 1;
  - p = 0, do mesmo modo;
  - D_s^σ(ξ²q/W*) = 0, logo q = 0.
- **As constantes.** δζ = c e δφ = c têm δω = 0, porque ω só vê D^σζ e D^σφ. Escritas como
  e^{sτ}(c e^{−sτ}), o perfil só é 4π-periódico se e^{4πs} = 1, o que Re s > 0 exclui.
- **Jatos.** Com δζ = e^{sτ}(z1 + τz0), D^σ do jato é e^{sτ}(h1 + τh0), com h0^σ = D_s^σz0 e
  h1^σ = D_s^σz1 + σz0. Isso dá exatamente T_s z1 = a(h1) − z0 e ∂ξz1 = b(h1), a recorrência do
  documento. Na convenção τʲ/j!, vale para comprimento qualquer, por indução.
  - A compatibilidade de ordem 1 é T_s b(h1) − ∂ξa(h1) + b(h0) = 0. É o coeficiente de τ⁰ de
    (δΩ3⁻ − δΩ3⁺)/(−2μξ) da perturbação secular. Conferido com funções genéricas.
  - O termo secular de **todos** os δΩ⁺ é dO(h1) + ∂_s dO(h0) + τ dO(h0) (L-afim, genérico).
  - Logo os resíduos do jato se anulam se e só se L(s0)h0 = 0, L(s0)h1 = −h0, C(s0)h0 = 0 e
    **C(s0)h1 + C1h0 = 0**. Esta última é a constraint de jato D h1 = 0 de
    `CONSTRAINT_ROOT_QUOTIENT.md`, e não C(s0)h1 = 0.
  - O documento diz só "uma cadeia h_j" (lacuna L4). q_j continua algébrico. "Por holomorfia em s"
    confere, porque ∂_s T_s⁻¹ = −T_s⁻².
- **O quociente.** O "núcleo da reconstrução trivial" e a injetividade do mapa (δQ, δζ, δφ) ↦ δω, nas
  Floquet com Re s > 0, dão a bijeção {perturbações Floquet no ansatz que resolvem as equações} ↔
  ker L(s) ∩ ker C(s): ida por R4, volta por V1–V6.

## V8. Consequência para T2

Seja h ≠ 0 em ker L(s) ∩ ker C(s), s ≈ 0,7332, com λ* ∈ [0,73316950; 0,73318982]:
1. **Espaço.** h está em D(Y+), pelo L-real (|s| <= 3,1). É real a menos de escala, porque L(s) é real
   para s real e a raiz é simples.
2. **Volta.** Dá δ = (δg, δφ), um modo físico de Floquet na classe H-reg (real-analítico em M̂, que
   contém o eixo e o cone), com δω = e^{sτ}h. δ é **não nulo**, porque δω ≠ 0.
3. **Não é de gauge.** Se δ = L_X(ḡ, φ*), com X de classe C¹ (ou contínuo com δω_X analítico, que é o
   caso aqui), então X é admissível, porque δ está no ansatz. Então δω_X = e^{sτ}h ≠ 0 é Floquet com
   Re s > 0, e o Lema G força e^{4π(s−μ)} = 1, isto é, s ≡ μ mod i/2. Mas Re s ≈ 0,7332 ≠ μ ≈ 0,1683.
   Contradição. O Lema G não supõe X Floquet; ele o deduz (Passo 6).
4. **A correspondência de `HREC_IDA.md` §1** é bijetiva nos autovetores, para todo Re s > 0:
   - **ida∘volta = id:** δ(h) já está no ansatz com δω = e^{sτ}h, e a ida com X = 0 devolve h.
   - **volta∘ida = id módulo gauge:** se δ' = δ + L_X ḡ está no ansatz com δω' = e^{sτ}h, então
     δ' − δ(h) está no ansatz, com δω = 0. Por R6/V7, (δζ, δQ, δφ) são constantes, logo nulas para
     Floquet com Re s > 0. Assim δ' = δ(h).
   - **Direção de gauge:** em s ≡ μ, a reconstrução de e^{(μ−s)τ}g* é L_{X0}(ḡ, φ*), pela unicidade.
     Os quocientes se correspondem.
   - **Ressalva para cadeias (L4):** "preserva o comprimento das cadeias" vale com cadeias de L com
     constraints de jato. Como as quatro raízes são simples, isso não afeta T2.
5. **Ressalva O1 (fora da volta):** a volta precisa de C(s)h = 0. Em T2 isso vem de K injetivo (S3a) e
   de KC = ML, que está em outra parte da cadeia (passos 11–12). Ver O1.

**Conclusão para T2:** com a volta, a contagem 1 é atingida. O vetor em 0,7332 dá uma perturbação
física, não nula e não de gauge, na classe H-reg, condicionado ao vetor estar em ker C, que vem dos
passos 11–12 de T2, não da volta.

---

## Lacunas (com gravidade e conserto)

- **L1 (moderada; de ligação, não de matemática).**
  - **Problema.** O documento não prova V1: assume "todos os δΩ^σ nulos" e não menciona L, C, S, S♯ nem
    Ξ⁻¹. A ligação existe em RT (Remark após 2DprobSeries) e em `SHARP_LINEARIZED_IDENTITY.md`
    (e = Ic + RF), mas não é citada. Sem ela, "a volta" não se aplica literalmente a ker L ∩ ker C.
  - **Conserto.** Um parágrafo "Lema V1" com os passos 1–4 de V1, citando
    `SHARP_LINEARIZED_IDENTITY.md`. Artefatos: `v123_identidades.py` (I5–I7) e `v1_L_C_pontual.py`.
- **L2 (moderada; de enunciado).**
  - **Problema.** O documento não diz em que espaço h vive (D(Y+) ou toro), nem enuncia que a
    perturbação reconstruída é um "modo físico de Floquet" no sentido de `HREC_IDA.md`: real-analítica
    em M̂, com (Θ²)*δg = e^{4πs}e^{−4K}δg. As conclusões são verdadeiras, mas estão espalhadas ou
    implícitas.
  - **Conserto.** Usar o enunciado de V0 acima, com as hipóteses (a)–(f).
- **L3 (menor; status desatualizado).**
  - **Problema.** O cabeçalho diz que "a identificação com o domínio espectral RT e todas as condições
    físicas de admissibilidade continuam pendentes", e o fim diz que "a contagem física permanece
    aberta". O §"Pendência geométrica precisa" descreve como aberta a fixação de gauge, que hoje é a
    ida (`HREC_IDA.md`). Isso contradiz `T2_ENUNCIADO.md`, que cita o documento como prova da volta.
  - **Conserto.** Atualizar o status. Separar o que é a volta (este documento) do que é a ida
    (`HREC_IDA.md`, Lema R).
- **L4 (menor; cadeias).**
  - **Problema.** Para cadeias, a volta exige que os resíduos do **jato** se anulem, isto é, a cadeia
    de L com C(s0)h0 = 0 e C(s0)h_j + C1h_{j−1} = 0. Não basta C(s0)h_j = 0 em cada vetor. O documento
    diz "uma cadeia h_j" sem condição. `HREC_IDA.md` §1 fala em "(ker L ∩ ker C) … preserva o
    comprimento das cadeias".
  - **Efeito.** Nenhum em T2, porque as raízes são simples.
  - **Conserto.** Escrever a condição de jato, como em `CONSTRAINT_ROOT_QUOTIENT.md`.
- **L5 (menor; natureza dos testes).**
  - **Problema.** Os testes do autor conferem identidades numa instância polinomial, com σ = + e W
    constante na reconstrução. Não são prova.
  - **Conserto.** Incorporar a verificação genérica com sympy (como `v123_identidades.py`,
    `v6_ricci_rt.py` e `v7_cadeias.py`) ou registrar que os testes são de regressão.
- **L6 (menor; redação).**
  - **Problema.** As constantes "2/(Re s)" e "2(1+b)/(Re s)" (com b não definido) são cotas corretas,
    mas não derivadas no texto. Não são usadas em T2.
  - **Conserto.** Escrever |1/(s + im/2)| <= 1/Re s, ‖·‖_C <= √2|·| e ‖ξ‖ <= b.

**Observação fora do escopo, que condiciona V8:**
- **O1 (a conferir; não é lacuna da volta).** C(s)h = 0 em 0,7332 vem de K injetivo e de KC = ML.
  - Por `SHARP_DOMAIN.md`, KC = ML vale com **perda de raio**: C(s)h ∈ D_K(a', b'), com a' < 65/64 e
    b' < 5/4. O próprio documento diz "Não provamos … injetividade sharp nos novos pesos".
  - O L-real (`T2_HIPOTESES.md` §3.1) transfere a injetividade de K do toro para Y+ "pela inclusão".
    Mas Y♯(a', b') **não** está contido no toro nem em Y♯+: c_n = b'⁻ⁿ/n² está em ℓ¹(b'ⁿ), mas não em
    ℓ²(κ2^{2n}).
  - Falta, nos documentos que pude ler, uma recuperação de raio para ker K(s), análoga a R5 para L
    (Q♯ = J♯⁻¹; existência forte e unicidade fraca). O argumento deve ser o mesmo, e
    `RADIUS_RECOVERY.md` menciona o "bloco K" no free_tail, mas não o encontrei escrito.
  - **Conserto:** escrever esse lema "R5♯", ou citar onde está (talvez `S3_CONDICIONAMENTO_SHARP.md`, que
    eu não podia ler).

## O que não foi verificado

- O rigor de W* > 0, a condição (d) de RT: só uma conferência numérica de margem (não rigorosa). A
  afirmação rigorosa é de RT (E1).
- A origem de C(s)h = 0 em 0,7332 (O1): a injetividade de K, KC = ML com perda de raio e S3a estão fora
  do escopo.
- O lema L-real, a ida (Lema R), o Lema G e os certificados (NK em 0,7332, Rouchés), que têm revisões
  próprias. Usei só os enunciados.
- As constantes de norma do documento (2/(Re s) etc.): conferi apenas que são cotas válidas, não seu uso.
- A ilustração numérica (`v3_reconstrucao_numerica.py`) é de Galerkin, com fundo RefA, e não é prova.
  Ela mostra que o papel de C na volta é a recuperação de h2.
- **Transparência.** A primeira execução de `v7_cadeias.py` teve falhas causadas por erros **meus** no
  script: `coeff(tau, 0)` em expressões com funções de τ, e uma chave de substituição que não batia
  depois de x2 = x3 = 0. Corrigido, tudo confere (`v7_cadeias_saida.txt`).

## Artefatos (nesta pasta)

Todos rodam da raiz do repositório. Os que usam sympy precisam de `PYTHONPATH=<sympy>`.

| script | saída | conteúdo |
|---|---|---|
| `v123_identidades.py` | `v123_identidades_saida.txt` | I1–I7 e controles: log W, compatibilidades, equivalência dos potenciais, simetria P, anulamento em ξ = 0, parte livre de C |
| `v1_L_C_pontual.py` | `v1_L_C_pontual_saida.txt` | L(s) e C(s) do código = S δΩ⁺ e S♯ δΩ⁺ do TeX, ponto a ponto, mais controles |
| `v6_ricci_rt.py` | `v6_ricci_rt_saida.txt` | fórmulas de curvatura de RT = Ricci 4D direto, com funções genéricas |
| `v7_cadeias.py` | `v7_cadeias_saida.txt` | jatos e cadeias; δg⁻¹; Floquet com e^{−4K} |
| `v4_W_refA.py` | `v4_W_refA_saida.txt` | W* > 0 e paridades do fundo (numérico) |
| `v3_reconstrucao_numerica.py` | `v3_reconstrucao_numerica_saida.txt` | reconstrução nos autovetores de Galerkin (μ*, s_K, 0,7332) |
| — | `testes_unittest.txt` | testes do autor: 5/5 e 3/3 |
| — | `manifest_inicio.txt`, `manifest_fim.txt` | conferência do snapshot |
