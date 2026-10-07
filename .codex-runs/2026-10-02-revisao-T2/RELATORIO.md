# Revisão independente de T2 (02/10/2026)

Revisor: agente Claude (Fable 5.1), sem acesso aos documentos proibidos pela TAREFA.

Snapshot: `.snapshots/2026-10-02-revisao-T2/MANIFEST.sha256` tem 445 arquivos, e todos batem com a árvore de
trabalho no início da revisão (commit c25d090). Os artefatos desta revisão estão em `art/` (scripts e saídas). No PC
do João: `~/revisao_T2/` (cópia de `scripts/`, saída em `~/revisao_T2/out/`).

**Resumo.** Não consegui quebrar o L-real, a redução B → A, a álgebra do quociente nem a identidade KC = ML. Em todas
as tentativas os números e as derivações fecharam, com folga. Achei 11 problemas (§3), nenhum de gravidade alta:
- 3 de gravidade média, todos de enunciado:
  - o item 1 do Teorema mistura L_A com o problema de período completo;
  - há hipóteses escondidas;
  - H-rec (ida) e H-gauge (global) estão classificadas como física, mas são matemática não feita;
- 5 de gravidade baixa, de citação ou redação, entre elas a injetividade de K na terceira raiz, os lemas de deflação e
  equivalência e "com cone fixo seria 2";
- 3 de gravidade mínima.

Nenhum deles muda a contagem 1 sob as hipóteses declaradas.

## 1. Tabela de veredictos

| item | veredicto | artefato | o que foi feito |
|---|---|---|---|
| R1 L-real, toro | **OK com ressalva** | `art/r1_toro.py/.txt/.json`, `art/r1_BL.py/.txt/.json`, `art/l_real_toro_notebook.*` | Desigualdade derivada (§2.1). As normas exatas de ‖J⁻¹_m P_{n≥4096}‖, calculadas por recursão triangular O(n) + Lanczos (truncagem até 16.096), ficam abaixo de `cauda_Rinv` (razão 0,836–0,841). Os modos inteiros \|m\| = 400…1023 ficam abaixo da cota uniforme (razão 0,23–0,43). ‖B_L‖ truncado no toro com sinal: 3,6935 (12×36), igual ao "3,69" de F3, sobe a 3,850 (40×120), 3,859 (60×180) e 3,863 (80×240), abaixo de 3,9641. `l_real_toro.py` refeito inteiro no notebook (o PC 3 saiu do ar no meio): JSON idêntico ao gravado. Ressalva: o "+10⁻⁶" do erro do fundo é margem não derivada, irrelevante pela folga (§2.1). |
| R2 L-real, RT | **OK** | `art/r2_transfer.txt`, `art/r2_colunas.py/.txt/.json` | `transfer_bounds(1024, 4096, 31/10)` reproduz `build/t2/radius_transfer_s3_1.json`. Fórmulas recalculadas à mão em racionais. Somas de coluna pesadas exatas de Q contra as cotas de `STRUCTURED_EXTERIOR` (razão ≤ 1, justa em m = 0). Constantes de RT conferidas no TeX (§2.2). Extensão a \|s\| ≤ 3,1 trivial (folga até \|s\| ≈ 40 no lado forte e ≈ 7 no fraco). |
| R3 L-real, argumento | **OK com ressalva** | prova escrita em §2.3 | A prova fecha: inclusão, caixa, J⁻¹ preserva G, B limitado em Y+, unicidade no toro, cadeias e multiplicidades. Ressalva: a passagem "Rouché no toro → multiplicidade de L" usa dois lemas que não estão em L-real: equivalência por Q0(s) e a deflação de Brauer exata (§2.3, item 7). A "norma ≤ 1" da inclusão vale, mas pela razão certa (η ≥ 243), não pela citação de S3 §9.1. |
| R4 L-afim | **OK** | §2.4; `art/r7_symbolic.txt` (dC/ds) | Derivado do TeX de RT: a única derivada em τ do sistema selecionado está em O_μ = ∂τ + μD_O, com coeficiente 1; Γ1 e Γ2 não têm derivadas. Em C, a parte em s é só −Ξh2 na 2ª linha ♯. Conferido no código (`J_livre`, `bloco_B` sem s, `bloco_C`) e por `independent_rt_symbolic` (dC/ds). |
| R5 L-iso, L-constr, Lema 5.1 | **OK** | `art/r5_iso.py/.txt/.json` | Setor B montado com o MESMO montador e as paridades trocadas. Bloco a bloco, J_B(s) = J_A(s + i/2), B_B = B_A e C_B(s) = C_A(s + i/2): diferença 0 exata. Os blocos de B só dependem de (d, componente): diferença 0. ‖Mx‖/‖x‖ = 1,015625 = κ1. i-periodicidade com erro 1,1·10⁻¹⁶. Lema 5.1 conferido, e as retas Im s = −1/4, 1/4, 3/4 e Re s = 0, 3 são contorno certificado (§2.5). |
| R6 quociente | **OK** | `art/r6_quociente.py/.txt/.json` | 100 sorteios exatos (racionais) de A, B, D entrelaçador, C1 livre, com KC = ML. Conferidos coeficiente a coeficiente: M1 = C1, M0 e BD = DA; Dh_j = C(s0)h_j + C1h_{j−1}; dim física (definida à parte, pela observabilidade de D e^{−At} em E) = dim E − rank(D\|E); K(s0) injetivo ⇒ D\|E = 0. O atalho "ker C(s0) em cada vetor" erra a dimensão (25/25 sorteios no caso injetivo com Jordan). |
| R7 KC = ML | **OK com ressalva** | `art/r7_symbolic.txt`, `art/r6r7_unittest.txt`, `art/r7_numerico.py/.txt/.json` | `independent_rt_symbolic.py`: TODAS OK. `test_root_constraints`, `test_sharp_propagation` e `test_gauge_action`: 14 OK. M não existe na base com sinal, então KC − ML não pode ser montado. Testei a consequência: em s_K e na raiz B (que violam as constraints), C(s0)h é colinear ao vetor singular inferior de K(s0) (cos = 1,000000) e ‖KCh‖/(‖K‖‖Ch‖) = 3·10⁻⁷ (16×48), contra 0,12–0,16 com h aleatório. Em 0,7332 e μ, ‖Ch‖/‖h‖ → 0 com a malha (10⁻⁴ → 6·10⁻⁶). Ressalva: a extensão aos domínios (`SHARP_DOMAIN.md`) só foi lida, não refeita. |
| R8 cadeia de T2 | **OK com ressalva** | §2.6; `art/r8_contagem_discos.txt` | A lógica fecha. Lacunas de citação: (i) o passo 11 deve citar S3a ∪ disco de S3b ∪ R♯ = 1,765, porque a raiz 0,7332 não está localizada; (ii) a realidade das raízes 0,7332 e B usa a simetria de conjugação, que só vale em RT (Y+), não no toro com sinal. Conferido à parte que o Rouché da faixa B fecha só com os discos da faixa B (máx t = 516,5, t·ε_B = 0,031). |
| R9 hipóteses | **LACUNA (enunciado)** | §3 | Faltam no enunciado: simetria esférica no título; "modo" = modo de Floquet analítico (não estabilidade de semigrupo); a identificação do fundo de RT com a solução de Choptuik. "H-rec ida" e "H-gauge global" são matemática não feita, não modelagem. "Com cone fixo seria 2" só vale numa das duas leituras de "cone fixo". O item 1 do Teorema mistura L_A (largura 1) com o problema de período completo (ERRO de redação). |

## 2. Detalhe por item

### 2.1 R1: L-real, lado do toro

**Espaço.** O toro com sinal é o L² da superfície {|w| = κ1} × {|z| = κ2}, com w = e^{iτ/2} e ξ = (z + 1/z)/2,
restrito a funções simétricas em z ↦ 1/z. Os coeficientes são os de Laurent ("simétricos"), x_{m,n} = x_{m,−n}, e
‖x‖² = Σ_{m∈ℤ} Σ_{n∈ℤ} κ1^{2m} κ2^{2n} |x_mn|² = Σ_m Σ_{n≥0} κ1^{2m} ω2(n) |x_mn|², com ω2(0) = 1. É a norma de
`S3_CONDICIONAMENTO_SHARP.md` §9.1 (F3) e a de `nk_L.pesos`. Conferi que a matriz de J no código está nessa
convenção. Do TeX de RT, O_A = BandI_{1,−1}(Band_{0, im/2+μ(n+1)} + Band_{1, −im/2+μn}) dá entradas acima da diagonal
μ(c+1) + μ(c−1) = 2μc em TODAS as linhas j < c, sem a parte imaginária. Na convenção de Chebyshev (a_n = 2c_n) a linha
0 teria μc, então a convenção uniforme de `J_livre` é a de Laurent, coerente com ω2(0) = 1. A componente 1 (O_B, passo
2) dá o mesmo 2μc nas linhas ímpares.

**Desigualdade.** H(s) − I = (B + s)Q, com Q = J_μ(0)⁻¹ fixo. T = I − G é projeção ortogonal (coordenadas são
ortogonais na norma pesada), ‖T‖ = 1. Logo ‖T(H − I)T‖ ≤ ‖B + s‖ ‖QT‖ ≤ (‖B_L‖ + |s|) ‖J⁻¹T‖. Correto.

**‖J⁻¹T‖.** J é bloco-diagonal em m, então J⁻¹T = ⊕_m J_m⁻¹ P_m, com P_m = projeção nas colunas de T do modo m (todas,
se |m| ≥ M; n ≥ N, se |m| < M). Logo ‖J⁻¹T‖ = sup_m ‖J_m⁻¹P_m‖, e o peso κ1^m cancela no bloco. A cota de modo inteiro
para |m| ≥ 400 domina a das colunas, e por isso M = 1024 não entra no lado do toro. Rederivei as duas cotas.
- Entradas da forma fechada: |(J+σ)⁻¹_{jc}| ≤ (1/|d_j|)(v_c/|d_c|)·Π|1 − v_k/d_k|.
- Com Re σ ≥ 0: v_c/|d_c| < 2, Π ≤ 1 (|d|² − |d − v|² = 2μn(2μ + 2Re σ) ≥ 0) e 1/|d_j| ≤ 1/|Im σ| = 2/|m|.
- Razão de pesos: sqrt(ω2(j)/ω2(c)) ≤ c_w κ2^{j−c}, inclusive em j = 0.
- O teste de Schur dá (1 + 2c_w/(κ2 − 1))/(|m|/2). Em m0 = 400 isso é 0,052496.

**Conferência independente** (`art/r1_toro.py`). O operador exato W J_m⁻¹ P W⁻¹ foi aplicado por recursão triangular
O(n), sem formar matriz e sem overflow: a recursão usa só razões de pesos. O adjunto confere a 10⁻¹⁵, e a recursão
inverte `J_livre` do código a 6·10⁻¹⁶. σ_max por Lanczos (svds).

| caso | truncagem | norma exata | cota | razão |
|---|---|---|---|---|
| colunas n ≥ 4096, m = 0, 1, 2, −7 (passo 1) | 7.096 e 16.096 | 0,012803 | 0,015317 (`cauda_Rinv`) | 0,836 |
| idem, m = 0 e 2, passo 2 (componente 1, n ímpar) | idem | 0,006486 | 0,007709 | 0,841 |
| idem, m = 40 / −150 (p. 2) / 399 | idem | 0,012793 / 0,006415 / 0,011849 | 0,015311 / 0,007669 / 0,014749 | 0,836 / 0,836 / 0,803 |
| modo inteiro m = ±400 (passo 1 / 2) | 3.000, 6.000, 12.000 | 0,0226930 / 0,0119160 | 0,052491 | 0,43 / 0,23 |
| modo inteiro m = 401, 700, 1023 (passo 1) | 6.000 (1023: até 12.000) | 0,02264 / 0,01299 / 0,0088933 | 0,05236 / 0,02999 / 0,02052 | 0,43 |

A norma truncada não muda de 7.096 para 16.096, nem de 3.000 a 12.000 nos modos inteiros (`art/r1_toro_convergencia.txt`): convergiu. Ela é cota inferior da infinita, porque o bloco líder de um
operador triangular superior é exato. Logo as cotas de `l_real_toro.py` são válidas e frouxas por um fator ~1,2 a 4.

**‖B_L‖ no mesmo espaço** (`art/r1_BL.py`). Montei W B W⁻¹ com os blocos do código (`signed_operator_L.bloco_B`) e
pesos κ1^m sqrt(ω2(n)), e calculei σ_max por Lanczos:

| truncagem | dim | ‖B_L‖ truncado |
|---|---|---|
| 12×36 | 1.422 | 3,6935 (F3, "verdadeiro 12×36": 3,69) |
| 24×72 | 5.868 | 3,8208 |
| 40×120 | 16.500 | 3,8501 |
| 60×180 | 37.350 | 3,8593 |
| 80×240 | 66.600 | 3,8626 |

O valor de 12×36 reproduz o de F3: F3 e o L-real usam a mesma norma. A sequência sobe devagar e fica abaixo de 3,9641.
(Para multiplicadores, a norma em L² é o sup do símbolo, e F3 cota esse sup em Arb.)

**Refazer `l_real_toro.py`.**
- PC 3 (`~/revisao_T2/out/`, scripts conferidos por sha256): só o primeiro μ terminou (máx cauda 0,015317 em
  (0, 1, μ−)), porque o PC saiu do ar no meio do segundo.
- Notebook (`art/l_real_toro_notebook.{txt,json}`, 2975 s): o script inteiro deu cauda 0,015317, uniforme 0,052496 e
  defeito 0,37084, CONTRAI. O JSON é idêntico, campo a campo, a `build/t2/l_real_toro.json`.

**Ressalvas.**
- (a) O "+10⁻⁶" do erro do fundo em ‖B_L‖ não está derivado no script. Como δB = 2SΞΓ2(ω* − RefA, ·) é multiplicação
  por campos com ‖·‖_{ℓ¹} ≤ ε = 2⁻²⁵ + 2⁻²⁷⁷ e coeficientes de tabela O(1), a ordem é 10⁻⁷. Mas a contração só exige
  ‖B_L‖ < 1/0,0525 − 3,1 ≈ 15,9, então isso não pesa.
- (b) A cauda é calculada só em μ_RefA ± 10⁻¹⁰, não no intervalo. A cota dominante (modo inteiro) não depende de μ, e a
  cauda (0,0153) está 3,4× abaixo dela, então isso também não pesa.
- (c) O argumento `--M` de `l_real_toro.py` não é usado (correto, mas engana quem lê).

### 2.2 R2: L-real, lado de RT

- `transfer_bounds(1024, 4096, 31/10)` reproduz exatamente `radius_transfer_s3_1.json` (`art/r2_transfer.txt`).
- Recalculado em racionais:
  - t_m = (3/2)/(b0(1 − r)) = 0,014648 (b0 = 512, r = 4/5);
  - t_n = 27/(2μ_min(N+1)) = 0,019771;
  - o `input_tail` antigo dá 0,035156, e o mínimo fica com o novo 0,019771;
  - defeito forte = (5191/500 + 31/10)·0,019771 = 0,26655;
  - lado fraco: cauda 0,037344, defeito 0,85470.
- Lado fraco (Y− = pesos (129/128, 9/8), usado em H-reg "entre os dois raios"): β− ≤ (162/85)β+.
  - A divisão regularizada por ξ vale 2κ2⁻¹/(1 − κ2⁻²): 40/9 em κ2 = 5/4 e 144/17 em 9/8.
  - As multiplicações não crescem ao diminuir o peso.
  - Cada termo de β+ tem no máximo um fator de divisão, então a escala vale termo a termo.
  - Com |s| ≤ 3,1, o defeito fraco 0,8547 < 1. Isso estende `RADIUS_RECOVERY.md` (escrito para |s| ≤ 5/4) à região
    toda.
- A contração vale até |s| ≈ 40 (forte) e |s| ≈ 7,0 (fraco). A extensão de 5/4 para 3,1 é legítima, porque Q = J_μ(0)⁻¹
  não depende de s e s só entra no fator (β + |s|).
- **Derivação conferida.**
  - Escrevi a coluna c da forma fechada como q_{c−k,c}·Π f_j, com f_j = (σ + μ(1−j))/(σ + μj). Vale
    |σ + μj|² − |σ + μ(1−j)|² = μ(2j−1)(μ + 2Re σ) ≥ 0, então |f_j| ≤ 1.
  - Daí a soma de coluna com pesos de RT (razão r^{c−j}; em j = 0, o peso 1 contra 2κ2^c só ajuda) é
    ≤ 1/|D_c| + |q_{c−k,c}| r^k/(1 − r^k), com as cotas 9/(μ(n+1)) (k = 1) e (1 + 32n/(9(n−1)))/(μ(n+1)) (k = 2).
  - A parte de Fourier usa 2b0x ≤ b0² + x².
  - O fator 3/2 dos blocos reais 2×2 é √2 ≤ 3/2.
- **Conferência numérica** (`art/r2_colunas.py`). Somas de coluna exatas nos modos m = 0, 1, 2, 7, 1023 e colunas
  4096…4495: razão máxima soma/cota = 0,9992 (k = 1) e 1,0000000000000 (k = 2, m = 0).
  - A igualdade é real: em σ = 0 o produto telescopa e todas as entradas fora da diagonal ficam iguais a |q_{c−2,c}|.
  - A cota é justa só no limite (falta a cauda geométrica r^n). O excesso de 3,8·10⁻¹⁴ é arredondamento.
  - Modos m = 1024, 1025, 3000: soma máxima 0,00922 contra a cota 0,00977.
- **Constantes que vêm de RT** (conferidas no TeX `choptuik.tex`):

  | constante | onde no TeX | leitura |
  |---|---|---|
  | κ = (65/64, 5/4) | l. 320 | ok |
  | norma SSPACE = Σ_i ‖v_i‖, ‖v‖ = Σ_{m,n≥0}(2−δ_{m0})(2−δ_{n0})κ1^m κ2^n(\|Re\|+\|Im\|) | l. 1543–1550 | ok: sem pesos η; os η = (916, 4096, 243, 1233) são do projeto, uma norma equivalente |
  | ‖RefB‖ ≤ 𝓛2 = 2⁻²⁵ | (dkjfhk2), l. 1722, 2184 | ok |
  | ‖ω* − (RefA + RefB)‖ ≤ 2⁻²⁷⁷ e \|μ* − (μ_A + μ_B)\| ≤ 2⁻²⁷⁷ | l. 2355–2368 | ok, logo ε = 2⁻²⁵ + 2⁻²⁷⁷ |
  | ‖SΞΓ2(v,w)‖ ≤ 3‖v‖‖w‖ | l. 1645 | ok: dá o termo 6ε·(max η/min η) |
  | 𝓛1 = 10,5 | l. 1721, 2183 | recalculei 10,136 ≤ 10,5 a partir de RefA |
  | μ* = 0,1683070789634… | l. 204 | ok, dentro de [1/6, 17/100] |

- β+ = 5191/500 não é constante publicada. É calculada pelo projeto a partir de RefA (dados de RT, sha fixado) mais o
  termo do erro. Recalculei com as funções de `independent_component_beta.py`, sem ler `build/spectrum`
  (`art/r2_beta.txt`): β_η = 10,38108 ≤ 10,382. Sem os pesos η a cota seria 21,78. Mesmo assim o lado forte contrairia
  (0,49), então a escolha de η não é crítica.

### 2.3 R3: L-real, o argumento (prova completa)

Notação: Y+ = complexificação dos DOFs reais de RT, com ‖x‖_+ = Σ_c η_c ‖x_c‖_{SSPACE-componente}; X = toro com sinal;
D+ = Q(Y+), D_X = Q(X); H(s) = L(s)Q = I + (B + s)Q; G = {|m| < 1024, n < 4096}, T = I − G.

1. **Mesma coordenada, mesmo operador.** A troca (a, b) ↦ (x_m, x_{−m}) = (a + ib, a − ib), com m = 0 só a, é uma
   bijeção C-linear entre DOFs reais complexificados e coeficientes de Laurent. É a matriz T de
   `signed_operator_L.validar`, e com ela L_sinal = T L_RT T⁻¹ entrada a entrada (6×18, 8×24). Os x_m são os
   coeficientes de Laurent da função: RT define ‖v‖ = Σ_{m∈ℤ}…, com v_{−m} = conj(v_m) para funções reais. Q, B e G
   são dados pela mesma matriz infinita nos dois espaços. Q é diagonal em m e B é convolução (Toeplitz em m, R5).
2. **Inclusão Y+ ⊂ X**, com constante explícita. sqrt(ω2(n)) ≤ (2 − δ_{n0})κ2^n. Para o par ±m,
   κ1^m|a + ib| + κ1^{−m}|a − ib| ≤ (2 − δ_{m0})κ1^m(|a| + |b|). Com ℓ² ≤ ℓ¹ por componente e
   Σ_c ‖x_c‖ ≤ ‖x‖_+/min η, sai ‖x‖_X ≤ ‖x‖_+/243. A afirmação "norma ≤ 1" é verdadeira. Mas a citação de
   `S3_CONDICIONAMENTO_SHARP.md` §9.1 não a prova: lá a inclusão é do ℓ¹ do próprio toro, não do SSPACE de RT com η.
   Também D+ ⊂ D_X, porque Q(Y+) ⊂ Q(X).
3. **Q preserva a caixa.** J é diagonal em m e triangular superior em n, então J⁻¹ também é: QG = GQG, logo TQG = 0. Daí
   THG = T(B + s)QG = TBQG.
4. **B limitado em Y+** (β_η ≤ 5191/500, conferido em §2.2) e em X (F3). Para x finito, QGx é finito, e BQGx ∈ Y+.
5. **Núcleo.** Seja x ∈ X com H(s)x = 0 e |s| ≤ 3,1. Aplicando T: THT(Tx) = −THGx =: f ∈ Y+.
   - Em Y+, ‖T(H − I)T‖_+ ≤ 0,2666 < 1. Então THT|_{Ran_+T} é invertível (Neumann), e existe y ∈ Ran_+T ⊂ Y+ com THTy = f.
   - THTy = f vale também em X. Os operadores de Y+ e de X coincidem nos vetores finitos (mesma matriz). Os vetores
     finitos são densos em Y+, e Y+ ↪ X é contínua. Então, para y ∈ Y+, A_X y = lim A y⁽ⁿ⁾ = A_+ y em X, com
     y⁽ⁿ⁾ = truncagens.
   - Em X, ‖T(H − I)T‖_X ≤ 0,371 < 1, então a solução em Ran_X T é única. Logo Tx = y, x = Gx + y ∈ Y+ e h = Qx ∈ D+.
   - A recíproca, ker_{Y+} ⊂ ker_X, é a inclusão.
6. **Cadeias e multiplicidade.** Como L(s) = L(0) + s (L-afim), uma cadeia é (L(0) + s0)h_j = −h_{j−1}. Com x_j = Q⁻¹h_j:
   THT(Tx_j) = −Th_{j−1} − THGx_j. Por indução, h_{j−1} ∈ D+ ⊂ Y+, e o passo 5 se repete. Logo os espaços de raízes
   E = ∪_k ker(L(0) + s0)^k coincidem como conjuntos. Todo vetor gerado é o último elo de uma cadeia que começa num
   autovetor. Como dim E é a multiplicidade algébrica (pencil linear com coeficiente identidade), as multiplicidades
   coincidem. Pela compacidade de Q nos dois espaços (as caudas → 0), H é Fredholm de índice 0 e E tem dimensão finita
   (F1).
7. **O que L-real não diz, e a contagem usa.** O Rouché do toro conta os zeros de H̃(s) = L̃(s)Q0(s), com
   Q0(s) = J(s)⁻¹ e L̃ = L deflacionado (fase e μ, Brauer). Para que isso seja a multiplicidade de L em X, são precisos
   dois lemas.
   - (a) **Equivalência.** Q0(s) é família holomorfa de isomorfismos X → D_X na norma de gráfico, com inversa holomorfa
     J(s) = J(0) + s. Ela é invertível em Re s > −μ, onde d_k = μ(k+1) + s + im/2 ≠ 0. Logo
     m(H̃; s0) = m(L̃; s0) (Gohberg–Sigal). Vale.
   - (b) **Brauer.** Com p* = ∂τω* e g* = (Z* + 1)ω* exatos (L(0)p* = 0 e L(μ*)g* = 0: identidades de covariância, que
     `independent_rt_symbolic` confere) e φ biortogonais, span(p*, g*) é invariante sob Ã, e Ã induz no quociente o
     mesmo operador que A. Então N(L; Ω_A) = N(L̃; Ω_A) + 1 (μ*), e s = 0 sai do contorno (vai para s = −1).

   (a) é padrão. (b) pertence à contagem já revisada, que trata a diferença entre os vetores de RefA e os exatos como
   "deflação verdadeira". Ambos deveriam estar citados no L-real ou no passo 6 de T2, porque o L-real sozinho transporta
   multiplicidades de L, não de L̃.
8. **Consequências do L-real.**
   - Os autovetores certificados no toro estão em Y+: certo.
   - "A injetividade de K no toro implica a injetividade em RT": certo, com a ressalva de §2.6, item 3.
   - Exclusão em Re s ≥ 2,93: basta a inclusão (núcleo em Y+ ⇒ núcleo em X), sem o L-real. Certo.

**Veredito R3.** O argumento está correto. As ressalvas são de completude da citação: itens 2 e 7.

### 2.4 R4: L-afim

- TeX de RT (`choptuik.tex`, definições de O_A, O_B, H_μ, Γ1, Γ2 e Ω⁺, l. 1093–1320):
  - Ω⁺(λ, v) = ΞH_λ v + λΓ1v + ΞΓ2(v, v);
  - O_μ = SΞH_μ = ∂τ + μD_O, com D_O = diag(D_B, D_A, D_A, D_A) e D_A = (1 + ξ)∂ξ + 1;
  - Γ1 e Γ2 são algébricos (convoluções, P, sem derivadas).
- Com μ fixo, a perturbação e^{sτ}h leva ∂τ → ∂τ + s só onde ∂τ age em h. O resultado é
  L(s)h = (∂τ + s + μD_O)h + μSΓ1h + 2SΞΓ2(ω*, h), ou seja, L(s) = L(0) + s·I exatamente. A fórmula do documento
  "(s + μD_O + B_*)h" omite o ∂τ: erro tipográfico.
- C(s) = S♯{ΞH_{μ,s}h + …}. Na 1ª linha ♯, ½(1 − P) mata Ξ∂τh1, que é par porque h1 é ímpar. Na 2ª, −P(−ΞO_A P h2)
  dá C1h = −Ξh2. Logo C(s) = C0 + sC1.
- O código bate:
  - `J_livre` usa σ = s + im/2;
  - `bloco_B` não depende de s;
  - em `bloco_C`, s só entra em −σΞ (componente 1);
  - `independent_rt_symbolic`: "dC/ds = (0, −ξh2)" OK.
- K(s) = B_K + s: `signed_operator.Jrad` usa σ = s + im/2, e a parte de fundo não depende de s.

### 2.5 R5: L-iso, L-constr, Lema 5.1

- **Conjugação.**
  - M = multiplicação por e^{iτ/2} satisfaz M⁻¹∂τM = ∂τ + i/2 e comuta com multiplicações por funções de τ, com P e
    com as operações em ξ.
  - Com a equivariância T∘L = L∘T (o que define os setores) e L-afim: M⁻¹L_B(s)M = L_A(0) + (s + i/2) = L_A(s + i/2).
  - O mesmo vale para C (a parte em s de C é a de ∂τ).
  - M leva D_A em D_B porque J_B(s)M = MJ_A(s + i/2).
- **No código** (`art/r5_iso.py`): construí o setor B com o próprio montador, trocando as paridades (componentes 0, 1 e 2
  nos m ímpares, componente 3 nos pares). Resultados:
  - max|J_B(s)(m+1) − J_A(s + i/2)(m)| = 0;
  - max|B_B(m+1+d, m+1) − B_A(m+d, m)| = 0;
  - max|C_B(s)(m+1 ← m_i+1) − C_A(s + i/2)(m ← m_i)| = 0;
  - os blocos de B para m de mesma paridade e mesmo d são idênticos (0 exato);
  - ‖Mx‖/‖x‖ = 1,015625 = κ1;
  - J_A(s + i)(m) = J_A(s)(m + 2) a 1,1·10⁻¹⁶.
- **Lema 5.1.** Q = {−1/4 < Im s ≤ 1/4} e Q + i/2 = {1/4 < Im s ≤ 3/4} são disjuntos e de união Q̃. Σ_B ∩ Q corresponde
  bijetivamente (s ↦ s + i/2), com multiplicidades, a Σ_A ∩ (Q + i/2). A prova está certa.
- As retas de fronteira não importam. Os contornos certificados (cobertura: faixa A, arestas Re s = 0, 3 e
  Im s = ±1/4; faixa B, Re s = 0, 3 e Im s = 1/4, 3/4) mostram que L̃ é invertível sobre elas. As duas faixas usam a
  mesma matriz deflacionada (`diagT_lab2.npy`). Na faixa B a deflação é inócua: os autovalores deflacionados ficam em
  s = −1 e as raízes 0 e μ não estão na faixa. Por
  i-periodicidade, Im s = −1/4 ≡ 3/4. Então "semiaberto" contra "fechado" não muda a contagem.
- `S4_GAUGE.md` usa a convenção {−1/4 ≤ Im s < 3/4} e `T2_HIPOTESES.md` usa {−1/4 < Im s ≤ 3/4}. A diferença é
  inofensiva pelo mesmo motivo.

### 2.5b R6: o quociente das constraints

**Derivação.** Com L = A + s, K = B + s, C = C0 + sC1 e M = M0 + sM1, igualando potências de s em
(B + s)(C0 + sC1) = (M0 + sM1)(A + s):
- s²: M1 = C1;
- s¹: C0 + BC1 = M0 + M1A, então M0 = C0 + BC1 − C1A;
- s⁰: BC0 = M0A.

Com D = C0 − C1A: BD = BC0 − BC1A = M0A − BC1A = (C0 − C1A)A = DA. Confere com o documento.

Reciprocamente, toda solução vem de um entrelaçador D e de C1 livre. É isso que o teste sorteia.

**Cadeias.** De (A + s0)h_j = −h_{j−1}: Dh_j = C0h_j − C1Ah_j = C(s0)h_j + C1h_{j−1}. Confere.

**Por que ker(D|E) é o espaço físico.**
- Uma solução generalizada de Floquet é u(τ) = e^{−Aτ}u0 com u0 ∈ E, isto é, e^{s0τ}·(polinômio em τ). Ela satisfaz as
  constraints em todo τ se e só se C(∂τ)u = C0u − C1Au = De^{−Aτ}u0 ≡ 0.
- Isso equivale a u0 pertencer ao maior subespaço A-invariante contido em ker D|E. Esse subespaço é o próprio ker(D|E),
  que já é A-invariante (BD = DA).
- Daí dim = dim E − rank(D|E). O teste confere essa dimensão pela matriz de observabilidade [D; DA; DA²; …] restrita
  a E, um cálculo independente da fórmula.

**Gauge.** g ∈ ker(D|E), com Ag = −μg, gera uma linha invariante. O quociente tem dimensão dim ker(D|E) − 1, sem
exigir complemento invariante. O exemplo da cadeia mista (A = [[−μ, 1], [0, −μ]], D = 0) confere:
`test_root_constraints` OK.

**Teste** (`art/r6_quociente.py`, racionais exatos, 4 configurações × 25 sorteios):
- KC = ML, a fórmula das cadeias, dim física = dim E − rank(D|E) e "K(s0) injetivo ⇒ D|E = 0": todos OK.
- O atalho "exigir C(s0)h = 0 em cada vetor de E" dá a dimensão errada em 25 de 25 sorteios do caso com K injetivo e
  bloco de Jordan 2, em 2 de 25 e em 1 de 25 nos outros. Isso confirma o alerta do documento.

### 2.6 R8: da contagem a "1 modo físico"

1. **μ\*.** A álgebra (R6) dá contribuição dim E − rank(D|E) − [gauge]. Com E = span(g*) (simplicidade, S4) e
   Dg* = C(μ*)g* = 0 (identidade de gauge no fundo exato), E_constr = E. A subtração da linha de gauge exige H-cone:
   cone móvel, g* admitido e quocientado.
   - A simplicidade é necessária. Sem ela, uma cadeia sobre g* sobreviveria no quociente (o exemplo do próprio
     documento; reproduzido em R6).
   - Dg* = 0 é automático pela identidade.
   - H-gauge (completude) NÃO é necessária para μ*. É necessária para 0,7332, para garantir que essa raiz não seja
     gauge.
2. **s_K e raiz B.** Para raiz simples, E = span(h), e Dh = C(s0)h (h_{−1} = 0). O testemunho ‖Π C h‖ > 0 dá
   rank(D|E) = 1 = dim E, então a contribuição é 0. Basta o testemunho mais a simplicidade (NK).
   - Consistência com KC = ML: se C(s0)h ≠ 0, então K(s0) é singular no mesmo s0. Em 16×48, em s_K e em
     0,041419 + 0,5i, σ_min(K)/‖K‖ = 3·10⁻¹² e 3·10⁻¹¹, e C(s0)h é colinear ao núcleo de K (R7). Isso fecha.
3. **0,7332.** BD = DA em E, e (B + s0)^k Dh = D(A + s0)^k h = 0. Com K(s0) injetivo, Dh = 0 em todo E (R6). A raiz é
   simples (a contagem 2 de L̃ menos s_K), então "todo o espaço de raízes" é só o autovetor. A injetividade no ponto
   certo exige atenção.
   - A raiz não está localizada rigorosamente. O passo 11 ("K injetivo em 0,7332 (região de S3a)") pressupõe a
     localização.
   - O argumento correto: K é injetivo em toda a faixa A fechada, Re s ≥ 0, menos {s_K}. A cobertura é S3a
     ([0; 1,765] × [−1/4, 1/4] fora do disco de S3b), mais o Rouché de K no disco (zero único s_K), mais R♯ = 1,765
     (sem espectro em Re s > 1,765, qualquer Im s). A terceira raiz é ≠ s_K, porque s_K é simples para L. Logo K é
     injetivo nela, onde quer que esteja.
   - Espaços: Dh ∈ Y♯(129/128, 9/8) (perda de raio, `SHARP_DOMAIN.md`). Y♯ ⊂ toro de K, de mesmos raios, pelo mesmo
     argumento de §2.3, item 2. A injetividade no toro de K passa para o subespaço. Atenção: o L² com sinal de raio
     65/64 NÃO está contido no de raio 129/128, por causa dos m → −∞. A passagem precisa do L-real (h ∈ Y+) antes de
     aplicar C. A cadeia está certa, mas a ordem importa.
4. **Realidade.** "O modo é real" (T2_HIPOTESES §1) e "raiz B ↔ modo real do setor B" usam a simetria de conjugação
   de Σ_A. Ela vale em RT (operador real nos DOFs reais, Y+ simétrico em m), não no toro com sinal (a conjugação leva
   |w| = κ1 em |w| = 1/κ1).
   - Com o L-real, Σ_A na região é o mesmo nas duas realizações e é simétrico.
   - Na faixa A (simétrica em Im s = 0), sobra uma raiz depois de μ* e s_K (ambas reais): ela é real.
   - Na faixa B (simétrica em Im s = 1/2), a raiz única fica sobre Im s = 1/2.
   - O argumento não está escrito. Também não é necessário para a contagem 1, porque S3a cobre |Im s| ≤ 1/4 inteiro.
5. **Periodicidade e fronteiras.** Conferidas em §2.5. Re s = 0 é aresta certificada.
   - Em s = 0 a fase foi deflacionada (para s = −1), e L̃(0) é invertível no contorno. Isso prova que a raiz de L em 0
     tem multiplicidade algébrica exatamente 1, e que não há outras raízes em Re s = 0 nas duas faixas.
   - A contagem é de fato em Re s > 0 aberto, e a raiz de fase fica de fora corretamente.
   - Os aliases de 0 e de μ (Im inteiro) não caem na faixa B.
   - A cobertura da faixa B foi refeita só com os discos de B (`art/r8_contagem_discos.txt`): máx t = 516,5,
     t·ε_B = 0,031, 1 zero. O log original cita um disco da faixa A como o pior, provavelmente porque foi gerado com
     A ∪ B (conservador). É questão de rastreabilidade, não de erro.
6. **A conta final.** 0 (μ*) + 0 (s_K) + 0 (raiz B) + 1 (0,7332) = 1. Por H-rec (volta: `GEOMETRIC_FLOQUET_RECONSTRUCTION`,
   válida em Re s > 0, inclusive cadeias) e H-gauge (0,7332 ∉ {μ(1 − n) + ik}), esse 1 é um modo físico. Por H-rec
   (ida), não há outros. A lógica está correta.

### 2.7 R7: KC = ML

- **Simbólico.** `independent_rt_symbolic.py` (transcrição independente do TeX) verifica
  K(s)C(s)h − M(s)L(s)h + G(h, e(w)) = 0. É uma identidade polinomial nos jatos de ordem ≤ 2, nos pontos (τ0, ±ξ0),
  para todo fundo w, todo s e todo μ.
  - Em cada grau homogêneo nos campos, a base τ^aξ^b usada realiza qualquer jato da ordem que aparece. As contagens
    batem: 12 condições para ordem 2 e grau 5 − 2a em ξ.
  - Isso dá a identidade pontual para campos C², portanto para os analíticos dos espaços usados.
  - Também conferidos: tabelas dos macros = TeX (0 divergências), dC/ds, covariâncias de gauge e de fase, e a
    reconstrução. Saída: TODAS OK.
- **Numérico** (`art/r7_numerico.py`, 12×36 e 16×48). Em s_K e na raiz B:
  - ‖KCh‖/(‖K‖‖Ch‖) cai de 6·10⁻⁶ para 3·10⁻⁷ com a malha (h aleatório: 0,12–0,16);
  - cos(Ch, vetor singular inferior de K) = 1,000000.

  Em 0,7332 e em μ, ‖Ch‖/‖h‖ cai de 1,1·10⁻⁴ / 1,5·10⁻⁴ para 5,7·10⁻⁶ / 8,9·10⁻⁶.
- **Por que não testei KC − ML diretamente.** M não está implementado na base com sinal. Além disso, num corte finito,
  KC − ML tem defeito de borda (derivadas acoplam n a n ± 1, e C, K saem do corte) e o defeito G(h, e(RefA)). Testar a
  consequência em ker L evita os dois problemas.

### 2.8 R9: as hipóteses

**Precisão do enunciado.**
- O item 1 do Teorema ("Contagem") não depende de nenhuma hipótese física. Depende só de E1 e dos certificados. As
  hipóteses H-* só entram no item 2. Vale separar em Teorema A (incondicional, módulo E1) e Teorema B (sob H-*).
- O item 1 fala de "L(s) … no setor de período completo" e de "classe de Floquet (faixa semiaberta de largura 1)", e
  lista "≈ 0,0414 + 0,5i (setor B)". Para o operador de período completo, o reticulado de Floquet é (i/2)ℤ. A classe é
  uma faixa de largura 1/2, e nela as raízes são μ*, s_K e 0,7332 (setor A) e s_B ≈ 0,0414, real (setor B). Numa faixa
  de largura 1, o operador de período completo tem 8 raízes. O que tem 4 raízes numa faixa de largura 1 é L_A, e
  "0,0414 + 0,5i" é a imagem da raiz B em L_A. É erro de redação: a contagem física está certa pelo Lema 5.1.
- "Gauge admissível" (H-rec, ida) não está definido. Precisa ser compatível com H-cone, isto é, com f(0) livre.

**Hipóteses escondidas** (usadas e não enunciadas em `T2_ENUNCIADO.md`):
1. **Simetria esférica.** O título de `T2_HIPOTESES.md` ("tem exatamente um modo instável físico") e o Teorema não
   dizem "esfericamente simétrico". Só H-rec diz. Perturbações não esféricas não são tratadas.
2. **Noção de "modo".** É um modo de Floquet e^{sτ}h, com h na classe analítica de H-reg, contado no espectro do pencil.
   O teorema é espectral. Não afirma nada sobre o semigrupo linear, completude dos modos, estabilidade linear de dados
   gerais nem a variedade estável não linear (`SPECTRAL_PROBLEM.md` §5 diz isso, T2 não).
3. **Identificação do fundo.** ω* de RT é "a" solução crítica de Choptuik, o atrator do limiar. RT provam a existência
   de uma solução DSS perto de RefA; a identificação com o atrator do colapso é física e numérica.
4. **Coordenadas e μ fixos.** A perturbação é descrita nas coordenadas de RT, com μ fixo (μ não é perturbado). Isso é
   coerente com H-rec, mas fica implícito.
5. **Direção do tempo.** "Instável" é Re s > 0 porque τ → +∞ é a aproximação da singularidade
   (u_σ = −(1 + σξ)e^{−μτ} → 0). Correto, mas implícito.

**Matemática não feita × modelagem.**

| item | natureza |
|---|---|
| H-rec, ida (pôr toda perturbação linearizada no ansatz de RT, com centro, cone e Floquet, por um gauge) | **matemática não feita** (teorema de fixação de gauge para EDP); depende da modelagem só pela escolha de classe e de cone |
| H-gauge, versão global e germes não analíticos | **matemática não feita** (classificação global dos difeomorfismos que preservam o ansatz) |
| H-rec, volta | feita no ansatz para Re s > 0 (`GEOMETRIC_FLOQUET_RECONSTRUCTION.md`; `test_geometric_reconstruction`: 5 OK) |
| H-cone (quocientar a mudança de T*) | **modelagem/convenção** |
| H-reg (classe analítica nos raios de RT) | **modelagem**. A independência entre Y+, Y− e o toro está provada (L-real, `RADIUS_RECOVERY`). Estender a C^∞ ou Sobolev é matemática não feita e, em geral, falsa no sentido ingênuo: as soluções homogêneas singulares em ξ = −1 ou 0 (`SHARP_DOMAIN.md`, passo 4) entram em classes de regularidade finita |
| simetria esférica; fundo = solução de Choptuik | **modelagem** |
| localização de 0,7332 (NK) | matemática não feita, mas desnecessária |

**"Com cone fixo seria 2".**
- Só é verdade na leitura "mesma realização espectral de RT, sem quocientar a translação nula": aí g* é um modo com
  constraints satisfeitas e conta. A frase de `T2_ENUNCIADO` diz isso.
- Numa realização genuína de cone fixo, `S4_GAUGE.md` diz que "g é excluído", e então μ* não seria raiz. Também
  `GEOMETRIC_FLOQUET_RECONSTRUCTION.md` mostra que preservar o cone impõe k₊₊|_{u₋=0} = 0, uma obstrução à fixação de
  gauge. É outro problema espectral, que o projeto não contou.
- Sugiro: "sem quocientar a translação nula das coordenadas, μ* contaria e a contagem seria 2".
- Nota relacionada: `S4_GAUGE.md` chama μ* de "modo de gauge da literatura (λ = 1 no tempo logarítmico)". Pela
  conversão do próprio projeto, λ = 2πs/K, s = μ* dá λ ≈ 0,614, não 1. Não há contradição: o expoente de um modo de
  gauge depende do gauge. Em RT o vértice se move por translação de u, e u não é tempo próprio
  (a menos de fator periódico, T* − t ∝ |u|^{K/(2πμ)} no centro, K/(2πμ) ≈ 1,63, porque Θ: τ ↦ τ + 2π reescala u por e^{−2πμ} e o tempo próprio por e^{−K}). Mas o texto deveria dizer isso em vez de "λ = 1".

## 3. Lacunas, gravidade e conserto

Gravidade: **alta** = pode mudar a conclusão; **média** = o enunciado não diz o que se prova; **baixa** = citação ou
redação; **mínima** = cosmética ou rastreabilidade. Nenhuma lacuna alta foi encontrada.

| # | onde | lacuna | gravidade | conserto proposto |
|---|---|---|---|---|
| 1 | `T2_ENUNCIADO` Teorema, item 1 | Mistura o operador de período completo com L_A: "faixa de largura 1" e "0,0414 + 0,5i (setor B)". Para período completo, a classe tem largura 1/2, e a raiz B é s_B ≈ 0,0414 real | média (erro de enunciado; a conta não muda) | Enunciar para L_A em Q̃ = {−1/4 < Im s ≤ 3/4} e traduzir pelo Lema 5.1: "no problema de período completo, por classe (largura 1/2): μ*, s_K, 0,7332 (A) e 0,0414 (B)" |
| 2 | `T2_ENUNCIADO`, `T2_HIPOTESES` §1 | Hipóteses escondidas: simetria esférica fora do título; "modo" = Floquet analítico, sem afirmação de estabilidade linear ou semigrupo; fundo de RT = solução de Choptuik; μ e coordenadas fixos | média | Pôr no Teorema: "perturbações esfericamente simétricas, modos de Floquet na classe H-reg" e um parágrafo "o que T2 não afirma" |
| 3 | `T2_ENUNCIADO` §hipóteses | H-rec (ida) e H-gauge (global) aparecem como "físicas", mas são teoremas de EDP não provados | média (classificação) | Separar "conjecturas matemáticas abertas" (H-rec ida, H-gauge global) de "convenções/modelagem" (H-cone, H-reg, simetria esférica) |
| 4 | passo 11 de `T2_ENUNCIADO`; `T2_HIPOTESES` §1, item 8 | "K injetivo em 0,7332 (região de S3a)" pressupõe a localização, que é só numérica | baixa | Citar: K injetivo em {Re s ≥ 0, \|Im s\| ≤ 1/4} ∖ {s_K} por S3a ∪ (Rouché de K no disco de S3b) ∪ (R♯ = 1,765); a terceira raiz é ≠ s_K porque s_K é simples |
| 5 | L-real ("Consequências") e passo 6 | A passagem do Rouché do toro (zeros de L̃Q0(s)) às multiplicidades de L usa a equivalência por Q0(s) e a deflação de Brauer exata (N(L) = N(L̃) + 1), não citadas no L-real | baixa | Acrescentar os dois lemas (§2.3, item 7) ou remeter ao lugar da contagem onde estão |
| 6 | L-real, "Inclusão" | "Norma ≤ 1" citada de S3 §9.1, que trata de outro ℓ¹. Vale por ‖x‖_toro ≤ ‖x‖_+/243 | baixa | Substituir pela conta de §2.3, item 2 |
| 7 | `T2_HIPOTESES` §1 ("o modo é real") | A realidade de 0,7332 e da raiz B usa a simetria de conjugação, que vale em Y+ e não no toro com sinal | baixa | Escrever: via L-real, Σ_A na região é o de RT, simétrico; uma raiz restante numa região simétrica é real (A) ou fica em Im s = 1/2 (B) |
| 8 | `S4_GAUGE` e "com cone fixo seria 2" | Ambíguo: vale só como "sem quocientar g* na mesma realização". Uma realização de cone fixo exclui g*. E "λ = 1" contradiz λ = 2πs/K (dá 0,614) | baixa | Reformular (§2.8) e explicar que o expoente de gauge depende do gauge |
| 9 | `scripts/l_real_toro.py` | `--M` não usado; "+10⁻⁶" do erro do fundo não derivado; cauda só em μ_RefA ± 10⁻¹⁰ | mínima (folga: ‖B‖ até 15,9 contrai) | Comentar `--M`; derivar δB ≤ 2·3·c·ε no toro, ou citar a folga |
| 10 | L-afim | A fórmula "(s + μD_O + B_*)h" omite ∂τ | mínima | Escrever (∂τ + s + μD_O + B_*)h |
| 11 | `build/faixaB/contagem_rouche.txt` | O log cita um disco da faixa A como o pior e não lista os discos usados | mínima (refeito só com os discos da faixa B: OK, `art/r8_contagem_discos.txt`) | Gravar a lista e o sha256 dos discos no log |

## 4. O que não verifiquei

- **Tomados como dados**, conforme a TAREFA:
  - N_A = 3 e N_B = 1 (discos, janelas, Schur), inclusive a "deflação verdadeira" (Brauer com vetores de RefA contra
    os exatos);
  - os NK de S3b, S4 e da faixa B, e os testemunhos 0,1685 e 0,3547;
  - o certificado em Arb de F3 (‖B_L‖ ≤ 3,9641, R_L ≤ 2,926) e o de R♯ = 1,765;
  - a cobertura de S3a.

  De F3 só conferi que a norma é a mesma (§2.1).
- `SHARP_DOMAIN.md` (D_min = D_max; perda de raio, com as constantes 8385 e 1440): li e achei coerente, não refiz as
  constantes.
- A origem geométrica de g* em `GAUGE_RT_ACTION.md` (derivada de Lie da métrica). Só rodei os testes simbólicos de
  covariância (OK). Também não refiz a classificação de Floquet do gauge além da leitura.
- As derivações dos majorantes `selected_majorant` e `gamma1_majorant` (β+). Reexecutei, não rederivei.
- `signed_operator_L.validar` e `signed_constraints.validar` contra o exportador de RT. Eles leem `build/spectrum/`,
  fora da lista de leitura permitida; usei a validação como dado.
- ‖B_L‖ truncado acima de 80×240, e o sup do símbolo (só o certificado de F3 o cota).
- As hipóteses físicas, por natureza.

## 5. Estado dos cálculos longos (fechado em 03/10, 01:30)

- ‖B_L‖ truncado (`art/r1_BL.txt`): **terminou**. 3,6935 (12×36), 3,8208 (24×72), 3,8501 (40×120), 3,8593 (60×180)
  e 3,8626 (80×240, dim 66.600). A sequência é crescente e converge abaixo da cota de F3, 3,9641.
- `l_real_toro.py`:
  - **notebook: terminou**, com JSON idêntico ao gravado;
  - **PC 3: incompleto** (só o primeiro μ, idêntico), porque o PC ficou offline no meio. O processo não existe
    mais lá; não sobrou nada rodando.
- Todos os outros scripts (`art/r1_toro`, `r1_toro_convergencia`, `r2_*`, `r5_iso`, `r6_quociente`, `r7_*`,
  `r8_contagem_discos`) terminaram, com as saídas gravadas.
- Nenhum arquivo do projeto fora desta pasta foi alterado: `MANIFEST.sha256` confere sem falhas no fim.
