# C1/T1 reavaliada à luz de S3b e S4 (26/09/2026)

**Pergunta.** A transferência da contagem finita para o operador infinito (C1, C1-G, C2, C3 e T1)
estava formulada, em 20/09, como "majorante do exterior < 1". Isso exigia β de 10,382 → 2,1 e caixas
de 556×3001, o que levou à discussão de HPC (`PRIORIDADES.md`, itens 1 a 4). A maquinaria de S3b e S4
muda esse quadro?

**Resposta curta: sim.** A contagem pode ser feita por **Rouché de operadores com referência finita**,
a mesma rota que certificou K em S3b (`rouche_K.py`), com o tratamento de cauda de L em S3b e S4
(bloco-Jacobi em janelas, Loewner, Perron de 3 níveis). As medidas abaixo estimam ~250 tiles e
~100–150 horas de máquina nas máquinas que já temos, sem HPC. Tudo aqui é diagnóstico em ponto
flutuante (`scripts/experiment_contour_L.py`), **não é prova**. O que falta decidir está no piloto
(§5).

## 1. A formulação

- **Rouché de operadores (Gohberg–Sigal).** Com P = Q0(c0) fixo, H(s) = L(s)P tem os mesmos zeros
  que L. Toma-se a referência A(s) = diag(Ĥ_ZZ(s), W, I), em que:
  - Ĥ_ZZ(s) = I + P_Z(B + s − c0)Q0(c0)P_Z;
  - W é bloco-diagonal, com os blocos das janelas de cauda num centro fixo, logo constante em s e
    sem zeros.

  Se ‖A⁻¹(H − A)‖ < 1 em todo o contorno, H e A têm o mesmo número de zeros dentro dele. Isso é
  exatamente o Perron de 3 níveis de `nk_L3_C.py`, sem a borda.
- **Zeros da referência = autovalores da truncagem.** Q0 é triangular em n, então
  P_Z Q0 P_Z = J_Z⁻¹ e Ĥ_ZZ(s) = L_Z(s) J_Z(c0)⁻¹ exatamente. Os zeros de A dentro do contorno são
  os autovalores de L_Z ali.
- **Contagem finita.** Tira-se a decomposição de Schur de L_Z uma vez e certifica-se o resíduo.
  Um Rouché trivial (resíduo ~10⁻¹⁰ contra o resolvente já certificado nos tiles) conta os
  autovalores pela diagonal. Isso dispensa o winding exato por CRT, que é viável só em 12×36.
- **Deflação das raízes de gauge (Brauer).** A raiz de fase s = 0 fica na borda Re s = 0, e a
  vizinha dela, μ, deixa o resolvente alto entre as duas. Com v = ∂τω* (autovetor exato de fase,
  controlado pelo lema de Corr de S4) e φ qualquer com φ*v ≠ 0, L̃(s) = L(s) + σvφ* satisfaz
  det(I + σL⁻¹vφ*) = (s + σφ*v)/s. Então N(L̃, Ω) = N(L, Ω) − [0 ∈ Ω] + [−σφ*v ∈ Ω]; as outras
  raízes e multiplicidades ficam intactas. Idem para μ, com g*. A diferença entre o vetor
  aproximado e o exato (~10⁻⁸) entra como perturbação, como o erro do fundo.
- **Contorno.** ∂([0, R] × [−1/4, 1/4]), com R = 3 (F3: nada em Re s >= 2,93). Só a metade
  superior é preciso, pela simetria de conjugação. Se N(L̃) = 2, então N(L; Re s > 0, faixa A) = 3,
  **sem faixa descoberta junto ao eixo imaginário**: a lacuna "0 < Re s < 1/8" de C4 desaparece.
  Isso une C1, C1-G, C2, C3 e T1 num único certificado e dispensa a comparação entre a resolvente
  comprimida e a inversa finita (C1-G).

## 2. Medidas

**Raízes da truncagem** (12×36 e 16×64 concordam em 5 casas):

- faixa A, Re s > 0: μ = 0,16831, s_K = 0,40102 e **0,73318** (a raiz física);
- na borda: 0 (fase);
- à esquerda: −0,168 e −0,175;
- faixa B: 0,04142 ± 0,5i.

**O nível Z converge cedo na norma ℓ² de S3b.** ‖L_Z(s)⁻¹‖ nos pesos de L:

| s | 8×24 | 12×36 | 16×64 | 24×96 |
|---|---:|---:|---:|---:|
| 0,08 | 361,1 | 348,6 | 348,6 | 348,6 |
| 0,08 + 0,25i | 69,5 | 68,1 | 68,1 | 68,1 |
| 0,55 + 0,25i | 22,7 | 22,6 | 22,6 | 22,6 |
| 1 | 9,91 | 9,91 | 9,91 | 9,91 |

Na medida antiga (norma ℓ¹ de coluna, `S3_CONDICIONAMENTO_SHARP.md`), ‖V_L(1)‖ **crescia com a
malha** (183 → 214), e o máximo em Γ_A era 1846. Na norma de Hilbert com pesos com sinal, que é a
que o NK de S3b usou, a mesma grandeza é 15,3 e está convergida. Boa parte do "mau condicionamento
de L" era artefato da norma.

**Perfil ao longo do contorno** (16×64): ‖L_Z⁻¹‖ vale 349 → 68 na borda esquerda (Re 0,08, de
Im 0 a 1/4), 68 → 7 na borda superior de Re 0,08 a 1, e 7 → 0,7 de 1 a 3.

**Com a deflação de Brauer** (16×64, vetores numéricos):

| s | sem deflação | fase | fase + μ |
|---|---:|---:|---:|
| 0,04 | 541 | 174 | 126 |
| 0,08 | 349 | 191 | 94 |
| 0,12 | 394 | 296 | 73 |
| 0,1i | 214 | 170 | 152 |

Deflacionar também −μ piora (o par −0,168/−0,175 é mal condicionado).

**O acoplamento Z←T não segue o condicionamento.** Em S3b (`build/s3b_escala_L.txt`), em
s0 = 0,276, sem borda e com ‖V‖ = 97, a vale 0,71, 0,58 e 0,50 em Z = 16×64, 20×80 e 24×96. Em
0,401, com borda e ‖V‖ = 21, a vale 0,48 em 24×96. Então a depende do tamanho de Z, não da distância
às raízes. Em 32×128 deve ficar em ~0,4 também no contorno, e θ(s) ~ 0,85–0,9 como nos NK.

## 3. Custo estimado

- **Raio de tile.** θ(s) varia com s pelo resolvente de L_Z, então basta |s − c| · ‖L_Z(c)⁻¹‖ <= ~0,15
  para manter θ < 1. O número de tiles é ~∫‖L_Z⁻¹‖ |ds| / 0,3. Na metade superior, com a deflação,
  dá **~200–250 tiles**. Cada tile pede o nível Z certificado (‖Ĥ⁻¹‖ e a por Loewner em 14.016,
  ~15–20 min na máquina do laboratório): ~70 horas de máquina.
- **Cauda.** Os blocos de janela ficam num centro c0 por trecho. A perturbação
  ‖V_J‖ · |s − c0| · ‖Q0 P_J‖ ≈ 0,3 · |s − c0| pede centros a cada ~0,3 ao longo do contorno
  (comprimento ~3,5): ~10 centros × ~4 horas de máquina (a etapa B de S4 levou 2 h em duas
  máquinas).
- **Total:** ~100–150 horas de máquina, ~2–3 dias no laboratório e no laboratorio2.

## 4. Riscos e pontos em aberto

- **θ(s) fora do eixo real ainda não foi medido.** Nos NK, θ = 0,84 e 0,88 com borda. No contorno a
  referência não tem borda. É o que o piloto decide.
- A tolerância real da cauda a |s − c0|, que define quantos centros.
- A certificação dos vetores de deflação: ∂τω* sai do lema de Corr de S4 (‖∂τCorr‖ é controlado);
  g* já está certificado.
- As bordas Im s = ±1/4 são compartilhadas com a faixa B (C4). O mesmo método se aplica lá, onde
  fica a raiz 0,0414 + 0,5i, que viola as constraints (S3 na faixa B).
- A contagem física continua pedindo S3 na região inteira: K injetivo em 1,765 < Re s < R (F3 de K)
  e na faixa B.

## 5. Piloto proposto

Adaptar `nk_L2_Z.py` (etapa A) e `nk_L3.py` (estimativa) para s complexo, sem borda e com deflação
de Brauer opcional. Rodar a estimativa de θ em três pontos do contorno:

- s = 0 com a fase deflacionada (o pior canto);
- s = 0,08 + 0,25i;
- s = 0,55 + 0,25i.

São ~2 h por ponto por máquina. Se θ < ~0,93 nos três, a rota está aberta e o dimensionamento acima
vale.

## 6. Resultado do piloto (26/09, manhã): a rota está aberta

Estimativa de θ do Rouché, sem borda, com deflação de Brauer da fase e de μ, em Z = 32×128 e
F = 200×400 (`build/piloto/`, uma máquina por ponto, ~2 h cada):

| ponto | máquina | ‖Ĥ_ZZ⁻¹‖ | a (Z'←T) | θ0 (pesos 1) |
|---|---|---:|---:|---:|
| s = 0 | laboratório | 188,8 | 0,413 | **0,886** |
| s = 0,25i | laboratorio2 | 107,8 | 0,419 | **0,888** |
| s = 0,55 + 0,25i | PC 3 | 34,8 | 0,366 | **0,811** |

- **Os três ficam abaixo de 0,93.** Nos NK de S3b e S4, o θ rigoroso ficou ~0,02 acima da
  estimativa (0,860 → 0,880), então o rigoroso deve dar ~0,90–0,91 no pior ponto.
- **a não segue ‖Ĥ⁻¹‖:** vale 0,41 com ‖Ĥ⁻¹‖ = 189. Confirma o §2.
- **A cauda quase não depende de s.** Na janela 32..47, ‖V‖ vale 1,78, 1,77 e 1,65, e o
  acoplamento a Z' vale 0,40, 0,39 e 0,34. Poucos centros de cauda devem bastar, menos que os ~10
  do §3.
- **Raio de tile.** Com margem ~0,1, o raio fica em |s − c| · ‖L_Z(c)⁻¹‖ <= ~0,1: ~6·10⁻⁴ em s = 0,
  ~10⁻³ em 0,25i e ~3·10⁻³ em 0,55 + 0,25i. Isso dá ~300 tiles na metade superior, ~90 horas de
  máquina no nível Z, ~30 h em três máquinas.

**A contabilidade da deflação, com 0 na borda.** Seja Ω = [0, R] × [−1/4, 1/4] e seja
L̃ = L + δ0 v φ0* + δ1 g φ1*, com v = ∂τω* e g = (Z* + 1)ω* autovetores direitos EXATOS em 0 e μ, e
φ0, φ1 quaisquer. O certificado usa v_A e g_A; a diferença, <= ~10⁻⁷ pelo lema de Corr e pela
etapa E de S4, entra como perturbação no Rouché. Indenta-se o contorno em torno de 0 por um
meio-círculo de raio ε pequeno, onde L e L̃ são invertíveis. Como L(s)⁻¹v = v/s e
L(s)⁻¹g = g/(s − μ), vale L̃ = L(I + L⁻¹E) com D(s) = det(I_2 + [δ_j φ_i* L⁻¹ x_j]) racional:
s(s − μ)D(s) = q(s) é um quadrático explícito (raízes ~ −1, fora de Ω). Então
N(L, Re s > 0, faixa A) = N(L̃, Ω) + 1 (o polo de D em μ) − #raízes de q em Ω. O Rouché certifica
N(L̃, Ω) = contagem de Schur da referência (esperado: 2, s_K e 0,733), logo N_A = 3. Um L̃(0)
invertível já certifica que 0 é raiz simples de L; se não fosse, o tile em s = 0 falharia.

**Próximo passo:** o certificado rigoroso (`rouche_L.py`). A partir das etapas A, B e C de S4,
são quatro peças:

1. etapa A sem borda e com a deflação certificada (vetores de RefA, Loewner) por centro de tile;
2. etapa B nuns poucos centros de cauda, com a perturbação |s − c0| explícita;
3. montagem por tile com o raio;
4. contagem de Schur da referência.

Em três máquinas, são ~2–3 dias.

## 7. Desenho do certificado (26/09, tarde)

**O problema do tile.** O Rouché precisa da desigualdade em TODO ponto do contorno. Qualquer cota
que multiplique normas perde a estrutura que faz a pequeno:

- a = ‖V_Z Y‖ vale 0,4, mas ‖V_Z‖ · ‖Y‖ ~ 95;
- em 12×48, s = 0: a variação verdadeira é Δa = 0,02 em |ε| = 0,01, e a quantidade estruturada
  de primeira ordem vale ‖(V − I)RY‖ = 4,8 (`experiment_tile_variacao.py`);
- a forma usada para K (linha Z dividida por 1 − e, e = |ε| ‖V_c Q_ZZ‖) dá raio ~3·10⁻⁴ em
  s = 0, ~1300 tiles de ~20 min: inviável;
- as duas metades do contorno precisam de certificado: a conjugação leva o espaço com sinal κ1^m
  ao de κ1^{−m}.

**A saída: tile barato, e não tile grande.**

1. **Referência global holomorfa.** H(s) = L̃(s)Q0(s) com Q0(s) = J(s)⁻¹, e
   A(s) = diag(Ĥ_ZZ(s), Ĥ_J(s) nas 26 janelas, I). Na forma de homotopia,
   H_t = A + t(H − A) invertível para t ∈ [0, 1] sai do mesmo Perron com os inversos aproximados
   diag(V_Z(s), V_J(c_k), I): as entradas fora da diagonal só crescem com t.
2. **Escala s-dependente no nível Z.** Com D_Z = J_Z(c_k) Q_ZZ(s), a entrada T←Z fica
   exatamente ‖V_J P_J B Q_ZZ(c_k)‖, independente de s. A entrada Z←T vira
   ‖J_Z(c_k)[R(s)Ỹ − Q_ZZ(s)J_ZT] Q0_TT(s)‖, com R(s) = L̃_Z(s)⁻¹, Ỹ = P_Z L̃(0) P_T e
   J_ZT = P_Z J P_T, ambos constantes (usa J_ZZ Q_ZT + J_ZT Q_TT = 0). Q0_TT(s) custa só o fator
   1/(1 − |s − c_k| q_TT).
3. **Nível Z rápido.** A decomposição de Schur de L̃_Z (float) é feita uma vez e serve só como
   solver. O certificado é a posteriori: Ŝ = U(T + s)⁻¹U* F, e o resíduo (L̃_Z + s)Ŝ − F sai de
   um produto com cota de arredondamento. O acoplamento vira posto baixo: o Gram de
   [[Ỹ Q0_TT(c_k)], [J_ZT Q0_TT(c_k)]] = F F* + resto, com o resto certificado uma vez por centro
   de cauda. Cada ponto do contorno custa resoluções triangulares com k colunas, segundos.
4. **Entre pontos.** Primeira ordem estruturada (‖J_Z R(p)·R(p)F‖, também por Schur) mais
   segunda ordem grosseira com a cota de resolvente. A cota de resolvente vem de uma grade grossa
   de σ_min certificados (Loewner), usando que σ_min(L̃_Z + s) é 1-Lipschitz em s.
5. **Cauda.** Centros c_k com as quantidades da etapa B, mais ρ_J(s) = ‖I − V_J(c_k)Ĥ_J(s)‖
   limitado por |s − c_k| vezes uma constante calculada no centro. O espaçamento sai do
   diagnóstico de sensibilidade.
6. **Contagem.** Os autovalores de L̃_Z em Ω vêm da diagonal de Schur: um Rouché finito trivial,
   com a cota de resolvente no contorno, cobre o resíduo. Depois vem a contabilidade de Brauer
   (§6).

**Diagnósticos em curso** (26/09): sensibilidade da cauda (PC 3), posto do acoplamento em
32×128 (laboratorio2) e custo de Schur em 14016 (laboratório).

## 8. Implementação (26/09, tarde)

**Deflação.** Os autovetores de fase e de gauge são quase paralelos (cos 0,87). Com φ = autovetor
esquerdo, ‖E‖ ~ 36; com φ = v separado, as duas deflações se acoplam e deixam uma raiz em
−0,056, colada em Re s = 0 (resolvente ~190). A deflação **conjunta biortogonal** resolve as duas
coisas: φ_i ∈ span{v, g} com φ_i* x_j = δ_ij, e as duas raízes vão exatamente a −1, com ‖E‖ ~ 1,3.
No teste em 12×48, o resolvente cai para 36 em s = 0, 42 em 0,1i, 49 em 0,25i e 25 em 0,08.
No certificado, v e g são os vetores explícitos de RefA na caixa Z:

- ‖P_Z(v* − v_A)‖ <= ‖∂τRefB‖ + ‖∂τCorr‖ = 1,46·10⁻⁹ (`gauge_refB.py`, lema de Corr);
- ‖g* − g_A‖ <= 2,5·10⁻⁸ (etapa E de S4);
- as caudas fora de Z entram nas linhas de T.

**Peças.**

| script | função |
|---|---|
| `rouche_L_Z.py` (`schur`, `fator`, `pontos`) | Schur da truncagem deflacionada, posto baixo do acoplamento por centro, estimativa por ponto |
| `rouche_L_rig.py` (`globais`, `ponto`) | nível Z rigoroso por ponto: t_p por Loewner, X_k = R(p)^k F_t por Schur com resíduo certificado, a0, b1, b2, ‖Y3‖ |
| `nk_L3_T.py --rouche` | cauda rigorosa em λ0 complexo, com K_J = ‖V_J P_J B Q0(c)² P_J‖ e ‖Q0‖ por janela |
| `rouche_L_disco.py` | montagem por disco: maior raio com θ <= alvo, com o Perron em racionais |

**Sensibilidade da cauda, rigorosa.**

- Janelas não centrais: ρ_J(s) <= ρ_J(c) + d·K_J/(1 − d·nQ_J).
- Janelas centrais: pela divisão Q0(s)P_J = descida para Z + bloco da janela e pelas identidades
  triangulares, aparecem só TZ_J, jq_J = ‖J_ZJ Q_JJ(c)‖ e q_ZJ(s) <= ‖Q_ZZ(p) J_ZJ Q_JJ(c)‖ com
  fatores de Neumann (d = |s − c|).

**Execução** (três máquinas, filas em andamento):

1. estimativa da cauda em 10 centros;
2. nível Z: Schur, fatores, Schur final e 72 pontos de estimativa;
3. cauda rigorosa nos 10 centros;
4. nível Z rigoroso nos pontos (~9 min por ponto);
5. montagem por disco e cobertura do contorno;
6. contagem de Schur e contabilidade de Brauer.

## 9. Colagem e contagem (30/09)

**Uma lacuna no desenho do §7, e a correção.** A referência do §7 usa os blocos de janela
**dependentes de s**, Ĥ_J(s) = P_J H(s) P_J, com os inversos aproximados V_J(c_k) fixos por trecho.
A homotopia H_t = A + t(H − A) dá N(H, Ω) = N(A, Ω), mas N(A, Ω) = N(Ĥ_ZZ, Ω) + Σ_J N(Ĥ_J, Ω), e o
Perron no contorno não diz nada sobre zeros de Ĥ_J(s) **no interior** de Ω. (No §1, com W constante,
a questão não existia.) A correção é o princípio do máximo, em `rouche_L_janelas.py`:

- Ω = [0, 3] × [−1/4, 1/4] é particionado em 10 retângulos Ω_k, um por centro de cauda c_k.
- Em Ω_k, F_k(s) = I − V_J(c_k)Ĥ_J(s) é analítica: Q0(s) é analítica em Re s > −μ, pois J_livre é
  triangular com diagonal μ(n + 1) + s + im/2. ‖F_k‖ é subharmônica, então
  sup_{Ω_k}‖F_k‖ <= max_{∂Ω_k}‖F_k‖, e ‖F_k‖ < 1 dá Ĥ_J(s) invertível em Ω_k.
- No bordo de Ω_k vale a mesma sensibilidade da montagem por disco:
  ρ_J(s) <= ρ_J(c) + d[TZ_J q_ZJ(s) + (K_J + TZ_J q_ZJ(c))(1 + d q_JJ(s))] nas janelas centrais, e
  ρ_J(c) + d K_J/(1 − d nQ_J) nas outras.
- q_ZJ(s) e q_JJ(s) são cotados em discos em torno de pontos de amostra do bordo, por Neumann no bloco
  triangular: Q_ZJ(s) = (I + εQ_ZZ(q))⁻¹Q_ZJ(q)(I + εQ_JJ(q))⁻¹.

**Resultado** (`build/rouche_L_rig/janelas.{txt,json}`, 89 pontos de amostra, 57 min no notebook):
todas as 10 regiões passam, com o máximo de ρ_J no bordo em **0,149** (centro 0, janela −15..0, em
0,3 + 0,1i) contra 1 exigido. À esquerda (Re <= 0,3) o máximo fica em 0,105–0,149; no meio, em 0,075;
à direita, em 0,028–0,059. Os blocos de janela estão longe de singulares em toda a faixa, logo
N(A, Ω) = N(Ĥ_ZZ, Ω).

**Contagem finita** (`rouche_L_contagem.py`). A decomposição de Schur em ponto flutuante, A ~ UTU*,
dá B = UTU⁻¹, cujos autovalores são exatamente a diagonal de T. Com
‖A − B‖ <= ‖AU − UT‖_F/sqrt(1 − ‖I − U*U‖_F) =: ε_B e ‖(A + s)⁻¹‖ <= t(s) no contorno (dos discos:
t_p/(1 − r t_p)), a homotopia A + t(B − A) é invertível no contorno se max t · ε_B < 1. Então A tem
em Ω exatamente #{i : −T_ii ∈ Ω} autovalores. No teste 12×48, ε_B = 2,4·10⁻⁷, e as duas raízes em Ω
são 0,401 e 0,733.

**Cobertura** (`rouche_L_cobertura.py`): em aritmética racional, a união dos discos certificados
cobre as quatro arestas de ∂Ω (a meia-corda de cada disco numa aresta é cotada por baixo).

**Contabilidade de Brauer.** Com φ_i* x_j^A = δ_ij (biortogonal), δ = (1, 1 + μ_A) e
e_ij = φ_i*(x_j − x_j^A), vale |e_ij| <= ‖φ_i‖ · ‖P_Z(x_j − x_j^A)‖ <= 5,5 · 2,5·10⁻⁸ (φ_i fica na
caixa Z). Então
q(s) = (s + 1 + e_00)(s − μ* + δ_1(1 + e_11)) − δ_0 δ_1 e_01 e_10,
e em Re s >= 0 os dois fatores têm módulo >= 1 − 10⁻⁶ (|μ* − μ_A| <= 3,9·10⁻¹¹, `DELTA_MU`). Logo q
não tem zeros em Ω, e N(L; 0 < Re s < 3, |Im s| < 1/4) = N(L̃, Ω) + 1.

A matriz A.npy é a mesma nas três máquinas (sha256 70791ce0…), então a referência é uma só em
todo o contorno. Os fatores F e as decomposições de Schur de cada máquina só servem de solver.

**Certificado do fator de posto baixo** (`rouche_L_Z.py certF`, 30/09). O `fator2` usava identidades
exatas sobre quantidades de ponto flutuante: S_J = F_J W* (SVD em float), [C_B; W*][C_B; W*]* = LL*
(eigh com autovalores cortados em 0) e o resto calculado com S em float, sem o erro de Q nem o
arredondamento do produto. São efeitos de ~10⁻¹³ relativos, mas não estavam certificados. O
certificado novo não depende disso. Para **qualquer** Φ vale ‖XS‖ <= ‖XF‖·‖Φ‖ + ‖X‖·ρ, com
ρ = ‖S − FΦ‖. Toma-se Φ = V f(Σ) U* S (Tikhonov sobre a SVD de F; Φ é só um dado), e ‖Φ‖ e ρ são
acumulados numa passada pelos modos de T, com S exata:
- Q0 com erro dq pelo resíduo do modo;
- arredondamentos elemento a elemento, γ‖|X||Y|‖_F;
- q_TT certificado nos mesmos modos.

Vários λ saem na mesma passada, e a montagem usa o melhor par por disco. Nos 10 centros:

| centro | ‖Φ‖ (λ ~ 3·10⁻⁷ σ_max) | ρ | resto_F antigo | q_TT |
|---|---:|---:|---:|---:|
| 0 | 1,00397 | 9,290·10⁻⁵ | 9,215·10⁻⁵ | 0,3396 |
| ±0,25i | 1,00398 / 1,00397 | 9,28·10⁻⁵ | 9,21·10⁻⁵ | 0,3395 |
| 0,55 ± 0,25i | 1,00399 / 1,00401 | 9,085·10⁻⁵ | 9,062·10⁻⁵ | 0,2878 |
| 1,0 | 1,00505 | 8,942·10⁻⁵ | 8,933·10⁻⁵ | 0,2558 |
| 1,5 | 1,00558 | 8,787·10⁻⁵ | 8,789·10⁻⁵ | 0,2275 |
| 2,0 | 1,00606 | 8,636·10⁻⁵ | 8,645·10⁻⁵ | 0,2047 |
| 2,5 | 1,00662 | 8,490·10⁻⁵ | 8,504·10⁻⁵ | 0,1861 |
| 3,0 | 1,00735 | 8,348·10⁻⁵ | 8,367·10⁻⁵ | 0,1707 |

O sha256 do F de cada certificado confere com o arquivo usado no nível Z da máquina de origem.

**Nível Z rigoroso.** Fase 1: Loewner com partida f0 = 1,005 + 0,13(est/45)². A margem de Rump cresce
com t², e t_p sai de 1,003 a 1,12 vezes a estimativa. Fase 2: Schur sem cópias n×n. A versão anterior
precisava de ~17 GB, mais do que os 14–15 GB das máquinas, e caía em swap. Os t_p certificados ficam
no log e podem ser reaproveitados, arredondados para cima, em outra máquina com a mesma A. Isso salvou
o centro −0,25i quando o PC 3 foi desligado.

**Resultado (01/10/2026).** Montagem por disco com o certificado do fator, alvo θ <= 0,999, Perron
em racionais (`build/rouche_L_rig/discos/`):

| centro | máquina | pontos | θ máximo | raios |
|---|---|---:|---:|---|
| 0 | laboratorio2 | 12 | 0,99915 | 0,0149–0,0227 |
| 0,25i | laboratorio2 | 18 | 0,99916 | 0,0156–0,0303 |
| −0,25i | PC 3 + laboratório | 16 | 0,99910 | 0,0156–0,0306 |
| 0,55 + 0,25i | laboratório | 15 | 0,95704 | 0,0360–0,0683 |
| 0,55 − 0,25i | PC 3 | 15 | 0,95010 | 0,0368–0,0700 |
| 1,0 | laboratório | 10 | 0,88950 | 0,081–0,221 |
| 1,5 | laboratório | 4 | 0,88606 | 0,31–0,46 |
| 2,0 | laboratório | 2 | 0,80671 | 0,4995 |
| 2,5 | laboratório | 2 | 0,75692 | 0,4995 |
| 3,0 | laboratório | 2 | 0,72136 | 0,4995 |

- **Cobertura** (`cobertura.txt`): os 96 discos cobrem as quatro arestas de ∂Ω.
- **Janelas** (`janelas.txt`): nenhum bloco de janela tem zero em Ω (máx ρ_J no bordo 0,149).
- **Contagem finita** (`contagem.json`, `contagem_rouche.txt`):
  - ‖AU − UT‖_F <= 5,97·10⁻⁵, ‖I − U*U‖_F <= 1,30·10⁻⁶ e ε_B = 5,97·10⁻⁵;
  - max t(s) no contorno <= 560,3, logo t·ε_B <= 0,033 < 1;
  - zeros de A + s em Ω: **0,401025 e 0,733180**;
  - o autovalor mais próximo do bordo, fora de Ω, é −0,1683, a 0,168 dele.
- **Brauer:** q(s) não tem zeros em Ω, então N(L) = N(L̃) + 1.

**Conclusão:** N(L; 0 < Re s < 3, |Im s| < 1/4) = **3**, contado com multiplicidade algébrica:
μ* (S4), s_K = 0,401 (S3b) e a raiz física 0,733. Com F3 (nada em Re s >= 2,93), isso vale para toda
a faixa A, Re s > 0 e |Im s| < 1/4, **sem a fita descoberta junto ao eixo imaginário**. Isso fecha C1,
C1-G, C2, C3 e T1 na faixa A, na realização de L com pesos com sinal e norma ℓ² do toro. É um
certificado computacional e ainda **precisa de revisão independente**, no protocolo de
`revisao-independente`.

Continuam abertas: a faixa B (C4, com a raiz 0,0414 + 0,5i), a parte física de S3 (constraints) e T2.

## 10. Revisão independente e resposta (01/10/2026)

**Revisão.** Sub-agente de contexto limpo, de outro modelo da família (Fable). Ficou proibido de ler
este documento e o ledger, recebeu as alegações A1–A9 enunciadas e entregou um artefato por alegação
(`.codex-runs/2026-10-01-revisao-c1/`: `TAREFA.md`, `RELATORIO.md`, scripts `rN_*` e saídas).
Snapshot em `.snapshots/2026-10-01-revisao-c1/`; o manifesto conferiu intacto antes e depois.

O que a revisão conferiu:
- **R1, lógica:** escreveu a prova de que θ < 1 implica H_t invertível; Gohberg–Sigal se aplica;
  rederivou a identidade de escala.
- **R2, LU densa:** 27 pontos s de 6 discos em 4 centros. Nenhuma violação: max ‖R‖/t(s) = 0,987 e
  max valor/E = 0,99999.
- **R3, certF próprio:** S montado do zero, Φ por mínimos quadrados; reproduz ‖Φ‖ e ρ.
- **R4, sensibilidade:** fórmulas rederivadas e testadas; pior razão verdadeiro/cota de 0,981.
- **R6, janelas:** pontos recalculados.
- **R7, cobertura:** verificador próprio com os cantos.
- **R8, contagem:** ε_B próprio de 6,15·10⁻⁶, contra 5,97·10⁻⁵.
- **R9, Brauer:** min|q| >= 0,99999997.
- **R10, reprodução.**

**Achados, e o que foi feito.** Conferi cada um antes de aceitar.

| | achado | gravidade | conserto |
|---|---|---|---|
| L1 | a deflação verdadeira (x_k exatos contra x_k^A de A) faltava nas linhas de T e do far: φ*Q0(s) não se anula nas colunas de T e do far (Q0 triangular superior), e (I − P_Z)x_k ≠ 0 | baixa (Δθ ~10⁻⁴) | `theta_disco`: dET·nQx em T←T e far←T, dET·Qf_s em T←far e far←far, dET·‖Q_ZZ(c)‖ em far←Z'; `rouche_L_janelas`: nV_J·dET·q_ZJ(s) nas janelas centrais |
| L2 | Qf_s = Qf(c)/(1 − d Qf(c)) só vale no bloco far←far; as linhas de Z/T de Q0(s)P_far trazem (I + (s − c)Q_NN)⁻¹ | média (cota grosseira válida daria θ > 1) | `closed_form.cauda_Rinv_disco` e `Qfar_disco`: sup no disco das formas fechadas com majorantes, \|d_k\| >= max(\|d_k(p)\| − r, μ(n+1)) e \|1 − v_k/d_k\| <= min(1, …) (o 1 vale em Re s >= 0); com r = 0 reproduz `cauda_Rinv` a 10⁻¹¹ |
| L3 | certF de 0 e 0,25i sem sha256 do F | baixa | refeitos no laboratorio2; mesmos valores, com o sha |
| L4 | 7 rigs no PC do laboratório sem sha de A registrado | baixa a média | rig v2 grava sha256 de A e de F; os 96 pontos refeitos no laboratorio2 e no PC 3 (t_p reaproveitados); a montagem confere sha(F) do rig = sha(F) do certF e sha(A) = 70791ce0… |
| L5 | arredondamentos: fl(A + p) antes da Loewner; Q_ZZ^k F_b com normas espectrais onde a cota elemento a elemento pede \|Q\| | cosmética | rig v2: t_p/(1 − t_p u max\|A_ii + p\|); γ_dim‖\|M\|\|X\|‖_F calculado nos produtos por blocos |
| L6 | 10⁻⁶⁰ para ‖∂τCorr‖ sem derivação | — | `corr_gauge_bound.py`: Corr = O⁻¹G, logo ‖∂τCorr‖_L <= √2 c_τ‖G‖ = 1,39·10⁻⁷⁸ |

**Achados meus na conferência** (do mesmo tipo das lacunas, todos consertados):
- ‖Q0(s)P_J‖ estava tomado no centro nos termos de perturbação global. Agora usa
  ‖Q0(s)P_J‖ <= (‖Q0(c)P_J‖ + d‖Q_ZZ(s)‖‖Q_ZJ(c)‖)fJ, pela inversa em blocos de I + (s − c)Q0(c).
- `eps_mu_L` já cota o operador inteiro δμ(J1 + SΓ1)Q0(s). Multiplicá-lo por ‖Q0 P_J‖ < 1 subestima o
  termo; nos termos novos o fator é max(1, ‖Q0(s)P_J‖). Na maquinaria de S3b/S4, o termo global
  (rt + ε_μ)·nQmax com nQmax ≈ 5 cobre a diferença, de ~10⁻⁹.
- J_Z(0) e J_Z(c) são calculadas em ponto flutuante: 4u‖J‖_F entra em ‖B̂‖ e nas cotas a0, b1 e b2.
- O arredondamento de J_Z(c)·D passou a ser cotado elemento a elemento.
- A entrada far←J não ganha o termo de Z da inversa em blocos: P_far B P_Z = 0 na banda (Z vai até
  n < 128 e |m| < 32; a banda em n é 100 e o far começa em n = 400 ou |m| = 200). O resto fora da banda
  já está na perturbação global.

**Resultado da versão 2** (`build/rouche_L_rig/discos_v2/`, `z2/`, `janelas_v2.*`, `cobertura_v2.txt`,
`contagem_rouche_v2.txt`):

| centro | θ máximo v2 | antes | raio mínimo |
|---|---:|---:|---:|
| 0 | 0,99797 | 0,99915 | 0,01485 |
| 0,25i | 0,99841 | 0,99916 | 0,01559 |
| −0,25i | 0,99837 | 0,99910 | 0,01564 |
| 0,55 + 0,25i | 0,95721 | 0,95704 | 0,03596 |
| 0,55 − 0,25i | 0,95027 | 0,95010 | 0,03677 |
| 1,0 | 0,87846 | 0,88950 | 0,08081 |
| 1,5 | 0,87800 | 0,88606 | 0,31339 |
| 2,0 / 2,5 / 3,0 | 0,80017 / 0,74063 / 0,70321 | 0,807 / 0,757 / 0,721 | 0,4995 |

- **Discos:** os 96 passam, com θ <= 0,99842. A cota do far no disco é mais justa que a antiga e
  compensa com folga os termos novos.
- **Cobertura:** confirmada pelo verificador do autor e pelo do revisor.
- **Janelas:** máx ρ_J = 0,1488 (os termos novos são de ~10⁻⁷).
- **Contagem finita:** max t <= 559,2 e t·ε_B <= 0,0334.

**Conclusão, revisada:** N(L; 0 < Re s < 3, |Im s| < 1/4) = 3, com multiplicidade algébrica.

**Não verificado pela revisão:**
- a maquinaria já revisada em S3b e S4 (Loewner, cota_norma, nk_L3_T, dB_L, formas fechadas, eps_mu_L);
- a bola de Corr de RT;
- a reexecução das caudas T3.

O revisor remontou θ em 7 dos 96 discos; a reexecução v2 recalculou todos. A localização do terceiro
zero (0,733) vem da contagem de Schur e não do Rouché; o Rouché garante só que são 3.

## 11. Faixa B (C4): fechada (02/10/2026)

**Formulação.** Pela redução do setor B ao A (`SECTOR_B.md`, Σ_B = Σ_A − i/2), a faixa B é
Ω_B = (0, 3) × (1/4, 3/4) para L_A. F3 (campo de valores, Re W(−B_L) <= 2,95) não depende de Im s.

**Contagem: o mesmo Rouché, com a mesma referência deflacionada.**
- A deflação de Brauer move 0 e μ para as raízes de q (≈ −1), e os polos de D(s) ficam em 0 e μ.
  Tudo isso está fora de Ω_B, logo N(L, Ω_B) = N(L̃, Ω_B).
- A referência A(s) = diag(Ĥ_ZZ(s), Ĥ_J(s), I) é uma só função analítica. Então os discos da faixa A
  certificam a aresta Im = 1/4.
- A diagonal de Schur de A (`build/rouche_L_rig/diagT_lab2.npy`) tem **um** autovalor em −Ω_B:
  0,0414194105 + 0,5i.
- Perto de Im = 3/4 ficam os transladados por i de μ, s_K e 0,733 (periodicidade m → m + 2), a 0,25 da
  aresta.

**Resolvente e a Loewner pré-condicionada.**
- Na aresta esquerda, ‖(A + s)⁻¹‖ chega a 201,9 em 0,5i (a raiz fica a 0,041), com o mesmo valor em
  12×48 e em 32×128. No canto 0,75i, 86.
- A Loewner direta perde margem como 2n·u·tr(t²AA*) e não fecha acima de t ~ 70.
- **Pré-condicionada** (`rouche_L_rig.tp_precond`): (A + p)⁻¹ = K⁻¹Q, com Q = Q_ZZ(p) exata e
  K = I + QB̂, e t²KK* − QQ* ⪰ 0. Em 32×128 sai t = 1,003 vezes a estimativa, na primeira tentativa,
  até t = 202.
- As fases 1 e 2 do nível Z rodam em processos separados (memória).

**Centros.**
- Novos: 0,5i, 0,75i e 0,55 + 0,75i. Os dois primeiros com cauda rigorosa no laboratorio2 e no PC do
  João, com checkpoints juntados; o terceiro no PC do laboratório.
- Reaproveitados da faixa A: 0,25i e 1,0 a 3,0. As caudas dos centros reais toleram d até ~1,25 na
  aresta Im = 3/4.

**Resultado** (`build/faixaB/`; rigs v2 com sha de A = 70791ce0… e de F conferidos contra o certF):

| centro | discos | θ máximo | raios | t_p máximo |
|---|---:|---:|---|---:|
| 0,25i | 6 | 0,99778 | 0,0118–0,0181 | 66 |
| 0,5i | 28 + 2 | 0,99907 | 0,0028–0,0090 | 202 |
| 0,75i | 23 | 0,99884 | 0,0077–0,0191 | 86 |
| 0,55 + 0,75i | 9 + 2 + 2 | 0,95810 | 0,030–0,080 | 30 |
| 1,0 | 1 | 0,98057 | 0,1227 | 7,3 |
| 1,5 | 2 | 0,97831 | 0,23–0,36 | 3,9 |
| 2,0 / 2,5 / 3,0 | 1 / 1 / 2 | 0,912 / 0,828 / 0,772 | 0,4995 | <= 1,3 |

- **Cobertura** de ∂Ω_B: 174 discos, sendo 78 da faixa B (79 entradas, uma com r = 0) e 96 da faixa A, na aresta Im = 1/4. Confirmada
  pelo verificador do autor e pelo do revisor adaptado, com os quatro cantos (`cobertura.txt`).
- **Janelas** (`janelas.txt`, 10 regiões, princípio do máximo): máx ρ_J = 0,1517 (centro 1,0, em
  1,1 + 0,75i).
- **Contagem finita** (`contagem_rouche.txt`): 1 zero de A + s em Ω_B, max t <= 559 e
  t·ε_B <= 0,0334.

Logo **N(L; 0 < Re s < 3, 1/4 < Im s < 3/4) = 1**, contado com multiplicidade. Com F3 e os bordos
certificados, vale para toda a faixa B com Re s > 0.

**Testemunho: a raiz da faixa B viola as constraints** (`build/faixaB/nk/`).
- **NK bordejado** de L em λ0 = 0,0414194105346554 + 0,5i. Esse λ0 é o autovalor da truncagem 32×128
  **sem** deflação, obtido por Newton bordejado na matriz pesada (resíduo 2·10⁻¹⁵).
- **Dois consertos** foram necessários para fechar:
  1. h0 por iteração inversa **na matriz pesada**. A sem pesos tinha erro absoluto ~u em todas as
     componentes, que os pesos √ω(n)·κ1^m (~10¹²) amplificavam: Y_Z caiu de 2·10⁻⁴ para 2,8·10⁻⁶.
  2. tQ >= ‖Mh⁻¹[Q_ZZ; 0]‖ certificado (`nk_L2_Z --so-tQ`) no lugar de t_V‖Q_ZZ‖: 140,9 contra 339,4.
- **Resultado:** θ = 0,9115, Z2 = 155, Y = 6,0·10⁻⁶, r = 7,8·10⁻⁵. Há um zero único (y*, λ*) com
  |λ* − λ0| <= 3,7·10⁻⁵ e DF(x*) invertível, logo λ* é algebricamente simples.
- **Etapa D:** ‖Π_F C(λ*) h*‖ >= 0,35519 − 4,95·10⁻⁴ = 0,35469 > 0.

Como N(L, Ω_B) = 1, λ* é **a** raiz da faixa B. Ela é simples, e o autovetor viola as constraints.
**Multiplicidade física na faixa B: 0** (C4).

**Faixas A e B juntas.**
- O espectro é i-periódico, e A ∪ B cobre uma faixa de largura 1 em Im s, com os bordos certificados.
- Em Re s > 0 ficam 4 raízes: μ (gauge, S4), s_K (viola as constraints, S3b), 0,733 e a raiz da
  faixa B (viola as constraints).
- Pela redução, a raiz 0,733 é física se K for injetivo em 0,733, o que S3a cobre (Re <= 1,765).
- Contagem física total: **1**. Ela é condicionada às hipóteses de S4 (H-rec, H-gauge, H-cone) e às da
  redução do setor B (`SECTOR_B.md` §6).
- Como não há outras raízes de L, a injetividade de K em 1,765 < Re s < R deixa de ser necessária.
- Revisão independente: §11.1.

### 11.1 Revisão independente da faixa B (02/10/2026)

**Revisão.** Sub-agente de contexto limpo, de outro modelo (Fable). Ficou proibido de ler este documento e o
ledger, recebeu as alegações B1–B9 enunciadas e entregou um artefato por alegação
(`.codex-runs/2026-10-02-revisao-B/`). Snapshot em `.snapshots/2026-10-02-revisao-B/`; os manifestos
conferem antes e depois (local 413/413, laboratorio2 60/60, João 21/21). Usou só o laboratorio2 e o PC
do João. **Não conseguiu quebrar o certificado.**

| R | veredicto | o que conferiu |
|---|---|---|
| R1 | OK com ressalva | i-periodicidade (erro 1,1·10⁻¹⁶); F3 não depende de Im s; formas fechadas até \|Im\| = 3; Brauer; mesma A e mesmas janelas nas duas faixas. Os números de F3 não estavam entre os arquivos permitidos. |
| R2 | OK | Loewner pré-condicionada contra LU densa e Lanczos em 3 pontos perto de 0,5i: verdadeiro 201,148, 176,557 e 150,045; t_p 201,855, 177,169 e 150,556 (+0,35%). |
| R3 | OK com ressalva | 12 arquivos de nível Z: sha de A e de F, versão 2, pontos e t_p conferidos. a0 por LU própria: 0,42524 contra 0,42532. |
| R4 | OK | As 8 janelas de 0,5i calculadas nas duas máquinas são iguais bit a bit; o T3 juntado tem as 26. ‖V_J‖ por SVD fica 0,13% abaixo da cota. |
| R5 | OK | 3 discos remontados em racionais com código próprio: θ = 0,99885, 0,98080 e 0,99860. |
| R6 | OK | Cotas das janelas reproduzidas; ‖I − V_J Ĥ_J(s)‖ direto em Im = 0,75: 0,11–0,12. |
| R7 | OK | Verificador racional próprio, com os cantos. Shift-invert acha exatamente 1 autovalor em −Ω_B. |
| R8 | OK com ressalva | λ0 a 2,7·10⁻¹¹ do autovalor de A; λ* em Ω_B; Z2 com tQ correto; c_ref = 0,355192 contra 0,355189. Pelas matrizes de constraints exportadas por RT: ‖Ch‖/‖h‖ ≈ 0,50 na raiz B e 0,32 em s_K, contra ~10⁻⁴ em 0,733, 0 e μ. |
| R9 | OK (condicional) | A lógica dim E − rank(D\|E) está correta; a contagem física 1 depende das hipóteses abaixo. |
| R10 | OK | Cobertura e contagem reproduzidas exatamente (174 discos, max t = 559,151, t·ε_B = 0,0334, 1 zero). |

**Lacunas, todas de rastreabilidade, e o que foi feito:**
- **L1:** são 79 entradas, mas uma (0,291688 + 0,75i, centro 0,75i) tem raio 0. São **78 discos
  efetivos**, e 78 + 96 = 174 na cobertura. Texto corrigido; a entrada fica, porque a cobertura descarta
  r = 0.
- **L2:** a fase 2 sobrescrevia o `tp_*.json` da fase 1, com `precond` falso. Corrigido em
  `rouche_L_rig.py`: a fase 2 com todos os t_p reaproveitados não grava. Os logs, com a origem de cada t_p,
  estão arquivados em `build/faixaB/logs/` com `MANIFEST_logs.sha256`.
- **L3:** o `--alvo` não ficava no JSON do disco. `rouche_L_disco.py` agora grava os argumentos. Todos os
  discos das faixas A (v2) e B foram montados com `--alvo 0.999`. Onde o raio bate no teto r·t_p < 0,9, o
  resultado não depende do alvo; daí a reprodução exata do revisor com 0,998.
- **L4, fechada:** `T3_0_55p0_75.json` não tinha manifesto remoto: foi calculado no PC do laboratório,
  que estava desligado na hora da conferência. A cauda foi **recalculada no laboratorio2** (3,5 h). Os
  dois T3 coincidem a 5·10⁻¹⁶ relativo, que é o arredondamento esperado com 3 threads de BLAS contra 4.
  - Os 13 discos de 0,55 + 0,75i foram remontados com o T3 do laboratorio2: raios idênticos, θ igual a
    6·10⁻¹⁶.
  - A região 4 das janelas foi refeita: 0,0800.
  - Cobertura e contagem foram reconferidas.
  - O T3 do PC do laboratório ficou como `T3_0_55p0_75_labpc.json`.

**Hipóteses de que depende a contagem física 1** (lista do revisor):
- **Redução do setor B:** H-afim, H-iso, H-dof e H-constr (`SECTOR_B.md` §6).
- **Identificação física:** ansatz de RT, gauge residual e cone móvel (H-rec, H-gauge e H-cone de S4).
- **S4:** μ é gauge.
- **Espaço:** o espectro no espaço com pesos do toro é o espectro físico.
- **Dados:** ε = 4,7·10⁻⁸, DELTA_MU e F3.
- **Faixa A:** K injetivo em 0,733 (S3a).

**Não verificado pela revisão:**
- os números de F3;
- uma reimplementação de bloco_B e bloco_C;
- nas caudas: N_ij, TZ, K_J, nVhi e Nf;
- os certF dos centros novos;
- no nível Z: b1, b2, ‖Y3‖ e ‖Q³F_b‖;
- as etapas A e B do NK, tQ, e as perdas c_Z e c_J da etapa D.

## 12. NK em 0,7332: a raiz física localizada (03/10/2026)

Mesma cadeia do testemunho da faixa B (§11), com λ0 = 0,7331796596606924, o autovalor da truncagem refinado
por Newton bordejado (resíduo 1,6·10⁻¹⁵, parte imaginária 2·10⁻¹⁶ descartada). O vetor inicial h0 vem de
iteração inversa na matriz pesada. Artefatos em `build/nk733/`, com `MANIFEST.sha256`.

| grandeza | valor |
|---|---|
| truncagem Z / F | 32×128 / 200×400, n = 14016 |
| t_V = ‖V_b‖ | 8,85 (58 na raiz B) |
| a (Loewner) | 0,3788 |
| tQ = ‖Mh⁻¹[Q_ZZ; 0]‖ | <= 8,11045 (Loewner), contra t_V‖Q_ZZ‖ = 10,48 |
| Perron de 3 níveis (racionais) | θ = 0,805721 |
| Z2 / Y | 10,32 / 4,19·10⁻⁶ |
| raio | r = 2,158·10⁻⁵ |

**Resultado:** há um único zero (y*, λ*) com ‖x* − x0‖ <= r e **|λ* − 0,7331796596606924| <= 1,016·10⁻⁵**,
e DF(x*) é invertível, logo λ* é algebricamente simples. Como N_A = 3 e as outras duas raízes (μ*, S4;
s_K, S3b) estão em discos disjuntos deste, λ* é a terceira raiz da faixa A.
- Pela simetria de conjugação na realização de RT (L-real), λ* é real: λ* ∈ [0,73316950; 0,73318982].
- **Conferência externa:** com Δ ≈ 3,4453, γ = Δ/(4πλ*) ≈ 0,3739.

A etapa B (26 janelas da cauda) foi dividida entre o laboratorio2 (em ordem direta) e o PC 3 (em
ordem reversa, `etapaB_joao.txt`), com junção dos checkpoints no laboratorio2. O PC 3 foi desligado
no meio e retomou do checkpoint. A etapa C usou o T3 juntado (26/26 janelas).

## 13. O eixo imaginário: só a fase, simples (03/10/2026)

**Corolário dos Rouchés das faixas A e B.** Não há cálculo novo.

**Os dados.**
- A cobertura de ∂Ω_A, com Ω_A = (0, 3) × (−1/4, 1/4) (`build/rouche_L_rig/cobertura_v2.txt`), e a de
  ∂Ω_B, com Ω_B = (0, 3) × (1/4, 3/4) (`build/faixaB/cobertura.txt`), incluem a **aresta esquerda
  Re s = 0**.
- Nelas, cada ponto está num disco em que H̃(s) = L̃(s)Q0(s) é invertível, pela cota t_p certificada e
  θ < 1. t·ε_B < 1 é a condição da homotopia da contagem finita.
- s = 0 está no interior de dois discos, com folgas 0,018 e 0,005. `disco_0.json` é o arquivo do centro de
  cauda c = 0, não um disco centrado em 0. (Revisão leve, L2 e L3.)
- As duas arestas juntas cobrem Re s = 0 com Im s ∈ [−1/4, 3/4], um período inteiro de L_A, porque Σ_A é
  invariante por s ↦ s + i (§11).

**A multiplicidade.**
- L̃ = L(I + L⁻¹E), com E de posto 2 e vetores exatos v = ∂_τω* (raiz 0) e g = (Z* + 1)ω* (raiz μ) (§9). Então
  det(I + L⁻¹E) = q(s)/(s(s − μ)), e pela fórmula de Gohberg–Sigal para perturbações de posto finito,
  mult_L̃(s) = mult_L(s) + ord_s(q/(s(s − μ))).
- Em Re s >= 0, q não se anula: os dois fatores têm módulo >= 1 − 10⁻⁶, pelo §9.
- Q0(s) é invertível e analítica em Re s > −μ.

**Conclusão.**
- No eixo, para s ≢ 0: mult_L(s) = mult_L̃(s) = 0.
- Em s = 0: mult_L(0) = mult_L̃(0) + 1 = 1. **A raiz de fase é algebricamente simples**, como o §9 já
  antecipava ("um L̃(0) invertível já certifica que 0 é raiz simples de L").

Logo, **entre as raízes de L**, o eixo imaginário só tem {0} (mod i/2), o modo de fase ∂_τω*, que é de gauge
(n = 1).

**Isto não exclui modos físicos neutros** (revisão leve, L1). A escala global δg = 2ḡ e a translação do
escalar δφ = const são perturbações físicas em s = 0, não são de gauge e têm δω = 0, logo L não as vê. São
simetrias exatas. Os Lemas G e R não se estendem a Re s = 0.

