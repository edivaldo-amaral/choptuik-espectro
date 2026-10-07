# Revisão independente de S3a (tiles de Perron) — 24/09/2026

Enunciado sob revisão: `K(s)` (sharp, setor A, Fourier com sinal) é injetivo no L² do toro de
raios (129/128, 9/8), em `0 <= Re s <= 1,765`, `0 <= Im s <= 1/4`, fora de `|s − 0,401| < 1/8`,
pelos 154 tiles de `build/cobertura_lab/grande-*.json`.

Li só o que foi indicado (S3 §10.4–10.13, S3B_DISCO, os scripts listados, Rump 2006 pp. 1–8).
Não li a auditoria anterior nem `test_tile_perron_blocks.py`. Todos os scripts e saídas estão nesta
pasta. Nenhum arquivo fora dela foi alterado.

## 1. Veredictos

| artefato | veredicto | evidência |
|---|---|---|
| R1 forma fechada e cotas `cauda_Rinv`, `Rinv_uniforme`, `descida` | **OK** | `r1_forma_fechada.py`, `r1_saida.txt`: fórmula contra `inv` com erro 2,5·10⁻¹⁵; 336 + 24 + 8 casos, pior verdade/cota 0,9714 (convergida na truncagem) |
| R2 enumeração dos 16 blocos | **OK com LACUNA** | tabela na §2; em (Z, far) e (T1, far) falta o termo de μR_x (L1) |
| R3 ataque às cotas | **OK** (nenhuma violação) | `r3_ataque.py`, `r3_saida.txt`, `r3_resultados.json`, `r3b_*`, `r3c_*`, `r3d_*`: 7 configurações, pior razão 0,9999998; far exatos ≤ 0,27 |
| R4 Rump, hipótese por hipótese | **OK** (com duas notas menores) | §3.3; `r4_rump.py`, `r4_saida.txt`: verdade rigorosa em Arb, nenhuma cota abaixo, 0 falsos positivos em 42 tentativas |
| R5 cobertura | **OK** | `r5_cobertura.py`, `r5_saida.txt`, `check_cover_saida.txt`: 1,2·10⁵ pontos racionais, 0 descobertos; `check_cover` COMPLETA |
| R6 estrutura (i)–(iv) | **OK** | §3.5 (demonstrações, com as hipóteses usadas) |
| Enunciado: μ fixo | **LACUNA a esclarecer** | L2: o tile usa só μ = 722873400/2³², mas o símbolo foi certificado para μ ∈ [1/6, 0,17] |
| Versões mistas do código | **LACUNA (rastreabilidade)** | L3: o JSON não guarda a versão do código; conferi à parte que as versões antigas também são válidas |

**Conclusão.** Não achei nenhum erro que invalide o certificado. Achei uma falha real de argumento:
a premissa da linha 475 de `tile_perron.py` é falsa. Numericamente ela não tem efeito, e a
completei com uma cota fechada e com a reverificação dos 154 tiles. Achei também uma questão de
enunciado (μ) que o autor precisa esclarecer. A região e a cobertura estão corretas.

## 2. Tabela R2: os 16 blocos de `E = I − A H`

Notação: `Q_Z = Q0(sZ)`, `Q_1 = Q0(s1)`, `Q_2 = Q0(s2)`; `δ_Z = s − sZ`, `δ_1 = s − s1`,
`δ_2 = s − s2`. A cota usada é `N0[i][j] + r_j N_j[i][j]`, com `N_j` = `NZm`, `N1m` e `N2m`
(este último para T2 e para far). `B_t` é o fundo truncado (|d| ≤ 20, coeficientes radiais ≤ 40),
e `dB = resto_banda + erro de RT + γ₈·Young_abs` cobre `B − B_t` e o arredondamento das entradas.
`nB = 1,64967 + dB`. As linhas citadas são de `scripts/tile_perron.py`.

| (linha ← coluna) | bloco de E | cota no código (linha) | por que vale / onde entra dB |
|---|---|---|---|
| Z ← Z | `(I − V H_ZZ(sZ)) − δ_Z V P_Z Q_Z P_Z` | `rho` (247) + rZ·`VQ_ZZ` (251) | `Q_Z` preserva Z. `dH` (243) soma o produto (γ), `dB·nQ` (B−B_t em Z×Z) e `nB·dQ` (erro de R). R3: 0,9997 |
| Z ← T1 | `−V P_Z (B + δ_1) Q_1 P_T1` | `a1 = ‖V H_Z,T1‖` (463) + r1·`VQ_Z1` (464) | `HZ1` são as linhas Z de `P_G+ B_t Q_1 P_T1` (323). `dH1` (321) inclui `dB·nQ1` e `nB·dQ1`. `P_Z Q_1 P_T1 ≠ 0` só nos modos de Z (334–335). R3: 0,983 |
| Z ← T2 | `−V P_Z (B + δ_2) Q_2 P_T2` | `a2` (470, sanduíche `V S_Z V*`) + r2·`VQ_Z2` (471) | `S_Z = Σ Y_Z Y_Z*` sobre todas as colunas de T2 com \|m\| < MP + 20 (389–400). `extra` (424) tem `dB·QT2` e `nB·dQ2`. R3: 0,996 |
| **Z ← far** | `−V P_Z (B + δ_2) Q_2 P_far` | `nV(nB·descida(NZ+41, Nc) + dB·Qf)` (478) + r2·`nV·descida(NZ, Nc)` (485) | **LACUNA L1.** O código supõe `P_Z B_t P_{n ≥ NZ+41} = 0`, mas μR_x (d = 0) leva c1 ímpar n → c2 par n' < n para todo n. O termo que falta tem cota fechada ≤ 1,9·10⁻¹⁴·nV no laboratório. O δ-termo está certo: `P_Z Q_2 P_far` só desce |
| T1 ← Z | `−P_T1 B Q_Z P_Z` (δ-termo nulo: `Q_Z` preserva Z) | `b1` (288) | Linhas de T1 alcançáveis de Z por B_t: n < NZ + 41 (261). `extra` (287) tem `dB·nQ` e `nB·dQ`. R3: 0,9999 |
| T1 ← T1 | `−P_T1 (B + δ_1) Q_1 P_T1` | `alfa11` (327) + r1·`Q_11` (333) | Bloco finito exato; `dH` (321). R3: 0,983 |
| T1 ← T2 | `−P_T1 (B + δ_2) Q_2 P_T2` | `c12` (426) + r2·`Q_12` (427) | Gram acumulado `S1`, `SQ1`. `extra` com `dB·QT2`. R3: 0,981 |
| **T1 ← far** | `−P_T1 (B + δ_2) Q_2 P_far` | `nB·descida(NP+41, Nc) + dB·Qf` (479) + r2·`descida(NP, Nc)` (485) | **LACUNA L1**, igual: falta μR_x de n ≥ NP + 41 até as linhas de T1 (≤ 2,3·10⁻¹⁰ no laboratório) |
| T2 ← Z | `−P_T2 B Q_Z P_Z` | `b2` (289) | B_t leva Z a n < NZ + 41 e \|m\| < MZ + 20. `extra` = `dB·nQ + nB·dQ`. R3: 0,998 |
| T2 ← T1 | `−P_T2 B Q_1 P_T1` (δ-termo nulo: `Q_1` preserva G+) | `c21` (355) | Linhas fora de G+ alcançáveis (`linhas_fora`); `dB·nQ1 + nB·dQ1`. R3: 0,9999 |
| T2 ← T2 | `−P_T2 (B + δ_2) Q_2 P_T2` | `alfa22` (428) + r2·`Q_22` (430) | Schur sobre d: `Σ_d sup_m`, exato por modo (\|d\| ≤ D); Young para D < \|d\| ≤ 20 (`resto`); `dB·QT2`; erro de R com `‖B_d‖ ≤ ‖B‖` (média de rotações unitárias em θ). R3: 0,44 (frouxa 2,3×) |
| T2 ← far | `−P_T2 (B + δ_2) Q_2 P_far` | `beta = nB·Qf` (436) + r2·`Qf` (485) | `‖P_T2 B Q_2 P_far‖ ≤ ‖B‖ ‖Q_2 P_far‖`. `Qf = max(cauda_Rinv(Nc), Rinv_uniforme(Mc))` (432–435), já que `Q_2 P_far` é bloco-diagonal em m. nB inclui dB |
| far ← Z | `−P_far B Q_Z P_Z` (δ nulo) | `dB·nQZ` (481) | B_t leva Z para dentro de F (Nc ≥ NP + 41 ≥ NZ + 41; Mc ≥ MP + 20); só B − B_t sai de F. R3 (exato): 0,0015 |
| far ← T1 | `−P_far B Q_1 P_T1` (δ nulo) | `dB·nQ1` (481) | Idem: B_t leva G+ a n ≤ NP + 39 < Nc. R3 (exato): 0,030 |
| far ← T2 | `−P_far B Q_2 P_T2` (δ nulo: `Q_2` preserva F) | `alfa_f2` (429) | Linhas far alcançáveis pela banda (n < Nc + 41, e o modo inteiro se \|m+d\| ≥ Mc), soma sobre d, `resto`, `dB·QT2`. R3d (exato, fundo completo): 0,18–0,27 |
| far ← far | `−P_far (B + δ_2) Q_2 P_far` | `beta` + r2·`Qf` | Como em T2 ← far |

Todos os 16 pares têm cota. As 14 cotas fora de (Z, far) e (T1, far) estão corretas.

## 3. Números

### 3.1 R1: forma fechada e cotas do inverso livre

Por substituição regressiva em `(J + σ) x = e_c`, com `S_j = Σ_{k≥j} v_k x_k` (logo
`S_j = S_{j+1}(1 − v_j/d_j)` e `S_c = v_c/d_c`):

```
(J+σ)⁻¹_{jc} = 1/d_c  (j = c);   −(v_c/(d_j d_c)) Π_{j<k<c} (1 − v_k/d_k)  (j < c);   0  (j > c),
d_k = μ(n_k+1) + σ,  v_k = 2μ n_k,  |1 − v_k/d_k| = |μ(1−n_k)+σ| / |μ(1+n_k)+σ| ≤ 1  se Re σ ≥ −μ.
```

Com pesos: `sqrt(ω(n_j)/ω(n_c)) ≤ c_w κ₂^{n_j−n_c}`, onde `c_w = sqrt(1+κ₂⁻⁴)`. Também
`|v_c/d_c| < 2` e `1/|d_j| ≤ 1/(μ(n_j+1))` para `Re σ ≥ 0`.

- **`Rinv_uniforme`**: `|d_k| ≥ h = |m|/2 − |Im s|`. Pelo teste de Schur, colunas e linhas ficam
  `≤ (1 + 2c_w Σ_{k≥1} κ₂^{−k})/h = (1 + 2c_w/(κ₂−1))/h`. Confere com o código.
- **`descida(J, Nc)`**: as entradas ficam `≤ (2c_w/(μ(n_j+1))) κ₂^{n_j−n_c}`. As colunas somam
  `≤ (2c_w/μ) κ₂^{J−Nc} H_J` e as linhas `≤ (2c_w/μ) κ₂^{J−Nc} κ₂/(κ₂−1)`. Confere, inclusive o
  fator 1,0001.
- **`cauda_Rinv`**: Schur com as mesmas majorantes, somadas até `Nbig = 3000` em log. Refiz as
  caudas para colunas `> Nbig`: dividindo `n_j < n_c/2` e `n_j ≥ n_c/2`, a soma é decrescente em
  `n_c`. As linhas acima de `Nbig` também conferem.

Testes (`r1_saida.txt`). Uso a norma pesada verdadeira com colunas `n ∈ [N, 1400)` e linhas completas;
como J é triangular superior, isso é uma cota **inferior** da verdade, e qualquer violação seria real.

| cota | casos | pior verdade/cota |
|---|---:|---:|
| fórmula fechada contra `np.linalg.inv` (n ≤ 40 e ímpares ≤ 81, 15 σ, passos 1 e 2) | 210 | erro relativo 2,5·10⁻¹⁵ |
| `cauda_Rinv` (s ∈ {0; 0,25i; 0,2+0,1i; 0,401+0,125i; 1+0,25i; 1,765}, m de −100 a 198, N ∈ {24, 80, 160, 400}) | 336 | **0,9714** (m = 198, s = 0,25i; o pico está em n ≈ 590 e a verdade converge: 600 → 0,0838, 1000 a 2600 → 0,08473 contra a cota 0,08722) |
| `Rinv_uniforme` (\|m\| de 20 a 202) | 24 | 0,393 |
| `descida` ((24,77) … (380,400)) | 8 | ≤ 1,5·10⁻³ (a mais frouxa: 10³×) |

### 3.2 R3: ataque às cotas, com montagem própria

Caixa Z 8×24, G+ 12×36, F 32×77, D = 4, banda 20,40, cotas de `tile_perron` no modo rigoroso.
Montei B com o fundo completo via `ft.B_rad`; confere com `signed_operator.operador − J` até
1,8·10⁻¹⁵. `Q0` vem da minha fórmula fechada. Uso o `V` do certificado, reordenado.

Em cada bloco testei três coisas: (a) `‖E0‖` contra `N0` e `‖C‖` contra `N_j` separadamente;
(b) s em 32 direções na borda do tile; (c) `|δ_j| = r_j` em 32 direções, o que é mais forte que o
tile.

Configurações `(sZ, rZ; s1, r1; s2, r2)`:
1. `(0,1i, 0,02; 0,02+0,1i, 0,05; 0,1+0,1i, 0,15)`, com Re sZ = 0;
2. `(0,25i, 0,02; 0,02+0,23i, 0,05; 0,125+0,125i, 0,2)`, com Re sZ = 0 e Im sZ = 1/4;
3. `(0, 0,015; 0,03+0,03i, 0,06; 0,125+0,125i, 0,2)`;
4. `(0,25+0,25i, 0,03; 0,23+0,22i, 0,07; 0,2+0,2i, 0,12)`;
5. `(0,401+0,13i, 0,02; 0,42+0,15i, 0,05; 0,375+0,125i, 0,1768)`;
6. `(0,276; 0,02)` com um centro só;
7. `(1+0,2i, 0,1; 1,05+0,15i, 0,18; 1,125+0,125i, 0,25)`.

Maior razão verdade/cota por par, sobre as 7 configurações:

| par | ‖E0‖/N0 | ‖C‖/N_j | borda do tile | \|δ_j\| = r_j |
|---|---:|---:|---:|---:|
| Z ← Z | 0,00016 | 0,999999 | 0,99974 | 0,99974 |
| Z ← T1 | 0,99992 | 0,9999998 | 0,98255 | 0,98255 |
| Z ← T2 | 0,99971 | 0,999999 | 0,99577 | 0,99577 |
| T1 ← Z | 0,99990 | — (C = 0) | 0,99990 | 0,99990 |
| T1 ← T1 | 0,99988 | 0,99977 | 0,98294 | 0,98294 |
| T1 ← T2 | 0,99997 | 0,999999 | 0,98089 | 0,98089 |
| T2 ← Z | 0,99831 | — | 0,99831 | 0,99831 |
| T2 ← T1 | 0,99991 | — | 0,99991 | 0,99991 |
| T2 ← T2 | 0,471 | 0,988 | 0,437 | 0,437 |
| far ← Z (exato, caixa 50×126) | 0,0015 | — | — | — |
| far ← T1 (exato, caixa 54×138) | 0,030 | — | — | — |
| far ← T2 (exato, `r3d_far_T2.py`, 3 valores de s2) | 0,18 / 0,20 / 0,27 | — | — | — |
| Z ← far, T1 ← far (parciais: colunas far até 48×260) | ≤ 10⁻⁴ | | | |

**Nenhuma violação.** Os blocos exatos ficam a menos de 10⁻⁴ da verdade, e os δ-coeficientes batem
no fator 1 + 10⁻⁶ da cota verificada.

Os níveis far na caixa pequena são dominados pela `descida`, que nela é enorme:
`descida(65, 77) = 24`. Por isso os blocos parciais Z ← far e T1 ← far não testam a lacuna L1, que
tratei à parte (§4). Para T2 ← far e far ← far a cota é `‖B‖·‖Q0 P_far‖`. `‖Q0 P_far‖` foi testado
no R1. Para `‖B‖` (entrada externa 1,64967), a compressão `‖P_F B P_F‖` dá 1,509 (40×60),
1,549 (24×120) e 1,557 (16×200), coerente (`r3c_saida.txt`).

Young por modo (`norma_Bd_young`, d ≠ 0, usada em `resto` e em `dB`): a verdade fica entre 0,42 e
0,50 da cota para d = ±2 … 40 (`r3b_saida.txt`). A cota é válida porque, em coeficientes
simétricos com ω(0) = 1, `mult` é a multiplicação de Laurent restrita, e
`‖mult(a)‖ ≤ sup_{|z|=κ₂}|a| ≤ l1(a)`.

### 3.3 R4: o critério de Rump em `verified_norms`

Hipóteses do Corolário 2.4 (com o Teorema 2.3 e a cota I da §3), e a linha que garante cada uma:

| hipótese | onde | conferência |
|---|---|---|
| Matriz fatorada exatamente simétrica | `cota_gram` 131–135: `fl(a+b) = fl(b+a)`, `fl(a−b) = −fl(b−a)`, ×0,5 exato; a diagonal (136, 141) não quebra a simetria. Além disso, `dpotrf` (66) só lê um triângulo | OK |
| `ã_ii ≤ a_ii − c` | 141 subtrai `_margem_rump`, que inclui `2u·max a_ii` (59): `fl(a−c'') ≤ a − c'' + u·a` | OK |
| `c ≥ ‖Δ(A)‖` pela cota I | 58: `γ_{n+1}/(1−γ_{n+1})·tr + n·3(2n+max a)·η`, que é exatamente a cota I (M = 3(2n + max a_νν)); o fator 1 + 10⁻⁶ cobre a soma | OK |
| dimensão em `γ_{n+1}` | 55–56: `n = len(d)` = dimensão da realificada (2n complexo) | OK |
| hermitianas complexas | realificação (130–133), PD se e só se a hermitiana for. O Teorema 2.3 já cobre hermitianas; a realificação só é mais conservadora | OK |
| erro de formar `X*X` | `cota_norma` 156: `γ_{4(k+1)}·‖Y‖_F² ≥ √2 γ_{k+2}·‖ \|Y\|*\|Y\| ‖`; entra em `extra` (138). A parte hermitiana de um erro não aumenta a norma | OK |
| arredondamento ao montar `t²I − G` | 138: `u(t² + max\|g_ii\|) + u‖G‖_F` | OK (nota N1) |
| arredondamento ao mais próximo, uma precisão, qualquer ordem | LAPACK/BLAS em binary64. Somas em árvore (blocos do `dpotrf`) satisfazem (2.4): cada termo leva no máximo k−1 fatores, e a razão leva no máximo k−2 | OK (nota N2) |

- **Nota N1.** A diagonal é `fl(fl(t·t) − g_ii)`. O que se prova é `‖X‖² < fl(t²)`, ou seja
  `‖X‖ < t(1 + u/2)`, mas a função devolve `t`. A diferença de 5,6·10⁻¹⁷ relativo é absorvida por
  `FATOR_PESO` e `ARRED`.
- **Nota N2.** Com FMA o produto sai exato, o que cabe no modelo com ε_i = 0. No subnormal, o
  arredondamento da soma fundida (≤ η/2) fica fora do modelo literal de Rump; o mesmo vale se o
  ambiente ligar FTZ/DAZ. Isso é irrelevante nas escalas usadas: `r4_saida.txt`, teste E, escalas
  10⁻²⁸⁰ a 10¹⁴⁰.

Testes (`r4_saida.txt`, 115 s). A verdade é rigorosa em Arb (python-flint, 200 bits; o mpmath não
está no .venv): `λmax ≤ (λ̃max + ‖F‖_F)/(1 − ‖D‖_F)`, com `F = U*GU − Λ` e `D = U*U − I` calculados
em bolas.

- **A/B.** Grams não hermitianos (um triângulo perturbado em ~10⁻¹⁵ relativo, com `erro` rigoroso),
  n = 10–400, nos casos quase singular (λmin/λmax = 10⁻¹⁰), topo aglomerado (metade dos autovalores
  iguais) e topo a 10⁻¹³: cota/verdade − 1 = 1,0·10⁻⁶ em todos. Nenhuma cota ficou abaixo.
- **C.** Falso positivo: o primeiro t com `t² = inf(λmax)(1 − ε)`, para ε de 10⁻³ a 0, n = 20, 100
  e 300. Em 42 tentativas o teste nunca aceitou.
- **D.** `cota_norma(X)` para 60×40, 200×120, 400×300, 300×400 e 150×150 de posto 1: todas acima da
  verdade (+10⁻⁶).
- **E.** Escalas 10⁻²⁸⁰, 10⁻¹⁵⁰ e 10¹⁴⁰: OK.

### 3.4 R5: cobertura

`r5_saida.txt`:

- São 156 tiles gravados, 154 com `ok`. Os 154 passam em `check_cover.reverifica` e na minha
  reimplementação, que também exige `Re` dos centros ≥ 0, sem nenhuma divergência. θ vai de 0,6127
  a **0,996115**; o pior tile é `sZ = 0,015625 + 0,234375i`. Todos os arquivos têm
  `rigoroso = True` e as caixas 20×80, 40×160, 200×400.
- **10⁵ pontos uniformes** (racionais de denominador 2⁴⁰): 0 descobertos.
- **9715 pontos** nas quatro bordas do retângulo e sobre a circunferência de D, estes racionais
  **exatamente** sobre o círculo pela parametrização racional, mais `0,401 ± 1/8` e `0,401 + i/8`:
  0 descobertos.
- **10⁴ pontos** nos 211 vértices das malhas de tiles, médios e grandes (o próprio vértice e
  deslocamentos ≤ 10⁻⁷ hZ): 0 descobertos.
- `check_cover.py build/cobertura_lab --disco 0.401,0.125` dá COBERTURA COMPLETA (1770 células).

Folga uniforme: todos os tiles continuam certificados com N multiplicado por até
`1/0,996115 = 1,0039`.

**Conjugação (§10.6).** Tomo `(Cx)_m = conj(x_{−m})`. Vale `C K(s) C = K(s̄)`, porque J é real,
`conj(σ(m)) = σ(−m)` em `s̄`, `modo(g, j, −d) = conj(modo(g, j, d))`, `combos` e `refl` são reais e
μR_x é real. `C` é isometria de Y (pesos `κ₁^{2m} + κ₁^{−2m}`), e Y ⊂ X com sinal, com
`Dom_Y = Q0(Y) ⊂ Q0(X)`. Então a injetividade em X para `Im s ∈ [0, 1/4]` dá a injetividade em Y
para `|Im s| ≤ 1/4`. **Não** dá a injetividade no espaço com sinal para `Im s < 0`; a frase
"logo também em −1/4 ≤ Im s ≤ 0" da §10.12 vale só para Y. A região e o disco aberto estão tratados
corretamente: a circunferência de D é coberta e as bordas são fechadas.

### 3.5 R6: estrutura

1. **(i) `K_livre` é bloco-triangular superior em [Z, T1, T2, far].** No modo m,
   `K_livre = J_rad + im/2` é triangular superior em n (a coluna c só tem entradas em j ≤ c), logo
   preserva `span{n < N}`. Z, G+ = Z ∪ T1 e F = G+ ∪ T2 são caixas em (m, n), logo invariantes.
   Daí as entradas (T1, Z), (T2, Z ∪ T1) e (far, F) de `K_livre` e de `Q0` são nulas.
2. **(ii) P é bijeção de X sobre `Dom K`.** Pela identidade do resolvente,
   `P = Q0(s2)·T`, com `T = I + (s2−sZ) Q0(sZ) P_Z + (s2−s1) Q0(s1) P_T1`. `T − I` tem posto finito
   e é bloco-triangular. Os blocos diagonais `P_Z (K_livre + s2) Q0(sZ) P_Z` e o análogo em T1 são
   invertíveis, porque a diagonal `μ(n+1) + s2 + im/2 ≠ 0` para `Re s2 ≥ 0`. Logo T é invertível e
   `P(X) = Q0(s2)X`. Esse é o **domínio maximal**: `K_livre + σ` é injetivo em X, porque de
   `T_j = T_{j+1}(1 − v_j/d_j)` e `Π_{k≥j}(1 − v_k/d_k) → 0` segue `T_j = 0`. A §10.6 diz
   "bijeção" sem esse argumento. Hipóteses usadas: `Re s_ℓ ≥ 0`, que `tile_perron` exige e eu
   conferi no JSON, e `Q0` limitado (R1).
3. **(iii) `K(s) P = I + Σ_ℓ (B + s − s_ℓ) Q0(s_ℓ) P_ℓ`.** Segue de
   `K(s) = (K_livre + s_ℓ) + (s − s_ℓ) + B` e de `Σ P_ℓ = I`, com o nível T2 + far no centro s2.
4. **(iv) Lema de Perron.** Suponha `‖P_i E P_j‖ ≤ N_ij`, `v > 0` e `N v ≤ θ v` com `θ < 1`. Se
   `Hx = 0`, então `Ex = x`. Tome `μ = max_i ‖P_i x‖/v_i`; então
   `‖P_i x‖ ≤ Σ_j N_ij v_j μ ≤ θ v_i μ`, logo μ = 0 e x = 0. Daí H é injetivo, e K(s) também,
   porque P cobre o domínio.
   - Usa: a decomposição ortogonal nos pesos diagonais, `N_ij` como **cotas superiores**, v > 0 e
     θ < 1.
   - **Não** usa V invertível, irredutibilidade nem o cálculo de ρ(N). O v é exibido e conferido
     em racionais. A frase da §10.4 "ρ(N) < 1 ⇒ ∃ v > 0 com Nv < v" é verdadeira (por
     perturbação), mas é dispensável.
   - Para todo s do tile, `|δ_j| ≤ r_j` pelas contenções, que `check_cover` confere em racionais.

## 4. Erros e lacunas, por gravidade

Nenhum ERRO afeta o resultado. As lacunas, da mais grave para a menos grave:

1. **L1. Premissa falsa nos blocos (Z, far) e (T1, far) (`tile_perron.py` 472–479).**
   - **O problema.** A cota `nV(nB·descida(NZ+41, Nc) + dB·Qf)` supõe que B_t não desce mais que
     BAND_N (comentário da linha 475). Mas `B_rad(d = 0)` contém μR_x: c1 ímpar n → c2 par n' < n,
     para **todo** n' (`fourier_tail_bound.py` 93–99). Exemplo explícito
     (`r3b_muRx_young.py` (a)): a entrada `B_t[(c2, 0), (c1, 401)] = 2μ = 0,336614`. Falta
     `T = P_Z μR_x P_{n ≥ NZ+41} Q0(s2) P_far`, e o análogo em T1.
   - **O tamanho.** Pela forma fechada, a entrada do produto é
     `≤ 4c_w² κ₂^{n'−n_c}(H(n_c+1)+1)`, e por Schur `‖T‖ ≤ 1,9·10⁻¹⁴` (Z) e `2,3·10⁻¹⁰` (T1) no
     laboratório (80/160/400). Medido: 1,1·10⁻¹⁶ e 1,0·10⁻¹².
   - **O conserto.** Somei `2·nV·1,9·10⁻¹⁴` a `N0[0][3]` e `2·2,3·10⁻¹⁰` a `N0[1][3]` e refiz o
     Perron dos 154 tiles em racionais com o v gravado: **todos passam, pior θ = 0,9961146191**
     (`r3b` (d)).
   - **Efeito numérico nulo.** O termo também é menor que a folga da própria `descida`
     (`nB·descida` = 9,3·10⁻¹³ e 1,2·10⁻⁸). Mas o argumento, como escrito, não é prova. Na caixa
     pequena o termo mede 1,7·10⁻³, e só não viola porque lá `descida` vale 24.
   - **Correção sugerida.** Somar esse termo fechado em `montar` e corrigir o comentário.
2. **L2. μ fixo no certificado dos tiles (a esclarecer pelo autor).**
   - **O problema.** Todos os tiles usam `μ = 722873400/2³²` ("MU_REFA, RT seção 5") em `J_rad` e
     em μR_x. Nenhum termo cobre `δμ`. Já o símbolo (§9.2) foi certificado para
     μ ∈ [1/6, 17/100] "pelos extremos", e `EPS_FUNDO` cobre só os campos.
   - **Quando pesa.** Se o μ do problema é por definição o valor diádico de RefA, não há lacuna. Se
     é conhecido só com erro, a perturbação `δμ(J₁ + R_x)` é ilimitada em `J_rad`. Ela é
     relativamente limitada (`J₁ Q0 = (I − σ' Q0)/μ`, de norma ≲ 130), mas precisa entrar
     explicitamente. Pela folga uniforme de 0,39% e `‖V‖ ≤ 25,5`, estimo que δμ ≲ 10⁻⁶ caberia.
     Isso é uma estimativa, não uma prova.
3. **L3. Rastreabilidade: tiles de versões mistas sem carimbo de versão.**
   - **O problema.** O JSON não registra o hash do código, e a §10.12 admite versões anteriores e
     posteriores a 35b118b e 9b40673.
   - **Conferi as diferenças:**
     - antes de 9b40673 usava-se `young[d]·dq`: vale, porque Young é cota para d ≠ 0 e
       `young[0] = 2,558 ≥ nB`;
     - o `cota_gram` antigo realificava sem simetrizar, mas o `dpotrf` lê um só triângulo. A
       matriz efetiva difere da exata por ≤ 2‖Δ‖_F em vez de ‖Δ‖₂, o que dá um defeito relativo
       ≤ ~10⁻⁷ nas cotas antigas;
     - o `(1+8u)` sobre c é coberto pela folga da cota I.
   - **Conclusão.** Tudo isso é coberto pela folga uniforme de 0,39%. A matemática fecha, mas a
     reprodução depende do código do laboratório, que não pude rodar (caixas grandes proibidas
     aqui).
4. **Entradas externas não auditadas aqui.** `‖B♯‖ ≤ 1,64967` (símbolo) e `EPS_FUNDO = 40·2⁻²⁵`
   (erro de RT). Há só uma checagem de coerência: `‖P_F B P_F‖ ≤ 1,557`.
5. **Notas menores.**
   - N1: `cota_gram` devolve t, e a cota provada é t(1 + u/2).
   - N2: FMA e FTZ/DAZ no modelo de Rump (irrelevante).
   - A §10.12 estende "também em −1/4 ≤ Im s ≤ 0" sem dizer que isso vale só no espaço simétrico Y.
   - A §10.6 afirma a bijeção de P sem o argumento do domínio maximal (dado aqui, §3.5).
6. **Frouxidões, que não são erros:** `descida` (10³×), `Rinv_uniforme` (2,5×), `alfa22` (2,3×),
   `alfa_f2` (4–5×) e Young (2×).

Resposta ao pedido: com L1 completada como acima e L2 esclarecida (μ exato por definição), os 154
tiles provam o enunciado. Não achei contraexemplo a nenhuma cota.
