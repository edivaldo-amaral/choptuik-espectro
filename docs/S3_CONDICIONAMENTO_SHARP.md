# S3: o sistema de constraints é 130× mais bem condicionado que L

23/09/2026. Diagnóstico em ponto flutuante, **não é prova**. Script:
`scripts/measure_sharp_vs_selected_V.py`.

## 1. A medida

Mesma grandeza nos dois operadores, cada um nos pesos do espaço em que o seu
certificado teria de valer:

```
||V(s)|| = || J (op + s I)^-1 ||_w       (norma de coluna l1 ponderada)
```

`L`: operador selecionado, 5 componentes, pesos RT. `K`: operador de propagação
das constraints (sistema sharp), 2 componentes, pesos `(129/128, 9/8)` — os que
`SHARP_DOMAIN.md` exige, por causa da perda de raio de `C(s)`.

Máximo sobre o contorno `Gamma_A = [1/8, 1] × [-1/4, 1/4]` (12 pontos por lado,
metade superior, por simetria de conjugação):

| malha | máx `||V_L||` | máx `||V_K||` | razão | `||V_L(1)||` | `||V_K(1)||` |
|---|---:|---:|---:|---:|---:|
| 4×12 | 201,5 | 10,3 | 20× | 79,1 | 2,9 |
| 6×18 | 3822,1 | 12,6 | 303× | 143,6 | 3,8 |
| 8×24 | 1636,2 | 13,5 | 121× | 183,5 | 4,3 |
| 12×36 | 1846,2 | 14,0 | 132× | 214,4 | 4,5 |

**Validação.** `||V_L(1)|| = 183,5` em 8×24 reproduz os 183,49 registrados; o
máximo 1636,2 reproduz a tabela de deflação de `ALTERNATIVE_ROUTES.md` §5.22; e
1846,2 em `s = 1/8` no 12×36 reproduz a varredura de campo de valores.

`||V_K||` converge (incrementos 2,3 → 0,9 → 0,5). `||V_L||` oscila conforme a raiz
de gauge encosta na borda, e mesmo longe de qualquer raiz, em `s = 1`, ainda
cresce com a malha. **Essa diferença de comportamento é ela própria um dado.**

Perto do eixo imaginário, `||V_K||` continua pequeno:

| `Re s` | máx `||V_K||` em `Im ∈ [0, 1/4]`, 12×36 | em `Im = 0` |
|---|---:|---:|
| 0 | 20,6 | 16,9 |
| 0,05 | 17,4 | 14,4 |
| 1/8 | 14,0 | 12,2 |

`K(0)` é invertível — onde `L` tem o modo de fase. É consistente: a translação de
fase satisfaz as constraints.

**Raízes de `K` em `Re s > 0`, `|Im| ≤ 1/4`:** uma só, em todas as malhas —
0,4013 → 0,4011 → 0,4010 → 0,4010.

## 2. O que isso faz com S3

Pela identidade `K(s)C(s) = M(s)L(s)` (`SHARP_LINEARIZED_IDENTITY.md`, auditada
como universal e exata), se `K(s0)` é injetivo e `L(s0)h = 0`, então
`C(s0)h = 0`. Com o argumento de jatos de `SHARP_LOCAL_EXCLUSION.md`, vale para
cadeias inteiras. Então S3 se divide em duas peças de natureza diferente:

- **S3a — exclusão de `K`:** `K` injetivo em `{Re s > 0} ∩ {|Im| ≤ 1/4}` fora de um
  disco `D` em torno de 0,401. Consequência: **todo** modo de `L` fora de `D` —
  onde quer que esteja — satisfaz as constraints, com cadeias. Não é preciso
  saber onde estão as raízes de `L`.
- **S3b — o disco `D`:** `L` tem ali uma raiz algebricamente simples cujo
  autovetor viola as constraints (`C(s)h ≠ 0`). É o testemunho já desenhado em
  `CONSTRAINT_WITNESS_NORMS.md` e na auditoria matemática, com
  `c_ref ≈ 0,011` na rota (a).

S3a é uma certificação **só de `K`**. É aí que o condicionamento importa.

## 3. Uma constante 104× frouxa na cadeia sharp existente

`sharp_guarded_bridge.py`, linha 51, usa

```
kappa = ||J|| · ||K^-1|| = 103,69 × 11,588 = 1201,5
```

como cota de `||H_FF^-1|| = ||J K^-1||`. O valor verdadeiro em `s = 0,733`, 12×36,
é **8,41** (e `||K^-1|| = 8,507` reproduz exatamente o `inverse_candidate_norm`
do relatório `sharp-disk-12x36-background.json`). O caminho direto já está no
código — a linha 57 toma `min(kappa, direct)` se o relatório trouxer
`preconditioned_inverse_bound` — mas o relatório do 12×36 não o traz. A
auditoria de 16/09 chegou a rodar `independent_sharp_disk.py --mode right`, que
calcula a quantidade certa.

**Nesta ponte isso quase não pesa:** `kappa` entra no defeito de Schur só por
`c·kappa = 1,82e-5 × 1201,5 = 0,022`. O corte `G = 160×960` é fixado pela cauda
sozinha, `d = (beta♯ + |s|)·t♯(G) < 1`, com `beta♯ ≤ 4,622`. A correção passaria
`delta` de 0,935 para ≈0,915 — marginal. Mas a mesma confusão entre
`||J K^-1||` e `||J||·||K^-1||` custaria duas ordens de grandeza em qualquer
dimensionamento no estilo de `TT_CEILING.md`, que é linear em `||V||`.

## 4. Onde S3a fica em custo

O que está **provado** (`SHARP_INTERMEDIATE_SHELL.md`, condicionado à bola RT):
no disco `|s − 0,733| ≤ 1/32`, a compressão a `F + T` é invertível com
`F = 12×36` e cauda a partir de `G = 160×960`. O que **falta** é o Schur finito
da casca intermediária, `E_U = H_UU − H_UZ H_ZZ^-1 H_ZU`, com 228.366 DOFs.

Comparação grosseira de tamanho, não de custo: a caixa sharp `160×960` tem
≈2,3·10⁵ DOFs; a caixa `107×384` de `L` que o orçamento original de HPC cobria
tem ≈1,35·10⁵. **S3a vive na escala de ~1,7× a caixa original**, contra os 9794×
que a contagem de `L` exige. Ainda assim, é preciso cobrir a região inteira
(contorno ou tiles), não só um disco: o número de tiles multiplica o custo, a
menos que a parametriz da casca seja reaproveitada entre eles — `H(s)` é afim em
`s`, então há esperança, mas nada foi verificado.

Falta também a exclusão para `Re s` grande: um argumento de campo de valores
para `K` análogo ao de F3. A parte livre de `K` é `diag(K_mu, J_mu) + sI`, e
`R_J` já está certificado para `J_mu`. Parece barato; não foi feito.

## 5. Ressalvas

- Ponto flutuante, truncado, fundo RefA. O contorno foi amostrado em 12 pontos
  por lado.
- `L` e `K` são medidos em normas diferentes. É a comparação certa — cada rota
  paga na sua própria norma —, mas não é a mesma quantidade.
- A truncação de `K` tem raízes com `Re s = +0,0267`, em `Im = ±2,48` (8×24) e
  `±4,48` (12×36). Andam exatamente `ΔM/2` com a malha: artefato da borda de
  Fourier, fora da faixa `|Im| ≤ 1/4`. O argumento infinito precisa excluí-las
  explicitamente, não por inspeção.
- "1,7× a caixa original" compara DOFs, não custo. O custo da casca `E_U` não foi
  derivado.
- **Nada aqui toca a contagem global de `L`**, que T2 continua exigindo e que
  continua a 10⁴× (`RELATORIO_SETEMBRO_2026.md` §3). S3/S4 decidem *quais* modos
  são físicos; não dizem *quantos* modos `L` tem.

## 6. Por que `L` é pior que `K` (23/09)

Script: `scripts/diagnose_L_vs_K_conditioning.py`. Mesma grandeza,
`||J(op+s)^-1||`, nas mesmas normas para os dois operadores, 12×36:

| | norma | `s = 1` | `s = 1/8` |
|---|---|---:|---:|
| `L` | forte `(65/64, 5/4)` | 214,4 | 1846 |
| `L` | fraca `(129/128, 9/8)` | 89,3 | 1686 |
| `L` | ℓ² | 19,3 | **3186** |
| `K` | forte | 9,0 | 17,5 |
| `K` | fraca | 4,5 | 12,2 |
| `K` | ℓ² | 2,66 | **13,2** |

Os pesos de componente `η` **não** entram: a norma própria de `L` é a forte pura.

São dois fenômenos distintos.

**(a) Na borda do contorno, é o operador.** Em `s = 1/8`, `L` é enorme em todas
as normas. É a raiz de gauge em `mu`, a 0,0433 da borda, com projeção espectral
muito não-normal (`||P|| ≈ 80`, §5.22 de `ALTERNATIVE_ROUTES.md`). `K` não
contém esse modo: modos de gauge satisfazem as constraints.

**(b) Longe das raízes, é o regime pré-assintótico.** Em `s = 1`, `L` é 7× pior
mesmo em ℓ², e piora com a malha:

- `sigma_min((A+s)J^-1)`: 0,080 → 0,052 (`K`: 0,48 → 0,38);
- a direção quase-nula vive em modos de Fourier numa **fração fixa** da truncação
  (`m_médio ≈ 0,67 M` nas duas malhas) e leva a saída para `n` radial alto;
- a pior coluna de `V_L` é sempre `(comp 1, m=0, n=0)`, e a norma cresce porque a
  massa se espalha pelos DOFs novos — cada um contribuindo pouco —, até em ℓ¹
  sem peso (22 → 33 → 62).

O fundo só fica pequeno diante da parte livre quando `mu n > beta` e `m/2 > beta`:

| | `beta` | limiar radial `beta/mu` | limiar de Fourier `2 beta` | 12×36 (`N=35`, `M=11`) |
|---|---:|---:|---:|---|
| `L` | 10,38 | **62** | **21** | fora do regime assintótico |
| `K` | 4,62 | **27,5** | **9,2** | dentro |

**Esses são exatamente os `n > 61` e `|m| > 21` que o certificado de `L` exige**
(`ALTERNATIVE_ROUTES.md` §5.7, linha do "bound provado (10,382)"): o certificado precisa que a parte
livre domine o fundo modo a modo. A verdade converge em `n0 = 16` porque o
acoplamento que o espectro sente é menor que a cota `beta`. **É isso que o vão
entre 16 e 61 mede.**

**Previsão registrada:** `||V_L(1)||` só satura em malhas ≳ 22×62.

**Consequência:** todo certificado de `L` — a contagem global, S3b, C1 — está
preso ao mesmo limiar `beta/mu`. A alavanca é o acoplamento efetivo, o que
reconecta ao item 1 de `PROBLEMAS_ABERTOS.md` (similaridade não-diagonal).

## 7. Quanto custa S3b (23/09)

Script: `scripts/measure_s3b_cost.py`. S3b é um certificado de `L` num disco em
volta de 0,401: raiz algebricamente simples, com autovetor enclausurado para o
testemunho.

| 12×36 | grandeza | valor | vs. `Gamma_A` (1846) |
|---|---|---:|---:|
| winding num círculo, `r = 0,10–0,15` | máx `||V_L||` | ≈ 950 | 1,9× |
| **Newton–Kantorovich** no sistema aumentado | `||DF^-1||` | **247** | **7,5×** |

Newton–Kantorovich em `F(x,s) = (H(s)x, l(x)-1)` é a rota certa: `DF`
invertível equivale à raiz ser algebricamente simples (Keldysh), então o mesmo
certificado dá simplicidade **e** o enclausuramento que o testemunho pede.

Mas `||DF^-1||` cresce com a malha (150 → 210 → 247) acompanhando `||V_L(1)||`
(143 → 183 → 214), na razão ≈1,1. **S3b fica no piso pré-assintótico de `L`.**
Pelo dimensionamento linear em `||V||`, na ordem de 10³× a caixa original: 10×
mais barato que a contagem global, e ainda fora de alcance.

## 8. `K` para `Re s` grande: estimativa e plano (23/09)

Em ponto flutuante, produto interno `D = diag(w)` com pesos `(129/128, 9/8)`:

| malha | `máx Re W(-K)` | parte livre | `||fundo||_2` |
|---|---:|---:|---:|
| 8×24 | 1,007 | 0,026 | 1,36 |
| 12×36 | **1,037** | 0,026 | 1,45 |

Parte livre por componente: `c1` (só `n` ímpar, acoplamento de mesma paridade)
−0,035; `c2` (o `J_mu` de F3) +0,026. Fundo por bloco (12×36): `c1←c1` 0,20,
`c1←c2` 0,85, `c2←c1` 1,17 (inclui `mu R_x`), `c2←c2` 0,88.

Uma coincidência útil: `D = diag(w)` com `(9/8)^n` é a convenção de `R_J` com
`kappa2 = (9/8)^2 = 81/64 ≈ 1,266`. Como `f(kappa2)` se anula em 1,3205, isso
prevê exatamente a parte livre pequena e positiva que se mede.

**Feito em 23/09: `R_K <= 5,06`, certificado exato** (`certify_numerical_range_l2.py`;
verdadeiro 1,006 na truncação). `K(s)` é injetivo para `Re s > 5,06`. O que resta
de S3a é a região `0 <= Re s <= 5,06` da faixa, fora do disco em 0,401.

**Plano original, para registro:** a extensão usa a mesma máquina de F3, e essa máquina tem os
dois passos injustificados de `ALTERNATIVE_ROUTES.md` §5.25. O certo é fazê-la
já numa convenção única (`D = diag(w)`), com a estrutura de blocos em ℓ² feita
como norma espectral da matriz de normas de blocos, e não como máximo. Pela
estimativa, `R♯` certificado deve cair em 2–3.

## 9. `R♯ <= 1,765`, pelo símbolo pontual (23–24/09)

### 9.1 O método

Em vez de majorar o fundo bloco a bloco, troca-se o espaço de Hilbert pelo **L² do
toro de raios `(κ1, κ2)`**, com coeficientes simétricos e índices com sinal:
`||x||² = Σ_{m,n∈ℤ} κ1^{2m} κ2^{2n} |x_mn|²`. Nele todos os operadores do fundo
ficam **pontuais**:

- multiplicação por um campo → multiplicação por `f(w, z)`, `|w| = κ1`, `|z| = κ2`;
- a reflexão `P` → `z ↦ −z`, isometria;
- a divisão por ξ → multiplicação por `1/ξ(z)`, `ξ = (z + 1/z)/2`.

Então `Re<−Bx, x>` é a integral de uma forma Hermitiana 3×3 nos valores
`(x1(z), x2(z), x2(−z))` (com `x1(−z) = −x1(z)`), e

```
max Re W(−B) <= sup_toro λ_max(N^-1/2 Herm(−S) N^-1/2),   ||B|| <= sup_toro ||M N^-1/2||,
```

com `N = diag(2,1,1)`. Isso captura **todos** os cancelamentos entre blocos e entre
fases de Fourier, sem parte de Hankel e sem desigualdade triangular. A inclusão
`l^1(w) ⊂ L²` continua com norma ≤ 1.

### 9.2 O certificado

`scripts/certify_symbol_numerical_range.py`, grade 1024×2048 sobre
`φ ∈ [0, π) × θ ∈ [0, 2π)`: **4.194.304 centros verificados em Arb** (python-flint),
nenhuma falha. O ponto flutuante só propõe o limiar; a positividade de `tI − H` é
verificada pelos menores principais em bolas. Folga por caixa por Weyl, com
Lipschitz exatos dos campos e cota local de `|d(1/ξ)/dφ|`: 0,054. `μ ∈ [1/6, 17/100]`
pelos extremos, porque `λ_max` e a norma são convexos em `μ`.

| | cota rigorosa | verdadeiro (truncação 12×36) |
|---|---:|---:|
| `max Re W(−B)` | **1,4396** | 1,29 |
| `||B♯||` | **1,6497** | 1,43 |
| parte livre na norma do toro (`c2`) | **0,325** | 0,322 |
| **`R♯`** | **1,765** | 1,03 |

**`K(s)` não tem espectro em `Re s > 1,765`** — contra 5,06 da cota por blocos.

Protótipo em ponto flutuante sem `μR_x` (só `Γ2♯`): `max Re W ≈ 0,99`,
`||Γ2♯|| ≈ 1,49`. **Tirar `μR_x` da norma quase não muda `||B||`**: o fundo é
dominado por `Γ2♯`, não pela divisão por ξ.

### 9.3 O que isso faz com a região finita de S3a — e uma correção

Eu tinha priorizado apertar `R♯` dizendo que isso encolheria a região finita em
~3×. **A premissa era fraca**: o custo da região finita não escala com a largura
dela. Para `Re s` grande, `K(s)` é quase o livre e os tiles saem baratos. O gargalo
é a caixa verificada perto de `Re s` pequeno, e ela é fixada pelo critério de
cauda, não por `R♯`.

**Restrição estrutural.** O bloco finito é invertido com `||H_FF^-1|| ≈ 15`, e esse
fator amplifica qualquer acoplamento de ordem 1 entre `F` e o resto. O fundo acopla
modos a distâncias de até 40 em `m` e 100 em `n`. Por isso a ponte existente usa
uma **faixa de guarda** entre `F` e a cauda (zerando `F → T` por suporte), e a
própria faixa precisa ser invertida. As constantes melhores **encolhem** a caixa
verificada; não a eliminam.

Caudas do inverso livre na norma do toro (`c2`, o pior; `c1` é ~2× melhor):

| | radial `N = 100 / 140 / 180` | Fourier `M = 20 / 30 / 40` |
|---|---|---|
| `s = 0` | 0,74 / 0,56 / 0,45 | 0,75 / 0,52 / 0,40 |
| `s = 0,4` | 0,57 / 0,46 / 0,38 | 0,58 / 0,43 / 0,35 |

O fundo decai depressa nos pesos sharp: truncá-lo em `m ≤ 20, n ≤ 40` deixa só
`5·10⁻⁵` em norma ℓ¹. A guarda pode encolher de (41, 101) para (~21, ~41).

| cota do fundo | caixa verificada por tile | DOFs |
|---|---|---:|
| ℓ¹ (antes) | 160×960 | ~228 mil |
| **símbolo, com guarda reduzida** | **~46–50 × 180–220** | **~25–33 mil** |

**Ressalvas.** Estimativa de ordem de grandeza, a partir de normas de cauda em ponto
flutuante. Uma matriz densa de 30 mil DOFs são ~7 GB em dupla precisão: não cabe na
máquina local (7 GB), cabe no Kaggle (30 GB). E o número de tiles ao longo do
contorno ainda não foi estimado.

## 10. O operador na memória e o tile-protótipo (24/09)

**Montador.** `scripts/sharp_operator_builder.py` monta `K = K_livre + μR_x + Γ2♯` numa
caixa qualquer direto dos campos de RefA, sem o exportador de texto. Contra
`rt-sharp-12x36.dat`: **diferença máxima 1,8·10⁻¹⁵**. Conferem de uma vez as
convenções de `R_x`, da tabela de `Γ2♯`, das paridades e da simetria conjugada em
`m`.

**Correção de dimensionamento.** O setor A só usa `m` par (os campos `v1, v2, v3`
têm só `m` par; `v4` só ímpar e não entra em `K`). Uma caixa `M × N` tem
`(M−1)·1,5N` DOFs reais — confere com os 594 do 12×36. As estimativas da §9.3
contavam `m` ímpar e estavam **2× acima**.

**Experimento de tile** (`scripts/tile_experiment.py`, ponto flutuante, **não é
prova**). `H(s0) = I + B·Q0`, `Q0 = (K_livre + s0)^-1`, norma do toro; bloco finito
`Z` com inverso `V`, e uma faixa finita além dele fazendo o papel da cauda. O
inverso aproximado `diag(V, I)` fecha se a norma espectral de
`[[ρ, ||V H_ZT||], [||H_TZ||, α]]` for < 1.

| `s0` | bloco 12×36 | bloco 20×80 |
|---|---:|---:|
| 1,765 | 0,48 | |
| 1,0 | 0,65 | |
| 0,6 | 0,79 | |
| 0,501 / 0,401+0,1i / 0,301 (círculo `r = 0,1`) | 0,84 / 0,89 / 0,95 | |
| 0,2 | 1,02 | **0,65** |
| 0,125 + 0,25i (canto) | 1,09 | **0,68** |
| 0 / 0,25i (borda esquerda) | 1,23 / 1,23 | **0,72 / 0,72** |

As normas **reais** dos blocos (`||V H_ZT|| ≈ 0,48`, `||H_TZ|| ≈ 0,22`) ficam muito
abaixo dos produtos de normas do dimensionamento. Com bloco 20×80, **todos os
pontos difíceis testados fecham**. `||V||` chega a 10, mas o inverso não amplifica o
acoplamento com a cauda.

**Como fica o certificado rigoroso.** Cauda dividida em faixa próxima (calculada
exatamente) e cauda distante, além de uma guarda de mais de 41 modos radiais.
Assim os acoplamentos distante↔bloco ficam geometricamente pequenos: o fundo
truncado não alcança, e o inverso livre só desce com peso `(9/8)^-Δn`. A cauda
distante é limitada pela forma fechada do inverso livre,
`(R^-1)_jn = −(1/d_j)(v_n/d_n) Π_{k=j+1}^{n−1}(1 − v_k/d_k)` (conferida a 10⁻¹³),
com teste de Schur a 10–25% da norma verdadeira. Isso pede `N_far ≈ 250–300`:
~18 mil DOFs complexos na faixa próxima, ~5 GB por matriz. Cabe no PC do
laboratório (15 GB) se o código não montar matrizes cheias desnecessárias.

### 10.1 Onde o tile rigoroso emperra (24/09)

O experimento em ponto flutuante mostra que o operador verdadeiro é bem
comportado: com bloco 20×80, `α = ||H_TT − I|| ≈ 0,5`. O rigor exige limitar a
cauda **infinita**, e a cota natural por produto é `||B||·||Q0 P_T||`. Com
`||B|| <= 1,65` e `||Q0 P_T|| ≈ 0,89` em `N = 80`, dá **1,4** — três vezes a verdade.
O produto de normas perde a estrutura: a direção que maximiza `Q0` não é a que
maximiza `B`.

Três tentativas de contornar isso, medidas:

1. **Precondicionar à esquerda** (`Q0·K`) para tirar a descida da cauda para os
   modos baixos. Não adianta: `||K_TT^-1|| ≈ ||Q0 P_T||` (0,83 contra 0,89). A
   própria cauda do inverso livre tem norma ~70/N.
2. **Coercividade radial.** `λ_min(Herm K_livre | n >= N) ≈ 0,065 μN` (`c2`) e
   `0,13 μN` (`c1`) supera a borda do fundo `R_B = 1,44` a partir de `N ≈ 125` e
   `N ≈ 65`. Isso **funciona na direção radial**: a cauda radial fica invertível sem
   produto de normas.
3. **Coercividade em Fourier.** Não existe: `∂τ = im/2` é antissimétrico e some da
   parte Hermitiana. Na cauda de Fourier só resta o produto, com cota rigorosa
   ~16/M·1,2, o que pede `M ≈ 60`.

**Estimativa do tile rigoroso:** bloco exato até `N ≈ 200` e `M ≈ 60`, mais faixas
de guarda. Isso dá ~20–40 mil DOFs **complexos** (5–25 GB): Kaggle, não o PC do
laboratório. É engenharia numérica séria — verificação em ponto flutuante com
cotas de erro de Rump num sistema desse porte, mais a montagem das caudas — e
ainda antes de contar os tiles.

**A pergunta que mudaria a escala:** uma cota rigorosa para a cauda de Fourier
melhor que o produto. Por exemplo, o símbolo de `B` combinado com o de `Q0` modo a
modo, já que `Q0` é diagonal em `m`.

### 10.2 A cauda de Fourier, modo a modo: ganho de 2× (24/09)

`scripts/measure_fourier_tail_modewise.py` (ponto flutuante, não é prova).

`Q0` é diagonal em `m`. O fundo se decompõe em modos de Fourier `B = Σ_d S_d B_d`,
e o bloco `(m+d, m)` de `B·Q0` é a composição **radial 1D** `B_d · R_m^-1`, com
`R_m = J_rad + s + im/2`. Pelo teste de Schur nos blocos de Fourier,

```
α_F <= Σ_d sup_{|m| >= M} || B_d R_m^-1 ||,
```

e cada termo é um problema 1D, onde a composição verdadeira pode ser calculada em
vez de majorada pelo produto de normas.

| `M = 30`, `s = 0` | composição | produto | ganho |
|---|---:|---:|---:|
| `d = 0` (com `μR_x`) | 0,178 | 0,701 | **3,9×** |
| `d = ±2` | 0,188 / 0,140 | 0,217 | 1,2–1,5× |
| `d = ±4` | 0,051 / 0,037 | 0,061 | 1,2–1,6× |
| **soma sobre `d`** | **0,626** | 1,304 | **2,1×** |

Em `M = 20`: 0,879 contra 1,879.

**Leitura.** O ganho grande de `d = 0` vem quase todo de `μR_x`: a divisão por ξ
se alinha mal com as direções em que `R_m^-1` é grande. As multiplicações genuínas
pelo fundo ganham pouco (1,2–1,5×), e agora dominam. A soma triangular sobre `d`
ainda perde ~1,8× das fases entre modos (o verdadeiro com bloco 20×80 é ≈0,5).

**Consequência.** A cauda de Fourier pede `M ≈ 40–50` em vez de ~60. Com as
composições 1D também na cauda radial e a estrutura de faixa de topo
(`||V P_topo||` calculável exatamente), o tile rigoroso fica em ~**15 mil DOFs
complexos** (~3,5 GB por matriz): cabe no PC do laboratório. A estimativa anterior,
20–40 mil, era escala Kaggle.

**O que ainda ganharia mais:** evitar a soma triangular sobre `d`, isto é, uma cota
que trate a convolução de Fourier com fase, como o símbolo fez para `B` sozinho.
`R_m^-1` varia ~50% entre `m` e `m+2`, então congelar `m` não funciona diretamente.

### 10.3 O tile-protótipo rigoroso: o que foi feito e onde parou (24/09)

**Implementado.**

- `scripts/tile_certificate.py`, **parte finita**. Três níveis (`G`, faixa
  `near`, `far`); todos os blocos entre `G` e `near` ficam dentro de `G+`, porque
  `Q0` preserva caixas, e são calculados exatamente. Normas por **cota superior
  verificada** (`t² I − X*X ≻ 0` por Cholesky, com margem para o arredondamento).
  Em `s0 = 1`, `G = 12×36`, `G+ = 20×60`: `||V|| <= 2,64`, `e_Gn <= 0,445`,
  `e_nG <= 0,168`, `e_nn <= 0,454` — os mesmos valores do experimento.
- `scripts/fourier_tail_bound.py`, **cauda infinita** `β_f = ||(B+ds) Q0 P_far||`.
  Composições radiais 1D modo a modo; teste de Schur nos blocos de Fourier (com o
  peso `κ1^|d|` do toro, que eu tinha esquecido na medida da §10.2); `||B_d||` por
  Young; `||R_m^-1 P_{>N}||` por Schur com majorantes da forma fechada (fatores
  `|1 − v_k/d_k| <= 1` quando `Re s0 >= 0`); cota uniforme
  `||R_m^-1|| <= (1 + 2 c_w/(κ−1))/h` para `|m|` grande.

**Onde parou.** Em `s0 = 1`, `M+ = 30`, `N+ = 120`, `N_r = 300`, `M_max = 60`:
**`β_f <= 3,4`**, contra ~0,6 verdadeiro e ~0,3 necessário. O excesso vem das
caudas das caudas, que ainda são produtos de normas:

| peça | contribuição |
|---|---:|
| `|m| > M_max`: `||B_d||_Young × 21,4/h` | ~1,7 |
| `n > N_r`: `||B_d||_Young × ||R^-1 P_{>N_r}||` | ~0,9 |
| composições finitas | ~0,5 |

Empurrá-las exige `N_r ≈ 3000` e `M_max ≈ 1000`: ~7.500 composições de ~4500²,
horas por tile. E ainda sobra o piso da soma triangular sobre `d` (0,63 em
`M+ = 30`, perda de fases), que força `M+ ≈ 60–70`. **O tile rigoroso, com as
ferramentas atuais, custa ~20–40 mil DOFs complexos por tile** — o mesmo que a
estimativa por produto da §10.1. O ganho de 2× da §10.2 é real, mas é comido
pelas caudas das caudas.

**O que mudaria isso.** Uma cota da composição `B·Q0` na cauda que preserve as
fases entre os modos de Fourier. O símbolo pontual fez exatamente isso para `B`
sozinho, e ganhou 3×. Estendê-lo à composição — símbolo matricial em `m` com
controle rigoroso do comutador — é o problema matemático que resta.

### 10.4 A parte matemática: o critério certo é o de Perron (24/09)

**O operador livre como campo no toro.** Com `ξ = (z + 1/z)/2`, vale
`(1+ξ)∂_ξ = z(z+1)/(z−1) ∂_z`, e em ângulos (`w = κ1 e^{iθ}`, `z = κ2 e^{iφ}`)

```
K_livre = ½∂_θ + A(φ)∂_φ + D(φ)(−i∂_φ) + μ + s0,
A = −2μκ sin φ / |z−1|²   (transporte),   D = μ(κ²−1)/|z−1|² > 0   (coercivo).
```

O símbolo `p = im/2 + iA n + D n + μ + s0` **explica as duas caudas medidas**:

- **radial:** a coercividade vem de `D(φ) n`, que é mínima em `φ = π` (junto ao
  cone de luz `ξ = −1`): `D = μ(κ−1)/(κ+1) ≈ 0,0099`, logo `1/|p| ~ 100/n` —
  os ~90/N medidos;
- **Fourier:** ressonância `m/2 = −A n`, pior em `sin φ = ±1`:
  `1/|p| ≈ 4κ/((κ²−1) m) ≈ 17/m` — os ~16/M medidos.

**Tentativa que não funcionou: peso pontual.** A ideia era fatorar
`||B Q0|| <= ||B W^-1|| · ||W Q0||` com `W = w(φ)`: o primeiro fator é um sup
pontual (pega as fases de Fourier) e o segundo é bloco-diagonal em `m` (sem soma
sobre `d`). Deu 0,98 contra 0,88 da soma modo a modo (`M = 20`). O motivo: o fundo
é **máximo** em `φ = π/2` (1,585), exatamente onde a ressonância de Fourier é pior.
Estão alinhados, e nenhum peso espacial os separa.
`scripts/measure_spatial_weight.py`.

**A verdade da cauda de Fourier** (`scripts/measure_fourier_tail_true.py`, matriz de
blocos `(m+d, m)` com `m` até `M+60` e 240 modos radiais):

| `M` | verdadeiro | soma modo a modo | produto |
|---|---:|---:|---:|
| 20 | **0,510** | 0,879 | 1,19 |
| 30 | **0,388** | 0,626 | 0,82 |

Há um ganho de fases de ~1,7× que a soma triangular sobre `d` perde.

**O que resolve: o critério de Perron, não a norma espectral.** Eu exigia
`||N||_2 < 1` para a matriz não negativa `N` de cotas dos blocos de `E = I − A H`.
Basta o **raio espectral** `ρ(N) < 1`: existe `v > 0` com `Nv < v`, e na norma
`max_i ||x_i||/v_i` o operador `E` é uma contração, logo `H` é injetivo. Para dois
níveis isso é `ρ + ab/(1 − α) < 1` (eliminação de Schur). Com números já medidos:

| ponto | `||N||_2` (antigo) | `ρ + ab/(1−α)` |
|---|---:|---:|
| `s0 = 0,2`, bloco 12×36 | 1,02 ✗ | **0,56 ✓** |
| `s0 = 0,301`, bloco 12×36 | 0,95 | **0,44** |
| `s0 = 0`, bloco 20×80 | 0,72 | **0,22** |

**Consequência para o tile rigoroso.** Cada cauda só precisa contrair (`α < 1`),
e as "caudas das caudas" viram níveis próprios da matriz de Perron:

- `G+` (caixa exata, invertida inteira);
- cauda de Fourier intermediária (`M+ <= |m| < M_max`, `n < N_c`): composições
  modo a modo, **radialmente finitas**;
- cauda de Fourier distante (`|m| >= M_max`): cota uniforme de `||R_m^-1||`;
- cauda radial distante (`n >= N_c`, todo `m`): produto `||B|| · ||Q0 P_R||`.

Estimativa grosseira com `G+ = 40×160` (~9 mil DOFs complexos), `M_max = 200`,
`N_c = 800`: `ρ(N) ≈ 0,65–0,8`. **Cabe no PC do laboratório**, e a matemática que
faltava era trocar o critério, não apertar as cotas.

### 10.5 O tile por Perron com casca exata fecha em ponto flutuante (24/09)

`scripts/signed_operator.py`, `scripts/tile_perron.py` (estimativa em ponto
flutuante, não é prova).

**Base certa.** O certificado passou a ser montado na base de Fourier **com sinal**,
isto é, na própria norma do toro em que o símbolo deu `||B♯|| <= 1,65`. O montador
real usa partes (cos, sen) com peso simétrico `ω1(|m|)`, e complexificado isso é
outra norma, equivalente mas diferente. O protótipo da §10.3 misturava as duas. O
operador com sinal foi validado contra o montador real pelo espectro, com
emparelhamento ótimo: mediana `10⁻¹²` em 8×24 e 12×36, e o `−0,401` aparece nos
dois.

**Bug de rigor corrigido em `cauda_Rinv`.** Para `passo = 2`, os pesos `ω(n)` com
`n` até 6001 estouravam, a razão virava `inf/inf = nan`, e o `max()` do chamador
descartava o `nan` em silêncio. O `passo = 1` domina (0,76 contra 0,41 em
`N = 120`), então o β_f da §10.3 saiu certo por acaso. Os pesos agora são
calculados em log, com asserção de finitude.

**Níveis.** `Z` (inverso denso `V`) ⊂ `G+` ⊂ `F`, mais o nível `far`, fora de `F`:

- `T1 = G+ \ Z`: **casca exata**. Todos os blocos que envolvem T1 são matrizes
  finitas calculadas inteiras (`Q0` preserva caixas). A casca **não é
  invertida**: entra só em produtos e normas.
- `T2 = F \ G+`: composições radiais modo a modo, exatas porque `F` é uma caixa.
  O bloco `T2 → T2` usa Schur sobre `d`; os blocos `T2 → Z` e `T2 → T1` são
  conjuntos (Gram acumulado).
- `far`: `||B|| = 1,65` vezes `||Q0 P_far||`. Com `Mc >= MP + 40`, a parte de
  Fourier distante não alcança `G+`, e a radial só alcança descendo de
  `n >= Nc`, o que custa `κ2^-(Nc − NP − 101)` pela forma fechada (~10⁻⁵).

`s0 = 0`, `r = 0`, `F = 200×400`:

| configuração | `Z ← T` | `T ← Z` | bloco de T | raio de Perron |
|---|---:|---:|---:|---:|
| `Z = 12×36`, sem casca | 0,97 | 0,30 | 1,82 | 2,22 |
| `Z = 20×80`, sem casca | 0,48 | 0,23 | 1,17 | 1,27 |
| `Z = 30×120`, sem casca | 0,33 | 0,16 | 0,87 | **0,950** |
| **`Z = 20×80`, casca `40×160`** | 0,48 | 0,23 | **0,52** (T1, exato) / 0,69 (T2) | **0,866** |

Com a casca, o inverso denso tem só **2280 DOFs** e a casca, 9360. O vetor de Perron
é `(0,55; 1; 0,87; 0,07)`. Folgas: α22, a única peça modo a modo, pode crescer
1,29×, α11 1,5×, e `a1` 3,4×.

**Onde o modo a modo perdia.** Sem casca, `d = 0` atinge o sup em `m = 0` (a cauda
radial do modo baixo, colada em Z) e `d = ±2` na borda de Fourier: a soma de sups
junta máximos de lugares diferentes. Um nível por modo, ou faixas de modos, **piora**
(2,55 contra 2,22 em 12×36): as linhas de Z e de `far` viram somas sobre modos, e
`sqrt(Σ a_m²) = 2,18` contra `a = 0,97` conjunto, porque `V` espalha. A casca exata
põe a fronteira do modo a modo longe de Z, onde as composições já são pequenas
(`d = 0`: 0,52 → 0,31).

**O custo que vai ditar o número de tiles é `ds`.** `N = N0 + r N1`; na linha de Z,
`r ||V Q0_ZZ||` entra como `1/(1 − r ||V Q0_ZZ||)`, e `||V Q0_ZZ||` é a norma do
resolvente finito: ~18 em `s0 = 0` (caixa pequena), inflada pela não normalidade,
com o autovalor a 0,401. Isso dá `r ≈ 0,03–0,04` perto do eixo imaginário.
O mapa `r_adm(s0)` sobre a região está sendo medido no PC do laboratório.

**O que falta para a versão rigorosa** (todas peças de engenharia, sem matemática
nova):

1. normas verificadas (Cholesky com margem de arredondamento) no lugar da
   iteração de potência, inclusive nos Gram acumulados;
2. erro de `fl(B Q)` e de `R_m^-1` calculado (resíduo `J R − I`), somado aos
   blocos;
3. erro do fundo de RT (`2^-25` na norma ℓ¹ sharp, que domina a do toro) como
   `ε ||Q0 P_j||` por coluna;
4. o vetor de Perron exibido, com `N v < v` verificado em aritmética intervalar;
5. a cobertura da região `0 < Re s <= 1,765`, `0 <= Im s <= 1/4` (K é real), fora
   do disco em 0,401, por discos `|s − s0| <= r_adm(s0)`.

### 10.6 Do tile à cobertura: modo rigoroso, banda truncada, três centros (24/09)

**Normas verificadas** (`scripts/verified_norms.py`). Cota superior de `||X||_2` em
ponto flutuante: `t` vale se `t² I − X*X` é definida positiva, verificada pelo
critério de Rump (BIT 2006): o Cholesky em ponto flutuante de `A − cI` termina, com
`c >= γ_{n+1}/(1 − γ_{n+1}) tr(A)` mais o termo de underflow. O erro de formar `X*X`
em ponto flutuante (`γ_k ||X||_F²`) entra somado a `c`. Para o caso complexo, uso
`γ` com `n → 4(n+1)`. A razão cota/verdade fica entre `1 + 10⁻⁶` e `1 + 7·10⁻⁴`.

**Contabilidade de erros do tile** (`tile_perron.py --rigoroso`):

- **pesos**: calculados com erro relativo `<= 4u`; a norma com eles difere da do
  toro por `<= 1 + 10u`, e todas as cotas levam o fator `1 + 10⁻¹²`;
- **fundo**: truncado em `|d| <= 20`, `|dn| <= 40` **para fora**. A parte de Hankel
  e `μR_x` só descem em `n`; o resto tem Young `<= 2,2·10⁻⁵`. Entram numa perturbação
  global `dB` o resto, o arredondamento das entradas (`γ_8 · Young_abs`, com campos
  em valor absoluto) e o erro de RT (`40·2⁻²⁵`, o mesmo que o símbolo usou).
  `μR_x` é exato (μ é diádico). `dB` entra em todas as colunas, inclusive
  `far ← G+` e `G+ ← far`, que antes eu tinha zerado: o fundo truncado não liga
  esses níveis, mas a perturbação global tem suporte completo;
- **`R_m⁻¹`**: pelo resíduo com pesos, `||R_exato − R|| <= ||R|| ||E||/(1 − ||E||)`;
- **produtos e Grams**: `||fl(AB) − AB|| <= γ_k ||A||_F ||B||_F`; Grams acumulados
  com `γ` da soma das dimensões internas;
- **Perron**: o vetor `v` é exibido e `N v <= θ v` é checado em racionais exatos.

No caso de teste, as cotas rigorosas diferem das estimadas em ~10⁻⁵ (dominadas pelo
resto da banda e pelo erro de RT).

**Três centros.** `K_livre` é bloco-triangular superior na decomposição
`[Z, T1, T2, far]`, porque preserva as caixas. Então

```
P = Q0(sZ) P_Z + Q0(s1) P_T1 + Q0(s2) P_(T2+far)
```

é bloco-triangular com diagonais invertíveis, uma bijeção sobre o domínio, e

```
H(s) = K(s) P = I + Σ_ℓ (B + s − s_ℓ) Q0(s_ℓ) P_ℓ.
```

Um tile `|s − sZ| <= rZ` exige `|s − s1| <= r1` e `|s − s2| <= r2` nele. **Só as
colunas de Z são sensíveis a `s`**: `||V Q0_ZZ|| ≈ 18` em `s0 = 0`, contra 0,83
(T1), 0,50 (T2) e 0,24 (far). T1 e T2 saem uma vez por tile médio ou grande, com
cache em disco. Por tile pequeno sobra o que envolve `V`: em 12×36 com caches, menos
de 1 s.

**Por que basta `Im s >= 0`.** A norma com sinal `Σ κ1^{2m}|x_m|²` **não** é
invariante pela conjugação `x_m → conj(x_{−m})`, que troca o toro de raio κ1 pelo de
raio 1/κ1. Mas o espaço simétrico `Y` (pesos `κ1^{2m} + κ1^{−2m}`, o natural para
S3) é invariante e está contido no espaço com sinal. Injetividade no espaço com sinal
para `Im s >= 0` dá injetividade em `Y` para `Im s >= 0`, e a conjugação em `Y` leva
isso a `Im s <= 0`.

**Mapa de `r_adm`** (um centro só, alvo: raio <= 0,95; `Z = 20×80`, `G+ = 40×160`,
`F = 200×400`; PC do laboratório):

| `s0` | raio (`r = 0`) | `||V||` | `r_adm` |
|---|---:|---:|---:|
| 0 | 0,866 | 7,3 | 0,028 |
| 0,25i | 0,864 | 10,4 | 0,019 |
| 0,2 | 0,810 | 5,4 | 0,061 |
| 0,2 + 0,25i | 0,808 | 6,2 | 0,059 |
| 0,401 + 0,15i | 0,761 | 7,1 | 0,067 |
| 0,6 | 0,721 | 5,8 | 0,109 |
| 1 + 0,25i | 0,655 | 2,5 | 0,38 |
| 1,765 | 0,569 | 1,8 | 0,88 |

A faixa difícil é estreita, junto ao eixo imaginário: com `Re s = 0,2`, o raio
admissível já triplica. A sensibilidade `||V Q_ZZ||` cai de 18 (eixo) para 11 a 0,15
do autovalor, e para 0,9 em `Re s = 1,765`. Quem domina é a não normalidade junto ao
eixo, não o polo em 0,401. Todas as peças melhoram com `Re s`: α11 0,52 → 0,31,
α22 0,69 → 0,46.

**Planejador** (`scripts/tile_cover.py`): quadtree em três escalas. Quadrados
grandes dão o centro de T2, os médios o de T1, e os pequenos são os tiles de Z,
subdivididos enquanto o Perron verificado não fecha. O disco circunscrito de um
subquadrado cabe no do quadrado-pai, então as condições de centro valem por
construção. Quando a própria base (`rZ = 0`) não fecha, o planejador para de
subdividir e acusa o quadrado médio.

### 10.7 Primeiro tile rigoroso (24/09)

`tile_perron.py --s0 0.2 --r 0.04 --Z 20x80 --Gp 40x160 --F 200x400 --rigoroso`
(`build/tiles_rig/tile-0_2-r0_04.json`), 27 min em 4 threads:

| bloco | cota |
|---|---:|
| `||V||` / `ρ` | 5,363 / 3,9·10⁻⁴ |
| `||V Q_ZZ||` | 10,55 |
| `Z ← T1` / `T1 ← Z` / α11 | 0,4297 / 0,1998 / 0,4851 |
| `T1 ← T2` / `T2 ← T1` / α22 | 0,2361 / 0,1153 / 0,6519 |
| `||Q0 P_far||` / β | 0,2333 / 0,3850 |
| raio de Perron, `r = 0` / tile | 0,8106 / 0,8752 |

Vetor exibido `v = (0,988; 1; 0,664; 0,053)`, com `max_i (N v)_i / v_i = 0,8752 < 1`
conferido em racionais: **`K(s)` é injetivo em `|s − 0,2| <= 0,04`**, no espaço com
sinal e, pela §10.6, no simétrico.

**Estado do rigor.** Ponto flutuante com cotas de erro a priori (Rump para
positividade, `γ_n` para produtos), não aritmética intervalar. As entradas
externas são `||B♯|| <= 1,64967` (símbolo, §5.26 das rotas alternativas) e as cotas
da forma fechada de `R_m⁻¹` (`cauda_Rinv`, `Rinv_uniforme`, `descida`), que somam
em ponto flutuante com margem relativa 10⁻⁴. A margem é argumentada (~10⁴ termos
positivos), não intervalar. **Antes de contar como prova, precisa de auditoria**:
releitura independente da contabilidade de erros, do lema de Perron e da
completude dos blocos, que conferi linha por linha, mas sozinho.

**Custo da cobertura.** Com os três centros, por tile pequeno sobra a parte Z
(~40 s) mais os produtos com `V`. As partes T1 (~6 min) e T2 (~25 min) saem uma vez
por quadrado médio e grande. Em `s0 = 0` (o pior ponto), `rZ = 0,015` e `r1 = 0,03`
ainda admitem `r2 ≈ 0,2`: quadrados grandes de lado 0,25 servem na região toda,
**~8 partes T2**.

**Correção no critério (24/09, mesma noite) — RETIFICADA pela auditoria.** *O fundamento abaixo está errado: o Teorema 2.3 de Rump cobre também hermitianas complexas, com as mesmas constantes; a versão complexa era válida e a realificação é só mais conservadora. Ver `AUDITORIA_S3A_TILES_24SET.md`.* Texto original: relendo o modo rigoroso: o critério de
Rump é um teorema para matrizes **reais** simétricas, e eu o tinha aplicado ao
Cholesky complexo com um `γ` inflado "com folga". Isso é argumento, não teorema.
Agora `verified_norms.definida_positiva` aplica o teorema à realificação
`[[Re A, −Im A], [Im A, Re A]]` (real simétrica `2n × 2n`, positiva se e só se `A`
é), como enunciado. A subtração `a_ii − c` entra na margem explicitamente. Os
produtos complexos continuam com `γ_{4(k+1)} >= √2 γ_{k+2}` (Higham, Lema 3.5). O
tile acima foi calculado antes da correção e está sendo refeito; a cobertura do
laboratório já roda com a versão corrigida.

### 10.8 Conferência independente das cotas de blocos (24/09)

`scripts/test_tile_perron_blocks.py`. Para linhas e colunas dentro de `F`, o bloco
`P_i E P_j` do operador infinito coincide com o da caixa `F`: `Q0` é triangular e
diagonal em `m`, e `F` é uma caixa. O teste monta na caixa, **com o fundo completo**
(sem truncar a banda), `H(s) = K(s) P` com os três centros e o mesmo `V` do
certificado, e compara a norma **verdadeira** (SVD) de cada bloco de `E = I − A H(s)`
com a cota do certificado, em 6 pontos de cada tile (centro e bordo).

| configuração (`Z 8×24`, `G+ 12×36`, `F 32×77`) | pior verdade/cota |
|---|---:|
| `sZ = 0,32`, `s1 = 0,3`, `s2 = 0,25` | 0,9999 |
| `sZ = 0,02 + 0,1i` (junto ao eixo), `s1`, `s2` deslocados | 0,9999 |
| `sZ = 1 + 0,2i`, `rZ = 0,1` | 0,9999 |

Nenhuma cota é violada. Os blocos calculados exatamente (`T1 ← Z`, `T2 ← Z`,
`T2 ← T1`) ficam a ~10⁻⁴ da verdade, o tamanho dos termos que cobrem o resto da banda.
`Z ← Z` na borda do tile é atingido com igualdade (`rZ ||V Q_ZZ||` é justo). O
modo a modo `T2 ← T2` perde 2,4× nesta caixa pequena. As entradas de `far` não cabem
na caixa: são as cotas da forma fechada e do símbolo, conferidas à parte.

**As entradas de `far`** (`scripts/test_far_bounds.py`), contra a norma verdadeira
numa truncagem radial de 1500:

| cota | verdade / cota (pior) |
|---|---|
| `cauda_Rinv(σ, N)` (`N = 160, 400`; `m = 0…−100`; três `s`) | ~0,90–0,95: justa |
| `Rinv_uniforme` (`|m| >= 60, 200`) | ~0,38: 2,5× folgada, e hoje a parte radial domina β |
| `descida(J, Nc)` | `10⁻¹³` contra `10⁻⁹`: folgada, e desprezível |
| `||B||` numa caixa 24×60 (fundo completo) contra o símbolo | 1,507 contra 1,650 |

Pior razão 0,965. Nenhuma cota violada.

### 10.9 Verificação independente da cobertura (24/09)

`scripts/check_cover.py` não confia no planejador nem na marca `ok` gravada pela
rodada. Para cada tile, a partir só do JSON, ele:

1. **reverifica o certificado em racionais exatos**: `N = N0 + rZ NZ + r1 N1 + r2 N2`
   (as matrizes gravadas) e o `v` gravado dão `N v <= θ v` com `θ < 1`; e as
   contenções `|sZ − s1| + rZ <= r1` e `|sZ − s2| + rZ <= r2` valem;
2. só então usa o disco. A cobertura é checada por uma quadtree própria, também em
   racionais: cada célula da região deve estar inteira no disco **aberto** D de
   S3b ou inteira num disco reverificado (os 4 cantos bastam, por convexidade).

**Uma lacuna real, microscópica, achada assim.** O planejador usava o disco
circunscrito exato de cada quadrado (`h/√2`). Os vértices da malha ficam **sobre**
os círculos, e `h/√2` arredondado para baixo os deixa de fora. Numa malha
sintética, o verificador reprova exatamente nos vértices. Agora todos os raios
(`rZ`, `r1`, `r2`) são inflados por `1 + 10⁻⁶`. Com o mesmo fator nos três, as
contenções continuam valendo, e a malha sintética passa. A cobertura do laboratório
foi reiniciada com isso; os caches de T1 e T2 não dependem do raio.

### 10.10 Perto do disco de S3b (24/09)

Parte Z em pontos a distância `d` do autovalor 0,401 (estimativa):

| `d` | `||V||` | `||V Q_ZZ||` | `rZ ~ 0,5/||V Q_ZZ||` |
|---:|---:|---:|---:|
| 0,15 | 6,5–7,4 | 9,8–12,7 | 0,039–0,051 |
| 0,10 | 9,5–10,6 | 15,0–17,7 | 0,028–0,033 |
| 0,07 | 13,6–14,8 | 21,9–24,5 | 0,020–0,023 |
| 0,05 | 19,2–20,4 | 31,0–33,6 | 0,015–0,016 |
| 0,04 | 24,1–25,3 | 39,0–41,6 | 0,012–0,013 |

Perto do polo, `||V Q_ZZ|| ≈ 1,6/d`, então `rZ ≈ 0,3 d`: o número de tiles junto a D
cresce só como `log(1/raio de D)`.

**A base não piora.** Temi que `Z ← T1 = ||V H_Z,T1||` crescesse com `||V|| ~ 1/d`, e
o produto `(Z ← T1)(T1 ← Z)` quebrasse o Perron perto de D, qualquer que fosse `rZ`.
Não acontece. Com T1 exato em cada ponto e T2 aproximado pelo de `0,401 + 0,15i`:

| `s` | `d` | `||V||` | `Z ← T1` | base (`r = 0`) |
|---|---:|---:|---:|---:|
| 0,551 | 0,15 | 7,35 | 0,364 | 0,749 |
| 0,501 | 0,10 | 10,59 | 0,372 | 0,753 |
| 0,451 | 0,05 | 20,38 | 0,381 | 0,758 |
| 0,401 + 0,05i | 0,05 | 19,85 | 0,390 | 0,762 |

A direção em que `V` explode (a do autovetor) é quase ortogonal à imagem do
acoplamento com a casca. Consequência para S3b: o raio de D pode ser escolhido pelo
que S3b precisar. S3a paga só `rZ ≈ 0,3 d` na borda.

### 10.11 Auditoria do modo rigoroso (24/09)

`AUDITORIA_S3A_TILES_24SET.md`: 16 alegações com veredicto. Achados: minha afirmação
errada sobre Rump (retificada acima); hipótese de simetria exata; margem de `c`; somas
de cotas em ponto flutuante; o termo de `μR_x` em `norma_Bd_young` não é cota provada
(retirado do certificado); vértices da malha. Todos corrigidos. Os tiles já calculados
continuam válidos pelas folgas existentes. Pendências: a folga de Lipschitz do símbolo
em Arb, o raio de D, a cobertura completa e uma auditoria por revisor independente.

### 10.12 A região finita de S3a está coberta (24/09)

Cobertura rigorosa no PC do laboratório: 8 quadrados grandes de lado 0,25, 30 quadrados
médios, 156 tiles calculados, **154 certificados** (os 2 que falharam, junto ao disco
D, foram subdivididos). Pior raio de Perron 0,9961 (junto ao eixo imaginário); menor
`rZ` 0,0147 (junto a D).

`scripts/check_cover.py build/cobertura_lab` reverifica cada tile em racionais exatos a
partir do JSON (`N v <= θ v` com `θ < 1`, contenções dos discos dos centros, com a
inflação de `10⁻⁹`) e confirma, por quadtree própria em racionais:

> **`K(s)` é injetivo em `{0 <= Re s <= 1,765, 0 <= Im s <= 1/4}` fora do disco aberto
> `|s − 0,401| < 0,05`**, no espaço com sinal. Pela §10.6, isso vale também no espaço
> simétrico `Y`, e **em `Y`** (não no espaço com sinal) também em `−1/4 <= Im s <= 0`.

Com `R♯ <= 1,765` (§9: sem espectro em `Re s > 1,765`), isso é **S3a inteira**, sujeita às
pendências da auditoria (`AUDITORIA_S3A_TILES_24SET.md` §3): o raio de D fixado com S3b
e uma revisão independente. (A folga de Lipschitz do símbolo, que também entra em `R♯`,
foi refeita em Arb em 24/09 e é coberta pela gravada.)

Os quadrados rodaram com versões do código anteriores e posteriores às correções da
auditoria. As anteriores continuam válidas pelas folgas descritas na auditoria (§2), e
o `check_cover` aplica a todos a mesma reverificação. Custo: ~11 h de CPU em 1 thread por
processo, 3 processos.

### 10.13 O disco de S3b: raio 1/8 (24/09)

`S3B_DISCO.md`: o raio que minimiza `max ||V_L||` em ∂D é `~1/8` (880 em 12×36, contra 1846 em
`Γ_A`), com folga de 0,108 até a raiz vizinha de L em 0,1683. A cobertura de S3a (feita com
0,05) já atende: `check_cover --disco 0.401,0.125` dá COMPLETA. Os scripts agora usam `1/8`
como padrão.

### 10.14 Revisão independente: duas lacunas, fechadas (24/09)

`AUDITORIA_S3A_TILES_24SET.md` §6. Um revisor com contexto limpo (outro modelo) achou:

1. o caminho de `μR_x` de `far` até G+. `μR_x` desce de qualquer `n`; a cota nova é
   `descida_rx`, de ~10⁻¹⁴ a 10⁻¹⁰;
2. **o μ verdadeiro de RT**, que difere do diádico de `RefA` em até `3,88·10⁻¹¹`. A cota nova é
   `eps_mu`, ~10⁻⁷ por coluna.

Somadas aos 154 tiles, o pior θ vai a 0,9961189, e a cobertura continua completa. O revisor
também deu o argumento de domínio maximal de `P`, que a §10.6 omitia: `P = Q0(s2)·T` com
`T − I` de posto finito e bloco-triangular com diagonais invertíveis, e `K_livre + σ` é
injetivo em X porque `Π_{k>=j}(1 − v_k/d_k) → 0`.
