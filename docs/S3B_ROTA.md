# S3b: a rota local (24/09/2026)

**O que S3b precisa garantir.** Entre as raízes de L em `Re s > 0`, contar exatamente quantas
violam as constraints, porque essas não são físicas. S3a já diz que fora do disco
`D = {|s − 0,401| < 1/8}` toda raiz de L satisfaz as constraints.

## 1. O argumento

Se `L(λ)h = 0` e `C(λ)h ≠ 0`, então `K(λ)C(λ)h = M(λ)L(λ)h = 0`: um autovetor de L que viola as
constraints **só pode estar numa raiz de K**. Daí:

1. **K tem uma única raiz em D, `s_K`, simples.** Newton–Kantorovich (NK) para K em `s_K`
   (K é bem condicionado: `‖V‖ ≈ 7–20`) dá unicidade num disco `D(s_K, ρ)`. Os tiles de S3a,
   estendidos de raio 0,05 até `ρ`, dão K injetivo no resto de D.
2. **L tem, perto de 0,401, um autovalor `λ_L` algebricamente simples, com autovetor
   enclausurado.** Um certificado NK do operador bordejado `[[L(λ), h],[ℓ, 0]]`. Como
   `L(s) = L_0 + s`, a derivada segunda é trivial.
3. **Testemunho:** `C(λ_L)h ≠ 0`, pelo enclausuramento de `h` (a constante `c_op` de H3,
   refeita nas normas do certificado).
4. **Fechamento:** por 3 e pelo argumento acima, `λ_L` é raiz de K, logo `λ_L = s_K`. Em `s_K`,
   L é algebricamente simples e viola as constraints. Qualquer outra raiz de L em D satisfaz
   as constraints (K é injetivo ali) e entra na contagem global de T2, como as de fora de D.

**Conclusão de S3b:** exatamente uma raiz de L em `Re s > 0` (com multiplicidade) viola as
constraints, e ela está em `s_K ≈ 0,401`. Isso **não** exige contar as raízes de L em D (a
antiga H1): exige um certificado de L num ponto, não ~400 tiles na borda.

## 2. Por que não contar L em D

Montador de L na base com sinal (`scripts/signed_operator_L.py`), validado entrada a entrada
contra o exportador de RT (≤ 1,8·10⁻¹⁵ em 4×12, 6×18, 8×24). Experimento de Perron
(`scripts/experiment_perron_L.py`, ponto flutuante) em `s = 0,276`, o pior ponto de ∂D, na
norma do toro (65/64, 5/4):

| Z | `‖V‖` | `‖V Q_ZZ‖` | `a` | `b` | `α` (cauda) | `ab/(1 − α)` |
|---|---:|---:|---:|---:|---:|---:|
| 12×48 | 96,9 | 159,9 | 0,94 | 0,84 | 1,31 | — |
| 16×64 | 96,9 | 159,9 | 0,71 | 0,67 | 1,08 | — |
| 20×80 | 96,9 | 159,9 | 0,58 | 0,56 | 0,93 | 4,4 |
| 24×96 | 96,9 | 159,9 | 0,50 | 0,48 | 0,81 | 1,26 |

- **`α ≈ 20/M`, e a cauda é de Fourier** (Z = 12×48: Fourier←Fourier 1,30, radial←radial 0,91).
  Em L, a componente 3 vive em `m` ímpar, e o fundo acopla modos a distância 1, não 2 como em
  K.
- **L fecha com Z ≈ 30–32 × 120–128**, um inverso denso de ~14 mil DOFs (em K eram 2.280).
- **`‖V Q_ZZ‖ ≈ 160`** dá tiles de raio ~0,002–0,003 em ∂D: ~300–400 fatorações desse
  tamanho, dias de máquina.
- Passar `μSΓ1` para a parte livre não ajuda (`α` 1,32): o acoplamento vem dos campos de `Γ2`.

## 3. Ordem de trabalho

1. Viabilidade (ponto flutuante) do certificado bordejado de L em 0,401: o raio de Perron
   do operador bordejado e o resíduo do autovetor aproximado. Daí o `ε_v` previsto, contra
   o necessário para o testemunho.
2. Porte da maquinaria rigorosa para L (4 componentes, κ = (65/64, 5/4), `‖B_L‖ <= 3,964`
   do símbolo em Arb), com o nível bordejado.
3. NK para K em `s_K` e extensão dos tiles de S3a até `ρ`.
4. Testemunho nas normas do certificado.

## 4. Viabilidade do certificado bordejado de L (24/09, noite)

`scripts/experiment_bordered_L.py` (ponto flutuante), no próprio autovalor
`λ0 = 0,4010247330`; cauda aproximada por uma banda G' além de Z:

| Z (G') | `‖V_b‖` | `a_b` | `b` | `α` | `ab/(1 − α)` | resíduo fora de Z |
|---|---:|---:|---:|---:|---:|---:|
| 12×48 (20×80) | 20,7 | 0,858 | 0,785 | 1,236 | — | 2,9·10⁻² |
| 20×80 (28×112) | 20,7 | 0,555 | 0,538 | 0,892 | 2,76 | 8,0·10⁻⁴ |
| 24×96 (32×128) | 20,7 | 0,477 | 0,465 | 0,783 | 1,02 | 1,2·10⁻⁴ |
| **28×112 (36×144)** | **20,7** | **0,419** | **0,409** | **0,699** | **0,57** | **1,8·10⁻⁵** |

- **O bordejamento remove a direção singular:** `‖V_b‖ = 20,7` no autovalor, contra
  `‖V‖ ≈ 97` na borda de D.
- **Com Z = 28×112** (~10.700 DOFs), o critério fecha com θ ≈ 0,89, e o raio NK previsto é
  1,6·10⁻⁴.
- **Testemunho** (`scripts/signed_constraints.py`, `C(s)` validado contra RT a 1,8·10⁻¹⁵):
  `‖C h‖/‖h‖ = 0,17174` nas normas do toro, estável desde 12×48. A queda de 17/09 era da
  normalização ℓ¹. A folga sobre o raio NK é de ~100×.
- **K e L têm o mesmo autovalor** até 5·10⁻¹¹ em 20×80 (6·10⁻⁸ → 10⁻⁹ → 5·10⁻¹¹).

**A ressalva.** A cauda aqui é uma banda finita. Na versão rigorosa, a cauda distante vai
por cota modo a modo, que em K perdeu ~2×. Para L, isso pede uma casca exata maior, dividida
em subníveis para caber na memória, ou blocos exatos por janelas de modos. É a próxima peça
de engenharia.

## 5. A cauda distante de L: modo a modo não serve; janelas servem (24/09, noite)

Na banda de `M = 20` (`n < 120`), no autovalor:

| cota da cauda de Fourier | valor | perda sobre a verdade |
|---|---:|---:|
| verdade (bloco inteiro da banda) | 0,893 | — |
| modo a modo (soma sobre `d` do sup em `m`), `experiment_modewise_L.py` | 2,36 | 2,65× |
| **janelas de 8 modos** (blocos exatos janela ← janela, Perron), `experiment_windows_L.py` | **1,26** | **1,4×** |

A perda do modo a modo vem dos `d` **ímpares**: `d = ±1` vale 0,42 cada, contra 0,54 de
`d = 0`. É o acoplamento da componente 3 (`m` ímpar) com as outras. Nas janelas, esses
acoplamentos ficam dentro de blocos exatos. O bloco da primeira janela sozinho vale 0,889
(quase a verdade), o vizinho ~0,40, e as janelas seguintes decaem como `1/M` (0,70, 0,54).

**Desenho do certificado rigoroso de L** (um ponto, NK bordejado):

- `Z'`: Z = 28×112 mais o bordejamento `λ` (~10.700 DOFs, denso);
- casca: blocos exatos em faixas, sem montar a matriz inteira;
- cauda de Fourier até `|m| ~ 200`: **janelas** de 8–16 modos, blocos exatos e Perron por
  janela (a linha de Z' usa o Gram conjunto, como em K);
- `far`: `‖B_L‖ <= 3,964043` (símbolo, em Arb) e a forma fechada com κ2 = 5/4;
- as correções já conhecidas: o caminho de `μSΓ1` (quatro termos de divisão por ξ), o μ
  verdadeiro de RT, o erro de RT e o resto da banda;
- NK: `Y0` (resíduo, ~2·10⁻⁵), `Z0` (Perron), `Z1` (trivial: L é afim em `λ`); testemunho
  `‖C Q0‖·ε <= ‖C h0‖` (0,172).

Custo estimado: 1–2 dias de porte, algumas horas de PC do laboratório por certificado,
e uma revisão independente com artefatos.

## 6. O testemunho, dimensionado (24/09, noite)

Com `h = Q0 y`: `‖Π_F C(λ)h‖ >= ‖Π_F C(λ0)h0‖ − ‖Π_F C(λ0)Q0‖·‖y − y0‖ − |λ − λ0|·‖ξ‖·‖h‖`,
com a saída na **norma de L** (a de K tem pesos maiores nos modos negativos, pelo fator
`(130/129)^{2|m|}`, e não seria controlada) e projetada numa caixa finita `F`.

- `‖C(λ0)Q0‖` **sem projeção** é grande e ainda cresce (31, 41, 47 em 8×32, 12×48, 16×64).
  O termo principal de `C` tem `(1 − ξ)ξ∂ξ`, enquanto a parte livre inverte `(1 + ξ)∂ξ`, e
  `|1 + ξ|` chega a 0,025 no toro de L, junto ao cone de luz.
- **Com projeção** (`|m| < 12`, `n < 48`): `‖Π_F C(λ0)Q0‖ = 43,7`, convergido no corte das
  colunas (`n < 96`, 160 e 240 dão o mesmo valor).
- `‖Π_F C h0‖/‖h0‖ = 0,184`, já com a caixa `|m| < 4`, `n < 16`.

O testemunho passa se `ε_y` for menor que ~4·10⁻³; o raio NK previsto é ~20× menor que isso.

## 7. A contagem local: tile bordejado com Rouché escalar (24/09, noite)

**O lema do bordejamento.** Se `M_b(λ) = [[K(λ), k0], [ℓ, 0]]` é invertível num disco, então
todo autovetor de `K(λ)` tem `ℓ(k) ≠ 0` (senão `M_b(λ)(k, 0) = 0`), o núcleo tem dimensão
`<= 1`, e os autovalores de K no disco são **exatamente os zeros**, com multiplicidade, da
função escalar analítica `t(λ) = e_λ^T M_b(λ)⁻¹ e_λ`. O mesmo vale para L.

**K é ideal para isso** (Z = 20×80, `Q0` centrado em `s_K`):

| `s` | `‖V_b‖` | `‖V_b Q0‖` | `t(s)` |
|---|---:|---:|---:|
| `s_K = 0,401024733` | 3,68 | 4,85 | ~10⁻¹⁴ |
| `s_K + 0,02` | 3,53 | 4,66 | −0,0200 |
| `s_K + 0,05`, `s_K + 0,05i`, `s_K − 0,05` | 3,3–4,1 | 4,4–5,4 | `−(s − s_K)` |

Sem bordejamento, `‖V Q0‖` era 11–30 no mesmo disco. `t` é praticamente **linear**.

**Enclausurar `t` sem a cota bruta de Neumann.** Com `A = diag(V_b, I, …)` e
`E = I − A H_b`, `t = e^T V_b e + e^T (I − E)⁻¹ E V_b e`. A cota bruta,
`θ‖V_b‖/(1 − θ) ≈ 4`, seria inútil contra `|t| ≈ 0,05`. Mas `E V_b e` é o acoplamento com a
cauda da solução aproximada do sistema bordejado, `−P_T H_b P_Z' V_b e`, da ordem do resíduo
do autovetor (~10⁻⁶). A parte `e^T V_b(λ) e` é um bloco finito, calculado nos pontos da
circunferência. Rouché escalar: `|t − t_aprox| < |t_aprox|` na borda dá exatamente um zero.

**S3b fica assim** (revisto em §8: K por Rouché de operadores, L só por NK):

1. **K:** tile(s) bordejado(s) de raio ~0,05 em torno de `s_K` e Rouché escalar dão
   exatamente uma raiz, simples. Com a cobertura de S3a (fora de `|s − 0,401| < 0,05`), é a
   **única raiz de K em D**.
2. **L:** NK em `λ_L` (existência, enclausuramento do autovetor), mais um tile bordejado
   pequeno em torno de `λ_L` com Rouché escalar, dão multiplicidade algébrica 1.
3. **Testemunho:** `Π_F C h ≠ 0`, logo `λ_L` é raiz de K, logo `λ_L = s_K`.
4. As demais raízes de L em D satisfazem as constraints (K é injetivo fora de `s_K`) e entram
   na contagem global de T2.

## 8. K: Rouché de operadores na circunferência (24/09, noite)

**Por que não o Rouché escalar de §7 sobre o operador infinito.** A correção de cauda de `t`
(o termo `e^T (I − E)⁻¹ E V_b e`) só é pequena **perto** do autovalor, onde `V_b e` é quase o
autovetor. Na circunferência, `V_b(λ) e` já não é, e a correção fica da ordem do próprio `t`.
O bordejamento serve para a contagem **finita**; a passagem do finito ao infinito é feita por
Rouché de operadores.

**Desenho** (`scripts/rouche_K.py`). `P` FIXO: `Q0(c)` nos três níveis, `c = 0,401`, e
`H(s) = K(s)P = I + (B + s − c) Q0(c)`. Referência

```
A(s) = diag(Ĥ_ZZ(s), I, I, I),   Ĥ_ZZ(s) = I + P_Z (B_t + s − c) Q0(c) P_Z,
```

com o fundo truncado `B_t` (dados de RT em aritmética exata) e `μ = μ_RefA` (diádico).
**Gohberg–Sigal:** `H = I + compacto` e `A = I + posto finito` são holomorfas; se
`‖A⁻¹(H − A)‖ < 1` em `|s − c| = ρ`, H e A têm o mesmo número de zeros no disco, com
multiplicidade. Os zeros de A são os de `det Ĥ_ZZ(s) = det K̂_Z(s)/det J_Z(c)`, isto é, os
autovalores da matriz finita `K̂_Z = J_Z + P_Z B_t P_Z` no disco. (`K_livre` é triangular
superior em `n`: `Q0 P_Z = P_Z Q0 P_Z` e `Q_ZZ = J_Z(c)⁻¹`.)

**A condição na circunferência sai do Perron de tiles.** Tile `j`: centro `p_j = c + dV_j`,
raio `rZ`, `V_j ≈ Ĥ_ZZ(p_j)⁻¹` (`tile_perron.parte_Z` ganhou o deslocamento `dV`: o
precondicionador fica em `c` e `V` anda). T1 e T2 saem uma vez, em `c`, com raio
`R = ρ + rZ`. Com `N` a matriz de Perron do tile e `e = N_ZZ ≥ ‖I − V Ĥ_ZZ(s)‖`:

```
‖Ĥ_ZZ⁻¹ H_Zk‖ <= N_Zk/(1 − e)   (k ≠ Z),
‖Ĥ_ZZ⁻¹ (H_ZZ − Ĥ_ZZ)‖ <= ‖V‖ (dB ‖Q_ZZ‖ + eps_mu(c))/(1 − e),
```

e as linhas de T1, T2 e far são as de `N`. Se `N v <= θ v`, a linha de Z de `N'` fica
`<= (θ − e)/(1 − e) v_Z <= θ v_Z`: **o Perron do tile já dá Rouché**, a menos da entrada
`(Z, Z)`, que é só a perturbação (~10⁻³). `N' v <= θ' v` é verificado em racionais.

**A norma com sinal não é simétrica pela conjugação** (§10.6 de S3_CONDICIONAMENTO_SHARP):
a circunferência inteira é coberta, não só a metade superior. A cobertura por `n` discos de
raio `rZ = 2ρ sen(π/2n)(1 + 10⁻⁶)` é conferida em Arb: `|ρ e^{iφ} − dV_j|²` é
quase-convexa em `φ`, logo o máximo em cada arco está nas pontas.

**Contagem finita, rigorosa** (`rouche_K.py contagem`, Z = 20×80, dimensão 2.280, local em
34 s). Bordejado com os autovetores aproximados direito e esquerdo,
`M(s) = [[K̂_Z(s), r], [l^H, 0]]`, `t(s) = det K̂_Z(s)/det M(s)`. Com `V_b ≈ M(c)⁻¹`,
`E(s) = E_c − (s − c) V_b D`:

| `‖V_b‖` | `ε_c` | `β = ‖V_b D‖` | `q = ε_c + ρβ` | `τ_0` | `τ_1` | `|τ_j|`, `j >= 2` | resto |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 5,215 | 1,2·10⁻⁴ | 5,130 | 0,257 | 4,35·10⁻⁵ | 1,7589 | ≤ 4·10⁻¹⁶ | 3,9·10⁻⁴ |

(`ε_c` é dominado pelo erro de montagem conservador `dB`, que inclui o resto da banda e o erro
de RT.) Rouché escalar contra `g(s) = τ_0 − (s − c)τ_1`: `3,9·10⁻⁴ < ρ|τ_1| − |τ_0| = 0,0879`.
**`K̂_Z` tem exatamente um autovalor no disco, simples, em ≈ 0,4010247311.** `t` é linear
porque `r` é quase o autovetor exato: `t(s) = −(s − s_K)/(l^H r)`.

**Circunferência, estimativa** (Z = 20×80, G+ = 40×160, F = 200×400; 16 tiles de raio
`rZ = 0,0098`; T1 e T2 em `c` com `R = 0,0598`; local, ~70 s por tile). T1: `α11` 0,452,
`‖Q_11‖` 0,633. T2: `α22` 0,618, `β` 0,368. Em todos os tiles, `‖V‖` 19,0–20,8,
`N_ZZ` 0,305–0,328, Perron do tile 0,816–0,818 e **Rouché `θ'` 0,810–0,812**. O `θ` quase não
varia na circunferência: quem manda são os níveis T, e o nível Z tem folga.

**Circunferência, RIGOROSA** (laboratório, 25/09, 23 min: T1 148 s, T2 760 s, ~38 s por
tile; `build/s3b/rouche_K_circ.{txt,json}`). T1 `α11` 0,4520; T2 `α22` 0,6176, `β` 0,3677;
`ρ` de cada tile ~9·10⁻⁴ (é `‖V‖·dH`); **os 16 tiles fecham, pior `θ'` = 0,812558**,
reverificado em racionais a partir do JSON. Com a contagem finita:

> **K tem exatamente uma raiz em `|s − 0,401| < 0,05`, e ela é simples** (em `X_κ1`, logo em
> `Y` pelo lema abaixo). Com S3a (K injetivo em `Y` no resto de D), **`s_K` é a única raiz de K
> em D.** O item 1 de S3b está fechado.

**Lema (regularidade: o espaço com sinal não cria raízes).** Autovetores e cadeias de Jordan
de K (e de L) no espaço com sinal `X_κ1` estão no espaço simétrico `Y`; raízes e
multiplicidades coincidem em `X_κ1` e em `Y`. *Prova.* `Y = X_κ1 ∩ X_{1/κ1}` e
`‖B‖_Y <= ‖B‖_{X_κ1}` (a conjugação troca `X_κ1` e `X_{1/κ1}` e `B` é real). Para `M` grande,
`‖P_{|m|>M} Q0 B‖ < 1` em `X_κ1` e em `Y` (`‖Q0 P_{|m|>M}‖ = O(1/M)`). Se
`x = −Q0 B x + f` com `f ∈ Y` (autovetor: `f = 0`; cadeia de Jordan: `f` do elo anterior),
`P_high x` é o ponto fixo, único em `X_κ1`, de
`y ↦ −P_high Q0 B (P_low x + y) + P_high f`. `P_low x` tem finitos `m`, e nelas as duas normas
são equivalentes, logo está em `Y`. O ponto fixo em `Y` também é ponto fixo em `X_κ1`, logo é
`P_high x`. ∎ Com isso, "K injetivo em `Y` fora de `D(c, ρ)`" (S3a, as duas metades por
conjugação) e "exatamente uma raiz em `D(c, ρ)` em `X_κ1`" (Rouché) falam das mesmas raízes.

**Para L, o NK basta (sem Rouché).** Se `DF(h, λ) = [[L(λ), h], [ℓ, 0]]` é invertível na solução
(o NK dá isso na bola inteira), `λ` é **algebricamente simples**: um segundo autovetor `h₂`
daria `DF(h₂ − ℓ(h₂)h, 0) = 0`, e uma cadeia de Jordan `L(λ)g = −h` daria
`DF(g − ℓ(g)h, 1) = 0`. Não é preciso contar raízes de L perto de `λ_L`: outra raiz `λ'` em
`D(c, ρ)` que violasse as constraints seria raiz de K, logo `λ' = s_K = λ_L`.

## 9. L em Z = 32×128: o NK não fecha com janelas; a perda é da agregação (25/09, madrugada)

Rodada completa de `nk_L.py` no laboratório (nível Z da rodada anterior, janelas W = 16,
2,7 h; `build/s3b_nk_32x128_janelas.txt`):

| `a` (Z'←T) | `b` (T←Z) | `α_T` (janelas) | far←T | `β` | `θ0` |
|---:|---:|---:|---:|---:|---:|
| 0,3735 | 0,3657 | 0,952 | 0,094 | 0,495 | **1,143** |

Com os demais fixos, `θ0` em função de `α_T`: 0,85 → 1,061; 0,78 → 1,007; 0,70 → 0,947;
0,65 → 0,911. É preciso `α_T <= 0,77`, e `<= 0,65` para ter margem.

**O `α_T` verdadeiro é 0,637** (`experiment_windows_loss_L.py`, potência com produtos implícitos
sobre todo T, F = 80×200 — a borda de Z, onde mora o singular, fica igual):

| cota | `α_T` | perda |
|---|---:|---:|
| verdadeira | 0,637 | — |
| janelas W = 16 | 0,917 | 1,44× |
| janelas W = 32 | 0,835 | 1,31× |

O maior bloco diagonal das janelas **sozinho** já dá 0,638: a perda é toda da soma dos
acoplamentos com as janelas vizinhas (~0,25), que a matriz de normas soma sem as fases. A
matriz de normas é tridiagonal: a banda efetiva do fundo é menor que 16 modos. O raio de Perron
dessa matriz é igual à norma (ela é quase simétrica): pesos por janela não ajudam.

**Onde mora o singular dominante:** 87% da massa em `32 <= |m| < 40`, `64 <= n < 128`, e 9% em
`128 <= n < 160` (entrada e saída). É uma região pequena, colada na borda de Z em `m`.

Duas saídas em teste:
1. **Deslocar a borda:** Z mais largo em `m` e mais curto em `n` (40×96, 40×112, 48×96), com a
   mesma dimensão. O `α` verdadeiro cai como ~20/M.
2. **Um segundo nível denso:** a região forte `T1 = {32 <= |m| < 48, n < 160}` (dois lados,
   desacoplados) ganha inverso aproximado próprio, `V_1 = (I + A_11)⁻¹`. O 0,637 sai da
   diagonal de Perron: sobram `‖V_1 A_12‖` e `‖A_21‖`, e as janelas só em `T2 = T \ T1`.
   É a arquitetura de S3a com dois blocos densos (V_b em Z', V_1 em T1).

**Resultado da saída 2** (`experiment_T1_L.py`, T1 = {32 <= |m| < 48, n < 160}, 4.480 por lado,
F = 80×200): `‖A_11‖` 0,637, mas **raio espectral 0,174** — `A` é muito não normal —, e
`‖V_1‖` 1,69. `‖V_1 A_12‖` 0,27, `‖A_21‖` 0,25; T2: verdadeiro 0,47, janelas 0,746 (perda
1,59×). Perron da parte T: 0,586 com `A_22` verdadeiro, **0,827 com janelas** (antes, 0,952).
Melhora, não basta. A não normalidade sugere o passo seguinte: **bloco-Jacobi nas janelas
fortes** (`V_J = (I + A_JJ)⁻¹` em cada janela de norma > 0,3): a diagonal de Perron zera nelas
e sobram acoplamentos ~0,25 pré-multiplicados por `V_J`.

**Resultado da saída 1** (laboratório): deslocar a borda **não ajuda**. `α_T` verdadeiro:
40×96 0,680; 40×112 0,606; 48×96 0,680 (contra 0,637 em 32×128); janelas perdem 1,54–1,56×.
Encurtar `n` expõe as caudas radiais dos modos de Z, que passam a limitar `α`. Descartada.

**Bloco-Jacobi nas janelas** (`experiment_blockjacobi_L.py`, W = 16, F = 80×200; densas as 8
janelas de norma > 0,3, `|m|` de 16 a 80): todas são muito não normais (norma 0,33–0,64, raio
espectral 0,12–0,16; `‖V_J‖` 1,27–1,69), e os acoplamentos pré-multiplicados quase não mudam
(0,262 contra 0,256). **Perron da parte T: 0,416** (sem precondicionar, 0,917).

Com Z e far (a = 0,3735, b = 0,3657, β = 0,495, far←T = 0,094), `θ0` estimado: 0,77 (sem
fatores), 0,83 (`V_J` no acoplamento com Z), 0,88 (`V_J` também no far — conservador demais:
as janelas densas ficam em `n < 200`, `|m| < 80`, e o far só as alcança por descidas ~κ2⁻²⁰⁰).
Raio NK previsto ~Y/(1 − θ0) ~ 10⁻⁵.

**Arquitetura do certificado de L:** níveis Z' (bordejado, `V_b`), D_1…D_8 (janelas fortes em
`n < 200`, `V_J` denso de ~5,6 mil), R (o resto de F, identidade, normas de blocos por janela)
e far. Linha de Z' com agregação ℓ² (`a·‖v_T‖_2`, pois `a` é a norma conjunta sobre T).


## 10. Bloco-Jacobi em F = 200×400: θ0 = 1,24 — a cauda radial e o far (25/09)

`nk_L2.py` no laboratório (Z = 32×128, F = 200×400, janelas W = 16 alinhadas com a borda de Z,
densas as 8 janelas com `16 <= |m| < 80` em `n < 200`; 85 min; `build/s3b/nk2_32x128.{txt,json}`):
**θ0 = 1,241** (linha de Z' com ℓ²), contra 0,77–0,88 previstos. O experimento de F = 80×200
cortava a cauda em `n < 200` e `|m| < 80`; aqui entram:

1. **As caudas radiais das janelas densas** (`n` de 200 a 400) e as janelas de fora, identidade
   como precondicionador: diagonal 0,24–0,37 (R 48..63: 0,374), vizinhas ~0,15.
2. **O far acoplado a todos os níveis R com β = 0,495 cada.** Na norma do máximo com 33 níveis, a
   linha do far soma ~0,03 de cada um e empurra `v_far` para ~1; isso volta a cada R como
   β·v_far/v_R ≈ 0,5–1,1. É a amplificação ~√(níveis) da norma do máximo quando um nível só
   se acopla a todos. A cota conjunta verdadeira é pequena (far←T = 0,094 na rodada de 3 níveis).

Tratando a cauda inteira como UM nível com norma ℓ² ponderada por dentro (T←T pela norma 2 da
matriz de normas ponderada; T←far e far←T conjuntos), o θ0 cai para 1,11–1,18, e o nível T
sozinho tem raio 0,74. Falta saber quanto disso é agregação e quanto é acoplamento de verdade:
`experiment_tail_true_L.py` mede a norma VERDADEIRA de `E_TT` com o bloco-Jacobi em F = 200×400
(potência implícita; Q0 em O(n) por modo, por substituição com soma acumulada).

**Resultado** (`build/s3b/cauda_verdadeira.txt`, laboratório, 40 min): em F = 200×400,
`‖P_T B Q0 P_T‖` verdadeiro = 0,6370 (igual a F = 80×200: o máximo mora na borda de Z) e
**`‖E_TT‖` verdadeiro com o bloco-Jacobi das 8 janelas = 0,388**. A cota por blocos dá 0,74: o
obstáculo é **só a agregação** (perda de 1,9×), não o acoplamento.

Outras medidas: a cauda radial `‖R_m⁻¹ P_{n>=400}‖` vale ~0,12 para todo `|m|` baixo (β ≈ 0,49
não depende de onde; é a descida radial com κ2 = 5/4); só Nc ~ 800 a reduziria à metade.

## 11. `nk_L3`: todas as janelas densas, cauda como um nível ℓ² (25/09)

Com a diagonal de blocos nula em TODAS as janelas (inversos densos das janelas inteiras,
`n < Nc`; a central dividida em duas), a agregação fica numa cadeia de diagonal nula: com os
números medidos, T←T ≈ 0,44–0,49 e θ0 ≈ 0,84–0,91 (modelo). Montagem em 3 níveis (Z', T, far)
com `‖x_T‖_W² = Σ w_J² ‖x_J‖²`: T←T = `‖W N W⁻¹‖₂`; T←far só pelas linhas altas de cada janela
(`n >= Nc − BAND_N − 60`; nas baixas, descida ~10⁻⁴), com `‖V_J P_alto‖`; far←T pela matriz de
blocos (grupos do far × janelas); Z'←T pela norma conjunta. Entradas densas por potência com a
LU de `Ĥ_J` (sem formar `V_J A_IJ`).

**Resultado de `nk_L3`** (laboratório, 2,3 h; `build/s3b/nk3_32x128.{txt,json}`; 26 janelas densas
de até 11.200, a central dividida): `‖V_J‖` 1,04–1,69, `‖V_J P_alto‖` 1,05–1,28, acoplamentos
vizinhos 0,04–0,25. Montagem de 3 níveis (pesos 1):

| T←T | T←Z' | Z'←T | T←far | far←T | β | **θ0** | raio NK |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0,396 (verdadeira: 0,388) | 0,559 | 0,373 | 0,636 | 0,094 | 0,495 | **0,828** | ~8·10⁻⁵ |

A cadeia de diagonal nula não perde quase nada (0,396 contra a norma verdadeira 0,388). **O
certificado NK de L fecha na estimativa**; segue o modo rigoroso: etapas A (`nk_L2_Z.py`, nível Z'),
B (`nk_L3_T.py`, cauda), C (`nk_L3_C.py`, montagem) e D (`nk_L3_D.py`, testemunho).

## 12. O certificado RIGOROSO de L e o fechamento de S3b (25/09)

Quatro etapas (`build/s3b/rig/`; Z = 32×128, F = 200×400, janelas W = 16):

**A — nível Z' bordejado** (`nk_L2_Z.py`, laboratório, 31 min). Inverso EXATO da matriz
calculada `M̂` (14.017), cotas por Loewner (Rump na realificação): `‖M̂⁻¹‖ <= 20,66748`,
`a <= 0,37889`, `ρ <= 4,8·10⁻⁵`, `Y_Z <= 2,47·10⁻⁵` (piso do erro publicado de RT),
Z'←far `<= 5,9·10⁻³`. O centro y0 é definido por `Q0 y0 = h0w` (`h0w.npy`).

**B — cauda com bloco-Jacobi em todas as janelas** (`nk_L3_T.py`; 26 janelas densas de até
11.200; PC gamer, 1,6 h + 16 min de far): por janela, `‖V_J‖ <= 1,692`, `‖V_J P_alto‖ <= 1,281`,
acoplamentos `<= 0,263`, diagonais `ρ_J <= 5,3·10⁻³`, T←Z' `<= 0,367`, resíduo `<= 6,0·10⁻⁶`;
far: `β = 0,4951`, descida `8,9·10⁻⁴`, far←janelas 0,041–0,063. Perturbação por bloco:
arredondamento + resto da banda (|d| > 16) + erro de RT + μ (`fundo_L.py`).

**C — montagem** (`nk_L3_C.py`, racionais exatos). Três níveis (Z', T em ℓ² ponderada, far),
com a perturbação GLOBAL (erro de RT e μ, que não respeita a banda) somada de novo em toda
entrada:

| T←T | T←Z' | Z'←T | T←far | far←T | β | **θ** | Z2 | Y | **r** |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0,3997 | 0,5739 | 0,3789 | 0,6398 | 0,0950 | 0,4951 | **0,839826** | 116,8 | 5,44·10⁻⁵ | **6,19·10⁻⁴** |

> **NK verificado:** existe um único zero `x* = (y*, λ*)` com `‖x* − x0‖ <= r`;
> **|λ* − 0,401024733| <= 2,81·10⁻⁴**; `DF(x*)` é invertível, logo **λ* é autovalor
> algebricamente simples de L** (§8).

**D — testemunho** (`nk_L3_D.py`; caixa F = 12×48 no espaço de K, pesos (129/128, 9/8)):
`c_ref = ‖Π_F C(λ0) h0‖ >= 0,17174`; perdas: nível Z' `2,82·10⁻³` (`‖Π_F C Q0 P_Z‖ <= 10,05`),
cauda `2,6·10⁻¹⁰` (as colunas de cauda começam em n = 128 e só descem até F por
(5/4)⁻⁸⁰), far `4·10⁻¹⁴`, λ `4,2·10⁻⁴`, perturbação global. **`‖Π_F C(λ*) h*‖ >= 0,16849 > 0`:
o autovetor de L em λ* VIOLA as constraints.**

**Fechamento de S3b.** Como `K C = M L`, λ* é raiz de K. Pelo Rouché de K (§8) e por S3a, K
tem uma única raiz em D, `s_K`, simples; λ* está no disco de raio 0,05 (a menos de 2,8·10⁻⁴ de
0,401), logo **λ* = s_K**. Qualquer outra raiz de L em D com autovetor que viole as constraints
seria raiz de K, logo igual a `s_K` — e em `s_K` L é algebricamente simples. Portanto:

> **Exatamente uma raiz de L em D (com multiplicidade) viola as constraints: `s_K ≈ 0,40102`,
> simples. Todas as outras raízes de L em D satisfazem as constraints.**

Condicionado, como todo o repositório, à bola publicada de RT (ε_ω = 2⁻²⁵ + 2⁻²⁷⁷) e ao μ
com 80 dígitos. Pendente: revisão independente (como a de S3a).

## 13. Depois da revisão independente (25/09, tarde)

`AUDITORIA_S3B_25SET.md`. O revisor (Fable, artefatos obrigatórios, sem acesso à narrativa)
não achou erro que derrube o certificado; reproduziu C e D bit a bit; 19 normas verdadeiras
abaixo das cotas. Achou cinco lacunas formais (L1–L5, corrigidas, sem efeito) e mostrou que a
margem do NK era mínima (discriminante 2,4·10⁻⁴). Com as normas RESTRITAS de Q0 no termo de 2ª
ordem (`nk_L3_q.py`: ‖P_Z Q0 P_T‖ = 0,169 em vez de 1,92), os números finais são:

| θ | Z2 | Y | r | \|λ* − 0,401024733\| | discriminante | testemunho |
|---:|---:|---:|---:|---:|---:|---:|
| 0,839936 | 42,96 | 5,44·10⁻⁵ | 3,78·10⁻⁴ | <= 1,71·10⁻⁴ | 0,0163 | >= 0,16976 |

A cota do símbolo `NORMA_BL` pode ser até 23% maior sem que o NK deixe de fechar. A ressalva
sobre cadeias de Jordan das outras raízes de L em D está respondida por
`CONSTRAINT_GAUGE_EQUIVALENCE.md` §1 (K(μ) invertível para μ ≠ s_K).
