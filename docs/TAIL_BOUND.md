# O bound de cauda `‖C − C_N‖ ≤ ε_N` (obrigação C1)

> Auditoria de 14/09/2026: a composição abaixo é condicional. P2 não decorre
> da analiticidade do fundo e falha já para `1+∂τ` na realização sem pesos.
> Além disso, a matriz de Galerkin da resolvente não é em geral a resolvente
> da matriz de Galerkin usada pelo produtor. A obrigação adicional C1-G e
> as demonstrações estão em [`ANALYTIC_AUDIT.md`](ANALYTIC_AUDIT.md).
> A alternativa de cauda algébrica com pesos fixos está em
> [`ALGEBRAIC_TAIL.md`](ALGEBRAIC_TAIL.md); a fórmula exponencial histórica
> abaixo não é usada nessa nova derivação.

O ledger [`PROOF_OBLIGATIONS.md`](PROOF_OBLIGATIONS.md) registra em C1 que
"faltam os três bounds PDE" sem nomeá-los. Este documento nomeia os três,
monta a cadeia quantitativa que os compõe em `ε_N`, diz **em que norma cada um
vive**, e separa o que existe na literatura de RT do que ainda não existe.

Nada aqui é uma demonstração de EDP. O que fica demonstrado é a **composição**:
dadas as três estimativas com constantes explícitas, a cauda do truncamento é
majorada pela fórmula final. A composição aritmética está verificada pelo
kernel do Lean em
[`formal/ChoptuikFormal/TailBound.lean`](../formal/ChoptuikFormal/TailBound.lean).

## 1. O objeto

Do lema de redução do ledger: `P(z) = P(z₀) + (z − z₀)J`, com `P(z₀)` bijetivo,
`R₀ = P(z₀)⁻¹` compacto e

```
C = J R₀ ,        C_N = truncamento de Galerkin Fourier×Chebyshev de C .
```

No pencil de Choptuik `z = s`, `P(s) = L_A(s) = L_A(0) + s·I`, de modo que
`J = I` — mas atenção, `I` aqui é a **inclusão do domínio no contradomínio**,
não a identidade de um espaço consigo mesmo (§3). Esse é o primeiro lugar onde
se troca de norma sem perceber.

O truncamento é o retângulo de índices

```
Π = Π_{N_F,N_C} :  mantém  |n| ≤ N_F  (Fourier em τ_R)  e  m ≤ N_C  (Chebyshev em ξ),
```

nas mesmas variáveis certificadas em que o código de RT já representa o fundo.
As malhas efetivamente calculadas no repositório são `N_F × N_C` = 4×12, 6×18,
8×24.

## 2. Duas escolhas de truncamento, e por que isso importa

| escolha | definição | erro a estimar | inversa em `I + λC_N` |
|---|---|---|---|
| unilateral | `C_N = Π C` | `(I − Π)C` | fórmula de Woodbury (§6) |
| bilateral | `C_N = Π C Π` | `(I−Π)C + ΠC(I−Π)` | `max(1, ‖(I+λM_N)⁻¹‖)` |

**Lema 2.1 (o determinante finito é o mesmo).** Seja `V_N = ran Π`,
`ι : V_N ↪ H` a inclusão e `B = ΠC : H → V_N`. Então `ΠC = ι B` e
`ΠCΠ = ι(BΠ)`. Pela identidade de Sylvester para operadores de posto finito
(`det(I+λAB) = det(I+λBA)` sempre que `AB` e `BA` são de traço),

```
det_H(I + λ ΠC) = det_{V_N}(I + λ B ι) = det_{V_N}(I + λ ΠCΠ|_{V_N}) = det_H(I + λ ΠCΠ) ,
```

e `B ι = ΠCΠ|_{V_N}` é a compressão da resolvente, que não deve ser
identificada com `(A_N+z0 I)^-1` sem a estimativa C1-G. Para essa compressão,
**o winding do determinante finito não distingue as duas
escolhas**; só o bound de erro e o bound da inversa as distinguem.

**Consequência.** A escolha unilateral consome uma hipótese analítica a menos
(§6), mas paga na norma da inversa; a §7.2 a adota. A escolha bilateral é a que
o `Tile.transferMajorant` de `RationalCertificate.lean` hoje pressupõe — a fórmula
`max(1, ‖V‖/(1−η))` só é um majorante porque `I + λΠCΠ` é diagonal por blocos
em `H = V_N ⊕ V_N^⊥` (age como identidade em `V_N^⊥`). Portanto, **o
`tailError` que aquele tile consome é o `ε_N` bilateral**, e o bilateral exige a
hipótese dual `P2*` da §6. Isso não está registrado no ledger e deveria estar.

## 3. Espaços e normas

Trocar de norma no meio da cadeia é o erro clássico que invalida esse tipo de
argumento. Fixamos quatro objetos e nunca mais mudamos:

| símbolo | espaço | norma | papel |
|---|---|---|---|
| `X` | domínio de `P(z₀)` | `‖·‖_X` (energia, com a regularidade em `∂_τ` do sistema selecionado) | onde vive a solução |
| `Y` | contradomínio de `P(z₀)` | `‖·‖_Y` | onde vive o dado; **é o espaço ambiente do argumento de Neumann/determinante**, pois `C : Y → Y` |
| `𝒜_{σ,ρ}` | classe analítica quantitativa | `‖u‖_{σ,ρ}` (§3.1) | onde a regularidade compra o decaimento |
| `V_N` | `ran Π` | herdada de `Y` | onde o determinante é finito |

Exigência de coerência, usada o tempo todo: `𝒜_{σ,ρ} ↪ Y` continuamente com
constante `1` na normalização adotada, e `Π` age diagonalmente nos coeficientes.

### 3.1 Duas realizações concretas

Escreva `u(τ_R, ξ) = Σ_{n∈ℤ} Σ_{m≥0} a_{n,m} e^{i n τ_R} T_m(ξ)`.

**(W) Realização de Wiener (ℓ¹ com pesos).**

```
‖u‖_𝒜 = Σ_{n,m} |a_{n,m}| ,        ‖u‖_{σ,ρ} = Σ_{n,m} |a_{n,m}| e^{σ|n|} ρ^m ,
```

com `Y = 𝒜 = 𝒜_{0,1}`. É uma álgebra de Banach para o produto pontual (os
coeficientes de Chebyshev satisfazem `T_j T_k = (T_{j+k} + T_{|j−k|})/2`, logo a
convolução tem norma `1`), o que importa porque `L_A(0)` contém multiplicação
pelo fundo. Nesta realização `‖Π‖ = 1` e `‖u‖_∞ ≤ ‖u‖_𝒜`.

**(S) Realização de norma do supremo (Cauchy + Bernstein).**

```
Y = C⁰ ,  ‖u‖_Y = sup |u| ;
𝒜_{σ,ρ} = {u holomorfa em {|Im τ| ≤ σ} × E_ρ} ,  ‖u‖_{σ,ρ} = sup_{faixa × E_ρ} |u| ,
```

onde `E_ρ` é a elipse de Bernstein de parâmetro `ρ > 1` (semi-eixos
`(ρ±ρ^{-1})/2`). É a realização mais próxima do enunciado clássico de
analiticidade, mas paga constantes maiores (§5.4) e tem uma armadilha: **em
`C⁰` a projeção de Fourier não tem norma `1`** (constante de Lebesgue `~ log N`).
Por isso (S) só serve para o truncamento unilateral, onde `‖Π‖` não aparece.

## 4. As três hipóteses de EDP

| ID | nome | enunciado | norma | constante |
|---|---|---|---|---|
| **P1** | invertibilidade coerciva no ponto-base | `P(z₀) : X → Y` é bijetivo e `‖R₀ f‖_X ≤ M₀ ‖f‖_Y` para todo `f ∈ Y` | `Y → X` | `M₀` |
| **P2** | regularidade analítica quantitativa | para todo `u ∈ ran R₀`: `u ∈ 𝒜_{σ,ρ}` e `‖u‖_{σ,ρ} ≤ A₀ ‖u‖_X`, com `σ > 0`, `ρ > 1` **e majorantes racionais explícitos** `q_F ≥ e^{−σ}`, `q_C ≥ 1/ρ`, `q_F, q_C < 1` | `X → 𝒜_{σ,ρ}` | `A₀, σ, ρ` |
| **P3** | limitação do coeficiente do pencil na escala analítica | `‖J v‖_{σ',ρ'} ≤ M_J ‖v‖_{σ,ρ}` com `0 < σ' ≤ σ`, `1 < ρ' ≤ ρ` | `𝒜_{σ,ρ} → 𝒜_{σ',ρ'}` | `M_J` (e a perda `σ→σ'`, `ρ→ρ'`) |

Três observações que não são decorativas:

1. **P2 é a priori sobre a solução, não sobre o dado.** É assim que estimativas
   de regularidade analítica realmente saem de EDP: limita-se a norma forte da
   solução pela norma fraca da própria solução. O fator `M₀` só entra depois,
   ao compor com P1 (§5.2). Se P2 fosse enunciada diretamente em termos do dado
   (`‖R₀f‖_{σ,ρ} ≤ A₀‖f‖_Y`), o `M₀` sumiria da fórmula final — e a fórmula
   continuaria correta, com outro significado de `A₀`. As duas versões não são
   intercambiáveis: a primeira é a que se demonstra.
2. **Fornecer `q_F` e `q_C` racionais faz parte de P2.** `e^{−σ}` e `1/ρ` são
   irracionais em geral, e o verificador exato não aceita "existe `σ > 0`". A
   obrigação analítica tem de terminar em um majorante racional explícito.
3. **`M_J = 1` não é de graça.** Como `L_A(s) = L_A(0) + s·I`, `J` é a inclusão
   `X ↪ Y`; `M_J = ‖ι_{X↪Y}‖` e vale `1` apenas se a norma de `X` dominar a de
   `Y` com constante `1`. Se `J` fosse um operador diferencial de primeira
   ordem com coeficientes analíticos — e é o que acontece em qualquer
   formulação em que a derivada `∂_τ` não venha com coeficiente identidade —
   então P3 traz perda estrita `σ' < σ` e `M_J` vira uma constante de Cauchy do
   tipo `C/(σ−σ')`. A forma afim `L_A(0) + s·I` é exatamente o que evita isso.

## 5. A derivação

### 5.1 Passo 0: reduzir ao lado do contradomínio

Com o truncamento unilateral, `C − C_N = (I − Π)C`. Nenhuma hipótese dual é
usada: só é preciso saber que **a imagem de `C` é analítica**. Pelo Lema 2.1,
essa escolha não altera o determinante finito.

### 5.2 Passo 1: compor P1 e P2

Para `f ∈ Y` com `‖f‖_Y ≤ 1`, seja `u = R₀f`. Por P1, `‖u‖_X ≤ M₀`. Por P2,

```
‖R₀ f‖_{σ,ρ} ≤ A₀ ‖R₀ f‖_X ≤ A₀ M₀ ‖f‖_Y .            (1)
```

### 5.3 Passo 2: aplicar P3

```
‖C f‖_{σ',ρ'} = ‖J R₀ f‖_{σ',ρ'} ≤ M_J ‖R₀ f‖_{σ,ρ} ≤ M_J A₀ M₀ ‖f‖_Y .   (2)
```

A partir daqui só se usa (2): `C` leva a bola unitária de `Y` na bola de raio
`M_J A₀ M₀` da classe analítica.

### 5.4 Passo 3: a cauda do truncamento

O complemento do retângulo de truncamento é coberto pela **união**
`{|n| > N_F} ∪ {m > N_C}`; a estimativa é subaditiva sobre essa união (não é
uma igualdade, e a interseção é contada duas vezes — a favor da desigualdade).

**Realização (W).** Para `v` com `‖v‖_{σ',ρ'} ≤ A`, escrevendo
`|a_{n,m}| = (|a_{n,m}| e^{σ'|n|} ρ'^m) · e^{−σ'|n|} ρ'^{−m}` e usando
`e^{−σ'|n|} ≤ q_F^{N_F+1}` quando `|n| > N_F` (com `ρ'^{−m} ≤ 1`), e
simetricamente na outra faixa:

```
‖(I − Π)v‖_𝒜 = Σ_{(n,m) ∉ caixa} |a_{n,m}| ≤ (q_F^{N_F+1} + q_C^{N_C+1}) ‖v‖_{σ',ρ'} .
```

Não sobra nenhum denominador geométrico: os pesos já somaram a cauda. Ou seja
`c_F = c_C = 1`.

**Realização (S).** Estimativa de Cauchy em `τ_R` (deslocar o contorno para
`Im τ = ∓σ'`) e desigualdade de Bernstein em `ξ` (para `f` holomorfa e limitada
por `A` em `E_ρ'`, os coeficientes de Chebyshev satisfazem `|a_m| ≤ 2A ρ'^{−m}`
para `m ≥ 1`, e `|a_0| ≤ A`) dão, no domínio produto,

```
|a_{n,m}| ≤ 2 A e^{−σ'|n|} ρ'^{−m} ≤ 2A q_F^{|n|} q_C^{m} .
```

Somando `|e^{inτ}| = 1`, `|T_m| ≤ 1` em `[−1,1]` e as duas séries geométricas
(`Σ_{n>N} q^n = q^{N+1}/(1−q)` e `Σ_{n∈ℤ} q^{|n|} = (1+q)/(1−q)`):

```
Σ_{|n|>N_F} |a_{n,m}| ≤ [ 4/((1−q_F)(1−q_C)) ] · A · q_F^{N_F+1} ,
Σ_{m>N_C}   |a_{n,m}| ≤ [ 2(1+q_F)/((1−q_F)(1−q_C)) ] · A · q_C^{N_C+1} ,
```

isto é `c_F = 4/((1−q_F)(1−q_C))` e `c_C = 2(1+q_F)/((1−q_F)(1−q_C))`.

### 5.5 Resultado

**Proposição 5.1 (composição).** Sob P1, P2 e P3 e com as convenções da §3,
para todos `N_F, N_C ≥ 0`,

```
‖C − Π_{N_F,N_C} C‖_{Y→Y}  ≤  ε_N(N_F,N_C)
  =  M_J · A₀ · M₀ · ( c_F · q_F^{N_F+1} + c_C · q_C^{N_C+1} ) ,
```

com `(c_F, c_C) = (1, 1)` na realização (W) e
`(c_F, c_C) = (4/((1−q_F)(1−q_C)), 2(1+q_F)/((1−q_F)(1−q_C)))` na realização (S).

Três honestidades sobre a forma:

* o expoente é `N+1`, não `N`: o truncamento guarda `|n| ≤ N_F`, logo o primeiro
  modo descartado é `|n| = N_F + 1`;
* **não há fator algébrico em `N`** em nenhuma das duas realizações. Fatores
  `N^k` apareceriam se P2 fosse trocada por regularidade de ordem finita
  (`H^k`), caso em que a cauda vira `N^{−k}` e toda a economia do certificado
  muda: com decaimento algébrico as malhas necessárias explodem;
* denominadores geométricos `1/(1−q)` aparecem em (S) e **não** em (W). A
  fórmula "bonita" sugerida pela intuição — `c·e^{−σN}` puro — é a de (W); a de
  (S) é pior por um fator `≈ 4/((1−q_F)(1−q_C))`, que para `q = 2/5` já é `11`.

## 6. A variante bilateral e a hipótese dual `P2*`

Com `C_N = ΠCΠ`,

```
‖C − ΠCΠ‖ ≤ ‖(I−Π)C‖ + ‖Π‖ · ‖C(I−Π)‖ .
```

O segundo termo **não segue de P2**. `‖C(I−Π)‖ = ‖(I−Π*)C*‖` exige a hipótese
dual

> **P2\*** (regularidade analítica do adjunto): `C*` leva a bola de `Y*` na
> classe analítica `𝒜_{σ*,ρ*}` com constante `A₀*`.

Em geral `P2*` não é consequência de `P2`: o operador adjunto de um sistema
hiperbólico com bordo característico tem domínio e condições de contorno
próprias. Além disso, em (S) o fator `‖Π‖` é a constante de Lebesgue (cresce
como `log N_F`), então a variante bilateral **só é admissível na realização
(W)**, onde `‖Π‖ = 1`.

Contrapartida, já registrada na §2 e decidida na §7.2: a bilateral dá a norma
da inversa de graça,

```
‖(I + λΠCΠ)⁻¹‖ = max(1, ‖(I + λM_N)⁻¹‖) ,
```

enquanto a unilateral exige Woodbury:
`(I + λιB)⁻¹ = I − λι(I + λBι)⁻¹B`, logo

```
‖(I + λΠC)⁻¹‖ ≤ 1 + |λ| · ‖(I + λM_N)⁻¹‖ · ‖C‖ ,   ‖C‖ ≤ M_J M₀ ,
```

que é estritamente pior. Resumo da escolha: **unilateral compra uma hipótese
analítica a menos ao preço de um majorante de inversa pior; bilateral faz o
contrário.** A infraestrutura existente em `RationalCertificate.lean` está escrita para a
bilateral; a §7.2 decide trocá-la, e
mede o preço em quatro modos de Fourier.

## 7. Duas decisões que a integração tem de respeitar

As duas primeiras entradas desta seção eram, na versão anterior deste
documento, um menu de opções. Estão decididas, com a razão de cada escolha.

### 7.1 N-coh — a cadeia inteira vive em `ℓ¹`, norma por **coluna**

**Decisão: realização (W), e o certificado passa a usar `RatMatrix.oneNorm`
(máxima soma de módulos por coluna) em todos os pontos onde hoje usa
`RatMatrix.infNorm`.** `RationalCertificate.lean` já expõe as duas, com
`oneNorm m = m.transpose.infNorm`; a troca é de uma linha em três lugares
(`inverseDefect` usa `infNorm` duas vezes, `inverseBound` uma).

Por que `ℓ¹` e não `ℓ^∞`, sendo que a estimativa de cauda fecha nos dois (em
`ℓ^∞` ela até fica marginalmente melhor, com `max` no lugar da soma):

1. **`ℓ¹` nos coeficientes controla a função; `ℓ^∞` não controla nada.**
   `‖u‖_∞ ≤ Σ|a_{n,m}| = ‖u‖_𝒜`. Com o supremo dos coeficientes não se conclui
   sequer que `u` é limitada. Como o certificado precisa avaliar o fundo no ball
   de RT e majorar produtos pontuais, isto sozinho decide.
2. **`ℓ¹` é álgebra de Banach com constante `1`** para o produto pontual
   (`T_jT_k = (T_{j+k}+T_{|j-k|})/2` dá convolução de norma `1`). `L_A(0)`
   contém multiplicação pelo fundo e o termo bilinear `Γ₂`; em `ℓ^∞` a
   multiplicação só é limitada **pagando** a norma `ℓ¹` do outro fator, o que
   reintroduz a norma `ℓ¹` pela porta dos fundos.
3. `‖Π‖ = 1` vale nas duas, e é o que o argumento de Woodbury da §6 usa.

O que **muda** na derivação: nada. As constantes de cauda continuam
`c_F = c_C = 1` (§5.4, realização (W)), porque a derivação já foi feita em `ℓ¹`
com pesos. O que muda é do lado do código: toda norma de matriz consumida pelo
tile tem de ser soma por **coluna**. Manter `infNorm` significaria refazer a §5
em `ℓ^∞` e perder os itens 1 e 2 acima.

### 7.2 P2\* — adotar o majorante **unilateral** e não carregar a hipótese

**Decisão: truncamento unilateral `C_N = Π C`. `P2*` não entra no ledger como
hipótese a demonstrar; entra como hipótese evitada.**

O preço está medido. Com o unilateral, o majorante de transferência é

```
q_N  =  sup_Γ |z-z₀| · ε_N · ‖(I + (z-z₀)Π C)⁻¹‖
     ≤  d · ε_N · ( 1 + d · β_N · ‖C‖ ) ,
     d = sup_Γ |z-z₀| ,  β_N ≥ ‖(I + (z-z₀)M_N)⁻¹‖ ,  ‖C‖ ≤ M_J · M₀ ,
```

contra o bilateral `d · ε_N · max(1, β_N)`. Com os valores de calibração do
projeto (`d = 3`, `β_N = 50`, `‖C‖ = M_J·M₀ = 10`), a razão é
`(1 + d β_N ‖C‖)/max(1,β_N) = 1501/50 = 30,02`: o unilateral exige `ε_N` trinta
vezes menor. Traduzido em malha, sob `H-opt` e `N_C = 24`, isso é
`N_F : 9 → 13`. **Quatro modos de Fourier é o preço inteiro de não precisar
demonstrar regularidade analítica do adjunto de um sistema hiperbólico com
bordo característico** — cujo próprio adjunto ninguém escreveu ainda. A troca é
boa e não é próxima.

Forma exata a implementar no lugar de `Tile.transferMajorant`, com um campo
novo `operatorBound` (que é `M_J · M₀`, já presente em `TailData`):

```
transferMajorant = maxDistance * tailError *
                     (1 + maxDistance * inverseBound * operatorBound)
```

A versão bilateral atual, `maxDistance * tailError * max(1, inverseBound)`,
**não** é majorante do erro unilateral e só se aplica a `Π C Π`. As duas estão
implementadas lado a lado, com as regressões que verificam a razão `1501/50` e o
salto `9 → 13`, em
[`formal/ChoptuikFormal/DeterminantBound.lean`](../formal/ChoptuikFormal/DeterminantBound.lean)
(`TransferData.unilateralMajorant`, `TransferData.bilateralMajorant`).

### 7.3 O que permanece como obrigação, sem escolha a fazer

* **Uniformidade em `z`.** `ε_N` acima é para o `C` fixo do ponto-base; o lema
  de transferência precisa dele uniformemente sobre todo o contorno `Γ`, o que
  é automático aqui porque `C` não depende de `z` (a dependência está toda no
  fator `(z - z₀)`), mas deixa de ser automático se o ponto-base for movido de
  tile em tile.
* **Uniformidade no ball de RT.** `M₀, A₀, M_J, σ, ρ` têm de valer para **todo**
  fundo no ball certificado de RT, não para um representante numérico.
* **`matrixVariation` não é obrigação de EDP.** Para o pencil afim vale
  `F_N(z) - F_N(z_j) = (z-z_j)C_N` exatamente, logo
  `matrixVariation = raio_do_tile · ‖C_N‖` é aritmética sobre a matriz racional.
  Ver [`WINDING_DATA.md`](WINDING_DATA.md) §3.6.

## 8. De onde cada hipótese teria de vir, e o que falta hoje

| ID | o que já existe na literatura de RT | o que falta |
|---|---|---|
| **P1** | RT provam invertibilidade dos sistemas linearizados usados no esquema de Newton, com certificados `INV` verificados exatamente (auditados em parte neste repositório). Reiterer, arXiv:1510.05310, dá o quadro de estabilidade de codimensão finita para equações hiperbólicas periódicas. | Os certificados `INV` são para o operador **no ponto fixo**, nas normas de RT, sem deslocamento espectral. É preciso um `z₀` fora do espectro (tipicamente `Re s` grande, via estimativa de energia — obrigação F3) e a constante `M₀` nessa norma. Nada disso está publicado como número. |
| **P2** | Teorema 1 de RT dá a solução **real-analítica** no domínio `M = {u₊<0, u₊<u₋<−u₊/81}`. Analiticidade é qualitativa no enunciado. | (i) Extrair `σ` e `ρ` explícitos: a faixa de analiticidade em `τ_R` e a maior elipse de Bernstein contida na imagem em `ξ` do domínio de holomorfia. O bordo `u₋ = −u₊/81` é o que limita `ρ`: qualquer `ρ` rigoroso tem de sair dessa razão, e o fundo **não** é analítico até o ponto singular `(0,0)` nem até o horizonte de Cauchy. (ii) Uma estimativa a priori analítica para a **solução do problema linearizado**, não para o fundo. Para sistemas hiperbólicos a ferramenta natural é uma escala de espaços de Banach à Ovsyannikov, que perde raio de analiticidade proporcionalmente ao tempo; aqui `τ_R` é periódico, de modo que a perda tem de ser recuperada pela periodicidade — é exatamente esse ponto fixo que ninguém escreveu. (iii) Uniformidade sobre o ball certificado, e uniformidade até o bordo característico `ξ = ±1` (cone passado e centro), onde o sistema degenera. (iv) Os majorantes **racionais** `q_F, q_C`. |
| **P3** | A estrutura afim `L_A(s) = L_A(0) + s·I` está derivada na §2 de [`SPECTRAL_PROBLEM.md`](SPECTRAL_PROBLEM.md) e implementada em `patches/rt-spectral.patch`. | A constante `M_J = ‖ι_{X↪Y}‖` depende da definição precisa de `X` e `Y`, que o repositório ainda não fixou (obrigação F1: "sistema simétrico hiperbólico e domínio fechado"). Sem `X` e `Y` fixados, `M_J` não tem valor, e P1/P2 também não. |

Em uma frase: **P1 e P3 são difíceis mas estão no caminho já aberto por RT;
P2 é de outra natureza** — é uma estimativa analítica quantitativa uniforme
para o linearizado de um sistema hiperbólico degenerado no bordo, e nem RT nem
Reiterer publicaram algo nessa forma.

## 9. O que a aritmética exata diz sobre as malhas de hoje

As constantes abaixo **não são certificadas**. São dois conjuntos batizados
para poderem ser citados sem se disfarçarem de teorema:

* **H-opt**: `M₀ = 10`, `A₀ = 4`, `M_J = 1`, `q_F = q_C = 2/5`, realização (W).
  O valor `2/5` foi escolhido para ser compatível com a razão de convergência
  observada no localizador em ponto flutuante (as diferenças sucessivas do
  candidato `s` caem por um fator `≈ 0,09` a cada `ΔN_F = 2`, e o resíduo de
  constraints por `≈ 0,18`). **Isso é diagnóstico, não é bound.**
* **H-pes**: `M₀ = 10`, `A₀ = 40`, `M_J = 1`, `q_F = 3/5`, `q_C = 1/2`,
  realização (W).

O alvo `1/150` vem de `sup_Γ |z − z₀| = 3` e `sup_Γ ‖(I+(z−z₀)C_N)⁻¹‖ = 50` no
majorante de transferência do ledger; também é hipotético.

| malha `N_F×N_C` | dimensão real | `ε_N` sob H-opt | passa em `1/150`? | `ε_N` sob H-opt, realização (S) |
|---|---:|---:|---|---:|
| 4×12 | 138 | `100065536/244140625 ≈ 4,10·10⁻¹` | não (61× acima) | `≈ 4,55` |
| 6×18 | 333 | `≈ 6,55·10⁻²` | não (9,8× acima) | `≈ 7,28·10⁻¹` |
| 8×24 | 612 | `≈ 1,05·10⁻²` | **não (1,57× acima)** | `≈ 1,17·10⁻¹` |
| 9×24 | — | `≈ 4,19·10⁻³` | sim | — |
| 12×36 | — | `≈ 2,68·10⁻⁴` | sim | `≈ 2,98·10⁻³` (sim) |

Sob **H-pes** nem 20×32 passa (`≈ 8,8·10⁻³`); seria preciso algo como 24×48.

Leitura honesta: mesmo com constantes escolhidas de forma favorável, **a melhor
malha já calculada fica fora por um fator da ordem de 1,6**, e basta trocar
`q_F` de `2/5` para `3/5` — uma mudança modesta na faixa de analiticidade —
para que o requisito salte de `N_F = 9` para `N_F ≈ 24`, isto é, uma matriz
várias vezes maior. **Toda a viabilidade do certificado depende do valor
numérico de `σ`, que hoje ninguém sabe.** Essa é a conclusão operacional mais
importante deste documento: a prioridade da obrigação C1 não é implementação, é
extrair `σ` e `ρ` do argumento de analiticidade de RT.

## 10. O que está formalizado

[`formal/ChoptuikFormal/TailBound.lean`](../formal/ChoptuikFormal/TailBound.lean)
(`namespace Choptuik.TailBound`) implementa a composição da Proposição 5.1 em
`Rat` exato:

* `TailData` guarda `M₀, A₀, M_J, q_F, q_C, c_F, c_C, N_F, N_C`;
* `TailData.epsilon` calcula `ε_N` exatamente, com potências racionais
  `^ (n : Nat)`;
* `TailData.admissible` recusa dados degenerados (`q ≥ 1`, constantes negativas);
* `TailCertificate.valid_sound` e `epsilon_le_target_of_valid` são teoremas de
  solidez no estilo `of_decide_eq_true` já usado em `RationalCertificate.lean`;
* todas as regressões fecham por **`decide`** (kernel), nenhuma por
  `native_decide` — ver [`KERNEL_TRUST.md`](KERNEL_TRUST.md). A aritmética `Rat`
  reduz no kernel depois de `unseal Rat.mul Rat.add Rat.sub Rat.inv`; as quatro
  são `@[irreducible]` no Lean 4.33, e `unseal` só relaxa transparência, sem
  acrescentar nenhuma hipótese.

[`formal/ChoptuikFormal/DeterminantBound.lean`](../formal/ChoptuikFormal/DeterminantBound.lean)
acrescenta, no mesmo estilo, o majorante de transferência unilateral da §7.2, a
estimativa de aresta por Jacobi de [`WINDING_DATA.md`](WINDING_DATA.md) §3, e
`matrixVariation` afim — com as regressões que medem o preço da decisão 7.2.

O que **não** está formalizado, e não pode ser dado como feito: P1, P2, P3,
a identidade de Sylvester do Lema 2.1, as estimativas de Cauchy e de
Bernstein da §5.4, e a própria afirmação de que `ε_N` limita `‖C − C_N‖` — que é
justamente o enunciado da Proposição 5.1, cuja demonstração vive neste
documento e não no Lean.
