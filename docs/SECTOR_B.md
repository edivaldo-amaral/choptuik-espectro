# Setor B por redução ao setor A (obrigações S2/C4)

> Atualização de 14/09/2026: H-conv foi conferida no artigo RT (§1.1,
> `u_sigma=-(1+sigma*xi) exp(-mu*tau)`, `Theta: tau→tau+2pi`), e a
> seleção H-dof foi inspecionada no código original; veja §6.3 abaixo.
> A redução continua condicionada aos domínios funcionais e à equivalência
> física. Ela não resolve as lacunas C1-G/P2 da auditoria analítica.
> O contorno reduzido 4×12 foi agora calculado: winding **0**, com
> 1946 arestas originais e 55 após união exata de segmentos colineares.
> O certificado das caixas foi aceito pelo kernel Lean. Reproduza com
> `make verify-reduced-b`. Isso fecha a regressão FINITA, não C4 infinito.

[`SPECTRAL_PROBLEM.md`](SPECTRAL_PROBLEM.md) §3 registra a conjugação
`M⁻¹ L_B(s) M = L_A(s + i/2)` com `M = ` multiplicação por `exp(i τ_R/2)`, e o
ledger marca S2 como "redução exata, implementação nativa aberta" e C4 como
dependente de "fornecer as arestas do setor B". Este documento (i) reverifica a
conjugação do zero, (ii) enuncia e demonstra o lema de redução que ela implica,
(iii) lista **todas** as hipóteses que a redução consome, incluindo as que não
são automáticas, e (iv) traduz C4 em um retângulo concreto do plano `s` no setor
A — eliminando a necessidade de implementar `B` nativamente.

**Veredito antecipado.** A conjugação está correta e a redução fecha *no nível
espectral*, sob quatro hipóteses nomeadas na §6, das quais uma (H-dof) não pode
ser confirmada sem ler o código-fonte de RT e outra (H-afim) já é pressuposta
pela própria derivação do pencil. A reação da redução com a realidade **quase
quebra**: a escolha ingênua de janela simétrica coloca exatamente os modos reais
do setor B em cima do contorno. A §7 mostra onde isso acontece e qual janela
evita o problema.

## 1. Normalização e a involução

Sejam as quatro componentes `h = (h₁,h₂,h₃,h₄)` da perturbação nas variáveis de
RT, funções de `(τ_R, ξ)`. A involução "meio eco + troca de sinal de `φ`" é

```
(T h)(τ_R, ξ) = σ · h(τ_R + 2π, ξ) ,        σ = diag(1,1,1,−1) ,
```

de modo que

```
setor A = ker(T − 1) :  h₁,h₂,h₃  2π-periódicas,  h₄ antiperiódica ;
setor B = ker(T + 1) :  paridades trocadas ;
T² h = h(· + 4π) = h   para perturbações de período completo (4π).
```

Isto é exatamente o §3 de `SPECTRAL_PROBLEM.md`: `Θ` (meio eco) é o
deslocamento de `2π` em `τ_R` junto com `φ ↦ −φ`, e o **eco completo é `4π`**.

**Conferência interna da normalização.** A normalização não é uma convenção
livre: ela está fixada pelo próprio deslocamento `i/2` documentado. Se o meio
eco fosse `τ_R ↦ τ_R + π`, o conjugador teria de satisfazer
`m(τ+π) = −m(τ)`, isto é `m = e^{iτ_R}`, e o deslocamento seria `i`, não `i/2`.
Com meio eco `= 2π`, `m(τ+2π) = −m(τ)` dá `m = e^{iτ_R/2}` e deslocamento
`i/2`. O `i/2` de `SPECTRAL_PROBLEM.md` §3 e o `s_B = s_A − i/2` impresso por
`scripts/analyze_spectrum.py` só são consistentes com a leitura `2π/4π`. Fica
registrado como hipótese **H-conv** (§6).

## 2. Verificação da conjugação

**Lema 2.1 (derivada).** Para toda `h` diferenciável em `τ`,

```
d/dτ ( e^{iτ/2} h ) = e^{iτ/2} ( d/dτ + i/2 ) h .
```

Regra do produto; nada mais. Equivalentemente `M⁻¹ ∂_τ M = ∂_τ + i/2`.

**Lema 2.2 (`M` troca os setores).** Com `M h := e^{iτ_R/2} h`,

```
T(M h) = σ e^{i(τ+2π)/2} h(τ+2π) = e^{iπ} e^{iτ/2} σ h(τ+2π) = − M(T h) ,
```

isto é `T M = − M T`. Logo `M` leva `ker(T−1)` bijetivamente em `ker(T+1)`:
`M : A → B` é isomorfismo linear, com inverso `M⁻¹ = ` multiplicação por
`e^{−iτ_R/2}`. As paridades declaradas no §3 de `SPECTRAL_PROBLEM.md` são
exatamente as de `ker(T∓1)` acima, e o cálculo usa apenas
`e^{i(τ+2π)/2} = −e^{iτ/2}` — vale igualmente para as quatro componentes, não há
componente privilegiada.

**Lema 2.3 (conjugação do pencil).** Suponha (hipótese **H-afim**) que o sistema
selecionado se escreve com a derivada `∂_τ` de coeficiente identidade, de modo
que a substituição `∂_τ ↦ ∂_τ + s` produz exatamente `L(s) = L(0) + s·I`
(é o quadro boxed do §2 de `SPECTRAL_PROBLEM.md`), e que `L` comuta com `T`
(equivariância, que é o que sustenta a própria decomposição em setores). Então,
sobre `A`,

```
M⁻¹ L_B(0) M = L_A(0) + (i/2)·I = L_A(i/2) ,
M⁻¹ L_B(s) M = M⁻¹ (L_B(0) + s·I) M = L_A(0) + (s + i/2)·I = L_A(s + i/2) .
```

A primeira linha é o Lema 2.1 aplicado à única ocorrência de `∂_τ`; a segunda
usa que `M⁻¹ (s·I) M = s·I` porque `M` é multiplicação e não depende de `s`.
Isto reproduz a fórmula do §3 de `SPECTRAL_PROBLEM.md`. ∎

**Conferência do sinal contra o código.** `analyze_spectrum.py --show-b` imprime
`s_B = s_A − i/2`, que é a leitura correta de `M⁻¹L_B(s)M = L_A(s+i/2)`:
`L_B(s)` é singular exatamente quando `L_A(s + i/2)` é, isto é
`s_B + i/2 = s_A`. Coerente.

Observe-se, porém, que **o sinal só é definido módulo `i`**: `M' = e^{−iτ_R/2}`
também satisfaz `T M' = −M' T` e dá `M'⁻¹ L_B(s) M' = L_A(s − i/2)`. As duas
afirmações são compatíveis porque `Σ_A` é invariante por `s ↦ s + i` (§4), de
modo que `Σ_A − i/2 = Σ_A + i/2`. Quem tentar usar o sinal do deslocamento para
distinguir alguma coisa estará usando informação que não existe.

## 3. Transferência de multiplicidades

`M` **não depende de `s`**. Portanto a família analítica `s ↦ M⁻¹ L_B(s) M` é
uma conjugação por um isomorfismo constante, e toda a estrutura de Fredholm
analítica é transportada sem correção: núcleos, co-núcleos, cadeias de Jordan do
pencil, multiplicidade algébrica no sentido de Gohberg–Sigal e o próprio winding
do determinante. Em particular, como divisores em `ℂ` (conjuntos com
multiplicidade),

```
Σ_B = Σ_A − i/2 .                                     (3.1)
```

Mais geralmente, uma equivalência por famílias holomorfas de isomorfismos
com inversas holomorfas também preserva a multiplicidade algébrica local.
A independência de `s` simplifica a demonstração aqui; não é condição
necessária. Uma mudança singular no parâmetro, por outro lado, exige cuidado.

## 4. Os dois reticulados de Floquet

Multiplicar `h` por `e^{i k τ_R}` preserva as paridades de `T` se e só se
`e^{2πik} = 1`, isto é `k ∈ ℤ`. Logo

```
Σ_A + i = Σ_A ,       Σ_B + i = Σ_B .                  (4.1)
```

Já no espaço de período completo (`4π`-periódico, sem decompor), multiplicar por
`e^{i k τ_R/2}` com `k ∈ ℤ` preserva a periodicidade, e portanto

```
Σ_A ∪ Σ_B  é invariante por  s ↦ s + i/2 .             (4.2)
```

(4.1) e (4.2) são consistentes com (3.1): `Σ_A − i/2 − i/2 = Σ_A − i = Σ_A`.

**Duas larguras, nenhuma contradição.** O §3 de `SPECTRAL_PROBLEM.md` fala em
"faixa fundamental de largura `1/2` em `Im(s)`" e este documento vai falar em
largura `1`. São reticulados diferentes:

* o problema de **período completo** tem reticulado `(i/2)ℤ` — domínio
  fundamental de largura `1/2`, e nele contam-se os **dois** setores;
* o setor **A sozinho** tem reticulado `iℤ` — domínio fundamental de largura
  `1`, e nele conta-se apenas `A`.

Os dois procedimentos produzem o mesmo número (Lema 5.1). É por isso que uma
única computação em `A` sobre uma faixa de largura `1` substitui duas
computações sobre faixas de largura `1/2`.

## 5. O lema de redução

Para uma região `R ⊆ ℂ` escreva `N_A(R)` e `N_B(R)` para o número de pontos de
`Σ_A` e `Σ_B` em `R`, contados com multiplicidade algébrica.

**Lema 5.1 (redução).** Seja `Q ⊂ ℂ` um domínio fundamental semiaberto do
reticulado `(i/2)ℤ` (largura `1/2` em `Im s`), e seja
`Q̃ = Q ∪ (Q + i/2)`, um domínio fundamental de `iℤ` (largura `1`). Então, para
qualquer região `Ω` invariante por translações imaginárias (por exemplo
`Ω = {Re s > 0}`):

```
N_A(Q ∩ Ω) + N_B(Q ∩ Ω)  =  N_A(Q̃ ∩ Ω) .
```

*Demonstração.* Por (3.1), `Σ_B ∩ Q = (Σ_A − i/2) ∩ Q`, e `x ∈ (Σ_A − i/2) ∩ Q`
se e só se `x + i/2 ∈ Σ_A ∩ (Q + i/2)`. A translação por `i/2` é uma bijeção que
preserva multiplicidades (§3) e preserva `Ω`. Logo
`N_B(Q ∩ Ω) = N_A((Q + i/2) ∩ Ω)`. Como `Q` é semiaberto de largura exatamente
`1/2`, `Q` e `Q + i/2` são disjuntos e `Q̃ = Q ⊔ (Q+i/2)`; somando, obtém-se a
identidade. ∎

**Corolário 5.2 (contagem física).** A multiplicidade instável física por classe
de Floquet do problema de período completo é

```
N_phys = N_A(Q̃ ∩ {Re s > 0}) ,
```

com `Q̃` qualquer faixa semiaberta de largura `1` em `Im(s)`. Em particular,
nenhum objeto do setor `B` precisa ser construído.

**Corolário 5.3 (C3 e C4 no mesmo operador).** Escolhendo
`Q = {−1/4 < Im s ≤ 1/4}` e portanto `Q̃ = {−1/4 < Im s ≤ 3/4}`:

```
C3  (winding 1 no setor A) :  N_A( {−1/4 < Im s ≤ 1/4} ∩ {Re s > 0} ) = 1 ,
C4  (winding 0 no setor B) :  N_A( { 1/4 < Im s ≤ 3/4} ∩ {Re s > 0} ) = 0 .
```

C4 deixa de exigir implementação nativa de `B`: vira o winding de
`det[I + (z−z₀)C_N]` **do mesmo operador `L_A`** sobre um segundo retângulo,
transladado de `i/2`.

## 6. Hipóteses consumidas

Nenhuma das quatro é automática; três delas o repositório ainda não pode
confirmar sozinho.

| ID | hipótese | estado |
|---|---|---|
| **H-conv** | o meio eco é `τ_R ↦ τ_R + 2π` com `σ = diag(1,1,1,−1)`, e o eco completo é `4π` na normalização usada pelo código | consistente com `SPECTRAL_PROBLEM.md` §3 e com o `i/2` impresso pelo diagnóstico; não reconferida contra o fonte de RT |
| **H-afim** | `L(s) = L(0) + s·I` com coeficiente identidade em `∂_τ`, e `L` comuta com `T` | já é a hipótese que sustenta o pencil (S1 do ledger); a redução não acrescenta nada aqui |
| **H-iso** | `M` é isomorfismo dos espaços funcionais relevantes: preserva domínio, regularidade e condições de contorno | argumentada abaixo; a parte sobre a seleção de graus de liberdade (**H-dof**) é a que falta |
| **H-constr** | a estrutura de constraints é levada em constraints: `M⁻¹ C_B(s) M = C_A(s + i/2)` | argumentada abaixo; redutível a `H-afim` + equivariância do sistema `♯` |

### 6.1 Por que `H-iso` é quase automática, e onde não é

`M` é multiplicação por uma função unimodular **só de `τ_R`**. Logo:

* preserva regularidade em `ξ` e todas as condições de contorno em `ξ = ±1`
  (cone passado e centro), que não dependem de `τ_R`;
* preserva o domínio do operador: `∂_τ M = M(∂_τ + i/2)` mostra que `M` leva o
  domínio de `∂_τ` nele mesmo, e a perturbação `i/2` é limitada;
* é isometria nas normas `L²` e nas normas de Wiener sem peso; nas normas com
  peso de analiticidade `e^{σ|n|}` ela **não** é isometria — desloca o índice de
  Fourier de `1/2` e multiplica a norma por no máximo `e^{σ/2}`. Isso é
  limitado com inverso limitado, portanto inofensivo para contagem, mas **as
  constantes do certificado de cauda (`docs/TAIL_BOUND.md`) mudam por esse
  fator** se alguém tentar usar as mesmas constantes nos dois setores. Como o
  Corolário 5.2 nunca sai do setor `A`, esse fator não chega a aparecer — mais
  um motivo para fazer a redução em vez de implementar `B`.

**H-dof (a parte que falta).** O código de RT seleciona graus de liberdade por
`IndexDOF_build_VBold_alloc`, e a extensão em `patches/rt-spectral.patch` exige
`sec.off_m == 0 && sec.off_n == 0`. Se essa seleção impuser **alguma condição
dependente do índice de Fourier `m`** — tipicamente um tratamento especial do
modo `m = 0`, que existe em `A` e não tem análogo em `B` (harmônicos
semi-inteiros nunca se anulam) — então o espaço que o código realiza não é
exatamente `ker(T−1)` e `M` não o leva no `ker(T+1)` natural. A redução
continuaria valendo para os espaços abstratos, mas a matriz calculada deixaria
de representá-los. Confirmar isto exige ler `IndexDOF_build_VBold_alloc` no
artefato de RT, o que este repositório só faz em `make rt-audit`.

### 6.2 Por que `H-constr` se reduz ao que já se assume

O resíduo de constraints do §4 de `SPECTRAL_PROBLEM.md` é
`C_A(s)h = S^♯{Ξ H_{μ,s}h + μΓ₁h + 2ΞΓ₂(ω_*,h)}`, e a extensão calcula-o
exatamente como **afim em `s`**: `patches/rt-spectral.patch` escreve um bloco de
grau `0` e um bloco de grau `1`, sendo o de grau `1` a derivada em `s` do
resíduo `♯` (que é `−Ξ h₂`, só na segunda linha `♯`). Ou seja `C(s) = C(0) + sD`
com `D = ∂C/∂s = ` coeficiente de `∂_τ` no resíduo `♯`. Mas essa é exatamente a
situação do Lema 2.3: a substituição `∂_τ ↦ ∂_τ + i/2` produz
`M⁻¹C_B(s)M = C_A(s + i/2)` pela mesma conta de uma linha. Consequência: **um
autovetor do setor B satisfaz as constraints se e só se seu transportado em A
satisfaz**, e a contagem "com constraints" transfere-se junto com a contagem
bruta.

Dois detalhes conferidos no patch, porque são onde uma conjugação costuma
quebrar:

* a única operação não trivial nas linhas `♯` é um sinal
  `(d == 1 && n par) ? −1 : +1`, isto é um operador de **paridade no índice de
  Chebyshev** (`P` aplicado à segunda linha). `M` age só em `τ_R`; logo `M` e
  `P` comutam trivialmente;
* a seleção de linhas `♯` é uma seleção de **componentes** (`d = 0` e `d = 2`),
  que comuta com `T = σ ∘ shift` porque `σ` é diagonal nas componentes.

O que resta como hipótese genuína é a equivariância `T`-por-`T` do próprio
sistema `♯` de RT — a mesma que já se usa para afirmar que `L` preserva os dois
setores. A obrigação S3 do ledger (propagação das constraints para `s` complexo)
**não** é resolvida por nada disto: ela continua aberta, e agora aberta uma vez
só, para `A`, em vez de duas.

### 6.3 Conferência do fonte em 14/09/2026

`IndexDOF.c`, função `IndexDOF_build_VBold_alloc`, chama a seleção `0`,
que preserva `m` par nas primeiras três componentes e `m` ímpar na quarta,
e exige `n` ímpar na primeira componente. Não aplica `VGauged` (isso ocorre
somente na seleção `1`). `Field_set_VBasic2` elimina a parte imaginária de
`m=0`, exatamente a redundância da representação de funções reais.

Essa condição em `m=0` não elimina o modo constante: sua coordenada real
permanece. Ao complexificar o espaço real de graus de liberdade, essa
coordenada também recebe amplitudes complexas. Portanto não se deve
interpretar a ausência do DOF imaginário redundante como uma condição
extra que impediria a conjugação. O deslocamento `exp(i*tau/2)` atua na
complexificação da série completa, não como simples deslocamento da tabela
de coeficientes positivos de uma função real.

H-conv está confirmada também pelo artigo original, §1.1, que fixa séries
de Fourier de período `4pi` e índice `m`. Assim a frequência é `m/2`;
o conjugador desloca o índice inteiro `m` por 1. Fonte primária:
[RT, §1.1](https://arxiv.org/pdf/1203.3766).

Esta inspeção resolve a suspeita de uma restrição oculta no seletor de DOFs,
mas não demonstra isomorfismo entre truncamentos finitos: um deslocamento
de frequência altera a borda da caixa de truncamento. A conjugação é uma
afirmação sobre os espaços completos; sua aplicação ao espectro infinito
ainda depende de H-iso/H-afim/H-constr e da transferência certificada.

`build_contour.py --imag-center 1/2 --height 1/4` implementa agora o
retângulo reduzido de B. Seu winding continua sendo o da matriz finita.

## 7. Realidade: onde este argumento costuma furar

O pencil do setor A é **real** na base de graus de liberdade de RT (a matriz
gravada por `SpectralMatrix` tem entradas diádicas reais; é o que
`analyze_spectrum.py` lê). Logo, com `κ` = conjugação complexa,

```
κ L_A(s) κ = L_A( conj s )      ⟹      Σ_A = conj(Σ_A) .      (7.1)
```

Combinando com (4.1) e (3.1): `conj(Σ_B) = conj(Σ_A) + i/2 = Σ_A + i/2 = Σ_B`.
Ou seja, `Σ_B` também é simétrico por conjugação — mas **por um motivo
transportado, e a estrutura real transportada não é a de `A`**. Explicitamente,
para `g ∈ A`,

```
M⁻¹ κ_B M g = M⁻¹ conj(e^{iτ/2} g) = e^{−iτ} conj(g) ,
```

isto é: lida em coordenadas de `A`, a realidade do setor `B` é
`κ` **composta com a translação de Floquet `−i`**. A simetria antilinear
correspondente no plano `s` é

```
s ↦ conj(s) − i ,      cujo conjunto fixo é  Im(s) = −1/2  (mod ℤ) .
```

**Consequência dura.** Enquanto os modos reais do setor `A` vivem na reta
`Im(s) = 0`, os modos reais do setor `B`, transportados para `A`, vivem na reta
`Im(s) = 1/2` (mod `ℤ`). Portanto:

* a janela "natural" e simétrica `Im(s) ∈ [−1/2, 1/2]`, de largura `1`, é a
  **pior escolha possível**: ela põe os modos reais de `B` exatamente sobre o
  bordo horizontal do contorno, que é justamente onde o certificado de winding
  exige ausência de espectro. Um modo instável real em `B` — o candidato mais
  provável a existir, se existir — cairia em cima da linha de integração;
* a janela deslocada `Q̃ = {−1/4 < Im s ≤ 3/4}` do Corolário 5.3 coloca as duas
  retas de realidade (`Im = 0` para `A`, `Im = 1/2` para a imagem de `B`) no
  **interior** das duas metades, e os bordos `Im = −1/4`, `Im = 3/4` em posição
  genérica. É essa a escolha correta, e a razão não é estética.

**Uso legítimo da realidade (economia opcional — e hoje desaconselhada).** De
(7.1) e (4.1), o divisor restrito a `Q̃` é simétrico pelas duas reflexões
`s ↦ conj(s)` (eixo `Im = 0`) e `s ↦ conj(s) + i` (eixo `Im = 1/2`). Logo
bastaria varrer `0 ≤ Im s ≤ 1/2` e duplicar, desde que se contem separadamente
as multiplicidades **sobre** as duas retas de simetria. Isso corta o contorno
pela metade e consome (7.1), que depende da matriz ser real na base usada.

**Correção posterior: essa economia colide com o verificador de winding.** O
contorno reduzido fecha *sobre* o eixo real, e por (7.1) `det` é real em todo o
segmento de fechamento. O classificador de
`formal/ChoptuikFormal/WindingInterval.lean` decide o incremento pelo cruzamento
do raio `ℝ_{<0}`: uma sequência inteira de arestas passaria a ter valores em
cima do raio, e todas devolveriam `none`. Ou se abandona a economia, ou se gira
o raio (multiplicar todas as caixas por `i`, exato em aritmética racional e sem
alterar o winding). A análise completa, com a tabela de tradução dos ramos, está
em [`WINDING_DATA.md`](WINDING_DATA.md) §2.7, que **recomenda abandonar a
economia**: no contorno completo da §8 a paridade de `Im det` ao longo dos lados
verticais é uma garantia de sinal que se perde ao girar o raio, e vale mais do
que a metade do trabalho economizada.

## 8. Consequência operacional exata

Fixado `R > 0` tal que não há espectro com `Re s ≥ R` (obrigação F3) e uma
condição de fase que trate o modo neutro em `s = 0` (§4 de
`SPECTRAL_PROBLEM.md`), o contorno a percorrer é **um único retângulo no plano
`s` para o operador `L_A`**:

```
Γ = ∂ { σ₀ ≤ Re s ≤ R ,  −1/4 ≤ Im s ≤ 3/4 } ,      σ₀ > 0 pequeno,
```

com a decomposição em dois sub-retângulos que separam C3 de C4:

```
Γ_A = ∂ { σ₀ ≤ Re s ≤ R ,  −1/4 ≤ Im s ≤ 1/4 }   →  winding alvo 1  (C3) ,
Γ_B = ∂ { σ₀ ≤ Re s ≤ R ,   1/4 ≤ Im s ≤ 3/4 }   →  winding alvo 0  (C4) .
```

`Γ_A` e `Γ_B` compartilham a aresta `Im s = 1/4`, percorrida em sentidos
opostos; a soma dos dois windings é o winding sobre `Γ` e é a contagem física
do Corolário 5.2. Na prática o produtor calcula as caixas dessa aresta uma só
vez e as reutiliza com o incremento trocado de sinal — o que é, de quebra, um
teste de consistência barato. Exigências adicionais que isto impõe, e que têm
de constar do certificado: ausência de espectro sobre as três retas
`Im s ∈ {−1/4, 1/4, 3/4}` no intervalo `σ₀ ≤ Re s ≤ R`, e sobre `Re s = σ₀` e
`Re s = R`.

**Duas regras de colocação de nós, que não são opcionais.** A reta `Im s = 0`
(realidade exata de `Γ_A`, por (7.1)) e a reta `Im s = 1/2` (realidade do setor
B transportada, aproximada no truncamento) cruzam os lados **verticais** dos
dois retângulos. Nesses pontos `det` é real e, quando é negativo, cai sobre o
raio `ℝ_{<0}` usado pelo classificador. Portanto:

1. nenhum nó do contorno sobre `Im s = 0` (em `Γ_A`) ou `Im s = 1/2` (em
   `Γ_B`) — use um número **ímpar** de subdivisões nos lados verticais, para que
   o cruzamento caia no interior de uma aresta, onde é classificado
   corretamente como `±1`;
2. nós verticais em posições simétricas `Im s = ±h` em `Γ_A`: por (7.1),
   `Im det(σ₀ + iy)` é **ímpar** em `y`, de modo que os sinais nos dois nós são
   automaticamente opostos e só um deles precisa ser resolvido por precisão.

Uma subdivisão uniforme com número par de segmentos viola a regra 1 e torna o
certificado irrecuperável por refinamento. Detalhes e justificativa em
[`WINDING_DATA.md`](WINDING_DATA.md) §2.

Em termos de implementação: **nada precisa ser complexificado**. Os harmônicos
semi-inteiros, que exigiriam reescrever a álgebra `MultiField`, nunca aparecem;
eles são absorvidos pela translação do argumento `s`, que o localizador já
trata (os autovalores complexos da matriz real de `L_A(0)` já são calculados de
uma vez). O custo de C4 passa a ser o de percorrer um segundo retângulo com a
mesma maquinaria de C3.

## 9. Diagnóstico numérico (ponto flutuante — não é prova)

Na matriz 4×12 já gravada em `build/spectrum/rt-A-4x12.dat` (dimensão real 138):

| região (`Re s > 0`) | autovalores |
|---|---:|
| todo o semiplano `Re s > 0`, sem faixa | 15 |
| faixa de A, `−1/4 < Im s ≤ 1/4` | 3 (todos reais: `0,7096`, `0,4387`, `0,2318`) |
| faixa da imagem de B, `1/4 < Im s ≤ 3/4` | **0** |

Três leituras, todas honestas:

1. a restrição à faixa fundamental **não é cosmética**: ela leva a contagem
   bruta de 15 para 3. Sem ela, os translados de Floquet e os artefatos de
   truncamento entram na conta;
2. o zero na faixa da imagem de `B` é o primeiro sinal numérico a favor de C4 —
   e é só isso: um sinal, em ponto flutuante, na malha mais grosseira, cujo
   resíduo de constraints é `1,1·10⁻¹`. Os dois candidatos reais extras
   (`0,4387`, `0,2318`) na faixa de `A` são exatamente o tipo de raiz que a
   obrigação S3 (constraints) e S4 (gauge) existem para eliminar ou explicar;
3. o truncamento **quebra** a invariância exata (4.1): o translado previsto de
   `0,7096` por `+i` aparece em `0,7143 + 1,0562 i`, com erro `≈ 0,06`. Isso é
   esperado — deslocar o índice de Fourier em `1` empurra modos para fora da
   caixa retida — e é mais um motivo para que a contagem em `B` seja feita por
   winding sobre `Γ_B`, e não por inspeção de autovalores deslocados.

## 10. Veredito

* A conjugação `M⁻¹L_B(s)M = L_A(s+i/2)` **está correta**, com a conta explícita
  nos Lemas 2.1–2.3, e o sinal impresso por `analyze_spectrum.py`
  (`s_B = s_A − i/2`) é o coerente com ela. O sinal, contudo, só está definido
  módulo `i`.
* A redução (Lema 5.1 e Corolários 5.2, 5.3) **fecha no nível espectral**: uma
  única computação em `A` sobre uma faixa de largura `1` em `Im(s)` cobre os dois
  setores, e C4 vira o winding de `L_A` sobre o retângulo transladado `Γ_B`.
  Isso muda a economia da obrigação S2: o item "implementação nativa de `B`"
  pode ser **retirado** do caminho crítico; o que fica é conferir H-dof.
* A interação com a realidade **não quebra a redução, mas dita a janela**: os
  modos reais de `B` vivem sobre `Im(s_A) = 1/2`, de modo que qualquer janela
  simétrica de largura `1` os coloca sobre o contorno. A janela `(−1/4, 3/4]` é
  a correção.
* Consequência operacional descoberta depois, ao revisar o verificador de
  winding: a simetria de realidade **força** `det` a ser real em quatro pontos
  do bordo dos dois retângulos, e é de um desses pontos que sai o `+1` de C3.
  Isso não invalida nada, mas impõe as duas regras de colocação de nós da §8 e
  desaconselha a economia da §7. Ver [`WINDING_DATA.md`](WINDING_DATA.md).
* O que a redução **não** faz: não demonstra S3 (propagação das constraints para
  `s` complexo), não demonstra F1–F3, não fornece o `ε_N` de C1, e não dispensa
  a hipótese H-dof sobre a seleção de graus de liberdade do código de RT.
  Continua valendo a regra do projeto: enquanto essas obrigações estiverem
  abertas, as raízes são candidatas, não modos.

## 11. Correção de 16/09/2026 (auditoria física)

A leitura 2 da §9 está **superada**. O zero na imagem da faixa de B em 4×12
é artefato da malha grossa. De 6×18 em diante, `L_A` tem exatamente uma raiz
na faixa, e ela converge: `0,063467+0,479488i` (6×18), `0,041834+0,500373i`
(8×24), `0,041373+0,500111i` (10×30), `0,041421+0,500008i` (12×36). Isso
corresponde a um modo real do setor B com `s_B≈0,0414`, fora do retângulo
`1/8<=Re s<=1`. A raiz tem a assinatura de violação de constraints:
‖Ch‖/‖h‖₂≈0,50 estável, colinear à raiz sharp `0,041419+0,5i`.

Consequência: o alvo C4 deixa de ser "winding 0" e passa a ser
multiplicidade física 0, isto é, winding de `L_A` em `Γ_B(σ0)` menos as
raízes com testemunho S3. A redução `Sigma_B = Sigma_A - i/2` e a reta de
realidade em `Im s = 1/2` foram confirmadas. Regiões descobertas e plano:
[`AUDITORIA_FISICA_16SET.md`](AUDITORIA_FISICA_16SET.md) §7 e §9.
