# O que um produtor de certificados de winding precisa entregar

`formal/ChoptuikFormal/WindingInterval.lean` recomputa o winding a partir de
três caixas racionais por aresta. Este documento responde às duas perguntas que
ficam abertas do lado de fora do kernel:

1. o raio `ℝ_{<0}` usado na classificação sobrevive à simetria de realidade dos
   contornos fixados em [`SECTOR_B.md`](SECTOR_B.md) §8?
2. quem produz o `edgeBox`, isto é, o enclosure de `f` sobre o **segmento
   inteiro**?

**Resumo das conclusões.** (1) O módulo está correto e o raio `ℝ_{<0}` é
utilizável nos dois retângulos de §8, mas **a simetria de realidade força `f` a
ser real em quatro pontos do bordo**, e daí saem duas regras de construção que
não são opcionais (§2.6). Existe um caso em que o raio `ℝ_{<0}` é
**inutilizável**: o contorno reduzido pela metade que a §7 de `SECTOR_B.md`
sugeria como economia. Lá a correção é girar o raio, o que se faz **sem tocar
no Lean** (§2.7). (2) Jacobi **basta** para o `edgeBox`, mas na forma
multiplicativa e com o traço *calculado*, não majorado por `n·‖·‖`; a estimativa
completa está na §3 e sua composição aritmética está verificada pelo kernel em
`formal/ChoptuikFormal/DeterminantBound.lean`. O que Jacobi **não** dá são os
valores nos nós; a §3.4 mostra como evitá-los por completo.

## 1. A função certificada

```
f(z) = det[I + (z - z₀) C_N] ,     C_N = truncamento de Galerkin de C = J R₀ .
```

Como `I + (z-z₀)R₀ = R₀ P(z)`, vale `f(z) = det(R₀)·det P(z)`: os zeros de `f`
são o espectro e `winding(f) = winding(det P)`, porque `det R₀ ≠ 0` é constante.

Disso decorre uma liberdade que a §3.4 usa: **o winding não muda se `f` for
multiplicado por qualquer constante complexa não nula**. Em particular pode-se
certificar `h(z) = f(z)/f(z₀-nó)` em vez de `f`, e `h(z₀-nó) = 1` exatamente.

## 2. O raio `ℝ_{<0}` contra a simetria de realidade

### 2.1 O que é forçado a ser real, e com que exatidão

`scripts/analyze_spectrum.py` lê a matriz de `SpectralMatrix` num array **real**
(`np.zeros(..., dtype=float)`, entradas diádicas exatas): `C_N` é uma matriz
real na base de graus de liberdade de RT, e a caixa de truncamento é estável por
conjugação (`n ↦ -n`). Logo, para `z₀` real,

```
f(conj z) = conj f(z)   exatamente, não aproximadamente.          (2.1)
```

Há uma segunda simetria, a do setor B transportada (`SECTOR_B.md` §7): a
involução antilinear `Λ' = (multiplicação por e^{-iτ_R}) ∘ conjugação` satisfaz
`Λ'⁻¹ L_A(s) Λ' = L_A(conj(s) + i)`, cuja reta fixa é `Im s = 1/2`. **Mas `Λ'`
envolve a translação inteira de Floquet, que desloca o índice de Fourier em `1`
e portanto não preserva a caixa `|n| ≤ N_F`.** Consequência:

```
f_N é exatamente real sobre Im s = 0 ;
f_N é apenas aproximadamente real sobre Im s = 1/2 .              (2.2)
```

O defeito é mensurável: em 4×12 o translado previsto de `s ≈ 0,7096` por `+i`
aparece em `0,7143 + 1,0562 i`, erro `≈ 0,06` (diagnóstico em ponto flutuante).

### 2.2 `Γ_A`: o cruzamento do raio é forçado, e é de onde vem o `+1`

`Γ_A = ∂{σ₀ ≤ Re s ≤ R, −1/4 ≤ Im s ≤ 1/4}` encontra o eixo real em exatamente
dois pontos, `(σ₀, 0)` e `(R, 0)`, ambos no **interior dos lados verticais**. Por
(2.1), `f_N` é real nesses dois pontos. Mais: sobre `[σ₀, R]` (que está no
*interior* de `Γ_A`, não no bordo) `f_N` é uma função real cujos zeros são os
autovalores reais; se o número deles, com multiplicidade, for **ímpar** — e a
contagem alvo é `1` —, então

```
sinal f_N(σ₀) = − sinal f_N(R) ,
```

ou seja **exatamente um dos dois lados verticais tem `f_N` real e negativo no
ponto em que cruza o eixo real**. Isso não é um caso patológico a evitar: é o
cruzamento do raio que produz o `+1` do winding. O ramo 3 de `edgeCrossing`
existe justamente para ele.

O que a simetria dá de brinde: por (2.1), `Im f_N(σ₀ + iy)` é **ímpar** em `y`.
Portanto, escolhidos nós em `y = ±h`, os sinais de `Im f_N` neles são
automaticamente opostos: basta resolver **um** sinal, e o outro vem de graça.
Esse é o ponto mais delicado de toda a computação e é o único em que a simetria
ajuda em vez de atrapalhar.

### 2.3 `Γ_B`: a faixa contém `Im s = 1/2`, mas o bordo só a toca em dois pontos

`Γ_B = ∂{σ₀ ≤ Re s ≤ R, 1/4 ≤ Im s ≤ 3/4}`. A reta `Im s = 1/2` está no
**interior** de `Γ_B`; o bordo a encontra apenas nos dois pontos `(σ₀, 1/2)` e
`(R, 1/2)`, de novo no interior dos lados verticais. Por (2.2), ali `f_N` é
*quase* real: `Im f_N` é pequeno, da ordem do defeito de truncamento, mas não é
zero e não é ímpar. Consequências:

* se `f_N` for quase-real **positivo** nesses pontos, não há problema algum: o
  ramo 1 (`edgeBox ⊆ {Re > 0}`) classifica com folga;
* se for quase-real **negativo**, há um cruzamento do raio ali, exatamente como
  em `Γ_A` — legítimo e classificável —, mas **sem** a garantia de paridade da
  §2.2: o produtor tem de resolver o sinal de um `Im f_N` pequeno por precisão
  bruta, nos dois nós vizinhos.

Note-se que `Γ_A` e `Γ_B` **compartilham a aresta `Im s = 1/4`**, percorrida em
sentidos opostos. O produtor deve calcular as caixas dessa aresta uma única vez
e reutilizá-las com o incremento trocado de sinal; a soma dos dois windings é o
winding sobre `Γ = ∂{σ₀ ≤ Re s ≤ R, −1/4 ≤ Im s ≤ 3/4}`, e a aresta partilhada
se cancela. É também um teste de consistência barato.

### 2.4 Lados horizontais e cantos: nada é forçado

Sobre `Im s = ±1/4` e `Im s = 3/4` nenhuma simetria força `f_N` a ser real. A
relação (2.1) liga o lado `Im s = -1/4` ao lado `Im s = +1/4` (valores
conjugados, percurso invertido) e os cantos `(σ₀, ±1/4)` entre si; isso torna
metade dos dados de `Γ_A` redundante, mas não coloca nenhum valor sobre o raio.

### 2.5 O que de fato pode quebrar

Só há um modo de falha, e ele é de **colocação de nós**:

> se um nó `z_k` cair sobre uma reta de realidade (`Im s = 0` em `Γ_A`,
> `Im s = 1/2` em `Γ_B`) num lado vertical onde `f_N` é negativa, então
> `f(z_k) ∈ ℝ_{<0}` está **sobre** o raio; nenhuma caixa que o contenha é
> estritamente `Im > 0` nem `Im < 0`, as duas arestas adjacentes devolvem
> `none`, e **refinar não adianta**, porque o nó continua onde estava.

Uma subdivisão uniforme do lado vertical `[−1/4, 1/4]` em um número **par** de
segmentos põe um nó exatamente em `Im s = 0`. É o erro natural a cometer.

### 2.6 Duas regras de construção, obrigatórias

1. **Nenhum nó sobre uma reta de realidade.** Nos lados verticais de `Γ_A`, use
   um número ímpar de subdivisões (ou um deslocamento explícito) de modo que
   `Im s = 0` caia no **interior de uma aresta**, não em um nó. O mesmo para
   `Im s = 1/2` nos lados verticais de `Γ_B`. Assim o cruzamento do raio é
   classificado pelo ramo 3, que é o comportamento correto.
2. **Nós simétricos em `Γ_A`.** Coloque os nós verticais em `y = ±h` de modo
   que a paridade de `Im f_N` (§2.2) force os sinais opostos; resolva por
   precisão apenas o sinal em `y = +h`.

Com essas duas regras, o raio `ℝ_{<0}` é seguro nos contornos de `SECTOR_B.md`
§8, e **não há motivo para mudar o módulo**.

### 2.7 Onde `ℝ_{<0}` é inutilizável, e o raio alternativo

Há um caso em que nenhuma colocação de nós salva: o **contorno reduzido pela
metade** que `SECTOR_B.md` §7 sugeria como economia (varrer só `0 ≤ Im s ≤ 1/2`
e usar a conjugação para dobrar). Esse contorno **fecha sobre o eixo real**, e
por (2.1) `f_N` é real em *todo* o segmento de fechamento. Onde esse segmento
passa por valores negativos, uma sequência inteira de arestas fica com
`edgeBox` cortando o raio e com nós exatamente sobre ele: todas devolvem `none`.
A economia da §7 e o raio `ℝ_{<0}` são **incompatíveis**.

O raio alternativo seguro é qualquer um transverso ao eixo real; o canônico é o
semieixo imaginário positivo `i·ℝ_{>0}`. A razão estrutural: os valores que a
simetria força a existir são **reais**, e um raio transverso ao eixo real nunca
os contém (a origem já está excluída).

**A troca não exige mudar `WindingInterval.lean`.** Multiplicar por `i` leva
retângulos racionais em retângulos racionais, exatamente:

```
i · ([a,b] × [c,d])  =  [−d,−c] × [a,b] ,
```

e `winding(i·f) = winding(f)`, porque `i` é constante não nula. Basta o produtor
emitir as caixas de `g = i·f`; o classificador de `ℝ_{<0}` aplicado a `g` é
exatamente o classificador de `i·ℝ_{>0}` aplicado a `f`. A tradução dos ramos,
para conferência:

| ramo em `g = i·f` | condição sobre `f` |
|---|---|
| `Re g > 0` | `Im f < 0` |
| `Im g > 0` / `Im g < 0` | `Re f > 0` / `Re f < 0` |
| `Re g < 0`, `Im g : + → −` | `Im f > 0`, `Re f : + → −` (incremento `+1`) |

**Recomendação.** Manter `ℝ_{<0}` e as duas regras da §2.6 no contorno completo
de §8 — porque lá a paridade de `Im f_N` é uma garantia de sinal que se perde ao
girar. Usar a rotação por `i` se e somente se o contorno reduzido pela metade
for adotado, e nesse caso perde-se a garantia de paridade em troca de poder
fechar sobre o eixo real. Como a §8 não depende da redução pela metade, o
caminho de menor risco é: **contorno completo, raio `ℝ_{<0}`, regras da §2.6**.

## 3. Quem produz o `edgeBox`

### 3.1 Jacobi dá um bound multiplicativo

`F(z) = I + (z-z₀)C_N` tem `F'(z) = C_N`, e a fórmula de Jacobi dá

```
f'(z) = f(z) · tr(F(z)⁻¹ C_N) .                                    (3.1)
```

O bound que sai daí é sobre `f'/f`, não sobre `f'`: é **multiplicativo**. Isso é
uma vantagem (a estimativa não depende da escala de `f`, que é astronômica em
dimensão 612) e obriga a converter para uma caixa aditiva na §3.3.

(3.1) só vale onde `F(z)` é invertível, e precisa de `‖F(z)⁻¹‖` **em toda a
aresta** — que é exatamente o que o tile de Neumann já entrega:
`β = ‖V‖/(1-η)` com `η = ‖I - V F(z_j)‖ + ‖V‖·δ_F`.

### 3.2 O traço tem de ser calculado, não majorado

A rota ingênua `|tr A| ≤ n·‖A‖` custa o fator `n = 612` e torna a malha
inviável (`K ≈ 6·10⁵` obrigaria arestas de comprimento `< 10⁻⁶`, isto é `~10⁷`
arestas). A rota correta separa o termo dominante:

```
tr(F⁻¹C_N) = tr(V C_N) + tr((F⁻¹ - V) C_N) ,
|tr((F⁻¹-V)C_N)| ≤ ‖F⁻¹-V‖_F ‖C_N‖_F ≤ ‖I - V F‖_F · β · ‖C_N‖_F ,
```

usando Cauchy--Schwarz de Frobenius (`|tr XY| ≤ ‖X‖_F‖Y‖_F`),
`F⁻¹ - V = (I - VF)F⁻¹` e `‖AB‖_F ≤ ‖A‖_F‖B‖`. O termo `tr(V C_N) =
Σ_{ij} V_{ij}(C_N)_{ji}` é uma soma barata, calculada exatamente em racionais.
Resultado:

```
K := |tr(V C_N)| + ‖I - V F‖_F · β · ‖C_N‖_F   ≥   sup_aresta |f'/f| .   (3.2)
```

`K` é agora da ordem da própria derivada logarítmica (`Σ_j 1/(z - s_j)` sobre os
autovalores), tipicamente dezenas, não `10⁵`. Como as normas de Frobenius
envolvem raiz quadrada, o produtor entrega **majorantes racionais** delas; é
esse majorante que entra no certificado.

### 3.3 Do bound multiplicativo para a caixa

Seja `L` o comprimento da aresta (racional: o contorno de §8 tem lados paralelos
aos eixos) e `m ≥ |f(z_k)|` um majorante racional do módulo no nó inicial — de
uma caixa racional basta `m = max|Re| + max|Im|`, que majora o módulo sem raiz.
Com `S = sup_aresta |f|`, de (3.2) e do teorema fundamental do cálculo,

```
|f(z)| ≤ m + K·L·S  para todo z na aresta  ⟹  S ≤ m/(1 − K·L)   (se K·L < 1),
|f(z) − f(z_k)| ≤ K·L·S ≤ K·L·m/(1 − K·L)  =:  δ .                (3.3)
```

E então

```
edgeBox := startBox inflado por δ em Re e em Im.                  (3.4)
```

Tudo racional e exato. A condição operacional é `K·L < 1/2`, que garante
`δ < m`; acima disso a caixa engole a origem e nenhuma aresta classifica.

**Verificação pelo kernel.** A composição (3.2)--(3.4) está implementada em
`formal/ChoptuikFormal/DeterminantBound.lean` (`LogDerivativeData.bound`,
`EdgeEnclosureData.radius`, `EdgeEnclosureData.valid`), com regressões `decide`.
Exemplo típico, com números plausíveis e não certificados: `|tr(VC_N)| = 12`,
`‖I-VF‖_F = 10⁻⁵`, `β = 50`, `‖C_N‖_F = 200` dão `K = 121/10`; com `L = 1/500` e
`m = 1`, a inflação é `δ = 121/4879 ≈ 2,5 %`. Para um contorno de perímetro
`≈ 12` isso dá `≈ 6000` arestas, e a classificação pede `|Im f| > 0,025` nos nós
— viável fora da vizinhança imediata do cruzamento forçado da §2.2.

### 3.4 O que Jacobi **não** dá: os valores nos nós

(3.3) enclausura `f` na aresta *em termos de* `m ≥ |f(z_k)|`. O valor no nó não
vem de Jacobi. Duas rotas:

* **Rota absoluta.** Calcular um enclosure verificado de `det F(z_k)` para cada
  nó. Para `n = 612` isto é caro e mal condicionado; a via usual
  (`det F = det(V)⁻¹ det(I - G)` com `‖G‖ = η`) exige `n·η ≲ 1`, isto é
  `η ≲ 1/612`, e ainda deixa `det V` para calcular exatamente. Não recomendada.
* **Rota relativa normalizada (recomendada).** Como o winding é invariante por
  divisão por constante não nula (§1), certifique
  `h(z) = f(z)/f(z_{nó 0})`. Então `h(nó 0) = 1` **exatamente**, e os nós
  seguintes saem por propagação:

  ```
  h(z_{k+1}) = h(z_k) · r_k ,   r_k = det[I + (z_{k+1}-z_k) C_N F(z_k)⁻¹] .
  ```

  Nenhum determinante absoluto é jamais calculado. O `m` da §3.3 sai da caixa
  propagada, e o valor `1` do primeiro nó é o que justifica `m = 1` no exemplo
  acima.

  **Atenção ao acúmulo.** Se `r_k` for enclausurado apenas por
  `|det(I+Δ)-1| ≤ (1+‖Δ‖)^n - 1 ≈ n‖Δ‖`, o erro relativo acumulado ao longo do
  contorno é `Σ_k n‖Δ_k‖ ≈ n‖C_N‖β·|Γ|`, que **não diminui ao refinar**: o termo
  de primeira ordem é informação, não erro. É preciso capturar essa primeira
  ordem exatamente — de novo o traço, `log r_k ≈ tr((z_{k+1}-z_k)C_N F(z_k)⁻¹)`
  — deixando um resto de segunda ordem cuja soma é `O(L)` e portanto **vai a
  zero linearmente sob refinamento**. É a mesma decisão da §3.2, e é o ponto em
  que um produtor ingênuo se quebra sem perceber.

### 3.5 Checklist do produtor

Por tile do contorno:

| dado | origem | já existe? |
|---|---|---|
| `V` (inverso aproximado) e `‖I - V F(z_j)‖` | tile de Neumann | sim, `RationalCertificate.Tile` |
| `β = ‖V‖/(1-η)` | tile de Neumann | sim |
| `matrixVariation = raio_do_tile · ‖C_N‖` | **aritmética**, não EDP (ver §3.6) | sim, calculável |
| `‖C_N‖_F`, `‖I - V F‖_F` (majorantes racionais) | soma de quadrados sobre a matriz | não |
| `tr(V C_N)` exato | soma `Σ V_{ij}(C_N)_{ji}` | não |

Por aresta: `L`, `m`, e daí `K`, `δ`, `edgeBox` por (3.2)--(3.4).

Fora disso, e continuando como hipótese analítica: `ε_N` (obrigação C1, ver
[`TAIL_BOUND.md`](TAIL_BOUND.md)) e a própria validade da transferência de
Neumann do truncado para o operador.

### 3.6 Um bound de EDP a menos

Para o pencil afim, `F_N(z) - F_N(z_j) = (z - z_j)C_N` **exatamente**, logo

```
matrixVariation  =  raio_do_tile · ‖C_N‖ ,
```

que é aritmética sobre a matriz racional. O ledger listava `matrixVariation` ao
lado de `tailError` como se ambos fossem bounds de EDP a demonstrar; só
`tailError` é. A composição está em `DeterminantBound.affineMatrixVariation`.

## 4. Revisão de `WindingInterval.lean`

**O módulo está correto.** Reverifiquei de forma independente, e sem encontrar
defeito, os quatro pontos em que esse tipo de argumento costuma falhar:

* **o argumento do piso.** Com `m_k = ⌊(θ_k - π)/2π⌋`, cada ramo mantém
  `⌊(θ-π)/2π⌋` constante ao longo da aresta (ramos 1 e 2) ou o faz saltar de
  exatamente `±1` (ramo 3). Recomputei os três: em `{Re>0}`,
  `(θ-π)/2π ∈ (m-3/4, m-1/4)`; em `{Im>0}`, `(m-1/2, m)`; em `{Im<0}`,
  `(m-1, m-1/2)`; em `{Re<0}`, `(m-1/4, m+1/4)`, com piso `m-1` sse `Im f > 0`.
  Confere com o cabeçalho do módulo;
* **robustez a cruzamentos múltiplos.** No ramo 3, `θ` permanece num único
  intervalo `(π/2+2πm, 3π/2+2πm)` por conexidade, de modo que `m` é constante na
  aresta mesmo que `Im f` oscile várias vezes. O incremento depende só dos dois
  nós. O módulo não anuncia isso, mas é verdade e é o que torna o critério
  robusto;
* **orientação.** `Im : + → −` no semiplano esquerdo vale `+1`; conferido em
  `f(z)=z` sobre o círculo anti-horário, onde o único cruzamento é em `θ = π`;
* **consistência nos nós.** A comparação de **caixas** em `chained` faz mais do
  que o módulo afirma: ela impede leituras de sinal contraditórias num nó
  partilhado. Um ciclo de duas arestas com incrementos `+1, +1`, por exemplo, é
  impossível, porque exigiria que a mesma caixa tivesse `imLo > 0` e
  `imHi < 0`.

Duas observações, nenhuma delas um defeito de correção:

1. **Incompletude conservadora do ramo 3.** Quando `edgeBox ⊆ {Re < 0}` e os
   **dois** nós têm `Im > 0` (ou os dois `Im < 0`), o incremento é
   demonstravelmente `0` — pelo mesmo argumento do módulo: `m` é constante e os
   dois pisos valem `m-1` (resp. `m`). Hoje `edgeCrossing` devolve `none` nesse
   caso e força refinamento. Acrescentar os dois casos é sólido e barato, e o
   ganho não é cosmético: é precisamente na vizinhança do cruzamento forçado da
   §2.2 que as arestas ficam inteiras no semiplano esquerdo sem cruzar o raio.
2. **A nota sobre `mkRat` está obsoleta.** O cabeçalho afirma que `n / d` trava
   a redução no kernel. Isso é verdade por omissão de `unseal`: com
   `unseal Rat.mul Rat.add Rat.sub Rat.inv` (como em `TailBound.lean` e
   `DeterminantBound.lean`), `n / d`, `+`, `-`, `*` e `^` reduzem normalmente no
   kernel. O produtor pode continuar emitindo `mkRat` — é inofensivo —, mas a
   restrição não é do Lean e não deve ser propagada como se fosse.

## 5. O que continua sendo hipótese

Nada neste documento demonstra que as caixas contêm `f`. A cadeia de
responsabilidade fica assim:

| afirmação | quem garante |
|---|---|
| `startBox ∋ f(z_k)`, `endBox ∋ f(z_{k+1})` | `scripts/build_contour.py`, exato (§6) |
| `edgeBox ∋ f(z)` em todo o segmento | `scripts/build_contour.py`, por Taylor exato (§6) |
| `K` majora `sup |f'/f|` | `scripts/build_contour.py`, exato via `p'/p` (§6.2) |
| a composição `K, L, m ↦ δ` | kernel (`DeterminantBound`) |
| forma, geometria, encadeamento, fechamento, `0 ∉ edgeBox`, soma | kernel (`WindingInterval`) |
| que `p` é mesmo `det(A + zI)` | verificado contra Bareiss independente (§6.1) |
| que `A` (truncamento 4×12) representa o operador infinito | **ninguém**: é C1, S3, S4, tudo aberto |

As três primeiras linhas eram, na versão anterior deste documento, hipóteses do
produtor. Desde a §6 elas são **calculadas**, para esta matriz, em aritmética
exata. A última linha é a que continua sendo o abismo: o certificado é sobre o
determinante do truncamento, e nada mais.


## 6. Primeira execução sobre dados reais (matriz 4×12)

`scripts/build_contour.py` liga a matriz de `build/spectrum/rt-A-4x12.dat` ao
verificador. Esta seção relata a primeira execução, com os números medidos.
**Nada aqui é prova de coisa alguma sobre o operador infinito**: o objeto
certificado é o determinante do truncamento de Galerkin, cuja origem é
numérica. O aviso está no cabeçalho dos arquivos gerados.

A matriz: `n = 138`, entradas diádicas exatas `A = M · 2^{-72}` com `M` inteira,
7528 entradas não nulas (40 % de densidade), maior coeficiente com 75 bits.

### 6.1 Arquitetura: a cadeia relativa, calculada pelo caminho barato

A prescrição da §3.4 é a cadeia `h(z_{k+1}) = h(z_k)·det[I + (z_{k+1}-z_k)C_N F(z_k)^{-1}]`.
No pencil afim, `C_N = (A + z₀I)^{-1}` e `F(z) = (A+zI)(A+z₀I)^{-1}`, logo o
fator relativo é `det(A+z_{k+1}I)/det(A+z_kI)` e a cadeia **telescopa
exatamente** para `h(z_k) = p(z_k)/p(z₀)` com `p(z) = det(A + zI)`. Os dois
caminhos produzem os **mesmos números**; o que muda é o custo. Medido:

| caminho | custo medido | por nó |
|---|---|---|
| cadeia nó a nó, determinante exato (Bareiss fraction-free) | **42 s por determinante** | 42 s |
| `p` inteiro por Hessenberg modular + CRT (358 primos de 31 bits) | **166 s uma vez** | ~2 ms (Horner) |

Para os 1946 nós do contorno completo isso é a diferença entre **23 h** e
**3 min**. A cota de Hadamard para os coeficientes é de 10734 bits, e o maior
coeficiente obtido tem 9994 bits — grandes, mas irrelevantes para inteiros de
Python.

**Verificação independente.** O polinômio foi conferido contra determinantes de
Bareiss exatos (algoritmo diferente, aritmética diferente) em `y = 1` e
`y = 2^71`: igualdade exata nos dois pontos, 44 s cada (`--verify`).

### 6.2 O traço, exato

Com `p` em mãos, `tr(F(z)^{-1}) = p'(z)/p(z)` **exatamente**. A decomposição da
§3.2 (`tr(VC_N)` calculado, resto por Frobenius) deixa de ser necessária: não há
inverso aproximado `V`, nem resíduo, nem o fator `n` que a §3.2 mostrou ser
proibitivo. O bound de Jacobi é montado a partir dos coeficientes de Taylor
exatos em cada nó, junto com o enclosure direto:

```
delta_taylor = Σ_{m≥1} |c_m| L^m ,
K = [Σ_{m≥1} m|c_m|L^{m-1}] / [|c_0| - Σ_{m≥1}|c_m|L^m] ,  K·L < 1/2 conferido.
```

### 6.3 O que não fechou na primeira tentativa

Com 16 subdivisões por lado horizontal e 9 por vertical (50 nós), **0 de 50
arestas classificaram**: `δ/|p|` da ordem de `7`, isto é `K·L ≈ 7,3`, contra o
`< 1/2` exigido. A causa é aritmética, não estrutural: `K = Σ_j 1/|z - s_j|`
sobre as 138 raízes vale até `134` sobre `Γ_A` (diagnóstico em ponto flutuante),
de modo que `L` precisa ficar abaixo de `~1,5·10^{-3}`. A correção foi refinar o
contorno; **nenhum parâmetro de tolerância foi afrouxado**.

### 6.4 Resultados

| contorno | retângulo | nós | `δ/|p|` máx | `K·L` máx | arestas que classificam | winding |
|---|---|---:|---:|---:|---:|---:|
| `Γ_A` completo (`SECTOR_B.md` §8) | `1/8 ≤ Re s ≤ 1`, `|Im s| ≤ 1/4` | 1946 (640+333 por lado) | 0,200 | 0,263 | **1946/1946** | **3** |
| caixa do candidato físico | `5/8 ≤ Re s ≤ 7/8`, `|Im s| ≤ 1/8` | 514 (128+129 por lado) | 0,274 | 0,401 | **514/514** | **1** |

As subdivisões verticais são ímpares (333 e 129), como manda a regra 1 da §2.6,
e os nós verticais são simétricos em `±h` (regra 2).

Os dois arquivos estão em `examples/contour-A-4x12-gamma-a-nodes.txt` e
`examples/contour-A-4x12-physical-nodes.txt`. Ambos passam por
`scripts/make_winding_certificate.py` sem erro, e o Lean emitido **fecha no
kernel** (`decide`, sem `native_decide`):

| certificado | arestas | tempo de `lean` | ajustes |
|---|---:|---:|---|
| caixa física | 514 | 144 s | `maxRecDepth 1000000` |
| `Γ_A` completo | 1946 | 673 s | `maxRecDepth 1000000`, `maxHeartbeats 0` |

O custo cresce essencialmente linear no número de arestas (3,8× arestas →
4,7× tempo).

### 6.5 O resultado que importa: `Γ_A` dá 3, não 1

O winding do truncamento sobre `Γ_A` é **3**, e as três raízes são reais:
`s ≈ 0,7096`, `0,4387`, `0,2318` (diagnóstico). A obrigação C3 pede **1**. A
diferença não é erro do verificador nem do contorno: são exatamente as raízes
que as obrigações **S3** (propagação das constraints) e **S4** (classificação de
gauge) existem para eliminar ou explicar, e que o resíduo de constraints
`1,14·10^{-1}` do 4×12 já sinalizava. Restringindo o retângulo ao entorno do
candidato físico, o winding é **1**, como C3 espera.

Dito com precisão: **a maquinaria de winding não é mais o gargalo**. Ela roda
ponta a ponta sobre a matriz real e fecha no kernel. O que falta para C3 é
analítico — C1, S3, S4 —, não computacional.

### 6.6 Correção à análise da §2: a normalização gira a reta de realidade

A §2.2 previa que o cruzamento do raio `ℝ_{<0}` seria **forçado** sobre um dos
lados verticais, no ponto em que o contorno cruza `Im s = 0`. A execução mostra
que isso **não aconteceu**: as arestas que contêm `Im s = 0` deram incremento
`0` nos dois contornos, e os cruzamentos apareceram em outros lugares (cinco no
caso físico, treze em `Γ_A`).

A previsão não estava errada sobre `p`; estava incompleta sobre **o que é
certificado**. Verificando em aritmética exata:

* `Im p(x) = 0` exatamente para `x` real, confirmando a simetria (2.1) no nível
  do truncamento;
* `p(1/8) < 0` e `p(1) > 0` — a troca de sinal existe, como manda o número
  ímpar (3) de raízes reais entre eles.

Só que o certificado não é sobre `p`, e sim sobre `h = p/p(z₀)`, e `z₀` é o
canto inferior esquerdo, um ponto **complexo**: `arg p(z₀) = −24,7°` em `Γ_A` e
`+102,4°` na caixa física. A normalização é uma rotação por constante não nula
— não muda o winding — mas **gira a reta de realidade para fora do raio**. Daí
os pontos forçados a serem reais em `p` caírem, em `h`, sobre uma reta oblíqua,
longe de `ℝ_{<0}`.

Consequências, ambas úteis:

1. o risco estrutural da §2.5 **é controlável pela escolha do nó base**:
   normalizar num nó complexo o dissolve; normalizar num nó real o reintroduz;
2. em compensação, a garantia de paridade da §2.2 (`Im p` ímpar em `y` dando os
   sinais opostos de graça) **também se perde** ao girar, e passa a ser
   acidental. Quem quiser a garantia tem de normalizar num ponto real e aceitar
   o cruzamento forçado.

A regra 1 da §2.6 (nenhum nó sobre a reta de realidade) continua valendo e
continua barata; ela foi respeitada nas duas execuções. O que muda é o enunciado:
ela é uma regra sobre **a função efetivamente certificada**, não sobre `det`.
