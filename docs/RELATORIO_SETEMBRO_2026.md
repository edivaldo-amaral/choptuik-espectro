# O problema espectral de Choptuik: viabilidade, mapa e um resultado

> **Atualização 03/10/2026:** relatório de 17–22/09, anterior à rota do Rouché e a T2. Ver `PANORAMA.md`.

Relatório consolidado, 17–22 de setembro de 2026. Reúne numa narrativa só o que
está disperso em `ALTERNATIVE_ROUTES.md`, `MAJORANT_ROUTE.md`, `TT_CEILING.md`,
`NAO_FAZER.md`, `METODO.md` e `ESTADO_DA_ARTE.md`. Cada afirmação aqui aponta
para o documento e o script que a sustentam.

---

## Sumário

Em cinco dias o projeto passou de *"inviável sem HPC, a ~57× do critério"* para:

1. **um resultado positivo**: a fronteira de exclusão de F3 cai de **414 para
   2,93**, certificado — fator 141 (24/09, pelo símbolo pontual no L² do toro; a
   versão de 23/09 dava 10,92 por majorantes de bloco). (A primeira versão dizia
   6,05, mas misturava convenções de norma; refeita em 23/09, §5.25 de
   `ALTERNATIVE_ROUTES.md`.)
2. **um dimensionamento honesto** da rota certificada: **~9800×** a caixa
   orçada, e **~360×** somando todas as alavancas nos seus pisos provados;
3. **um mapa de nove rotas fechadas com o piso de cada uma**, mais cinco rotas
   de certificado e três de infraestrutura;
4. **uma nota de método** que explica por que sete tentativas renderam 1–3× e
   quatro renderam ordens de grandeza.

Nada disso conta os modos instáveis. `T2` segue aberta.

---

## 1. O resultado positivo: F3 de 414 para 6,05

> **Correção de 23/09:** a cadeia abaixo misturava convenções de norma e levava a estrutura de blocos de ℓ¹ para ℓ². Foi refeita num único produto interno, com cauda uniforme e sem ponto flutuante: **`R <= 10,92`, exato** (fator 38 contra 414). O número 6,05 abaixo está retirado; a narrativa fica como registro. Ver `ALTERNATIVE_ROUTES.md` §5.25.5.

`F3` pede: *ausência de espectro para `Re s >= R`*, com evidência exigida
"estimativa de energia quantitativa". Estava **aberta**, e o `R` em uso era
**414**, obtido por série de Neumann grosseira (`ALGEBRAIC_TAIL.md`).

### 1.1 A ideia

`spec(T) ⊂ fecho(W(T))`: a borda direita do **campo de valores** dá um `R`
direto. E campo de valores **é** uma estimativa de energia, `<Lu,u>` — ou seja,
exatamente a evidência que a linha pedia.

A dificuldade era que o projeto vive em `l^1` ponderado, herdado de RT, e campo
de valores exige Hilbert. Isso se resolve numa linha: para pesos positivos,
`sum (w_j|x_j|)^2 <= (sum w_j|x_j|)^2`, logo `l^1(w) ⊂ l^2(w)` com norma 1.
Toda autofunção em `l^1` é autofunção em `l^2`, e **excluir espectro em `l^2`
exclui em `l^1`** — a direção que F3 precisa. A recíproca não vale e não é
necessária.

### 1.2 A cadeia

```
R  <=  R_J  +  ||B||_2  <=  0,08885  +  [ 1,511  +  2 · 2,2273 ]  =  6,05
```

| elo | valor | estado |
|---|---:|---|
| `R_J`, parte livre | < 0,08885 | **certificado**, racional, com cauda |
| parte linear `S Gamma1` | 1,511 | RT, publicado |
| máx. das 7 combinações de `Gamma2`, em `l^2` | 2,2273 | grade + Lipschitz exato + cota de arredondamento |

### 1.3 `R_J`: de constante medida a constante derivada

Quatro passos, cada um destravando o seguinte:

- **redução a 1-D**: o autovetor de `lambda_min(Herm(J_mu))` vive **inteiramente
  na componente 3**, com perfil radial fixo; o índice `m` só acompanha a malha.
  Dimensão 1422 vira dimensão efetiva **18**.
- **o pencil**: conjugar por `sqrt(w)` traz irracionais, mas o mesmo número é o
  menor `lambda` de `Herm(WJ)y = lambda W y`, e **`WJ` e `W` são racionais**.
  Prova-se `lambda_min >= lambda0` por `LDL^T` exato de pivôs positivos.
- **a estrutura**: `Herm(W J_mu)` é **semisseparável de posto 1 mais diagonal**,
  o que colapsa o `LDL^T` numa recursão **escalar** de duas linhas, verificada
  exatamente contra a original.
- **a cauda**: essa recursão tem ponto fixo `1/kappa2`, donde
  `pivo_k/(w_k mu(k+1)) -> 1 - 1/kappa2 = 1/5`, e o intervalo `[0,79; 0,80]` é
  invariante para `k>=15`, com checagem exata abaixo disso.

Resultado: `R_J = mu · f(kappa2)`, com `f` monótona, `f -> -1`, zero em
`kappa2 = 1,3204632636`. Script: `certify_free_numerical_range.py`.

### 1.4 `||B||_2`: a substituição fiel

O `choptuik.tex` está no tarball do `.cache`. A estimativa de `Gamma2` de RT
**não é a norma do fundo**: é um **máximo sobre sete combinações explícitas**
das componentes, cada uma na norma `SSPACE`. Sobre `RefA` o máximo é **8,605**,
e o `10,5` usado em `ALGEBRAIC_TAIL` é cota conservadora disso.

A substituição correta troca **cada `Space(·)` da fórmula de RT** pelo bound
`l^2` da multiplicação daquela combinação. Esse bound vem da decomposição
exata, na base de Chebyshev,

```
M_f = Toeplitz + Resto,
Toeplitz[i,j] = f_{|i-j|}/2 + delta_ij f_0/2   (símbolo = f na elipse de Bernstein)
Resto[i,j]    = f_{i+j}/2   - delta_i0 f_j/2   (primeira linha NULA)
```

com `||Resto||_F` em **forma fechada**:
`(1/4)(k2^2/(k2^2-1)) sum_u f_u^2 (k2^u - k2^-u)`, exata sobre os coeficientes
diádicos. Razão `l^2/l^1` **estável em 0,18–0,30** nas sete combinações — sinal
de que a substituição é estrutural. Máximo: **2,2273**.

Scripts: `certify_background_multiplier.py`, `certify_f3_chain.py`.

### 1.5 A estimativa de energia

> Para todo `u` no domínio, `Re <(J_mu+B)u, u> >= -R ||u||^2` com `R <= 6,05`.
> Logo, para `Re s > R`, `Re <L(s)u,u> >= (Re s - R)||u||^2 > 0`, donde
> `||L(s)u|| >= (Re s - R)||u||`: injetividade e cota de resolvente de uma vez.

### 1.6 A ressalva que resta

Tudo é exato menos a avaliação da grade, em `float` com **cota de erro
explícita** de `3,97e-12`, contra folga de Lipschitz de `5,68e-2`. Trocar por
aritmética intervalar muda a décima segunda casa. É a diferença entre erro
provado *por cota* e erro provado *por construção*.

---

## 2. O mapa negativo

### 2.1 Nove alavancas, todas entre 1,0× e 2,7×

| alavanca | ganho | piso que a fecha |
|---|---:|---|
| pré-condicionador de casca | 1,93× | teto 29–32× contra 57× necessários |
| rota analítica k=0 | — | piso de Perron 9,5–12,3 contra alvo 0,6 |
| reponderação por componentes | 1,0009× | já a 0,09% de `rho(M) = 10,3717` |
| família produto `k1^m k2^n` | 1,01× | o 5/4 de RT já é quase ótimo |
| ponto-base `z0` | pior | a cauda piora mais rápido que `||V||` melhora |
| borda do contorno `sigma0` | 1,22× | perfil em U, sem janela limpa |
| Combes–Thomas | — | taxa `theta = d/(2L||B||)`, **proporcional a `d`** |
| reponderação livre por DOF | 2,22× | piso `rho(|B|) = 3,635`; o alvo está **abaixo** dele |
| campos de valores não convexos | — | Cassini exige `||B||·||C|| < 0,0052`; o menor é 8,31. **Falta 1608×** |

E a deflação do modo de gauge (S4), medida corretamente: **2,62×** — dentro da
mesma faixa, não a exceção que se supôs.

### 2.2 Cinco rotas de certificado para a contagem

| rota | por que falha |
|---|---|
| bound uniforme `1/(min|D|-C)` | usa o mínimo num bloco onde `|D|` cresce; **vazio** |
| grade diagonal | `J_mu` **não é diagonalmente dominante** |
| determinante de Fredholm | `Q0 B` não é de classe traço, `sigma_k ~ k^-0,66` |
| `det_2` / Hilbert–Schmidt | cauda `~N^-1/2`, pior que `1/N` |
| grade na fatorização | já em uso — é o `81/(N+1)` |

### 2.3 O que une as falhas

As de enclausuramento falham pela **mesma geometria**: o contorno passa a
**0,0414** da raiz de gauge, e essa distância entra multiplicativamente em
qualquer critério. O critério que sobreviveu a todas as rodadas:

> **Não é "precisar da resolvente perto do contorno". É deixar a resolvente
> multiplicar uma cauda algébrica.**

Detalhe em `NAO_FAZER.md`.

---

## 3. O dimensionamento

A rota certificada (majorante + inversão certificada) foi dimensionada com
constantes medidas, não estimadas.

| cenário | DOFs vs 107×384 |
|---|---:|
| hoje | **9794×** |
| S4+S3 deflacionados (2,71× medido) | 3616× |
| + `beta` no piso de Perron | 540× |
| + absorver `B` | **360×** |

**A rota não tem teto matemático** — na matriz majorante de blocos, `e_TS` não
depende da caixa enquanto `e_ST` e `e_TT` vão a zero com ela. O que existe é
uma caixa exigida, e ela é de 10² a 10⁴ vezes a que já era cara demais.

E a verdade converge muito antes: contando autovalores em truncamentos do
12×36, **a contagem estabiliza em 3 a partir de `n0 = 16`**, contra `n > 61` que
o certificado exige. **O obstáculo não é o objeto, é o certificado.**

---

## 4. O método

Ver `METODO.md`. O padrão, em uma frase:

> Uma sequência de ganhos de 1 a 3× não é falta de engenho: é o sinal de que a
> formulação está no seu ótimo e o ganho tem de vir de **trocá-la**.

Os quatro movimentos que renderam ordens de grandeza mudaram **o que** se
calcula: trocar a representação (pencil racional), checar se o obstáculo é real
(`l^1 ⊂ l^2` grátis), perguntar **onde** o extremo vive (1422 → 18), e olhar a
**forma** da matriz (semisseparável → recursão escalar → `1/5`).

E o que pegou os erros: **nunca** a revisão do raciocínio; sempre uma
verificação numérica barata que o resultado deveria ter passado. Cinco
resultados meus foram retratados em cinco dias — o bound que saiu **menor** que
o valor real, o peso cuja forma denunciou degenerescência, o pico que sumiu ao
refinar. A prática que fica: **depois de um resultado favorável, procurar a
medida mais barata que o refutaria.**

---

## 5. Estado das obrigações

| id | obrigação | estado |
|---|---|---|
| F3 | ausência de espectro à direita | **`R <= 2,93`, certificado** (§5.26) |
| F2 | coercividade de alta frequência | sondada: limiar `n0 ~ 431` com bound global; apertar exige malha indisponível |
| C1 | aproximação compacta | aberta; a rota substituta pede caixas ≳176×944 |
| C3 | winding 1 no setor A | maquinaria completa; dá **3**, e a redução a 1 depende de S3/S4 |
| S3, S4 | constraints e gauge | parciais; o gerador de gauge é **analítico**, `(Z+1)omega` em `s=mu` |
| T2 | exatamente um modo físico instável | **aberta**, depende de todas as acima |

---

## 6. O que fazer a seguir

1. **Aritmética intervalar na grade de F3** — duas horas, remove a última
   ressalva.
2. **Transcrever `R_J` para o Lean** — o `LDL^T` racional é o que o kernel faz
   bem; é a primeira peça do projeto pronta para isso.
3. **S3 e S4**, sabendo que **não são atalho de custo** (2,7×) e sim obrigação
   para `T2`.
4. **Não** voltar a nenhuma das nove alavancas, nem a ajuste de peso, norma,
   ponto-base, contorno ou enclausuramento.

---

## 7. O que este relatório NÃO afirma

- Não contamos os modos instáveis.
- Não fechamos F3 — temos uma cadeia com um elo certificado e três somas
  finitas, e uma ressalva de ponto flutuante.
- Não há nada em Lean.
- **Não há técnica nova.** A decomposição Toeplitz+Hankel+posto 1 é conhecida
  (Olver & Townsend); elipse de Bernstein, campo de valores, semisseparável,
  seção finita e Perron são clássicos. O que é nosso são as **constantes para
  este operador** e o **mapa**. Ver `ESTADO_DA_ARTE.md`.
- O dimensionamento supõe custo proporcional aos DOFs, o que subestima; e
  `||V||` foi extrapolado de quatro pontos, com `p = 0,42`.
