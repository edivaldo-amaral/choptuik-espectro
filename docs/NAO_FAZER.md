# Rotas fechadas: o que não tentar, e por qual constante

Lista consolidada em 22/09/2026. Cada entrada foi **medida**, não estimada, e
traz o piso que a fecha. Ler antes de propor qualquer trabalho novo.

> **Regra geral: nada que ajuste peso, norma, ponto-base, contorno ou
> enclausuramento.** Nove rotas fechadas com piso medido, a última delas por
> fator 1608.

## As nove

| # | rota | ganho medido | piso que a fecha |
|---|---|---:|---|
| 1 | pré-condicionador de casca | 1,93× | teto 29–32× contra 57× necessários; `SHELL_PRECONDITIONER.md` |
| 2 | rota analítica k=0 | — | piso de Perron `rho(|X|)` = 9,5–12,3 contra alvo 0,6 |
| 3 | reponderação por **componentes** | **1,0009×** | os pesos em uso estão a 0,09% de `rho(M) = 10,3717` |
| 4 | família produto `kappa1^m kappa2^n` | **1,01×** | captura só 19% do ganho disponível; o 5/4 de RT já é quase ótimo |
| 5 | ponto-base `z0` | **pior** | em `z0` grande `||V||` cai, mas a cauda piora mais rápido (`|delta| ~ z0`) |
| 6 | borda do contorno `sigma0` | **1,22×** | perfil em U com mínimo em `s ~ 0,100`; sem janela limpa |
| 7 | Combes–Thomas | — | a taxa é `theta = d/(2L||B||)`, **proporcional a `d`**; na borda daria `N ~ 5·10^5` |
| 8 | reponderação **livre por DOF** | **2,22×** | piso de Perron `rho(|B|) = 3,635`; o alvo `beta <= 2,064` está **abaixo do piso** |
| 9 | campos de valores **não convexos** (quadratic / block / pseudo numerical range) | — | Cassini exige `||B||·||C|| < 0,0052`; o menor medido é 8,31. **Falta 1608×** |

## Rotas de certificado para a contagem, também fechadas

| rota | por que falha |
|---|---|
| bound uniforme `1/(min|D| - C)` | usa o mínimo de `|D|` num bloco onde ele cresce; **vazio** em todos os cortes |
| grade diagonal `D(I + D^-1 P)` | `J_mu` **não é diagonalmente dominante** (a coluna `n` tem `n` entradas de tamanho `2 mu n`); `||D^-1P||` fica em 8–9 |
| determinante de Fredholm | `Q0 B` **não é de classe traço**: `sigma_k ~ k^-0,66` |
| `det_2` / Hilbert–Schmidt | cauda HS `~ N^-1/2`, **pior** que o `1/N` da norma de operador |
| grade na fatorização | já está em uso — é o `81/(N+1)` |

## Infraestrutura

| rota | por quê |
|---|---|
| cluster PXE no laboratório | ciclo de trabalho: 20 PCs usados 2 h por vez rendem menos que **uma** máquina ligada 24/7 |
| malhas acima de 12×36 | `RefA.dat` não sustenta o setor expandido; abortam em `Field_shrink` |
| comprar HPC antes de fechar a analítica | o dimensionamento mudou de 332× para 9794× em dois dias; qualquer orçamento feito antes disso nasce errado |

## O padrão que as une

Ver `METODO.md`. Todas as nove tentam melhorar uma **constante dentro de uma
formulação fixa**, e por isso rendem de 1 a 2,7×: a formulação já estava perto
do seu ótimo, porque quem a montou também otimizou parâmetros.

E as de enclausuramento (5, 6, 7, 9) falham todas pela **mesma geometria**: o
contorno passa a `0,0414` da raiz de gauge, e essa distância entra
multiplicativamente em qualquer critério. **A única forma de mexer nisso é
deflacionar o modo de gauge — que é S4.**

## O que continua aberto

| frente | ganho medido | estado |
|---|---:|---|
| **S4**, deflação do modo de gauge | **6,3×** | parcial no ledger; alvo atual |
| **S3**, junto com S4 | **13,5×** | parcial |
| `beta` ao piso de Perron | 3,7× | 2c morto, 2b vale 1,24× |
| cadeia de F3 | — | `R_J` certificado; três somas finitas e três identificações |
| rota 5, decaimento a priori | — | não verificável com os dados publicados; atacável no papel |

---

## Acréscimos de 22–23/09/2026

### Não tentar encolher `B` restringindo às altas frequências

Parece natural: F2 bate `lambda_min(Herm(J_mu)|_{n>N0}) ≈ 0,061·mu·(N0+2)`
contra o bound **global** `||B||_2 = 4,4546`, e usar um bound restrito à cauda
deveria baixar `N0 ≈ 431`. **Não baixa.**

A parte de alta frequência de um operador de multiplicação é a sua parte de
**Toeplitz**, que é invariante por deslocamento: restringi-la a `n > N` dá um
operador unitariamente equivalente a ela mesma. A norma não tende a zero — tende
ao `sup` do símbolo na elipse de Bernstein. O que decai é só o acoplamento de
Hankel, que é a parte pequena.

Medido em `||B[n>N, n>N]||_2`, e o número confirma: a queda é lenta e **não
converge em malha**, porque o que se está medindo é o truncamento e não a cauda.

| `N` | 8×24 | 12×36 |
|---:|---:|---:|
| 0 | 5,97 | 8,80 |
| 12 | 3,36 | 5,97 |
| 16 | 2,61 | 5,07 |
| 20 | 1,68 | 4,20 |

Em `N = 16` o valor quase dobra entre as duas malhas. Qualquer ganho lido aqui é
artefato de ficar sem malha.

**O que sobra para F2:** apertar a constante `4,4546` ou o coeficiente `0,061`.
Ambos movem `N0` proporcionalmente. E nenhum dos dois se resolve medindo, porque
a malha acessível (12×36) está uma ordem de grandeza abaixo de `N0 ≈ 431`.

### Correção: o piso de Perron não fecha a direção inteira

A entrada sobre reponderação diz que a direção está fechada porque
`rho(|B|) = 4,49` está acima do alvo `2,064`. Isso fecha a reponderação
**diagonal**, que é o que Perron–Frobenius cobre. O piso para similaridade geral
é `rho(B) = 1,4482`, **abaixo** do alvo.

Não é convite para voltar a ajustar pesos — essa parte continua fechada, e
continua na lista. É o registro de que a direção não-diagonal nunca foi
examinada. Ver item 1 de [`PROBLEMAS_ABERTOS.md`](PROBLEMAS_ABERTOS.md).
