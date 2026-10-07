# A rota do majorante: estado de fechamento

Documento de fechamento, 20/09/2026. Detalhe rodada a rodada em
[`TT_CEILING.md`](TT_CEILING.md); sequenciamento em
[`PRIORIDADES.md`](PRIORIDADES.md).

## Veredito em uma linha

**A rota não tem teto matemático, e por isso não está refutada — mas todas as
alavancas baratas foram medidas e fechadas, e a caixa exigida com as constantes
reais é da ordem de 10⁴ vezes a que já era cara demais.**

## O que se tentou fechar

C1, F2 e F3 se reduzem a um problema só: tornar o majorante do exterior menor
que 1. O critério de blocos de `FREQUENCY_BLOCKS` (existe `r>0` com
`sum_i r_i e_ij < r_j` por coluna) é, para `E` não negativa, exatamente
`rho(E)<1`. A pergunta é portanto de Perron–Frobenius sobre a matriz majorante
de **blocos**.

## Por que não há teto

Na matriz de blocos, `e_TS = beta*t(P)` **não depende da caixa**, enquanto
`e_ST` e `e_TT` são proporcionais a `t(G) -> 0`. A condição passa para toda
caixa suficientemente grande. É a diferença estrutural para a rota k=0, onde a
reponderação era similaridade diagonal sobre uma matriz FIXA e o ínfimo era
`rho(|X|)` — um piso genuíno.

Consequência: nenhum argumento de reponderação fecha esta rota. O que existe é
uma **caixa exigida**, e a questão passou a ser o tamanho dela.

## As sete alavancas, todas medidas

| alavanca | estado | valor |
|---|---|---:|
| teto matemático da rota | **não existe** | — |
| reponderação por componentes | fechada, **negativa** | 1,0009× |
| ponto-base `z0` | fechada, **negativa** | pior |
| borda do contorno `sigma0` | fechada, **negativa** | 1,6× |
| absorver `B` numa parametriz exterior | quantificada | 1,24× |
| separação núcleo/caixa aditiva | colhida | ~1,7× |
| razão de aspecto 5,4 | colhida | ~1,5× |

**Reponderação por componentes.** `beta_ref` é a norma de coluna ponderada de
`M = A_ref + A_linear`, e o mínimo sobre todos os pesos positivos é `rho(M)`.
Medido: pesos atuais 10,3811, piso 10,3717. Os pesos que o autovetor numérico
propôs já estavam a **0,09% do ótimo**.

**Ponto-base.** Em `z0` grande a série de Neumann converge e `||V||` ganha bound
analítico pequeno (1,34 em `z0=200` contra 80 medido), mas a cauda piora mais
rápido porque `|delta| ~ z0`. `z0=0` ganha, e continua ganhando mesmo supondo a
cauda estruturada estendida a `z0>0`.

**Borda do contorno.** Ver a seção seguinte: `||V||` não é monótono na distância
à raiz, e não há janela limpa.

## O termo que domina: `||V||`

O registro usava `||V|| = 80`, medido em `s = 1`. **`s = 1` é o ponto mais longe
de qualquer raiz** — era o melhor caso possível, em dois eixos ao mesmo tempo:

| centro `s` | 4×12 (138 DOFs) | 6×18 (333 DOFs) |
|---|---:|---:|
| 1 — onde os 80 foram medidos | 79,13 | 143,65 |
| **1/8 — a borda de `Gamma_A`** | **201,48** | **3822,06** |

`||V||` cresce com a proximidade da raiz **e** com a malha, e os dois efeitos se
compõem. Isso não é patologia do método: `V` é a inversa do operador, e o
contorno passa perto de onde o operador é singular **por construção**. A raiz
responsável converge para `s = mu = 0,168311` — o candidato a modo de gauge de
S4 — e não sai do contorno com o refinamento.

Caixa exigida com os valores medidos, `beta = 10,372`:

| `||V||` | origem | DOFs vs 107×384 |
|---:|---|---:|
| 80 | `s=1`, 4×12 — o registro | 332× |
| 2399 | `sigma0=1/10`, o melhor medido | **9127×** |
| 3822 | `sigma0=1/8`, a borda atual | 14607× |

E `||V||` ainda cresce com a malha, então estes são pisos.

## Achado espectral colateral, REFUTADO em 21/09

Varrendo `sigma0` na malha 6×18, `||V||` dá **39183 em `sigma0 = 1/32`** — dez
vezes o valor da borda atual, e fora de qualquer tendência monótona. Há uma
quase-singularidade perto de `s ~ 0,03`.

Essa faixa estava listada em `ALGEBRAIC_TAIL.md` e em C4 como **descoberta**, e
nunca tinha sido sondada. Agora há evidência numérica de que não está vazia.
É um dado sobre o espectro, independente da questão de custo, e afeta C3, C4 e
F3: qualquer argumento que precise excluir espectro ali tem de lidar com isto.

**A ressalva disparou: era artefato.** Em 8×24 o pico desaparece e a curva volta
a ser monótona (3568 em `Re s=0,030`, contra 26295 em 6×18). Em 4×12 também não
havia feature. A faixa continua **descoberta**, não povoada. Ver a seção
corrigida em `PROOF_OBLIGATIONS.md`.

## O que NÃO está fechado

- A rota não está refutada. Não há teto; há custo.
- O agravante de norma não foi medido: os `||V||` estão na norma com pesos RT e
  o `beta` na norma `eta` por componentes. O transporte custa até
  `4096/243 ~ 16,9`, e o recálculo direto de `V` na norma `eta` foi iniciado e
  **interrompido** (`COMPONENT_EXTERIOR.md`), sem resultado. Se o transporte for
  necessário, some mais uma ordem de grandeza.
- `a_fin` continua validado só em `s = 1`.
- O dimensionamento supõe custo proporcional aos DOFs, o que subestima.

## Regra que fica

Continuar por esta rota exige uma ideia **estruturalmente diferente** — não
ajuste de peso, de norma, de ponto-base ou de contorno. Os quatro foram medidos
e fechados em 20/09/2026.

Reprodução:

```bash
.venv/bin/python scripts/tt_ceiling.py --out build/spectrum/tt-ceiling.json
.venv/bin/python scripts/tt_dimensioning.py --out build/spectrum/tt-dimensioning.json
.venv/bin/python scripts/measure_parametrix_norm.py 4x12 6x18
.venv/bin/python scripts/measure_contour_border.py
```
