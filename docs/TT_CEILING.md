# Análise de teto do defeito TT (20/09/2026)

Pergunta: a rota do majorante tem piso, como tiveram o pré-condicionador de
casca (teto 29–32× contra 57× necessários) e a rota k=0 (piso de Perron
9,5–12,3 contra alvo 0,6)?

**Resposta: não tem.** O que existe é uma caixa exigida, e ela é muito
sensível a duas alavancas. Script: `scripts/tt_ceiling.py`; relatório
`build/spectrum/tt-ceiling.json`.

## 1. Por que a pergunta é de Perron, e sobre qual matriz

O critério de blocos de `FREQUENCY_BLOCKS` — existe `r>0` com
`sum_i r_i e_ij < r_j` para cada coluna `j` — é, para `E` não negativa,
exatamente `rho(E)<1`. Então a pergunta do teto É uma pergunta de
Perron–Frobenius, só que sobre a **matriz majorante de blocos**, e não sobre
a matriz infinita de DOFs. Isso é o que diferencia este caso do k=0: lá a
reponderação era uma similaridade diagonal sobre uma matriz fixa, e o ínfimo
sobre todos os pesos era `rho(|X|)`, um piso. Aqui um dos blocos **depende da
caixa** e vai a zero com ela.

## 2. Por que não há piso

Com a parametriz acoplada `V` sobre `G=P+S`, os acoplamentos com o exterior
são os da §7 de `FREQUENCY_BLOCKS`:

```
e_TS = beta * t(P)                (NÃO depende de G)
e_ST = ||V|| * (beta+rho) * t(G)
e_TT =         (beta+rho) * t(G)
```

Reduzindo ao critério de dois blocos com `F=P+S`, a condição é
`e_TT<1` e `e_TS*e_ST < (1-a_fin)(1-e_TT)`. Como `e_TS` é constante em `G` e
os outros dois são proporcionais a `t(G) -> 0`, **a condição passa para toda
caixa suficientemente grande.** Nenhum argumento de reponderação fecha esta
rota, porque a quantidade que precisa encolher de fato encolhe.

## 3. A caixa exigida, e como ela escala

Com `||V||=1` (otimista) e `a_fin=1e-9`, núcleo `P=107×384`, cauda
estruturada:

| `beta` | caixa G mínima | DOFs relativos a 107×384 |
|---:|---:|---:|
| 23,000 (uniforme) | 2124×11470 | 592,9× |
| **10,382 (componentes, hoje)** | **556×3001** | **40,6×** |
| 5,000 | 193×1039 | 4,9× |
| 3,941 | 214×768 | 4,0× |
| 3,000 | 104×562 | 1,4× |
| 2,064 | 107×384 | 1,0× |

**DOFs escalam como `beta^2,8`.** É por isso que a redução de 23 para 10,382,
obtida pela norma adaptada por componentes, valeu sozinha um fator 14,6 em
graus de liberdade. A alavanca é extraordinariamente sensível.

O núcleo é a segunda alavanca, mais fraca: `DOFs ~ t(P)^2`, e dobrar o núcleo
de 107×384 para 214×768 vale 2,3×.

## 4. O alvo, agora numérico

| objetivo | `beta` necessário | fator a ganhar sobre hoje |
|---|---:|---:|
| caixa 214×768 (4× o custo já orçado) | <= 3,941 | **2,6×** |
| caixa 107×384 (**custo zero adicional**) | <= 2,064 | **5,0×** |

Ou seja: o item 2 de `PRIORIDADES.md` deixou de ser "melhorar o majorante" e
passou a ser **"reduzir `beta` de 10,382 para 2,1"**. Se isso for conseguido,
a etapa computacional cabe na caixa que já foi medida e orçada, e o custo
extra é zero.

Para calibrar a plausibilidade: `beta` já caiu de 23 para 10,382 (fator 2,2)
com um único movimento, a norma adaptada por componentes. E há folga
conhecida na própria constante — a norma ponderada exata de `B_6x18`, medida,
vale **8,078**, abaixo do bound 10,382 em uso. Parte do fator 5 pode ser
apenas contabilidade mais apertada.

## 5. O que este diagnóstico NÃO diz

- `||V||=1` é otimista; `||V||>=1` na prática, e as caixas acima são piso.
- `a_fin=1e-9` vem do experimento acoplado em 6×18 com núcleo 4×12, validado
  **somente em s=1**. Uniformidade no contorno não está estabelecida.
- Os acoplamentos `P<->T` foram omitidos por serem pequenos pelo argumento de
  cauda do fundo (§4 de `FREQUENCY_BLOCKS`), mas isso exige verificar a
  separação para cada par `(P,G)` usado.
- A tradução DOFs -> horas-núcleo supõe custo proporcional aos DOFs, o que é
  otimista: o custo por coluna também cresce com a caixa.
- Nada aqui é certificado. É dimensionamento, como a §7 de
  `G_PARAMETRIX_FEASIBILITY`, não prova.

## 6. Consequência para a sequência

O item 1 de `PRIORIDADES.md` fecha **positivo**: a rota não tem teto, e o
item 2 está liberado com alvo numérico. Pela primeira vez desde o
pré-condicionador, o trabalho analítico seguinte tem uma meta que se mede —
`beta <= 2,1` — e um retorno conhecido: ele zera o custo computacional em vez
de apenas reduzi-lo.

Reprodução:

```bash
.venv/bin/python scripts/tt_ceiling.py --out build/spectrum/tt-ceiling.json
```

---

# Rodada 2 (20/09/2026): o piso do item 2c, e o núcleo como alavanca livre

Três resultados, dois negativos e um positivo. Script auxiliar embutido em
`scripts/tt_ceiling.py` (`min_box_self_consistent`).

## 7. O item 2c está fechado, NEGATIVO

O bound `beta_ref = max_j sum_i eta_i M_ij / eta_j` é a norma de coluna
ponderada de `M = A_ref + A_linear`, e `min` sobre todos os pesos positivos
`eta` é `rho(M)` por Perron–Frobenius. Medido:

| quantidade | valor |
|---|---:|
| `beta` com os pesos atuais `(916,4096,243,1233)` | 10,3811 |
| `rho(M)` — piso sobre TODA reponderação por componentes | **10,3717** |
| folga disponível | **1,0009×** |

Os pesos que o autovetor numérico propôs já estavam **a 0,09% do ótimo**.
Não há nada a colher em 2c.

Correção de um registro anterior: a frase "há 1,29× de folga só em
contabilidade", baseada na norma ponderada exata de `B_6x18` (8,078) contra
o bound 10,382, estava frouxa. Aquele 8,078 é de uma truncagem FINITA e não
limita o operador infinito — o próprio `COMPONENT_EXTERIOR.md` diz isso. A
folga real da reponderação é 0,09%, não 29%.

## 8. O item 2b vale menos do que se supunha

Absorver `A_linear` **inteiro** na parametriz triangular por componentes — o
melhor caso concebível de 2b, com inversa certificada e sem perda:

| | `beta` | ganho |
|---|---:|---:|
| hoje, `rho(A_ref + A_linear)` | 10,3717 | — |
| 2b completo, `rho(A_ref)` | **8,3780** | 1,24× |

Um fator 1,24 em `beta`, que pela escala `beta^2,8` vale ~1,5× em DOFs. É
real, mas não é decisivo, e custa certificar a inversa de um principal novo.

**O tamanho verdadeiro do pedido.** Para chegar a `beta <= 2,064` seria
preciso absorver cerca de **93% do termo não linear** `2 S Xi Gamma2(RefA,·)`:

| fração de `A_ref` absorvida | `beta` resultante |
|---:|---:|
| 0% | 10,372 |
| 50% | 6,099 |
| 75% | 3,905 |
| 90% | 2,522 |
| ~93% | 2,064 |

A massa está concentrada nas entradas `h2` e `h4` (somas de coluna 21,78 e
17,25 contra 4,98 e 2,75). Zerar o bloco de `h2` sozinho levaria `rho` a
4,265 — é o maior alavanca isolada, e indica onde uma parametriz exterior
deveria morder primeiro.

## 9. Positivo: o núcleo nunca tinha sido otimizado

A tabela da §3 fixou `P=107×384` porque era a caixa já orçada. Isso
**subestima o núcleo e infla a caixa**: `e_TS = beta*t(P)` cai quando `P`
cresce, e `P` está dentro de `G` de qualquer forma. Otimizando `P` junto,
com a folga `G >= 2P` que o argumento de cauda do fundo exige:

| `beta` | núcleo P | caixa G | DOFs vs 107×384 |
|---:|---:|---:|---:|
| 23 (uniforme) | 559×2006 | 1118×4012 | 109,2× |
| **10,382 (hoje)** | **259×929** | **518×1858** | **23,4×** |
| 8,378 (2b completo) | 211×757 | 422×1514 | 15,5× |
| 3,000 | 89×319 | 178×638 | 2,8× |

**A caixa exigida hoje é 23,4× e não 40,6×** — um fator 1,7 de graça, obtido
escolhendo o núcleo em vez de herdá-lo. Não é matemática nova; é uma escolha
que ninguém tinha feito.

Correção: numa primeira passagem desta rodada foi reportado 14,4×. Estava
errado — ali a restrição que mordia era a separação, e foi tomada a caixa do
critério em vez do máximo com `k*P`. O valor correto é 23,4×.

## 10. Estado e próximo passo

| item | estado |
|---|---|
| 1 — teto | fechado POSITIVO: não há piso |
| 2c — reponderar componentes | fechado NEGATIVO: 0,09% de folga |
| 2b — absorver B | quantificado: 1,24× no melhor caso; ~93% de `A_ref` para o alvo |
| novo — otimizar o núcleo | colhido: 40,6× → 23,4× |

O que falta para tornar este dimensionamento honesto, antes de qualquer
decisão de compra: incluir os acoplamentos `P<->T` com os dados reais de
`background-couplings.json` em função da separação (em `k=2` o bound do §4
de `FREQUENCY_BLOCKS` é 0,140, que NÃO é desprezível nas margens em que
estamos), e substituir `||V||=1` por um majorante real da parametriz finita.
Ambos empurram a caixa para CIMA. As 23,4× são piso.

---

# Rodada 3 (20/09/2026): dimensionamento fechado, e o termo dominante trocou

`scripts/tt_dimensioning.py`, `build/spectrum/tt-dimensioning.json`. Fecha os
dois buracos da rodada 2: os acoplamentos `P<->T` passam a usar a cauda medida
do fundo, e `||V||` deixa de ser 1.

## 11. Dois ganhos livres a mais

**A coroa é ADITIVA, não multiplicativa.** `RefA` tem suporte finito: a cauda
`tail_ref(g_m,g_n)` é exatamente 0 em `g=(42,102)`, e além dali só resta o
`epsilon` publicado, dando `e_TP <= 1,93e-5`. Logo a coroa precisa de largura
~(42,102) **somada** ao núcleo, não de um fator sobre ele. A restrição
`G >= 2P` da rodada 2 era conservadora demais.

**A razão de aspecto de 107×384 não é a ótima.** `t=max(15/M, 81/(N+1))` é
minimizada quando `N+1 = 5,4 M`; a razão 3,59 de 107×384 desperdiça 1,5× no
termo radial.

## 12. `||V||` medido, e o dimensionamento fechado

O que faltava era a norma da parametriz finita. Ela está medida:

| fonte | `||V||` |
|---|---:|
| `free-preconditioner-4x12.json`, `finite_inverse_bound` | 79,13 |
| `frequency-blocks-6x18-core4x12.json`, parametriz acoplada | 80,17 e 63,48 |

**`||V|| ~ 80`, não 1.** Com esse valor e `beta=10,372`:

| `beta` | `||V||=1` | `||V||=20` | `||V||=40` | **`||V||=80`** |
|---:|---:|---:|---:|---:|
| 23 (uniforme) | 46,9× | 426,0× | 793,0× | 1508,9× |
| **10,372 (hoje)** | 11,1× | 95,0× | 175,4× | **332,0×** |
| 8,378 (2b completo) | 7,7× | 64,5× | 118,7× | 223,6× |
| 5,000 | 3,4× | 26,1× | 47,5× | 88,7× |
| 3,000 | 1,6× | 11,4× | 20,3× | 37,5× |

**Resposta fechada: ~332× os DOFs de 107×384** — núcleo 1559×8418, caixa
1601×8520, `rho(E)=0,9995`.

**O termo dominante trocou.** `||V||` de 80 para 1 valeria **30×** em DOFs;
`beta` de 10,372 para 5 vale 3,7×. A prioridade analítica não é mais `beta`.

## 13. O ponto-base é uma alavanca morta

Em `z0` grande a série de Neumann converge e `||V||` ganha bound analítico
pequeno (1,34 em `z0=200`, contra 80 medido), mas a cauda piora mais rápido:
`e_TT=(beta+|delta|)t(G)` com `|delta| ~ z0`. Varredura em `z0`:

| `z0` | `||V||` | DOFs |
|---:|---:|---:|
| 0 | 80 (medido) | **melhor** |
| 200 | 1,34 | 3,6× pior |
| 828 | 0,028 | 50× pior |

E o mesmo vale **mesmo supondo** que o argumento de cancelamento de
`STRUCTURED_EXTERIOR` fosse estendido a `z0>0` — testado hipoteticamente,
`z0=0` continua ganhando. Fechado, negativo: não vale a pena reabrir o
ponto-base, nem estender a cauda estruturada com esse objetivo.

## 14. Consequência para o orçamento

332× muda tudo o que foi dito sobre custo. Se o custo escalar com os DOFs — e
ele escala pior —, a conta de nuvem sai de US$ 626–1.066 para a ordem de
**US$ 200 mil a 350 mil**, e a alocação no SDumont sairia de 24–41 mil UA para
**8 a 14 milhões de UA**. Nenhum dos dois é viável.

**O orçamento de R$ 6.900 de `~/possivel_projeto.docx` está obsoleto e não
deve ser submetido como está.** Ele foi calculado para a caixa 107×384, que
o dimensionamento fechado agora mostra ser insuficiente para fechar C1.

## 15. Riscos conhecidos que empurram para PIOR

- `||V||=80` foi medido em `s=1`, no centro. O contorno `Gamma_A` passa a
  0,043 da raiz em `s=0,168311`, e `||H^-1||` cresce perto de raiz. O `||V||`
  **uniforme no contorno** pode ser muito maior que 80.
- `||V||=80` foi medido em caixas minúsculas (4×12, 6×18). Se crescer com a
  caixa, o problema é autorreferente.
- `a_fin` continua validado só em `s=1`.

## 16. Onde atacar agora

Em ordem de retorno medido:

1. **`||V||`** — 30× de prêmio. Duas perguntas baratas antes de qualquer
   teoria: `||V||` cresce com a caixa? (medir em 12×36, que já existe) e
   quanto ele vale perto da raiz em `s=0,168311`?
2. **`beta`** — 3,7× indo a 5, mas 2c está morto e 2b vale só 1,24×; exigiria
   absorver ~75% do termo não linear.
3. Ponto-base, reponderação por componentes, razão de aspecto, separação:
   **todos fechados**, os dois primeiros negativos, os dois últimos colhidos.

---

# Rodada 4 (20/09/2026): `||V||` medido no contorno — a rota morre no custo

`scripts/measure_parametrix_norm.py`. A pergunta era se o `||V||=80` usado no
dimensionamento sobrevive perto da raiz. **Não sobrevive.**

## 17. A medida

`||V||` é `finite_inverse_bound`: a norma da inversa diádica validada, na norma
com pesos RT. Varrendo o centro `s`:

| centro `s` | 4×12 (dim 138) | 6×18 (dim 333) |
|---|---:|---:|
| 1 — onde os 80 foram medidos | **79,13** | **143,65** |
| 1/2 | 634,30 | 674,28 |
| 1/4 | 626,00 | 555,49 |
| **1/8 — a borda de `Gamma_A`** | **201,48** | **3822,06** |
| 0,15 | 193,76 | — |
| 0,168311 (= `mu`) | 216,57 | — |

Dois efeitos, ambos para pior, e eles se compõem:

1. **Proximidade de raiz.** O `s=1` do registro original é o ponto mais LONGE
   de qualquer raiz. Na borda do contorno o valor é 2,5× maior em 4×12 e
   **27× maior em 6×18**, porque ali a raiz do truncamento está a 0,019 da
   borda.
2. **Crescimento com a malha.** No mesmo `s=1`, `||V||` vai de 79,1 (138 DOFs)
   para 143,7 (333 DOFs). Não há sinal de que estabilize.

Isso não é patologia: `V` é a inversa do operador, e procuramos exatamente onde
o operador é singular. A raiz em questão converge para **`s = mu`**
(`mu = 722873400/2^32 = 0,168311`), a 0,0433 da borda `sigma0 = 1/8` — é o
candidato a modo de gauge de S4, e ela **não sai do contorno** com o
refinamento. A geometria é estrutural, não um artefato.

## 18. O que isso faz com o dimensionamento

Mapa `||V||` → caixa, com `beta = 10,372`:

| `||V||` | caixa exigida | DOFs vs 107×384 |
|---:|---:|---:|
| 80 (o que a rodada 3 usou) | 1601×8520 | 332× |
| 144 (`s=1`, 6×18) | 2108×11258 | 578× |
| 201 (borda, 4×12) | 2470×13213 | 794× |
| 626 | 4283×23003 | 2398× |
| 1000 | 5395×29008 | 3809× |

Extrapolando `||V|| ~ C/d` para a geometria convergida (`d = 0,0433`) com a
constante medida em 6×18 (`C = 3822 × 0,019 = 72,6`), dá `||V|| ~ 1700` **e
ainda crescendo com a malha**. A caixa correspondente passa de 5000×29000.

**A rota, como formulada hoje, está morta no custo.** Não por teto matemático
— a rodada 1 mostrou que não há piso —, mas porque a constante que governa o
tamanho da caixa é uma ordem de grandeza pior do que o registro sugeria.

## 19. Agravante de norma, ainda não medido

Os `||V||` acima estão na norma com pesos **RT**. O `beta = 10,372` está na
norma **eta** por componentes. Transportar de uma para a outra custa até
`max(eta)/min(eta) = 4096/243 ~ 16,9`, e o recálculo direto de `V` na norma
`eta` foi iniciado e **interrompido** (`COMPONENT_EXTERIOR.md`), sem resultado.
Se o transporte for necessário em vez do recálculo, some mais uma ordem de
grandeza. Isso precisa ser medido antes de qualquer conclusão final.

## 20. A alavanca que sobra: `sigma0` é uma escolha

A explosão é governada pela distância da borda à raiz em `s = mu`, e
`sigma0 = 1/8` **não é obrigatório** — é a borda esquerda escolhida para
`Gamma_A`. Movê-la para a esquerda afasta o contorno de `mu` (de 0,0433 para
0,15 em `sigma0 = 1/64`) e deveria derrubar `||V||` proporcionalmente.

O preço é analítico, não computacional: a faixa `0 < Re s < 1/8` está
explicitamente listada como **não coberta**, e baixar `sigma0` traz parte dela
para dentro do contorno, exigindo excluir espectro ali. Ou seja, a alavanca
**troca custo por obrigação analítica** — e essa obrigação já está na lista
(é da família de F2/F3).

É a próxima medida, e é barata: `||V||` em `sigma0 = 1/10, 1/16, 1/32, 1/64,
1/256` na malha 6×18, ~40 s por ponto.

---

# Rodada 5 (20/09/2026): a alavanca `sigma0` também está morta

`scripts/measure_contour_border.py`. A hipótese da §20 era que a explosão de
`||V||` fosse governada pela distância da borda à raiz em `s=mu`, e que baixar
`sigma0` a resolvesse. **Não é, e não resolve.**

## 21. A medida, na malha 6×18

| `sigma0` | `||V||` | distância a `mu` |
|---|---:|---:|
| 1/8 (a borda atual) | 3822,1 | 0,0433 |
| **1/10** | **2399,2** | 0,0683 |
| 1/16 | 3847,1 | 0,1058 |
| **1/32** | **39183,2** | 0,1371 |
| 1/64 | 5503,4 | 0,1527 |
| 1/256 | 3413,4 | 0,1644 |

**`||V||` não é monótono na distância a `mu`.** Afastar a borda piora tanto
quanto melhora.

**Correção de 21/09 sobre a explicação, não sobre a conclusão.** Eu li o pico em
`sigma0=1/32` como quase-singularidade real e escrevi que a faixa
`0 < Re s < 1/8` seria povoada. O refinamento refutou: em 8×24 o pico some. O
perfil convergido é um **U limpo**, com mínimo em `s ~ 0,100`:

| `Re s` | 0,020 | 0,050 | **0,100** | 0,125 | 0,140 | 0,155 |
|---|---:|---:|---:|---:|---:|---:|
| `||V||`, 8×24 | 5444 | 2138 | **1325** | 1636 | 2356 | 5185 |

`||V||` sobe para a esquerda (rumo ao eixo imaginário) e para a direita (rumo à
raiz de gauge em 0,1664). A alavanca continua morta, mas por outro motivo: de
`sigma0=1/8` para `1/10` ganha-se **1,24×**, não os 1,6× que a malha 6×18
sugeria, e não porque haja raiz no caminho.

O melhor ponto de toda a varredura é `sigma0 = 1/10` — **1,24× melhor que a
borda atual** na malha 8×24, contra as ordens de grandeza que seriam
necessárias.

## 22. O dimensionamento com os valores realmente medidos

| `||V||` | origem | caixa exigida | DOFs vs 107×384 |
|---:|---|---:|---:|
| 80 | `s=1`, 4×12 — o registro | 1601×8520 | 332× |
| 2399 | `sigma0=1/10`, o melhor medido | 8345×44938 | **9127×** |
| 3822 | `sigma0=1/8`, a borda atual | 10554×56866 | 14607× |
| 39183 | `sigma0=1/32`, o pior | 35368×190862 | 164291× |

E `||V||` ainda cresce com a malha, então estes são pisos.

## 23. Estado das alavancas, fechado

| alavanca | estado | valor |
|---|---|---:|
| teto matemático da rota | **não existe** | — |
| reponderação por componentes (2c) | fechada, negativa | 1,0009× |
| ponto-base `z0` | fechada, negativa | pior |
| borda do contorno `sigma0` | **fechada, negativa** | 1,6× |
| absorver `B` (2b) | quantificada | 1,24× |
| separação aditiva | colhida | ~1,7× |
| razão de aspecto 5,4 | colhida | ~1,5× |

**Conclusão honesta.** A rota do majorante não tem teto matemático, e é por
isso que não se pode declará-la refutada. Mas todas as alavancas baratas estão
exauridas, e a caixa exigida com as constantes medidas é da ordem de **10⁴
vezes** a que já era cara demais. A distância até a viabilidade não é de um
fator 5 em `beta`, como a rodada 1 sugeriu: é de três a quatro ordens de
grandeza, e o termo que domina é `||V||`, a condicionamento do operador perto
do próprio espectro que se quer contar.

Qualquer continuação por aqui exige uma ideia **estruturalmente diferente** —
não um ajuste de parâmetro, de peso, de norma ou de contorno. Todos os quatro
foram medidos e fechados.

---

# Rodada 6 (21/09/2026): duas correções de números meus

## 24. O `||V||` da borda era 2,3× pessimista

O dimensionamento das rodadas 3 e 4 usou `||V|| = 3822`, medido em **6×18**. Ali
a raiz de gauge do truncamento está em 0,1440, a **0,019** da borda. Na malha
8×24 ela já está em 0,1664, a **0,0414** — e o valor convergido é 0,168311, a
0,0433. Ou seja, 6×18 media uma geometria que **não é a convergida**.

| malha | raiz de gauge | distância à borda | `||V||` em `s=1/8` |
|---|---:|---:|---:|
| 4×12 | 0,2318 | 0,1068 | 201 |
| 6×18 | 0,1440 | 0,0190 | 3822 |
| **8×24** | **0,1664** | **0,0414** | **1636** |
| (convergido) | 0,1683 | 0,0433 | ~1570 estimado |

A razão 3822/1636 = 2,34 bate com a razão de distâncias 0,0414/0,019 = 2,18,
confirmando o modelo `||V|| ~ C/d` e que a diferença é geométrica, não de malha.

**O `||V||` a usar é ~1600, não 3822.** Isso melhora o dimensionamento da rodada
4 por um fator ~2,3 — e não muda a conclusão: continua em milhares de vezes a
caixa orçada.

## 25. A faixa não é povoada; a alavanca `sigma0` morre por outro motivo

Ver §21 corrigida. O perfil convergido é um U limpo com mínimo em `s ~ 0,100`.
A alavanca vale 1,24×, não 1,6×, e a explicação de 20/09 (faixa povoada)
estava errada.

Registro de método: as duas correções vêm da mesma causa — **tratar 6×18 como
se fosse a geometria convergida**. A malha tinha a raiz de gauge 2,2× mais perto
da borda do que o limite, e isso contaminou tanto o `||V||` quanto a leitura do
pico. A lição vale para as rodadas futuras: antes de extrair um número que
dependa da distância a uma raiz, conferir onde a raiz está **naquela** malha.

## 26. `||V||` satura com a malha — a última incerteza estrutural, e ela é favorável

Em `s=1`, ponto mais longe de qualquer raiz, portanto livre do efeito
geométrico da §24:

| malha | DOFs | `||V||` | incremento |
|---|---:|---:|---:|
| 4×12 | 138 | 79,13 | — |
| 6×18 | 333 | 143,65 | +64,52 |
| 8×24 | 612 | 183,49 | +39,84 |

Os incrementos caem com razão **0,617** — geométrico. Extrapolando,
`V_inf ~ 248`, ou **1,35×** o valor de 8×24.

Isso **encerra a preocupação registrada na §15**, de que `||V||` pudesse crescer
com a caixa e tornar o problema autorreferente. Ele não cresce sem limite: a
sequência converge.

Aplicando o mesmo fator ao valor da borda, `||V||` convergido ~**2209**.
Dimensionamento final desta linha de investigação:

| cenário | `||V||` | DOFs vs 107×384 |
|---|---:|---:|
| como estava na rodada 4 (6×18) | 3822 | 14607× |
| medido em 8×24 | 1636 | 6218× |
| **extrapolado a malha infinita** | **2209** | **8402×** |
| + deflação ao vão (item 1) | 648 | 2480× |
| + reponderação no piso (item 2) | 648 | 372× |

**O número honesto para a rota da resolvente é ~8400×**, com todas as alavancas
conhecidas somadas levando a ~370×. As rodadas 3 e 4 diziam 332× e 14607×; as
duas estavam erradas, por motivos opostos, e o valor correto está entre elas.

Ressalva: três pontos são pouco para extrapolar. A desaceleração é clara e
consistente, mas `10x30` e `12x36` confirmariam — a `12x36` custaria ~40 min em
`s=1` na máquina atual.

---

# Rodada 7 (21/09/2026): quatro pontos em 12×36, e a deflação encolhe

Tarefa noturna `scripts/run_overnight_norms.sh`, log
`build/spectrum/normas-20260921.log`. `||V||` na maior malha existente
(dim 1422), quatro centros, ~90 min cada por causa da inversa exata.

| centro `s` | `||V||` |
|---|---:|
| 1 (longe de raiz) | 214,44 |
| **1/8 — a borda** | **1846,19** |
| 0,568 — meio do vão | 692,78 |
| 1/10 | 1518,40 |

## 27. A saturação era 21% otimista

Com o quarto ponto, o ajuste `V = V_inf - c·dim^-p` dá **`V_inf = 299,1`**,
`p = 0,416`, resíduo máximo 3,10:

| malha | dim | `||V||` | incremento |
|---|---:|---:|---:|
| 4×12 | 138 | 79,13 | — |
| 6×18 | 333 | 143,65 | 64,52 |
| 8×24 | 612 | 183,49 | 39,84 |
| **12×36** | **1422** | **214,44** | **30,95** |

A extrapolação de ontem, feita com três pontos, dava 248. **O valor com quatro
é 299 — 21% maior.** A convergência é lenta (`dim^-0,42`), então 299 é ele
próprio um piso: um quinto ponto provavelmente subiria de novo. Não há quinto
ponto possível — 12×36 é o teto dos dados publicados (ver
`ALTERNATIVE_ROUTES.md` §5.6).

`||V||` de borda convergido: **~2575**, contra os 2209 estimados na rodada 6.

## 28. A correção que importa: a deflação vale 2,66×, não 8×

O valor do item 1 é a razão entre o `||V||` da borda e o do meio do vão:

| malha | borda | meio do vão | razão |
|---|---:|---:|---:|
| 6×18 | 3822,1 | 480,5 | **7,95×** |
| **12×36** | **1846,2** | **692,8** | **2,66×** |

A diferença é inteiramente a mesma armadilha da §24: em 6×18 a raiz de gauge
está a 0,019 da borda, contra 0,0433 no limite, o que inflava o numerador.
Com a geometria convergida, **a deflação analítica vale 2,66×** — não as oito
vezes que a rodada de ontem registrou.

Isso rebaixa o item 1 de "a primeira alavanca com sinal positivo" para "mais
uma na faixa de 1 a 3×", junto com as outras.

## 29. Dimensionamento com tudo convergido

| cenário | `beta` | `||V||` | DOFs vs 107×384 |
|---|---:|---:|---:|
| rodada 6, como estava | 10,372 | 2209 | 8402× |
| **corrigido com o 4º ponto** | 10,372 | 2575 | **9800×** |
| + deflação ao vão (item 1) | 10,372 | 966 | 3680× |
| + reponderação no piso (item 2) | 3,635 | 966 | 550× |
| + 2b completo | 2,862 | 966 | **366×** |

O número honesto para a rota da resolvente é **~9800×**, e o melhor caso
concebível somando todas as alavancas nos seus pisos provados é **366×** —
contra os 278× que a rodada de ontem estimava.

## 30. O que se confirmou

A alavanca `sigma0` deu **1,22×** em 12×36, contra 1,24× em 8×24: estável, e
continua morta. O perfil da faixa continua sem raiz.
