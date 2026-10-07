# Problemas em aberto

> **Atualização 03/10/2026:** lista de 22–23/09, anterior a T2. Vários itens foram resolvidos por outras rotas; ver `PANORAMA.md` §5 para o que falta hoje.

Este documento existe para que alguém de fora — ou nós daqui a três meses —
consiga escolher um problema e atacá-lo sem reconstruir o histórico. Cada
abertura traz **o que se sabe**, **o que bloqueia**, **o que contaria como
progresso** e **nossa aposta honesta**.

Antes de começar qualquer coisa, leia [`NAO_FAZER.md`](NAO_FAZER.md). São nove
alavancas e cinco rotas de certificado já fechadas, cada uma com piso medido. O
mapa vale mais que as aberturas: ele economiza meses.

Visão geral do projeto: [`RELATORIO_SETEMBRO_2026.md`](RELATORIO_SETEMBRO_2026.md).

---

## 1. Similaridade não-diagonal — a lacuna no argumento que fecha tudo

**Status: aberta. Descoberta em 22/09/2026; alcance corrigido em 23/09 — leia
"Quanto isso vale" abaixo antes de investir aqui.**

### O que se sabe

Nove alavancas de reponderação estão fechadas por um único argumento: a norma de
coluna ponderada é `||D B D^-1||_1` com `D` diagonal, e por Perron–Frobenius

```
inf  sobre pesos diagonais w>0  de  ||D B D^-1||_1   =   rho(|B|).
```

Como `rho(|B|) = 4,49 > 2,064` (o alvo da rodada 1 de `TT_CEILING.md`), concluiu-se
que "viabilidade por reponderação é impossível, não difícil".

**Essa conclusão está correta e continua valendo — para similaridades
diagonais.** O piso para similaridade **geral** não é `rho(|B|)`, é `rho(B)`:

```
inf  sobre S invertível  de  ||S B S^-1||   =   rho(B),
```

e as duas coincidem apenas se não houver cancelamento. Medimos
(`scripts/measure_similarity_floor.py`):

| malha | `rho(B)` | `rho(|B|)` | razão |
|---|---:|---:|---:|
| 4×12 | 1,2472 | 2,9629 | 2,38× |
| 6×18 | 1,3608 | 3,6352 | 2,67× |
| 8×24 | 1,4091 | 4,0452 | 2,87× |
| 10×30 | 1,4339 | 4,3113 | 3,01× |
| 12×36 | 1,4482 | 4,4931 | 3,10× |

(`rho(|B|)` reproduz os 3,6352 já registrados em 6×18 — a maquinaria confere.)

**O alvo 2,064 está entre os dois pisos.** Acima do diagonal, abaixo do geral.
A direção não está fechada; está fechada só na metade que foi examinada.

### A objeção óbvia, e por que ela não morde

O ínfimo `rho(B)` é aproximado por similaridades cada vez pior condicionadas, e
este projeto já foi mordido exatamente assim uma vez: um peso balanceado por
Perron "melhorou" a norma e tinha condicionamento `10^151`. Então a medida é
obrigatória.

Há um limite rigoroso para **toda** similaridade, via pseudoespectro: se
`||S B S^-1|| = r`, então `rho_eps(B) <= r + eps·cond(S)`, donde
`r >= rho_eps(B) − eps·cond(S)`. Medido:

| `eps` | `rho_eps(B)` em 8×24 | piso que isso impõe |
|---|---:|---|
| 1e-1 | 1,4254 | nenhum |
| 1e-2 | 1,4108 | nenhum |
| 1e-4 | 1,4091 | nenhum |

`rho_eps ≈ rho + 0,2·eps` — o pseudoespectro mal se afasta do espectro. Para um
operador normal seria `rho + eps` exatamente, então **B é bem-comportado nesse
aspecto específico**, apesar de `||P||≈71` na projeção espectral. Nenhum valor de
`eps` produz piso acima de 2,064. **O condicionamento não é a obstrução.**

### O que bloqueia de verdade

Não é matemática, é implementabilidade — e por isso pode ser contornável:

1. **Aritmética exata.** Todo o projeto roda em racionais. `S` precisa ser
   racional e de denominadores controlados, ou o certificado morre.
2. **Estrutura.** `S` densa destrói a esparsidade que torna 1422 DOFs viáveis em
   exato, e transforma `J_mu` — hoje quase diagonal, e a fonte da recursão
   semisseparável que certificou `R_J` — em algo denso.
3. **O majorante inteiro tem de ser rederivado em `S`.** O alvo 2,064 foi
   calculado na formulação fixa, variando só os pesos. Ele não se transplanta de
   graça.

### O que contaria como progresso

- Uma família **estruturada** de similaridades não-diagonais — em blocos, banda
  estreita, ou diagonal mais posto baixo — que seja racional, barata de aplicar,
  e que desça `||S B S^-1||` abaixo de 2,064 em 12×36.
- Ou o contrário: um argumento de que qualquer `S` que preserve a estrutura
  necessária está **de volta** ao piso diagonal. Isso fecharia a abertura de
  verdade, e seria tão valioso quanto abri-la.

### Quanto isso vale — correção de 23/09

O alvo `2,064` foi calculado na rodada 1 de `TT_CEILING.md`, com `||V|| = 1`,
**antes** de `||V||` ser medido. Ali `beta` era o termo dominante e `2,064`
zerava o custo extra. No dimensionamento corrigido (`ALTERNATIVE_ROUTES.md`
§5.22.4) o termo dominante é `||V|| = 949`, e `beta` entra com expoente ≈1,8
(implícito na tabela: 3616× → 540× quando `beta` vai de 10,38 a 3,635).

Levar `beta` até `rho(B) ≈ 1,45` daria então **540× → ~100×** a caixa
107×384 — extrapolado de dois pontos, com expoente que ninguém verificou tão
longe. É a maior alavanca de constante já encontrada no projeto (~5×, contra
1,0–2,7× das outras nove), mas **não abre a rota**: continua sendo uma alavanca
de constante dentro da formulação da resolvente, e `METODO.md` explica por que
esse é o lugar errado para procurar.

Onde ela pode valer mais é fora da rota da resolvente: numa estimativa de
energia, escolher `S` é escolher o produto interno, e o campo de valores **não**
é invariante por similaridade. Em F2 e F3, `||B||` entra diretamente; uma `S`
com `||S B S^-1|| ≈ 1,45` no lugar de `4,45` levaria o limiar de F2 de
`N0 ≈ 431` para algo como `~140` — se `Herm(S J_mu S^-1)` mantiver o
crescimento linear, o que não foi verificado. Ainda acima de 61.

### Aposta honesta

`rho(B)` cresce com a malha (1,247 → 1,448) com incrementos que caem por
aproximadamente metade a cada refinamento, extrapolando para ≈1,47. Isso é
extrapolação, não prova, e para o operador infinito o piso pode ser maior — a
mesma ressalva que `rho(|B|)` carrega. Mas a folga até 2,064 é de 40%, o que
sobrevive a bastante erro de extrapolação.

O ponto 2 é o que nos faria apostar contra: a estrutura semisseparável que deu o
único resultado certificado do projeto é exatamente o que uma `S` densa destrói.
Uma `S` que preserve essa estrutura é um problema bem-posto, e não sabemos
resolvê-lo.

---

## 2. O vão entre 16 e 61 — qual instrumento de certificação o fecha?

**Status: aberta. É a questão mais profunda do projeto.**

Contando autovalores em truncamentos sucessivos, **a contagem estabiliza em 3 a
partir de `n0 = 16`**. O certificado, para garantir o mesmo fato, exige `n > 61`.
O objeto matemático converge quatro vezes antes do instrumento que o certifica.

Todo o custo do projeto — o dimensionamento de 9794× — mora nesse vão. Fechá-lo
não é otimizar; é trocar o instrumento.

**O que já se sabe sobre a forma do vão.** Toda rota fechada compartilha a mesma
estrutura: *inverter algo perto do espectro, e depois limitar uma cauda por série
de Neumann*. O `1/d` da distância ao contorno (`d = 0,0414`) entra
multiplicativamente em todas. A única rota que funcionou — F3 — funcionou porque
**trocou inversão por positividade**.

**O que contaria como progresso.** Um instrumento de contagem que não passe por
inversão. A analogia natural é a lei de inércia de Sylvester: para um problema
auto-adjunto, contam-se autovalores num intervalo pela assinatura dos pivôs de
uma fatoração `LDL^T` — e nós **já temos** uma `LDL^T` exata e certificada,
construída para `R_J`. A obstrução é que o pencil não é auto-adjunto.

A pergunta precisa: *existe uma estrutura (J-auto-adjunção, simetria de
Krein, ou a estrutura hamiltoniana que uma linearização em torno de uma solução
auto-similar costuma ter) que reduza esta contagem a uma contagem de inércia?*

Quem souber responder isso provavelmente não precisa saber nada sobre Choptuik.

---

## 3. Quais outras obrigações admitem a troca inversão → positividade?

**Status: aberta. F3 é a prova de conceito.**

(**Correção de 23/09:** a primeira versão, `R <= 6,05`, misturava convenções. Refeita: **`R <= 2,93`, certificado** (24/09, §5.26; antes 10,92), em `ALTERNATIVE_ROUTES.md` §5.25.5.)

F3 saiu de `R = 414` para `R <= 6,05` porque foi reescrita como estimativa de
energia. O obstáculo aparente — o projeto vive em `l^1` ponderado e campo de
valores exige Hilbert — custou uma linha: para pesos positivos
`sum (w_j|x_j|)^2 <= (sum w_j|x_j|)^2`, logo `l^1(w) ⊂ l^2(w)` com norma 1, e
excluir espectro em `l^2` exclui em `l^1`.

F2 (coercividade) **já é** positividade, e é por isso que cedeu à mesma
maquinaria. A pergunta aberta é sobre **C1** (aproximação compacta), hoje
enunciada como obrigação de inversão/aproximação, e sobre a própria contagem
(item 2).

Detalhes do método em [`METODO.md`](METODO.md).

---

## 4. F2: dois números decidem, e ambos são atacáveis

**Status: sondada, não fechada.**

A sondagem deu `lambda_min(Herm(J_mu)|_{n>N0}) ≈ 0,061·mu·(N0+2)`, que cruza o
bound `||B||_2 = 4,4546` em **`N0 ≈ 431`**. Como a dependência é linear, `N0`
escala inversamente com o coeficiente e diretamente com a constante: apertar
qualquer um dos dois move `N0` proporcionalmente.

**O que já se tentou e não funciona.** Restringir `B` às altas frequências para
obter constante menor. Medimos, e a norma cai devagar e **não converge em malha**
(em `N=16`: 2,61 no 8×24 contra 5,07 no 12×36). O motivo é estrutural: a parte de
alta frequência de um operador de multiplicação é a sua parte de **Toeplitz**,
que é invariante por deslocamento — a norma tende ao sup do símbolo na elipse de
Bernstein, não a zero. Registrado em [`NAO_FAZER.md`](NAO_FAZER.md).

**O que contaria como progresso.** Apertar `4,4546` (é a substituição fiel em
`l^2`, e as sete combinações de RT têm razão `l^2/l^1` entre 0,18 e 0,30 —
pode haver folga), ou melhorar o coeficiente `0,061`, que veio de medida e não de
derivação. A malha acessível (12×36) está muito abaixo de `N0 ≈ 431`, então isto
**não se resolve medindo**: precisa de um argumento.

---

## 5. S3 e S4 — propagação das restrições e classificação de gauge

**Status: parciais. É aqui que mora o risco do projeto.**

São elas que decidem se a contagem finita cai de 3 para 1. As três raízes têm
assinaturas distintas e consistentes com "exatamente um modo físico" — e
consistência não é prova. Ver [`MODE_DISCRIMINATION.md`](MODE_DISCRIMINATION.md).

O gerador de gauge é **analítico**, `(Z+1)omega` em `s = mu`, o que é mais do que
se costuma ter. Ver [`GAUGE_RT_ACTION.md`](GAUGE_RT_ACTION.md) e
[`CONSTRAINT_GAUGE_EQUIVALENCE.md`](CONSTRAINT_GAUGE_EQUIVALENCE.md).

**Aviso de custo, para não se repetir um erro nosso.** Deflacionar o modo de
gauge parecia ser também um atalho de custo. Não é: a raiz de gauge está a
`d = 0,0414` da borda e a seguinte a `0,276`, o que prevê `6,7×` por `1/d` — mas
o ganho **medido** é `2,71×`. A diferença é a não-normalidade da projeção
espectral, `||P|| = 71,3`, que é intrínseca e não some com deflação exata. S3 e
S4 são obrigação para o teorema, não economia.

**Achado de 23/09 que muda o custo de S3.** O operador de constraints `K` é ~130×
mais bem condicionado que `L` no mesmo contorno, e tem uma única raiz em
`Re s > 0`. Excluí-lo fora de um disco em 0,401 prova que todo modo de `L` fora
desse disco satisfaz as constraints — e essa exclusão vive na escala de ~1,7× a
caixa original, não 10⁴×. Ver [`S3_CONDICIONAMENTO_SHARP.md`](S3_CONDICIONAMENTO_SHARP.md).

Esta é a abertura com maior distância entre a nossa competência e o que ela pede:
é análise de EDP e relatividade geral clássicas, não aritmética certificada.

---

## 6. Fechar a ressalva de F3 — a melhor primeira contribuição

**Status: mecânica, cerca de duas horas.**

**Resolvido em 23/09, por outro caminho.** A cadeia foi refeita com Young no ℓ²
simétrico, o que dispensa a avaliação em grade: o certificado novo
(`certify_numerical_range_l2.py`) é inteiramente exato, sem ponto flutuante.
**Feito em 24/09:** `R <= 2,93` pelo símbolo pontual (verdadeiro 2,28). O que
sobrava como abertura era **apertar** a cota (10,92 contra 2,35 verdadeiro),
refazendo o sup na elipse corretamente em 2D. O texto original segue.

A cadeia `R <= 6,05` tem três somas finitas avaliadas em ponto flutuante, com
cota de erro explícita de `3,97e-12` contra folga de `5,68e-2`. Trocar por
aritmética intervalar muda a décima segunda casa decimal e remove a única
ressalva do único resultado original do projeto.

É trabalho bem delimitado, verificável, e não exige entender o resto. Ponto de
entrada: `scripts/certify_f3_chain.py`.

---

## 7. Transcrever `R_J` para o Lean 4

**Status: pronta para começar.**

A fatoração `LDL^T` racional com a recursão escalar de duas linhas é a primeira
peça do projeto pronta para formalização: aritmética exata, recursão finita, e
cauda fechada por intervalo invariante a partir de `k >= 15`. É exatamente o que
o kernel faz bem. Ponto de entrada: `scripts/certify_free_numerical_range.py` e
`formal/`.

---

## Como verificar o que afirmamos

```
make verify-reduced-b          # certificado de winding aceito pelo kernel
python3 scripts/certify_free_numerical_range.py   # R_J < 0,08885, exato
python3 scripts/certify_f3_chain.py               # a cadeia inteira
python3 scripts/measure_similarity_floor.py       # os dois pisos do item 1
```

Cinco resultados nossos foram retratados entre 17 e 22 de setembro, todos por
verificação numérica barata e nenhum por releitura do raciocínio. As retratações
estão no registro, não apagadas — ver [`RETOMADA.md`](RETOMADA.md) e
[`METODO.md`](METODO.md). Se você encontrar a sexta, é contribuição.
