# As constraints discriminam as três raízes instáveis?

`docs/WINDING_DATA.md` §6.5 deixou o programa com uma pergunta precisa: o
winding do truncamento 4×12 sobre `Γ_A` é **3**, enquanto a obrigação C3 pede
**1**. As três raízes são reais:

| `s` | `lambda` |
|---|---|
| 0,709626658 | 2,588174 |
| 0,438698175 | 1,600035 |
| 0,231792504 | 0,845401 |

A hipótese natural é que `0,7096` é o modo físico e que as outras duas violam as
constraints (ou são gauge) — é para isso que S3 e S4 existem. Este documento
testa essa hipótese com o 4×12 que o repositório já tinha (§§2–6) e com um
6×18 gerado para este fim (§8), porque um único truncamento não decide.

**NOTA DE NAVEGAÇÃO.** Este documento cresceu em rodadas (4×12 e 6×18; depois
8×24; depois o teste de gauge; depois o 10×30), e por isso os números **8** e
**9** aparecem duas vezes. As referências externas — o ledger e
`scripts/gauge_projection.py` — apontam para a **segunda** ocorrência: §8 = o
terceiro truncamento (8×24), §9 = o teste de gauge. Dentro da §9, as §§9.6–9.8
são a rodada do **quarto** truncamento (10×30, dimensão 975): §9.6 pré-registra,
§9.7 mede, §9.8 é o veredito. Renumerar exigiria editar o ledger, que está fora
do alcance deste sub-agente.

**AVISO.** Tudo aqui é diagnóstico em ponto flutuante. Nada disto é prova, e
nada disto autoriza descartar raiz alguma da contagem. O entregável é evidência
para orientar S3/S4, com a honestidade de dizer onde ela não decide.

## 1. As três medidas

Produzidas por `scripts/discriminate_modes.py`, todas para o **mesmo**
truncamento (logo comparáveis entre si):

1. **`residuo`** — `‖C(s)v‖/‖v‖` para o autovetor `v` do pencil truncado. É a
   quantidade que `analyze_spectrum.py` imprime (`1,14·10⁻¹` para o candidato no
   4×12).
2. **`q`** — o mesmo resíduo dividido pelo que um vetor genérico produziria,
   `‖C(s)‖_F/√n`. Adimensional: `q ≈ 1` significa "viola tanto quanto um vetor
   aleatório", `q ≪ 1` significa "o autovetor é especial".
3. **`gap_S3`** — `min ‖(A+sI)v‖/‖A+sI‖` sobre `v` unitário **no núcleo de
   `C(s)`**. É a pergunta de S3 posta diretamente: existe vetor que anula o
   pencil *e* satisfaz as constraints? Não depende do condicionamento do
   autovetor calculado. O núcleo tem dimensão exatamente 84 = 138 − 54.

Também se reporta `sigma_min` do sistema empilhado `[L_A(s); C(s)]` e o número
de condicionamento do autovalor — sem este último, um resíduo grande num
autovalor mal condicionado não significaria nada.

## 2. Resultado no 4×12

```
                         s       lambda      residuo          q   sigma_min      gap_S3      cond
+0.709626658-0.000000000i    +2.588174   1.1427e-01     0.0297  7.8695e-04  2.0259e-03  5.99e+00  real
+0.438698175-0.000000000i    +1.600035   3.2321e-01     0.0845  2.2699e-03  3.3027e-03  8.68e+00  real
+0.231792504-0.000000000i    +0.845401   1.1639e-01     0.0306  1.5473e-03  2.5748e-03  1.42e+01  real
+0.714264522-1.056239957i    +2.605090   1.0022e-01     0.0259  1.3937e-03  2.4776e-03  4.34e+00  complexa
+0.369712761-0.861449136i    +1.348429   3.4104e-01     0.0890  2.2657e-03  3.4070e-03  6.47e+00  complexa
-0.019014125-0.000000000i    -0.069349   4.4308e-01     0.1169  4.7970e-04  8.8921e-04  7.03e+01  FASE
```

Linhas de base em `s` genérico, medidas sobre **27 pontos reais a mais de
`0,05` de qualquer raiz** (esse afastamento é necessário: perto de uma raiz o
bloco do pencil já é quase singular sozinho, e a comparação ficaria viciada):

| medida | mínimo | mediana | máximo |
|---|---:|---:|---:|
| `sigma_min` | 8,50·10⁻⁴ | 1,85·10⁻³ | 4,68·10⁻³ |
| `gap_S3` | 1,85·10⁻³ | 2,77·10⁻³ | 5,84·10⁻³ |

Distribuição de `q` sobre os 138 modos: mínimo 0,0225, mediana 0,1314, máximo
0,8246.

## 3. O controle interno que decide a questão metodológica

A última linha da tabela é o ponto central deste documento. O modo de **fase**
(`s ≈ −0,019`, que no problema exato está em `s = 0` e corresponde à fase da
solução DSS) é um modo de gauge **genuíno e conhecido**: ele satisfaz as
constraints por construção. No 4×12 ele tem:

* o **pior resíduo de todos**, `4,43·10⁻¹`, quatro vezes o do candidato físico;
* `q = 0,1169`, no **percentil 44,9** — pior que 45 % de todos os 138 modos.

Um teste que classifica o modo bom conhecido como o pior de todos não é um
teste. Isso estabelece, sem apelo a opinião, que **no 4×12 o resíduo das
constraints não tem poder de discriminação**. A causa é visível na própria
tabela: o modo de fase tem o maior condicionamento (70) e está deslocado de
`0,019` em `s` — ele mesmo mal resolvido nesse truncamento.

Há uma explicação coerente para o resíduo falhar justamente aí: ele é calculado
sobre o autovetor *computado*, e quando o autovalor é mal condicionado esse
vetor não é o objeto certo. As outras duas medidas, que minimizam sobre um
subespaço em vez de usar um vetor específico, não têm esse defeito — e de fato
colocam o modo de fase em primeiro lugar, não em último.

## 4. O que as medidas robustas dizem

Ordenando pelo `gap_S3` (a medida mais fina, que é literalmente o enunciado de
S3) e comparando com a faixa genérica `1,85·10⁻³ – 5,84·10⁻³` (mediana
`2,77·10⁻³`):

| modo | `gap_S3` | posição relativa à faixa genérica |
|---|---:|---|
| fase (gauge conhecido) | 8,89·10⁻⁴ | **abaixo do mínimo dos 27 pontos**, por fator 2,1 |
| `0,7096` (candidato físico) | 2,03·10⁻³ | **dentro** da faixa, perto do limite inferior |
| `0,7143 ± 1,0562i` (translado de Floquet) | 2,48·10⁻³ | **dentro** |
| `0,2318` | 2,57·10⁻³ | **dentro**, junto da mediana |
| `0,4387` | 3,30·10⁻³ | **dentro**, acima da mediana |
| `0,3697 ± 0,8614i` | 3,41·10⁻³ | **dentro** |

Duas leituras, e a segunda é a que importa.

**A medida não é ruído puro: ela detecta o modo que sabemos ser especial.** O
modo de fase — o único da lista que é, com certeza estrutural, compatível com as
constraints — é o único que sai da faixa genérica, por um fator 2,1 abaixo do
mínimo de 27 pontos. Isso calibra a resolução do teste: ele enxerga uma
diferença de fator ~2.

**E, tendo essa resolução, ele não separa as três raízes de `Γ_A`.** As três
estão **dentro** da faixa genérica, e diferem entre si por fatores de 1,3 a 1,6
— abaixo da resolução que o próprio teste demonstrou ter. O candidato físico é o
menor dos três, mas não sai da faixa; `0,2318` cai praticamente na mediana dos
pontos genéricos.

As três medidas ainda **discordam entre si** sobre a ordenação: o modo de fase é
o pior pelo resíduo e o melhor pelo `gap_S3`; o translado de Floquet é o melhor
pelo resíduo e o terceiro pelo `gap_S3`. Discordância entre medidas que deveriam
medir a mesma coisa é a assinatura de que se está no ruído.

## 5. Veredito no 4×12 (parcial — ver §9 da rodada 1, e §§8–9 da rodada 2)

Das três possibilidades que a pergunta admitia:

* **não é** o caso de as constraints separarem limpamente os modos;
* **não está** estabelecido que a hipótese é falsa;
* **é** o caso de que o truncamento 4×12 é grosseiro demais para decidir.

E isso é afirmável com número, não por impressão: o resíduo do próprio candidato
físico cai de `1,14·10⁻¹` (4×12) para `3,61·10⁻³` (8×24), documentado em
`SPECTRAL_PROBLEM.md` §6 — um fator **32**. Ou seja, no 4×12 o resíduo de um
modo genuíno ainda está dominado por erro de truncamento, no mesmo patamar
`10⁻¹` em que estão todos os outros. Um discriminador cujo sinal e cujo ruído
têm a mesma ordem de grandeza não discrimina.

Uma ressalva que corta nos dois sentidos, e que não pode ser omitida: o
`gap_S3` tem resolução demonstrada de fator ~2 (§4), e as três raízes diferem
entre si por menos que isso. Logo o dado é compatível com "as três são
igualmente compatíveis com as constraints" **e** com "há uma diferença menor que
a resolução atual". O 4×12 não separa os dois casos.

## 6. A tensão lógica que o experimento revela

Há um problema com a hipótese antes mesmo dos números, e vale registrá-lo
porque ele reorienta S3 e S4.

**Se S3 for verdadeira, o teste do resíduo é vazio por teorema.** A obrigação S3
afirma exatamente que as equações selecionadas *implicam* as constraints, com o
deslocamento `s` incluído — é a versão para o pencil do argumento que RT provam
no problema não linear periódico. Se isso vale, então **toda** raiz de `L_A`
carrega um autovetor que satisfaz as constraints, e nenhum resíduo pode separar
coisa alguma. Ou seja, as duas afirmações

* "as raízes extras violam as constraints" (a hipótese testada aqui), e
* "as equações selecionadas implicam as constraints" (a obrigação S3),

não podem ser ambas verdadeiras. Os cenários são mutuamente exclusivos:

| cenário | consequência para o teste | consequência para a contagem |
|---|---|---|
| **(a)** S3 vale | resíduo/`gap_S3` nunca discriminam | a distância entre 3 e 1 tem de ser fechada por **S4** (gauge) ou pela região de contagem, não pelas constraints |
| **(b)** S3 falha para o pencil deslocado | as constraints são condição extra genuína | C3 precisa ser reenunciado como winding **restrito** (§6.1) |

A evidência do 4×12 — nenhuma discriminação, com o modo de gauge conhecido
saindo pior que todos no resíduo — é fraca evidência a favor de **(a)**. Não a
estabelece: um truncamento grosseiro produziria o mesmo padrão sob (b). A §8.3
retoma esta tabela com o 6×18, que aponta na direção oposta.

### 6.1 O que um winding restrito exigiria formalmente

Se o cenário for (b), a formulação ingênua não funciona, e é importante dizer
por quê. "Restringir `L_A(s)` a `ker C(s)`" não fecha: `ker C(s)` tem dimensão
`84` e `L_A(s)` o leva no espaço inteiro de dimensão `138`. A não-injetividade
de um mapa `84 → 138` não é condição de codimensão um; o conjunto de `s` onde
ela ocorre seria vazio genericamente, e não um conjunto discreto de autovalores.
O que se precisa é de um **pencil reduzido quadrado**:

1. uma família analítica em `s` de projeções `Π_c(s)` sobre o subespaço
   constrangido, e um operador reduzido
   `L̃(s) = Π_c(s) ∘ L_A(s)|_{ker C(s)}` cujos domínio e contradomínio tenham a
   **mesma** dimensão e dependam analiticamente de `s`;
2. **posto constante**: `dim ker C(s)` tem de ser constante em toda a região.
   Um salto de posto é exatamente onde a redução quebra, e certificar
   constância de posto é obrigação nova, não herdada. (No 4×12 o posto numérico
   de `C(s)` foi `54` em todos os pontos testados, com núcleo de dimensão `84`;
   é evidência fraca e só para este truncamento.)
3. refazer **toda** a cadeia C1–C3 para `L̃`, não para `L_A`: o bound de cauda
   `ε_N` teria de majorar a cauda do operador **reduzido**, cuja regularidade
   analítica não é consequência de P2 para `L_A` — a projeção `Π_c(s)` entra na
   estimativa e precisa da sua própria hipótese de regularidade;
4. só então o winding de `det L̃(s)` conta modos físicos, e o pipeline de
   `build_contour.py` passaria a ser alimentado por `L̃`, com o mesmo bound de
   Jacobi, mas com `K` e `δ` recalculados para o operador reduzido.

Nenhum desses quatro itens é barato, e o item 3 é o que efetivamente duplica o
trabalho analítico de C1.

### 6.2 Se o cenário for (a), o problema é de gauge, e ele não fecha sozinho

`SPECTRAL_PROBLEM.md` §4 nomeia dois modos não físicos esperados: a fase da
solução DSS (`s = 0`, removível por condição de fase) e a variação do tempo de
acumulação, que "corresponde geometricamente a `lambda = 1`", isto é
`s = K/(2π) = 0,27418`.

Isso **não basta** para fechar a conta:

* o modo de fase está em `s ≈ −0,019`, **fora** de `Γ_A` (cujo lado esquerdo é
  `Re s = 1/8`), logo não é uma das três;
* a terceira raiz está em `s = 0,2318`, isto é `lambda = 0,8454`, contra `1` do
  modo de acumulação — 15,5 % abaixo. No mesmo truncamento, o candidato físico
  erra 3,2 % (`2,588` contra `2,674`). A identificação é plausível mas não é
  imediata; ela tem uma previsão falsificável (§7);
* **mesmo que a identificação valha, sobra `0,4387` (`lambda = 1,600`) sem
  explicação.** Nenhuma combinação dos modos de gauge já nomeados explica as
  duas raízes extras.

Portanto, olhando só o 4×12, nem (a) nem (b) fecham a distância entre 3 e 1, e o
repositório não tem explicação para a raiz `lambda ≈ 1,600`. A §8 muda isso:
essa é justamente a raiz que o 6×18 reprova nas constraints. O que continua sem
explicação depois da §8 é a **outra**, a de `lambda ≈ 0,5`, que passa no teste
das constraints e mesmo assim não converge de posição.

## 7. Previsões registradas antes de olhar o 6×18

O teste que decide não é o valor absoluto num truncamento, e sim como cada
quantidade **escala** com o refinamento. Estas quatro previsões foram escritas
antes de o 6×18 existir, justamente para que o resultado não possa ser lido de
forma acomodatícia:

* **P1 (calibração).** O candidato físico deve ir de `s = 0,7096` para `0,7312`
  e seu resíduo de `1,14·10⁻¹` para `1,94·10⁻²` — ambos já documentados em
  `SPECTRAL_PROBLEM.md` §6. Se isso não se reproduzir, o pipeline está errado e
  nada mais na tabela vale.
* **P2 (a hipótese em teste).** Se as constraints discriminam (cenário (b)), os
  resíduos de `0,4387` e `0,2318` **não** devem cair pelo mesmo fator ~6 do
  candidato, e seus `gap_S3` devem continuar dentro da faixa genérica enquanto o
  do candidato desce.
* **P3 (identificação do modo de acumulação).** Se `0,2318` é o modo do tempo de
  acumulação, ele deve **subir** na direção de `s = 0,27418` (`lambda → 1`),
  movendo-se ~18 %, muito mais do que os ~3 % do candidato físico.
* **P4 (artefato).** Se `0,4387` é artefato de truncamento, ele deve se mover
  muito, sair do eixo real, ou desaparecer do semiplano instável.

Nota de robustez: as quatro raízes analisadas estão bem isoladas no espectro do
4×12 (vizinho mais próximo a `0,207`, `0,207`, `0,271` e `0,090`), de modo que a
contaminação de autovetores por degenerescência **não** é a explicação do
resíduo ruim do modo de fase. O que resta é erro de truncamento, que é a leitura
da §5.

## 8. O 6×18, e a resposta (rodada 1)

O 6×18 não existia no repositório; foi gerado para este teste
(`./scripts/generate_spectrum.sh 6 18`, dimensão real 333, constraints
135×333). Com dois truncamentos o quadro muda de figura.

| modo | `s` 4×12 | `s` 6×18 | deslocamento | resíduo 4×12 | resíduo 6×18 | fator | `gap_S3` 4×12 | `gap_S3` 6×18 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| físico (`λ≈2,67`) | 0,709627 | 0,731184 | +3,0 % | 1,143·10⁻¹ | 1,937·10⁻² | **5,90** | 2,03·10⁻³ | 2,84·10⁻⁴ |
| segundo (`λ≈1,5`) | 0,438698 | 0,404385 | −7,8 % | 3,232·10⁻¹ | 3,195·10⁻¹ | **1,01** | 3,30·10⁻³ | 9,33·10⁻⁴ |
| terceiro (`λ≈0,5`) | 0,231793 | 0,144011 | **−37,9 %** | 1,164·10⁻¹ | 3,267·10⁻² | **3,56** | 2,57·10⁻³ | 3,11·10⁻⁴ |
| fase (gauge) | −0,019014 | +0,033811 | — | 4,431·10⁻¹ | 4,257·10⁻² | **10,41** | 8,89·10⁻⁴ | 1,45·10⁻⁴ |

Faixa genérica do `gap_S3`: `1,85·10⁻³ – 5,84·10⁻³` no 4×12;
`4,15·10⁻⁴ – 2,50·10⁻³` (mediana `8,63·10⁻⁴`) no 6×18.

**Calibração (P1).** O candidato físico vai a `s = 0,731184` e resíduo
`1,9369·10⁻²`, contra os `0,7311839` e `1,94·10⁻²` documentados em
`SPECTRAL_PROBLEM.md` §6. Reproduz ao dígito: o pipeline está calibrado e o
resto da tabela pode ser lido.

**P2 confirmada, para uma das duas raízes extras.** O modo `λ≈1,5` é o único da
tabela cujo resíduo **não cai**: fator `1,01` contra `3,56`, `5,90` e `10,41`
dos outros três. E é o único cujo `gap_S3` permanece **dentro** da faixa
genérica — em `9,33·10⁻⁴`, praticamente sobre a mediana `8,63·10⁻⁴` —, enquanto
os outros três descem abaixo do mínimo dos 25 pontos genéricos. Essa é a
assinatura de um modo que **viola as constraints**: refinar a malha melhora tudo
menos ele.

**P3 refutada.** O terceiro modo não vai para `λ = 1` (`s = 0,27418`): ele anda
na direção **oposta**, de `0,2318` para `0,1440`, e de `λ = 0,845` para
`λ = 0,525`. A identificação com o modo do tempo de acumulação está **errada**,
ou pelo menos não é o que estes dois truncamentos mostram.

**P4, parcialmente.** O modo `λ≈1,5` não sumiu nem saiu do eixo real; moveu-se
8 %. Não é um artefato que evapora — é uma raiz persistente do pencil que viola
as constraints. A distinção importa: ele não desaparece refinando, ele é
**excluído** pelas constraints.

### 8.1 A decomposição da contagem

As três raízes que o winding conta em `Γ_A` não são três coisas do mesmo tipo:

| raiz | satisfaz constraints? | posição estável? | de quem é o problema |
|---|---|---|---|
| `λ ≈ 2,67` | sim (resíduo cai 5,9×, `gap` abaixo da faixa) | sim (+3,0 %) | é o candidato físico |
| `λ ≈ 1,5` | **não** (resíduo estagnado, `gap` na mediana genérica) | razoável (−7,8 %) | **S3** |
| `λ ≈ 0,5` | sim (resíduo cai 3,6×, `gap` abaixo da faixa) | **não** (−37,9 %) | **S4** ou artefato |

Ou seja, a distância entre 3 e 1 se decompõe em dois problemas **diferentes**, e
nenhum dos dois é do verificador de winding:

* uma raiz que as constraints excluem — exatamente o que S3 existe para
  formalizar;
* uma raiz que as constraints **aceitam** e que ainda assim não parece física,
  porque sua posição não converge. Essa é do território de S4 (gauge) ou do
  controle de artefatos, e continua sem explicação.

### 8.2 Um aviso operacional para C3

O terceiro modo está em `s = 0,1440` no 6×18 e o lado esquerdo de `Γ_A` está em
`σ₀ = 1/8 = 0,125`. Ele está **dentro por 0,019**, e movendo-se para a esquerda
a ~38 % por passo de refinamento. No 8×24 ele plausivelmente sai do contorno.

Isto é: **o winding sobre `Γ_A` não é estável sob refinamento**, não por defeito
do certificado, mas porque a escolha `σ₀ = 1/8` corta uma região por onde uma
raiz está migrando. Qualquer enunciado de C3 tem de fixar `σ₀` com uma
justificativa que o torne imune a isso — tipicamente uma cota inferior provada
para o espectro físico não nulo, que hoje não existe. Enquanto não existir,
"winding sobre `Γ_A`" é um número que depende de onde se pôs a borda.

### 8.3 De volta aos cenários da §6

A §6 opôs dois cenários: **(a)** S3 vale, e então nenhuma raiz do pencil pode
violar as constraints; **(b)** S3 falha para o pencil deslocado, e as
constraints cortam o espectro. O 4×12 sugeria fracamente (a). O 6×18 exibe uma
raiz cujo resíduo **não cai** — o que sob (a) seria impossível para uma raiz
exata.

Mas isso ainda não escolhe entre os cenários, e é importante dizer por quê. Um
resíduo estagnado é compatível com duas coisas diferentes:

* **(i)** o modo `λ≈1,5` é raiz do pencil exato e viola as constraints — então
  S3, como enunciada, é falsa para o pencil, e C3 tem de virar o winding
  restrito da §6.1;
* **(ii)** o modo `λ≈1,5` **não é** raiz do pencil exato: é uma raiz espúria do
  truncamento, e um autovetor espúrio não tem razão alguma para satisfazer as
  constraints. Nesse caso S3 continua de pé e o que falha é a aproximação.

O deslocamento de posição não separa as duas: `−7,8 %` é mais que os `+3,0 %` do
candidato físico e muito menos que os `−37,9 %` do terceiro modo. **O
diagnóstico numérico chega até aqui e para.** Separar (i) de (ii) exige ou mais
truncamentos (o 8×24 diria se a posição converge) ou a prova de S3 — isto é, o
argumento, não o experimento.

Vale notar a consequência prática, que é a mesma nos dois casos: a raiz `λ≈1,5`
não é contada como modo físico. Muda quem tem de justificar isso — S3 no caso
(i), o bound de cauda C1 no caso (ii).

## 9. Veredito da rodada 1 (4×12 e 6×18)

À pergunta "as constraints discriminam as três raízes?", a resposta honesta em
duas partes:

1. **No 4×12, não** — e isso é demonstrável, não opinião: o modo de gauge
   conhecido sai como o pior de todos no resíduo (§3). Quem tivesse rodado só o
   4×12 e concluído alguma coisa teria concluído errado.
2. **No 6×18, sim, e para exatamente uma das duas raízes extras.** O modo
   `λ ≈ 1,5` é reprovado pelas constraints com clareza (resíduo estagnado em
   `3,2·10⁻¹` enquanto os demais caem por fatores de 3,6 a 10,4). O modo
   `λ ≈ 0,5` **passa** no teste das constraints e, ainda assim, não converge de
   posição.

O que **não** está estabelecido, e não deve ser insinuado:

* que o modo `λ ≈ 0,5` é espúrio ou gauge. Ele satisfaz as constraints tão bem
  quanto o candidato físico; o único sinal contra ele é a deriva de posição, e
  dois truncamentos não determinam um limite. Um 8×24 decidiria.
* que o modo `λ ≈ 1,5` está eliminado. Um resíduo que não cai é evidência forte
  de violação das constraints, mas a eliminação da contagem exige S3 provada —
  o enunciado, não o diagnóstico.
* qualquer coisa sobre o operador infinito. Tudo aqui é o truncamento.

O ganho concreto para o programa é que a distância entre 3 e 1 deixou de ser um
buraco único e passou a ser dois problemas nomeados, um por raiz, com evidência
numérica apontando qual obrigação ataca qual (§8.1).

## 10. Previsões pré-registradas para o 8×24

Escritas e commitadas **antes** de a matriz 8×24 existir, pelo mesmo motivo da
§7: impedir leitura acomodatícia depois do fato. As três primeiras são as
perguntas que dois truncamentos não fecham.

* **Q1 — estagnação do `λ≈1,5`.** Se a estagnação do resíduo é padrão e não
  coincidência, o resíduo desse modo fica em `≈3,2·10⁻¹` (fator entre 1,0 e 1,3
  sobre o 6×18) enquanto o físico cai para o `3,61·10⁻³` documentado (fator
  ~5,4). **Refutada** se o resíduo do `λ≈1,5` cair por fator ≳ 3: nesse caso a
  leitura de dois pontos era ruído e a §8 tem de ser reescrita.
* **Q2 — a terceira raiz sai de `Γ_A`.** Ela foi de `0,2318` para `0,1440`,
  razão `0,621` por passo de refinamento. Extrapolando, prevejo
  `s(8×24) ≈ 0,089`, portanto **abaixo** da borda `σ₀ = 1/8 = 0,125`: a raiz sai
  do contorno e o winding de `Γ_A` cai de **3 para 2**. **Refutada** se
  `s > 0,125`.
* **Q3 — `gap_S3`.** O `λ≈1,5` continua **dentro** da faixa genérica enquanto
  físico, terceiro modo e fase descem abaixo do mínimo dos pontos genéricos.
* **Q4 — calibração.** O candidato físico vai a `s ≈ 0,7331` e resíduo
  `3,61·10⁻³`, ambos documentados em `SPECTRAL_PROBLEM.md` §6. Se não bater, o
  pipeline está errado e o resto da tabela não vale.

Previsão adicional, de menor peso: o modo de fase deve continuar oscilando em
torno de `s = 0` com amplitude decrescente (foi `−0,019` e depois `+0,034`).

---

## 8. O terceiro truncamento (8×24, dimensão 612)

> Seção acrescentada pelo orquestrador. A malha `8×24` foi gerada antes de o
> sub-agente que conduzia esta investigação ser interrompido por limite de API;
> as medidas abaixo foram produzidas com o `scripts/discriminate_modes.py` dele,
> sem alteração. Calibração: o resíduo do candidato físico dá
> `3,612412792401e-03`, contra o `3,61e-3` já documentado em
> [`SPECTRAL_PROBLEM.md`](SPECTRAL_PROBLEM.md).

### 8.1 Resíduo das constraints nos três truncamentos

| modo | 4×12 | 6×18 | 8×24 | fatores de queda |
|---|---:|---:|---:|---|
| físico `λ≈2,67` | 1,14e-1 | 1,94e-2 | 3,61e-3 | 5,88 · 5,37 |
| segundo `λ≈1,46` | 3,23e-1 | 3,19e-1 | 3,25e-1 | **1,01 · 0,98** |
| terceiro `λ≈0,61` | 1,16e-1 | 3,27e-2 | 4,24e-3 | 3,55 · **7,71** |
| fase (gauge) | 4,43e-1 | 4,26e-2 | 4,92e-3 | 10,40 · 8,66 |

### 8.2 As previsões pré-registradas

| previsão | resultado |
|---|---|
| o resíduo do `λ≈1,46` continua estagnado | **confirmada.** Três truncamentos, fatores 1,01 e 0,98. Deixou de ser coincidência de dois pontos |
| a terceira raiz sai de `Gamma_A` | **refutada.** Foi `0,2318 -> 0,1440 -> 0,1664`: desceu e **voltou**, e continua dentro por `0,041` |
| o `gap_S3` do `λ≈1,46` fica dentro da faixa genérica | **confirmada.** `4,58e-4` contra mediana genérica `4,52e-4` — praticamente em cima. Os outros três ficam abaixo da mediana; a fase fica abaixo até do mínimo dos 26 pontos genéricos |

O percentil de `q` corrobora: `λ≈1,46` está no percentil **29,2** dos 612 modos,
enquanto físico, terceiro e fase estão em 0,0, 0,2 e 0,3. A separação é de ordem
de grandeza, não de ajuste fino.

### 8.3 O que isto diz, e é incômodo

As constraints discriminam **uma** raiz, não duas. O `λ≈1,46` não as satisfaz e é
o que S3 teria de excluir. Mas o `λ≈0,61` as satisfaz **cada vez melhor** — fator
7,71 no último refinamento, o melhor da tabela depois da fase.

Logo, depois de S3 fazer o seu trabalho, sobram **duas** raízes que satisfazem as
constraints dentro de `Gamma_A`: `λ≈2,67` e `λ≈0,61`. O teorema-alvo diz uma.

Três leituras, e o diagnóstico numérico não escolhe entre elas:

1. **`λ≈0,61` é modo de gauge** — território de S4. A favor: a fase também
   satisfaz as constraints e também é gauge, então satisfazê-las não é atestado
   de ser físico.
2. **`λ≈0,61` é artefato de truncamento**, apesar do resíduo convergente. A favor,
   e é a anomalia mais forte da tabela: **a posição não converge**
   (`0,2318 -> 0,1440 -> 0,1664`, não monótona), enquanto o físico converge
   limpo (`0,7096 -> 0,7312 -> 0,7331`). Um autovalor genuíno deveria convergir
   em posição, e este não converge.
3. **A contagem-alvo precisa de revisão.**

O que **não** se pode fazer é escolher a leitura conveniente. O resíduo
convergente do `λ≈0,61` e a sua posição não convergente estão em tensão, e três
truncamentos não a resolvem. O passo natural é um quarto truncamento, olhando
especificamente para a posição.

## 9. O `λ≈0,61` é modo de gauge? Previsões pré-registradas

Escritas **antes** de qualquer projeção ser calculada, pelo mesmo motivo das §7
e §10.

O teste é uma projeção sobre a direção de gauge da **fase** (translação em
`τ`), construída de forma independente a partir do fundo `ω*` — não do espectro.
Isso é possível porque o artefato de RT traz `RefA.dat`, e o código-fonte fixa a
convenção sem ambiguidade:

* `Field_U_0_m_times_im_unit` multiplica o coeficiente de índice de Fourier `m`
  por `i·m`: é o `∂_τ` do código;
* `MultiField_in_VBold` fixa o espaço de perturbações: componentes 0,1,2 em
  `VTauTwoPiPeriodic` (só `m` **par**), componente 3 em
  `VTauTwoPiAntiPeriodic` (só `m` **ímpar**), componente 0 também `VXiOdd` (só
  `n` ímpar), e `VBasic2` zera a parte imaginária de `m = 0`. Isso confirma, no
  fonte, as paridades que `SPECTRAL_PROBLEM.md` §3 declara.

Logo a direção de fase no espaço de DOFs é `(i·m)·ω*_{d,m,n}`, com `ω*` lido de
`RefA.dat`. Um fator global não importa: o teste é de sobreposição normalizada.

**Previsões.**

* **G1 (calibração, elimina o teste se falhar).** A direção construída tem
  sobreposição `|cos|` ≥ 0,9 com o autovetor medido do modo de fase
  (`s = 1,39·10⁻⁴` no 8×24). Se ficar abaixo de 0,5, a construção está errada e
  **nada** do resto da seção vale. Referência de escala: dois vetores genéricos
  em dimensão 612 têm `|cos| ≈ 1/√612 ≈ 0,040`.
* **G2 (a pergunta).** Se o `λ≈0,61` é modo de gauge do tipo fase, sua
  sobreposição com a direção de fase é **alta e crescente** com o refinamento.
  Se é físico, **baixa** (ordem de `0,04`) e decrescente.
* **G3 (controle).** O candidato físico `λ≈2,67` tem sobreposição baixa e
  decrescente.
* **G4.** Sobre o `λ≈1,46` (o que viola as constraints) a teoria de gauge não
  prevê nada; reporta-se o número sem interpretá-lo.

**Limitação registrada antes do resultado, e ela é séria.** A direção de gauge
da fase **em expoente `s ≠ 0`** não é `∂_τω*`: é a imagem de um campo vetorial
de forma de Floquet `e^{sτ}∂_τ`, que produz `∂_τω*` **mais** termos de ordem `s`
vindos dos pesos tensoriais, e esses termos exigem a ação geométrica dos
difeomorfismos sobre as variáveis de RT, que o repositório não tem escrita. Em
`s = 0` a identificação é exata — por isso a calibração G1 é legítima. Para
`s ≠ 0` a sobreposição com `∂_τω*` é um **proxy**: alta sobreposição é evidência
de caráter gauge; baixa sobreposição é evidência mais fraca contra ele. Se o
resultado cair na faixa ambígua, a conclusão correta é "não decide", e é isso
que será escrito.

### 9.1 A construção é válida (G1 confirmada, por dois caminhos independentes)

**Caminho 1, estrutural, sem tocar no espectro.** `Field_set_VGauged` de RT zera
**exatamente um** DOF: a parte `Re` da componente 3 em `(m=1, n=0)`. Se a direção
construída é mesmo a direção de gauge que essa condição remove, ela tem de ter
peso grande ali. Tem: `0,683`, `0,667` e `0,665` da norma nos três truncamentos
— dois terços da direção concentrados num único DOF entre 138, 333 e 612.
Como subproduto, isso diz que **no arcabouço de RT a liberdade de gauge deste
espaço é unidimensional**: uma condição a remove.

**Caminho 2, espectral.** Sobreposição da direção construída com o autovetor
medido do modo de fase:

| truncamento | `s` do modo de fase | `|cos|` com `∂_τω*` |
|---|---:|---:|
| 4×12 | −0,019014 | 0,3860 |
| 6×18 | +0,033811 | 0,9929 |
| 8×24 | +0,000139 | **0,9999** |

A sobreposição converge para 1 exatamente na mesma medida em que o autovalor
converge para `s = 0`. **G1 pedia ≥ 0,9 e obteve 0,9999.** A construção está
certa e o resto da seção pode ser lido.

Registre-se o que isso implica para trás: no 4×12 o modo rotulado "fase" tinha
`|cos| = 0,386` com a própria direção de fase, enquanto o modo `λ≈0,85` tinha
`0,82`. **No truncamento mais grosseiro o rótulo estava trocado** — mais uma
razão, independente da §3, para não se concluir nada do 4×12 sozinho.

### 9.2 Resultado

Sobreposição `|cos|` de cada autovetor com a direção de gauge da fase. O nível
genérico (dois vetores quaisquer) é `1/√n`: `0,085`, `0,055` e `0,040`.

| modo | 4×12 | 6×18 | 8×24 | × genérico (8×24) | cond (8×24) |
|---|---:|---:|---:|---:|---:|
| físico `λ≈2,67` | 0,1120 | 0,1177 | 0,1171 | 2,9 | 5,7 |
| `λ≈1,46` (viola constraints) | 0,0113 | 0,0159 | 0,0278 | 0,7 | 11,1 |
| **`λ≈0,61`** | **0,8198** | **0,9349** | **0,9176** | **22,7** | 71,3 |
| fase (gauge conhecido) | 0,3860 | 0,9929 | 0,9999 | 24,7 | 299 |

* **G2: o ramo "alta" está satisfeito.** O `λ≈0,61` está `0,92` alinhado com a
  direção de gauge, 23 vezes o nível genérico, estável entre os dois
  truncamentos bons. "Crescente" vale de 4×12 para 6×18 e fica plano depois.
* **G3: parcialmente errada, e o erro é meu.** Previ que o candidato físico
  teria sobreposição "baixa e **decrescente**". Ela é baixa (`0,117`) mas
  **não decresce**: fica plana em valor absoluto, e só cresce relativamente
  porque o nível genérico encolhe com a dimensão. A previsão estava mal
  formulada — o que um modo físico deve ter é sobreposição de ordem genérica,
  não necessariamente decrescente. `0,117` contra `0,040` é ~3×, pequeno mas não
  nulo.
* **G4: o modo que viola as constraints é também o menos gauge de todos**
  (`0,028`, abaixo do nível genérico). Os dois testes, constraints e gauge,
  apontam para modos **diferentes** — o que é uma consistência não trivial: se
  ambos apontassem para o mesmo, um deles seria redundante.

### 9.3 O que nenhuma previsão cobria

Três fatos apareceram sozinhos, e juntos eles mudam a leitura:

1. **Os autovetores da fase e do `λ≈0,61` são quase paralelos entre si**:
   `|cos|` mútuo de `0,157` (4×12), `0,961` (6×18), `0,915` (8×24).
2. **Os dois são, de longe, os pior condicionados**: `cond` `299` e `71` no
   8×24, contra `5,7` e `11` dos outros dois. E o condicionamento do modo de
   fase **cresce monotonicamente** com o refinamento: `70 → 194 → 299`.
   Autovetores que coalescem e condicionamento que diverge é a assinatura de um
   autovalor **defectivo** (bloco de Jordan) sendo aproximado — e
   `SPECTRAL_PROBLEM.md` §4 já avisava que a variação do tempo de acumulação
   "pode surgir como **polo** em vez de autovetor do BVP".
3. **Como direções, os quatro modos são estáveis sob refinamento.** Restringindo
   aos DOFs comuns, `|cos|` entre 6×18 e 8×24: `0,9998` (físico), `0,9995`
   (`λ≈1,46`), `0,9962` (`λ≈0,61`), `0,9923` (fase). Ou seja: o `λ≈0,61` tem
   **direção estável e posição instável** (`0,2318 → 0,1440 → 0,1664`), que é o
   oposto do comportamento de um artefato numérico comum.
4. **A parte não-gauge do `λ≈0,61` também é estável, e isso decide um ponto.**
   Removendo de cada autovetor a componente ao longo da direção de gauge e
   comparando o que sobra entre 6×18 e 8×24:

   | modo | autovetor inteiro | parte ortogonal ao gauge |
   |---|---:|---:|
   | físico `λ≈2,67` | 0,9998 | 0,9998 |
   | `λ≈1,46` | 0,9995 | 0,9995 |
   | **`λ≈0,61`** | 0,9962 | **0,9812** |
   | fase | 0,9923 | **0,3899** |

   O modo de fase é o controle e se comporta como deve: ele é gauge quase puro
   (`0,9999`), então o que resta depois de remover o gauge tem norma `0,014` e
   direção **instável** (`0,39`) — é resíduo numérico, como tinha de ser. Já a
   componente não-gauge do `λ≈0,61`, de norma `0,40`, é uma direção **bem
   definida e convergente** (`0,98`). Ele não é o modo de fase ressurgido com
   ruído em cima.

### 9.4 Veredito (da rodada do 8×24 — **corrigido** pelas §§9.7.1 e 9.7.5)

> Duas afirmações desta subseção não sobreviveram ao 10×30 e à medida dos
> autovetores à esquerda: a leitura "92 % dele é gauge" superestima o que um
> cosseno entre autovetores de operador **não normal** significa (§9.7.5), e a
> posição, aqui dita instável, **converge** com quatro truncamentos (§9.7.1).
> O veredito vigente é o da §9.8.

**O `λ≈0,61` não é uma direção física independente.** `92 %` dele (no sentido
do cosseno) é a direção de gauge da fase, construída do fundo e validada duas
vezes de forma independente. Um modo físico não tem razão para estar 23 vezes
acima do nível genérico de alinhamento com uma direção de gauge — o candidato
físico está em 2,9 e o modo que viola as constraints, em 0,7.

Das três leituras que estavam na mesa, isso as repesa assim:

| leitura | como fica |
|---|---|
| **gauge** | fortemente favorecida pelo alinhamento de `0,92`, estável em dois truncamentos |
| **artefato** | **desfavorecida na forma simples** e apenas parcialmente sustentada na forma "desdobramento de estrutura defectiva em `s = 0`": o condicionamento crescente (`70 → 194 → 299`) e os autovetores coalescendo (`0,915`) apontam nessa direção, mas a componente não-gauge do `λ≈0,61` é estável a `0,98` (§9.3, item 4), o que um simples desdobramento numérico do modo de fase não produziria |
| **contagem-alvo errada** (o `λ≈0,61` seria físico) | fortemente desfavorecida |

As duas primeiras não seriam alternativas excludentes se o `λ≈0,61` fosse o
modo de fase redividido — seria o mesmo fenômeno visto de dois lados. Mas o item
4 da §9.3 separa-as: a parte não-gauge do `λ≈0,61` tem norma `0,40` e direção
**convergente**, enquanto a do modo de fase é ruído. Ou seja, o `λ≈0,61` tem
estrutura própria além do gauge. A leitura que sobra em pé é a **gauge**, na
forma "direção de reparametrização em `s ≠ 0`, com correções de ordem `s`
estruturadas" — que é exatamente o que o proxy pré-registrado não consegue
separar de um modo físico com forte componente de gauge.

**E há um limite que este teste não pode ultrapassar, por construção.** Um campo
vetorial de reparametrização `e^{sτ}∂_τ` existe para **todo** `s`. Portanto
"estar alinhado com a direção de gauge" diz que o modo vive no setor de
reparametrização, mas **não** diz que `s = 0,166` é um valor admissível para um
modo de gauge. Decidir quais difeomorfismos crescentes preservam o domínio, o
cone e a fixação de gauge é exatamente o conteúdo de **S4**, e nenhum
diagnóstico numérico o substitui. `SPECTRAL_PROBLEM.md` §4 já dizia que essa
classe "não deve ser subtraída da contagem por decreto"; esta seção fornece a
evidência de que ela **precisa** ser classificada, e indica em qual modo o
problema se materializa.

**O que fica estabelecido** (como diagnóstico, não como prova):

* a direção de gauge da fase é construtível a partir de `RefA.dat` e está
  validada a `0,9999`;
* o `λ≈0,61` está `0,92` alinhado com ela, com direção estável e posição
  instável;
* o `λ≈1,46` é simultaneamente o que mais viola as constraints e o menos
  alinhado com o gauge — os dois testes são independentes e apontam para modos
  diferentes;
* portanto a contagem `3` em `Γ_A` decompõe-se em **1 físico + 1 de S3 + 1 de
  S4**, com um mecanismo nomeado para cada um.

**O que não fica estabelecido**: que o `λ≈0,61` seja gauge *em vez de* artefato
defectivo; que a componente de `0,40` ortogonal à direção de fase seja apenas
correção de ordem `s` (o proxy pré-registrado); e nada sobre o operador
infinito. Um quarto truncamento decide a primeira dessas três — §9.5.

### 9.5 Quarto truncamento (10×30): previsões pré-registradas

Lançado antes de existir, pelas mesmas razões de sempre. A dimensão é **975**
(contra 612 do 8×24). A pergunta é uma só: **a posição do `λ≈0,61` converge?**

* **H1 (gauge genuíno em `s ≠ 0`).** A posição estabiliza perto de `0,16–0,17`;
  a sequência `0,2318 → 0,1440 → 0,1664` estaria oscilando em torno de um limite
  não nulo.
* **H2 (desdobramento de estrutura defectiva em `s = 0`).** A posição **desce**
  na direção de `0`, e o condicionamento do par continua crescendo.
* **Controle.** O candidato físico fica em `0,7331 ± 0,001` e o modo de fase em
  `|s| < 10⁻³`. Se o controle falhar, o ponto não vale.
* **Acessório.** O alinhamento de gauge do `λ≈0,61` permanece em `≈0,92` e o
  condicionamento do modo de fase continua a crescer (foi `70 → 194 → 299`).

`H1` e `H2` são exclusivas e a medida as separa. Se a posição fizer algo que
nenhuma das duas prevê — subir, por exemplo —, o correto é dizer que não
decidiu, não escolher a leitura mais conveniente depois do fato.

### 9.6 Previsões adicionais pré-registradas para o 10×30

A §9.5 foi escrita **antes** de o 10×30 existir: `MODE_DISCRIMINATION.md` tem
mtime `09:52` e `rt-A-10x30.dat`, `10:02` do mesmo dia. `H1`, `H2` e os
controles valem como estão escritos e não serão reescritos depois do fato.

Esta subseção acrescenta o que a §9.5 não cobre, e foi escrita antes de
qualquer medida ser tomada sobre a matriz de dimensão **975**.

#### 9.6.1 Uma observação estrutural, anterior a qualquer medida

O briefing desta rodada pedia a construção de uma **segunda** direção de gauge,
a da variação do tempo de acumulação (`SPECTRAL_PROBLEM.md` §4, `lambda = 1`,
isto é `s = K/(2π) = 0,27418`). Ela não é uma direção nova, e isso é
demonstrável antes de rodar qualquer coisa.

Uma reparametrização `tau_R -> tau_R + ε f(tau_R)` agindo sobre o fundo produz

```
   δω_d = ε [ f(tau_R) ∂_tau ω*_d + f'(tau_R) (W ω*)_d ],
```

onde `W` é a matriz (constante, desconhecida neste repositório) dos pesos
tensoriais das quatro componentes de RT. Com `f = e^{s tau_R}` — que é
exatamente a forma de Floquet com perfil periódico constante — isso dá

```
   δω = ε e^{s tau_R} [ ∂_tau ω* + s · W ω* ],   ou seja   h = ∂_tau ω* + s·W ω*.
```

Três consequências, todas anteriores ao experimento:

1. **A fase (`s = 0`) e o tempo de acumulação (`s = 0,274`) são dois membros da
   mesma família**, e o termo de ordem zero é **o mesmo vetor** `∂_tau ω*` nos
   dois casos. Não há segunda direção independente a construir; há uma família
   a um parâmetro cujo termo dominante não depende de `s`.
2. **Portanto a projeção sobre `∂_tau ω*` é cega a `s` por construção.** O
   alinhamento de `0,92` do `λ≈0,61` medido na §9.2 diz que o modo vive no
   setor de reparametrização; ele **não pode**, nem em princípio, dizer que
   `s = 0,166` é um valor admissível. Isso reforça — não enfraquece — a
   ressalva já registrada na §9.4: quais difeomorfismos crescentes preservam
   domínio, cone e fixação de gauge é o conteúdo de **S4**.
3. **O teste honesto da leitura "gauge com correções de ordem `s`" é projetar
   sobre o espaço inteiro da família**, não sobre um vetor. Como `W` é
   desconhecida, o espaço que contém **toda** escolha possível de `W` é

   ```
      G = span{ ∂_tau ω* } + span{ E_{d<-e} ω*_e },
   ```

   onde `E_{d<-e} ω*_e` é a componente `e` do fundo colocada no slot `d` do
   espaço de DOFs (as combinações proibidas por paridade caem sozinhas, porque
   os endereços correspondentes não existem na matriz). A versão **diagonal**
   (`W` diagonal, um peso por componente, que é o caso tensorial usual) tem
   `k = 5`; a versão com mistura arbitrária é uma **cota superior** com `k`
   maior. Isto ataca diretamente o item que a §9.4 deixou explicitamente em
   aberto: se a componente de norma `0,40` do `λ≈0,61`, ortogonal a
   `∂_tau ω*`, é ou não "apenas correção de ordem `s`".

#### 9.6.2 Previsões

* **E1 (identificação).** Os quatro modos continuam identificáveis no 975:
  físico em `s ≈ 0,733`, o das constraints em `s ≈ 0,40`, o `λ≈0,61` em algum
  lugar de `[0, 0,25]`, e a fase em `|s| < 10⁻²`. **Refutada** se aparecer uma
  raiz real nova em `Re s > -0,05` ou se uma das quatro sumir — nesse caso o
  pareamento por ordenação usado na §9.2 é inseguro e a tabela tem de ser
  refeita pelo pareamento por sobreposição (que é o que o script passa a fazer).
* **E2 (subespaço de reparametrização).** Se o `λ≈0,61` é direção de
  reparametrização com correções de ordem `s`, sua sobreposição com `G` sobe de
  `0,92` para **≥ 0,98** e a componente residual cai correspondentemente.
  **Refutada** se ficar em `≈0,92`, isto é, se o acréscimo de `G` sobre
  `∂_tau ω*` não explicar praticamente nada da componente de `0,40`.
* **E2' (condição de vacuidade, registrada antes).** `G` tem dimensão pequena
  mas não nula, e o nível genérico sobe de `1/√n = 0,032` para `√(k/n)`
  (`0,072` com `k = 5`). Se o **candidato físico** também subir acima de `0,5`
  na projeção sobre `G`, o teste de subespaço é **vazio** e nada se conclui
  dele — isso será dito, não contornado.
* **E3 (resíduo das constraints, quarto ponto).** O `λ≈1,46` continua estagnado
  (fator entre `0,8` e `1,3` sobre o 8×24) e o físico continua caindo (fator
  `≥ 3`). Uma quebra aqui invalida a §8.1.
* **E4 (estabilidade da parte não-gauge).** Entre 8×24 e 10×30, a parte
  ortogonal a `∂_tau ω*` do `λ≈0,61` mantém `|cos| ≥ 0,95` e a do modo de fase
  permanece ruído (`≤ 0,7`). Se a do `λ≈0,61` desabar, o item 4 da §9.3 —
  que é o que hoje separa a leitura "gauge" da leitura "artefato defectivo" —
  cai junto.

#### 9.6.3 Regra de decisão, fixada antes de olhar

Sobre a posição do `λ≈0,61`, que é a pergunta de `H1` contra `H2`:

| `s(10×30)` | leitura |
|---|---|
| dentro de `[0,15; 0,18]` | `H1`: limite não nulo; compatível com gauge em `s ≠ 0` |
| `< 0,12` e abaixo do valor do 8×24 | `H2`: descida para `s = 0`, desdobramento defectivo |
| qualquer outra coisa — inclusive subir acima de `0,19`, ou cair na faixa `[0,12; 0,15]` | **não decidiu**: registra-se a sequência e diz-se que quatro truncamentos não fixam a posição |

A terceira linha não é uma saída retórica: é a única leitura admissível se a
sequência continuar a oscilar, e está escrita aqui exatamente para que não seja
substituída depois por uma das outras duas.

#### 9.6.4 Dois testes acrescentados no meio da rodada (e o que neles é pós-hoc)

Esta subseção é escrita **depois** de medir a posição no 10×30 e **antes** de
rodar os dois testes que ela descreve. A distinção é registrada porque importa:
a observação que motiva o primeiro teste é **pós-hoc**, e é assim que ela será
tratada.

**A observação pós-hoc.** A posição do `λ≈0,61` no 10×30 é
`s = 0,168310887`. O parâmetro `mu` de RT, lido como diádico **exato** do
cabeçalho `DyadicQ` de `RefA.dat` (`coeff 722873400`, `TwoExp -32`), é

```
   mu_RT = 722873400 · 2^-32 = 0,16830707900226116   (exato em binário)
```

Isto **não** foi previsto por ninguém; foi notado ao olhar o número. Vale
portanto como pista, não como confirmação de previsão, e é assim que entra no
veredito.

* **F1 (o que a pista prevê, e que é falsificável).** Se a identificação
  `s* = mu_RT` for real, um truncamento `12×36` tem de dar `|s - mu_RT| < 10⁻⁷`,
  e a multiplicidade do autovalor perto de `-mu` em `L_A(0)` tem de ser **um**
  — se houver um bloco inteiro de autovalores em `-mu`, a coincidência é
  estrutural do termo `mu·Gamma_1` e não diz nada sobre este modo em
  particular. A contagem de autovalores a menos de `10⁻³` de `-mu` é medida
  agora; o `12×36` fica registrado como previsão para quem rodar.
* **F2 (a fixação de gauge de RT, com controle).** `Field_set_VGauged` zera
  exatamente um DOF: `Re` da componente 3 em `(m=1, n=0)`. O teste é apagar a
  linha **e** a coluna desse DOF da matriz de Galerkin e reexaminar as raízes
  reais.
  * **Controle, que decide se o teste vale**: o modo de **fase** tem de sair do
    espectro (ou afastar-se muito de `s = 0`). Se ele ficar onde está, apagar
    linha e coluna não realiza a fixação de gauge e **nada** se conclui do
    teste.
  * **Previsão**: se o `λ≈0,61` for gauge do tipo que essa condição remove, ele
    sai junto. Se permanecer perto de `0,1683`, então a condição de gauge de RT
    **não** o remove — o que não o torna físico, mas diz que a fixação de gauge
    que S4 precisaria é **outra**, e nomeia o trabalho que falta.
  * **Limitação registrada antes**: apagar linha e coluna é uma escolha
    particular de redução (espaço de teste restrito do mesmo modo que o espaço
    de tentativa). Ela é legítima como diagnóstico **apenas** se o controle
    passar, e mesmo então não é a única redução possível.

### 9.7 O 10×30 (dimensão 975): resultado

Tudo abaixo é **diagnóstico em ponto flutuante**; nada disto é prova nem
autoriza descartar raiz alguma da contagem. Produzido por
`scripts/gauge_projection.py` (agora com `--constraints` e `--gauge-fixed`), em
45 s para os quatro truncamentos.

**Calibração obrigatória, antes de qualquer leitura.** O script reproduz, ao
dígito, os números já publicados nas §8.1 e §9.2: resíduos
`1,143e-1 / 3,232e-1 / 1,164e-1 / 4,431e-1` no 4×12 e
`1,937e-2 / 3,195e-1 / 3,267e-2 / 4,257e-2` no 6×18; sobreposições de fase
`0,1120 / 0,0113 / 0,8198 / 0,3860` e `0,1177 / 0,0159 / 0,9349 / 0,9929`. A
direção de fase construída de `RefA.dat` tem `|cos| = 1,0000` com o autovetor
medido do modo de fase no 10×30 (era `0,9999` no 8×24): **G1 continua
confirmada, agora no quarto truncamento**, e a construção segue validada.

#### 9.7.1 Posição: `H1` confirmada, `H2` refutada

| modo | 4×12 | 6×18 | 8×24 | 10×30 | passo final |
|---|---:|---:|---:|---:|---:|
| físico `λ≈2,67` | 0,709627 | 0,731184 | 0,733091 | **0,733178** | +8,7e-5 |
| `λ≈1,46` | 0,438698 | 0,404385 | 0,401248 | **0,401024** | −2,2e-4 |
| `λ≈0,61` | 0,231793 | 0,144011 | 0,166389 | **0,168311** | +1,9e-3 |
| fase | −0,019014 | +0,033811 | +0,000139 | **−0,000304** | — |

Os controles da §9.5 passam: o físico fica em `0,733178` (pedido:
`0,7331 ± 0,001`; e `lambda = 2,674071` contra `2,674 ± 0,009` de Gundlach) e o
modo de fase em `|s| = 3,0e-4 < 10⁻³`. O ponto vale.

Pela regra de decisão fixada na §9.6.3, `s = 0,168311` cai **dentro** de
`[0,15; 0,18]`: **`H1` confirmada, `H2` refutada**. E a refutação é mais forte
do que a regra pedia — a posição não só ficou na faixa como **convergiu**: os
passos são `−8,8e-2`, `+2,2e-2`, `+1,9e-3`, isto é, redução por fatores `4,0` e
`11,8`, o mesmo padrão do modo físico (`2,2e-2`, `1,9e-3`, `8,7e-5`).

**Isto retira da mesa a leitura "artefato de truncamento" na forma em que ela
estava.** A §8.3 a sustentava num único argumento — "a posição não converge
(`0,2318 → 0,1440 → 0,1664`, não monótona)". Com quatro pontos, a sequência é
não monótona **e convergente**; a não monotonicidade era o 6×18 passando do
ponto. O que era a anomalia mais forte da tabela deixou de existir.

#### 9.7.2 As previsões `E1`–`E4`

| previsão | resultado |
|---|---|
| **E1** (identificação) | **confirmada.** Exatamente quatro raízes reais em `Re s > −0,05` nos quatro truncamentos; o pareamento por sobreposição (não por ordenação) casa `8×24 → 10×30` com `|cos|` `1,0000 / 1,0000 / 0,9999 / 0,9999`, contra segundo-melhor `0,69 / 0,69 / 0,92 / 0,91`. Nenhuma raiz nova, nenhuma perdida |
| **E2** (subespaço `G`) | **sem veredito, por vacuidade pré-registrada.** O `λ≈0,61` vai a `0,9406` (`G_diag`) e `0,9809` (`G_full`), mas o **candidato físico** vai a `0,7109` e `0,7776` — muito acima do limiar `0,5` que a §9.6.2 fixou como condição de vacuidade. O piso do teste é alto (o fundo e os autovetores são ambos objetos suaves de ordem baixa, e se sobrepõem por isso, não por gauge). **O teste não separa e é assim que ele é reportado** |
| **E3** (constraints) | **confirmada.** Físico `3,612e-3 → 5,609e-4` (fator `6,44`); `λ≈1,46` `3,250e-1 → 3,243e-1` (fator `1,002` — **quarto** ponto de estagnação); `λ≈0,61` `4,244e-3 → 8,275e-4` (fator `5,13`); fase `4,922e-3 → 8,387e-4` (fator `5,87`). A §8.1 sobrevive intacta ao quarto truncamento |
| **E4** (parte não-gauge) | **confirmada nos dois ramos.** De 8×24 para 10×30, a parte ortogonal à direção de fase tem `|cos| = 0,9997` no `λ≈0,61` (pedido: `≥ 0,95`) e `0,5197` no modo de fase (pedido: `≤ 0,7`, isto é, ruído). O item 4 da §9.3 fica de pé: o `λ≈0,61` tem estrutura própria além do gauge, e não é o modo de fase com ruído em cima |

Acessórios da §9.5, também confirmados: o alinhamento de gauge do `λ≈0,61`
permanece em `0,9160` (foi `0,9176`), e o condicionamento do modo de fase
continua a crescer, `70 → 194 → 299 → **430**`. O do `λ≈0,61`, ao contrário,
**estabiliza**: `14 → 105 → 71 → 75`. Se os dois estivessem coalescendo num
bloco de Jordan, os dois condicionamentos divergiriam; só um diverge. E a
sobreposição mútua dos dois autovetores fica parada em `0,9154 → 0,9157` — não
tende a `1`. A leitura "desdobramento de estrutura defectiva" (`H2`) perde
também esse apoio.

#### 9.7.3 A observação pós-hoc: `s* = mu_RT`

Registrada como pista na §9.6.4, **não** como previsão. A posição do `λ≈0,61`
não converge para um número qualquer:

| truncamento | `s` | `\|s − mu_RT\|` |
|---|---:|---:|
| 4×12 | 0,231793 | 6,35e-2 |
| 6×18 | 0,144011 | 2,43e-2 |
| 8×24 | 0,166389 | 1,92e-3 |
| 10×30 | **0,168311** | **3,81e-6** |

com `mu_RT = 722873400 · 2^-32 = 0,16830707900226116`, lido como diádico exato
do cabeçalho de `RefA.dat` (não é um valor ajustado: é o parâmetro que RT
gravaram no artefato). A distância cai por fatores `2,6`, `12,7` e `504` —
convergência do mesmo tipo espectral que a do modo físico, e o resíduo final
`3,8e-6` é da ordem do erro de truncamento esperado nesse nível (o passo final
do próprio modo físico é `8,7e-5`).

Dois contra-testes, para não confundir coincidência com estrutura trivial:

* **multiplicidade.** Exatamente **um** autovalor de `L_A(0)` a menos de `10⁻³`
  de `−mu` no 10×30 (zero nos truncamentos anteriores, onde o modo ainda não
  havia chegado lá). Não é um bloco degenerado: o critério `F1` de
  multiplicidade um está satisfeito.
* **não é entrada diagonal.** Nenhuma entrada diagonal da matriz está a menos
  de `10⁻⁴` de `−mu` (a mais próxima está a `0,337`), de modo que o autovalor
  não é um DOF desacoplado. Em compensação — e isto é honesto reportar —
  **há nove entradas exatamente iguais a `+mu`**, todas na diagonal, nos DOFs
  `(comp 2, ·, m par, n = 0)`. Ou seja, `mu` **está** explicitamente no
  operador (via o termo `mu·Gamma_1`), com o sinal oposto ao do autovalor
  encontrado. Um autovalor ligado a `mu` é portanto estruturalmente plausível,
  não milagroso — e continua sem derivação neste repositório.

`SPECTRAL_PROBLEM.md` §1 diz o que `mu_RT` é: o parâmetro que **controla a
reescala dos nulos `u_±`** — uma liberdade de parametrização da construção, não
um dado físico. Um modo instável cujo expoente é exatamente esse parâmetro é
exatamente o que se esperaria de uma direção de reparametrização; mas
"esperaria" não é "provou". A previsão falsificável fica registrada: num
`12×36`, `|s − mu_RT| < 10⁻⁷`.

#### 9.7.4 `F2`: a fixação de gauge de RT **não** valida o teste

Apagando linha e coluna do único DOF que `Field_set_VGauged` zera:

| truncamento | físico | `λ≈1,46` | `λ≈0,61` | fase |
|---|---:|---:|---:|---|
| 4×12 | 0,765394 | 0,433674 | 0,156628 | sai |
| 6×18 | 0,778102 | 0,404610 | 0,212784 | sai |
| 8×24 | 0,779320 | 0,401411 | 0,222024 | sai |
| 10×30 | **0,779375** | **0,401175** | **0,223237** | **sai** |

O controle pré-registrado passa: o modo de fase **sai** do espectro nos quatro
truncamentos (não sobra raiz real perto de `0`). Mas um segundo controle, que
eu não havia registrado e que é decisivo, **falha**: o modo **físico** também se
move, de `0,733178` para `0,779375` — um deslocamento de `4,6e-2`, que é
**convergente** (`5,6e-2 → 4,7e-2 → 4,62e-2 → 4,62e-2`) e portanto não é erro de
truncamento: é outro problema. Uma fixação de gauge legítima não pode mover o
autovalor físico em 6 %.

Logo **`F2` não decide nada sobre o `λ≈0,61`** (que também se move, para
`0,2232`), e o que ele de fato exibe é a advertência da §6.1 em forma numérica:
restringir o espaço de tentativa sem a projeção correspondente no
contradomínio **muda o espectro**, inclusive o do modo que deveria ser
invariante. Apagar uma linha é uma escolha arbitrária de qual equação
descartar.

Dois subprodutos que valem para S4, e que não dependem da validade do teste:

* o problema reduzido tem espectro próprio e **convergente** (`0,7654 →
  0,7781 → 0,7793 → 0,77938`) — é um objeto bem definido, só não é o objeto
  certo;
* mesmo se fosse legítimo, **ele não fecharia a distância entre 3 e 1**: as
  três raízes reduzidas (`0,779`, `0,401`, `0,223`) continuam todas dentro de
  `Re s > 1/8`. A condição de gauge de RT remove o modo de fase, que já estava
  fora de `Gamma_A`.

#### 9.7.5 Uma ressalva que enfraquece a §9.4, e que tem de ser dita

A §9.4 escreveu que "`92 %` dele (no sentido do cosseno) é a direção de gauge da
fase". A medida está certa, mas a força que lhe foi atribuída estava
superestimada, por uma razão que só apareceu agora ao medir os autovetores à
esquerda:

```
   w_i · v_j  (autovetor à esquerda i, à direita j), 10×30:
   fora da diagonal: 2e-17 a 3e-15   -- zero de máquina
```

Isto é biortogonalidade, e vale exatamente porque os quatro autovalores são
distintos. Consequência: **no sentido espectral, o modo `λ≈0,61` contém zero do
modo de fase** — a decomposição de um autovetor na base de autovetores é ele
mesmo. Os `0,916` são o **ângulo geométrico** entre duas direções invariantes
de um operador não normal (24°), não uma "fração de gauge contida no modo". E
ângulos pequenos entre autovetores são também a assinatura genérica de um par
mal condicionado — e este par é justamente o dos dois piores
condicionamentos da tabela (`430` e `75`).

O que sobrevive dessa medida, e é real: o `λ≈0,61` está a 24° da direção de
gauge construída do fundo, enquanto o candidato físico está a 83° e o
`λ≈1,46` a 88°. É uma assimetria grande e estável em três truncamentos. Mas ela
**não** é o argumento quantitativo de "92 % de gauge" que a §9.4 sugeria, e o
veredito abaixo é reescrito com essa correção embutida.

## 9.8 Veredito: o `λ≈0,61` é modo de gauge?

**Não decidiu** — e o "não decidiu" aqui é mais informativo do que era, porque
duas das três leituras mudaram de peso e a pergunta ficou sharp.

**O que ficou estabelecido (como diagnóstico em ponto flutuante, jamais como
prova):**

1. **O `λ≈0,61` é uma raiz convergente do pencil truncado**, não um artefato
   que evapora. Posição converge (`0,2318 → 0,1440 → 0,1664 → 0,168311`, passos
   caindo por fatores 4 e 12), direção converge (`|cos| = 0,9999` entre os dois
   últimos), condicionamento estabiliza (`71 → 75`), resíduo das constraints
   continua caindo (fator `5,1` no último passo). A leitura **"artefato de
   truncamento" está refutada na forma em que era sustentada**, porque o único
   fato que a sustentava — posição não convergente — era um efeito de ter
   apenas três pontos.
2. **A posição limite coincide com `mu_RT` a `3,8e-6`**, isto é, com o
   parâmetro que na construção de RT controla a reescala dos nulos `u_±`. É uma
   coincidência pós-hoc, com multiplicidade um e sem entrada diagonal que a
   explique trivialmente, e com uma previsão falsificável registrada para o
   `12×36`.
3. **O `λ≈1,46` continua sendo o único que viola as constraints**, agora com
   quatro truncamentos (fatores `1,01 · 0,98 · 1,00`), e continua sendo o menos
   alinhado com o gauge. Os dois testes seguem apontando para modos diferentes.

**O que não ficou estabelecido, e por quê — esta é a parte que importa:**

* **Que o `λ≈0,61` seja gauge.** A evidência a favor é de dois tipos e nenhum
  fecha. O ângulo de 24° com `∂_tau ω*` é real e estável, mas a §9.7.5 mostra
  que ele não significa "92 % de gauge": é um ângulo entre autovetores de um
  operador não normal, e o par em questão é o mal condicionado. A coincidência
  com `mu_RT` é mais forte por ser independente da geometria de autovetores,
  mas é pós-hoc e **não tem derivação**: este repositório não contém a ação
  geométrica que diria por que a liberdade parametrizada por `mu` produz um
  autovetor em `s = mu`.
* **Que ele seja físico.** Nada aqui favorece essa leitura, mas nada a refuta:
  um modo físico cujo expoente calhasse de ser `mu_RT` seria um acidente
  notável, e é só isso que se pode dizer.
* **Os dois testes que poderiam ter decidido, não decidiram, e cada um falhou
  num controle pré-registrado:** o teste do subespaço de reparametrização
  `G` é **vazio** (o candidato físico marca `0,71`, acima do limiar `0,5` que a
  §9.6.2 fixou); e a fixação de gauge de RT (`F2`) é **inválida** como redução,
  porque move o autovalor físico em `4,6e-2` de forma convergente (§9.7.4).
  Registrar isso é o ponto: os dois desenhos experimentais tinham controle, os
  controles dispararam, e as conclusões que eles produziriam foram descartadas.

**Consequência para o programa.** A distância entre 3 e 1 em `Gamma_A` continua
a decompor-se em `1 físico + 1 de S3 + 1 de S4`, como a §9.4 dizia, mas o
terceiro item agora tem um **enunciado matemático preciso** em vez de uma
suspeita: *a linearização em torno da solução de RT tem um autovalor em
`s = mu`, e a questão é se ele é gerado pela liberdade que `mu` parametriza (a
reescala dos nulos `u_±`)*. Isso é decidível analiticamente, é exatamente o
conteúdo de **S4**, e é falsificável numericamente (previsão: `12×36` dá
`|s − mu_RT| < 10⁻⁷`).

E a advertência da §8.2 fica **pior**, não melhor: o `λ≈0,61` converge para
`s = 0,1683`, e a borda esquerda de `Gamma_A` está em `sigma_0 = 1/8 = 0,125`.
A raiz **não** sai do contorno com o refinamento — ela se instala a `0,043` de
dentro. O winding sobre `Gamma_A` é **3** no limite, não `3` por acidente de
malha grossa. Qualquer enunciado de C3 tem de lidar com essa raiz por S4, não
esperando que ela vá embora.


---

## 10. O quinto truncamento refuta a previsão registrada

> Seção do orquestrador. A malha `12x36` (dimensão real **1422**) foi gerada
> (2859 s) e analisada com o mesmo caminho das anteriores.

A §9 registrou, antes de olhar, a previsão falsificável: *um `12x36` dá
`|s - mu_RT| < 10^-7`*. **Refutada.**

| truncamento | dim | `s` do modo | `\|s - mu_RT\|` | fator |
|---|---:|---|---:|---:|
| 4x12 | 138 | 0,231793 | 6,35e-2 | — |
| 6x18 | 333 | 0,144011 | 2,43e-2 | 2,6 |
| 8x24 | 612 | 0,166389 | 1,92e-3 | 12,7 |
| 10x30 | 975 | 0,168311 | 3,81e-6 | 504 |
| **12x36** | **1422** | **0,168309292** | **2,21e-6** | **1,7** |

A convergência **estagnou**: depois de um salto de fator 504, o passo seguinte
rende apenas 1,7. O modo está a `2,2·10^-6` de `mu_RT` e não parece estar mais se
aproximando no ritmo anterior.

### O que isso faz com a identificação

**Enfraquece, não mata.** Continua verdade que o modo se instala a poucos
`10^-6` de `mu_RT`, o que é uma coincidência numérica forte e não tem explicação
alternativa no repositório. Mas a leitura "converge **para** `mu_RT`" era
baseada em quatro pontos com convergência acelerando, e o quinto ponto não a
sustenta. Duas leituras compatíveis com os dados, e o diagnóstico não escolhe:

1. o limite é `mu_RT` e o que trava a convergência é outra coisa — precisão da
   diagonalização em dimensão 1422, ou um termo subdominante lento;
2. o limite **não** é `mu_RT`, mas um valor a ~`2·10^-6` dele, e a proximidade é
   aproximada em vez de exata.

Distinguir exigiria ou um sexto truncamento (dimensão ~1950, custo crescendo
rápido) ou, o que seria muito melhor, a **derivação** que liga o autovalor à
liberdade que `mu` parametriza — que continua sendo o que falta desde o começo.

### O modo físico não estagnou

No mesmo truncamento, o candidato físico segue convergindo limpo:
`s = 0,7331796783`, `lambda = 2,6740777376`, resíduo de constraints
`9,17e-05`. A sequência de resíduos é `1,14e-1 -> 1,94e-2 -> 3,61e-3 ->
5,61e-4 -> 9,17e-5`, caindo por ~6x a cada refinamento, sem sinal de piso. O
contraste com a estagnação do modo de `mu` é ele próprio um dado.

## 11. Revisão de 16/09/2026: a estagnação em 2,21e-6 é truncamento

A §10 concluiu que a identificação com `mu_RT` ficava enfraquecida. A
auditoria física de 16/09 a revê. A distância é exatamente
`s_N−mu_RefA = −y^H r/(y^H v)`, e em 12×36 ela se decompõe em truncamento
do gerador exato `(Z+1)ω` (+2,2126e-6) e defeito do fundo (+1,1e-10). Além
disso, \|cos\|((Z+1)ω, autovetor)→1. O "fator 504" do passo anterior era uma
troca de sinal, não convergência.

A identificação numérica com `s=mu` fica **fortalecida**, mas não provada.
A subtração da raiz continua dependendo de S4 (domínio, cone móvel,
multiplicidade). Leitura física: `s=mu` (`lambda≈0,614`) é o representante,
no gauge analítico de RT, da translação do ponto singular; o `lambda=1` do
tempo próprio central não é analítico no cone nessas coordenadas. Ver
[`AUDITORIA_FISICA_16SET.md`](AUDITORIA_FISICA_16SET.md) §4.
