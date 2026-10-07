# Rotas alternativas: levantamento de 20/09/2026

Escrito depois de `MAJORANT_ROUTE.md` fechar as sete alavancas da rota atual.
Não é uma lista de ideias soltas: cada item diz o que muda estruturalmente, que
evidência existe, e qual é a ressalva.

## O diagnóstico estrutural

O custo vem de **uma** coisa: o método precisa de um bound uniforme de
`||H(s)^-1||` num contorno que, por construção, passa perto do espectro que se
quer contar. As sete alavancas fechadas hoje tentaram melhorar constantes sem
tocar nessa estrutura — e por isso todas deram entre 1,0× e 1,7×.

Uma alternativa de verdade tem de mudar um destes três:

- **(A)** o que o contorno precisa cercar;
- **(B)** se é preciso inverter;
- **(C)** o que se conta.

## 1. Deflação analítica dos modos não físicos ANTES da contagem — (A)

**Hoje:** contar 3 em `Gamma_A`, depois subtrair 2 por S3 e S4.
**Alternativa:** remover os dois modos não físicos do pencil primeiro, e contar
**1** num contorno colocado em região limpa.

**Por que é possível.** O modo em `s=mu` tem gerador **explícito**,
`e^{mu tau}(Z+1)omega`, com a identidade
`DOmega(w)[e^{mu tau}(Z+1)w] = e^{mu tau}(Z+1)Omega(w)` verificada
universalmente e `|cos|` com o autovetor numérico tendendo a 1
(`GAUGE_RT_ACTION.md`, `MODE_DISCRIMINATION.md`). Não é um vetor calculado
numericamente: é conhecido em forma fechada.

**Por que ajuda.** As três raízes convergem para 0,1683 (gauge, S4), 0,4010
(violação de constraint, S3, resíduo estagnado) e 0,7332 (física, resíduo cai).
Removidas as duas primeiras, o contorno só precisa cercar 0,733, e a borda pode
ficar no vão entre 0,404 e 0,731.

**Medido (6×18), não modelado** (`scripts/measure_gap_conditioning.py`):

| borda `s` | `d` à singularidade mais próxima | `||V||` |
|---|---:|---:|
| 1/8 — a borda atual | 0,019 (modo de gauge) | **3822,1** |
| 0,45 | 0,046 (raiz de S3) | 1288,9 |
| 0,50 | 0,096 | 674,3 |
| 0,55 | 0,146 | 506,7 |
| **0,568 — meio do vão** | 0,163 | **480,5** |
| 0,60 | 0,131 (raiz física) | 538,8 |
| 0,65 | 0,081 (raiz física) | 772,8 |

**Ganho medido: 8,0×, e é o máximo possível desta alavanca.** A curva é um U com
mínimo **exatamente no meio do vão** — 480,5 em 0,568, subindo para os dois
lados conforme se aproxima da raiz de S3 (à esquerda) ou da física (à direita).
Isso valida o modelo `||V|| ~ C/d` com as duas raízes vizinhas governando o
condicionamento, e **limita o ganho**: não adianta procurar um ponto melhor no
vão, porque o ótimo já foi medido.

Uma ordem de grandeza é o que a geometria oferece. As três a quatro que faltam
têm de vir de outro lugar.

**Ressalva séria.** Deflação como **pré-condicionador** já foi tentada e
fracassou (`SHELL_PRECONDITIONER.md`: transientes de 10² a 10⁴ pela
não-normalidade). Mas ali o vetor era calculado e servia para acelerar iteração;
aqui é analítico e serve para mudar o pencil. São situações diferentes — isso
precisa ser **verificado, não assumido**.

**Custo.** O trabalho analítico é o mesmo de S3 e S4. Só muda a ordem: antes da
contagem em vez de depois. Não é trabalho extra.

## 2. Otimizar a realização para o problema espectral — constantes

**Correção honesta do que fechei hoje.** A reponderação que está a 0,09% do piso
de Perron é a de **componentes** (quatro pesos). Os pesos por DOF,
`kappa1=65/64` e `kappa2=5/4`, são **herdados da prova de existência de RT**,
escolhidos para aquele problema, e nunca foram otimizados para este. Essa
família é muito maior e continua **aberta**.

**Barato de sondar:** varredura em `(kappa1, kappa2)` medindo simultaneamente a
norma ponderada de `B` e o `||V||` da parametriz. As duas podem mover juntas.

**Ressalva.** Mudar os pesos invalida bounds herdados de RT que dependem deles —
(22), (34)–(38) em `ALGEBRAIC_TAIL.md`. O ganho tem de ser grande o bastante
para pagar a rederivação.

### 2.1 Medido em 20/09/2026 — e o item está fechado

`scripts/measure_weight_floor.py`, `build/spectrum/weight-floor-6x18.json`.
Validação da maquinaria: a norma `eta` calculada aqui dá **8,077804** contra os
**8,077805** registrados em `COMPONENT_EXTERIOR.md`.

A norma de coluna ponderada de `B` sob pesos diagonais `w` é `||D B D^-1||_1`
com `D=diag(w)`; por Perron–Frobenius o ínfimo sobre **todos** os `w>0` é
`rho(|B|)`. Isso limita de uma só vez `kappa1`, `kappa2`, os pesos de componente
e quaisquer pesos livres por DOF.

| quantidade (malha 6×18) | valor |
|---|---:|
| norma `eta` com os pesos RT atuais | 8,0778 |
| melhor da família produto `kappa1^m kappa2^n` | 6,9407 (`kappa1=0,861`, `kappa2=1,284`) |
| **`rho(|B|)` — piso sobre QUALQUER reponderação diagonal** | **3,6352** |

**Dois resultados, e o segundo é o que importa:**

1. **O item 2 como enunciado vale 1,16×.** A família produto captura só **19%**
   do ganho disponível. Otimizar `kappa1` e `kappa2` não é o caminho — e o ótimo
   tem `kappa1=0,861 < 1`, ou seja, despesar os modos de Fourier altos, o que
   ninguém teria adivinhado mas quase não ajuda.
2. **A reponderação livre por DOF vale 2,22×, e esse é o teto de toda a
   direção.** O alvo `beta <= 2,064` que zeraria o custo extra está **abaixo do
   piso**: viabilidade por reponderação é impossível, não difícil.

Ressalvas: `rho(|B|)` foi medido na truncagem 6×18, e para o operador infinito o
piso é **maior** (acrescentar linhas e colunas não negativas não diminui `rho`).
Além disso, parte dos 2,22× exige apertar o bound analítico de 10,372 até o
valor medido de 8,078, que é trabalho separado.

> **Correção de alcance, 22/09/2026.** O piso de Perron fecha a reponderação
> **diagonal**, e só ela. Para similaridade geral o ínfimo de `||S B S^-1||` é
> `rho(B)`, não `rho(|B|)`, e medimos `rho(B) = 1,4482` contra
> `rho(|B|) = 4,4931` em 12×36 — fator 3,10×, crescente com a malha. **O alvo
> 2,064 está entre os dois pisos.** A frase "viabilidade por reponderação é
> impossível, não difícil" continua correta; a frase "e isso fecha a direção"
> não. Ver item 1 de [`PROBLEMAS_ABERTOS.md`](PROBLEMAS_ABERTOS.md) e
> `scripts/measure_similarity_floor.py`.

## 3. Contagem por inércia / positividade em vez de winding — (B)

Já **nomeada** em F2 ("matrizes corretoras e certificado LDL* intervalar") e
nunca desenvolvida. Um certificado de inércia não usa inversa e é estável perto
de raiz — é imune exatamente ao que matou a rota atual.

**Obstáculo:** o pencil não é auto-adjunto, então Sylvester não se aplica
diretamente ao problema de autovalores.

**Uso imediato mesmo assim:** aplicado a `H*H`, que é sempre hermitiana, dá
certificado de **exclusão** — "não há espectro nesta região" —, que é
precisamente o que F3 e a faixa `0 < Re s < 1/8` precisam. E a sondagem de hoje
mostrou que essa faixa **não pode ser suposta vazia**.

**Vantagem de kernel:** positividade em aritmética racional é o que o Lean faz
bem. Contorno, não.

## 4. Exibir o modo físico e excluir o resto — (C)

Em vez de contar: (i) enclosure certificado do modo em 0,733 — já contemplado na
linha S3 ("exige enclosure certificado do autovetor E multiplicidade algébrica
1") — e (ii) excluir todo o resto por argumento global de coercividade ou
energia, que é o conteúdo de F2 e F3. Troca **inversão** por **positividade**.

## O que NÃO tentar

Qualquer ajuste de peso de componente, de norma, de ponto-base ou de contorno
dentro do arcabouço atual. Medidos e fechados em 20/09/2026, todos entre 1,0× e
1,7×. Ver `MAJORANT_ROUTE.md`.

## O que sobra somando tudo, cada alavanca no seu piso

| cenário | `beta` | `||V||` | caixa | DOFs vs 107×384 |
|---|---:|---:|---:|---:|
| hoje | 10,372 | 3822 | 10554×56866 | 14607× |
| só deflação (item 1) | 10,372 | 480 | 3761×20184 | 1848× |
| só reponderação livre (item 2, no piso) | 3,635 | 3822 | 4067×21837 | 2161× |
| **os dois, ambos no melhor caso** | 3,635 | 480 | 1466×7791 | **278×** |
| os dois + 2b completo | 2,862 | 480 | 1200×6355 | 186× |

**Empilhando todas as alavancas conhecidas, cada uma no seu piso provado,
chega-se a ~200–280×** — um ganho de 50 a 80 vezes sobre hoje, e ainda **duas
ordens de grandeza** da caixa já orçada.

O cálculo é otimista: trata as duas otimizações como independentes, quando
mudar os pesos também muda `||V||`.

## Avaliação honesta

**Nenhuma destas fecha a distância sozinha.** A rota atual está três a quatro
ordens de grandeza da viabilidade, e a deflação analítica vale ~10×.

O valor delas é outro: os itens 1, 3 e 4 **compõem** com o trabalho analítico
que já é obrigatório (S3, S4, F2, F3) em vez de competir com ele — fazem a mesma
prova em ordem diferente, e de quebra melhoram o condicionamento. O item 2 é
puro reconhecimento de que fechei uma sub-família e chamei de família.

Ordem sugerida, por custo crescente de descoberta: **2** (varredura de dois
parâmetros, dias), **1** (maior retorno estrutural, e o trabalho já é
obrigatório), **3** (maior retorno de longo prazo, e é o que o kernel prefere).

---

# Item 3 sondado (20–21/09/2026)

`scripts/scan_strip_sigma_min.py`.

## 3.1 A premissa da técnica não se sustenta como esperado

`sigma_min(H) = 1/||V||`. Um certificado de positividade mede **exatamente a
mesma grandeza** que a inversa: trocar winding por inércia **não melhora o
condicionamento**. E o caminho do complemento de Schur para a cauda dá
`lambda_min(finito) > ||acoplamento||^2`, isto é, `||V||*acoplamento < 1` — a
mesma escala da condição de blocos. Não há ganho estrutural aí.

O que sobra de vantagem real é secundário mas não nulo: **positividade em
aritmética racional é o que o kernel do Lean faz bem**, e contorno não. Trocar
`RationalCertificate` por um `LDL*` de `H*H` é um caminho mais curto até o
kernel, não até a viabilidade.

## 3.2 A aplicação pretendida está bloqueada: a faixa não é vazia

A ideia era usar positividade para **excluir** espectro nas regiões descobertas,
começando por `0 < Re s < 1/8`. Um certificado de exclusão só funciona se a
região for de fato vazia. Perfil de `sigma_min` na malha 6×18:

| `Re s` | 0,010 | 0,020 | **0,030** | 0,040 | 0,050 | 0,100 | 0,120 |
|---|---:|---:|---:|---:|---:|---:|---:|
| `||V||` | 4234 | 7229 | **26295** | 16453 | 6466 | 2399 | 3190 |
| `sigma_min` | 2,4e-4 | 1,4e-4 | **3,8e-5** | 6,1e-5 | 1,5e-4 | 4,2e-4 | 3,1e-4 |

Há um mergulho agudo em **`s ~ 0,031`**. O melhor condicionamento de toda a
faixa está em `s = 0,100` (`sigma_min = 4,2e-4`), entre esse mergulho e a raiz
de gauge em 0,144 — o que explica, a posteriori, por que a varredura de
`sigma0` tinha achado 1/10 como melhor ponto.

**Consequência, se o mergulho for uma raiz:** a faixa teria de ser **contada**,
não excluída, e o item 3 perderia sua aplicação principal.

**Não é. Refutado em 21/09 pelo refinamento** — ver 3.3. Em 8×24 o pico
desaparece e a curva volta a ser monótona. A faixa continua **descoberta**, e a
aplicação do item 3 continua **disponível**: um certificado de exclusão ali não
está bloqueado por raiz conhecida. O que permanece é que ninguém provou que a
faixa é vazia.

## 3.3 O teste de convergência, e por que ele importa muito

Na malha 4×12 **não há feature nenhuma**: `||V||` cai monotonicamente de 450,7
em `Re s=0,010` para 205,1 em `0,120`. O mergulho **aparece com o refinamento**.

O projeto já viu esse padrão duas vezes, com desfechos opostos: o "winding 0" do
setor B em 4×12 era artefato de malha (de 6×18 em diante aparece raiz), e a
terceira raiz do setor A apareceu e **ficou**. Só refinar decide.

**Por que importa muito.** O contorno `Gamma_A` tem borda esquerda em
`sigma0=1/8=0,125`, então um modo em `s ~ 0,031` está **fora da região
contada** — e tem `Re s > 0`, ou seja, é candidato a modo instável não contado.
Se sobreviver ao refinamento, isso não é um problema de custo: é uma questão
sobre o próprio enunciado de C3.

**Resultado do teste (21/09): artefato.**

| `Re s` | 4×12 | 6×18 | **8×24** |
|---|---:|---:|---:|
| 0,020 | 432 | 7229 | 5444 |
| **0,030** | 400 | **26295** | **3568** |
| 0,040 | 367 | 16453 | 2660 |
| 0,080 | 260 | 2699 | 1452 |

Em 8×24 a curva volta a ser monótona decrescente. O pico existia **só** em
6×18 — é artefato daquela malha, o desfecho oposto ao da terceira raiz do setor
A. Não há modo instável não contado em `s ~ 0,031`, e nada deve ser construído
sobre o achado de 20/09. `PROOF_OBLIGATIONS.md` e `MAJORANT_ROUTE.md` foram
corrigidos.

Fica a lição de método, que o projeto já tinha aprendido uma vez: **uma feature
que aparece em exatamente uma malha não é um achado, é um candidato.** A
ressalva registrada em 20/09 foi o que permitiu corrigir em um dia.

## 3.4 Uma ferramenta que o repositório nunca menciona

Procurando por "Combes–Thomas", "decaimento exponencial" e equivalentes em
`docs/`: **nenhuma ocorrência**.

A cauda em uso, `t(G) ~ 81/(N+1)`, é **algébrica**, e vem só do crescimento
diagonal da parte livre `J_mu`. Ela não explora que **`B` é bandada** — `RefA`
tem suporte finito, com separação (42,102). Operador com decaimento fora da
diagonal mais parte livre com diagonal crescente é exatamente o cenário das
estimativas tipo Combes–Thomas, que dão decaimento **exponencial** dos blocos
fora da diagonal da resolvente, com taxa governada pela distância ao espectro.

Se aplicável, trocaria `N ~ 10^4` por `N ~ log`, e é a **única coisa vista até
agora capaz de mover ordens de grandeza** — todas as alavancas medidas ficaram
entre 1,0× e 8×.

**Ressalva já registrada no repositório** (`AUDITORIA_MATEMATICA_16SET.md`):
só `RefA` tem suporte finito; a bola `epsilon` do fundo exato **não** tem, e
"cresce ~100 índices radiais por potência". Um argumento Combes–Thomas teria de
tratar a bola à parte. Isso não o inviabiliza — perturbação pequena e não
bandada sobre operador bandado ainda admite decaimento —, mas é onde ele seria
atacado primeiro.

**Estado:** não sondado. É a candidata mais promissora da lista, e a única que
não foi medida hoje.

---

# Combes–Thomas sondado (21/09/2026) — fechado, negativo, e explica o resto

`scripts/measure_resolvent_decay.py`. Três testes, em ordem crescente de
fidelidade ao argumento real.

## CT.1 A premissa: a resolvente tem decaimento fora da diagonal?

Medindo `max |V_ij|` por distância de índice radial, na norma conjugada pelos
pesos RT, malha 6×18, `s=1`: o valor **oscila entre 0,4 e 2 em todas as
distâncias de 0 a 17**. Não há decaimento nenhum.

Isso não refuta a premissa — **precisa**. Os pesos RT `kappa2^n = (5/4)^n` já
extraem exatamente o decaimento que existe. Qualquer ganho de Combes–Thomas
teria de vir de um peso **além** de 5/4.

## CT.2 A versão ingênua: empurrar `kappa2`

`||B||` explode: 18,79 em `kappa2=1,25`, 32,8 em 1,4, 992 em 2,0, 4,2·10⁷ em
4,0 — crescimento compatível com `kappa2^(espalhamento radial)`.

Mas `kappa2` também **melhora** o prefator da cauda,
`P(kappa2)=2(1+r)/(1-r)` com `r=1/kappa2`, que cai de 18 para 6 entre 1,25 e
2,0. Contando os dois efeitos:

| `kappa2` | `beta` | `P` | DOFs vs 107×384 |
|---:|---:|---:|---:|
| 1,150 | 7,98 | 28,67 | 129856× |
| **1,250 (RT)** | 8,08 | 18,00 | **52415×** |
| 1,284 | 9,04 | 16,08 | 51729× |
| 1,400 | 13,89 | 12,00 | 65025× |
| 2,000 | 991,9 | 6,00 | — |

**O ótimo é `kappa2=1,284`, e vale 1,01× sobre o 5/4 de RT.** Os dois efeitos
se cancelam quase exatamente: a escolha de peso de RT, feita para a prova de
existência, já está no ótimo do problema espectral.

*Correção ao item 2.1:* lá eu otimizei só `beta` e **esqueci** que `kappa2`
também entra no prefator da cauda. Corrigido aqui, e o resultado não muda de
sinal — continua sendo ~1×.

## CT.3 O argumento de verdade: o comutador

Combes–Thomas não usa `||H_theta||`, usa `||H_theta - H||`, que para `theta`
pequeno é o comutador e é muito menor. Medido:

| `theta` | `||B_theta - B||` | `/theta` |
|---:|---:|---:|
| 1e-4 | 0,00781 | 78,11 |
| 1e-3 | 0,07835 | 78,35 |
| 1e-2 | 0,80862 | 80,86 |

Linear, com inclinação **78** — contra `L·||B|| = 66` previsto, com
`L=3,52` o espalhamento radial médio de `B`. A premissa estrutural está certa.

A condição do argumento é `||B_theta - B|| < d/2`, com `d` a distância ao
espectro. Logo **`theta_max = d/(2 L ||B||)`: a taxa é PROPORCIONAL a `d`.**

| ponto | `d = sigma_min` | `theta_max` | `N` para um fator `e` |
|---|---:|---:|---:|
| `s=1`, longe de raiz | 7,0e-3 | 5,3e-5 | 18 982 |
| **borda do contorno, `s=1/8`** | **2,6e-4** | **2,0e-6** | **505 066** |
| meio do vão, `s=0,568` | 2,1e-3 | 1,6e-5 | 63 495 |

**Na borda do contorno seriam necessários `N ~ 5·10^5` para ganhar um único
fator `e`** — pior que a cauda algébrica já entrega. Fechado, negativo.

## CT.4 Por que isto explica todas as outras falhas

A taxa de Combes–Thomas é proporcional a `d`; o `||V||` é `1/d`; a caixa
exigida escala com `||V||`. **É a mesma obstrução vestida de três maneiras.**

Na borda do contorno, `d = sigma_min = 2,6·10^-4`. Esse número sozinho explica:
o `||V||=3822`, o `theta_max=2·10^-6`, e as ~10^4 caixas. Qualquer método que
precise controlar a resolvente perto do contorno paga `1/d`, e `d` é pequeno
**por construção** — o contorno tem de passar perto do espectro que conta.

Isso dá o critério para julgar qualquer rota futura, e é o resultado mais útil
do levantamento:

> **Uma rota só tem chance se não precisar da resolvente perto do contorno.**

As quatro que foram sondadas — reponderação, ponto-base, borda, Combes–Thomas —
todas precisam, e todas ficaram entre 1,0× e 8×. Deflação analítica (item 1)
ataca `d` diretamente, aumentando a distância ao remover a raiz mais próxima, e
foi a única a render 8×; mas o piso de condicionamento do vão a limita aí.

---

# CT.5 O que a sondagem de Combes–Thomas revelou sem querer: as autofunções decaem

`scripts/measure_resolvent_decay.py` mede a resolvente. Medindo as
**autofunções** em vez dela, na malha 6×18 (`max |u_n|` por índice radial,
normalizado):

| modo | razão típica por índice (n=8..15) |
|---|---:|
| gauge, `s~0,144` | 0,63 |
| S3, `s~0,404` | 0,62 |
| física, `s~0,731` | 0,64 |

Decaimento **geométrico**, taxa `-ln(0,63) ~ 0,46` por índice radial — ou
`~0,24` na norma com pesos RT, já descontado o `(5/4)^n`.

**Compare com a taxa que Combes–Thomas consegue CERTIFICAR na borda do
contorno: `theta_max = 2·10^-6`.** São **cinco ordens de grandeza**. O
arcabouço de norma de operador está jogando fora praticamente todo o decaimento
que os objetos de fato têm.

E a razão é estrutural: a taxa de decaimento da autofunção **não depende de
`d`**. Ela vem do crescimento diagonal `mu(n+1)` vencer o acoplamento, não da
distância ao espectro. Para `n` além do suporte do fundo, a própria equação dá
`u_n = -(J_mu+s)^-1 (Bu)_n`, com coeficiente `~||B||/(mu*n) -> 0`: a recursão
passa a ser contrativa a partir de `n ~ ||B||/mu ~ 138`, **independentemente de
onde `s` esteja**.

Isso sugere uma quinta rota, que não estava na lista e que é a única compatível
com o critério da seção seguinte:

> **5. Estimativa a priori de decaimento das autofunções, em vez de bound de
> norma da resolvente.** Em vez de controlar `(H-s)^-1` no contorno, provar que
> QUALQUER autofunção com autovalor na região decai a uma taxa dada, e usar isso
> para limitar o erro de truncagem. O `1/d` não aparece.

Estado: **não sondada**, e é agora a candidata mais promissora. Ressalvas: a
medida é de uma truncagem 6×18 e os últimos índices têm efeito de borda; e
decaimento de autofunção controla o erro de truncagem de um modo **conhecido**,
enquanto contar exige excluir os **desconhecidos** — é aí que a rota teria de
provar que se sustenta.

---

# Rota 5 sondada (21/09/2026) — a primeira com sinal positivo

`scripts/measure_apriori_decay.py`, `build/spectrum/apriori-decay-12x36.json`.

## 5.1 O limiar de contração, e que ele não depende de `s`

Pela própria equação de autovalor, `u = -(J_mu+s)^-1 B u`. O fator de contração
local é `C/|D(m,n)|`, com `C` a soma de coluna ponderada de `B` e
`|D| >= max(mu(n+1), |m|/2)`. Medido na norma `eta`, em duas malhas:

| | 8×24 | 12×36 |
|---|---:|---:|
| `C` (maior soma de coluna ponderada) | 8,608 | 8,830 |
| limiar radial, `C/mu - 1` | 50 | **51** |
| limiar de Fourier, `2C` | 17 | **18** |

**Esses limiares não dependem de `s`.** Vêm do crescimento diagonal vencer o
acoplamento, não da distância ao espectro — ao contrário de `||V||`, da taxa de
Combes–Thomas e de tudo o mais medido nestas rodadas.

## 5.2 `B` é exponencialmente localizada, não bandada larga

A objeção óbvia seria que `B` acopla longe, e então a contração valeria por
banda e não por passo. Medido em 12×36:

| distância radial `L` | fração da massa ponderada |
|---:|---:|
| 1 | 0,339 |
| 3 | 0,645 |
| 8 | 0,912 |
| 12 | 0,972 |
| 20 | 0,997 |

E as entradas decaem com **razão 0,786 por índice radial**. A objeção não
procede: a localização efetiva é ~12, não os ~100 do suporte bruto de `RefA`.

## 5.3 A cauda a priori

Com `C=8,83` e contração a partir de `n>51`, o produto dos fatores dá:

| `N` | cauda relativa |
|---:|---:|
| 80 | 1,1e-3 |
| 100 | 1,6e-8 |
| 120 | 4,7e-15 |
| 250 | 1,8e-85 |
| **384 (a caixa já orçada)** | **5,4e-190** |

Na direção de Fourier, `M=107` dá 4,9e-46.

## 5.4 O ponto que muda o quadro: `||V||` deixa de doer

A rota 5 **não elimina** `||V||`. Converter "a autofunção está concentrada" em
"o autovalor finito está perto do infinito" ainda passa por um bound da
parametriz finita. A diferença é o que ele multiplica:

- **rota da resolvente:** `||V||` multiplica uma cauda **algébrica**, `~1/N`.
  Com `||V||=3822`, exige `N ~ 10^4` a `10^5`.
- **rota 5:** `||V||` multiplica uma cauda **fatorial**. Com `||V||=3822` e a
  caixa já orçada, o produto é `3822 x 5,4e-190`.

O critério que eu havia escrito estava, de novo, um pouco errado. A versão que
sobrevive a esta rodada:

> **Não é "não precisar da resolvente perto do contorno". É não deixar a
> resolvente multiplicar uma cauda algébrica.**

## 5.5 O que falta, e a ressalva que domina tudo

**O limiar de contração é `n ~ 51`, e a maior malha existente vai até `n = 35`.
Nenhuma computação deste projeto jamais entrou no regime em que a estimativa
morde.** Tudo em 5.3 é extrapolação a partir da estrutura, não medida.

Isso define o próximo passo, e ele é concreto e barato comparado a tudo que se
discutiu: **gerar uma malha com `n >= 60..80`** e verificar (a) que `C` continua
saturado em ~8,8 em vez de crescer, e (b) que o decaimento observado das
autofunções entra no regime previsto. Para calibrar o custo, `12x36` (dimensão
1422) levou 2859 s para ser gerada.

Outras ressalvas, menores:
- `C` é medido na truncagem; para o operador infinito o bound analítico na
  norma `eta` é 10,372, que moveria o limiar de 51 para 62. Não muda a ordem.
- A bola `epsilon` do fundo exato acrescenta `~3e-6` a `C`. Desprezível.
- Decaimento de autofunção controla a truncagem de um modo **conhecido**;
  contar exige excluir os **desconhecidos**. O argumento completo precisa ligar
  as duas coisas, e é aí que ele tem de ser provado, não medido.

## 5.6 A verificação está bloqueada pelos dados, não pela CPU (21/09/2026)

Tentativa de gerar malhas que cruzassem o limiar `n ~ 51`. Resultado: **12×36 é
o teto**. Toda malha acima aborta em `Field_shrink`, na asserção
`Sector_subset(sec, a->sec)`, em 5–6 s:

| malha | resultado |
|---|---|
| 7×21, 11×33, 12×36 | geram normalmente |
| 13×39, 14×42, 16×48, 18×54, 20×60, 24×72, 30×90 | **abortam** |

`RefA.dat`, o fundo de referência publicado por RT, declara `num_m 41` e
`num_n 101`, e o gerador precisa do fundo num setor **expandido** em relação ao
pedido. O teto efetivo fica em 12×36 — `n` até 35.

**Consequência dura: o limiar da rota 5 é `n ~ 51` e o dado publicado chega a
`n = 35`.** A verificação direta não é uma questão de CPU nem de dinheiro: não
há como fazê-la com o fundo disponível. Só seria possível regenerando o próprio
fundo com truncamento maior, que é a computação de RT, não a nossa.

Isso não refuta a rota 5 — mantém a §5.3 como extrapolação estrutural, e agora
com a informação de que ela **não é verificável com os dados publicados**. É uma
ressalva mais forte do que a registrada ontem.

### Duas armadilhas de método, pagas e registradas

1. **Corrida de compilação.** `build_spectral.sh` recompila o binário. Se um
   gerador anterior ainda o estiver executando, a recompilação corrompe a
   execução e o processo aborta **na mesma asserção** — o que parece limite de
   truncamento e não é. A primeira cadeia da noite morreu assim, e por isso
   uma sondagem intermediária deu "aceita" para 20×60 e 30×90.
2. **Status de saída trocado por encanamento.** `gerador | grep` devolve o
   status do `grep`: uma falha vira sucesso, e o log registrou "PRONTA em 6s
   (0 MB)". O `scripts/run_overnight_meshes.sh` final blinda as duas coisas, e
   com ele o teto de 12×36 se confirma de forma limpa.

## 5.7 Correção (21/09, tarde): a rota 5 NÃO depende da malha bloqueada

A §5.6 concluiu que a verificação estava bloqueada pelos dados. Verdade — mas
eu tratei isso como se bloqueasse a **rota**, e não bloqueia.

O `C` usado em 5.1 e 5.3 foi o **medido** (8,83 em 12×36). Mas existe bound
**analítico provado** para a mesma quantidade: `beta <= 10,382`, demonstrado em
`COMPONENT_EXTERIOR.md` e condicionado à bola publicada de RT. O limiar de
contração é consequência direta dele:

| origem de `C` | limiar radial | limiar de Fourier |
|---|---:|---:|
| medido em 12×36 (8,830) | n > 52 | \|m\| > 18 |
| **bound provado (10,382)** | **n > 61** | **\|m\| > 21** |

Cauda a priori usando **apenas o bound provado**:

| `N` | 80 | 100 | 120 | 150 | **384** |
|---|---:|---:|---:|---:|---:|
| cauda | 6,8e-2 | 3,2e-5 | 2,9e-10 | 1,9e-20 | **1,6e-165** |

E `M=107` dá 8,0e-40.

**Consequência.** A malha de `n >= 61` serviria para mostrar que o bound é
folgado — o medido é 15% menor que o provado. Folga não atrapalha uma prova:
usa-se 10,382, o limiar sobe de 52 para 61, e a cauda na caixa já orçada
continua absurdamente pequena.

**A rota 5 é atacável hoje, no papel, com o que já está provado.** O que falta
não é dado nem CPU:

1. **O lema da cauda, feito com rigor.** A recursão `u_n = -(J_mu+s)^-1 (Bu)_n`
   soma sobre TODOS os `n'`, não só a vizinhança; o produto de contrações da
   §5.3 é heurístico. A localização exponencial de `B` medida em 5.2 (razão
   0,786 por índice, 97% da massa em `|dn|<=12`) é o que torna o lema
   plausível, mas ela também precisa virar bound provado, não medida.
2. **Da concentração ao autovalor.** Converter "a autofunção está concentrada"
   em "o autovalor finito está perto do infinito", com constantes certificadas.
   É onde `||V||` reaparece — e onde ele deixa de doer, porque multiplica
   1,6e-165 em vez de 1/N.
3. **A contagem.** Decaimento controla a truncagem de um modo **conhecido**;
   contar exige excluir os **desconhecidos**. Este é o passo que pode falhar, e
   é onde o esforço deve ir primeiro, não por último.

Os itens 1 e 2 são análise comum. O item 3 é a pesquisa de verdade.

---

# 5.8 Ataque à contagem (21/09/2026): o passo que pode matar a rota

Seguindo a ordem decidida — atacar primeiro o passo que pode falhar. Três
medidas, uma positiva, uma negativa e um diagnóstico.

## 5.8.1 A seção finita converge, e muito antes do certificado

Contando autovalores dentro de `Gamma_A` (`Re s` em `(1/8,1]`, `|Im s| <= 1/4`)
em truncamentos do próprio 12×36:

| corte radial `n0` | 8 | 12 | **16** | 20 | 24 | 28 | 35 |
|---|---:|---:|---:|---:|---:|---:|---:|
| contagem | 3 | 2 | **3** | 3 | 3 | 3 | 3 |
| raiz de gauge | — | — | 0,16658 | 0,16778 | 0,16867 | 0,16837 | 0,16831 |

A contagem estabiliza em **3 a partir de `n0=16`**, e as raízes convergem a
`1e-5` por `n0=28`. Cortando o Fourier em vez do radial, estabiliza por
`m0=8`.

**A verdade converge em `n0 ~ 20`. O bound certificado exige `n > 61`, e a rota
da resolvente exige 9800× a caixa 107×384.** Todo o problema do projeto é essa
distância entre o que é verdade e o que é certificável — agora medida.

## 5.8.2 O bloco exterior se inverte bem, e o bound grosseiro é vazio

`||L_TT^-1||` na norma ponderada, com `s=1`, cortando em `n > n0`:

| `n0` | 5 | 10 | 15 | 20 | 25 | 28 | 30 |
|---|---:|---:|---:|---:|---:|---:|---:|
| `||L_TT^-1||` real | 6,94 | 3,57 | 2,45 | 1,74 | 1,29 | 1,06 | 0,90 |
| bound `1/(min|D|-C)` | vazio | vazio | vazio | vazio | vazio | vazio | vazio |

O bound uniforme é **vazio em todos os cortes** (`min|D| < C = 10,382`),
enquanto o valor real é pequeno e decrescente. A não-normalidade, que matou a
deflação como pré-condicionador, aqui **não atrapalha** — o comportamento real é
ordens de grandeza melhor que o certificado.

## 5.8.3 Por que a estimativa graduada falha: `J_mu` não é diagonalmente dominante

A perda do bound uniforme é usar o **mínimo** de `|D|` sobre um bloco em que
`|D|` cresce sem limite. A correção natural é graduar: escrever
`L = D(I + D^-1 P)` com `D` a diagonal e usar o `|D|` de cada linha.

Medido, com a decomposição correta (`P` inclui a parte **fora da diagonal de
`J_mu`**, não só `B`):

| `n0` | 5 | 15 | 24 | 30 | 32 |
|---|---:|---:|---:|---:|---:|
| `||D^-1 P||` | 9,13 | 8,53 | 7,50 | 5,07 | 3,45 |

Nunca chega perto de 1. O motivo é estrutural: a diagonal de `J_mu` é
`mu(n+1)`, mas há cerca de `n/2` termos de tamanho `2 mu n` abaixo do grau de
entrada, somando `~ mu n^2`. **`J_mu` não é diagonalmente dominante, e a
dominância nunca se instala** por mais alto que seja `n`.

*Registro de erro:* numa primeira passagem usei `P = B` apenas, esquecendo os
acoplamentos `±m/2` e `2 mu n` de `J_mu`. O resultado parecia funcionar —
`||D^-1 B||` caía abaixo de 1 em `n0=25` — mas produzia "bounds" MENORES que o
valor real, o que é impossível e denunciou o erro.

## 5.8.4 Onde isso deixa a contagem

- **Positivo:** a contagem finita é correta e estável muito antes do que
  qualquer bound exige; e o bloco exterior é bem-comportado.
- **Negativo:** nenhuma rota certificada foi encontrada. A dominância diagonal
  está estruturalmente barrada, e voltar à fatorização `V[k,f] D^-1 (I-U_k)`
  traz de volta a cauda algébrica `t(G) ~ 81/(N+1)`, que é o ponto de partida.
- **Próxima ideia, não sondada:** graduar **nas coordenadas da fatorização**,
  não nas coordenadas originais. A fatorização de RT resolve exatamente a
  estrutura que quebra a dominância diagonal; uma estimativa graduada ali pode
  funcionar onde a ingênua falha. É o caminho natural, e é análise de papel.

A conclusão de método continua valendo, e hoje ficou quantificada: **o obstáculo
não é o objeto, é o certificado.** O operador se comporta bem; o que não se
consegue é provar que ele se comporta bem sem pagar 10³ em tamanho de caixa.

## 5.8.5 Mais duas rotas de certificado, ambas negativas

**Determinante de Fredholm.** Como `Q0` é compacta, a contagem é o winding de
`det(I + Q0 B)`, e o erro de truncagem seria governado pela norma **traço**, que
vê o decaimento dos valores singulares. Medido em 12×36, `s=1`:

| `k` | 1 | 50 | 200 | 400 | 800 |
|---|---:|---:|---:|---:|---:|
| `sigma_k` | 2,13 | 1,12 | 0,589 | 0,283 | 8,4e-3 |

O decaimento é **algébrico**, `sigma_k ~ k^-0,66` na faixa útil (o colapso
depois de `k~800` é artefato da matriz finita). Norma traço medida: 321,7, e
`sum k^-0,66` **diverge**. **`Q0 B` não é de classe traço**, e o determinante de
Fredholm usual não está definido.

**Determinante regularizado `det_2`.** Para Hilbert–Schmidt — que `Q0 B` é,
porque `sum k^-1,32` converge — existe `det_2(I+K) = det((I+K)e^{-K})`, que
anula exatamente onde `I+K` é singular. O erro de truncagem passa a ser
governado pela norma HS. Medido:

| corte `n0` | 4 | 12 | 20 | 28 | 32 |
|---|---:|---:|---:|---:|---:|
| `||K-K_n0||_HS`, `s=1` | 14,71 | 11,70 | 8,79 | 5,65 | 3,64 |
| `s=1/8` | 23,78 | 17,81 | 12,86 | 8,02 | 5,07 |

Cai devagar, e **estruturalmente é pior que a norma de operador**: as entradas
de `K` na casca `n` escalam como `1/(mu n)` e o número de pares por casca é
constante, logo `||K-K_N||_HS ~ N^-1/2`, contra o `1/N` da cauda de operador.
Norma HS é maior **e** decai mais devagar. `det_2` é um retrocesso, não um
avanço.

**E a grade na fatorização já estava sendo usada.** A ideia anunciada em 5.8.4
era graduar nas coordenadas de `V[k,f] D^-1 (I-U_k)` em vez das originais. Ao
escrever a estimativa vê-se que `||Q0` fora da casca `N|| <= 81/(N+1)` **é**
exatamente essa grade — ela usa o índice da casca, não um mínimo global. Dá
`N > 840` para o bloco exterior, contra `n0 ~ 20` da verdade. Não havia grade
inexplorada.

## 5.8.6 Saldo do ataque à contagem

Quatro tentativas de certificado, quatro negativas, cada uma por um motivo
preciso:

| tentativa | resultado | motivo |
|---|---|---|
| bound uniforme `1/(min|D|-C)` | vazio | usa o mínimo de `|D|` num bloco onde ele cresce |
| grade diagonal `D(I+D^-1P)` | `||D^-1P||` fica em 8–9 | `J_mu` não é diagonalmente dominante |
| determinante de Fredholm | indefinido | `Q0B` não é de classe traço, `sigma_k ~ k^-0,66` |
| `det_2` / Hilbert–Schmidt | pior | cauda HS `~N^-1/2` contra `1/N` do operador |
| grade na fatorização | já em uso | é o `81/(N+1)`, dá `N>840` |

E o fato empírico continua de pé, sem explicação certificável: **a contagem
finita já está correta em `n0 ~ 20`**, e nenhum argumento disponível prova isso
com menos de `N ~ 840` no bloco exterior — e muito mais no critério acoplado.

O obstáculo tem agora nome e tamanho. Não é o operador: é que toda ferramenta
de certificação disponível mede o objeto por uma norma que perde a estrutura
que o faz funcionar.

---

# 5.9 Revisão de completude do mapa (21/09/2026)

O mapa de 5.8.6 tinha quatro becos. A pergunta certa é se são só quatro.
Varrendo `docs/` por famílias de técnica padrão para este tipo de problema:

| técnica | ocorrências em todo o projeto |
|---|---:|
| campo de valores (numerical range) | **0** |
| realização em `l^2` / espaço de Hilbert | **0** — tudo é `l^1` com pesos, herdado de RT |
| operadores limite / teoria da seção finita | **0** |
| pseudoespectro | **0** — apesar de a não-normalidade ser tema central |
| renormalização / escala auto-similar como ferramenta | **0** |
| formas quadráticas / setorial | **0** |
| tipo Gershgorin | **0** |

**O mapa não estava completo.** Nenhuma dessas famílias tinha sido considerada.

## 5.9.1 Campo de valores: falha no contorno, por convexidade

Para `T` num espaço de Hilbert, `||(T-s)^-1|| <= 1/dist(s, W(T))`, um bound
**certificado** que não passa por série de Neumann. Medido em 12×36, na norma
`l^2` ponderada:

| ponto `s` | dentro de `W`? | `||V||` real |
|---|---|---:|
| 1,000 | dentro, margem 1,63 | 214,44 |
| 0,125 (a borda) | dentro, margem 2,50 | 1846,19 |
| 0,568 | dentro, margem 2,06 | 692,78 |

**Inútil no contorno, e o motivo é estrutural.** `W` é **convexo** e contém todo
o espectro. Como `J_mu` tem autovalores `mu(n+1) -> infinito`, o espectro de `L`
se estende a `Re s -> -infinito`; um convexo que contém pontos muito à esquerda
**e** as raízes em 0,168–0,733 cobre necessariamente a borda em `Re s = 1/8`.
Nenhuma escolha de contorno que separe raízes escapa disso.

*Registro de erro:* na primeira leitura tomei a distância positiva como "fora de
`W`". É o contrário — o mínimo sobre os semiplanos de suporte é positivo
justamente quando o ponto está **dentro**.

## 5.9.2 Mas o subproduto é grande: F3 de 414 para ~12

`spec(T) subset W(T)`, então a borda **direita** de `W` dá um `R` certificado
com "não há espectro para `Re s > R`". Medido:

| malha | `R = max Re W(-A)` | só `J_mu` | `||B||_2` |
|---|---:|---:|---:|
| 4×12 | 1,9998 | **0,0888** | 2,766 |
| 6×18 | 2,3219 | **0,0888** | 3,282 |
| 8×24 | 2,4875 | **0,0888** | 3,605 |
| 12×36 | **2,6290** | **0,0888** | 3,976 |

`R_J = 0,0888` é **constante nas quatro malhas** — o operador livre quase não se
estende à direita, e isso tem cara de propriedade exata, não de coincidência
numérica. E `R <= R_J + ||B||_2`.

Para certificar `||B||_2` sem medir, o teste de Schur:
`||B||_2 <= sqrt(||B||_col · ||B||_lin)`. Em 12×36 dá 12,20 contra 3,98 reais —
folgado 3×, mas é um bound. Com o `beta <= 10,382` já provado na coluna,

> **`R <= 0,0888 + sqrt(10,382 · ||B||_linha) ~ 12 a 15`, contra os 414 de hoje.**

`F3` está no ledger como **aberta**, com evidência exigida "estimativa de
energia quantitativa" — e campo de valores **é** uma estimativa de energia, `
<Lu,u>`. Isso é um caminho concreto para fechá-la, e encolhe a região descoberta
de `1 < Re s < 414` para `1 < Re s < ~12`.

Ressalvas, todas sérias: (a) isto é em `l^2` ponderado e o projeto vive em `l^1`
ponderado — o espectro depende do espaço, e relacionar as duas realizações é
trabalho a fazer; (b) a norma de **linha** de `B` não tem bound provado, só a de
coluna; (c) `R_J = 0,0888` é medido em quatro malhas, não derivado.

## 5.9.3 O que continua sem ter sido tentado

- **Operadores limite / teoria da seção finita** (Rabinovich–Roch–Silbermann).
  É literalmente a teoria da pergunta "a seção finita converge, e quão rápido".
  Para operador com diagonal crescente ela deve responder que sim; a versão
  **quantitativa com constantes explícitas** é o que falta. Zero ocorrências no
  projeto.
- **A realização `l^2` como um todo.** Todo o arcabouço está em `l^1` porque os
  bounds de RT estão lá. O problema espectral pode ser muito melhor comportado
  em `l^2` — 5.9.2 é a primeira evidência disso.
- **Estrutura auto-similar / renormalização.** O fundo é DSS; se houver
  covariância de escala entre cascas, a cauda poderia sair por ponto fixo em vez
  de série. Especulativo, e as coordenadas de RT já absorvem parte da
  auto-similaridade.

**Pseudoespectro**, embora ausente, não é novidade: a condição "`Gamma` fora do
`epsilon`-pseudoespectro de `L_N`" é `||Delta L|| · ||(L_N-s)^-1|| < 1`, a mesma
já em uso. É reembalagem, não ferramenta nova.

## 5.10 Aprofundando o F3 pelo campo de valores (21/09/2026)

`scripts/numerical_range_f3.py`, `build/spectrum/numerical-range-f3.json`.
Três buracos a fechar; um fecha sozinho, um fica derivável, um continua aberto.

### 5.10.1 A ponte `l^1` → `l^2` é gratuita, e na direção que importa

A ressalva mais séria era que o projeto vive em `l^1` ponderado e o campo de
valores exige Hilbert. Ela se resolve numa linha: para pesos positivos,

```
||x||_{2,w}^2 = sum (w_j |x_j|)^2  <=  (sum w_j |x_j|)^2 = ||x||_{1,w}^2,
```

logo `l^1(w) subset l^2(w)` com norma 1. **Toda autofunção em `l^1` é
autofunção em `l^2`**, e o domínio acompanha: `L u = 0` em `l^1(w)` dá
`J_mu u = -(B+s)u` em `l^1(w) subset l^2(w)`.

Portanto **excluir espectro em `l^2` exclui em `l^1`** — exatamente a direção
que F3 pede. A recíproca não vale, e não é necessária.

### 5.10.2 `R_J` é uma constante 1-D, derivável

O autovetor de `lambda_min(Herm(J_mu))` vive **inteiramente na componente 3**,
com perfil radial fixo (`n=0`: 0,86; `n=1`: 1,00; `n=2`: 0,44), e o índice `m`
apenas acompanha a malha. Reconstruindo o problema 1-D correspondente —
diagonal `mu(n+1)`, acoplamento `2 mu n` para **todo** `n' < n`, pesos
`sqrt((2-[n=0]) kappa2^n)` — obtém-se

| `N` do modelo 1-D | 8 | 12 | 20 | 40 | 400 |
|---|---:|---:|---:|---:|---:|
| `R_J` | 0,0846210 | 0,0887652 | **0,0887710** | 0,0887710 | 0,0887710 |

**Reproduz exatamente o valor da matriz completa (0,088770993) e estabiliza em
`N = 20`.** É portanto certificável: aritmética exata ou intervalar num bloco
20×20, mais um argumento de cauda para o resto — onde a diagonal cresce e os
acoplamentos decaem com `kappa2^{-(n-n')/2}`.

### 5.10.3 O buraco que fica: `||B||_2` não convergiu

| malha | `R = max Re W(-A)` | `||B||_2` | `rho(|B|)` |
|---|---:|---:|---:|
| 4×12 | 1,9998 | 2,766 | 2,963 |
| 6×18 | 2,3219 | 3,282 | 3,635 |
| 8×24 | 2,4875 | 3,605 | 4,045 |
| 12×36 | **2,6290** | **3,976** | 4,493 |

`rho(|B|)` também cresce, então **não é artefato de peso** — é do operador. A
norma de **linha** de `B` cresce e seu máximo migra para a borda da malha, o
que é esperado: os pesos de RT foram escolhidos para limitar a **coluna**, que é
o que `l^1` pede, e não a linha.

Ajustando `X = X_inf - c·dim^-p`:

| quantidade | limite | `p` | resíduo |
|---|---:|---:|---:|
| `R` medido | **2,914** | 0,503 | 0,006 |
| `||B||_2` | 6,138 | 0,191 | 0,005 |

O `R` medido converge bem (`p ~ 1/2`, resíduo 0,006). O bound solto
`R <= R_J + ||B||_2` daria `~6,2`; o valor direto extrapola para **~2,9**.

### 5.10.4 Saldo

> **`R ~ 3` (medido, extrapolado) ou `~6` (por bound solto), contra os 414
> certificados hoje. Uma melhora de 70× a 140× na fronteira de exclusão de F3.**

A região descoberta encolhe de `1 < Re s < 414` para `1 < Re s < ~6`.

O que falta para virar certificado: (a) `R_J` em aritmética exata no bloco 1-D
de dimensão 20 mais o argumento de cauda — trabalho delimitado; (b) um bound
**provado** para `||B||_2` no infinito, que hoje não existe nem na medida (ela
não convergiu); (c) escrever a estimativa de energia diretamente, sem passar
pela matriz finita — que é o que a linha F3 do ledger pede.

E a ressalva de método: tudo em 5.10.3 é medida em malhas até 12×36, o teto dos
dados publicados. O `p = 0,503` do ajuste de `R` é tranquilizador; o
`p = 0,191` do `||B||_2` não é.

## 5.11 `R_J` caracterizado por completo (21/09/2026)

### 5.11.1 Homogeneidade e forma

A matriz do modelo 1-D é `mu·[diag(n+1) + M]`, com `M[n',n] = n` para `n' < n`,
conjugada pelos pesos. Logo, **exatamente**,

```
R_J = mu · f(kappa2),     f independe de mu.
```

Em `kappa2 = 5/4`: `f = 0,527434696`, e `mu · f = 0,088770993` — o valor medido
na matriz completa, a nove dígitos.

### 5.11.2 `f` é monótona decrescente, e cruza zero

| `kappa2` | 1,00 | 1,10 | 1,20 | **1,25 (RT)** | 1,30 | 1,50 | 2 | 8 | 4096 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `f` | ~2·10⁴ | 7,184 | 1,267 | **0,5274** | 0,1171 | −0,5025 | −0,7655 | −0,9409 | −0,9999 |

- **`f -> -1` quando `kappa2 -> infinito`**, e o motivo é claro: os acoplamentos
  `2 mu n` são multiplicados por `kappa2^{-(n-n')/2}` e somem, sobrando
  `Herm(J_mu) -> diag(mu(n+1))`, cujo mínimo é `mu`. Então `R_J -> -mu`.
- **`f` cruza zero em `kappa2 = 1,3204632636`.** Abaixo disso `J_mu` não é
  acretivo; acima, é — e deixa de contribuir para `R`.
- `f` diverge quando `kappa2 -> 1`: sem peso o operador livre é catastrófico.

### 5.11.3 Dimensão efetiva: 18

Em `kappa2 = 5/4`, `f` converge a nove dígitos já em **`N = 18`**:

| `N` | 10 | 14 | 16 | **18** | 20 | 25 |
|---|---:|---:|---:|---:|---:|---:|
| `f` | 0,526136 | 0,5274341 | 0,52743469 | **0,527434696** | 0,527434696 | 0,527434696 |

**`R_J` é certificável por um cálculo exato 18×18 mais um argumento de cauda.**
É a peça mais fácil de toda a agenda analítica do projeto.

### 5.11.4 Mas `R_J` é o termo pequeno — e isso encerra a linha

A liberdade de escolher o peso da realização `l^2` é real: a inclusão
`l^1(w_RT) subset l^2(w')` pede `w' <= C·w_RT`, o que sozinho só permitiria
**diminuir** `kappa2` — e isso piora `R_J`. Mas o decaimento a priori da rota 5
põe toda autofunção em **qualquer** espaço com peso, liberando `kappa2` para
cima. As duas rotas se encaixam.

Otimizando, em 8×24:

| `kappa2` | `R_J` | `||B||_2` | `R` total |
|---|---:|---:|---:|
| 1,25 (RT) | 0,0888 | 3,605 | 2,4875 |
| 1,40 | −0,0509 | 2,977 | 2,1515 |
| **1,60** | **−0,1030** | 2,895 | **2,0397** |
| 2,00 | −0,1288 | 4,274 | 2,4067 |

**Ganho total: 1,22×.** `R_J` vai de +0,089 a −0,103 — uma variação de 0,19 —
enquanto `||B||_2` fica em torno de 3. **Quem domina `R` é `B`, não o operador
livre.**

### 5.11.5 Saldo

`R_J` está encerrado: forma fechada em `mu·f(kappa2)`, monotonia e limites
conhecidos, ponto crítico localizado, dimensão efetiva 18, certificação barata.
E, justamente por estar encerrado, ficou claro que **ele não é onde está a
dificuldade**. Para o `R` de F3, tudo depende de um bound provado para
`||B||_2` no infinito — que hoje não existe, e cuja medida não convergiu
(§5.10.3).

## 5.12 `||B||_2`: um artefato retratado e uma rota provável (21/09/2026)

### 5.12.1 Retratação: o "peso balanceado" era degenerescência

Tentei balancear coluna e linha pelo peso de Perron `D = diag(sqrt(v/u))`, com
`u, v` os vetores de Perron de `|B|`. O resultado parecia ótimo — `||B||_2`
caía de 3,98 para 2,54 e passava a convergir bem (limite 2,875 com `p=0,453`
contra 6,138 com `p=0,191`).

**É artefato.** `|B|` é **redutível**: o vetor de Perron à direita tem 612
entradas nulas de 1422 em 12×36 (43%), e o peso resultante tem condicionamento
`10^151`. O "balanceamento" estava simplesmente zerando blocos, não
equilibrando normas. O ganho não existe.

| malha | entradas nulas de `v` | condicionamento do peso |
|---|---:|---:|
| 4×12 | 60 de 138 | 3,8·10¹⁵⁰ |
| 6×18 | 144 de 333 | 5,0·10¹⁵⁰ |
| 8×24 | 264 de 612 | 6,5·10¹⁵⁰ |
| 12×36 | 612 de 1422 | 1,2·10¹⁵¹ |

E a família produto `kappa1^m kappa2^n` **não consegue balancear**: o melhor dá
coluna 7,87 contra linha 13,77, ainda 1,75× de desequilíbrio, e o `R` resultante
(2,68) é **pior** que o de RT (2,49).

### 5.12.2 A rota provável: Young, não balanceamento

O caminho certo não é escolher peso, é usar a **estrutura**. `B` é
essencialmente multiplicação pelo fundo, e para convolução com peso
submultiplicativo (`w_{i+j} <= w_i w_j`, que `kappa^n` satisfaz) vale a
desigualdade de Young:

```
||f * g||_{l^2(w)}  <=  ||f||_{l^1(w)} · ||g||_{l^2(w)}
```

ou seja, **`||B||_{l^2(w)} <= ||B||_{coluna, l^1(w)}`** — a mesma quantidade que
já tem bound provado. E é a mesma submultiplicatividade da norma
Fourier–Chebyshev que `COMPONENT_EXTERIOR.md` já usa para limitar a convolução
pela norma do coeficiente.

Consistente em todas as malhas, com folga estável em torno de 0,42:

| malha | `||B||_2` | coluna `l^1` | razão |
|---|---:|---:|---:|
| 4×12 | 2,766 | 6,591 | 0,420 |
| 6×18 | 3,282 | 8,078 | 0,406 |
| 8×24 | 3,605 | 8,608 | 0,419 |
| 12×36 | 3,976 | 8,830 | 0,450 |

### 5.12.3 O que isso daria

Com `beta <= 10,382` **já provado** em `COMPONENT_EXTERIOR.md`:

> **`R <= R_J + ||B||_2 <= 0,0888 + 10,382 = 10,47`**, contra os **414**
> certificados hoje — **fator 40** na fronteira de exclusão de F3.

E a região descoberta de F3 encolhe de `1 < Re s < 414` para `1 < Re s < 10,5`.

### 5.12.4 O que falta escrever

1. **A desigualdade de Young na realização exata.** `B` não é uma convolução
   pura: tem a parte linear `S Gamma1` (que baixa grau e preserva `m`), as
   projeções de paridade e o `R_x`. Cada peça precisa do seu argumento; a
   submultiplicatividade já registrada cobre a parte bilinear do fundo.
2. **`R_J` em aritmética exata**, bloco 1-D 18×18 mais cauda (§5.11.3).
3. **A estimativa de energia escrita diretamente**, sem passar pela matriz
   finita — que é a evidência que a linha F3 do ledger pede.

Nenhum dos três é pesquisa; os três são redação com constantes. É a agenda mais
concreta que apareceu em três dias.

## 5.13 `R_J` certificado em aritmética exata (21/09/2026)

`scripts/certify_free_numerical_range.py`,
`build/spectrum/free-numerical-range.json`. Primeiro item da agenda de 5.12.4,
feito.

### 5.13.1 O truque que torna o problema racional

Conjugar por `sqrt(w)` introduz irracionais (`sqrt(4/5)`), o que estraga a
aritmética exata. Evita-se reformulando como **pencil generalizado**:

```
lambda_min( Herm(D J D^-1) ),  D = diag(sqrt(w))
   =  menor lambda de   Herm(W J) y = lambda W y,   W = diag(w),
```

e **`W J` e `W` são ambos racionais**. A prova de `lambda_min >= lambda0` é
então exibir um `lambda0` racional e verificar

```
Herm(W J) - lambda0 · W  definida positiva,
```

por `LDL^T` exato com todos os pivôs positivos. Nada de intervalo, nada de
ponto flutuante.

### 5.13.2 O certificado do bloco finito

Com `lambda0 = -1777/20000 = -0,08885`:

| `N` | definida positiva? | menor pivô |
|---:|---|---:|
| 18 | **sim** | 0,2571571 |
| 20 | **sim** | 0,2571571 |
| 30 | **sim** | 0,2571571 |
| 110 | **sim** | 0,2571571 |

> **`R_J < 0,08885`, certificado em aritmética racional.**

O valor medido é 0,088770993, então o bound está 0,09% acima — apertado. E o
menor pivô **não muda** com `N`: ele vive no bloco inicial, e acrescentar
cascas só acrescenta pivôs maiores.

### 5.13.3 A cauda: os pivôs crescem, com assintótica limpa

Dominância diagonal **não** vale — a coluna `n` tem `n` entradas de tamanho
`2 mu n` contra diagonal `mu(n+1)`. O que vale é a positividade dos pivôs do
Schur. Calculados exatamente até `N = 110`, todos positivos, e normalizados
pela diagonal crua `w_n mu (n+1)`:

| `n` | 20 | 40 | 60 | 80 | 100 |
|---|---:|---:|---:|---:|---:|
| `pivo / (w_n mu (n+1))` | 0,26954 | 0,23874 | 0,22683 | 0,22053 | 0,21662 |

Ajuste `r(n) = a + b/n`: **`a = 0,2043`, `b = 1,331`**, resíduo máximo
`1,3e-3`. A razão converge a **0,204 > 0**, não a zero.

> **A cauda é benigna: `pivo_n >= 0,20 · w_n · mu (n+1)` para `n >= 20`.**

Isso dá a forma do lema que falta escrever: uma indução mostrando que a razão
do complemento de Schur permanece acima de `1/5`. É análise elementar sobre uma
recursão explícita, não pesquisa.

### 5.13.4 Estado da agenda de F3

| item | estado |
|---|---|
| `R_J` no bloco finito | **feito**, exato, `< 0,08885` |
| cauda de `R_J` | forma do lema identificada; falta escrever a indução |
| Young para `||B||_2` | falta (§5.12.2) |
| estimativa de energia direta | falta |

Com `R_J < 0,08885` certificado e `||B||_2 <= 10,382` pela rota de Young,
`R <= 10,47` contra os 414 de hoje.

## 5.14 A cauda de `R_J`, fechada em forma fechada (21/09/2026)

### 5.14.1 A matriz é semisseparável mais diagonal

Verificado exatamente:

```
Herm(W J_mu)_ij = u_min(i,j) · v_max(i,j)  (i != j),
Herm(W J_mu)_ii = u_i v_i + w_i mu,
com  u_i = w_i = (2-[i=0]) kappa2^i   e   v_j = mu j.
```

Matrizes dessa forma — semisseparável de posto 1 mais diagonal — são as
inversas de tridiagonais, e seu `LDL^T` colapsa numa recursão **escalar**.

### 5.14.2 A recursão escalar

Escrevendo o complemento de Schur no passo `k` como `u^{(k)}_i = u_i - alpha_k v_i`:

```
p_k       = (w_k - alpha_k mu k)·mu k + w_k(mu - lambda0)
alpha_{k+1} = alpha_k + (w_k - alpha_k mu k)^2 / p_k
```

Duas linhas, um escalar. **Reproduz o `LDL^T` completo exatamente** — os 30
primeiros pivôs batem como racionais, dígito a dígito.

### 5.14.3 O ponto fixo, e o `1/5`

Com `gamma_k = alpha_k mu k / w_k`, a recursão vira, para `k` grande e
`lambda0 -> 0`,

```
gamma·kappa2 = gamma + (1 - gamma)   =>   gamma = 1/kappa2.
```

E como `p_k = w_k mu[k(1-gamma_k) + 1 - lambda0/mu]`, o pivô normalizado tende a

> **`p_k / (w_k mu (k+1))  ->  1 - 1/kappa2 = 1/5`.**

Confirmado numericamente até `k = 200`, em aritmética exata:

| `k` | 20 | 50 | 100 | 150 | 199 | limite |
|---|---:|---:|---:|---:|---:|---:|
| `gamma_k` | 0,79338 | 0,79422 | 0,79649 | 0,79751 | 0,79806 | **4/5** |
| `r_k` | 0,26954 | 0,23170 | 0,21662 | 0,21127 | 0,20857 | **1/5** |

O ajuste numérico da §5.13.3 dava `a = 0,2043`; era **`1/5`**, com a correção
`1/n` que o ajuste absorveu em `b = 1,331`.

**A cauda está fechada:** todos os pivôs são positivos e crescem como
`(1/5)·w_k·mu(k+1)`, com constante explícita e dependente só de `kappa2`. O que
falta é redigir a atratividade do ponto fixo — a recursão é escalar e
monotônica na região relevante, o que é exercício, não pesquisa.

### 5.14.4 Estado de `R_J`

| peça | estado |
|---|---|
| bloco finito | **certificado exato**, `R_J < 0,08885` |
| cauda | **forma fechada**, `r_k -> 1 - 1/kappa2`, pivôs positivos verificados a `k=200` |
| atratividade do ponto fixo | exercício de redação |

`R_J` deixou de ser uma constante medida e virou uma constante **derivada**,
com a dependência em `kappa2` explícita. É a primeira peça do projeto que
atravessou inteira de "medido" para "provado".

## 5.15 `R_J` completo, incluindo a cauda (21/09/2026)

A §5.14 deu a forma fechada; falta a atratividade do ponto fixo. Ela fecha em
três passos, todos verificáveis.

### 5.15.1 A recursão, em forma fechada

```
gamma_{k+1} = ((k+1)/kappa2) · [ gamma_k/k + (1-gamma_k)^2 / delta_k ],
delta_k     = k(1-gamma_k) + c,        c = 1 - lambda0/mu = 1,527904...
p_k         = w_k · mu · delta_k.
```

Verificada **exatamente** contra a recursão original para `k = 1..58`.

O invariante certo **não** é `gamma_k < 1` — de fato `gamma_5 = 1,1018` — e sim
`delta_k > 0`, que é o que torna o pivô positivo.

### 5.15.2 Os três passos

**(a) Checagem finita, `k <= 15`.** Em aritmética racional exata, `delta_k > 0`
em todos, com mínimo `0,74983` em `k = 8`.

**(b) Intervalo invariante, `k >= 15`.** O mapa `T` leva `[0,79 ; 0,80]` em si
mesmo, e é **monótono em `gamma`** no intervalo — então basta verificar os
extremos:

| `k` | 15 | 20 | 30 | 100 | 1000 | 10⁵ |
|---|---:|---:|---:|---:|---:|---:|
| `T(0,79)` | 0,794803 | 0,792946 | 0,792782 | 0,796492 | 0,799585 | 0,7999958 |
| `T(0,80)` | 0,795743 | 0,793565 | 0,793110 | 0,796531 | 0,799586 | 0,7999958 |

Ambos dentro de `[0,79 ; 0,80]` em todos os `k` testados, e convergindo a
`4/5 = 1/kappa2`.

**(c) Dentro do invariante**, `delta_k = k(1-gamma_k) + c >= 0,20 k + 1,528 > 0`.

> **Logo todos os pivôs são positivos, `Herm(WJ) - lambda0 W` é definida
> positiva no espaço infinito, e `R_J < 0,08885` está certificado incluindo a
> cauda.**

### 5.15.3 Estado

| peça de `R_J` | estado |
|---|---|
| redução ao problema 1-D | **feita** (§5.11) |
| bloco finito, exato | **feito** (§5.13) |
| forma fechada da cauda | **feita** (§5.14) |
| atratividade do ponto fixo | **feita** (§5.15) |

`R_J` está **completo**. Falta transcrever para o Lean, que é o trabalho de
sempre, não descoberta.

## 5.16 `||B||_2`: a estrutura do argumento de Young

A desigualdade `||B||_{l^2(w)} <= ||B||_{coluna, l^1(w)}` **não vale para
matriz qualquer** — vale pela estrutura. Duas peças:

**Parte bilinear (multiplicação pelo fundo).** No índice de Fourier é
convolução, e Young com peso submultiplicativo dá o resultado direto. No índice
de Chebyshev, `T_m T_n = (T_{m+n} + T_{|m-n|})/2`, e o espaço de coeficientes
com peso `kappa2^n` é o espaço de funções analíticas numa elipse; multiplicação
por `f` tem norma `<= sup_elipse |f| <= sum |f_n| kappa2^n`, que é exatamente a
norma `l^1` ponderada. **É a mesma submultiplicatividade que
`COMPONENT_EXTERIOR.md` já usa**, agora lida em `l^2`.

**Parte linear (`S Gamma1`, isto é `R_x` e reflexões).** `R_x f = (f - f(0))/xi`
baixa o grau e preserva `m`; precisa do seu próprio bound, e
`COMPONENT_EXTERIOR.md` já registra `||R_x|| <= 40/9` na norma `l^1`.

Numericamente a desigualdade vale com folga em **todas** as combinações de peso
testadas (8×24):

| `kappa2` | 1,15 | 1,25 | 1,40 | 1,60 |
|---|---:|---:|---:|---:|
| `||B||_2 / ||B||_col` | 0,52 | 0,43 | 0,19 | 0,05 |

Nunca chega perto de 1 — a folga cresce quando o peso cresce, como a leitura
via elipse prevê (o `sup` na elipse cresce muito mais devagar que a soma `l^1`).

**O que falta:** escrever as duas peças com as constantes. É redação sobre
estimativas que o projeto já tem, não técnica nova.

## 5.17 O fundo tem 6,7× de folga entre a norma `l^1` e o sup na elipse

A leitura via elipse (§5.16) não é só uma justificativa para
`||B||_2 <= ||B||_{l^1}` — ela é **muito mais apertada**, e a folga é
mensurável direto no fundo publicado.

Lendo os quatro campos de `RefA.dat` (41×101 cada, 4141 coeficientes por campo)
e comparando, para os mesmos pesos, a soma `l^1` ponderada com o supremo no
domínio complexificado (anel `|z|=kappa1` em Fourier, elipse de Bernstein de
parâmetro `kappa2` em Chebyshev):

| `kappa1` | `kappa2` | `sup / l^1` |
|---:|---:|---:|
| 65/64 | 1,15 | 0,1741 |
| **65/64** | **1,25 (RT)** | **0,1491** |
| 65/64 | 1,40 | 0,1307 |
| 1,00 | 1,25 | 0,1494 |
| 1,05 | 1,25 | 0,1484 |

**Estável em 0,13–0,17, e 0,149 no peso de RT.** A soma `l^1` superestima o
objeto que de fato limita o operador em `l^2` por um fator **6,7**.

### 5.17.1 O que isso faria com `beta`

`ALGEBRAIC_TAIL.md` decompõe o bound de `B` na norma RT:

```
||B|| <= (17/100)(80/9)  +  2·10,5  +  6(2^-25 + 2^-277)  <  23
          `- S Gamma1        `- bilinear no fundo     `- correcao
           = 1,511             = 21
```

O termo de 21 é a contribuição bilinear do fundo — exatamente a quantidade cuja
leitura via elipse vale 0,149 do valor `l^1`. Recuperando essa folga:

```
||B||_2  <~  1,511 + 21·0,149 + eps  ~  4,6   (contra 23)
```

E `4,6` é notavelmente próximo do `||B||_2` **medido**, que extrapola para
`~6,1` com os pesos de RT e vale 3,98 em 12×36. O bound deixaria de ser
folgado por 2,6× e passaria a ser quase justo.

### 5.17.2 O que isso faria com F3

| rota | `R` |
|---|---:|
| Neumann grosseiro (hoje, no ledger) | **414** |
| campo de valores + `beta<=10,382` (`l^1`) | ~10,5 |
| campo de valores + leitura via elipse | **~4,7** |
| medido, extrapolado | ~2,9 |

> **`R ~ 4,7` contra 414 é um fator 88**, e põe o certificado a 1,6× do valor
> medido.

### 5.17.3 Ressalva

O fator 0,149 é do **fundo**, medido em `RefA.dat`. Transportá-lo para o bound
de `B` exige que a estimativa bilinear de RT seja reescrita com o supremo no
lugar da soma `l^1` — e é aí que ela pode perder, porque o bound de RT usa a
estrutura de `Gamma2`, não só a norma do fundo. **A folga está medida; a
transferência ainda não está escrita.**

## 5.18 A transferência para a elipse: vale, mas com constante (21/09/2026)

A §5.17 mediu 6,7× de folga entre a soma `l^1` do fundo e seu supremo no
domínio complexificado, e sugeriu trocar uma pela outra no bound de `B`. Fui
testar a transferência antes de usá-la, e ela **não é exata**.

### 5.18.1 O contraexemplo

O espaço de **coeficientes de Chebyshev** com peso `rho^n` **não é um espaço de
Hardy**: a base `T_n` não é ortogonal no produto interno natural da elipse.
Logo `||M_f||_{l^2} <= sup_elipse |f|` não vale automaticamente. Testado num
modelo controlado, com `rho = 5/4`:

| `f` | `||M_f||_2` | `sup` na elipse | soma `l^1` | vale? |
|---|---:|---:|---:|---|
| `T_0 = 1` (identidade) | 1,0000 | 1,0000 | 1,0000 | sim |
| `T_1` | 1,2315 | 1,0250 | 2,5000 | **não, +20%** |
| `T_5` | 1,7779 | 1,6897 | 6,1035 | não, +5% |
| `soma 0,7^k T_k` | 3,3756 | 3,9276 | 10,1895 | sim |

### 5.18.2 A constante

Varrendo 400 coeficientes aleatórios com decaimento variado, mais os `T_j`
isolados:

| `rho` | pior `||M_f||_2 / sup_elipse` |
|---:|---:|
| 1,15 | 1,270 |
| **1,25 (RT)** | **1,210** |
| 1,40 | 1,211 |

> **A transferência vale com `C ~ 1,2`, não com `C = 1`.**

O pior caso é sempre um `T_j` isolado de grau baixo — coeficiente concentrado,
que é o oposto do fundo, cujos coeficientes decaem. Para funções com
decaimento, como `RefA`, a razão fica **abaixo de 1**.

### 5.18.3 As contas refeitas, agora consistentes

| quantidade, norma RT | valor |
|---|---:|
| bound `l^1` de `ALGEBRAIC_TAIL` | **23** |
| elipse com `C = 1` (o que eu ia afirmar) | 4,64 |
| **elipse com `C = 1,2`** | **5,6** |
| `||B||_2` medido, extrapolado | 4,90 |

Com `C = 1` a "previsão" ficava **abaixo** do medido — um bound não pode fazer
isso, e foi esse desacordo que denunciou o problema. Com `C = 1,2` a ordem fica
certa: `5,6 > 4,90`, e o bound está a 14% do valor medido.

### 5.18.4 Saldo para F3

| rota | `R` |
|---|---:|
| Neumann grosseiro (ledger hoje) | **414** |
| campo de valores + `beta <= 10,382` | ~10,5 |
| **campo de valores + elipse com constante** | **~5,7** |
| medido, extrapolado | ~3,0 |

> **`R ~ 5,7` contra 414: fator 73**, e o certificado a menos de 2× do medido.

**O que falta escrever:** a constante `C` da transferência precisa de prova, não
de varredura numérica — é uma desigualdade entre a norma de operador de
multiplicação em `l^2(rho^n)` de coeficientes de Chebyshev e o supremo na elipse
de Bernstein. É um resultado de análise clássica, e a varredura diz que a
constante certa está perto de `6/5`.

## 5.19 A estrutura da constante: Toeplitz mais resto (22/09/2026)

A §5.18 mediu a constante `C ~ 1,2` por varredura. Aqui ela ganha estrutura,
que é o que uma prova precisa.

### 5.19.1 A decomposição exata

Da fórmula do produto `T_j T_n = (T_{j+n} + T_{|j-n|})/2`, a matriz de
multiplicação por `f = sum_k f_k T_k` na base de Chebyshev é

```
M_f = Toeplitz + Resto,
Toeplitz[i,j] = f_{|i-j|}/2 + delta_{ij} f_0/2,
Resto[i,j]    = f_{i+j}/2   - delta_{i0}  f_j/2.
```

As duas correções não são cosméticas:

- a **diagonal** `delta_{ij} f_0/2` existe porque em `i=j` os dois termos do
  produto coincidem;
- a **primeira linha** `-delta_{i0} f_j/2` corrige o *aliasing do índice 0*:
  para `i=0` os casos `k = n+i` e `k = n-i` são o **mesmo** `k`, e a contagem
  ingênua os soma duas vezes.

Sem essas duas, a decomposição erra — e errava: verificada agora, bate
**exatamente** em todos os casos testados.

### 5.19.2 O Toeplitz tem símbolo `f` na elipse

Com a diagonal incluída, o símbolo da parte Toeplitz conjugada pelo peso é

```
sigma(z) = (1/2) sum_d f_{|d|} (rho z)^d + f_0/2  =  f(xi(z)),
xi(z) = (rho z + 1/(rho z))/2,
```

isto é, **exatamente `f` avaliada na elipse de Bernstein**. Logo
`||Toeplitz|| <= sup_elipse |f|`, e a medida confirma: a razão
`||Toeplitz||/sup` fica em `0,68–0,99`, nunca acima de 1.

### 5.19.3 Toda a constante está no Resto

| `f` | `||Toep||/sup` | `||Resto||/sup` | `||M||/sup` |
|---|---:|---:|---:|
| `T_1` | 0,9814 | **0,5454** | 1,2015 |
| `T_2` | 0,9299 | **0,5675** | 1,2101 |
| `T_5` | 0,6834 | 0,5169 | 1,0522 |
| `1+0,5T_1` | 0,9937 | 0,1848 | 1,0023 |
| `0,7^k, k<=8` | 0,8549 | 0,2009 | 0,8596 |
| `0,4^k, k<=12` | 0,9696 | 0,1543 | 0,9734 |
| alternado `0,6` | 0,9071 | 0,1903 | 0,9114 |

> **`||M_f|| <= (1 + 0,57) · sup_elipse|f|`** pela desigualdade triangular, e
> `<= 1,21 · sup` observado.

E o padrão é nítido: o Resto pesa `~0,55` para um `T_j` **isolado** e `~0,18`
para coeficientes que **decaem**. O fundo `RefA` é do segundo tipo — seus
coeficientes caem por várias ordens de grandeza. A constante efetiva para ele
deve ficar perto de `1,2`, não de `1,57`.

### 5.19.4 O que falta

Provar um bound para `||Resto||`. Ele é Hankel mais uma correção de posto 1, e
operadores de Hankel têm teoria própria (Nehari, Hartman) que dá a norma pela
distância a `H^infinito`. Para `f` analítica numa elipse maior, essa distância
é pequena — o que explica por que o Resto encolhe quando os coeficientes
decaem. **É o último passo da cadeia, e é o único que ainda não tem forma.**

A cadeia completa, para registro:

```
R  <=  R_J + ||B||_2
     <=  0,08885  +  C · [ sup_elipse do termo linear + 2 · sup_elipse do fundo ]
```

com `R_J` **certificado** (§5.15), os supremos **medidos** (§5.17) e `C`
**estruturado mas não provado** (aqui).

## 5.20 O último elo: o Resto tem forma fechada (22/09/2026)

`scripts/certify_background_multiplier.py`,
`build/spectrum/background-multiplier.json`.

### 5.20.1 A primeira linha do Resto é nula

Da §5.19, `Resto[i,j] = f_{i+j}/2 - delta_{i0} f_j/2`. Em `i = 0` isso é
`f_j/2 - f_j/2 = 0`. **O Resto é um Hankel puro com a primeira linha
removida** — verificado a `0,0e+00` em todos os casos.

### 5.20.2 E a norma de Frobenius tem forma fechada

```
||Resto||_F^2  =  (1/4) · (rho^2/(rho^2 - 1)) · sum_u f_u^2 (rho^u - rho^-u),
```

com `rho^2/(rho^2-1) = 25/9` em `rho = 5/4`. Verificada a nove dígitos contra o
cálculo direto, em quatro famílias de coeficientes, inclusive aleatórios.

E `||.||_2 <= ||.||_F` é apertado justamente onde importa: a razão vale `1,000`
e `1,005` para coeficientes que decaem, contra `1,58` para um `T_j` isolado.

> **Não é preciso uma constante universal `C`.** Basta o bound para *este* `f`,
> o fundo — e ele é uma soma finita sobre coeficientes publicados.

### 5.20.3 Aplicado ao fundo

| quantidade (peso RT, 4 campos de `RefA`) | valor |
|---|---:|
| `sup` no domínio complexificado | 1,8487 |
| `||Resto||_F` pela forma fechada | 0,6837 |
| **`||M_fundo||_2 <= sup + Resto`** | **2,5324** |
| norma `l^1` ponderada (o que se usa hoje) | 12,4024 |
| **ganho** | **4,90×** |

### 5.20.4 A cadeia completa de F3

```
R  <=  R_J  +  ||B||_2
    <=  0,08885  +  [ 1,511  +  2 · 2,5324 ]
    =   0,08885  +  6,576
    =   6,66
```

| elemento | origem | estado |
|---|---|---|
| `R_J < 0,08885` | 1-D, `LDL^T` exato, recursão escalar, invariante | **certificado** |
| `1,511` (parte linear `S Gamma1`) | `ALGEBRAIC_TAIL.md` | já provado |
| `sup = 1,8487` | soma finita sobre `RefA` | certificável |
| `||Resto||_F = 0,6837` | forma fechada, soma finita | certificável |

> **`R <= 6,66` contra os 414 do ledger: fator 62.** E o valor medido extrapola
> para `~3`, então o certificado está a `2,2×` do real.

A região descoberta de F3 encolhe de `1 < Re s < 414` para `1 < Re s < 6,7`.

### 5.20.5 As três ressalvas que restam

1. **A identificação do termo bilinear.** `ALGEBRAIC_TAIL` limita a contribuição
   bilinear por `2 · 10,5`, e eu troquei `10,5` pelo bound `l^2` do multiplicador
   (2,5324). Isso supõe que o `10,5` de RT é a norma do multiplicador e não outra
   coisa — a norma `l^1` medida do fundo é 12,40, próxima mas não igual. **A
   identificação precisa ser conferida contra a derivação de RT.**
2. **Mistura de normas entre direções.** Somei em Fourier com peso `l^1` e usei
   Frobenius em Chebyshev. É conservador, mas é uma escolha, não um teorema.
3. **Tudo em ponto flutuante.** As duas somas são finitas e explícitas; refazê-las
   em racionais é mecânico, e é o que falta para dizer "certificado" sem aspas.

Nenhuma das três é pesquisa. A cadeia inteira de F3 passou de "estimativa de
energia quantitativa: **aberta**" para quatro elementos explícitos, um deles já
certificado e três que são somas finitas.

## 5.21 Continuando de onde a literatura está: campos de valores não convexos (22/09/2026)

O levantamento de `ESTADO_DA_ARTE.md` mostrou que nosso campo de valores falha
no contorno **por convexidade** (§5.9.1), e que existe uma família inteira de
enclausuramentos **não convexos** construída exatamente contra essa limitação:
*quadratic numerical range* e *block numerical range* para matrizes de
operadores em blocos (Tretter e colaboradores), e *pseudo numerical range* para
funções de operador (Gerhat & Tretter, 2022).

Era a continuação natural: nosso operador tem estrutura de blocos (quatro
componentes, e qualquer corte por índice), e o QNR **pode ter componentes
separadas** — o contorno poderia passar entre elas.

### 5.21.1 O teste

Para um bloco `[[A,B],[C,D]]`, o enclausuramento do QNR é uma região de
**Cassini**: `s` está fora se

```
dist(s, W(A)) · dist(s, W(D))  >  ||B|| · ||C||.
```

Testado em 8×24, partindo por índice radial em `n0 = 2, 4, 6, 8, 10, 12, 16`:
em **todos**, `s = 1/8` está **dentro** de `W(A)` e de `W(D)` — as distâncias
saem negativas, e o critério nem chega a ser avaliado.

### 5.21.2 Por que não podia funcionar, e por quanto

O espectro em torno da borda, medido:

| lado | distância de `s = 1/8` |
|---|---:|
| direita (raiz de gauge, em `mu`) | **0,04139** |
| esquerda | 0,12486 |
| **produto** | **0,005168** |

O critério de Cassini exige `||B||·||C|| < 0,005168`. O menor valor medido
entre todos os cortes foi **8,31**.

> **Falta um fator 1608.**

E isso **não depende do corte**: a exigência é sempre o produto das distâncias,
e a distância à direita está fixada pela geometria — a raiz de gauge está em
`s = mu` e a borda em `1/8`, a 0,0414 (0,0433 no limite). Qualquer
enclausuramento do tipo Cassini herda essa exigência.

### 5.21.3 O que isso acrescenta ao mapa

É a **mesma obstrução de sempre**, agora vista numa terceira forma. A §5.9
mostrou que o campo de valores convexo cobre o contorno; aqui, os não convexos
também cobrem, e o déficit é quantificado em 1608×.

Entra no mapa de onde não ir, com constante:

| enclausuramento | por que falha |
|---|---|
| campo de valores clássico | convexo, e o espectro está dos dois lados |
| quadratic / block numerical range | não convexo, mas exige `||B||·||C|| < 0,0052`; o menor medido é 8,31 |
| pseudo numerical range | herda a mesma exigência de produto de distâncias |

**A geometria manda:** o contorno passa a 0,0414 de uma raiz, e essa distância
entra multiplicativamente em qualquer critério de enclausuramento. É o mesmo
`1/d` do critério da §5.8.6, vestido de Cassini.

## 5.22 S4 medido: a deflação vale 2,6×, não 6,3× (22/09/2026)

Ataque a S4 começado pela pergunta prática: quanto vale, em condicionamento,
deflacionar o modo de gauge? A resposta corrige uma projeção minha.

### 5.22.1 O modo de gauge é analítico, não numérico

`GAUGE_RT_ACTION.md` estabelece mais do que eu supunha. O gerador é
`h = (Z+1)omega` com `s = mu`, vindo da **ação do difeomorfismo** (translação
nula), e vale a identidade, para qualquer fundo regular `w`:

```
D Omega(w)[ e^{mu tau}(Z+1) w ] = e^{mu tau}(Z+1) Omega(w).
```

Em solução, `Omega(w) = 0`, logo `(Z+1)w` está no núcleo do linearizado em
`s = mu`. **O autovetor de deflação é conhecido em forma fechada.** Como o
documento diz: "não se trata de identificar uma raiz por proximidade numérica:
é a ação induzida pelo difeomorfismo."

### 5.22.2 A medida

Projeção de Riesz do modo de gauge, operador `(A+s)Q0`, pesos RT, borda `1/8`:

| malha | raiz de gauge | `d` até 1/8 | `||V||` | deflacionado | ganho | `||P||` |
|---|---:|---:|---:|---:|---:|---:|
| 6×18 | 0,144011 | 0,0190 | 3822,1 | 910,7 | 4,20× | 105,1 |
| 8×24 | 0,166389 | 0,0414 | 1636,2 | 603,3 | 2,71× | 71,3 |
| **12×36** | **0,168309** | **0,0433** | **1846,2** | **703,4** | **2,62×** | **79,6** |

E deflacionando também a raiz de S3, em 12×36:

| cenário | `||V||` | ganho |
|---|---:|---:|
| nada | 1846,2 | 1,00× |
| só gauge (S4) | 703,4 | **2,62×** |
| só S3 | 1986,5 | **0,93×** — piora |
| gauge + S3 | 681,1 | **2,71×** |

### 5.22.3 Retratação da projeção

Eu havia projetado **6,3×** para S4 e **13,5×** para S4+S3, a partir do modelo
`||V|| ~ C/d` com a distância à raiz mais próxima. **Errado, por dois motivos:**

1. o modelo supõe **normalidade**, e os modos são fortemente não-normais —
   `||P_gauge|| = 79,6`, de modo que a deflação não remove o polo de forma
   limpa;
2. sob não-normalidade, `||V||` **não é governado pela raiz mais próxima**: a
   raiz de S3, a `0,276` da borda, contribui praticamente nada, e deflacioná-la
   sozinha **piora**.

Um controle confirma o mecanismo: em `s = 1`, longe de qualquer raiz, a
deflação dá ganho `0,58×` e `0,85×` — **piora**, como tem de ser quando não há
polo a remover.

### 5.22.4 Dimensionamento corrigido

| cenário | `||V||` | DOFs vs 107×384 |
|---|---:|---:|
| hoje | 2573 | 9794× |
| S4+S3 deflacionados (2,71× medido) | 949 | **3616×** |
| + `beta` no piso de Perron | 949 | 540× |
| + absorver `B` | 949 | 360× |

contra os `726×` e `111×` que eu havia projetado.

### 5.22.5 O que isto significa

S4 **entra na mesma faixa de 1,0× a 2,7×** de todas as outras alavancas de
constante. Não é a exceção que eu tinha anunciado.

E isso é mais uma confirmação da regra de `METODO.md`: deflacionar o modo de
gauge, do ponto de vista do condicionamento, é **melhorar uma constante dentro
da formulação fixa** — e cai na mesma faixa que reponderar, mover a borda ou
trocar o ponto-base.

**O valor real de S4 continua sendo outro**: ele é obrigação para `T2`, e sem
ele a contagem de 3 não vira 2. Isso não mudou. O que mudou é que ele **não é
atalho de custo**, e não deve ser priorizado por esse motivo.

## 5.23 A cadeia de F3, fechada computacionalmente (22/09/2026)

`scripts/certify_f3_chain.py`, `build/spectrum/f3-chain.json`.

### 5.23.1 As duas somas, em forma certificável

**(a) `||Resto||_F`: exata.** A forma fechada da §5.20 sobre coeficientes
diádicos (`TwoExp -70`) sai em aritmética racional, sem aproximação:

```
||Resto||_F = 0,683719
```

**(b) `sup` no domínio complexificado: grade mais Lipschitz.** Um supremo sobre
um contínuo não sai de amostragem. A cota rigorosa usa

```
sup <= max_na_grade + L_theta·(dtheta/2) + L_phi·(dphi/2),
L_theta <= sum |c_mn| m kappa1^m (kappa2^n + kappa2^-n)/2   = 16,6
L_phi   <= sum |c_mn| kappa1^m n (kappa2^n + kappa2^-n)/2   = 87,6
```

com as duas constantes de Lipschitz calculadas **exatamente** dos coeficientes.
Refinando:

| grade | máx. na grade | folga | `sup <=` | `R <=` |
|---:|---:|---:|---:|---:|
| 360 | 1,849868 | 0,9086 | 2,7585 | 8,48 |
| 1440 | 1,850108 | 0,2272 | 2,0773 | 7,12 |
| **5760** | **1,850137** | **0,0568** | **1,9069** | **6,78** |

O máximo na grade estabiliza em `1,85014` — a folga restante é grossura de
grade, e cai linearmente.

### 5.23.2 A cadeia

```
R  <=  R_J  +  [ 1,511  +  2 · (sup + ||Resto||_F) ]
    <=  0,08885  +  [ 1,511  +  2 · 2,5906 ]
    =   6,78
```

| elemento | valor | origem |
|---|---:|---|
| `R_J` | < 0,08885 | **certificado exato** (§5.13–5.15) |
| parte linear `S Gamma1` | 1,511 | `ALGEBRAIC_TAIL`, já provado |
| `sup` do fundo | <= 1,9069 | grade + Lipschitz exato |
| `||Resto||_F` | 0,683719 | forma fechada exata |

> **`R <= 6,78` contra os 414 do ledger: fator 61.**

### 5.23.3 Teste de sanidade

A cadeia dá `||B||_2 <= 1,511 + 2·2,5906 = 6,69`, contra os **4,90** medidos e
extrapolados na mesma norma (§5.18.3). **O bound está acima do valor real**, com
36% de folga — que é a condição necessária que a tentativa de `C = 1` havia
violado.

### 5.23.4 O que falta para dizer "certificado" sem aspas

1. **Aritmética intervalar na grade.** As constantes de Lipschitz e o
   `||Resto||_F` já são exatos; a avaliação nos 5760² pontos é em ponto
   flutuante. Trocar por intervalos é mecânico e não muda a estrutura.
2. **A identificação do termo bilinear.** `ALGEBRAIC_TAIL` limita a contribuição
   bilinear por `2·10,5`, e aqui `10,5` foi substituído pelo bound `l^2` da
   multiplicação, `2,5906`. Isso supõe que o `10,5` de RT é a norma desse
   multiplicador. Evidência a favor: a norma `l^1` medida do fundo é **12,40**,
   e `10,5 < 12,40` é compatível com um bound refinado do mesmo objeto. Mas os
   fatores `S` e `Xi` que acompanham `Gamma2` têm normas próprias — `Xi` é
   multiplicação por `T_1`, que em `l^2(rho^n)` vale 1,23 — e um acerto de
   contas completo precisa passar por (34)–(38) de RT.
3. **A redação como estimativa de energia**, que é a forma que a linha F3 do
   ledger pede.

Nenhum dos três é pesquisa. O primeiro é mecânico, o segundo é leitura atenta
de uma derivação publicada, o terceiro é redação.

## 5.24 Os três itens, feitos (22/09/2026)

### 5.24.1 Item 2 primeiro: a identificação, resolvida no fonte de RT

O `arXiv-1203.3766v1.tar.gz` no `.cache` contém `choptuik.tex`. A estimativa
para `Gamma2` está lá, e **não é a norma do fundo**: é um

```
||S Xi Gamma2(v, ·)||  <=  max sobre SETE combinações explícitas
                            das componentes de v, cada uma em ||·||_SSPACE,
||v||_SSPACE = ||v1|| + (1/2)||v2±Pv2|| + (1/2)||v3±Pv3|| + (1/2)||v4±Pv4||.
```

Avaliando as sete sobre `RefA`:

| combinação | 1 | 2 | 3 | 4 | 5 | 6 | 7 | **máx** |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `SSPACE` | 1,630 | **8,605** | 4,014 | 0,870 | 1,134 | 5,648 | 5,480 | **8,605** |

**O `10,5` de `ALGEBRAIC_TAIL` é uma cota conservadora do `8,605` que a fórmula
de RT dá sobre `RefA`** — 22% de folga. A identificação está confirmada, e
agora sabe-se exatamente o que substituir.

### 5.24.2 A substituição fiel

Não se troca "a norma do fundo" por um número: troca-se **cada `Space(·)` da
fórmula de RT** pelo bound `l^2` da multiplicação daquela combinação, e toma-se
o mesmo máximo. A derivação de RT reduz o bilinear a convoluções, e o que muda
é só a norma em que cada convolução é medida.

| combinação | `SSPACE` (`l^1`) | análogo `l^2` | razão |
|---|---:|---:|---:|
| 1 | 1,630 | 0,482 | 0,296 |
| **2** | **8,605** | **2,227** | 0,259 |
| 3 | 4,014 | 0,714 | 0,178 |
| 6 | 5,648 | 1,353 | 0,240 |
| 7 | 5,480 | 1,335 | 0,244 |
| **máximo** | **8,605** | **2,227** | **0,259** |

A razão é **estável em 0,18–0,30** nas sete, o que é o sinal de que a
substituição é estrutural e não um acidente de uma combinação.

> **`R <= 0,08885 + 1,511 + 2 · 2,2273 = 6,05`**

### 5.24.3 Item 1: o arredondamento é irrelevante, com cota

A avaliação na grade é em ponto flutuante. Cota explícita do erro acumulado,
com `u = 2^-53` e `n_ops = nm + nn + 2 = 144`:

```
erro <= 20 · n_ops · u · (soma dos módulos ponderados)
      = 20 · 144 · 1,110e-16 · 12,4024  =  3,97e-12
```

O fator 20 é generoso para cobrir `cos`, `arccos`, `exp` e acúmulo. Comparado
com a folga de Lipschitz em `G = 5760`, que é `5,68e-2`, o arredondamento é
**10¹⁰ vezes menor**. Aritmética intervalar mudaria a décima segunda casa.

### 5.24.4 Item 3: a estimativa de energia

A linha F3 do ledger pede "estimativa de energia quantitativa". Ei-la, no
produto interno `l^2` com os pesos:

> Para todo `u` no domínio,
> ```
> Re <(J_mu + B) u, u>  >=  -R ||u||^2,     R = R_J + ||B||_2  <=  6,05.
> ```
> Logo, para `Re s > R`,
> ```
> Re <L(s) u, u>  =  Re <(J_mu+B)u, u>  +  (Re s)||u||^2  >=  (Re s - R)||u||^2  >  0,
> ```
> donde `L(s)u ≠ 0` para `u ≠ 0` e, por Cauchy–Schwarz,
> `||L(s)u|| >= (Re s - R)||u||` — **injetividade e cota de resolvente de uma vez**.

E a passagem a `l^1`, que é onde o projeto vive, é gratuita na direção
necessária: `l^1(w) ⊂ l^2(w)` com norma 1, então toda autofunção em `l^1` é
autofunção em `l^2`, e excluir espectro em `l^2` exclui em `l^1`.

### 5.24.5 Estado

| elemento | valor | estado |
|---|---:|---|
| `R_J` | < 0,08885 | exato, racional, com cauda |
| parte linear `S Gamma1` | 1,511 | RT, publicado |
| máx. das 7 combinações, `l^2` | 2,2273 | grade + Lipschitz exato + cota de arredondamento |
| **`R`** | **<= 6,05** | **contra 414 no ledger: fator 68** |

As aspas de "certificado" caem em tudo menos num ponto: a avaliação da grade
usa `float` com **cota de erro explícita** em vez de intervalos. É a diferença
entre `3,97e-12` provado por cota e `0` provado por construção, num número cuja
folga de Lipschitz é `5,68e-2`.

## 5.25 A cadeia de F3 mistura convenções de norma (23/09/2026)

Achado ao reutilizar a máquina de F3 para o operador de constraints `K`, não
ao reler a derivação. **`R <= 6,05` tem dois passos injustificados.** A cota é
plausivelmente verdadeira — os valores medidos ficam bem abaixo —, mas não está
provada.

### 5.25.1 Primeiro passo: as parcelas estão em produtos internos diferentes

`R <= R_J + ||B||_2` só vale se as duas parcelas estiverem no mesmo produto
interno. Com `w_n = (2-[n=0]) kappa2^n`:

| parcela | convenção |
|---|---|
| `R_J = 0,0888` | `D = diag(sqrt(w))` — o pencil `Herm(WJ)y = lambda W y` |
| `||Resto||_F`, forma fechada | também `sqrt(w)`: com `D = diag(rho^n)`, a soma `sum_{i=1..u} rho^{2(2i-u)}` reproduz a fórmula **exatamente** quando `rho^2 = kappa2` |
| sup do fundo | elipse de parâmetro `kappa2` — a de `sqrt(w)` é a de `sqrt(kappa2)`, menor; pelo máximo do módulo, é **conservador**, válido |
| **parte linear `1,511`** | **`w` inteiro** (e ℓ¹): é `(17/100)(80/9)` de RT |

`Gamma1` não divide por ξ — é combinação de `I` e da reflexão `P`, isometria
em Chebyshev. A divisão vem de `S`: as componentes 1, 2, 4, 5 de `Gamma1 v`
estão na imagem de `(1-P)` (RT, observação após a eq. do sistema), e `S` as
divide por ξ. Logo `mu S Gamma1` é do tipo `mu · 2 R_x`, e `80/9 = 2 · 40/9`.

`R_x` é Toeplitz triangular com entradas `±2`; o sinal alternado alinha todos os
termos do símbolo em `theta = pi/2`, então a norma ℓ² é **igual** à soma
geométrica `2 rho^-1/(1 - rho^-2)`:

| `rho` | `||R_x||` | `mu_max · 2 ||R_x||` |
|---|---:|---:|
| `5/4` (w inteiro) | 40/9 = 4,44 | **1,511** |
| `sqrt(5/4)` (convenção de `R_J`) | 8,94 | **3,04** |

Na convenção em que `R_J` e o `Resto` foram calculados, a parte linear é ~3,0,
não 1,511.

### 5.25.2 Segundo passo: a estrutura de blocos de ℓ¹ não se transplanta

A "substituição fiel" (§5.24.2) trocou cada `Space(·)` da fórmula de RT por um
bound ℓ² **e manteve o máximo sobre as sete combinações**. Esse máximo é a
estrutura de ℓ¹, em que a norma da entrada é a **soma** das normas das
componentes. Em ℓ² a norma da entrada é `sqrt(sum ||v_k||^2)`, e o análogo
correto é a norma espectral da matriz de normas de blocos — que pode exceder o
máximo. Se as combinações agem em entradas distintas, o análogo é da ordem de
`sqrt(sum c_k^2) ≈ 3,1`, contra 2,227. Não foi refeito.

### 5.25.3 A checagem barata

Valores verdadeiros na truncação, ponto flutuante:

| malha | convenção | `máx Re W(-A)` | `||D B D^-1||_2` | `R_J` livre |
|---|---|---:|---:|---:|
| 12×36 | `sqrt(w)` | **2,54** | **3,77** | 0,0888 |
| 12×36 | `w` inteiro | 2,63 | 4,14 | −0,121 |

A cadeia usou `||B||_2 <= 5,97`; o valor na truncação é 3,77. Então a cota
**provavelmente vale** para o operador infinito, com folga. O que caiu foi a
prova.

### 5.25.4 Estado e reparo

- `R_J < 0,08885` continua certificado, na sua convenção.
- `R <= 6,05` passa de "cadeia com uma peça certificada" a **"cota plausível com
  dois passos injustificados"**. Uma correção grosseira dá ordem de 9 —
  **ainda não rigorosa**, porque o passo 2 não foi refeito. O fator contra 414
  cai de 68 para ~45; o resultado qualitativo sobrevive.

**Reparo recomendado: tudo em `D = diag(w)`, pesos inteiros.** Nessa convenção
o sup na elipse `kappa2` é exato, a parte linear 1,511 é exata (norma de Toeplitz
= soma), e `R_J` é **negativo** (−0,121 em ponto flutuante; certifica-se pela
mesma recursão com `kappa2 -> kappa2^2 = 25/16`). Falta:

1. refazer `||Resto||_F` para `rho = kappa2`:
   `||Resto||_F^2 = sum_u (f_u^2/4) (kappa^{2u+4} - kappa^{4-2u})/(kappa^4 - 1)`;
2. trocar o máximo sobre as sete combinações pela norma espectral da matriz de
   normas de blocos, com a estrutura de componentes de `S Xi Gamma2` e de `S Gamma1`.

**A lição, para `METODO.md`:** ao trocar ℓ¹ por ℓ², não basta trocar cada norma.
Toda estrutura que dependia da **soma** sobre componentes — máximos sobre
combinações, somas de colunas — precisa ser refeita.

### 5.25.5 Reparo feito (23/09/2026): `R_L <= 10,92` e `R_K <= 5,06`, exatos

`scripts/certify_numerical_range_l2.py`, `build/spectrum/numerical-range-l2.json`.
Substitui `certify_f3_chain.py`.

**O espaço.** Os coeficientes de RT são os **simétricos** (índice `n ∈ ℤ`) — por
isso o operador livre acopla `2μn` uniformemente, inclusive a linha `n = 0`. O
ℓ² certo é o simétrico, `W(k) = κ1^|k1| κ2^|k2|`, que em coordenadas reais é
`D = diag(sqrt(c_m c_n) κ1^m κ2^n)`. Nele: `l^1(w) ⊂ l^2(W)` com norma ≤ 1;
multiplicação é convolução em ℤ² com peso submultiplicativo, e **Young** dá
`||M_f|| <= field_norm(f)`, a norma já auditada; `||R_x|| <= 2κ^-1/(1-κ^-2)`
(40/9 e 144/17), com somas de linha e de coluna iguais; e os pedaços
(componente, paridade) são ortogonais, então os blocos se combinam por norma
**espectral**, não por máximo.

**A parte livre, com cauda uniforme.** Pencil `Herm(W̃J) - λ0 W̃`,
`W̃ = diag(c_n κ2^(2n))`, racional. Recursão escalar exata no bloco finito e,
para a cauda, um argumento que vale para **todo** `k`:

```
gamma_{k+1} = rho_k [gamma_k + q_k(1-gamma_k)^2/(q_k(1-gamma_k)+1)]  <  rho_k   se gamma_k < 1,
rho_k = (1 + s/n_k) κ2^(-2s)  decrescente,  < 1 a partir de N0.
```

O colchete é `γ + (1-γ)t` com `t ∈ (0,1)`, logo menor que 1. Então `γ_k < 1` para
sempre e todo pivô `e_k(q_k(1-γ_k)+1)` é positivo. **Isso também fecha uma lacuna
da §5.15**: lá a invariância de `[0,79; 0,80]` foi verificada em `k` amostrados,
não para todo `k`.

`λ_min` escala linearmente com `μ` (a derivada temporal é antissimétrica e some
no Hermitiano), então certifica-se `μ = 1` e aplica-se `μ ∈ [1/6, 17/100]`.

| modelo | `λ0` (μ=1) | `γ_N0` | `ρ_N0` | `R_livre` |
|---|---:|---:|---:|---:|
| `L`, passo 1 (`v2,v3,v4`) | +0,57787 | 0,650 | 0,656 | −0,0963 |
| `L`, passo 2 (`v1`) | +0,83056 | 0,419 | 0,420 | −0,1384 |
| `K`, passo 1 (`c2`) | −0,37515 | 0,786 | 0,810 | +0,0638 |
| `K`, passo 2 (`c1`) | +0,20630 | 0,634 | 0,640 | −0,0344 |

**O resultado**, com `R <= λ_max(diag(r) + Sym(η C η^-1)) + correção`, pesos de
componente `η` otimizados e a desigualdade final por `LDL^T` exato:

> **`R_L <= 10,92`** (F3, contra 414: fator 38) e **`R_K <= 5,06`**.

**Validação na mesma convenção**, truncação 12×36: `R` verdadeiro 2,35 (`L`) e
1,006 (`K`); cada bloco verdadeiro fica abaixo do majorante (pior razão 0,38 e
0,51; nenhuma violação); e a parte livre verdadeira de `L` vale −0,0973, que é
exatamente `0,57787 × 0,1683` — **cinco dígitos de concordância, confirmando a
convenção**.

**O certificado é inteiramente exato.** O ponto flutuante só propõe `λ0` e `η`;
tudo o que se afirma é verificado em racionais. A ressalva de ponto flutuante
da cadeia antiga desaparece, porque Young dispensa a avaliação em grade.

**O preço é a folga:** 10,92 contra 2,35 verdadeiro. A cota antiga usava o sup na
elipse, 0,18–0,30 da norma ℓ¹, mas com convenções misturadas. Refazê-la
corretamente em 2D — Toeplitz mais Hankel, com o resto de Fourier por
desigualdade triangular e não por média — é refinamento opcional. Nenhum uso
atual de F3 depende do valor exato de `R`: ele só fecha a região a contar.

## 5.26 O espaço que torna o fundo pontual (24/09/2026)

`S3_CONDICIONAMENTO_SHARP.md` §9. No L² do toro de raios `(κ1, κ2)`, multiplicação,
reflexão `ξ ↦ −ξ` e divisão por ξ são todas pontuais, e a borda do campo de valores
do fundo vira um `sup` de `λ_max` de uma matriz pequena. Para `K`: `R♯` de 5,06
(majorantes por bloco) para **1,765** (símbolo), certificado em Arb, contra 1,03
verdadeiro.

**Aplicada a `L` (F3), 24/09.** A forma pontual é 7×7 nos valores
`(x1(z), x2(±z), x3(±z), x4(±z))`, com `B_L = μSΓ1 + 2SΞΓ2(ω,·)` lido do TeX de RT
(fórmula `skdfhkfdjhdkhjdgfdjh33hk`) e `μSΓ1 x = μ(0, (x1 + 2Ox2)/ξ, −x1/ξ, 2Ox4/ξ)`
— que reproduz o `gamma1_majorant` auditado. Toro de raios `(65/64, 5/4)`, grade
2048×4096 em quatro faixas no PC do laboratório: **16.777.216 centros verificados
em Arb, sem falha**, folga 0,077.

| | cota | verdadeiro (truncação 12×36) |
|---|---:|---:|
| `max Re W(−B_L)` | 2,9492 | 2,72 |
| `||B_L||` | 3,9641 | 3,69 |
| parte livre (toro) | −0,0232 | −0,0234 |
| **`R_L`** | **<= 2,926** | 2,28 |

> **F3: `R <= 2,93`**, contra 414 do ledger — fator 141 —, certificado.

A cota fica 5% acima da verdade na truncação, e a verdade ainda sobe com a malha.

**Natureza do rigor.** Os valores nos centros são bolas de Arb, a positividade é
verificada por `LDL^T` em bolas, e as constantes de Lipschitz dos campos são
racionais exatos. A folga final — norma de Frobenius de matrizes 7×7 montadas a
partir dessas constantes — é calculada em ponto flutuante e multiplicada por
1,001: a margem fica ~10¹² vezes acima do erro de arredondamento. É prova por cota
a priori, não por construção; passar essa etapa para Arb é mecânico.
**Feito em 24/09** (`scripts/symbol_slack_arb.py`): a folga refeita em bolas fica abaixo
da gravada, para L (0,077120/0,090042 contra 0,077197/0,090132) e para K
(0,053525/0,053616 contra 0,053578/0,053670). Cotas em Arb: `R_fundo,L <= 2,949121`,
`||B_L|| <= 3,964043`, `max Re W(−B_K) <= 1,439525`, `||B_K|| <= 1,649617`.
`scripts/symbol_forms.py`, `scripts/certify_symbol_numerical_range.py`,
`build/spectrum/symbol-L-2048x4096.json`.

**A lição de método:** quando os operadores de um problema são todos
"multiplicação, reflexão e divisão", procurar o espaço em que eles comutam com
a avaliação pontual. Majorar bloco a bloco perde todos os cancelamentos, e aqui
isso custava um fator 3.
