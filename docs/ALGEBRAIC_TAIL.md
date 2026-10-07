# Uma rota de cauda com pesos fixos

Esta derivação substitui a tentativa de ganho exponencial P2. É uma prova
em papel para a realização por fechamento do operador selecionado, sob as
hipóteses explicitadas abaixo. Não demonstra a equivalência desse domínio
com o quociente físico nem certifica os contornos existentes.

## Dados de partida

Use a complexificação dos graus de liberdade REAIS de RT, com norma
`sum_j w_j |v_j|`, agora com amplitudes complexas, pesos RT
`κ1=65/64`, `κ2=5/4`, e soma das quatro componentes. A norma de uma
matriz real nessa complexificação é a mesma máxima soma de colunas
ponderadas; os bounds reais de B se estendem sem fator adicional.
O operador selecionado é `L(s)=J_mu+B+sI`, onde
`J_mu=diag(K_mu,J_mu,J_mu,J_mu)` e
`B=mu S Gamma1+2 S Xi Gamma2(omega,·)`.

As fórmulas (22), (34)–(38) de
[RT](https://arxiv.org/pdf/1203.3766) fornecem o ponto de partida algébrico
e os bounds de multiplicação. Seus parâmetros publicados permitem tomar
`1/6 <= mu <= 17/100` e `||B|| <= 23`: a contribuição de `S Gamma1` é
no máximo `(17/100)(80/9)`, a contribuição bilinear no fundo `refA` é
no máximo `2·10,5`, e a correção é no máximo
`6(2^-25+2^-277)`. Essa soma é estritamente menor que 23.
O bound da correção usa (39b) e (60); sua auditoria computacional integral
permanece parcial neste repositório. As constantes acima dependem do
teorema publicado de RT, não apenas da inspeção de `RefA.dat`.

## Inversa livre deslocada

Escreva `a=Re z >=0`, `b=Im z`. Nos blocos escalares, somar `zI`
às fatorizações de RT substitui a diagonal por

```
D_z(m,n) = im/2 + mu(n+1) + z.
```

Os fatores de deslocamento têm razões

```
f_J(m,n) = (-im/2 + mu*n - z) / D_z(m,n),
f_K(m,n) = (-im/2 + mu*(n+1) - z) / D_z(m,n).
```

Ambas têm módulo no máximo 1: as partes imaginárias do numerador e
denominador têm módulos iguais, e a diferença dos quadrados das partes
reais é não negativa quando `a>=0`. Para J, essa diferença é
`mu*(2*n+1)*(mu+2*a)`; para K, é `4*mu*(n+1)*a`.
Como `|D_z| >= mu+a`, não há denominadores nulos.

As inversas são `V[k,f] D_z^-1 (I-U_k)`, com `k=1` ou `k=2`, onde
`U_k` baixa o índice de Chebyshev em k e preserva o índice de Fourier.
Para estimar Q, passe temporariamente aos coeficientes complexos completos
de Fourier, com índices positivos e negativos independentes e módulo
complexo usual. Com `r=κ2^-1`, a série de deslocamentos de V tem norma
no máximo `1/(1-r^k)` e a diagonal, `1/(mu+a)`.
A mudança `(v_cos,v_sin)→(v_cos+i*v_sin,v_cos-i*v_sin)` tem produto
das normas de ida e volta no máximo 2 (o modo zero tem constante 1).
Esse fator é aplicado uma única vez ao produto completo. Consequentemente,
para `Q(z)=(J_mu+zI)^-1`,

```
||Q(z)|| <= 2*(1+r)/(1-r)/(mu+a) = 18/(mu+a).
```

As fórmulas são identidades em séries de suporte finito e os inversos
preservam suporte finito. Defina o domínio como o fechamento do gráfico
de J nos polinômios admissíveis, na norma ambiente com pesos fixos.
As identidades de inversa e o bound acima se estendem a esse fechamento
por densidade; a bijetividade segue aproximando cada dado por polinômios.
Não se presume que esse fechamento seja automaticamente o domínio físico.

## Cauda de saída, sem hipótese sobre o adjunto

Fixe agora um ponto-base REAL `z0>=0`. Seja `Π` o corte
`|m|<M`, `n<N`, com `M,N>=1`. Os fatores V somente baixam n:
`(I-Π)VΠ=0`. É possível, portanto, inserir `I-Π` entre V e a diagonal
ao estimar `(I-Π)Q(z0)`. Na região descartada,

```
1/|D_z0(m,n)| <= max(2/M, 1/(z0+mu*(N+1))).
```

O primeiro caso é `|m|>=M`; o segundo, `n>=N`.
Aplicando os bounds dos três fatores e a mesma conversão de norma,

```
t(M,N) = 18 * max(2/M, 1/(z0+mu_min*(N+1))),
||(I-Π)Q(z0)|| <= t(M,N).
```

O índice `N+1` aqui é próprio da cauda de SAÍDA. Uma cauda de ENTRADA
`Q(I-Π)` atravessa `I-U_k` primeiro e pode perder até dois índices;
não se deve copiar esta fórmula para essa cauda.

Em particular `t(M,N)→0`. Q é limite em norma de operadores de posto
finito `ΠQ`, portanto compacto nessa realização com pesos fixos.

## Perturbação pelo fundo

Suponha `||B||<=beta`. Se `q0*beta<1`, com
`q0=18/(mu_min+z0)`, então

```
R0 = (L(0)+z0 I)^-1 = Q(z0) (I+B Q(z0))^-1,
||R0|| <= q0/(1-beta*q0),
||(I-Π)R0|| <= t(M,N)/(1-beta*q0).
```

Isso dá um bound uniforme de cauda unilateral do operador infinito
selecionado, condicionado aos bounds uniformes de mu e B, sem P2.
Também dá compactação de R0. A família afim é Fredholm de índice zero
pela fatorização por `I+(s-z0)R0`, compacto perturbando a identidade.
Esta última passagem usa teoria funcional padrão; não está formalizada
no kernel local, que hoje trabalha com matrizes racionais.

Com `mu_min=1/6`, `beta=23`, um ponto-base seguro é `z0=828`:
`q0=108/4969`, `beta*q0=2484/4969<1/2` e
`||R0||<=108/2485`. Para todo `Re s>=414`, a mesma estimativa livre
dá `beta*||Q(s)||<=2484/2485<1`, logo não há espectro nessa região
da realização selecionada. É uma borda muito mais conservadora que
`Re s=1` usada nos diagnósticos existentes.

## Limites desta conclusão

O bound está na norma COM PESOS RT. Os tiles históricos usam a norma de
colunas sem esses pesos. A opção nova `build_tile.py --rt-weights` conjuga
as matrizes pela diagonal dos pesos extraídos dos endereços dos DOFs,
permitindo que a norma de colunas emitida seja a norma ponderada correta.
Isso não vincula automaticamente os campos de bounds fornecidos pelo usuário
às constantes desta derivação.
A cauda unilateral acima é para `ΠR0`, cuja compressão não é a inversa
finita. C1-G continua exigindo uma ponte até o determinante implementado.
Também é preciso cobrir `1<Re s<414` e a região entre o eixo imaginário e
`sigma0=1/8`, ou obter estimativas melhores para excluí-las.
S3/S4 e a equivalência do domínio permanecem abertas.

`AlgebraicTail.lean` confere os valores racionais desta seção, não as
identidades de operadores, a compactação, a teoria de Fredholm ou as
hipóteses do teorema RT.

## Ponte alternativa para o determinante finito

A dificuldade C1-G da resolvente COMPLETA pode ser evitada se o
precondicionador for a inversa LIVRE `Q0=(J_mu+z0)^-1`.
Isto é outra fatorização e requer outro gerador de parametrizes;
`build_tile.py` ainda usa a resolvente finita completa.

Defina `delta=s-z0`, `D(s)=B+delta I` e

```
H(s) = L(s) Q0 = I + D(s) Q0,
H_N(s) = I + D(s) Q0 Π.
```

As séries de Q0 só baixam n e preservam m, portanto `Q0 ι = ι Q0_N` e
`Q0_N=(J_mu,N+z0 I)^-1` EXATAMENTE. O determinante de posto finito
de H_N é, por Sylvester,

```
det_V(I + Π D(s) Q0 ι)
 = det_V((A_N+sI) Q0_N)
 = det(A_N+sI) / det(J_mu,N+z0 I).
```

O denominador é constante e não nulo. Assim, este caminho reproduz o
winding do mesmo pencil finito sem identificar a inversa completa
com a sua compressão. Aqui `A_N` tem de ser o truncamento do fundo
EXATO: a diferença para o fundo numérico exportado ainda precisa ser
incluída em um enclosure ou em uma segunda homotopia.

Para a cauda de entrada, `I-U_k` pode baixar n em até 2. Com `N>=2`,
as mesmas fatorizações dão o bound seguro

```
t_in(M,N) = 18 * max(2/M, 1/(z0+mu_min*(N-1))),
||Q0(I-Π)|| <= t_in(M,N).
```

Se `d >= sup_Gamma |s-z0|` e `b=beta+d`, então
`||H-H_N|| <= b*t_in`. Seja `g` um bound uniforme da inversa da matriz
finita `I+ΠD(s)Q0ι`. Pela identidade de inversa de posto finito,

```
H_N^-1 = I - D(s)Q0ι (I+ΠD(s)Q0ι)^-1 Π,
||H_N^-1|| <= 1 + b*q0*g.
```

Portanto a condição quantitativa suficiente para transferência é

```
b*t_in(M,N) * (1+b*q0*g) < 1.
```

Ela garante que a homotopia entre H_N e H não encontre singularidades
no contorno. Essa condição ainda não foi satisfeita para as malhas
existentes. A cauda conservadora, o grande ponto-base e a distância d
podem torná-la impraticável; será necessário apertar os bounds, usar
uma decomposição em blocos ou aumentar os truncamentos.

## Experimento exato de 15/09/2026

`scripts/check_free_preconditioner.py` agora constrói `J_mu,N` diretamente
nas coordenadas reais exportadas. O bloco espacial tem diagonal `mu(n+1)`
e entradas `2*mu*n` nos graus inferiores; o bloco K preserva paridade.
A derivada temporal é o bloco real de multiplicação por `im/2`.
O código foi comparado com diferenciação de polinômios de Chebyshev,
testa a restrição exata da inversa livre e a identidade de determinantes.

Na matriz 4×12 inteira (138 coordenadas), com `mu=mu_refA`, `z0=0` e
centro real `s=1`, constrói `F_N=(A_N+I)Q0_N`, conjuga pelos pesos RT,
obtém uma inversa aproximada diádica e recalcula seu defeito com inteiros
e racionais. O defeito é menor que 1; há um bound finito da inversa.
O relatório exato está em `build/spectrum/free-preconditioner-4x12.json`.

Entretanto, `b*t_in=2592/11>1` antes mesmo de multiplicar por
`1+b*q0*g`. Todas as cinco malhas disponíveis falham nesse primeiro fator;
na maior, 12×36, ele ainda é `2592/35>1`. Essas contas estão verificadas
no kernel em `AlgebraicTail.lean`.

Com esse majorante, base 0 e distância 1, passar apenas no primeiro fator
exige `M>864` e `N>2593`: 865×2594 é a primeira escolha inteira que passa.
Isso NÃO garante o majorante completo, cujo fator de inversa é maior que 1.
Portanto não há justificativa para aumentar a malha às cegas: o gargalo
identificado é o majorante global. O próximo passo analítico é separar
blocos de frequência e estimar os acoplamentos, como a estratégia de RT,
em vez de substituir todos por `||B||<=23`.

Reprodução:

```bash
OPENBLAS_NUM_THREADS=2 .venv/bin/python scripts/check_free_preconditioner.py \
  build/spectrum/rt-A-4x12.dat \
  --out build/spectrum/free-preconditioner-4x12.json
```

O relatório não aceita um certificado físico: verifica apenas o centro
finito e mostra a falha da condição de transferência. O erro do fundo exato
e o restante do contorno continuam separados e não certificados.
