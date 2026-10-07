# Estimativas por blocos e o fechamento físico

Atualização: `STRUCTURED_EXTERIOR.md` melhora a cauda livre preservando
cancelamentos nas colunas. Em G=64×192, o majorante TT uniforme cai de
2619/191 para 7857/772, ainda maior que um. As estimativas abaixo são
mantidas como comparação conservadora.

## 1. Separação correta dos objetos

Trabalhamos com `H(s)=L(s)Q0=I+(B+(s-z0)I)Q0` na norma de coeficientes
com pesos RT. Sejam P um núcleo finito, G uma caixa maior,
S=G-P a coroa intermediária e T=I-G o exterior infinito.
Uma matriz exportada sobre G fornece blocos entre P e S.
Não fornece THT, PHT, THP nem os acoplamentos S/T.
Declarar que o bloco restante da matriz é a cauda infinita perderia esses
termos e produziria uma conclusão sem fundamento.

O núcleo e a coroa devem ser extraídos da MESMA matriz externa. Comparar
diretamente duas exportações não garante a mesma aproximação do fundo:
o gerador de RT corta o fundo conforme a dimensão solicitada.

## 2. Critério que mantém os acoplamentos separados

Para uma divisão finito/exterior e uma parametriz `R=diag(V,I)`, ponha
`E=I-RH`. Suponha bounds não negativos

```
||E_FF|| <= a,    ||E_FT|| <= b,
||E_TF|| <= c,    ||E_TT|| <= d.
```

Na norma `||x_F|| + r||x_T||`, a norma de E é majorada por
`max(a+r*c, d+b/r)`. Logo basta fornecer um racional positivo r com

```
a+r*c < 1,    d+b/r < 1.
```

Essa possibilidade equivale, para a,b,c,d não negativos, a

```
a<1,    d<1,    b*c < (1-a)*(1-d).
```

Com c>0, escolha r no intervalo aberto
`b/(1-d) < r < (1-a)/c`. Com c=0, basta r>b/(1-d).
O critério não exige b e c individualmente pequenos: importa seu produto.
Por exemplo, a=1/10, b=4, c=1/100, d=1/5 e r=10 dão contração 3/5,
enquanto r=1 falha. Essas desigualdades são conferidas no kernel em
`FrequencyBlocks.lean`.

Se V é invertível no bloco finito, a série de Neumann para I-E prova a
invertibilidade de H. No caso V é uma matriz quadrada candidata, o próprio
bound `||I-V H_FF||<1` prova sua invertibilidade finita.
Aplicar o critério em todo o contorno ainda exige bounds uniformes nos
segmentos e transferência da contagem por homotopia. Para uma homotopia
que desliga os acoplamentos e o termo de cauda, os mesmos bounds funcionam
desde que sejam estabelecidos uniformemente ao longo dessa homotopia.

Para três ou mais blocos, procure pesos r_i>0 tais que, para cada coluna j,
`sum_i r_i e_ij < r_j`. Trata-se de normas por COLUNA, compatíveis com ℓ¹.
Os e_ij devem incluir as incertezas do fundo e da variação no contorno.

## 3. Experimento: núcleo 4×12 dentro da caixa 6×18

Foi calculado `H(1)=(A_6x18+I)Q0_6x18`, com z0=0, mu de RefA e pesos RT.
Os 138 DOFs do núcleo e os 195 DOFs da coroa são selecionados por seus
índices, não por uma posição presumida na matriz.
As inversas propostas em ponto flutuante são arredondadas para diádicos;
os quatro blocos de I-RH são recalculados exatamente com inteiros.

| Parametriz | a | b: coroa→núcleo | c: núcleo→coroa | d |
|---|---:|---:|---:|---:|
| inversa do núcleo; identidade na coroa | 5,25e-10 | 4,407 | 17,126 | 8,958 |
| inversas separadas | 5,25e-10 | 4,407 | 22,660 | 7,22e-11 |
| inversa acoplada de núcleo+coroa | 5,64e-10 | 6,11e-11 | 7,87e-10 | 8,06e-11 |

Na segunda linha o problema não são as inversas diagonais: o produto
b*c está próximo de 100, contra margem próxima de 1. Nenhuma escolha de
peso faz esse majorante passar. Isso não prova singularidade da matriz;
prova que essa separação perde informação demais dos acoplamentos.
O programa também testa uma inversa acoplada do conjunto núcleo+coroa,
conferindo novamente os quatro blocos do defeito. Essa terceira opção
passou: até o majorante sem reponderação é menor que `1,36e-9`.
É uma validação racional finita no centro s=1, não uma certificação da cauda.

## 4. Núcleo→exterior: usar a cauda do fundo, não sua norma inteira

Se o núcleo mantém |m|<M_P,n<N_P e o exterior começa fora de
|m|<M_G,n<N_G, defina os limites de separação

```
g_m=M_G-M_P+1,    g_n=N_G-N_P+1.
```

Uma multiplicação de Fourier–Chebyshev só pode sair de G partindo de P
se o coeficiente do fundo tiver |m|>=g_m ou n>=g_n. Isso decorre dos
índices soma/diferença no produto de Chebyshev e da convolução de Fourier.
Reflexões e projeções de paridade não aumentam os índices.
O termo S Gamma1 somente baixa n e preserva m, logo `T S Gamma1 P=0`.
Usando o bound bilinear de RT já registrado em `ALGEBRAIC_TAIL.md`,

```
||T B P|| <= 6 ||omega_fora_da_caixa(g_m,g_n)||
          <= 6 (||RefA_fora_da_caixa(g_m,g_n)|| + epsilon_background).
```

Aqui `epsilon_background=2^-25+2^-277` é o bound PUBLICADO relativo a RefA,
não algo que o produtor deduz de um arquivo. A norma da cauda de RefA é
somada exatamente, com fatores de multiplicidade e pesos 65/64 e 5/4.
Como Q0 preserva P, isso implica
`||T H P|| <= 6(tail_ref+epsilon_background) ||Q0 P||`.

| Núcleo P | Caixa G | Bound de ||T B P|| |
|---|---|---:|
| 4×12 | 6×18 | 36,652 (pode ser substituído pelo global 23) |
| 4×12 | 12×36 | 1,227 |
| 4×12 | 24×72 | 0,001427 |
| 4×12 | 45×113 | <2e-7 |
| 12×36 | 24×72 | 0,140 |
| 12×36 | 53×137 | <2e-7 |

Os valores exibidos foram arredondados para cima. Os racionais completos
estão em `build/spectrum/background-couplings.json`.
Nas linhas em que a cauda de RefA zera, a incerteza do fundo exato permanece:
suporte finito da aproximação não implica suporte finito da solução.

## 5. O sentido contrário tem um termo adicional

Não se pode transpor o bound anterior: Q0 baixa graus de Chebyshev, então
pode levar dados do exterior para frequências baixas antes de B atuar.
Escolha uma caixa intermediária U entre P e G. Para T=I-G,

```
P D Q0 T = P D U Q0 T + P B (I-U) Q0 T,
D=B+(s-z0)I.
```

O segundo termo não contém delta pois `P(I-U)=0`.
O primeiro é limitado por `(beta+|s-z0|)||U Q0 T||`.
O segundo pode ser limitado por `||P B(I-U)|| ||Q0 T||`.
Para a parte bilinear de `P B(I-U)`, vale o mesmo argumento de separação
de frequências do fundo; para S Gamma1 há uma série descendente adicional.

Em z0 real não negativo, usando as séries de deslocamentos de Q0,
uma cota conservadora para a primeira passagem é

```
||U Q0 T|| <= 4/(mu_min+z0) *
              (4/5)^(N_G-N_U+1) / (1-4/5).
```

Somente a cauda em n contribui aqui, pois Q0 preserva m e M_U<=M_G.
A cota resulta do bound 4/(mu_min+z0) dos coeficientes da inversa,
somado com a razão dos pesos radiais; não é um bound para Q0 T inteiro.
O termo S Gamma1, por sua série `4 U_1 V_2`, tem cota
`4*mu_max*(4/5)^(N_U-N_P+1)/(1-(4/5)^2)` nesse acoplamento descendente.
Assim, uma estimativa utilizável é

```
||P D Q0 T|| <= (beta+|delta|)*q_UT +
  [6*(tail_ref(P,U)+epsilon_background) +
   4*mu_max*(4/5)^(N_U-N_P+1)/(1-(4/5)^2)] * t_in(G).
```

Essa cota foi avaliada com base 0 e |delta|<=5/4, suficiente para o
retângulo atual de A (Re s<=1, |Im s|<=1/4), embora a parametriz finita
acima só tenha sido verificada no centro real s=1:

| P | U | G | Bound de ||P D Q0 T|| |
|---|---|---|---:|
| 4×12 | 12×36 | 24×72 | 2,633 |
| 4×12 | 24×72 | 45×113 | 0,249 |
| 4×12 | 24×72 | 64×192 | 0,000808 |
| 12×36 | 32×96 | 64×192 | 0,000809 |

Os valores exibidos foram arredondados para cima. Ainda é necessário
multiplicar pelos blocos relevantes da parametriz finita para obter b.
Em particular, não se deve usar o bound direto de P→T também para T→P.

## 6. Consequência para a prova física

O resultado útil é uma escolha de estratégia: resolver núcleo e coroa
próxima de forma ACOPLADA, e usar a separação de índices para o exterior
distante. Ainda faltam bounds para S/T e T/T; ter THP pequeno não basta,
pois a coroa encosta no exterior. O critério por três blocos deve incluir
todos esses caminhos, além da variação complexa ao longo do contorno.

Mesmo uma transferência espectral completa contará o operador selecionado.
Para concluir um único modo físico, S3 deve excluir os modos que violam as
constraints, e S4 deve identificar o subespaço de gauge no domínio correto.
Os testes de resíduos e sobreposições existentes não substituem essas provas.

## 7. Controle uniforme completo da coroa e exterior

Se `|s-z0|<=rho` em todo o contorno, escreva `t_P=t_in(P)` e
`t_G=t_in(G)`. As projeções coordenadas têm norma um. Como
`S=(I-P)S`, `T=(I-G)T` e Q0 preserva G, obtemos diretamente

```
||T H(s) S|| = ||T B Q0 S|| <= 23 t_P,
||S H(s) T|| <= (23+rho) t_G,
||T(H(s)-I)T|| <= (23+rho) t_G.
```

Na primeira linha, tanto `TS` como `T Q0 S` são zero; por isso o termo
rho desaparece. Nas demais ele permanece. São bounds de operadores
infinitos, condicionais às estimativas analíticas de ALGEBRAIC_TAIL,
não extrapolações da matriz finita. Valem também quando cada acoplamento
e o defeito exterior são multiplicados por parâmetros em [0,1].

Com P=12×36, z0=0 e rho=5/4 (contorno A atual), os valores exatos são:

| G | S→T | T→S | defeito T→T |
|---|---:|---:|---:|
| 24×72 | 2484/35 | 2619/71 | 2619/71 |
| 64×192 | 2484/35 | 2619/191 | 2619/191 |
| 874×2621 | 2484/35 | 2619/2620 | 2619/2620 |

A última caixa faz apenas o defeito exterior ficar abaixo de um;
NÃO fecha o critério acoplado. As duas primeiras falham já na diagonal.
Assim, agora há majorantes para todos os caminhos, mas ainda não controle
contrativo. Reponderar blocos não conserta um bound diagonal maior que um.
O produtor `crown_exterior` registra essa distinção explicitamente.

Para usar uma inversa acoplada V de G, não basta multiplicar o bound
núcleo←T por ||V|| e esquecer a coroa: para i=P,S,
`||E_iT|| <= ||V_iP|| ||PHT|| + ||V_iS|| ||SHT||`.
Além disso `||E_TP||=||THP||`, `||E_TS||=||THS||`.
Essas expressões fornecem a matriz majorante de três blocos completa.
V deve ser certificado por segmento; uma inversa validada somente em s=1
não fornece uniformidade no contorno.

## Reprodução

```bash
OPENBLAS_NUM_THREADS=2 .venv/bin/python scripts/check_frequency_blocks.py \
  build/spectrum/rt-A-6x18.dat --core-m 4 --core-n 12 \
  --out build/spectrum/frequency-blocks-6x18-core4x12.json
.venv/bin/python scripts/background_couplings.py \
  .cache/rt-1203.3766v1/sourcecode/RefA.dat \
  --out build/spectrum/background-couplings.json
```
