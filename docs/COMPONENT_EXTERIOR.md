# Exterior com estrutura entre componentes

## Parte descendente exata

Escreva `R_x f=(f-f(0))/xi` e P para a reflexão espacial. Nas quatro
componentes selecionadas, a inspeção de Gamma1 fornece EXATAMENTE

```
S Gamma1 h = (0,
              R_x h1 + R_x(1-P)h2,
             -R_x h1,
              R_x(1-P)h4).
```

A primeira componente é zero porque h1 é ímpar, R_x h1 é par e a seleção
nessa componente retém a parte ímpar. Não é necessário estimá-la por um
número positivo. Todos esses termos baixam o grau radial e preservam m.
O bound `||R_x||<=40/9`, com mu<=17/100, dá a matriz majorante

```
A_linear = [[0,     0,     0, 0],
            [34/45, 68/45, 0, 0],
            [34/45, 0,     0, 0],
            [0,     0,     0, 68/45]].
```

Ela pode ser absorvida em um principal triangular por componentes:
diag(K_mu, J_mu+mu R_x(1-P), J_mu, J_mu+mu R_x(1-P)), com fontes
mu R_x h1 nas componentes 2 e 3. A inversa desse novo principal ainda
não foi certificada; aqui usamos a matriz acima, sem presumir seu ganho.

## Multiplicação pelo fundo com cancelamentos preservados

O programa `component_exterior.py` lê os quatro campos de RefA sem perder
os sinais dos coeficientes e interpreta as entradas de `GAMMA2_macro.c`.
Para cada componente de entrada e sua paridade radial, combina primeiro
os campos do fundo conforme o macro e SOMENTE depois toma valores absolutos.
Assim, por exemplo, diferenças entre campos correlacionados não são
substituídas imediatamente pela soma das duas normas.

O operador é `2 S Xi Gamma2(RefA,h)`. Na entrada h1 ímpar, o fator do macro
é 1/2; nas outras entradas, os fatores 2 vindos de h_j±Ph_j deixam fator 1.
Na primeira saída, a projeção ímpar seleciona a parte do coeficiente cuja
paridade é complementar à da entrada. Para cada bloco i,j, toma-se o máximo
entre as duas paridades de entrada permitidas. A submultiplicatividade da
norma Fourier–Chebyshev limita a convolução pela norma do coeficiente.

Isso dá uma matriz não negativa A_ref tal que
`||(B_ref h)_i|| <= sum_j (A_ref)_ij ||h_j||`.
Os racionais exatos estão em `build/spectrum/component-exterior.json`.
A matriz da parte não linear, apenas para visualização aproximada, é

```
[[1.029826, 5.516999, 2.753932, 6.494014],
 [1.217055, 5.507865, 0,        0       ],
 [1.217055, 2.753932, 0,        6.494014],
 [0,        6.494014, 0,        2.753932]].
```

## Norma adaptada e incerteza do fundo

Use `||h||_eta=sum_i eta_i ||h_i||`, com pesos exatos proporcionais a
`(916,4096,243,1233)`. O autovetor numérico apenas propôs esses pesos;
o bound final foi recalculado com racionais:

```
beta_ref = max_j sum_i eta_i (A_ref+A_linear)_ij / eta_j.
beta_true <= beta_ref + 6 epsilon_background max(eta)/min(eta)
          < 5191/500 = 10.382.
```

Aqui `epsilon_background=2^-25+2^-277`, condicionado à bola publicada RT.
O erro do fundo NÃO foi estimado na nova norma sem conversão: o fator
4096/243 está explicitamente incluído. Os hashes do fundo e do macro são
registrados no relatório; hashes identificam os dados, não provam a bola RT.

Q0 é diagonal nas quatro componentes, de modo que seu bound de
STRUCTURED_EXTERIOR continua válido na nova norma sem perda adicional.
O termo delta I também mantém norma |delta|. Logo, uniformemente em
|s|<=5/4, para H(s)=I+(B+sI)Q0,

```
||T(H(s)-I)T||_eta <= (5191/500+5/4) t_structured(G).
```

Em G=64×192 isso é menor que 4.882: o majorante original era 13.712.
Em G=175×942 o bound é `117774/117875 <1`.
Isso certifica SOMENTE a desigualdade numérica do bound diagonal exterior,
sob as hipóteses analíticas; não certifica a matriz completa nem um winding.

Na mesma norma, também valem
`||THS||_eta <= (5191/500) t_structured(P)` e
`||SHT||_eta <= (5191/500+5/4) t_structured(G)`.
O termo s desaparece da primeira desigualdade porque TQ0S=0.

## Não misturar normas ao fechar os blocos

Os defeitos e inversas finitos anteriores foram calculados com eta_i=1.
Para usar o novo beta é necessário recalculá-los conjugados pelos pesos
eta, ou majorá-los pela razão 4096/243. O mesmo vale para bounds de
acoplamento baseados na cauda do fundo. Usar beta=10.382 junto com os
blocos antigos sem conversão seria uma lacuna. A parametriz exterior
estruturada, a validação uniforme da parte finita, constraints e gauge
continuam pendentes.

Foi conferido o transporte da parametriz acoplada existente de 6×18:
com W a diagonal dos novos pesos, tome V_eta=W V W^-1 e
H_eta=W H W^-1. Então I-V_eta H_eta=W(I-VH)W^-1. Os quatro bounds
anteriores podem ser multiplicados por 4096/243. O novo majorante sem
reponderação adicional núcleo/coroa é estritamente menor que 2.28e-8.
Portanto essa parametriz FINITA permanece validada na nova norma no
centro s=1. `transport_defect` implementa essa passagem racional.

O recálculo direto de outra parametriz diádica com `--component-weights`
foi interrompido depois de validar esse transporte; não há resultado
concluído daquele recálculo. A opção permanece disponível para tentar
bounds melhores. O transporte usa o relatório existente
`frequency-blocks-6x18-core4x12.json`, terceira parametriz, sem alterar o fundo.

Reprodução:

```bash
.venv/bin/python scripts/component_exterior.py \
  .cache/rt-1203.3766v1/sourcecode/RefA.dat \
  --out build/spectrum/component-exterior.json
OPENBLAS_NUM_THREADS=2 .venv/bin/python scripts/check_frequency_blocks.py \
  build/spectrum/rt-A-6x18.dat --core-m 4 --core-n 12 \
  --component-weights 916,4096,243,1233 \
  --out build/spectrum/frequency-blocks-6x18-component-weights.json
```

Controle independente finito: subtraindo o J_mu construído analiticamente
da matriz A_6x18 exportada, a norma ponderada exata de B_6x18 tem valor
decimal aproximado 8.077805 (arredondado para cima), abaixo de 10.382.
É um teste de consistência; não foi usado para limitar o exterior infinito.
