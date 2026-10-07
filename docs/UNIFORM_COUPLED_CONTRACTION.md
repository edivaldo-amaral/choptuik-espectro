# Critério acoplado uniforme, com homotopia explícita

Fixe uma única caixa finita F, T=I-F e uma única norma de componentes.
Para H(s)=I+(B+(s-z0)I)Q0, escolha uma parametriz finita V_j por tile.
Ela fica fixa dentro daquele tile. Suponha, em TODO o tile,

```
a >= ||I-V_j H_FF(s)||,
b >= ||V_j H_FT(s)||,
c >= ||H_TF(s)||,
d >= ||H_TT(s)-I||.
```

Se d<1, elimine o exterior: o complemento de Schur é
S_F=H_FF-H_FT H_TT^-1 H_TF. Temos

```
||I-V_j S_F|| <= a + b*c/(1-d).
```

Assim basta `d<1` e `a+b*c/(1-d)<1`. Isso equivale ao critério
ponderado de dois blocos já implementado. Não é suficiente d<1 sozinho.
O bound também prova que V_j é invertível, pois a<1 implica que
V_j H_FF é invertível e ambas as matrizes finitas são quadradas.

## Variação dentro de um tile

Com centro s_j e |s-s_j|<=r_j, um bound suficiente é

```
a_j = a_center + r_j ||V_j F Q0 F|| + error_center.
```

error_center inclui as diferenças entre o fundo exato e o exportado,
arredondamento e, quando aplicável, a mudança de mu no principal e em Q0.
O erro de mu NÃO é apenas uma perturbação multiplicativa pelo fundo:
ele multiplica a derivada espacial. Não basta inserir somente
6 epsilon_background em seu lugar. A avaliação posterior, com Q_ref
fixo, está em DAMPED_CROWN_PARAMETRIX.md: erro total <0.0003257,
condicionado à bola RT. Isso não avalia os acoplamentos infinitos.

Se b,c,d são uniformes nesse tile, o raio deve satisfazer estritamente

```
r_j < [1-a_center-error_center-b*c/(1-d)] / ||V_j F Q0 F||.
```

Quando o numerador é não positivo, subdividir o tile não resolve esse
majorante. O script uniform_coupling.py confere exatamente a margem e
a estrita desigualdade; não infere a proveniência dos bounds.

## Homotopia e cobertura

Use, com 0<=t<=1,

```
H_t(s) = [[H_FF(s), t H_FT(s)],
          [t H_TF(s), I+t(H_TT(s)-I)]].
```

O defeito de Schur é majorado por `a+t²bc/(1-td)`, no máximo
`a+bc/(1-d)`. Logo o mesmo certificado serve para toda a homotopia.
Em t=0, o operador é diag(H_FF,I). Para transportar a contagem, ainda
se exige a estrutura Fredholm/compacta apropriada e que a união dos tiles
cubra o contorno sem lacunas; invertibilidade em um centro não é cobertura.

## Aplicação honesta aos dados atuais

Na norma eta, o novo bound d em 64×192 é aproximadamente 4.882.
Ele falha mesmo se a=b=c=0. Em 175×942, d<1, mas sua margem é apenas
101/117875: o termo bc é amplificado por 117875/101, cerca de 1167.
Além disso, não existe uma parametriz certificada para essa caixa.
A inversa finita disponível é de 6×18: não se pode combiná-la com o
bound exterior de 175×942 omitindo a coroa entre as duas caixas.

Portanto a contração uniforme completa ainda não foi obtida. O trabalho
restante está localizado: uma parametriz da coroa/exterior que absorva
parte de B, ou bounds de acoplamento suficientemente fortes, mais tiles
finitos e a incerteza completa do centro. Nenhum novo winding físico
decorre somente dos bounds desta nota.
