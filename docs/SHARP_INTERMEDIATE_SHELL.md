# Faixa sharp incorporada e redução do restante a um Schur finito

Os resultados desta nota valem no disco |s-0.733|<=1/32, nos pesos
(129/128,9/8), condicionados à bola RT. Não são uma exclusão do pencil
sharp infinito completo.

## Primeira faixa incorporada sem descartar os acoplamentos

Foi auditada a matriz sharp 12×36 inteira, dimensão 594. Sua compressão
para os 252 DOFs de 8×24 coincide EXATAMENTE, em diádicos, com a matriz
anterior. Os 342 novos DOFs e ambos os acoplamentos entram na inversa
completa; não foram invertidos separadamente e depois justapostos.

Incluindo a incerteza de mu e do fundo, o defeito uniforme é
<0.265850 e a inversa do pencil não precondicionado é <11.588.
Relatório: sharp-disk-12x36-background.json.

## Acoplamento validado com uma cauda realmente infinita

Ponha F=12×36, T=I-G e H=K_true(s) Q_ref, com Q_ref sharp FIXO.
Como Q_ref preserva F, a inversa finita precondicionada satisfaz

```
H_FF^-1=J_ref,F K_true,FF^-1,
||H_FF^-1|| <= kappa <1201.510.
```

Essa constante inclui a norma racional de J_ref,F: não reutilizamos
11.588 como se fosse a inversa do pencil precondicionado.

O principal, a identidade e Gamma1sharp preservam F. Portanto T H F
vem somente de Gamma2sharp do fundo. A filtragem de suporte por
G-F+1 seguida do majorante por componentes dá

```
||H_TF|| <= c <0.000018240,
||H_FT|| <= d,
||H_TT-I|| <= d,
d=(4.62237+4/5)t_sharp(G)+210 epsilon_mu.
```

Para c foi usado ||Q_ref||<=204, válido nos pesos menores, e o erro
3 epsilon_omega. RefA tem suporte contido em 41×101: para as caixas
abaixo, o termo de referência filtrado é zero, mas o erro do fundo NÃO.

O Schur em T, R_T=H_TT-H_TF H_FF^-1 H_FT, tem defeito
`delta=d+c*kappa*d`. As verificações racionais dão:

| G | bound de delta | DOFs restantes em U=G-F |
|---|---:|---:|
| 256×1536 | <0.585 | 586926 |
| 160×960 | <0.936 | 228366 |

Assim H comprimido a F+T é invertível em todo o DISCO, não só em sua
fronteira. Todos os blocos são holomorfos ali. Os relatórios são
sharp-guarded-bridge.json e sharp-guarded-bridge-160x960.json.

## O que permanece, sem omitir a faixa intermediária

Se Z=F+T e U=G-F, o problema sharp infinito no disco equivale ao
problema finito, ainda NÃO auditado,

```
E_U(s)=H_UU-H_UZ (H_ZZ)^-1 H_ZU.
```

A equivalência segue da fatorização triangular por operadores limitados,
holomorfos e invertíveis. Portanto kernels e multiplicidades analíticas
são preservados. Como H_ZZ é invertível em TODO o disco, ele não fornece
raízes ocultas ali. Não afirmamos isso para a antiga coroa selecionada,
cuja cobertura era apenas no contorno.

Para construir E_U, não basta usar a matriz crua H_UU. A inversa em Z é

```
[[A + A H_FT R_T^-1 H_TF A, -A H_FT R_T^-1],
 [-R_T^-1 H_TF A,           R_T^-1]],       A=H_FF^-1.
```

A parte infinita R_T^-1 tem aproximação quantitativa de Neumann. Se
R_T=I+D, ||D||<=delta, truncar sum_0^N(-D)^k deixa erro
<=delta^(N+1)/(1-delta). Para a caixa menor, delta<=117/125 e N=512
dão erro <2.873e-14; para a maior, delta<=117/200 e N=64 dão <1.767e-15.
Esses são erros da inversa R_T: ainda devem ser multiplicados pelos
fatores dos acoplamentos para controlar E_U. Nenhuma dessas potências
ou a matriz E_U foi construída nesta rodada.

O gargalo ficou explícito: validar a parametriz estruturada desse Schur
finito de 228366 DOFs, incluindo as correções, ou reduzir o corte com
bounds mais fortes. A matriz densa dessa dimensão não é a estratégia
computacional adotada. O bound de cauda sozinho não fecha essa etapa.

## Reprodução

```bash
OPENBLAS_NUM_THREADS=2 .venv/bin/python scripts/certify_sharp_disk.py \
  build/spectrum/rt-sharp-12x36.dat --include-background \
  --inner-matrix build/spectrum/rt-sharp-8x24.dat \
  --out build/spectrum/sharp-disk-12x36-background.json
.venv/bin/python scripts/sharp_guarded_bridge.py --outer-m 160 --outer-n 960 \
  --out build/spectrum/sharp-guarded-bridge-160x960.json
```

Os hashes vinculam os relatórios às matrizes e a RefA. O checker escalar
não refaz os resíduos matriciais do produtor anterior. Os resultados de
operadores continuam condicionais aos lemas analíticos RT e não foram
formalizados em Lean.

## Continuação: o erro exterior já propagado até a faixa

O relatório `sharp-guarded-bridge-148x840.json` incorpora as auditorias
mais recentes: kappa<11.461 para a inversa FINITA precondicionada,
||Q_ref||<5.942 e delta<0.989223. A faixa restante tem 184626 DOFs.
Essas constantes substituem os bounds anteriores nessa caixa, sem
alterar as hipóteses analíticas nem eliminar U.

É possível calcular a tolerância que uma parametriz de U precisará
consumir, incluindo ambos os acoplamentos. Com A=H_FF^-1, escreva

```
B_U = H_UU-H_UF A H_FU,
L_U = H_UT-H_UF A H_FT,
R_U = H_TU-H_TF A H_FU,
E_U = B_U-L_U R_T^-1 R_U.
```

Esta forma é a mesma eliminação em Z=F+T, com a ordem de todos os
fatores preservada. Se R_T=I+D e P_N=sum_(j=0)^N(-D)^j, então

```
E_(U,N)=B_U-L_U P_N R_U,
||E_U-E_(U,N)|| <= ell*r*delta^(N+1)/(1-delta),
ell >= ||H_UT||+||H_UF||*kappa*||H_FT||,
r   >= ||H_TU||+||H_TF||*kappa*||H_FU||.
```

Em particular, o erro de R_T^-1 sozinho não é a tolerância do Schur.
Para explicitar bounds válidos sem omitir nenhuma banda, ponha

```
d_F=(4.62237+4/5)t_sharp(12,36)+210 epsilon_mu,
e=4.62237 ||Q_ref||.
||H_FU||, ||H_TU|| <= d_F,
||H_UF|| <= e,
||H_UT||, ||H_FT|| <= d,
||H_TF|| <= c.
```

As duas primeiras colunas usam U fora de F. Na coluna F, o principal,
sQ_ref e a parte descendente preservam F, e beta_sharp também é um
majorante conservador da multiplicação. Todos os termos do erro do
fundo usados pela ponte permanecem incluídos. Portanto podemos tomar
ell=d(1+e*kappa) e r=d_F(1+c*kappa).

`sharp_shell_remainder.py` calcula e arredonda esses bounds para fora
antes de elevar delta a potências. Na caixa 148×840, obtém

```
ell <= 312.334548, r <= 22.42237, delta <= 0.989223,
N=3785 (3786 termos),
||E_U-E_(U,N)|| <= 9.92351e-13 < 10^-12.
```

O grau é o menor que passa para esses majorantes arredondados. A busca
e a desigualdade final usam racionais exatos. O relatório é
`build/spectrum/sharp-shell-remainder-148x840.json`, ligado ao relatório
de entrada por SHA-256. Não recalcula seus resíduos matriciais.

Reexecutando a ponte 160×960 com os MESMOS relatórios melhorados da
inversa finita e de Q_ref, o defeito cai para <0.915031. O relatório
`sharp-guarded-bridge-160x960-preconditioned.json` mantém 228366 DOFs
intermediários; `sharp-shell-remainder-160x960-preconditioned.json`
exige N=437 (438 termos) para erro propagado <=9.7964e-13. Isso explicita
o custo de aproximar a inversa exterior: reduzir a dimensão da faixa
até 148×840 piora muito a convergência da série com esses majorantes.

Isso fornece um orçamento de erro, NÃO uma parametriz de U. P_N e
E_(U,N) ainda não foram construídos. Embora a compressão final seja
finita, seus fatores usam operadores do fundo exato: implementá-los a
partir de RefA exige também orçar os erros dessa implementação. Os
3786 termos não são anunciados como um cálculo já realizado.

O critério suficiente restante é explícito. Se uma matriz candidata V
satisfaz, uniformemente no disco,

```
||I-V Ehat_(U,N)|| <= eta, ||V|| <= v,
||Ehat_(U,N)-E_(U,N)|| <= epsilon_impl,
eta+v*(epsilon_impl+ell*r*delta^(N+1)/(1-delta)) < 1,
```

então o Schur é invertível e a ponte dá exclusão sharp INFINITA no
disco. Nenhum valor validado de eta, v e epsilon_impl foi produzido.
`close_shell` apenas verifica esse critério escalar quando esses dados
são fornecidos; não os inventa. Os testes com eliminação racional de
matrizes de três blocos verificam a fórmula e a propagação do erro.

```bash
.venv/bin/python scripts/sharp_shell_remainder.py \
  build/spectrum/sharp-guarded-bridge-148x840.json \
  --out build/spectrum/sharp-shell-remainder-148x840.json
```
