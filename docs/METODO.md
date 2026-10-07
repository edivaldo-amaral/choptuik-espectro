# O que funcionou e o que não funcionou, como método

Nota de método, 21/09/2026. Escrita depois de três dias em que sete alavancas
foram fechadas com medida e quatro peças atravessaram para o lado certificado.
O objetivo é que a próxima retomada não repita o padrão que falhou.

## O padrão que falhou

Sete tentativas, todas entre 1,0× e 2,7×:

| alavanca | ganho |
|---|---:|
| reponderação por componentes | 1,0009× |
| família produto `kappa1^m kappa2^n` | 1,01× |
| ponto-base `z0` | pior |
| borda do contorno `sigma0` | 1,22× |
| absorver `B` na parametriz | 1,24× |
| razão de aspecto | ~1,5× |
| separação aditiva núcleo/caixa | ~1,7× |
| deflação analítica dos modos não físicos | 2,66× |

**O que elas têm em comum:** todas tentam melhorar uma *constante* dentro de uma
formulação fixa. Escolhem um parâmetro melhor — um peso, um ponto-base, uma
borda, um corte — sem mudar o que está sendo calculado.

Há um teto natural para isso, e ele apareceu: a formulação já estava perto do
seu ótimo, porque quem a montou também tinha otimizado parâmetros. Os pesos de
RT, por exemplo, estão a 0,09% do piso de Perron e a 1% do ótimo da família
produto — escolhidos para a prova de existência, e ainda assim quase ótimos
para o problema espectral.

## O padrão que funcionou

Quatro movimentos, todos com ganho qualitativo:

- **Trocar a representação da mesma quantidade.** `lambda_min` com peso `sqrt(w)`
  tem irracionais; o *mesmo* número é o menor autovalor do pencil generalizado
  `Herm(WJ) y = lambda W y`, e ali tudo é racional. Certificado exato, de graça.
- **Checar se o obstáculo é real antes de contorná-lo.** A ressalva `l^1` vs
  `l^2` parecia bloquear o campo de valores. Uma linha —
  `sum (w_j|x_j|)^2 <= (sum w_j|x_j|)^2` — mostra que a inclusão é gratuita na
  direção necessária.
- **Perguntar onde o extremo vive, não quanto ele vale.** O autovetor de
  `lambda_min(Herm(J_mu))` está inteiramente na componente 3, com perfil radial
  fixo. Isso reduziu um problema de dimensão 1422 a um 1-D de dimensão efetiva
  18.
- **Olhar a forma da matriz, não as suas normas.** `Herm(W J_mu)` é
  semisseparável de posto 1 mais diagonal. Isso colapsa o `LDL^T` numa recursão
  escalar de duas linhas, cujo ponto fixo é `1/kappa2` — e a constante `0,2043`
  que eu tinha ajustado numericamente era `1/5`.

**O que eles têm em comum:** todos mudam *o que* está sendo calculado, não
*quão bem*. E todos vieram da mesma pergunta — "o que é este objeto,
exatamente?" — e não de "como limito isto melhor?".

## A regra prática

> Quando uma sequência de tentativas rende ganhos de 1 a 3×, isso não é falta de
> engenho: é o sinal de que a formulação está perto do seu ótimo e o ganho tem
> de vir de **mudar a formulação**, não de ajustá-la.

E o corolário, que custou caro aprender: um levantamento de alternativas precisa
ser auditado quanto à completude. O mapa de becos sem saída de 21/09 tinha
quatro entradas; varrendo o projeto por famílias de técnica padrão, **sete
famílias não apareciam em lugar nenhum** — e a primeira que testei (campo de
valores) rendeu um fator 40 numa obrigação aberta.

## Os erros, e o que os pegou

Cinco resultados meus foram retratados em três dias. Vale registrar o que os
detectou, porque foi sempre a mesma coisa:

| erro | o que o pegou |
|---|---|
| "1,29× de folga em contabilidade" | comparar com o piso de Perron |
| caixa de 14,4× | conferir qual restrição mordia |
| `||V||` monótono na distância à raiz | medir em vez de supor |
| grade diagonal que "funcionava" | o bound saiu **menor** que o valor real |
| peso balanceado de Perron | olhar a *forma* do peso, não só o resultado |

Nenhum foi pego por revisão do raciocínio. Todos foram pegos por **uma
verificação numérica barata que o resultado deveria ter passado e não passou**.
A prática que ficou: depois de um resultado favorável, procurar a medida mais
barata que o refutaria.

---

## Acréscimo de 23/09: reutilizar é uma forma de verificar

A cadeia de F3 (`R <= 6,05`) tinha dois passos injustificados, e nenhum deles foi
achado relendo a derivação. Foram achados ao **reutilizar** a máquina para outro
operador (`K`, o sistema de constraints): para adaptá-la foi preciso perguntar em
que produto interno cada peça vive, e as respostas não batiam.

1. **Convenções misturadas.** `R_J` e o `Resto` estavam em `D = diag(sqrt(w))`;
   a parte linear 1,511 estava em `w` inteiro. Numa convenção a divisão por ξ vale
   40/9; na outra, 8,94.
2. **Estrutura de ℓ¹ levada a ℓ².** O máximo sobre sete combinações é a forma que
   a norma ℓ¹ produz, porque nela a entrada é a *soma* das componentes. Em ℓ² é
   outra estrutura.

A lição geral: **ao trocar ℓ¹ por ℓ², não basta trocar cada norma.** Toda
estrutura que dependia da soma sobre componentes precisa ser refeita.

E a lição de método: um resultado que nunca foi reutilizado nunca foi testado
nas suas convenções. Reutilizar a máquina em outro problema é uma verificação
barata, e deveria ser feita de propósito.
