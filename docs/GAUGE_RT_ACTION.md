# Ação RT da translação nula e completude residual local

**Atualização 03/10/2026:** a completude GLOBAL do gauge residual, para campos de gauge C¹ em todo o
domínio de RT e sem analiticidade dos germes, está provada em `GAUGE_GLOBAL.md` (Lema G).

Continuação: GEOMETRIC_FLOQUET_RECONSTRUCTION.md estabelece a reconstrução
global e recíproca DENTRO do ansatz para Re s>0, incluindo jets. Ainda
não prova que toda perturbação física possa ser posta nesse ansatz por
um gauge global admissível; a condição de cone permanece explícita.

Convenção: variação ativa por Lie, com mu fixo. Trocar para a convenção
passiva muda o sinal do gerador, não seu expoente. As definições geométricas
usadas são as equações `coordinatesxitau`, `ekhfkjfhkfhfkhf` e `ekjehee`
do artigo RT, disponíveis no arquivo-fonte local `choptuik.tex`.

## Gerador que move o cone

Para delta u_-=delta u_+=1, escreva

```
X = exp(mu tau) Z,       Z = mu^-1 partial_tau + xi partial_xi.
W = mu + xi² Q,          F = zeta + log W.
```

O bloco radial da métrica, em coordenadas nulas, tem fator proporcional a
`exp(-2 zeta)/(u_++u_-)²`. Sua derivada de Lie dá
`delta zeta = exp(mu tau)(Z zeta-1)`. O raio de área
`R=exp(-zeta)xi/W` é um escalar. Comparando delta log R com X log R,
obtemos `delta W=exp(mu tau) ZW`, ou

```
delta Q    = exp(mu tau)(ZQ+2Q),
delta F    = exp(mu tau)(ZF-1),
delta phi  = exp(mu tau) Zphi.
```

Usando `D^sigma=sigma partial_tau+mu(1+sigma xi)partial_xi`,
vale `[D^sigma,X]=exp(mu tau)D^sigma`. Substituição nas definições
de omega2, omega3 e omega4 dá o mesmo perfil para todas as componentes:

```
h_i = (Z+1)omega_i,       s=mu.
```

Para omega1, usa-se `omega1=2xi Q+2mu partial_xi F`; os termos extras
também somam exatamente omega1. Não se trata de identificar uma raiz
por proximidade numérica: é a ação induzida pelo difeomorfismo.

## Identidade fora das soluções

Para qualquer fundo regular w, com primeira componente ímpar,

```
D Omega(w)[exp(mu tau)(Z+1)w]
    = exp(mu tau)(Z+1)Omega(w).
```

Isso segue pelo comutador acima, pela homogeneidade quadrática de Gamma2
e pelo fator xi nos termos diferenciais e quadráticos. Em particular,
depois de aplicar a seleção que divide por xi,

```
L_w(mu)(Z+1)w = (Z+2) F_w,
C_w(mu)(Z+1)w = (Z+1) c_w.
```

`test_gauge_action.py` confere as cinco componentes da identidade completa
e a seleção/constraints com polinômios racionais, usando as tabelas RT.
Testa também a covariância por fase em s=0. Para o fundo EXATO, F_w=c_w=0,
logo h é solução das equações completas, com constraints nulas.
Para RefA, os lados direitos NÃO podem ser descartados.

## Não nulidade e domínio

Atualização: RADIUS_RECOVERY.md resolve a perda de raio para este vetor
usando a equação L(mu)h=0. O argumento de derivação direta abaixo só o
coloca inicialmente nos pesos menores; a recuperação o leva ao domínio
RT forte sem supor derivadas limitadas nesse mesmo espaço.

Se `(Z+1)omega4=0`, cada coeficiente Fourier de omega4, desenvolvido em
Taylor em xi, satisfaz `(1+n+im/(2mu))a_mn=0`. O fator tem parte real
positiva; portanto omega4 seria identicamente zero, excluído no fundo RT.
Assim o gerador é não nulo. Preserva as paridades e periodicidades RT.

As derivadas em Z são limitadas com perda de raio conforme SHARP_DOMAIN.
Logo o perfil pertence aos pesos estritamente menores. A equação L(mu)h=0
o coloca no domínio máximo nesse espaço e, pela igualdade dos domínios,
no fechamento mínimo. A derivação direta NÃO prova pertinência aos pesos
fortes; o argumento posterior de recuperação de raio completa essa etapa.
Isso ainda não permite subtrair uma raiz de uma contagem finita antiga.

## Completude local, dentro do ansatz fixado

Considere vetores esfericamente simétricos, analíticos, que preservam a
forma du_+ du_- da métrica radial, com mu fixo. Os componentes du_+² e
du_-² de sua derivada de Lie precisam se anular. Como g_+- é não nulo,
isso força `partial_+ X^-=0` e `partial_- X^+=0`, portanto
`X^-=f_-(u_-)`, `X^+=f_+(u_+)`. Preservar o centro força f_-=f_+.
Assim a classe de GAUGE_DOMAIN_CLASSIFICATION é localmente completa
para ESTE gauge residual do ansatz, com a regularidade especificada.
Sua classificação de Floquet implica que o único expoente real positivo
é mu, com gerador acima, se for permitido mover o cone.

Fixar o cone exige f(0)=0 e exclui esse gerador, pois X(u_-)=1 no cone.
Nesse caso não há gauge positivo nesta classe. Não confundir essa
condição geométrica com simplesmente exigir um perfil RT analítico.

## O que falta para equivalência física completa

Faltam a fixação do gauge de uma perturbação geométrica arbitrária nesse
ansatz, com estimativas e condições globais, e a escolha consistente de
cone fixo/móvel na realização espectral. A reconstrução também precisa
fixar constantes invisíveis aos omega: delta zeta constante e delta phi
constante não alteram esses campos. A primeira é uma escala conforme,
não automaticamente um difeomorfismo; a segunda é removida pela condição
antiperiódica adotada para phi. As normalizações devem ser mantidas na
equivalência linearizada. Finalmente, a classificação de autovetores
não prova ausência de cadeias generalizadas gauge. Nenhuma multiplicidade
algébrica foi subtraída nesta etapa. CONSTRAINT_ROOT_QUOTIENT.md esclarece
a subtração condicional da linha gauge e o tratamento de cadeias mistas;
a equivalência geométrica global continua necessária.
