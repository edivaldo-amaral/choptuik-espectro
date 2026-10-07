# Gauge residual, cone fixo e seção de fase

**Atualização 03/10/2026:** a completude GLOBAL do gauge residual, para campos de gauge C¹ em todo o
domínio de RT e sem analiticidade dos germes, está provada em `GAUGE_GLOBAL.md` (Lema G).

Continuação: GAUGE_RT_ACTION.md deriva a ação do gerador n=0 nos campos
RT e prova completude LOCAL do gauge residual dentro do ansatz com mu
fixo. As limitações históricas abaixo sobre a ação foram resolvidas ali;
a equivalência física global e as multiplicidades continuam abertas.

## Classe de transformações estudada

Esta nota classifica uma classe PRECISA, não todos os difeomorfismos físicos:
reparametrizações infinitesimais das coordenadas nulas de RT
`delta u_-=f_-(u_-)`, `delta u_+=f_+(u_+)`, preservando o centro u_-=u_+,
com germes analíticos em u=0. Preservar o centro força f_-=f_+=f no
intervalo comum. A condição adicional de preservar o cone u_-=0 é f(0)=0.
Resta provar que essa classe representa todo o gauge do problema físico
após fixar o ansatz; essa completude não é presumida aqui.

Use `u_sigma=-(1+sigma*xi) exp(-mu*tau)`, sigma=±1. Resolvendo as duas
equações para o vetor em coordenadas (tau,xi), obtemos

```
delta tau = exp(mu*tau) [f(u_+)+f(u_-)]/(2 mu),
delta xi  = exp(mu*tau) [(1+xi)f(u_-)-(1-xi)f(u_+)]/2.
```

Para f(u)=u^n, n inteiro não negativo, esse vetor tem a forma
`exp(s_n*tau) (A_n(xi),B_n(xi))`, com s_n=mu(1-n), e

```
A_n = (-1)^n [(1+xi)^n+(1-xi)^n]/(2 mu),
B_n = (-1)^n [(1+xi)(1-xi)^n-(1-xi)(1+xi)^n]/2.
```

Os perfis são polinomiais. `gauge_null_profiles.py` e os testes verificam
exatamente o pushforward para ambas as coordenadas nulas.

## Classificação de um modo de Floquet

Se o vetor tem um único multiplicador de Floquet exp(sT), então
`f(exp(-mu*T)u)=exp((s-mu)T)f(u)`. Desenvolva o germe
f(u)=sum a_n u^n. Para qualquer coeficiente não nulo,
`exp(-n*mu*T)=exp((s-mu)T)`, logo Re s=mu(1-n).
Como mu>0, dois graus diferentes não podem compartilhar esse multiplicador.
Isso classifica os modos, módulo os aliases imaginários de Floquet.

- n=0: Re s=mu>0, vetor exp(mu*tau)(1/mu,xi). Ele move o cone:
  delta u_-|_(u_-=0)=1. É uma translação conjunta das coordenadas nulas;
  não se identificou aqui com a translação do tempo próprio físico.
- n=1: s=0 módulo aliases, vetor (-1/mu,0), a translação da fase tau.
- n>=2: Re s<0, e o cone é preservado.

Portanto, nesta classe analítica, **não há gauge instável que preserve o
cone fixo**. Se o domínio físico permite mover o cone, o modo n=0 precisa
ser incluído numa formulação com domínio móvel e depois quocientado.
Não se pode alternar entre essas duas realizações dentro da mesma contagem.
A proximidade de uma raiz numérica a mu não decide qual realização ela
representa, nem fornece sua ação nos campos RT.

## Seção de fase: transversalidade quantitativa

O código RT fixa `b(w)=Re w4_(m=1,n=0)=0` (`Field_in_VGauged`). Para
p=partial_tau w*, a convenção Fourier im/2 dá
`b(p)=-Im w4_(1,0)/2`. Em RefA,

```
Im RefA4_(1,0) = 138672388040959354547 / 590295810358705651712.
```

O peso desse coeficiente é 65/32. Pela bola publicada de raio
epsilon=2^-25+2^-277, segue o bound racional

```
|b(p)| >= (|Im RefA4_(1,0)| - epsilon/(65/32))/2 > 117/1000.
```

Assim a seção é transversal à fase. No espaço analítico onde p pertence
ao domínio (por exemplo com a perda de raio de SHARP_DOMAIN), a projeção
`Pi h=h-p*b(h)/b(p)` identifica o quociente pelo gauge de fase com ker b,
no kernel NEUTRO, pois L(0)p=0. Isso verifica a hipótese BG invertível do
lema de seção para esse subgrupo unidimensional, não para o gauge inteiro.

## Por que não aplicar essa projeção aos modos positivos

Se L(s)=L(0)+sI e L(0)p=0, então L(s)p=s p. Para h em ker L(s),

```
L(s) Pi h = -s p*b(h)/b(p).
```

Logo, para s diferente de zero, a projeção de fase em geral NÃO preserva
o kernel. Cortar uma coordenada de fase no pencil positivo pode mover o
autovalor físico; não é uma operação de quociente justificada nesse s.
A eliminação de gauge em outros expoentes exige um gerador G(s) que seja
ele próprio uma solução das equações completas no mesmo domínio.

## Estado

Estão demonstradas a classificação desta classe de coordenadas nulas e a
transversalidade da fase, condicionada ao erro publicado do fundo.
GAUGE_RT_ACTION.md acrescenta a ação induzida do gerador positivo em omega
e nas variáveis geométricas, e a completude residual local dentro do ansatz.
Continuam abertas a fixação global do gauge, a equivalência do domínio
físico com a realização espectral e o controle das cadeias algébricas.
Não foi subtraída nenhuma raiz positiva.
