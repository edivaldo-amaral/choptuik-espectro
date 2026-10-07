# Transporte global do gauge nulo com condição no cone

Esta nota resolve o transporte necessário à remoção das componentes
radiais diagonais para uma classe analítica explícita, com perda de
vizinhança radial. Não identifica essa classe com toda a realização
espectral RT. O resultado é analítico em papel; os testes verificam
identidades racionais, não substituem a prova funcional abaixo.

## Da métrica ao transporte periódico

Use u_-=-(1-xi)exp(-mu tau), u_+=-(1+xi)exp(-mu tau), mu>0,
e g_+- não nulo. Para uma perturbação esfericamente simétrica k, ponha

```
X^- = exp((s-mu)tau) y_-(tau,xi),
X^+ = exp((s-mu)tau) y_+(tau,xi),
f_- = -mu exp(-s tau) k_++/g_+-,
f_+ = -mu exp(-s tau) k_--/g_+-.
```

Assumimos que f_± são 4pi-periódicas em tau e analíticas no domínio
descrito abaixo. A regra da cadeia dá

```
partial_+ = exp(mu tau)/(2mu) [partial_tau+mu(xi-1)partial_xi],
partial_- = exp(mu tau)/(2mu) [partial_tau+mu(xi+1)partial_xi].
```

Assim k_+++2g_+- partial_+ X^-=0 equivale a T_1(s)y_-=f_-,
e a outra componente equivale a T_-1(s)y_+=f_+, onde

```
T_c(s)=partial_tau+mu(xi-c)partial_xi+s-mu.
```

## Inversa no cilindro inteiro

Escolha uma vizinhança complexa convexa e limitada Omega de [-1,1],
com fecho contido numa vizinhança maior de analiticidade das fontes.
Uma elipse de Bernstein menor serve. Trabalhe na norma

```
||f|| = sum_m exp(sigma |m|) sup_(xi in Omega) |f_m(xi)|,
f(tau,xi)=sum_m f_m(xi) exp(i m tau/2).
```

A hipótese na vizinhança maior garante L=||partial_xi f|| finito
na menor, por Cauchy. Não afirmamos o mesmo bound de derivada sem
perda de domínio. Translações REAIS de tau são isometrias nessa norma
de coeficientes complexos. A conversão para os DOFs reais RT exigirá
seus fatores próprios.

Escreva f_c(tau)=f(tau,c), g=f-f_c e
R_c(s)=(partial_tau+s-mu)^-1 sobre funções periódicas de tau.
Se Re s>0 e s não pertence a {mu-i m/2 : m inteiro}, então

```
y(tau,xi) = R_c(s) f_c(tau)
  + integral_0^infty exp(-(s-mu)t)
      g(tau-t, c+exp(-mu t)(xi-c)) dt.                 (*)
```

O caminho radial permanece em Omega por convexidade. Se
D_c=sup_Omega |xi-c|, o integrando em norma é no máximo
D_c L exp(-Re s t), porque g se anula em c. Portanto (*) converge
também quando 0<Re s<mu; a integral sem subtrair o traço não teria
essa propriedade. Diferenciando o semigrupo sob a integral, ou primeiro
em polinômios Fourier e passando ao limite em domínios menores, obtém-se
T_c(s)y=f. O termo de fronteira no infinito é zero pelo mesmo bound.

Defina d(s)=inf_m |s-mu+i m/2|. A inversa temporal é o multiplicador
1/(s-mu+i m/2), e temos

```
||y|| <= ||f_c||/d(s) + D_c L/(Re s),
||partial_xi y|| <= L/(Re s).
```

Os bounds são locais uniformes em s fora das ressonâncias, e a solução
é holomorfa em s ali. A equação controla a derivada temporal com perda
adicional de faixa, se necessária.

Para unicidade, desenvolva uma solução homogênea em Taylor em xi=c.
O coeficiente de grau n satisfaz
(partial_tau+s+mu(n-1))y_n=0. Para n>=1 sua parte real é positiva;
para n=0 a hipótese de não ressonância dá invertibilidade. Todos os
coeficientes se anulam. A identidade analítica estende isso a Omega.

## Ressonância e cadeias

Em s0=mu-i m0/2, a condição necessária e suficiente é

```
(f_c)_(m0)=0.
```

A necessidade vem de avaliar a equação em xi=c. Para suficiência,
use (*) com o inverso temporal no complemento do modo m0 e fixe em
zero o coeficiente livre desse modo. A solução é única módulo
y0=exp(i m0 tau/2), constante em xi. No complemento, o inverso temporal
tem norma <=2, pois os outros denominadores são i(m-m0)/2.

A equação de cadeia T_c(s0)y1=-y0 não admite solução periódica:
seu traço tem coeficiente ressonante -1. Portanto esse núcleo do
TRANSPORTE não possui cadeia generalizada. Isso não decide cadeias
mistas físicas/gauge do pencil RT, nem permite subtrair multiplicidades
de uma contagem truncada.

Fora das ressonâncias, jets são reconstruídos por
T_c(s0)y_j=f_j-y_(j-1), com y_-1=0. Nas ressonâncias, é preciso impor
a compatibilidade do traço dessa fonte efetiva em cada ordem.
Se um setor RT restringe as frequências temporais, só as ressonâncias
permitidas nesse setor entram; não se pressupõe aqui essa identificação.

## Centro e cone

Para uma perturbação com simetria radial suave, suponha
f_+(tau,xi)=f_-(tau,-xi). A reflexão conjuga T_1 a T_-1.
Pela unicidade (ou pela mesma normalização ressonante),
y_+(tau,xi)=y_-(tau,-xi). Os perfis do vetor em coordenadas RT são

```
A=(y_++y_-)/(2mu),
B=((1+xi)y_--(1-xi)y_+)/2,
X=exp(s tau)(A partial_tau+B partial_xi).
```

A é par e B é ímpar, de modo que o centro é preservado e B/xi é
analítico ali. A construção vale no cilindro inteiro, não apenas num
germe em um ponto temporal. A periodicidade segue diretamente de (*).

Preservar o cone xi=1 equivale a y_-(tau,1)=0. Pelo transporte,
isso exige f_-(tau,1)=0 como FUNÇÃO inteira de tau, e não apenas a
anulação de um coeficiente ressonante. Essa condição também é suficiente:
em (*) o termo de traço é zero; na ressonância escolhemos zero para a
constante livre. Pela reflexão, a condição correspondente vale em xi=-1.

No cone móvel, fora das ressonâncias, o valor no cone é determinado
por R_1(s)f_-(tau,1). Na ressonância a compatibilidade de um coeficiente
ainda é necessária, mesmo permitindo mover o cone. Um gauge secular
em tau resolveria certas fontes incompatíveis, mas sairia da classe de
perfis Floquet ordinários e não foi admitido silenciosamente.

## Alcance e pendências

As hipóteses acima bastam para colocar o BLOCO RADIAL em coordenadas
nulas, com centro e periodicidade globais. A condição de cone foi
explicitada. Não bastam, por si sós, para afirmar equivalência física
com o espaço RT forte: faltam controlar a passagem da métrica e raio
de área às variáveis RT nos pesos exigidos, as normalizações e a classe
de perturbações admissíveis no cone. A reconstrução inversa dentro do
ansatz continua em GEOMETRIC_FLOQUET_RECONSTRUCTION.md.

`inverse_periodic_mode` em `scripts/global_null_gauge.py` verifica essa
solução em cada modo Fourier com polinômios radiais racionais. Os testes
substituem o resultado no operador diferencial, cobrem expoentes em
ambos os lados de mu, aliases ressonantes, obstrução de cadeia, reflexão
e cone fixo. Os antigos testes em polinômios de tau continuam locais.
