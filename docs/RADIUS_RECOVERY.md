# Recuperação do raio nos espaços de raízes selecionados

Este argumento é para a realização RT fechada, não uma equivalência
global com todas as perturbações geométricas. Usa a existência e os bounds
publicados do fundo, os lemas de inversa livre e D_min=D_max já registrados.

## Duas realizações compatíveis

Sejam Y+ os pesos (65/64,5/4) e Y- os pesos (129/128,9/8), ambos com
eta=(916,4096,243,1233). Temos inclusão contínua Y+ em Y-.
Nesta prova use Q=J_mu_true^-1, NÃO o precondicionador de RefA usado nas
auditorias numéricas. Ele é limitado em ambos os espaços, coincide na
interseção e preserva retângulos. A diferença de mu não é uma perturbação
aqui: mu é o mesmo parâmetro exato nas duas realizações.

O bound de B nos pesos fortes é beta+<=5191/500. Nos pesos menores,
a norma das multiplicações não aumenta. A norma de divisão regularizada
por xi muda de 40/9 para 144/17, razão 162/85. Portanto, reponderando
as mesmas componentes, basta beta-<=(162/85)beta+.

Para r=8/9, a prova por colunas de STRUCTURED_EXTERIOR dá

```
||Q(I-G)||_- <= max(27/M,153/(N+1)).
```

A fórmula geral, incluindo o bloco K e seu fator N/(N-1), é implementada
em sharp_exterior.free_tail; vale igualmente para os quatro blocos
selecionados. Não depende do peso Fourier porque Q preserva m.
Assim, para H(s)=L(s)Q=I+(B+s)Q e |s|<=5/4, G=1024×4096 dá

```
||T(H-I)T||_+ <0.230,
||T(H-I)T||_- <0.786,     T=I-G.
```

São bounds de cauda analíticos, não inversões de uma matriz dessa dimensão.
`radius_transfer.py` verifica os racionais; a faixa finita inteira é
mantida como dado na prova seguinte, sem hipótese de invertibilidade.

## Autovetores: existência forte e unicidade fraca

Se h pertence ao kernel de L(s) em D-, ponha x=J_mu h em Y-.
Então H(s)x=0. A parte Gx é finita, portanto pertence a Y+.
No espaço forte resolva, por Neumann,

```
y_T=-(H_TT,+)^-1 H_TG Gx.
```

O lado direito pertence a Y+: Gx é finito e H é limitado em Y+.
Pela compatibilidade das realizações, y_T também satisfaz a mesma
equação no espaço fraco. Mas Tx é a solução fraca, que é ÚNICA porque
H_TT,- é invertível. Logo Tx=y_T, x pertence a Y+, e h=Qx pertence a D+.
A inclusão inversa é imediata. Os kernels coincidem como conjuntos.

## Cadeias completas

Para L(s0)h_j=-h_(j-1), suponha por indução que h_(j-1) está em Y+.
Com x_j=J_mu h_j, a equação exterior é

```
H_TT T x_j = -T h_(j-1)-H_TG Gx_j.
```

O lado direito é forte. Existência forte e unicidade fraca repetem o
argumento. Portanto TODA cadeia finita recupera o raio forte, com os
domínios correspondentes. Não precisamos supor que derivar um fundo em
Y+ preserve Y+: usamos a equação satisfeita pelo vetor derivado.

## Consequência para o gauge

GAUGE_RT_ACTION constrói o vetor não nulo
h=(1+mu^-1 partial_tau+xi partial_xi)omega* primeiro em Y- e prova
L(mu)h=0, C(mu)h=0. Como mu<5/4, o argumento acima coloca h em D+.
Fecha-se essa pendência de admissibilidade na realização RT forte.
O vetor move o cone: sua inclusão ou exclusão no domínio GEOMÉTRICO
continua dependendo de cone móvel versus fixo. Não se subtrai uma raiz
de um certificado finito, nem se declara a equivalência física completa.

## Fredholm e limite da conclusão

Os bounds de cauda tendem a zero: Q é limite em norma dos operadores
de posto finito QG, portanto compacto em cada realização. H=I+(B+s)Q
é uma família analítica Fredholm de índice zero. Via a isomorfia
Q:Y->D(J), a afirmação passa para L:D->Y. A invertibilidade em algum
s real suficientemente grande vem do bound livre deslocado e Neumann;
assim não estamos no caso de uma família singular em todo ponto.
As raízes são isoladas e têm multiplicidade algébrica finita.

Com a coincidência das cadeias, as multiplicidades selecionadas coincidem
nos dois raios no disco estudado. Nada disso fornece seus valores, a
ausência de outras raízes, ou a invertibilidade sharp. A prova é analítica
em papel, apoiada em aritmética racional Python, não formalizada em Lean.
