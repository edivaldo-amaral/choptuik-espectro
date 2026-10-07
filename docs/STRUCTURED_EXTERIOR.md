# Estrutura das colunas da inversa livre

Continuação: `COMPONENT_EXTERIOR.md` também reduz o bound de B por uma
norma adaptada às componentes. A combinação dá TT<4.882 em 64×192;
é obrigatório converter os demais blocos para essa mesma norma.

## Cancelar antes de estimar

Mantemos Q0=(J_mu+z0)^-1, z0 real não negativo, e os pesos RT.
Para um índice Fourier assinado, ponha b=m/2 e
`D_j=mu(j+1)+z0+ib`. As fatorizações de ALGEBRAIC_TAIL dão
`Q0=V[k,f] D^-1(I-U_k)`, k=1 para J e k=2 para K.
Para uma coluna de grau n, a diagonal é 1/D_n. O primeiro coeficiente
fora da diagonal é EXATAMENTE

```
q_(n-k,n) = -f_(n-k)/D_n - 1/D_(n-k)
           = -2 mu n/(D_n D_(n-k)).
```

Os demais coeficientes descendentes multiplicam esse coeficiente por
produtos de -f_j. Como |f_j|<=1, a soma ponderada da coluna é majorada por

```
1/|D_n| + 2 mu n/(|D_n| |D_(n-k)|) * r^k/(1-r^k),  r=4/5.
```

Se n<k, há somente a diagonal. A multiplicidade do grau zero só diminui
a razão dos pesos. Essa estimativa mantém o cancelamento entre os dois
termos; estimar V, D^-1 e I-U separadamente o perdia.

Como z0 é real, cada coeficiente complexo q é representado nos DOFs reais
pelo bloco [[Re q,-Im q],[Im q,Re q]]. Sua norma de colunas é
`|Re q|+|Im q|<=sqrt(2)|q|<=3|q|/2`, inclusive após complexificar os DOFs.
Não se usa aqui a constante 2 da mudança geral de coordenadas. Para base
complexa esse argumento não foi estabelecido e a função não a aceita.

## Supremo sobre TODO o exterior

Na cauda Fourier |m|>=M, escreva b0=M/2. Para k=1,
`x=mu*n+z0`, temos mu*n<=x e |D_n|>=|D_(n-1)|.
A desigualdade `2 b0 x<=b0^2+x^2` implica o bound
`3/[2 b0(1-r)]` para a coluna.
Para k=2 use `x=mu(n-1)+z0` e mu*n<=x+mu_max. Resulta

```
t_m = (3/2) max(1/[b0(1-r)],
  1/[b0(1-r²)] + 2 mu_max r²/[b0²(1-r²)]).
```

Para n>=N>=2, ignore b e z0 nos denominadores. A parte J é
`(3/2)*9/[mu(n+1)]`. Na parte K, o fator adicional é
`1+32n/[9(n-1)]<=9`, pois n>=2. Logo

```
t_n = 27/[2 mu_min(N+1)],
||Q0(I-G)|| <= max(t_m,t_n).
```

O supremo foi dividido em duas regiões que cobrem o exterior infinito;
não é um máximo tomado somente sobre colunas amostradas.
O mínimo desse bound com o anterior continua válido. Para os parâmetros
RT e M>=1, t_m=15/M e t_n=81/(N+1).

## Efeito nos blocos uniformes

Substituir t_in pelo novo bound nas três desigualdades da seção 7 de
FREQUENCY_BLOCKS é legítimo. Nenhum outro fator muda: a melhora é uniforme
em |s|<=5/4 com Q0 fixo e z0=0.

| G | defeito TT anterior | novo defeito TT |
|---|---:|---:|
| 24×72 | 2619/71 | 7857/292 |
| 64×192 | 2619/191 | 7857/772 |
| 364×1964 | >1 | 7857/7860 <1 |

Em 64×192, o majorante cai de aproximadamente 13,712 para 10,177.
A primeira caixa que satisfaz ambas as desigualdades estritas desse
majorante tem M=364, N=1964 (sem otimizar paridades). Só o bloco TT passa;
os acoplamentos e a parametriz finita ainda impedem concluir contração.
Esse resultado melhora a constante, não resolve a dificuldade estrutural
remanescente de B. A próxima redução deve absorver uma parte significativa
de B em uma parametriz exterior, com resto certificado; a norma global
23 ainda descarta os sinais e a estrutura entre componentes de B.

Implementação: `scripts/structured_exterior.py`. Prova analítica em papel,
testes racionais em colunas finitas; não formalizado em Lean.
