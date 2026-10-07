# Bounds sharp próprios, nos pesos de propagação

Continuação: SHARP_INTERMEDIATE_SHELL.md incorpora a faixa 8×24->12×36,
valida seu acoplamento à cauda infinita distante e identifica o Schur
finito restante. A exclusão sharp completa continua pendente.

Todos os bounds usam (a,b)=(129/128,9/8), sem reponderação entre as
duas componentes sharp. Não usamos beta do operador selecionado aqui.

## Norma do termo limitado e incerteza

A tabela GAMMA2_SHARP, restrita às constraints, usa somente k=1,4,5.
Para cada saída i, a coluna da primeira entrada tem majorante
||L_i1(w)||/2; a segunda tem max(||L_i4(w)||,||L_i5(w)||)/2,
pois as partes par e ímpar da entrada têm suportes disjuntos na norma l1.
Somamos ainda mu*(144/17) na entrada (segunda saída, primeira entrada),
oriunda da divisão regularizada de uma função ímpar por xi.

Essa construção mantém os cancelamentos de L_ik ANTES de tomar módulos.
`sharp_exterior.py` calcula o majorante racional e acrescenta 3 epsilon_omega,
dando beta_sharp<=4.62237. O fator de Gamma2sharp é UM, não dois.
O Lipschitz 3 segue diretamente das colunas da tabela: seus majorantes
são no máximo 3 vezes a soma das normas dos campos do erro. A passagem
dos pesos fortes aos fracos não aumenta essa soma.

Na compressão 8×24, se j_N é a norma ponderada da derivada espacial do
principal em relação a mu,

```
||K_true,N-K_ref,N|| <= (j_N+144/17)epsilon_mu+3 epsilon_omega
                     <1.040e-7.
```

O produtor calcula j_N exatamente; remove o bloco temporal real/imaginário,
que não depende de mu. O corte do fundo em 2N do exportador não altera
a compressão, porque Gamma2sharp é multiplicação sem divisão do produto.
Não estendemos j_N ao principal infinito: lá usamos erro RELATIVO.

## Disco finito agora com fundo incerto

`sharp-disk-8x24-background.json` usa o mesmo disco |s-.733|<=1/32:
o defeito é <0.263453 e a inversa <11.447. Condicionado à bola RT,
a compressão do pencil sharp do fundo exato é invertível nesse disco.
A cauda infinita e seus acoplamentos NÃO estão incluídos nesse certificado.

## Cauda do precondicionador livre sharp

Para Qsharp=Jsharp_muRef^-1, a prova por colunas vale com r=8/9.
Em geral, com b0=M/2, o termo Fourier é

```
(3/2) max(1/[b0(1-r)],
          1/[b0(1-r²)]+2 mu_max r²/[b0²(1-r²)]).
```

O termo radial é (3/2)/(mu_min(N+1)) vezes
max(1+2r/(1-r), 1+2N r²/((N-1)(1-r²))).
Nos parâmetros presentes o máximo final é max(27/M,153/(N+1)).
São bounds para TODA a cauda, não uma amostragem finita.

Com ||partial_tau Qsharp||<=34 e muRef>=1/6,
||Jslash Qsharp||<=(1+34)/muRef<=210. Assim no disco estudado, |s|<4/5,

```
||T(Hsharp-I)T|| <= (4.62237+4/5)t_sharp(G)+210 epsilon_mu.
```

O relatório sharp-exterior.json dá 33.185 em 8×24, 4.299 em 64×192,
1.144 em 128×768 e <0.572 em 256×1536 (valores arredondados para cima).
Só o último satisfaz o teste diagonal listado. NÃO combinamos a inversa
8×24 com a cauda fora de 256×1536 omitindo a região intermediária.
A contração sharp acoplada continua aberta.

```bash
.venv/bin/python scripts/sharp_exterior.py \
  .cache/rt-1203.3766v1/sourcecode/RefA.dat --out build/spectrum/sharp-exterior.json
OPENBLAS_NUM_THREADS=2 .venv/bin/python scripts/certify_sharp_disk.py \
  build/spectrum/rt-sharp-8x24.dat --include-background \
  --out build/spectrum/sharp-disk-8x24-background.json
```

Hipóteses: bola RT e lemas analíticos de operadores. Os produtos finitos
são auditados em inteiros/racionais Python, não em Lean.
