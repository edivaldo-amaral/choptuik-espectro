# Exclusão sharp local: componente finita e extensão necessária

Continuação: SHARP_EXTERIOR_BOUNDS.md incorpora o erro do fundo no mesmo
disco finito e constrói os bounds próprios da cauda sharp. A exclusão
infinita acoplada ainda não foi obtida. CONSTRAINT_ROOT_QUOTIENT.md dá
o mapa D=C0-C1 A para contar constraints nas raízes generalizadas.

## Disco finito auditado nos pesos corretos

`certify_sharp_disk.py` auditou o pencil sharp 8×24 (dimensão 252), nos
pesos (129/128,9/8), sobre TODO o disco complexo

```
|s-733/1000| <= 1/32.
```

A matriz vem de RefA. A inversa candidata no centro foi arredondada
para diádicos de 48 bits. Depois de conjugar pelos pesos racionais,
o arredondamento de cada entrada da matriz tem erro <=2^-49.
O produto residual é calculado em inteiros exatos. A cota no disco é
`q=||I-V Khat(s0)||+||V||[252*2^-49+1/32] <0.264`, e
`||K_N(s)^-1||<=||V||/(1-q)<11.446`.

O relatório `build/spectrum/sharp-disk-8x24.json` preserva os racionais,
o hash e os marcadores de alcance. Essa exclusão cobre o candidato
físico NUMÉRICO próximo de 0.73318, mas não prova que o autovalor do
operador selecionado infinito esteja no disco. Tampouco inclui erro do
fundo ou cauda sharp: não é ainda injetividade do operador sharp infinito.

```bash
OPENBLAS_NUM_THREADS=2 .venv/bin/python scripts/certify_sharp_disk.py \
  build/spectrum/rt-sharp-8x24.dat --out build/spectrum/sharp-disk-8x24.json
```

## Condição de extensão, sem trocar o espaço

Precondicione o sistema sharp por uma inversa livre FIXA e decomponha-o
em F e T nos mesmos pesos menores. Para transportar esta exclusão, é
necessário validar o defeito finito da parametriz PRECONDICIONADA e os
blocos b,c,d de UNIFORM_COUPLED_CONTRACTION, incluindo a mudança do fundo.
O bound 11.446 da inversa de K_N não é, por si, o bound de uma inversa
de K_N Q_N: esta última é Q_N^-1 K_N^-1. Confundir essas duas matrizes
produziria uma constante incorreta. A desigualdade a+bc/(1-d)<1 ainda
precisa ser verificada, com bounds sharp e pesos menores, não os de A.

## Constraints de cadeias, não só de autovetores

Escreva L(s)=L0+(s-s0)I, K(s)=K0+(s-s0)I e
C(s)=C0+(s-s0)C1 no fundo exato, com KC=ML. Para uma cadeia
L0 h_j+h_(j-1)=0 (h_-1=0), defina

```
c_j=C0 h_j+C1 h_(j-1).
```

Comparando coeficientes das séries, a identidade completa dá
`K0 c_j+c_(j-1)=0`. Se K0 for injetivo na realização de raio menor,
indução implica c_j=0 para todos os j. Não se deve exigir apenas
`C0 h_j=0`: para j>0 existe o termo C1 h_(j-1).
O argumento vale para jets finitos, com a compatibilidade de domínios
de SHARP_DOMAIN aplicada a cada coeficiente e às derivadas do pencil.

Se K0 tiver kernel, a identidade só define um mapa entre cadeias;
não prova que esse mapa seja sobrejetivo. Assim não é legítimo subtrair
a contagem TOTAL sharp da contagem selecionada. É preciso certificar
a imagem efetiva das constraints nas cadeias selecionadas, ou excluir
sharp na região de interesse e tratar separadamente os discos onde ele
pode ser singular. O disco acima inicia a segunda rota; a raiz sharp
perto de 0.401 continua exigindo uma análise distinta.
