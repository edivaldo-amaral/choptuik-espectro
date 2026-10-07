# Localizador espectral sharp e testemunho de violação

Continuação: SHARP_LOCAL_EXCLUSION.md acrescenta uma exclusão racional
para o pencil sharp finito 8×24 perto de s=0.733, nos pesos menores.
Não transforma os diagnósticos desta nota em afirmações sobre o espectro infinito.

## Implementação independente do exportador selecionado

`rt_sharp_matrix.c` exporta o pencil K(s)=K(0)+sI sobre os DOFs sharp de RT.
Usa as três rotinas da implementação original:
JBold_SHARP, XiDiv_GAMMA1_SHARP e GAMMA2_SHARP. O último termo tem fator
UM, não o fator dois do pencil selecionado. Os coeficientes exportados
são diádicos exatos, mas fundo e truncamento continuam aproximados.

O formato tem uma assinatura própria CHOPTUIK_RT_SHARP_MATRIX_V1 para
evitar confundi-lo com uma matriz selecionada. A comparação exige que a
ordem de endereços coincida exatamente com as linhas do arquivo de constraints.
O exportador recusa sobrescrever arquivos existentes.

Os testes verificam o principal no fundo zero, o acoplamento mu R_x c1,
o fator um num fundo constante e a recusa de sobrescrita. A interpretação
polinomial independente da identidade completa continua nos testes anteriores.

## Resultados numéricos

| caixa | dimensão sharp | raiz selecionada candidata a violação | raiz sharp real | distância |
|---|---:|---:|---:|---:|
| 4×12 | 54 | 0.438698175 | 0.401299570 | 3.740e-2 |
| 6×18 | 135 | 0.404385467 | 0.401099535 | 3.286e-3 |
| 8×24 | 252 | 0.401248171 | 0.401028542 | 2.197e-4 |
| 12×36 | 594 | 0.401024411 | 0.401024731 | 3.199e-7 |

São autovalores calculados em ponto flutuante, NÃO intervalos certificados.
No último truncamento, para esse modo, o resíduo relativo
`||(K_N+sI) C_N(s)h|| / ||C_N(s)h||` é aproximadamente 0.0592 na norma
de pesos (129/128,9/8); em 4×12 era aproximadamente 0.734.
Ele ainda não é desprezível. A identidade infinita não comuta com o
truncamento, e o fundo de referência tem resíduo não nulo.

O candidato físico em s≈0.73318 permanece separado da raiz sharp real
por aproximadamente 0.332. A proximidade entre as duas raízes próximas
de 0.401 reforça a hipótese de violação de constraints; NÃO prova que
elas sejam o mesmo modo no limite nem exclui o modo do espaço físico.

Também aparecem raízes sharp complexas com parte real positiva abaixo
de 1/8 e aliases próximos de deslocamentos inteiros imaginários. Os
relatórios preservam essas raízes: não é legítimo declarar que o sistema
sharp tenha só uma raiz instável em todo o semiplano.

Relatórios: `build/spectrum/sharp-comparison-{4x12,6x18,8x24,12x36}.json`.
Eles identificam matrizes por hashes e distinguem normas fraca/fraca e,
na versão atual do produtor, forte/fraca para compatibilidade com SHARP_DOMAIN.

## Como transformar a observação num filtro provado

Não basta subtrair o winding sharp do winding selecionado. Uma alternativa
é encerrar o autovetor selecionado numa bola e provar que sua constraint
não zera. Fixe Y forte=(65/64,5/4), Ysharp fraco=(129/128,9/8), um vetor
aproximado v e expoente s0. Suponha dados CERTIFICADOS

```
||h-v||_Y <= epsilon_v,   |s-s0|<=epsilon_s,   ||v||_Y<=v_norm,
||P C_ref(s0)v|| >= c_ref,
||P(C_true(s0)-C_ref(s0))v|| <= epsilon_applied,
||C_true(s0)||_(Ysharp<-Y) <= c_op.
```

P aqui é uma projeção coordenada finita das constraints, de norma um.
Como C1 h=(0,-Xi h2), seu bound forte→fraco é no máximo 9/8. A desigualdade
triangular implica

```
||C_true(s)h|| >= c_ref - epsilon_applied - c_op epsilon_v
                 - (9/8) epsilon_s (v_norm+epsilon_v).
```

Se o lado direito é estritamente positivo, h não é um modo físico com
constraints nulas. O checker constraint_witness.py implementa essa conta,
mas os enclosures infinitos de h e s ainda NÃO foram fornecidos. Erros de
mu, fundo, truncamento e arredondamento devem entrar nas parcelas corretas.
Os resíduos numéricos da tabela não preenchem essas hipóteses automaticamente.

## Reprodução

```bash
bash scripts/build_sharp.sh
OMP_NUM_THREADS=1 build/sharp-matrix \
  .cache/rt-1203.3766v1/sourcecode/RefA.dat 12 36 build/spectrum/rt-sharp-12x36.dat
OPENBLAS_NUM_THREADS=2 .venv/bin/python scripts/analyze_sharp.py \
  build/spectrum/rt-A-12x36.dat build/spectrum/rt-sharp-12x36.dat \
  --out build/spectrum/sharp-comparison-12x36.json
```

Se a matriz já existe, pule a etapa de exportação: ela não sobrescreve dados.
