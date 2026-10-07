# Acoplamento com o exterior infinito usando uma faixa de separação

Fixe F=6×18, P=4×12 e S=F-P. Para uma caixa G contendo F, ponha
U=G-F e T=I-G. Usamos os pesos RT fortes, eta=(916,4096,243,1233),
Q_ref fixo e |s|<=5/4.

## Acoplamento de saída por componentes

Q_ref preserva F. Os blocos T J_true Q_ref F, T s Q_ref F e
T mu S Gamma1 Q_ref F são zero: esses operadores preservam retângulos.
Resta a multiplicação pelo fundo. Para F=(M,N) e G=(M',N'), apenas
coeficientes do fundo com m>=M'-M+1 OU n>=N'-N+1 podem levar F a T.
Isso segue do suporte de produtos Fourier–Chebyshev; reflexão e projeção
de paridade não aumentam os índices.

Filtramos primeiro esses coeficientes e depois construímos o majorante
4×4 de COMPONENT_EXTERIOR, preservando cancelamentos entre campos antes
do módulo. Se beta_guard é sua norma eta,

```
||T H F|| <= [beta_guard + 6(4096/243)epsilon_omega] ||Q_ref F||.
```

A auditoria livre finita dá ||Q_ref F||<5.942; usá-la aqui é legítimo
porque F é invariante. Não substituímos ||Q_ref T|| por essa norma finita.
O erro do fundo não é filtrado: sua localização espectral é desconhecida.

Há ainda uma melhora GLOBAL para a inversa livre: a cauda de entrada
fora de F tem norma <=max(15/6,81/19)=81/19, menor que o bound finito
5.942. Como a norma é l1 ponderada e a divisão é por colunas,

```
||Q_ref|| <= max(||Q_ref F||, ||Q_ref(I-F)||) < 5.942.
```

Aqui o máximo, em vez da soma, é essencial: as normas das duas partes
do vetor de entrada somam sua norma total. Isso não pressupõe que Q_ref
preserve o exterior. Com esse bound no orçamento do fundo,

```
||H_true-H_ref|| <=114 epsilon_mu
  +||Q_ref||[g1_eta epsilon_mu+6(4096/243)epsilon_omega]
  <0.000017918.
```

É cerca de 18 vezes menor que o orçamento antigo. Os tiles anteriores
permanecem válidos com seus bounds conservadores; não foram reescritos
nem melhorados por subtração de erros já arredondados.

| G | majorante de ||T H F||, arredondado para cima |
|---|---:|
| 6×18 | 59.947 |
| 12×36 | 2.531 |
| 24×72 | 0.003579 |
| 64×192 | 0.000017909 |
| 176×944 | 0.000017909 |

São bounds uniformes com saída INFINITA, condicionados à bola RT,
não máximos tomados somente sobre colunas amostradas.

## Subsistema coroa mais exterior distante

Sejam e=114 epsilon_mu, d=(5191/500+5/4)t_structured(G)+e,
c o bound acima e kappa=669414234798400/18478486838353.
Temos ||H_TT-I||<=d, ||H_ST||<=d, ||H_TS||<=c e
||H_SS^-1||<=kappa. O termo e mantém J_true-J_ref separado de B_true.
O complemento de Schur em T satisfaz

```
||H_TT-H_TS H_SS^-1 H_ST-I|| <= d+kappa*c*d.
```

Para G=176×944, o lado direito é aproximadamente 0.997675403828,
racionalmente menor que 0.997676. Logo o complemento é invertível por
Neumann, com inversa <431. A fatorização triangular prova invertibilidade
de H comprimido a S+T em TODO o contorno. Reduzir ambos os acoplamentos
S/T por t em [0,1] também preserva a desigualdade (a correção recebe t²).

Isso controla um subsistema infinito, NÃO H inteiro. O complemento
espacial restante é P+U, com U=(176×944)-(6×18). Falta construir/controlar
seu operador de Schur e seus acoplamentos; U nunca foi considerado vazio.
A contribuição espectral de S tampouco foi descartada.

## Reprodução e confiança

```bash
.venv/bin/python scripts/guarded_exterior.py \
  .cache/rt-1203.3766v1/sourcecode/RefA.dat \
  build/spectrum/crown-contour-6x18-background.json \
  --out build/spectrum/guarded-exterior.json
```

O script verifica hashes de RefA/matriz e a cobertura escalar da coroa.
Os resíduos matriciais desse relatório dependem da auditoria anterior.
As estimativas analíticas e a bola RT permanecem hipóteses externas ao
verificador racional Python; não há novo teorema Lean de operadores.
