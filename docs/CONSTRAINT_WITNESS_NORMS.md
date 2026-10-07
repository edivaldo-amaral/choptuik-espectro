# Normas do testemunho de constraints: a hipótese H3 (c_op racional) e c_ref

Sub-agente FÍSICO, rodada 3 (17/09/2026). O que é exato está marcado como **EXATO**;
o que é ponto flutuante, como **diagnóstico**.

Código: `scripts/independent_constraint_operator_bound.py`. Testes (0,2 s):
`scripts/test_independent_constraint_operator_bound.py`.

```bash
OPENBLAS_NUM_THREADS=2 .venv/bin/python scripts/independent_constraint_operator_bound.py --out <json>
.venv/bin/python -m unittest discover -s scripts -p 'test_independent_constraint_operator_bound.py'
```

## 1. O testemunho e as normas

Para uma raiz s0 de L e um vetor aproximado v, com ‖v‖=1:

```
‖Π_F C_true(s) h‖ >= c_ref - ε_apl - c_op ε_v - ‖C1‖ |s-s0| (1+ε_v),
c_ref = ‖Π_F C_ref(s0) v‖,   ‖h-v‖ <= ε_v,
```

com Π_F a projeção coordenada numa caixa de saída F (norma 1).

- **Entrada Y+:** pesos fortes (65/64, 5/4), **sem η**. É a norma RT de SSPACE, soma
  sobre componentes.
- **Saída Y♯-:** pesos sharp menores (129/128, 9/8).
- **Soma:** normas l1 com soma sobre ℤ², como em RT.
- **Complexificação:** módulo em cada coordenada real. Assim ‖s x‖ = |s|‖x‖, e a
  norma de operador é o supremo das normas l1 das colunas reais.

**Por que projetar a saída.** O testemunho só precisa de C h ≠ 0, e para isso basta
uma cota inferior em QUALQUER norma. Projetar numa caixa finita F elimina a perda de
raio de Fourier (∂τ) fora de F. Isso reduz c_op sem exigir injetividade sharp nos
pesos menores. (A injetividade de K nesses pesos só é necessária na rota oposta,
"K injetivo ⇒ constraints satisfeitas".)

## 2. Decomposição exata de C(s)

**EXATO, universal:** `test_C_equals_D_plus_A_universally` usa o método de jatos de
`independent_rt_symbolic.py`, a partir do TeX de RT; a mutação que remove o termo
ξ²∂ξ falha.

```
C(s)h = D(s)h + A_ω h,
D(s)h = ( μ(ξ∂ξ+2)h1 ,  -ξ(∂τ+s+μ)h2 + μ(ξ∂ξ - ξ²∂ξ)h2 + 2μ(h2-h3) ),
A_ω h = ( 2 O(ξ Γ2_1(ω,h)) ,  -2 P(ξ Γ2_3(ω,h)) ),
```

onde O=(1-P)/2 e Γ2_1, Γ2_3 são as componentes (Ω1, Ω2♯) de Γ2 do TeX. Logo
`C(s) = C(0) + s C1`, com `C1 h = (0, -ξh2)` e `‖C1‖_{Y♯-←Y+} = 9/8`. O máximo vem
da coluna n=0, porque ξT_0=T_1.

## 3. Cota racional de ‖C_true(0)‖ (**EXATO**)

Coluna unitária e_{k,p,m,n}, componente de entrada k:

- **Parte D (sem fundo).** É diagonal em m. As razões de coluna são calculadas em
  `Fraction` para n<=200, com majorantes decrescentes para n>=201 (testados). Os
  fatores são (129/130)^m <= 1 e (9/10)^n.
  - A derivada em τ contribui `sup_m (m/2)(129/130)^m`. Sem projeção, o máximo é em
    m<=130 e vale 23,8199; é exato, pois a sequência decresce para m>=130.
  - A parte em μ usa μ* = μ_RefA + ε_μ, que majora μ_ref e μ_true.
  - Para k=2 a parte D vale exatamente
    `(129/130)^m [ (m/2)‖ξe_n‖ + μ ‖(ξ∂ξ - ξ²∂ξ + 2 - ξ)e_n‖ ] / ‖e_n‖` (em s=0).
- **Parte A_ω.** Cota uniforme via ‖f*g‖ <= ‖f‖‖g‖ em Y♯-, ‖O‖, ‖P‖ <= 1,
  ‖ξ‖_{Y♯-} = 9/8 e ‖e‖_{Y♯-} <= ‖e‖_{Y+}.
  - As normas das combinações do fundo (por exemplo 3Pω2-2Pω3+ω1-ω2) são somas
    racionais exatas sobre RefA.
  - A bola RT soma (Σ|coef|)·ε_ω, com ε_ω = 2^-25+2^-277 por componente.
  - Com isso A_ω cobre C_true; a diferença C_true-C_RefA entra em ε_apl (último item desta seção).

| caixa de saída | sup_m (m/2)(129/130)^m | D: k=1 / k=2 / k=3 | A_ω: k=1 / 2 / 3 / 4 | c0 >= ‖C_true(0)‖ | **c_op(\|s\|<=5/4) = c0 + (5/4)(9/8)** |
|---|---:|---|---|---:|---:|
| sem projeção | 23,8199 | 3,810 / 27,323 / 0,337 | 1,826 / 11,776 / 7,021 / 15,462 | 4887429/125000 = 39,0994 | **20252841/500000 = 40,5057** |
| m<12 | 5,0521 | 3,810 / 8,997 / 0,337 | idem | 10386527/500000 = 20,7731 | **2772413/125000 = 22,1793** |
| m<4 | 1,4657 | 3,810 / 7,979 / 0,337 | idem | 19755347/1000000 = 19,7553 | **21161597/1000000 = 21,1616** |

- **Faixa de s.** O retângulo B (Re s<=1, |Im s|<=3/4) cabe em |s|<=5/4: o canto
  1+0,75i tem módulo exatamente 5/4. A faixa A e o candidato B
  (|0,0414+0,5i| = 0,502) também cabem. Para |s| maior vale `c_op(s) = c0 + (9/8)|s|`.
  Nas raízes: c_op(0,4010) = c0 + 0,4512 e c_op(0,0414+0,5i) = c0 + 0,5645.
- **Folga (diagnóstico).** A maior norma de coluna de Π C0 na caixa 12×36 (entradas
  exportadas) é 8,719. Ela majora por baixo o verdadeiro c0 (m<12) e está na coluna
  (h2, m=10, n=10). Portanto a cota exata 20,77 tem folga <= 2,4×. A parcela
  dominante é A_ω, estimada pela desigualdade triangular.
- **ε_apl (EXATO):** ‖Π_F (C_true - C_RefA)(s0) v‖ <= (81/8)ε_ω + 45,19·ε_μ < 3,04e-7
  para ‖v‖=1. É desprezível.

## 4. c_ref no truncamento 12×36 e a razão c_ref/c_op (**diagnóstico**)

v é o autovetor da matriz 12×36, com ‖v‖_{Y+}=1 (sem η) e c = (C0 + s0 C1)v.

| raiz | c_ref, F = 12×36 | c_ref, F = 4×12 | c_ref/c_op (sem projeção) | c_ref/c_op (m<12) | c_ref/c_op (m<4, F=4×12) |
|---|---:|---:|---:|---:|---:|
| A: 0,401024411 | 0,02631 | 0,02277 | 6,7e-4 | **1,24e-3** | 1,13e-3 |
| B: 0,041420731+0,500008241i | 0,11192 | 0,09756 | 2,8e-3 | **5,25e-3** | 4,80e-3 |

`c_ref/c_op` é o ε_v máximo em Y+, relativo a ‖v‖=1, com ε_apl e |s-s0| desprezados.
O valor 0,02631 coincide com o `constraint_ratio_strong_to_weak` do
`sharp-comparison-12x36.json`.

## 5. O testemunho é alcançável? Não com a máquina da rodada 2

Dois números decidem, ambos diagnósticos em float.

1. **Erro do próprio centro v.** Comparei os autovetores 10×30 e 12×36, com fase
   alinhada, nos DOFs comuns mais a massa fora deles, na norma Y+:
   - 0,401: ‖v12 - v10‖ = 3,8e-2; massa de v12 na casca m>=10 ou n>=30 = 2,2e-2;
   - raiz B: 5,6e-2 e 2,4e-2.

   O erro de truncamento do centro já é 30× (A) e 10× (B) maior que o ε_v tolerável.
2. **Amplificação do enclausuramento.** Um par autovalor/autovetor do operador
   infinito (H2) sai de Newton–Kantorovich/Krawczyk com ε_v ≈ ‖H^{-1}‖·resíduo. Com as
   margens de Schur da rodada 2 (1-0,99876 ≈ 1,2e-3, κ≈31), ‖H^{-1}‖ é da ordem de
   10⁴. Seria preciso um resíduo <= 1e-7, contra um erro de centro ~1e-2 no 12×36. O resíduo do centro no operador infinito não foi medido.

A **rota (b)**, com h e C h nos pesos fracos e saída projetada, não resolve sozinha:

- c_ref sobe (0,0559 para A; 0,230 para B);
- o erro do centro cai para 1,3e-2 e 2,3e-2;
- mas as normas de coluna de Π C0 de fraco para fraco chegam a 102, porque sem perda
  de raio ∂τ e ξ∂ξ não são amortecidos;
- a razão fica em ~5e-4, sem ganho.

**Conclusão.** H3 está resolvida com constante racional explícita: c_op <= 40,51
(saída sem projeção) ou <= 22,18 (saída m<12), uniforme em |s|<=5/4 e na bola RT.
O testemunho é possível em princípio, porque a margem c_ref é robusta e não cai com
o refinamento (AUDITORIA_MATEMATICA §(d)). Mas exige:

- ε_v <= ~1e-3 (raiz A) ou ~5e-3 (raiz B) em Y+;
- um centro ≥30× melhor que o autovetor 12×36;
- uma inversa do problema de autovalor com amplificação muito menor que a das margens
  de 10⁻³.

**Nenhum desses requisitos é atendido hoje.**

## 6. O que falta, por hipótese

| hipótese | estado |
|---|---|
| **H1:** multiplicidade 1 da família infinita num disco | aberta; exige winding 1 com cauda e bola RT em torno de 0,401 e de 0,0414+0,5i |
| **H2:** enclausuramento do par autovalor/autovetor infinito, com ‖v‖>ε_v | aberta; é o gargalo quantitativo da §5 |
| **H3:** C(s) limitado com perda de raio, c_op explícito | **resolvida nesta nota** (EXATO); folga <= 2,4× |
| **H4:** identificação com a multiplicidade física, módulo gauge e cone | aberta (S4); para 0,401 e para a raiz B não há gauge a subtrair (B6.6) |

## 7. Próximo passo sugerido, se esta rota continuar

1. **Reduzir c_op.** Trocar a cota triangular de A_ω por normas de coluna exatas (como
   os β do Matemático), com ganho esperado de até 2,4×. Escolher a caixa de saída F
   onde c_ref/c_op é máximo.
2. **Obter um centro melhor.** Autovetor em caixa maior (24×72 ou via
   iteração inversa com o operador de colunas exatas) com ‖v - h‖_{Y+} <= 1e-4.
   Critério: a diferença entre duas caixas sucessivas deve cair abaixo de 1e-4.
3. **Enclausurar só depois das margens.** O enclausuramento só é realista com d*(G)
   de ordem 0,5, não 0,998. Portanto depende da parametriz da faixa U ou de caixas
   maiores, e isso é trabalho do Matemático.
