# O disco D de S3b (24/09/2026)

**Escolha: `D = {|s − 0,401| < 1/8}`.**

S3a exclui `K` fora de D; S3b precisa mostrar que, dentro de D, `L` tem **uma única** raiz,
algebricamente simples, e que o autovetor dela viola as constraints (obrigação H1 de
`CONSTRAINT_WITNESS_NORMS.md` e testemunho). O raio decide duas coisas: quantas raízes de
`L` ficam dentro de D, e quão bem condicionado `L` é em ∂D, onde o *winding* de S3b vai ser
certificado. `scripts/measure_s3b_disk_radius.py` (ponto flutuante, diagnóstico).

## 1. As raízes de L perto de 0,401

Estáveis de 8×24 a 12×36:

| raiz de L | distância a 0,401 |
|---|---:|
| **0,4010** | 0 |
| 0,1683 | 0,233 |
| 0,7332 | 0,332 |
| 0,0000 (modo de fase) | 0,401 |

O raio precisa ficar bem abaixo de 0,233. (0,1683 e 0,7332 estão fora de D; por S3a,
satisfazem as constraints, e classificá-las como gauge ou físicas é tarefa de S4.)

## 2. O condicionamento de L em ∂D

Máximo de `||V(s)|| = ||J (op + s)⁻¹||_w` na circunferência `|s − 0,401| = r` (a norma que
dimensionou a rota certificada de L; 49 pontos na metade superior):

| `r` | L 10×30 | L 11×33 | L 12×36 | K 12×36 |
|---:|---:|---:|---:|---:|
| 0,02 | — | 4192 | 4253 | 134 |
| 0,10 | — | 945 | 959 | 28 |
| 0,110 | 857 | 891 | 906 | 25 |
| **0,125** | **833** | **866** | **880** | **22** |
| 0,135 | 840 | 873 | 888 | 21 |
| 0,15 | — | 926 | 942 | 19 |

Até ~0,12 domina o polo (`~1/r`); acima disso, a raiz vizinha 0,1683. O máximo em `r = 1/8`
cai no eixo real, em `s = 0,276`, o lado voltado para 0,1683. O ótimo é largo
(0,11–0,14 fica a 3% do mínimo), e `1/8` é diádico, conveniente para aritmética exata no
contorno. Para comparação, o contorno `Γ_A = [1/8, 1] × [−1/4, 1/4]` tinha máximo 1846: o
círculo é ~2× mais bem condicionado.

Folgas: 0,108 até a raiz vizinha de L; 0,207 até 0,7332; `K` bem condicionado em ∂D (22).

## 3. Compatibilidade com S3a

A cobertura de S3a foi feita com um disco provisório de raio 0,05. Uma cobertura da região
menos um disco menor cobre também a região menos o disco maior, e
`check_cover.py build/cobertura_lab --disco 0.401,0.125` dá **COBERTURA COMPLETA**
(154 discos). Nada a recalcular.

## 4. O que fica para S3b

- *Winding* 1 de `L` em ∂D, certificado (o custo escala com `max ||V_L|| ≈ 880–900`).
- Testemunho: `C(0,401) h ≠ 0` para o autovetor, com enclosure certificado.
- Multiplicidade algébrica 1 (segue do *winding* 1).
