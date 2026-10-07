# Pré-condicionador de casca: a rota está fechada nesta máquina

Rodada 6, 17/09/2026. DIAGNÓSTICO em float na caixa reduzida 12×36, nos três
centros do protótipo. Não é prova; é a medida que decide para onde vai o
esforço. Scripts: `scripts/shell_preconditioner_probe.py`,
`scripts/shell_balance_floor.py`. Relatórios:
`build/spectrum/shell-preconditioner-probe.json`,
`build/spectrum/shell-balance-floor.json`.

Contexto: a §7 de `G_PARAMETRIX_FEASIBILITY.md` deixou a arquitetura (c) a
**~57×** do critério de 2 h por região, e apontou o pré-condicionador de casca
como a maior alavanca restante. Esta nota mede essa alavanca.

## 1. O grau k=32 não era desperdício

O custo por região é `colunas × k × custo_por_aplicação`, e o alvo do defeito
da casca **não** é o 0,9 do pré-registro: o raio do disco de validade uniforme
escala com `(1-α)`, então baixar o alvo encolhe a região e multiplica o número
de regiões (que já tem teto de 200). O ponto de operação real é α≈10⁻⁵.

Medindo `k_min(α)` na norma de coluna, com o amortecimento ω ótimo de cada
candidato (pior centro, núcleo 6×18):

| candidato | ρ | pico | k(α=0,9) | k(α=10⁻⁵) | ganho |
|---|---:|---:|---:|---:|---:|
| Neumann amortecida (rota atual) | 0,58 | 5,9 | 9 | **29** | 1,00× |
| Jacobi | 0,56 | 5,9 | 9 | 27 | 1,07× |
| bloco-n | 0,54 | 5,8 | 9 | 26 | 1,12× |
| deflação espectral, posto 8 | 0,47 | **185** | 9 | 26 | 1,12× |
| deflação espectral, posto 32 | 0,37 | **1813** | 10 | 23 | 1,26× |
| **bloco-m** | **0,34** | 5,3 | 5 | **15** | **1,93×** |

O k=29 medido confirma o k=32 do documento: a escolha estava certa.

**Deflação é ativamente ruim.** Ela baixa o raio espectral, mas os autovetores
de H_OO são fortemente não ortogonais e o projetor espectral é malcondicionado:
o transiente sobe a 10²–10⁴ antes de cair. Como o que se certifica é a norma de
coluna, e não o raio, o pico é que manda. Isso também explica por que ρ é um
guia enganoso aqui: ρ≈0,5 convive com normas de coluna ≈58.

## 2. O teto da família inteira fica abaixo do necessário

Mesmo um pré-condicionador **perfeito e de graça** (P=H_OO⁻¹ exato, k=1) só
divide o custo por 29–32. Isso dá 228/32 ≈ 7,1 h-núcleo por região, ou **3,6 h
de relógio em 2 núcleos: ainda 1,8× acima do critério de 2 h**, e sem contar o
custo de construir o pré-condicionador.

O melhor candidato real, bloco-m, dá 1,93×, e ainda exigiria inverter de forma
certificada 107 blocos de ~3072 DOFs por região na caixa G — custo de montagem
que provavelmente come o ganho.

**Contra os 57× necessários, a família toda não serve.**

## 3. A rota k=0 também está fechada, e por um piso exato

Restava a cota analítica estilo RT (k=0, nenhum trabalho por coluna), que
baixaria o custo por região para ~0,5 h. Ela exige colunas de X=(B+s)Q_ref com
norma ≲0,6; hoje são ≈58. A pergunta natural é se outra família de pesos
resolveria — afinal os pesos η_d·(65/64)^m·(5/4)^n são uma escolha.

A resposta é **não, e o limite é exato**. A norma do projeto é
`max_j Σ_i (|Re X_ij| + |Im X_ij|)`, ou seja a norma 1 da matriz não negativa
`Bm = |Re X| + |Im X|`. Trocar pesos por DOF é exatamente uma similaridade
diagonal positiva, que leva `Bm → D⁻¹ Bm D` entrada a entrada. Por
Perron-Frobenius,

```
min sobre D>0 diagonal de  ||D⁻¹ Bm D||₁  =  ρ(Bm).
```

Nenhuma escolha de pesos por DOF fura esse piso. Medido (raiz de Perron pelo
espectro denso; os dois quocientes de Collatz-Wielandt coincidem em 4 casas,
então o valor é firme):

| centro | ‖X‖ hoje | piso ρ(Bm) | ganho máximo de repesagem | falta para 0,6 |
|---|---:|---:|---:|---:|
| s=1/8+i/4 | 54,69 | **9,48** | 5,8× | 15,8× |
| s=1 | 58,41 | **12,27** | 4,8× | 20,5× |
| s=0,733 | 56,82 | **10,91** | 5,2× | 18,2× |

A melhor repesagem concebível tira 58 para ~10. O alvo é 0,6. Faltam 16–20×, e
não há para onde ir dentro da família.

## 4. Conclusão

As duas rotas locais estão fechadas:

- **pré-condicionador**: teto de 29–32×, melhor candidato real 1,93×, contra 57×;
- **cota analítica k=0**: piso de Perron 9,5–12,3, contra o alvo 0,6.

A decisão **INVIÁVEL SEM HPC** da §5 de `G_PARAMETRIX_FEASIBILITY.md` não só se
mantém como passa a ter, do lado do pré-condicionador, um argumento e não só
uma estimativa. Restam duas saídas, ambas fora desta máquina: uma alocação de
HPC (3,4–5,8·10⁴ h-núcleo, ver §7), ou uma arquitetura que não certifique a
casca coluna a coluna.

**Ressalvas honestas.** Tudo aqui é float na caixa reduzida 12×36; o piso de
Perron é um fato de álgebra linear, mas calculado em ponto flutuante e nessa
truncagem. O piso tende a CRESCER com a caixa (mais termos somados por coluna),
então a conclusão deve endurecer, não afrouxar, em 107×384 — mas isso não foi
medido. E o argumento cobre famílias de peso por DOF; um pré-condicionador
não diagonal e não coberto pelos quatro candidatos testados continua, a rigor,
fora do alcance desta nota.
