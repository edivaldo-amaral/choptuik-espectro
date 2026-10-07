# Auditoria de S4: revisão independente (25–26/09/2026)

Revisor com contexto limpo, de outro modelo da família (Fable, não o Opus que escreveu o
certificado), proibido de ler a narrativa do autor (`S4_GAUGE.md`, `RETOMADA.md`,
`PROOF_OBLIGATIONS.md`, auditorias e a revisão de S3b). A tarefa pedia um artefato por
alegação: `.codex-runs/2026-09-25-revisao-s4/TAREFA.md`. O relatório está em
`.codex-runs/2026-09-25-revisao-s4/RELATORIO.md`, com os scripts `rN_*.py` e as saídas
`rN_saida.txt`.

- **Snapshot antes:** `.snapshots/2026-09-25-revisao-s4/MANIFEST.sha256`, 242 arquivos. Conferido
  intacto depois.
- **Execução:** ~32 min. O revisor usou o PC do laboratório e o laboratorio2, em pastas separadas.

## Veredicto

**Nenhum erro numérico ou lógico no certificado de simplicidade e identificação.** Tudo se
reproduz bit a bit (`C3.json` e `E.json` regerados numa cópia), e as normas verdadeiras ficaram
abaixo das cotas:

| Grandeza | Verdadeira (ponto flutuante) | Cota |
|---|---|---|
| ‖M̂⁻¹‖ (etapa A, 14.017, laboratório) | 52,0404 | 52,14972 |
| a (Z'←T) | 0,39526 | 0,40461 |
| Y_Z sem margens | 5,1·10⁻⁷ | 3,0·10⁻⁶ |
| janela 32..47: ‖V‖ | 1,74117 | 1,74134 |
| janela 32..47: T←Z' | 0,38047 | 0,39412 |
| janela 32..47: ←48..63 | 0,25844 | 0,25919 |
| janela 32..47: ←16..31 | 0,26343 | 0,26538 |

O revisor também conferiu cinco pontos:

- **RefB por leitor próprio:**
  - ‖RefB‖_L = 2,350632·10⁻¹⁰;
  - ‖(Z_A+1)RefB‖_L = 2,454216·10⁻⁸;
  - Young ‖B_L(RefB)‖ <= 4,670165·10⁻⁸, com a tabela do artigo transcrita de forma independente.
- **Bola de RT:** o artigo centra a bola em RefA + RefB, com raio 2⁻²⁷⁷ na soma das quatro
  componentes, o que é mais forte do que o "por componente" do autor. ‖B_L(Corr)‖ <= 7·10⁻⁸³.
- **Gauge:** a identidade L_w(μ)(Z+1)w = (Z+2)F_w vale nos operadores do código com erro
  1,4·10⁻¹⁶, e os controles (sinal de ∂τ trocado, ou (Z+1)F) falham em 10⁻². A identidade das
  constraints também vale, o que dá D g* = 0.
- **Q.json:** recalculado e idêntico, com λ0 = μ_A. A junção dos checkpoints é exata.
- **Unicidade:** θ + 2 Z2 ρ <= 1 é conservador para todo x, porque o Z2 do código já é a
  constante de Lipschitz inteira.

Conferi por amostragem, antes de aceitar: o manifesto (intacto), as quatro normas da janela
32..47 (`r4_janela_saida.txt`), a periodicidade no código (`r6_saida.txt`) e a igualdade de
`C3.json` e `E.json` reproduzidos.

## Lacunas e respostas

| # | Gravidade | Lacuna | Resposta |
|---|---|---|---|
| L1 | média (enunciado) | Σ_A é invariante por s ↦ s + i: μ + ij são raízes de gauge, e "N_A = raízes em Re s > 0" sem faixa é infinito. | **Corrigido no texto** (`S4_GAUGE.md`, "A faixa fundamental"). A convenção do projeto já existia (`SECTOR_B.md`, Lema 5.1: L_A em −1/4 <= Im s < 3/4 = faixa de Γ_A ∪ faixa B), mas S4 não a enunciava. Em Q̃ só μ aparece, e N_A passa a ser explicitamente a faixa A. |
| L2 | baixa–média (rigor escrito) | ‖(Z*+1)Corr‖_L <= 10⁻⁶⁰ é afirmado sem argumento. A bola de RT controla Corr só em ℓ¹, e Z é derivada. | **Fechado** (`corr_gauge_bound.py`): O_{μ*}Corr = G pela definição de C em RT, ‖G‖ <= 5,5·10⁻⁸⁰ com as constantes publicadas; ‖∂τO⁻¹‖ <= 18 e ‖n Corr‖ <= 9‖(1+ξ)∂ξ Corr‖ pelas cotas de banda de RT, sem divisão por (1 + ξ). Resultado: ‖(Z*+1)Corr‖_L <= 1,1·10⁻⁷⁵. |
| L3 | baixa | O domínio forte de g* é citado de `RADIUS_RECOVERY.md`, escrito noutra realização. | **Fechado** (no mesmo lema): g* ∈ Y; (J + λ0)g_n ∈ Y coeficiente a coeficiente; J + λ0 é injetivo no domínio máximo (as soluções homogêneas são singulares na elipse), logo y_g = Q0⁻¹g_n ∈ Y. |
| L4 | cosmética | ‖B_*‖ omitia δμ‖SΓ1‖ ≈ 3,4·10⁻¹⁰. | **Corrigido** em `nk_L3_E.py`. δ muda em 10⁻¹⁶, e `E.json` foi refeito (dist = 1,8155·10⁻⁵). |
| L5 | baixa (rastreabilidade) | `nk_L3_T.py` lia `--lam0` da linha de comando (padrão 0,401). T3, Q e nk3_mu não gravavam λ0, e o 4,7·10⁻⁸ vinha só do ambiente. | **Corrigido:** `nk_L3_T` lê λ0 do `Z.json`, ou confere se foi dado, e grava `lam0` e `eps_fundo_L` no T3; `nk_L3_q` grava `lam0`; `nk_L3_C` confere que λ0, `rt` e `eps_mu` são os mesmos em Z, T3 e Q. Os `C3.json` de S3b e de S4 saem idênticos com os asserts novos. |
| L6 | hipótese (conceitual) | Faltavam H-reg, a mesma realização para N_A, rank(D\|E) nas outras raízes e em s_K, e a faixa. | **Incorporado** em `S4_GAUGE.md`: (iv) rank(D\|E_{s0}) = 0 onde K é injetivo, também ao longo de cadeias (derivando K C = M L), coberto por S3a e S3b, com o resto da faixa em F3 e C4; (v) mesma realização, a cargo de C3 e T1; (H-reg); completude do gauge só para germes analíticos, dentro de H-gauge. |

## O que o revisor não verificou

- A contagem N_A (C3) e a certificação de s_K (S3b, revisada à parte).
- A validação de L e C contra o exportador de RT.
- A reexecução das cotas rigorosas de Loewner/Rump: comparou com estimativas em ponto flutuante.
- NORMA_BL e a prova de RT.
- H-rec e H-gauge globais, que são hipóteses.

## Estado

A parte computacional de S4 (μ* simples, autoespaço span(g*), D g* = 0) está confirmada, e o lema de
L2/L3 está escrito (26/09): todas as lacunas da revisão estão respondidas. A
contagem "N_físico = N_A − 2 na faixa A" depende, além de S4, das obrigações abertas C3, C4, F3
e T1 e das hipóteses H-rec, H-gauge, H-cone e H-reg.
