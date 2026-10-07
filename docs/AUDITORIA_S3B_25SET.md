# Auditoria de S3b: revisão independente (25/09/2026)

Revisor com contexto limpo, de outro modelo da família (Fable, não o Opus que escreveu o
código), proibido de ler a narrativa do autor (`S3B_ROTA.md`, `RETOMADA.md`,
`PROOF_OBLIGATIONS.md`, auditorias). As alegações foram dadas de forma neutra, com um artefato
obrigatório por alegação: `.codex-runs/2026-09-25-revisao-s3b/TAREFA.md`. Relatório:
`.codex-runs/2026-09-25-revisao-s3b/RELATORIO.md` (scripts `rN_*.py`, saídas `rN_saida.txt`).
Snapshot antes: `.snapshots/2026-09-25-revisao-s3b/MANIFEST.sha256` (239 arquivos); conferido
intacto depois. O revisor rodou por ~35 min e usou o PC do laboratório (numa pasta separada).

## 1. Veredictos

| R | alegação | veredicto | o que fez |
|---|---|---|---|
| R1 | Rouché de K | OK com ressalva | derivou a N' do código; cobertura da circunferência em Arb, método próprio (sobreposição mínima 3,9·10⁻⁷ rad: rigorosa, justa); montou H e A no pior tile (caixa 40×200): todos os blocos abaixo de N'; Perron numérico 0,763 contra θ' = 0,8126 |
| R2 | contagem finita de K | OK | montador independente (base real, truncamento próprio): uma raiz simples em 0,4010247311, a próxima a 0,616 do centro |
| R3 | lema X_κ1 / Y | OK | prova própria por ponto fixo na cauda de Fourier; hipóteses: Young de B nos dois pesos e Re λ >= 0 |
| R4 | estrutura do NK | OK com ressalva | derivou o NK bordejado, a simplicidade algébrica, Z2 e todas as entradas de 3 níveis; achou L1–L5 |
| R5 | ataque numérico | OK | 19 normas verdadeiras (potência com LU, no laboratório), todas abaixo de `T3.json`; as mais justas: `‖V‖` da janela 32..47, 1,69163 contra 1,69180; `N[32..47 ← 16..31]`, 0,2508 contra 0,2525 |
| R6 | cotas verificadas | OK | Loewner/Rump contra Arb em 20 casos: nenhuma cota falsa (mal condicionados falham sem dar cota errada); Young de `fundo_L` acima dos blocos (razão máxima 0,58); 40 ε_ω seguro: pela norma ℓ¹ de RT, Young dá 5,01 ε_ω |
| R7 | testemunho | OK | c_ref com montagem própria 0,1717424; os c_J pequenos são reais (fundo ≤ 6·10⁻⁸ no modo d = 21); `‖Π_F C1‖` ≈ 0,96 < 1,5 usado |
| R8 | fechamento e espaços | OK com ressalva | D no retângulo de S3a; float(0,401) coberto pelos tiles; X_L com sinal NÃO está em X_K, então o lema A3 é indispensável; ressalva sobre cadeias de Jordan (§3) |
| R9 | reprodução | OK | C3.json e D.json reproduzidos campo a campo |

## 2. Lacunas e o que fiz com cada uma (conferido por mim)

| # | gravidade | o quê | conserto |
|---|---|---|---|
| L1 | baixa (formal) | em T←far, o far de Fourier (\|m\| >= 200, todo n) alcança pela banda \|d\| <= 40 as linhas BAIXAS das janelas com \|m\| >= 160: não é descida radial | `nk_L3_C.py`: nessas janelas usa `‖V_J‖` em vez de `‖V_J P_alto‖` (sem efeito: 1,089 < 1,281) |
| L2 | desprezível | linhas do far em 216 <= \|m\| < 240 recebem pelo resto da banda e não estavam em grupo | `nk_L3_C.py`: far←T += resto·`‖Q0 P‖`(janelas \|m\| >= 176)/min w |
| L3 | desprezível | componente far do resíduo e a parte de RT nas janelas não alcançadas | `nk_L3_C.py`: Y = max(Y_Z/v0, (Y_T + max w‖V‖ rt)/v_T, rt/v_f) |
| L4 | desprezível | a linha l Q0 P_far em Z'←far | `nk_L3_C.py`: epsZ += t_V·descida(NZ, Nc) (~10⁻²⁵) |
| L5 | baixa | a_extra supunha `‖Q0‖` decrescente em \|m\| | `nk_L3_q.py`: máximo CALCULADO nos 48 modos (0,17199, justamente o da borda) |

**A margem do NK (o achado mais importante).** O revisor mostrou que θ só podia subir 7,5·10⁻⁴
e Y·Z2 só 0,94%: o discriminante (1 − θ)² − 4 Z2 Y valia 2,4·10⁻⁴. A causa era Z2 = 116,8, com
uma cota frouxa: a linha de Z' usava `‖Q0‖` do modo inteiro (1,92) onde só entra
`‖P_Z Q0 P_T‖`, a descida das colunas de cauda (n >= 128) para Z. `nk_L3_q.py` calcula as
normas RESTRITAS (Rump, com o erro de Q0): **q_ZT = 0,169**, q_TT = 0,300, q_Zf ~ 10⁻²⁵.
Com elas e L1–L5:

| θ | Z2 | Y | r | \|λ* − 0,401024733\| | discriminante |
|---:|---:|---:|---:|---:|---:|
| 0,839936 | **42,96** | 5,44·10⁻⁵ | **3,78·10⁻⁴** | **<= 1,71·10⁻⁴** | **0,0163** (antes 2,4·10⁻⁴) |

θ agora pode subir até ~0,903, e **`NORMA_BL` pode ser até 23% maior (4,89 em vez de 3,964)**
antes de o NK falhar: a dependência crítica que o revisor apontou deixou de ser frágil.

## 3. A ressalva das cadeias de Jordan (R8)

O revisor observou que A7 fala de autovetores, e que uma raiz μ ≠ s_K de L em D poderia ter
autovetor que satisfaz as constraints e vetor generalizado g com C(μ)g ≠ 0. Isso já está
resolvido em `CONSTRAINT_GAUGE_EQUIVALENCE.md` §1: para μ ≠ s_K, K(μ) é injetivo e Fredholm de
índice 0 (K P = I + compacto), logo invertível, e K(s)⁻¹ é analítico perto de μ. Para uma
função-raiz h(s) = h0 + (s − μ)h1 (L(s)h(s) = O((s − μ)²)),
`C(s)h(s) = K(s)⁻¹ M(s) L(s) h(s) = O((s − μ)²)`: C(μ)h0 = 0 e C(μ)h1 + C'(μ)h0 = 0. A cadeia
satisfaz as constraints na ordem certa — é a condição para a solução e^{μt}(h1 + t h0) do
problema linearizado; C(μ)g = 0 sozinho não é a condição física.

## 4. O que ficou sem reconferência independente

- **`NORMA_BL` = 3,964043** (símbolo em Arb, `certify_symbol_numerical_range.py` e
  `symbol_slack_arb.py`). Mitigado: a margem agora tolera 23% de erro; uma cota inferior
  independente em ponto flutuante (`experiment_norma_BL.py`, compressão P B P) confere a
  consistência (§5).
- A norma verdadeira de Z'←T (a = 0,37889): só a montagem foi revisada.
- As formas fechadas de cauda para κ2 = 5/4 além da regressão contra as de K (auditadas em S3a);
  Qfar conferido por seção finita (0,1167 <= 0,1249).
- As validações `--validar` de `signed_operator_L` e `signed_constraints` contra o exportador de
  RT (feitas pelo autor, ≤ 1,8·10⁻¹⁵; o revisor não reexecutou).
- O μ de 80 dígitos (constante de RT; já tratado na revisão de S3a).

## 5. Checagem independente de `NORMA_BL` (25/09, tarde)

`scripts/experiment_norma_BL.py` (laboratório, 1,8 h; `build/s3b/norma_BL.txt`): cota INFERIOR
de ‖B_L‖ pela compressão P B_L P numa caixa |m| < 60, n < 200 (41.500 incógnitas; ‖P B P‖ <= ‖B‖),
por potência: **3,8544** após 80 iterações (3,818 → 3,844 → 3,851 → 3,854, ainda subindo
devagar). Fica 2,8% abaixo da cota do símbolo em Arb, 3,964043: consistente (a cota superior
está acima do que a compressão alcança e perto dela). Com a margem atual (até 4,89), um erro
da cota do símbolo teria de passar de 23% para derrubar o NK, e a compressão mostra que
‖B_L‖ >= 3,85.
