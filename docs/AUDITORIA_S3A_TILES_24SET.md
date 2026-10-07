# Auditoria do modo rigoroso dos tiles de S3a (24/09/2026)

Autor: Claude, **o mesmo que escreveu o código auditado**. Isso é uma limitação: a
auditoria de 16/09 foi feita por um sub-agente separado. Ver "O que esta auditoria não
é", no fim.

Escopo: o certificado de injetividade de `K(s)` por tiles (`scripts/tile_perron.py
--rigoroso`), as normas verificadas (`scripts/verified_norms.py`), o planejador
(`scripts/tile_cover.py`) e o verificador de cobertura (`scripts/check_cover.py`).
Seções 10.4–10.10 de `S3_CONDICIONAMENTO_SHARP.md`.

Regra: **CONFIRMADO** exige rederivação escrita aqui ou recomputo independente (que não
importa o código auditado, salvo onde indicado). **PARCIAL** quando há uma etapa aceita
por argumento, não por construção. Ponto flutuante só como proposta, ou como teste
numérico rotulado.

## 1. Veredictos

| # | Alegação | Como foi checada | Veredicto |
|---|---|---|---|
| A1 | Espaço com sinal `X_κ1` (`Σ κ1^{2m} ω2(n)|x_mn|²`); `Im s >= 0` basta | Rederivação: o espaço simétrico `Y` é invariante pela conjugação e está contido em `X_κ1` (§10.6) | CONFIRMADO |
| A2 | `K_livre` é bloco-triangular em `[Z, T1, T2, far]`; `P` dos três centros é bijeção; `H(s) = I + Σ (B + s − s_ℓ) Q0(s_ℓ) P_ℓ` | Rederivação: caixas preservadas (diagonal em `m`, triangular em `n`); diagonais `(K_ℓℓ + s_ℓ)⁻¹` | CONFIRMADO |
| A3 | Lema de Perron: `N v <= θ v`, `θ < 1` ⇒ `E` contrai na norma `max ||P_i x||/v_i` ⇒ `H` injetivo | Rederivação | CONFIRMADO |
| A4 | Lista de blocos completa, e cada cota majora o bloco | Releitura linha a linha (§10.6) e teste com a **verdade** na caixa `F`, fundo completo, três configurações (`test_tile_perron_blocks.py`) | CONFIRMADO (pior razão 0,9999) |
| B1 | Operador com sinal = montador real | `K_sinal = T K_real T⁻¹` **entrada a entrada**, com `T` explícito (`test_signed_similarity.py`) | CONFIRMADO (`<= 3,6·10⁻¹⁵` em 8×24, 12×36, 16×50) |
| B2 | Montador real = exportador de RT | Auditoria de 16/09 (1,8·10⁻¹⁵) | CONFIRMADO antes |
| C1 | `||B♯|| <= 1,64967` na norma do certificado | Leitura do script e do JSON: toro com sinal, raios (129/128, 9/8), `μ ∈ [1/6; 0,17]` por convexidade, `μR_x` incluído; conferido contra 1,507 numa caixa 24×60 | CONFIRMADO: centros em Arb e, desde 24/09, também a folga de Lipschitz em Arb (`symbol_slack_arb.py`): 0,053616 contra 0,053670 gravado; cota refeita `||B♯|| <= 1,649617` (usada: 1,64967) |
| C2 | `cauda_Rinv`, `Rinv_uniforme`, `descida` | Rederivação (forma fechada de `R_m⁻¹`, Schur) e verdade numa truncagem de 1500 (`test_far_bounds.py`) | CONFIRMADO (pior razão 0,965), com margem 1,0001 a priori |
| C3 | Young dos campos (`norma_Bd_young` para `d ≠ 0`, `young_abs`, `resto_banda`) | Rederivação: estrutura 2×3 com as metades par/ímpar de `c2` ortogonais; ℓ¹ pesado majora o sup no toro | CONFIRMADO |
| C3' | Termo de `μR_x` em `norma_Bd_young` (`d = 0`) | Verdade 1,4253 contra 1,4257 da fórmula; Schur rigoroso dá 1,5226 | **NÃO É COTA PROVADA**; retirado do certificado (ver §2) |
| C4 | Erro do fundo de RT: `40·2⁻²⁵` em `dB` | Estrutura de `Γ2♯` pede ~5–10ε; o símbolo já soma 20ε; o tile soma outros 40ε | CONFIRMADO (folga dupla) |
| D1 | Critério de Rump | **Fonte primária** (BIT 46:433–452, Teorema 2.3, cota I da §3, Corolário 2.4) | CONFIRMADO após correções (§2) |
| D2 | As normas verificadas majoram a verdade | Cota independente **em bolas** (Arb, por congruência com base de autovetores), 5 matrizes, incluindo blocos reais do certificado e condicionamento 10¹² (`test_verified_norms_arb.py`) | CONFIRMADO (pior Arb/Rump = 0,999999) |
| D3 | Contabilidade de erros (pesos, entradas do fundo, resíduo de `R_m⁻¹`, produtos, Grams, sanduíche com `V`) | Releitura com majorantes entrada a entrada; cada erro de produto vem de `|fl(AB) − AB| <= γ |A||B|` | CONFIRMADO após correções (§2) |
| D4 | Verificação final de Perron | Racionais exatos no tile e, de novo, no `check_cover` a partir do JSON | CONFIRMADO |
| E1 | A união dos discos cobre a região menos D | `check_cover.py`: quadtree própria em racionais, sem usar o planejador | CONFIRMADO após correção (§2); faixa `0,25 <= Re s <= 0,75` (quadrados 1 e 2, incluindo a volta do disco D) completa, 46 discos. A quadtree do verificador precisa descer até a folga dos raios (`~10⁻⁶ r`) nos vértices da malha: profundidade ~26 |

## 2. O que a auditoria achou

1. **Afirmação errada minha sobre Rump.** Ontem à noite escrevi que o critério "é um
   teorema para matrizes reais simétricas" e troquei o Cholesky complexo pela
   realificação, chamando a versão anterior de "argumento, não teorema". A fonte diz
   outra coisa: o Teorema 2.3 cobre `A* = A ∈ M_n(F + iF)` com as mesmas constantes.
   As duas versões eram válidas; a realificação é só mais conservadora (`n → 2n`).
   Fica retificado aqui e na §10.7.
2. **Hipótese de simetria exata.** O teorema pede `A` exatamente simétrica
   (hermitiana). Um Gram calculado em ponto flutuante não sai assim, e o LAPACK lê só
   um triângulo; a simetrização implícita pode aumentar a norma de um erro por até
   `√2` em Frobenius. Antes isso estava coberto apenas pela folga do `γ`
   (`γ_{4(k+1)} >= 2 γ_{k+2}`). Agora usa-se a parte hermitiana `(G + G*)/2`, que em
   ponto flutuante é exatamente simétrica, com o arredondamento dela na margem.
3. **`c` precisa majorar a cota.** A soma `tr(A)`, com ~1,4·10⁴ termos, arredonda
   mais do que o fator `1 + 8u` que eu aplicava. Antes a sobra vinha do termo
   `2u·max a_ii` (o necessário é `u·max a_ii`). Agora `c` é inflado por `1 + 10⁻⁶`.
4. **Somas de cotas em ponto flutuante.** Cada entrada de `N` é soma de poucas dezenas
   de cotas positivas e pode arredondar para baixo em ~10⁻¹⁴ relativo. Agora é
   inflada por `1 + 10⁻⁹`, na montagem e na reverificação do `check_cover`.
5. **Termo de `μR_x` em `norma_Bd_young`.** Não é cota provada (acima). Era usado só
   como `||B_0||` no erro das composições (`× dq ~ 10⁻¹²`). Trocado por
   `||B_d|| <= ||B||`: a extração do modo `d` é média de conjugações unitárias por
   rotação em θ.
6. **Vértices da malha de tiles** (achado pelo verificador; §10.9): raios inflados por
   `1 + 10⁻⁶`.

**Efeito sobre a cobertura em curso.** Os quadrados que já rodavam usam a versão
anterior às correções 2–5. Eles continuam válidos:

- a correção 2 estava coberta pela folga do `γ` (fator ≥ 1,89 > √2);
- a 3, pela sobra do termo `2u·max a_ii`;
- a 4 e a 5, pela inflação de `10⁻⁹` que o `check_cover` aplica ao reverificar
  (≥ 3·10⁻¹⁰ absoluto, contra os ~10⁻¹² em jogo).

Os quadrados que ainda estão na fila rodam com a versão corrigida.

## 3. Pendências antes de chamar de prova

1. ~~Folga de Lipschitz do símbolo em Arb~~ — **feita** (C1): a gravada cobre a de Arb, para K e para L.
2. **Raio do disco D de S3b** fixado junto com S3b (hoje 0,05, provisório; §10.10).
3. ~~Cobertura completa~~ — **feita** (§10.12 de `S3_CONDICIONAMENTO_SHARP.md`): 154 tiles, `check_cover` COMPLETA.
4. **Auditoria independente** (abaixo).

## 4. O que esta auditoria não é

É uma releitura crítica feita pelo próprio autor, com checagens independentes
(testes que comparam contra a verdade numérica e contra aritmética de bolas, a fonte
primária de Rump, e um verificador de cobertura que não usa o planejador). Ela pegou
cinco problemas. Não substitui um revisor sem o contexto de quem escreveu, que é o
padrão das auditorias de 16/09. Recomendo uma rodada assim quando a cobertura fechar,
com escopo: A2–A4, D1–D3 e o `check_cover`.

## Arquivos

`scripts/test_tile_perron_blocks.py`, `scripts/test_far_bounds.py`,
`scripts/test_signed_similarity.py`, `scripts/test_verified_norms_arb.py`,
`scripts/check_cover.py`; correções em `scripts/verified_norms.py` e
`scripts/tile_perron.py` (commits de 24/09, madrugada).

## 5. Revisão pelo Codex (24/09, tarde)

`.codex-runs/2026-09-24-revisao-s3a/` (`gpt-5.6-luna`, esforço médio, 131 mil tokens, 8 min).
Protocolo da memória: snapshot com manifesto de 185 arquivos (conferiu intacto na volta),
lista fechada de 7 perguntas, e todo achado exigindo evidência executável.

**Resultado: os 7 itens OK, nenhum ERRO ou LACUNA.** Minha avaliação: **fraca como evidência
independente.**

- O revisor rodou os três testes existentes e um teste próprio de `cota_norma` (600 matrizes
  aleatórias, pior razão 1,0010, nenhuma cota abaixo da verdade).
- As perguntas 1, 3 e 5 foram respondidas descrevendo o código, sem conta ou teste
  independente.
- Na pergunta 6, a forma fechada que ele escreve está **errada**:
  `(1/d_j) Π_{k=j}^{c−1}(−v_k/d_k)`, em vez de `−(1/d_j)(v_c/d_c) Π_{k=j+1}^{c−1}(1 − v_k/d_k)`.
  Mesmo assim, ele afirma que o código "implementa exatamente esse produto". O OK dessa
  pergunta, portanto, não vem de uma rederivação.

Conclusão: não achou nada, mas também não mostrou ter procurado a fundo. Uma segunda rodada
precisa exigir **artefatos conferíveis por pergunta**: a forma fechada checada contra
`np.linalg.inv` num exemplo; a enumeração dos 16 pares de blocos `(i, j)`, cada um casado com
a linha do código; e uma tentativa ativa de violar cada cota numa caixa pequena.

## 6. Revisão independente por agente (24/09, tarde)

`.codex-runs/2026-09-24-revisao-agente/RELATORIO.md`. Revisor com contexto limpo, de outro
modelo da família (Fable 5.1, não o Opus 5.5 que escreveu o código), proibido de ler esta
auditoria e o relatório do Codex. A tarefa exigia artefatos conferíveis. O manifesto (187
arquivos) conferiu intacto. O revisor rodou por ~1 h.

| artefato | veredicto do revisor | o que fez |
|---|---|---|
| R1 forma fechada | OK | fórmula própria contra `inv` (2,5·10⁻¹⁵); `cauda_Rinv`, `Rinv_uniforme`, `descida` contra a verdade (pior 0,9714) |
| R2 16 blocos | OK com lacuna | cada par casado com a linha do código; faltava um termo em `(Z, far)` e `(T1, far)` |
| R3 ataque | OK | montagem **própria** com fundo completo, 7 configurações de centros, borda do tile: pior 0,9999998 |
| R4 Rump | OK | verdade em Arb; nenhuma cota abaixo dela |
| R5 cobertura | OK | 1,2·10⁵ pontos racionais, incluindo ∂D e vértices: 0 descobertos |
| R6 estrutura | OK | deu o argumento de domínio maximal de `P`, que a §10.6 omitia |

**Lacunas que ele achou, e o que fiz com cada uma (conferido por mim):**

1. **L1 — `μR_x` desce de qualquer `n`.** Real: no montador, a coluna `c1` de `n` ímpar recebe
   `2μ(−1)^{k−j}` em todas as linhas `2j < n`. Meu comentário em `montar` supunha que o fundo
   truncado só desce `BAND_N`. Deduzi a cota `descida_rx(J, Nc)` (composição `μR_x·Q0`, Schur)
   e a testei contra a verdade (folgada ~60–90×). Ela vale `1,3·10⁻¹⁴` (Z) e `1,6·10⁻¹⁰`
   (T1). Está somada em `montar` e, para os tiles já calculados, na reverificação do
   `check_cover`.
2. **L2 — μ fixo.** Real, e mais séria do que o revisor pôde ver. RT provam
   `|μ − 0,168307078963449969510…| <= 10⁻⁸⁰` (80 dígitos); o `μ_RefA` diádico difere disso
   em **3,88·10⁻¹¹**. Os tiles certificavam `K(μ_RefA)`. Conserto:
   `K(μ) = K(μ_RefA) + δμ(J₁ + R_x)`, com `J₁ Q0 = (I − σ_m R_m⁻¹)/μ` modo a modo, e cotas da
   forma fechada. Isso dá `ε_μ(s) ~ 5·10⁻⁸–1,1·10⁻⁷` por coluna (vezes `||V||` na linha de Z).
   Está em `eps_mu()`, somado em `tile()` e na reverificação.
3. **L3 — versões mistas.** A análise do revisor (as versões anteriores são válidas a menos de
   ~10⁻⁷ relativo, coberto pela folga de 0,39%) confere com a minha (§2).

**Com L1 e L2 somadas aos 154 tiles:** pior θ = **0,9961189** (antes 0,9961146); `check_cover`
dá COBERTURA COMPLETA.

Notas menores aceitas: a §10.12 só vale para `Im s <= 0` no espaço simétrico `Y` (corrigido
lá); `cota_gram` prova `t(1 + u/2)`, absorvido pelas inflações; FMA e FTZ ficam fora do
modelo literal de Rump e são irrelevantes (o numpy não liga FTZ).

**Leitura do processo.** A revisão pelo Codex (Luna, esforço médio, perguntas abertas) deu
"7 OK" em 8 minutos e não achou nada. Este revisor, com tarefa de artefatos obrigatórios,
achou duas lacunas reais, uma delas (μ) de princípio. O que fez a diferença foi o formato da
tarefa, não só o modelo.
