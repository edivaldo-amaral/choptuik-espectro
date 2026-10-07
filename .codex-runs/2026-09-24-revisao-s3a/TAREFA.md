# Revisão independente: certificado rigoroso de S3a por tiles de Perron

Você é um revisor matemático independente. Outro modelo escreveu e auditou (ele mesmo) um
certificado assistido por computador. Seu trabalho é procurar erros. Escreva em português.

## Contexto em 5 linhas

Projeto: contagem rigorosa de modos instáveis em torno da solução crítica de Choptuik
(continuação de Reiterer–Trubowitz). A obrigação S3a diz: o operador de constraints `K(s)`
é injetivo em `0 <= Re s <= 1,765`, `|Im s| <= 1/4`, fora do disco aberto `|s − 0,401| < 0,05`.
Isso foi "provado" cobrindo a região com 154 discos (tiles), cada um certificado por um
critério de Perron sobre cotas de blocos de `E = I − A·K(s)·P`. Leia primeiro
`docs/S3_CONDICIONAMENTO_SHARP.md` §10.4 a §10.12 e `docs/AUDITORIA_S3A_TILES_24SET.md`.

## Arquivos a revisar (SÓ LEITURA)

- `scripts/tile_perron.py` — o certificado de um tile (funções `Ctx`, `parte_Z`, `parte_T1`,
  `parte_T2`, `montar`, `tile`, `Rw_cert`, `descida`, `verifica_perron`, `young_abs`,
  `resto_banda`, `truncar`).
- `scripts/verified_norms.py` — cotas superiores de normas em ponto flutuante (critério de
  Rump). O artigo está em `.codex-runs/2026-09-24-revisao-s3a/Rump2006_verification_positive_definiteness.pdf`
  (Teorema 2.3, cota I da §3, Corolário 2.4).
- `scripts/check_cover.py` — reverifica cada tile em racionais e checa a cobertura.
- `scripts/tile_cover.py` — o planejador (quadtree em três escalas).
- `scripts/signed_operator.py` e, em `scripts/fourier_tail_bound.py`, as funções
  `modos_campos, modo, refl, combos, mult, omega, B_rad, norma_Bd_young, livre_inv,
  cauda_Rinv, Rinv_uniforme`.
- Resultados: `build/cobertura_lab/grande-*.json`.

## As 7 perguntas (responda TODAS, uma seção por pergunta)

1. **Estrutura (A2).** `K_livre` é mesmo bloco-triangular superior na decomposição
   `[Z, T1, T2, far]`? O operador `P = Q0(sZ)P_Z + Q0(s1)P_T1 + Q0(s2)P_(T2+far)` é bijeção sobre
   o domínio? A identidade `K(s)P = I + Σ_ℓ (B + s − s_ℓ) Q0(s_ℓ) P_ℓ` está certa?
2. **Lema de Perron (A3).** Se `N` (matriz 4×4 de cotas dos blocos de `E`) tem `v > 0` com
   `N v <= θ v`, `θ < 1`, então `I − E` é injetivo? Confira as hipóteses exatas usadas.
3. **Completude e correção dos blocos (A4).** Em `montar()` e `parte_*`, cada entrada
   `N[i][j]` majora `||P_i E P_j||` para TODO `s` do tile, com os três centros e os termos
   `(s − s_ℓ)`? Atenção especial: entradas envolvendo `far`; onde entra `dB` (perturbação
   global do fundo: resto da banda truncada, erro de RT, arredondamento); as linhas de Z
   (multiplicadas por `V`); os Grams acumulados.
4. **Critério de Rump (D1).** `verified_norms.py` satisfaz as hipóteses do Corolário 2.4?
   Simetria exata da matriz fatorada, `ã_ii <= a_ii − c`, `c >= ||Δ(A)||`, realificação de
   hermitianas, o erro de formar `X*X` e Grams em ponto flutuante, o uso do LAPACK `dpotrf`.
5. **Contabilidade de erros (D3).** `err_mm`, `err_gram`, `Rw_cert` (inverso por resíduo),
   `FATOR_PESO` (pesos calculados), o sanduíche `V S V*`, `cota_norma(X, erro_X)`: cada termo
   de erro majora o que diz majorar?
6. **Cotas da forma fechada (C2).** Rederive `cauda_Rinv`, `Rinv_uniforme` e `descida`
   (em `fourier_tail_bound.py` e `tile_perron.py`) a partir da forma fechada de
   `(J + σ)⁻¹`, com `J` triangular superior. Estão corretas? As margens bastam?
7. **Cobertura (E1).** `check_cover.py` prova mesmo que a união dos discos reverificados
   cobre a região menos o disco ABERTO D? O argumento de conjugação (§10.6) está certo? A
   região `0 <= Re s <= 1,765`, `0 <= Im s <= 1/4` é a pedida?

## Regras

- **NÃO edite nenhum arquivo existente.** Só pode criar arquivos DENTRO de
  `.codex-runs/2026-09-24-revisao-s3a/`.
- Para cada pergunta, dê um veredicto: **OK**, **LACUNA** (argumento incompleto, mas
  provavelmente conserta) ou **ERRO** (algo falso).
- **Todo LACUNA ou ERRO precisa de evidência verificável:** a linha exata do código e um
  contraexemplo matemático explícito, OU um script em `.codex-runs/2026-09-24-revisao-s3a/`
  que, rodado com `.venv/bin/python`, exibe o problema (por exemplo, uma cota menor que a
  norma verdadeira numa caixa pequena). Achado sem evidência não conta.
- Pode rodar scripts pequenos para testar (caixas como `--Z 8x24 --Gp 12x36 --F 32x77 --D 4`;
  os testes existentes `scripts/test_tile_perron_blocks.py --rigoroso`,
  `scripts/test_far_bounds.py`, `scripts/test_verified_norms_arb.py`). NÃO rode a cobertura
  inteira nem caixas grandes (`Z 20x80` e acima): leva horas e usa gigabytes.
- Use `.venv/bin/python`. Seja específico e curto; não reescreva o documento do autor.

## Saída

`.codex-runs/2026-09-24-revisao-s3a/RELATORIO.md`, com:
1. uma tabela `pergunta | veredicto | evidência (arquivo:linha ou script)`;
2. uma seção por pergunta, com o raciocínio;
3. uma lista final dos ERROS e LACUNAS, em ordem de gravidade.
