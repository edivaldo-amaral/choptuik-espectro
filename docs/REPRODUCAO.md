# Reprodução dos certificados de T2

Guia para conferir o teorema T2 (`T2_ENUNCIADO.md`) a partir deste repositório. O mapa completo, com
cada alegação, seus arquivos, sha256 e a máquina em que estão, é gerado por
`scripts/inventario_certificados.py` em `docs/REPRODUCAO_MAPA.md` e `build/reproducao/inventario.json`.

## Dois níveis

| nível | o que refaz | custo | o que prova |
|---|---|---|---|
| **1** | as verificações finais, em aritmética racional, a partir dos certificados gravados (JSONs pequenos) | ~2 h em 4 núcleos (discos) + < 1 min (o resto) | que os números gravados implicam T2 com o código do repositório |
| **2** | os certificados intermediários, a partir dos dados de RT: matriz de referência, cotas de Loewner, fatores, caudas, janelas, Arb dos símbolos, tiles | dias, em 3–4 máquinas | que os números gravados vêm dos dados de RT |

O nível 1 é o que um leitor deve rodar primeiro. O nível 2 é a reprodução completa.

## Ambiente

- Python 3.14, com `numpy==2.5.3`, `scipy==1.18.1` e `python-flint==0.9.0` (Arb), conforme
  `requirements.txt`. O sympy só é usado em conferências simbólicas auxiliares e nas revisões.
- **Variável obrigatória:** `CHOPTUIK_EPS_FUNDO_L=4.7e-08`. É o erro do fundo de RT em L,
  ‖B_L(RefB)‖ <= 4,67·10⁻⁸ (`S4_GAUGE.md`). Todos os certificados de L foram produzidos com ela.
  `scripts/verifica_T2.py` a define sozinho. Sem ela, `fundo_L.py` usa a cota genérica 40ε_ω, maior.
- Os θ racionais dos discos são construídos a partir de cotas em ponto flutuante, e o arredondamento das
  entradas depende do BLAS e do número de threads.
  - Por isso a reprodução de um disco em outra máquina dá o mesmo θ a ~10⁻¹⁶, não a igualdade racional.
  - O critério do nível 1 é θ < 1 em todos os pontos, com diferença na escala do arredondamento. Uma
    diferença maior indica uma entrada diferente.

## Dados de RT

`bash scripts/baixa_dados_rt.sh` baixa o fonte do arXiv:1203.3766v1, confere o sha256 do tarball e de
`RefA.dat` e `RefAplusB.dat`, e extrai em `.cache/rt-1203.3766v1/`. O repositório não redistribui esses
arquivos. O artigo é Reiterer & Trubowitz, *Choptuik's critical spacetime exists*, Comm. Math. Phys. 368
(2019).

## Nível 1: `scripts/verifica_T2.py`

```
.venv/bin/python scripts/verifica_T2.py              # tudo (os discos levam ~6 h num processo)
.venv/bin/python scripts/verifica_T2.py --so s3a,nk_s4
.venv/bin/python scripts/verifica_T2.py --disco build/faixaB/discos/disco_0_5.json   # um disco (para paralelizar)
```

| verificação | passo de T2 | o que refaz | tempo |
|---|---|---|---|
| `discos_A`, `discos_B` | 6, 7 | o θ racional de **cada ponto** de cada disco (Perron de 3 níveis), com `rouche_L_disco.contexto`/`theta_disco`, no (p, r) gravado: θ < 1, e a diferença para o gravado na escala do arredondamento | ~100 s por ponto + ~150 s por centro; 22 discos, ~175 pontos |
| `cobertura_A`, `cobertura_B` | 6, 7 | a união dos discos cobre ∂Ω_A e ∂Ω_B (racionais) | < 1 s |
| `contagem_A`, `contagem_B` | 6, 7 | t·ε_B < 1 no contorno e a contagem de Schur: 2 zeros na faixa A (s_K e 0,7332; com μ*, N_A = 3) e 1 na faixa B (s_B) | < 1 s |
| `nk_s3b`, `nk_s4`, `nk_raizB`, `nk_0733` | 9, 8, 10, 11a | o Perron de 3 níveis do NK bordejado, em racionais, a partir de Z.json, Q.json e T3.json; θ, r e o erro de λ idênticos aos gravados | < 1 s cada |
| `s3a` | 11 | o Perron de cada um dos 154 tiles, em racionais, e a cobertura da região de S3a (`check_cover.py`) | ~5 s |
| `rouche_K` | 9, 11 | o Rouché de K no disco de S3b: os 16 tiles gravados como certificados (o recálculo é nível 2) | — |
| `simbolo_L`, `simbolo_K` | 5, 11 | F3: recombina as faixas em Arb gravadas (R_fundo = 2,9492, ‖B_L‖ = 3,9641; com a parte livre −0,0232, R_L <= 2,926); R♯ (1,4396) gravado | < 1 s |
| `L_real_toro`, `radius_transfer` | 4 | a contração da cauda no toro (0,371) e em RT (0,267) | minutos |
| `F1`, `R5` | 11b, 11c | ω4 ≢ 0 no centro e no cone; a recuperação de raio de qualquer peso | ~30 s |

**Resultado (03–07/10/2026): nível 1 completo, TUDO OK.** Nenhuma verificação falhou. Os relatórios estão em
`build/reproducao/`.

| verificação | resultado | onde rodou |
|---|---|---|
| discos A e B | 22 discos, **175 pontos**, θ < 1 em todos (pior 0,999074); diferença para o gravado <= 7,8·10⁻¹⁶ | laboratorio2 (4 processos) |
| coberturas | ∂Ω_A coberto por 96 discos; ∂Ω_B por 174 | notebook |
| contagens | A: 2 zeros (s_K e 0,7332), mais μ* = 3; B: 1 zero (s_B); t·ε_B <= 0,0334 | notebook |
| janelas A e B | 20 regiões, max ρ_J <= 0,152, iguais ao gravado | laboratorio2, PC 3 |
| NK (4) | θ = 0,839936 (S3b), 0,880172 (S4), 0,911526 (raiz B), 0,805721 (0,7332), **idênticos** | notebook |
| testemunhos D e identificação E | c_ref = 0,17174 (S3b) e 0,35519 (raiz B); dist = 1,82·10⁻⁵ (S4), **idênticos** | notebook |
| S3a | 154 tiles reverificados em racionais (pior θ 0,9961), cobertura completa | notebook |
| símbolos, Rouché de K | R_fundo,L = 2,9492, ‖B_L‖ = 3,9641 (recombinados); R♯ = 1,4396; 16 tiles de K | notebook |
| F1/F1′, R5, transferência de raio | idênticos | notebook |
| L-real no toro | contração 0,371 da cauda no toro, cauda_max, qT e defeito idênticos | laboratorio2 |

**Nível 2 (amostras), 07/10/2026:**
- **Matriz de referência A,** reconstruída dos dados de RT com os parâmetros do `schur.json`: **idêntica
  bit a bit** (sha256 `70791ce0…`, diferença máxima 0). Os dados da deflação também batem exatamente
  (`build/reproducao/nivel2_A.json`, no PC do laboratório).
- **Resíduo de Schur** (ε_R, ε_G, ε_B e os zeros em Ω), refeito a partir de A, T e U: **idêntico**
  (`build/reproducao/nivel2_contagem.json`, no laboratorio2). T e U são dados; o certificado é o
  resíduo.
- **Cauda T3 do centro 3,0,** refeita do zero com `nk_L3_T.py --rouche`: em curso no PC do laboratório.

**O que só o nível 2 recalcula:**
- os símbolos em Arb (as faixas);
- os tiles do Rouché de K;
- os níveis Z dos discos (Loewner);
- os fatores F e as caudas T3;
- os NK nas etapas A e B.

O nível 1 confere os certificados que esses cálculos produziram.

## Nível 2: a cadeia completa

Cada certificado tem o seu pipeline documentado. Os mais longos são:

| certificado | pipeline | onde está descrito | custo (aprox.) |
|---|---|---|---|
| matriz de referência A (Z = 32×128) e Schur (T, U) | `rouche_L_Z.py`; `rouche_L_contagem.py residuo` | `C1_REAVALIACAO.md` §8–9 | horas; A, T e U têm 3,1 GB cada |
| nível Z dos discos (t_p por Loewner pré-condicionada) | `rouche_L_rig.py` (versão 2) | §9–10 | horas por centro |
| fator de posto baixo (certF) | `rouche_L_Z.py certF` | §9 | ~1 h por centro; F tem 594 MB |
| caudas T3 (janelas fortes) | `nk_L3_T.py --rouche` | §9, §11 | 3–6 h por centro |
| janelas | `rouche_L_janelas.py` | §9 | horas |
| NK (S3b, S4, raiz B, 0,7332): etapas A, tQ, q, B, C, D, E | `nk_L2_Z.py`, `nk_L3*.py` | `S3B_ROTA.md`, `S4_GAUGE.md`, §11–12 | 3–6 h cada |
| S3a (154 tiles) | `tile_certificate.py`, `tile_cover.py` | `S3_CONDICIONAMENTO_SHARP.md` §10 | ~1 dia em 3 máquinas |
| símbolos (F3, R♯) em Arb | `certify_symbol_numerical_range.py` | `ALTERNATIVE_ROUTES.md` §5.25–5.26 | horas em 4 faixas |

**Dados pesados** (no Zenodo, DOI reservado **10.5281/zenodo.23222363**, rascunho até o artigo; manifesto em
`build/reproducao/pacote_dados.json`):
- A.npy (sha256 `70791ce0…`);
- os 13 fatores F_*.npy;
- T.npy e U.npy.

Hoje estão no notebook e no laboratorio2 (ver o mapa).
