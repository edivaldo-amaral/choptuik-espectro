# Checagem independente do NK de L em 0,7332 (09/10/2026)

Você é um revisor independente de um certificado assistido por computador. O autor é outro modelo. Seu
trabalho é **tentar quebrar** o certificado. Um "OK" sem artefato conferível (script mais saída, ou derivação
escrita passo a passo) não conta. Escreva o relatório em português do Brasil.

## O objeto

O certificado de Newton–Kantorovich (NK) bordejado de L em λ0 = 0,7331796596606924, em `build/nk733/`
(Z.json, Q.json, T3.json, C3.json, h0w.npy, os logs etapa*.txt e MANIFEST.sha256). Ele sustenta, no artigo
(`paper/sec/roots.tex`, Proposição 6.2, linha de s*), a afirmação:

> existe um único zero (y*, λ*) do sistema bordejado na bola B(x0, r), com r <= 2,16·10⁻⁵ e
> |λ* − λ0| <= 1,0163·10⁻⁵, e DF(x*) é invertível; logo λ* é raiz algebricamente simples de L, e
> s* ∈ [0,73316949; 0,73318983].

O mesmo código já foi checado por revisões anteriores nas raízes s_B (`build/faixaB/nk/`), s_K
(`build/s3b/rig/`) e μ (`build/s4/rig/`). Esta instância nunca foi checada à parte. Foque no que é próprio
dela: os dados de entrada, a montagem e os números.

## Regras

- **Pode ler:**
  - `scripts/*.py`, especialmente `nk_L2_Z.py`, `nk_L3_T.py`, `nk_L3_C.py`, `nk_L3_q.py`, `fundo_L.py`,
    `closed_form.py`, `verified_inverse.py` e `verifica_T2.py::check_nk`;
  - `build/nk733/**`, e, para comparação, `build/faixaB/nk/**`, `build/s3b/rig/**` e `build/s4/rig/**`;
  - os dados de RT em `.cache/rt-1203.3766v1/`;
  - `paper/sec/roots.tex` e `paper/sec/counting.tex`, só para o enunciado e o método.
- **NÃO pode ler:** `docs/C1_REAVALIACAO.md` e os outros documentos de narrativa ou de conclusões em `docs/`;
  nada em `.codex-runs/` além desta pasta; a memória do autor.
- **Não altere nenhum arquivo existente.** Escreva só nesta pasta e, no laboratorio2, só em `~/revisao_nk733/`.
- Ambiente: `.venv/bin/python` e `CHOPTUIK_EPS_FUNDO_L=4.7e-08`. No notebook (7 GB, 4 núcleos) use
  `OPENBLAS_NUM_THREADS=2`. Para cálculo pesado, use o **laboratorio2**
  (`ssh -o BatchMode=yes <usuario>@<host>`, projeto em `~/Choptuik`, 14 GB, 4 núcleos). Não use outras
  máquinas. Trabalhos longos: `setsid nohup ... &` e esperas com `until ...; do sleep 30; done` (os processos em
  segundo plano do harness morrem em ~30 min). Mate processos pelo PID.
- Snapshot: `.snapshots/2026-10-09-revisao-nk733/MANIFEST.sha256` (commit em `COMMIT.txt`). Confira-o no
  início e no fim.
- Grave resultados parciais no `RELATORIO.md`. Orçamento: até ~3 h.

## Artefatos obrigatórios

- **N1 (integridade).** Confira `build/nk733/MANIFEST.sha256`. Confira que Z.json, Q.json e T3.json foram
  produzidos para o mesmo λ0, as mesmas caixas (Z 32×128, F 200×400, janelas W = 16) e o mesmo h0w, e que o
  T3.json juntado tem as 26 janelas.
- **N2 (montagem).** Refaça, com código próprio em racionais exatos, a montagem de Perron de 3 níveis de
  C3.json a partir de Z.json, Q.json e T3.json: θ, Z2, Y, o discriminante, r e o erro de λ. Compare com o
  gravado e com o que `verifica_T2.py --so nk_0733` reproduz.
- **N3 (fórmulas).** Derive a forma de Z2 (sistema bordejado quadrático em (y, λ)), a cota de Y (resíduo do
  centro, com o erro do fundo e de μ) e a fórmula de r. Diga se o erro de λ é corretamente extraído de r e dos
  pesos v.
- **N4 (dados de entrada, por amostragem).** Escolha pelo menos:
  - 2 janelas de T3.json (uma central e uma periférica), e refaça no laboratorio2 as quantidades da janela
    (‖V_J‖, os acoplamentos, ρ_J, o resíduo) com `nk_L3_T.py` ou código próprio;
  - o nível Z (Z.json): refaça a cota de Loewner de ‖M̂⁻¹‖, ou um subconjunto conferível, e o resíduo Y_Z.

  Compare com o gravado.
- **N5 (enunciado do artigo).** Confira que a linha de s* da tabela da Proposição 6.2 e o intervalo
  [0,73316949; 0,73318983] do Teorema Principal saem dos números certificados e estão arredondados para o lado
  seguro.
- **N6 (síntese).** Uma lista de achados com gravidade (erro, lacuna, imprecisão, cosmético), local, evidência e
  correção sugerida, e um veredito: o certificado sustenta a afirmação?

## Entrega

`RELATORIO.md` nesta pasta, com N1–N6, mais os scripts e saídas `nN_*`, e a conferência do manifesto no fim.
