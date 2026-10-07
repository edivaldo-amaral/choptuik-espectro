# Revisão independente de S3b (25/09/2026)

Você é o revisor independente de um certificado assistido por computador. O autor é outro
modelo. Seu trabalho é **tentar quebrar** o certificado, não confirmá-lo. Um "OK" sem artefato
conferível não conta. Escreva em português do Brasil.

## Regras

- **Pode ler:** `scripts/*.py`, `build/s3b/**` (JSON, txt, npy), `build/s3b/rig/**`,
  `docs/S3_CONDICIONAMENTO_SHARP.md` §10.1–10.3 (definições do espaço com sinal e de S3a),
  `docs/CONSTRAINT_WITNESS_NORMS.md` §1–2 (definição de C(s)), `docs/SHARP_LINEARIZED_IDENTITY.md`
  (identidade K C = M L), os PDFs em `.codex-runs/2026-09-24-revisao-agente/` e o artigo de RT
  em `.cache/rt-1203.3766v1/`.
- **NÃO pode ler:** `docs/S3B_ROTA.md`, `docs/RETOMADA.md`, `docs/PROOF_OBLIGATIONS.md`,
  `docs/AUDITORIA_*.md`, `.codex-runs/2026-09-24-*/RELATORIO.md`, a memória do autor. São a
  narrativa e as conclusões do autor; as alegações estão enunciadas abaixo.
- **Não altere nenhum arquivo existente.** Escreva SÓ em `.codex-runs/2026-09-25-revisao-s3b/`
  (scripts `rN_*.py`, saídas `rN_saida.txt`, e o `RELATORIO.md`).
- Máquinas: o notebook local tem 7 GB de RAM (rode com `OPENBLAS_NUM_THREADS=2`,
  `.venv/bin/python`). Para contas pesadas (blocos densos de 11.200), pode usar o PC do
  laboratório: `ssh <usuario>@<host>` (chave já instalada), `~/Choptuik/.venv/bin/python`,
  15 GB de RAM, 4 núcleos; trabalhe SÓ numa pasta nova `~/revisao_s3b/` lá (copie o que
  precisar; não altere `~/Choptuik`). Rode coisas longas com `setsid nohup ... &` e acompanhe.
- Orçamento: até ~2 h de trabalho. Priorize o que pode estar errado.

## O que é alegado

Notação: K(s) e L(s) são os operadores de constraints e selecionado de RT, na base de Fourier
COM SINAL (`scripts/signed_operator.py`, `scripts/signed_operator_L.py`), com normas do toro
`||x||² = Σ κ1^{2m} ω2(n) |x_mn|²`, ω2(n) = κ2^{2n} + κ2^{-2n}: K em (129/128, 9/8), L em
(65/64, 5/4). Q0 = inverso da parte livre (diagonal em m, triangular superior em n). O disco é
D = {|s − 0,401| < 1/8}; S3a (já revisada, fora do escopo) certifica K injetivo na região
0 <= Re s <= 1,765, |Im s| <= 1/4, fora de |s − 0,401| < 0,05.

**A1 (K, Rouché).** Com P = Q0(c) fixo (c = 0,401), H(s) = K(s)P e a referência
A(s) = diag(Ĥ_ZZ(s), I, I, I), Ĥ_ZZ(s) = I + P_Z(B_t + s − c)Q0(c)P_Z (fundo truncado B_t,
μ diádico), vale ||A⁻¹(H − A)|| < 1 em |s − c| = 0,05, na norma max_i ||x_i||/v_i. Por
Gohberg–Sigal, H e A têm o mesmo número de zeros no disco. Os zeros de A são os autovalores
(com multiplicidade) da matriz finita K̂_Z = J_Z + P_Z B_t P_Z no disco. Código:
`scripts/rouche_K.py circ`, resultado `build/s3b/rouche_K_circ.json` (16 tiles, pior θ'
0,812558). A passagem da matriz de Perron N do tile para N' (linha Z dividida por 1 − e,
entrada (Z,Z) = ||V||(dB ||Q_ZZ|| + eps_mu)/(1 − e)) está em `circ()`.

**A2 (K, contagem finita).** K̂_Z (Z = 20×80) tem exatamente um autovalor no disco, simples.
Código: `rouche_K.py contagem`, `build/s3b/rouche_K_contagem_20x80.json` (bordejamento,
Taylor de t(s) = e^T M(s)⁻¹ e, Rouché escalar).

**A3 (lema).** Autovetores e cadeias de Jordan de K (e de L) no espaço com sinal X_κ1 estão
no espaço simétrico Y = X_κ1 ∩ X_{1/κ1}; raízes e multiplicidades coincidem.

**A4 (L, NK bordejado).** x = (y, λ), h = Q0 y, F(x) = (L(λ)h, ⟨h0w, h⟩ − 1), com
Q0 y0 = h0w (`build/s3b/rig/h0w.npy`, Z = 32×128). Aproximação do inverso
A = diag(M̂⁻¹, V_J, I) (M̂ bordejado de Z; V_J = Ĥ_J⁻¹ densos nas 26 janelas W = 16 da cauda
F \ Z, F = 200×400; far = complemento de F). Norma: max(||x_Z'||/v0, ||x_T||_W/vT, ||x_f||/vf),
||x_T||_W² = Σ w_J² ||x_J||². Etapas: A = `scripts/nk_L2_Z.py` → `build/s3b/rig/Z.json`;
B = `scripts/nk_L3_T.py` → `T3.json`; C = `scripts/nk_L3_C.py` → `C3.json`. Alegação:
θ = 0,839826 < 1, Z2 = 116,8, Y = 5,44·10⁻⁵, r = 6,19·10⁻⁴ satisfazem Y + θr + Z2 r² <= r e
θ + 2 Z2 r < 1; logo há um único zero x* na bola, DF(x*) é invertível e λ* é autovalor
algebricamente simples de L, com |λ* − 0,401024733| <= 2,81·10⁻⁴.

**A5 (cotas verificadas).** `scripts/verified_inverse.py` (||H⁻¹X|| por Loewner:
t²HH* − XX* ⪰ 0, Rump na realificação) e `scripts/verified_norms.py` dão cotas SUPERIORES
válidas; `scripts/fundo_L.py` dá a perturbação do fundo de L (arredondamento, resto da banda
|d| > 16, erro de RT = 40 ε_ω, μ).

**A6 (testemunho).** `scripts/nk_L3_D.py` → `build/s3b/rig/D.json`: ||Π_F C(λ*) h*|| >=
0,16849 > 0 (F = 12×48 no espaço de K). Obs.: as perdas das janelas de cauda saíram ~10⁻¹³ a
3,5·10⁻⁷ ("c_J"); o autor diz que isso é real (colunas de cauda em n >= 128, caixa em n < 48).

**A7 (fechamento).** Como K C = M L, um autovetor de L que viole as constraints só existe em
raiz de K; com A1–A6 e S3a, exatamente uma raiz de L em D viola as constraints: s_K, simples.

## Artefatos obrigatórios (um por alegação; script + saída + veredicto)

- **R1 (A1, estrutura):** derive você mesmo por que o Perron do tile implica
  ||A⁻¹(H − A)|| < 1 com a N' do código, e por que zeros de A = autovalores de K̂_Z (confira a
  triangularidade de K_livre nas caixas no código). Monte, em ponto flutuante e com SEU
  código, H(s) e A(s) numa caixa moderada num ponto da circunferência e estime
  ||A⁻¹(H − A)|| na norma dos pesos v do JSON. Confira a cobertura da circunferência pelos 16
  discos com código próprio (Arb ou racionais).
- **R2 (A2):** recalcule, com montagem própria, os autovalores de K̂_Z (20×80) no disco; confira
  a lógica da desigualdade do Rouché escalar no código.
- **R3 (A3):** escreva a prova (ou a refutação) do lema, com as hipóteses usadas.
- **R4 (A4, estrutura):** derive as condições de NK para F bordejado; confira que DF(x*)
  invertível implica simplicidade algébrica; confira a fórmula de Z2 em `nk_L3_C.py` contra
  DF(x) − DF(x0) e a estrutura de Q0 (quais P_l Q0 P_c são nulos). Confira as desigualdades da
  montagem de 3 níveis (T←T = ||W N W⁻¹||₂, T←far, far←T, Z'←T, T←Z') e a linha de Z'.
- **R5 (A4, ataque numérico):** com montagem PRÓPRIA (pode reusar `signed_operator_L` para B e
  J), calcule a norma verdadeira de alguns blocos de E = I − A DF(x0) (por exemplo, 2–3 pares
  de janelas vizinhas perto da borda de Z, a entrada T←Z' da janela 32..47, e ||V_J||) e
  verifique que as cotas de `T3.json` estão acima. Tente achar um termo esquecido (acoplamento
  não contado, banda, far, perturbação global, pesos).
- **R6 (A5):** teste `verified_inverse` contra valores exatos (Arb) em matrizes pequenas,
  incluindo casos mal condicionados; confira a contabilidade de erros. Confira `fundo_L`:
  as cotas de Young majoram as normas dos blocos montados; o 40 ε_ω tem justificativa?
- **R7 (A6):** recalcule c_ref com montagem própria de C (ou reuse `signed_constraints`, que é
  validado contra o exportador de RT); verifique diretamente os c_J pequenos para 2 janelas;
  confira as perdas.
- **R8 (A7 e espaços):** confira a lógica final, incluindo em que espaço cada afirmação vale
  (X_κ1, Y, L-space → K-space) e se a região de S3a + disco de A1 cobre D.
- **R9 (reprodução):** rode `nk_L3_C.py` e `nk_L3_D.py` com `--dir` apontando para uma CÓPIA de
  `build/s3b/rig` dentro da sua pasta e confira que reproduzem θ, r e a folga.

## Relatório

`.codex-runs/2026-09-25-revisao-s3b/RELATORIO.md`: tabela R1–R9 (veredicto: OK / OK com
ressalva / LACUNA / ERRO; artefato; o que fez), a lista de lacunas com gravidade e um conserto
proposto, e o que você não conseguiu verificar.
