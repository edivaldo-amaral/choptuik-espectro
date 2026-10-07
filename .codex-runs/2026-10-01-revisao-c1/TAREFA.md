# Revisão independente do Rouché de L na faixa A (01/10/2026)

Você é o revisor independente de um certificado assistido por computador. O autor é outro
modelo. Seu trabalho é **tentar quebrar** o certificado, não confirmá-lo. Um "OK" sem artefato
conferível não conta. Escreva em português do Brasil.

## Regras

- **Pode ler:**
  - `scripts/*.py`;
  - `build/rouche_L_rig/**` (JSON, txt, npy) e `build/rouche_L_est/**`;
  - `build/s4/rig/E.json`;
  - `docs/S3_CONDICIONAMENTO_SHARP.md` §10.1–10.3 (definição do espaço com sinal);
  - o artigo e os dados de RT em `.cache/rt-1203.3766v1/`;
  - o PDF de Rump em `.codex-runs/2026-09-24-revisao-agente/`.
- **NÃO pode ler** (são a narrativa e as conclusões do autor, ou de revisões anteriores; as
  alegações estão enunciadas abaixo):
  - `docs/C1_REAVALIACAO.md`, `docs/PROOF_OBLIGATIONS.md`, `docs/RETOMADA.md`, `docs/PRIORIDADES.md`;
  - `docs/S4_GAUGE.md`, `docs/S3B_ROTA.md`, `docs/AUDITORIA_*.md`;
  - nada em `.codex-runs/` além desta pasta e do PDF de Rump;
  - nem a memória do autor.
- **Não altere nenhum arquivo existente.** Escreva SÓ em `.codex-runs/2026-10-01-revisao-c1/`:
  scripts `rN_*.py`, saídas `rN_saida.txt` e o `RELATORIO.md`. Os scripts que gravam em `--dir`
  devem receber uma CÓPIA dentro desta pasta.
- Toda conta usa `CHOPTUIK_EPS_FUNDO_L=4.7e-08` no ambiente. Use `.venv/bin/python` (o Python do
  sistema não tem scipy) e `OPENBLAS_NUM_THREADS=2` no notebook.
- O snapshot do estado revisado está em `.snapshots/2026-10-01-revisao-c1/`:
  - `MANIFEST.sha256` (local, commit em `COMMIT.txt`);
  - `MANIFEST_lab2_remoto.sha256` (os arquivos grandes no laboratorio2).
- Máquinas: o notebook local tem 7 GB de RAM e 4 núcleos. Para contas com as matrizes densas de
  n = 14.016 (A, T, U, cada uma com 3,1 GB complexa), use:
  - o **laboratorio2**: `<usuario>@<host>`, 14 GB, 4 núcleos, disponível 24 h. Acesso por
    senha via
    `SSH_ASKPASS=<scratch> SSH_ASKPASS_REQUIRE=force DISPLAY=:0 ssh -o PubkeyAuthentication=no <usuario>@<host> ... < /dev/null`
    (não copie a senha para lugar nenhum). Os dados estão em `~/Choptuik/build/rouche_L/`: A.npy,
    T.npy, U.npy e defl.npy, os fatores F dos centros 0 e ±0,25i, e `schur.json` e `globais.json`.
    Os fatores dos centros 0,55 ± 0,25i, 1,0 … 3,0 que estão em `~/Choptuik/build/rouche_L/` são
    **OBSOLETOS** (26/09, de outra versão); os corretos estão em `~/Choptuik/build/certF_extra/`
    e, no notebook, em `build/rouche_L_rig/F/` (todos os 10).
  - o **PC 3**: `<usuario>@<host>`, 15 GB, 4 núcleos, via
    `SSH_ASKPASS=<scratch> SSH_ASKPASS_REQUIRE=force DISPLAY=:0 ssh -o PubkeyAuthentication=no <usuario>@<host> ... < /dev/null`.
    Os dados estão em `~/choptuik_trabalho/` (mesma A, sha256 70791ce0…). **É o PC pessoal de
    outra pessoa:** nada de root, nada de mudanças no sistema e nada de instalar pacotes.
  - **NÃO use** o PC do laboratório (<host>): ele está em uso por outras pessoas.

  Nas máquinas remotas, trabalhe SÓ numa pasta nova `~/revisao_c1/` (copie o que precisar; não
  altere `~/Choptuik` nem `~/choptuik_trabalho`). Uma conta com A, T e U juntas usa ~10 GB: no
  laboratorio2 rode uma por vez. Rode coisas longas com `setsid nohup ... &` e acompanhe; ao matar
  processos, use o PID (`pkill -f` casa com a própria linha de comando).
- Orçamento: até ~3 h de trabalho. Priorize o que pode estar errado e o que é novo.

## Contexto

L(s) é o operador de Reiterer–Trubowitz (RT) linearizado em torno da solução crítica, na base de
Fourier COM SINAL (`scripts/signed_operator_L.py`), com a norma do toro
`||x||² = Σ κ1^{2m} ω2(n) |x_mn|²`, (κ1, κ2) = (65/64, 5/4), ω2(n) = κ2^{2n} + κ2^{-2n}.
L(s) = J_μ + B + s, com J_μ a parte livre (diagonal em m, triangular superior em n) e B o fundo.
Q0(s) = (J_μ + s)⁻¹ em forma fechada por modo. A maquinaria de cauda (`nk_L3_T.py`: janelas
bloco-Jacobi com V_J = Ĥ_J⁻¹, Loewner, `verified_inverse.py`, `fundo_L.py`) já foi revisada em
duas aplicações anteriores, em s real. **Foque no que é novo aqui:**
- o Rouché no contorno com referência dependente de s;
- o nível Z por Schur;
- o fator de posto baixo;
- a sensibilidade em s;
- os blocos de janela no interior;
- a contagem finita;
- a deflação.

**Alegação final:** N(L; 0 < Re s < 3, |Im s| < 1/4) = 3, com multiplicidade algébrica: μ ≈ 0,1683,
0,4010 e 0,7332.

## O que é alegado

**A1 (formulação e homotopia).** Ω = (0, 3) × (−1/4, 1/4).
- **Deflação de Brauer.** L̃(s) = L(s) + Σ_{k=0,1} δ_k x_k φ_k*, onde:
  - x_0 = ∂τω* é o autovetor exato de fase (L(0)x_0 = 0) e x_1 = g* = (Z* + 1)ω* o de gauge
    (L(μ*)x_1 = 0, L(s) = L(0) + s);
  - φ_k ∈ span{x_0^A, x_1^A}, com φ_i* x_j^A = δ_ij (biortogonal com os vetores aproximados de RefA
    na caixa Z, `rouche_L_Z.py schur --vetores refA --phi biortogonal`);
  - δ = (1, 1 + μ_A).
- **Operador e referência.** H(s) = L̃(s) Q0(s) tem os mesmos zeros que L̃ em Ω. A referência é
  A(s) = diag(Ĥ_ZZ(s), Ĥ_J(s) nas 26 janelas, I no far), com Ĥ_ZZ(s) = P_Z H(s) P_Z, tomada como a
  matriz de ponto flutuante A + s conjugada (A = W L̃_Z(0) W⁻¹ em `A.npy`, exata como dado; a
  diferença para a verdadeira entra como perturbação), e Ĥ_J(s) = P_J H(s) P_J.
- **Rouché.** Com os inversos aproximados diag(V_Z(s), V_J(c_k), I), V_Z implícito na escala
  D_Z = J_Z(c) Q_ZZ(s), um Perron de 3 níveis com θ < 1 em cada ponto do contorno dá
  H_t = A + t(H − A) invertível para todo t ∈ [0, 1]. Então, por Gohberg–Sigal,
  N(H, Ω) = N(A, Ω) = N(Ĥ_ZZ, Ω) + Σ_J N(Ĥ_J, Ω).
- **Escala no nível Z.** A entrada T←Z vira ‖V_J P_J B Q_ZZ(c)‖, independente de s. A entrada Z←T
  vira ‖J_Z(c)[R(s)Ỹ − Q_ZZ(s)J_ZT] Q0_TT(s)‖, com:
  - R(s) = (A + s)⁻¹, Ỹ = P_Z L̃(0) P_T e J_ZT = P_Z J P_T;
  - a identidade J_ZZ Q_ZT + J_ZT Q_TT = 0.

**A2 (nível Z por ponto, `rouche_L_rig.py ponto`).** Por ponto p do contorno:
- t_p >= ‖(A + p)⁻¹‖ por Loewner (`verified_inverse.cota_loewner`, cota de arredondamento de
  A A* elemento a elemento). Valores de t_p de execuções anteriores, impressos no log, foram
  reaproveitados (`--tps`, arredondados para cima pela casa impressa).
- X_k ≈ R(p)^k F_t por Schur, com resíduo certificado, e as cotas a0, b1, b2 = ‖J_Z(c)[Y_{j+1} −
  Q_ZZ(p)^{j+1} F_b]‖, ‖Y_3‖ e ‖Q_ZZ(p)³ F_b‖.
- Para |s − p| <= r, pela identidade do resolvente:
  ‖J_Z(c)[R(s)F_t − Q_ZZ(s)F_b]‖ <= E = a0 + r b1 + r² b2 + r³(τ‖Y_3‖ + JQ‖Q³F_b‖), com:
  - τ = 1 + (|c − s| + ‖B̂‖) t_p/(1 − r t_p), B̂ = A − J_Z(0), ‖B̂‖ em `globais.json`;
  - JQ = 1 + d nQZ/(1 − r nQZ).

**A3 (fator de posto baixo, `rouche_L_Z.py certF`).** S = [[(B + J)_ZT Q0_TT(c)], [J_ZT Q0_TT(c)]]
(exata; B são os blocos calculados `ctx.B`, tomados como dados).
- Para qualquer Φ: ‖[M1, M2] S‖ <= ‖[M1, M2] F‖ ‖Φ‖ + ‖[M1, M2]‖ ρ, com ρ >= ‖S − FΦ‖.
- Φ_m = V f(Σ) U* S_m (Tikhonov sobre a SVD de F) por modo m; ‖Φ‖ e ρ são acumulados com
  arredondamentos elemento a elemento e o erro de Q0 pelo resíduo do modo. q_TT é certificado.
- Saídas: `build/rouche_L_rig/F/certF_*.json`, com o sha256 do F usado.

**A4 (cauda e sensibilidade, `nk_L3_T.py --rouche` → `build/rouche_L_rig/cauda/T3_*.json`).**
Em cada centro c, N(c), TZ_J, K_J = ‖V_J P_J B Q0(c)² P_J‖, nQ por janela e o far. Para |s − c| <= d,
com fJ = 1/(1 − d nQ_JJ):
- janela central: ρ_J(s) <= ρ_J(c) + d[TZ_J q_ZJ(s) + (K_J + TZ_J q_ZJ(c)) fJ];
- janela não central: ρ_J(s) <= ρ_J(c) + d K_J fJ;
- fora da diagonal: N_ij(s) <= N_ij(c) + d[TZ_i q_Zj(s) + N_ij(c) nQ_jj fJ_j].
Os termos restantes estão na docstring de `scripts/rouche_L_disco.py`.

**A5 (montagem por disco, `rouche_L_disco.py`).** Para cada um dos 96 pontos (10 centros,
`build/rouche_L_rig/discos/disco_*.json`, com o rig em `build/rouche_L_rig/z/<máquina>/` e a cauda em
`build/rouche_L_rig/cauda/`), o maior r tal que o Perron 3×3 (Z', T, far), com Perron em racionais,
fica < 1 para todo s do contorno com |s − p| <= r. Todos os 96 discos dão θ <= 0,99916.

**A6 (blocos de janela sem zeros em Ω, `rouche_L_janelas.py` → `build/rouche_L_rig/janelas.*`).**
- Ω é particionado em 10 retângulos, um por centro.
- F_k(s) = I − V_J(c_k)Ĥ_J(s) é analítica, Q0 é analítica em Re s > −μ e ‖F_k‖ é subharmônica.
  Logo a cota da A4, verificada no bordo de cada retângulo (amostras mais Neumann estruturado),
  vale no interior.
- Máximo: 0,149 < 1. Então N(Ĥ_J, Ω) = 0.

**A7 (cobertura, `rouche_L_cobertura.py` → `cobertura.txt`).** Os 96 discos cobrem ∂Ω (racionais).

**A8 (contagem finita, `rouche_L_contagem.py` → `contagem.json`, `contagem_rouche.txt`).**
- ‖A − UTU⁻¹‖ <= ε_B = 5,97·10⁻⁵ (Frobenius, por blocos).
- No contorno, ‖(A + s)⁻¹‖ <= max t_p/(1 − r t_p) <= 560,3, e t·ε_B <= 0,033 < 1.
- Logo N(Ĥ_ZZ, Ω) = #{i : −T_ii ∈ Ω} = 2 (0,401025 e 0,733180).

**A9 (Brauer).** D(s) = det(I + L⁻¹ Σ δ_k x_k φ_k*) é racional, e s(s − μ*)D(s) = q(s), com
q(s) = (s + δ_0(1 + e_00))(s − μ* + δ_1(1 + e_11)) − δ_0 δ_1 e_01 e_10, e_ij = φ_i*(x_j − x_j^A).
q não tem zeros em Ω̄, e o polo de D em μ* está em Ω. Com o contorno indentado em 0 (L(0) singular,
L̃(0) invertível): N(L, Ω) = N(L̃, Ω) + 1 = 3. As cotas de ‖x_j − x_j^A‖ vêm de `gauge_refB.py`
(1,46·10⁻⁹) e de `build/s4/rig/E.json` (2,5·10⁻⁸).

## Artefatos obrigatórios (um por alegação; script + saída + veredicto)

- **R1 (A1, lógica).**
  - Escreva a prova, ou ache a falha, de que θ < 1 (na norma max ponderada do Perron) implica H_t
    invertível para todo t.
  - Mostre que Gohberg–Sigal se aplica: H(s) − I é compacto e analítico? Q0 é analítica perto de Ω?
  - Mostre que N(A, Ω) se decompõe como alegado, com a referência dependente de s e os V_J(c_k)
    trocando de centro ao longo do contorno.
  - Confira a identidade da escala s-dependente.
- **R2 (A2, numérico).**
  - No laboratorio2, para um ponto de um centro difícil (0,25i ou 0) e 3 a 5 s sorteados no disco
    (no contorno), calcule por LU densa o valor verdadeiro de ‖(A + s)⁻¹‖ e de
    ‖J_Z(c)[R(s)F_t − Q_ZZ(s)F_b]‖. Confira que ficam abaixo de t(s) e de E.
  - Re-derive a expansão de Taylor e os papéis de τ e JQ.
- **R3 (A3).**
  - Derive a desigualdade.
  - Recalcule ‖Φ‖ e ρ para um centro com código próprio (pode ser outro Φ, por exemplo mínimos
    quadrados puros), montando S por conta própria a partir de `ctx.B` e Q0.
  - Confira o sha256 do F e que esse F é o mesmo que o rig usou.
- **R4 (A4).**
  - Derive as três fórmulas de sensibilidade, inclusive o papel de q_ZJ(s) e de fJ.
  - Teste numericamente, numa truncagem menor (a cauda menor que couber no notebook) ou numa
    janela da truncagem cheia: calcule ρ_J(s) verdadeiro para s a uma distância d do centro e
    compare com a cota.
- **R5 (A5).**
  - Para 3 discos (um de 0,25i, um de 0 e um de 1,0), recalcule com código próprio, a partir dos
    JSONs, a matriz 3×3 e θ em racionais.
  - Confira cada termo contra a docstring e contra a sua própria derivação. Procure termos
    faltando: arredondamentos, a perturbação do fundo, a deflação verdadeira contra a aproximada
    (dEZ, dET, nE), e μ.
- **R6 (A6).**
  - Confira a analiticidade (inclusive em Re s = 0) e que a partição cobre Ω.
  - Recalcule a cota em 2 pontos de amostra do bordo.
  - Diga se o argumento do princípio do máximo está certo para operadores (norma de função
    analítica matricial).
- **R7 (A7).** Escreva o seu próprio verificador de cobertura (racional), incluindo os cantos, e
  rode-o sobre os JSONs.
- **R8 (A8).**
  - No laboratorio2, recalcule ‖AU − UT‖ e ‖I − U*U‖ com código próprio.
  - Confira o argumento (autovalores de B = UTU⁻¹ = diagonal de T, homotopia).
  - Confira o max t(s) dos discos.
  - Se der tempo, confira por outro caminho que A tem 2 autovalores em −Ω: por exemplo,
    `scipy.linalg.eigvals` em A (só como estimativa) ou o argumento de det por LU num contorno grosso.
- **R9 (A9).**
  - Derive D(s) e q(s), e confira os valores de δ e de e_ij.
  - Confira que a indentação em 0 e o polo em μ dão exatamente +1.
  - Diga se a contagem com multiplicidade algébrica é a que Gohberg–Sigal entrega aqui.
- **R10 (consistência e reprodução).**
  - Confira que:
    - o λ0 de cada T3 bate com o centro;
    - os pontos de cada rig são os do plano (`build/rouche_L_est/plano_rig.*`);
    - todo disco usou certF;
    - os sha256 de A nas máquinas são iguais.
  - Rode `rouche_L_disco.py` para um centro, `rouche_L_cobertura.py` e
    `rouche_L_contagem.py rouche` numa cópia, e confira que reproduzem.
  - Grave os sha256 dos originais antes e depois.

## Relatório

`.codex-runs/2026-10-01-revisao-c1/RELATORIO.md` deve trazer:
- uma tabela R1–R10 com veredicto (OK / OK com ressalva / LACUNA / ERRO), artefato e o que foi
  feito;
- a lista de lacunas, com gravidade e um conserto proposto;
- o que você não conseguiu verificar.
