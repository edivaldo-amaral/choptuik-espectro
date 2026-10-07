# Revisão independente da faixa B (C4) do Rouché de L (02/10/2026)

Você é o revisor independente de um certificado assistido por computador. O autor é outro modelo. Seu
trabalho é **tentar quebrar** o certificado, não confirmá-lo. Um "OK" sem artefato conferível (script +
saída) não conta. Escreva em português do Brasil.

## Regras

- **Pode ler:**
  - `scripts/*.py` e `build/faixaB/**`;
  - `build/rouche_L_rig/**`: o certificado da faixa A, já revisado; use para comparar e para os discos
    da aresta Im = 1/4;
  - `build/s4/rig/E.json` e `build/s4/corr_gauge.json`;
  - `docs/SECTOR_B.md`, `docs/SPECTRAL_PROBLEM.md` e `docs/CONSTRAINT_ROOT_QUOTIENT.md`. Foram escritos
    por outro modelo: trate-os como **alegações a conferir**, não como premissas;
  - `docs/S3_CONDICIONAMENTO_SHARP.md` §10.1–10.3 (espaço com sinal);
  - o artigo e os dados de RT em `.cache/rt-1203.3766v1/`;
  - o PDF de Rump em `.codex-runs/2026-09-24-revisao-agente/`.
- **NÃO pode ler** (são a narrativa e as conclusões do autor, ou de revisões anteriores):
  - `docs/C1_REAVALIACAO.md`, `docs/PROOF_OBLIGATIONS.md`, `docs/RETOMADA.md`, `docs/PRIORIDADES.md`;
  - `docs/AUDITORIA_*.md`, `docs/S4_GAUGE.md`, `docs/S3B_ROTA.md`;
  - nada em `.codex-runs/` além desta pasta e do PDF de Rump;
  - a memória do autor.
- **Não altere nenhum arquivo existente.** Escreva SÓ em `.codex-runs/2026-10-02-revisao-B/` (`rN_*.py`,
  `rN_saida.txt`, `rN_resultado.json` e `RELATORIO.md`) e, nas máquinas remotas, SÓ em `~/revisao_B/`.
  Scripts que gravam em `--dir` recebem uma CÓPIA.
- Ambiente: `CHOPTUIK_EPS_FUNDO_L=4.7e-08` em toda conta; `.venv/bin/python` (o Python do sistema não tem
  scipy). No notebook (7 GB de RAM, 4 núcleos), use `OPENBLAS_NUM_THREADS=2`.
- Snapshot do estado revisado: `.snapshots/2026-10-02-revisao-B/`:
  - `MANIFEST.sha256` (local), commit em `COMMIT.txt`;
  - `MANIFEST_lab2_remoto.sha256` e `MANIFEST_joao_remoto.sha256`.
- **Máquinas** (ssh direto, sem senha: `ssh -o BatchMode=yes ...`):
  - **laboratorio2**: `<usuario>@<host>`, 14 GB, 4 núcleos, disponível 24 h.
    - `~/Choptuik/build/rouche_L_v2/`: A.npy, T.npy e U.npy (links; A com sha 70791ce0…), os fatores
      F_*.npy, os níveis Z `rig_<c>_B*.json` e os `tp_<c>_B*.json`.
    - `~/Choptuik/build/faixaB/nk/`: o NK do testemunho (Z.json, Q.json, T3.json, C3.json, D.json,
      h0w.npy).
    - `~/Choptuik/build/caudaB_rig/0_5/`: os checkpoints T3.lab2.json e T3.joao.json e o T3.json juntado.
    - `~/Choptuik/build/caudaB_rig/0_75/`.
  - **PC 3**: `<usuario>@<host>`, 15 GB, 4 núcleos.
    - `~/choptuik_trabalho/build/rouche_L/` (A.npy com o mesmo sha) e `rouche_L_v2/`.
    - `~/choptuik_trabalho/build/faixaB/nk/` (etapas A e B do NK).
    - `~/choptuik_trabalho/build/caudaB_rig/0_5/T3.parcial.json`.
    - **É o PC pessoal de outra pessoa:** nada de root, nada de mudanças no sistema e nada de instalar
      pacotes.
  - **NÃO use** o PC do laboratório (<host>): está em uso por outras pessoas.
  - Uma conta com A, T e U juntas usa ~10 GB; a Loewner pré-condicionada, ~12 GB. Rode uma de cada vez
    por máquina, com `setsid nohup ... &`, e acompanhe com `until ...; do sleep 30; done` (sleep isolado
    longo é bloqueado). Mate processos pelo PID (`pkill -f` casa com a própria linha de comando).
  - O notebook pode desligar: grave saídas parciais em arquivo com frequência.
- Orçamento: até ~3 h de trabalho. Priorize o que é novo e o que pode estar errado.

## Contexto

- **O operador.** L(s) = J_μ + B + s, na base de Fourier COM SINAL do setor A (`scripts/signed_operator_L.py`),
  com a norma do toro ‖x‖² = Σ κ1^{2m} ω2(n) |x_mn|². Q0(s) = (J_μ + s)⁻¹ tem forma fechada por modo.
- **A faixa A, já certificada e revisada.**
  - N(L; 0 < Re s < 3, |Im s| < 1/4) = 3: Rouché de operadores com referência
    A(s) = diag(Ĥ_ZZ(s), Ĥ_J(s), I).
  - Deflação de Brauer da fase e de μ (L̃ = L + Σ δ_k x_k φ_k*).
  - Discos com Perron de 3 níveis e janelas sem zeros pelo princípio do máximo.
  - Contagem de Schur com ‖A − UTU⁻¹‖ <= ε_B.
  - A maquinaria (rouche_L_rig v2, rouche_L_disco, rouche_L_janelas, certF, nk_L3_T, cobertura e contagem)
    foi revisada na faixa A. **Foque no que é novo na faixa B.**

## O que é alegado

**B1 (formulação).**
- Pela redução do setor B ao A (Σ_B = Σ_A − i/2, `SECTOR_B.md`), a faixa B é
  Ω_B = (0, 3) × (1/4, 3/4) para L_A.
- O espectro de L_A é i-periódico (m → m + 2), então A ∪ B cobre um período em Im s.
- F3 (`R_L <= 2,926`: Re W(−B_L) <= 2,9492 e parte livre <= −0,0232 no toro) vale para todo Im s.
- Brauer: N(L, Ω_B) = N(L̃, Ω_B), porque os zeros de q (~ −1) e os polos 0 e μ de D(s) estão fora de Ω_B.
- A referência é a mesma da faixa A (a mesma A.npy, uma só função analítica de s). Então os discos da
  faixa A em Im = 1/4 (`build/rouche_L_rig/discos_v2/`) valem como parte de ∂Ω_B.

**B2 (Loewner pré-condicionada, `rouche_L_rig.tp_precond`).**
- (A + p)⁻¹ = K⁻¹Q, com Q = Q_ZZ(p) = J_Z(p)⁻¹ exata, B̂ = A − J_Z(0) e K = I + QB̂.
- Loewner: t²K_cK_c* − Q_cQ_c* ⪰ 0 ⇒ ‖K_c⁻¹Q_c‖ <= t.
- Passagem para os exatos, com ‖K_e − K_c‖ <= dK e ‖Q_e − Q_c‖ <= dq:
  ‖(A + p)⁻¹‖ <= t/(1 − kc·dK) + dq·kc/(1 − kc·dK), onde kc = t·nJ/(1 − nJ·dq).
- Usada quando a estimativa passa de 40: até t_p = 201,9 em 0,5i.

**B3 (nível Z da faixa B).**
- Os `rig_<c>_B*.json` (v2: sha256 de A e de F registrados) foram feitos com a fase 1 e a fase 2 em
  processos separados; a fase 2 reaproveita os t_p da fase 1 (`--tps`, do `tp_<c>_B*.json`).
- Pontos novos nos centros 0,25i, 0,5i, 0,75i, 0,55 + 0,75i e 2,5. Para 1,0, 1,5, 2,0 e 3,0, as caudas
  reais da faixa A foram reaproveitadas na aresta Im = 3/4 (`build/faixaB/z/rig_*_Bteste.json`).

**B4 (caudas novas).**
- 0,5i, 0,75i e 0,55 + 0,75i, com `nk_L3_T --rouche`.
- A de 0,5i foi feita em duas máquinas e juntada: checkpoints T3.parcial.json (João, ordem direta) e
  T3.lab2.json (laboratorio2, ordem reversa), unidos com `--juntar` numa passada final no laboratorio2.
  Mesmo código (md5) e mesma estimativa `nk3_0_5.json`, esta feita no PC do laboratório.
- 8 janelas foram calculadas nas duas máquinas.

**B5 (montagem por disco).** `build/faixaB/discos/*.json`: 79 discos, todos com θ < 1 (máx 0,99907) e
Perron em racionais. As fórmulas de `rouche_L_disco.py` são as da faixa A (versão 2).

**B6 (janelas, `rouche_L_janelas.py --faixa B`).** Ω_B particionado em 10 retângulos (`REGIOES_B`). Máximo
de ρ_J no bordo: 0,1517. Logo N(Ĥ_J, Ω_B) = 0.

**B7 (cobertura e contagem).**
- Os 79 discos da faixa B mais os da faixa A em Im = 1/4 cobrem ∂Ω_B (`build/faixaB/cobertura.txt`).
- Contagem: −diag(T) (`build/rouche_L_rig/diagT_lab2.npy`, do T.npy do laboratorio2) tem um elemento em
  Ω_B, 0,0414194105 + 0,5i. Com ε_B de `build/rouche_L_rig/contagem.json` e max t <= 559,
  t·ε_B <= 0,0334 (`build/faixaB/contagem_rouche.txt`).
- Logo N(L, Ω_B) = 1.

**B8 (NK do testemunho, `build/faixaB/nk/`).**
- λ0 = 0,0414194105346554 + 0,5i: autovalor da truncagem 32×128 de L sem deflação (Newton bordejado na
  matriz pesada).
- h0 por iteração inversa na matriz PESADA (`nk_L2_Z.py`).
- tQ >= ‖Mh⁻¹[Q_ZZ; 0]‖ por Loewner (`nk_L2_Z.py --so-tQ`, com o mesmo h0w.npy) substitui t_V·nQZ em Z2
  (`nk_L3_C.py`).
- Resultado: θ = 0,9115, Z2 = 155, Y = 6,0·10⁻⁶ e r = 7,8·10⁻⁵, logo há um zero único e λ* é
  algebricamente simples, com |λ* − λ0| <= 3,7·10⁻⁵.
- Etapa D (`nk_L3_D.py`): ‖Π_F C(λ*)h*‖ >= 0,35469 > 0, logo o autovetor viola as constraints.

**B9 (conclusão).**
- N(L, Ω_B) = 1, λ* é essa raiz, simples e violadora: a multiplicidade física na faixa B é 0.
- Com a faixa A (μ é gauge, s_K = 0,401 viola as constraints e 0,733 é física porque K é injetivo ali),
  a contagem física total em Re s > 0 é 1.
- A injetividade de K em 1,765 < Re s < R deixa de ser necessária.

## Artefatos obrigatórios (um por alegação; script + saída + veredicto)

- **R1 (B1).**
  - Confira a redução e a i-periodicidade no código (paridades por componente em `signed_operator_L`).
  - Confira que F3 não depende de Im s.
  - Confira a contabilidade de Brauer em Ω_B.
  - Confira que os discos de Im = 1/4 da faixa A certificam a mesma referência que a faixa B usa.
  - Procure qualquer dependência "da faixa A" nas fórmulas da montagem: Re s >= 0, o eps_mu_L, as formas
    fechadas (cauda_Rinv_disco, Rinv_uniforme com |Im s| até 0,75 + r) e a deflação.
- **R2 (B2).**
  - Derive a cota.
  - No laboratorio2 ou no PC 3, para 2 a 3 pontos perto de 0,5i (t ~ 150 a 200), calcule
    ‖(A + p)⁻¹‖ verdadeiro (LU mais Lanczos ou SVD parcial) e confira que t_p fica acima.
  - Procure erros que faltem em dK: J_Z(0) calculada, Q_cert, o produto Q·B̂ por blocos e a soma com I.
- **R3 (B3).**
  - Confira em todos os rig_*_B*.json: sha de A = 70791ce0…, sha de F = o do certF do centro, pontos =
    plano, e t_p = o do tp_*.json correspondente (fase 1).
  - Para 1 ponto, refaça a0 com LU densa e compare.
- **R4 (B4).**
  - Compare as 8 janelas em comum entre T3.parcial.json (João) e T3.lab2.json: os valores devem
    coincidir (mesmo código, mesma entrada).
  - Confira que o T3.json juntado tem as 26 janelas e que cada entrada vem de um dos checkpoints.
  - Para uma janela de 0,5i, recalcule ‖V_J‖ por conta própria.
- **R5 (B5).** Para 3 discos (o pior de 0,5i perto de 0,49i, o de 1,0 + 0,75i com cauda real e d ~ 0,87,
  e um de 0,75i no canto), remonte a matriz 3×3 e θ em racionais com código próprio, a partir dos JSON.
- **R6 (B6).**
  - Confira que a partição cobre Ω_B.
  - Recalcule a cota em 2 pontos de amostra, inclusive da região do centro 1,0, que vai até
    Im = 0,75 com d ~ 0,8.
- **R7 (B7).**
  - Verificador de cobertura próprio para ∂Ω_B (racional, com os cantos).
  - Confira por outro caminho que A tem exatamente 1 autovalor em −Ω_B: argumento de det por LU num
    contorno grosso, ou eigvals em parte, como estimativa.
- **R8 (B8).**
  - Confira λ0, a iteração pesada e a Loewner de tQ, e se Z2 com tQ está correto (DF(x) − DF(x0) e as
    componentes).
  - Recalcule c_ref = ‖Π_F C(λ0)h0‖ com código próprio.
  - Confira que λ* está em Ω_B e que a unicidade do NK mais N = 1 identificam λ* com a raiz da contagem.
- **R9 (B9).** Examine a lógica física e diga se a contagem física 1 está bem fundamentada e o que é
  hipótese.
- **R10 (reprodução).**
  - Numa cópia, rode `rouche_L_disco.py` para um centro, `rouche_L_cobertura.py --regiao 0,3,1/4,3/4` e
    `rouche_L_contagem.py rouche --regiao 0,3,0.25,0.75 --diag ...`.
  - Grave os sha256 dos originais antes e depois.

## Relatório

`.codex-runs/2026-10-02-revisao-B/RELATORIO.md`:
- uma tabela R1–R10 com veredicto (OK / OK com ressalva / LACUNA / ERRO), artefato e o que foi feito;
- a lista de lacunas, com gravidade e conserto proposto;
- o que você não conseguiu verificar.
