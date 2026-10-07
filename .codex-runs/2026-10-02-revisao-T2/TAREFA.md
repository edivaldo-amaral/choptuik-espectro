# Revisão independente de T2: lemas novos, passos em papel e cadeia lógica (02/10/2026)

Você é o revisor independente de um resultado assistido por computador. O autor é outro modelo. Seu
trabalho é **tentar quebrar** o resultado, não confirmá-lo. Um "OK" sem artefato conferível (script +
saída, ou derivação escrita passo a passo) não conta. Escreva em português do Brasil.

## O objeto da revisão

`docs/T2_ENUNCIADO.md` (o teorema condicional e a cadeia da prova) e `docs/T2_HIPOTESES.md` (o inventário
e os lemas L-real, L-afim, L-iso e L-constr). As **contagens** (N_A = 3 e N_B = 1) e os certificados NK
(S3b, S4 e a faixa B) já foram revisados por outras revisões independentes: tome-os como dados e não os
refaça. **Foque no que não foi revisado:**
- os lemas novos;
- os passos em papel 1, 2, 3, 4 e 12 da tabela de `T2_ENUNCIADO.md`;
- a lógica que leva das contagens a "1 modo instável físico".

## Regras

- **Pode ler:**
  - `docs/T2_ENUNCIADO.md` e `docs/T2_HIPOTESES.md`;
  - `docs/SECTOR_B.md`, `docs/SPECTRAL_PROBLEM.md`, `docs/CONSTRAINT_ROOT_QUOTIENT.md`,
    `docs/RADIUS_RECOVERY.md`, `docs/SHARP_LINEARIZED_IDENTITY.md`, `docs/SHARP_DOMAIN.md`,
    `docs/STRUCTURED_EXTERIOR.md`, `docs/GAUGE_RT_ACTION.md`, `docs/GAUGE_DOMAIN_CLASSIFICATION.md` e
    `docs/GEOMETRIC_FLOQUET_RECONSTRUCTION.md`;
  - `docs/S3_CONDICIONAMENTO_SHARP.md` §9–10, `docs/ALTERNATIVE_ROUTES.md` §5.26 (F3) e
    `docs/S4_GAUGE.md`.

  Vários desses documentos foram escritos por outro modelo (14–16/09) e auditados em 16/09: trate tudo
  como **alegações a conferir**. Também pode ler `scripts/*.py`, `build/t2/**`, `build/faixaB/**`,
  `build/rouche_L_rig/**`, `build/s4/**` e o artigo, o código e os dados de RT em `.cache/rt-1203.3766v1/`.
- **NÃO pode ler** (são narrativa e conclusões do autor ou de revisões anteriores):
  - `docs/C1_REAVALIACAO.md`, `docs/PROOF_OBLIGATIONS.md`, `docs/RETOMADA.md`, `docs/PRIORIDADES.md`;
  - `docs/AUDITORIA_*.md`;
  - nada em `.codex-runs/` além desta pasta;
  - a memória do autor.
- **Não altere nenhum arquivo existente.** Escreva só em `.codex-runs/2026-10-02-revisao-T2/` e, no PC do
  João, só em `~/revisao_T2/`.
- Ambiente: `.venv/bin/python` e `CHOPTUIK_EPS_FUNDO_L=4.7e-08`. No notebook (7 GB, 4 núcleos), use
  `OPENBLAS_NUM_THREADS=2`.
- **Máquinas:** o notebook e o **PC 3**, `ssh -o BatchMode=yes <usuario>@<host>`, sem senha,
  15 GB, 4 núcleos, com o projeto em `~/choptuik_trabalho/`. É o PC pessoal de outra pessoa: nada de root,
  mudanças no sistema ou instalação de pacotes.
  - **NÃO use** o laboratorio2, que roda outro cálculo, nem o PC do laboratório.
  - Trabalhos longos: `setsid nohup ... &`, acompanhados com `until ...; do sleep 30; done`. Os processos
    em segundo plano do harness morrem em ~30 min. Mate pelo PID.
- Snapshot: `.snapshots/2026-10-02-revisao-T2/MANIFEST.sha256` (commit em `COMMIT.txt`).
- O notebook pode desligar: grave resultados parciais e vá atualizando o RELATORIO.md.
- Orçamento: até ~3 h.

## Artefatos obrigatórios

- **R1 (L-real, lado do toro).**
  - Confira a desigualdade ‖T(H − I)T‖_toro <= (‖B_L‖_toro + \|s\|)·‖J⁻¹T‖_toro e a escolha da caixa.
  - Confira que `closed_form.cauda_Rinv` e a cota do modo inteiro (1 + 2c_w/(κ2 − 1))/(\|m\|/2) cotam a
    norma de J⁻¹ restrita às colunas T no **toro com sinal**. O peso κ1^m cancela dentro de um modo? A
    convenção de componentes e coeficientes simétricos é a mesma?
  - Para 3 a 5 modos, recalcule ‖J⁻¹_m P_{n >= N}‖ por SVD de uma truncagem grande e compare.
  - Confira que ‖B_L‖ <= 3,9641 (F3) é a norma no MESMO espaço (pesos ω2(n) contra κ2^{2n} com n ∈ ℤ).
  - Refaça `scripts/l_real_toro.py` ou um subconjunto, e confira `build/t2/l_real_toro.json`.
- **R2 (L-real, lado de RT).**
  - Confira `radius_transfer.transfer_bounds`, `structured_exterior.structured_input_tail`,
    `sharp_exterior.free_tail` e β+ = 5191/500, e a extensão a \|s\| <= 3,1 com Q = J_μ⁻¹ independente
    de s.
  - Diga quais constantes vêm da bola publicada de RT e se a leitura delas é correta (confira no artigo
    o que puder).
- **R3 (L-real, o argumento).**
  - Escreva a prova completa, ou ache a falha. Pontos a conferir:
    - a inclusão Y+ ⊂ toro, com índices, coeficientes e a norma;
    - a mesma caixa G nos dois espaços;
    - J⁻¹ preserva a caixa;
    - B é limitado em Y+;
    - a unicidade no toro;
    - as cadeias de Jordan.
  - Diga se a conclusão "contagens de Rouché no toro = multiplicidades algébricas na realização de RT"
    segue, e o que mais seria preciso. Por exemplo: o Rouché conta zeros de H̃ = L̃Q0(s) no toro;
    isso é a multiplicidade algébrica de L no toro?
- **R4 (L-afim).** Derive L(s) = L(0) + s·I e C(s) = C0 + sC1 a partir das equações de RT (artigo e
  `SPECTRAL_PROBLEM.md`), e confira no código (`signed_operator_L.J_livre`, `signed_constraints`).
- **R5 (L-iso, L-constr e o Lema 5.1 de `SECTOR_B.md`).**
  - Derive a conjugação M⁻¹L_B(s)M = L_A(s + i/2) e M⁻¹C_B(s)M = C_A(s + i/2).
  - Confira numericamente, nos blocos do código, que J_{m+1}(s) = J_m(s + i/2) por componente e que os
    blocos de B e de C só dependem de d = m_out − m_in e da componente.
  - Confira a norma de M no toro com sinal e o Lema 5.1, inclusive o domínio fundamental semiaberto e as
    retas Im s = ±1/4.
- **R6 (o quociente das constraints).**
  - Confira `CONSTRAINT_ROOT_QUOTIENT.md`: M1 = C1, M0 = C0 + BC1 − C1A, BD = DA, D = C0 − C1A,
    dim E_constraints = dim E − rank(D|E) e a subtração do gauge.
  - Faça um teste com matrizes aleatórias que satisfazem KC = ML.
- **R7 (KC = ML).**
  - Confira a identidade de `SHARP_LINEARIZED_IDENTITY.md` no que for possível.
  - Pelo menos teste numericamente, nos blocos do código, K(s)C(s) ≈ M(s)L(s) numa truncagem, ou
    explique por que não é testável assim. Veja também `test_root_constraints.py` e
    `independent_rt_symbolic.py`.
- **R8 (a cadeia de T2).** Confira a lógica de `T2_ENUNCIADO.md`, da contagem a "1 modo físico":
  - μ*: precisa de simplicidade E Dg = 0 E H-cone/H-gauge?
  - s_K e a raiz B: rank(D|E) = dim E = 1 só com o testemunho?
  - 0,7332: a injetividade de K (região de S3a) cobre 0,7332? Isso basta para rank(D|E) = 0 em todo o
    espaço de raízes?
  - A periodicidade em Im s, o domínio fundamental e as retas de fronteira estão bem tratadas?
  - Re s = 0, com a raiz de fase em 0, é fronteira?
- **R9 (as hipóteses).**
  - Estão enunciadas com precisão? Há hipóteses escondidas?
  - Separe o que é matemática não feita do que é modelagem física.
  - Diga se "com cone fixo seria 2" está certo.

## Relatório

`.codex-runs/2026-10-02-revisao-T2/RELATORIO.md`:
- uma tabela R1–R9 com veredicto (OK / OK com ressalva / LACUNA / ERRO), artefato e o que foi feito;
- a lista de lacunas, com gravidade e conserto proposto;
- o que você não conseguiu verificar.
