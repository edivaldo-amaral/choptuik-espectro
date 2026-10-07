# Revisão independente da H-rec (ida): Lema R (03/10/2026)

Você é o revisor independente de um lema matemático, quase todo em papel. O autor é outro modelo. Seu
trabalho é **tentar quebrar** o lema, não confirmá-lo. Um "OK" sem artefato conferível não conta: o
artefato é um script com saída, ou uma derivação escrita passo a passo no relatório.
Escreva em português do Brasil.

## O objeto da revisão

`docs/HREC_IDA.md`: o Lema R, que diz que todo modo físico de Floquet com Re s > 0 entra no ansatz de RT
por um gauge analítico, com as variáveis de RT no espaço espectral Y+. A prova tem seis peças, R1–R6.

**O que está em jogo.** O teorema condicional T2 (`docs/T2_ENUNCIADO.md`) conta modos instáveis físicos
como (ker L(s) ∩ ker C(s)) módulo gauge. A H-rec (ida) é a hipótese de que todo modo físico aparece
nessa conta. Se o Lema R estiver certo, a hipótese vira teorema, salvo a modelagem H-reg e H-cone.

**Pontos de maior risco:**
- **R2, a ressonância s ≡ μ.** O autor argumenta que, se a compatibilidade falhasse, o gauge secular κτX0
  produziria uma cadeia de Jordan de L em μ*, e que isso contradiz a simplicidade certificada em S4.
  Confira cada passo: a decomposição δ' = A + τB, a forma de δω', as equações da cadeia, o espaço em que
  a cadeia vive e a identificação de μ* com s = μ.
- **R5, a recuperação de raio a partir de qualquer peso > 1.** Confira a validade de `free_tail` para
  peso radial genérico, a reponderação β(κ2') = β+·D(κ2')/D(5/4) e a cota explícita.
- **R1 e R4,** que vêm de provas do interlúdio de 14–16/09 em `GLOBAL_NULL_GAUGE_TRANSPORT.md` e
  `GEOMETRIC_FLOQUET_RECONSTRUCTION.md`. Não foram revisadas de forma independente: trate como
  alegações.
- **R3,** a regularidade de δQ no eixo, que está em esboço.

## Regras

- **Pode ler:**
  - `docs/HREC_IDA.md`, `docs/GLOBAL_NULL_GAUGE_TRANSPORT.md`, `docs/GEOMETRIC_FLOQUET_RECONSTRUCTION.md`,
    `docs/RADIUS_RECOVERY.md`, `docs/STRUCTURED_EXTERIOR.md`, `docs/SHARP_EXTERIOR_BOUNDS.md`,
    `docs/COMPONENT_EXTERIOR.md`, `docs/GAUGE_GLOBAL.md`, `docs/GAUGE_RT_ACTION.md`,
    `docs/GAUGE_DOMAIN_CLASSIFICATION.md`, `docs/S4_GAUGE.md`, `docs/T2_ENUNCIADO.md`,
    `docs/T2_HIPOTESES.md`, `docs/SPECTRAL_PROBLEM.md`, `docs/SECTOR_B.md`, `docs/SHARP_DOMAIN.md` e
    `docs/CONSTRAINT_ROOT_QUOTIENT.md`;
  - os scripts em `scripts/*.py` (por exemplo `radius_recovery_geral.py`, `sharp_exterior.py`,
    `structured_exterior.py`, `radius_transfer.py`, `global_null_gauge.py`, `test_global_null_gauge.py`,
    `test_geometric_reconstruction.py` e `gauge_null_profiles.py`);
  - `build/t2/**` e os dados de RT em `.cache/rt-1203.3766v1/`.
- **O artigo de RT** está em `.cache/rt-1203.3766v1/arXiv-1203.3766v1.tar.gz`, com `choptuik.tex` dentro.
  Extraia para o scratchpad da sua sessão (**não** para esta pasta: tem 168 MB).
- **NÃO pode ler:**
  - `docs/C1_REAVALIACAO.md`, `docs/PROOF_OBLIGATIONS.md`, `docs/RETOMADA.md`, `docs/PRIORIDADES.md`;
  - `docs/AUDITORIA_*.md`;
  - nada em `.codex-runs/` além desta pasta;
  - a memória do autor.
- **Não altere nenhum arquivo existente.** Escreva só em `.codex-runs/2026-10-03-revisao-R/`.
- **Ambiente:** `.venv/bin/python`.
  - O `.venv` **não** tem sympy nem mpmath. Se precisar, instale com `pip install --target <scratchpad>`
    e use `PYTHONPATH`.
  - Os testes rodam com `cd scripts && ../.venv/bin/python -m unittest -v <modulo>`, sem pytest.
- **Use só o notebook,** com `OPENBLAS_NUM_THREADS=2` e sem ssh. Não deve haver cálculo longo.
- **Snapshot:** `.snapshots/2026-10-03-revisao-R/MANIFEST.sha256`, com o commit em `COMMIT.txt`.
  Confira-o no início e no fim.
- **Orçamento:** até ~2,5 h. Grave o RELATORIO.md aos poucos.

## Alegações e artefatos obrigatórios

- **R0 (enunciado e definições).**
  - Diga se a definição de "modo físico de Floquet" é razoável e se esconde hipóteses. Por exemplo:
    - "Floquet num gauge analítico";
    - o fator e^{−4K};
    - analiticidade em M̂ inteiro, incluindo eixo e cone.
  - Diga se a classe de gauge (C¹ ou G′) está declarada e é coerente com `GAUGE_GLOBAL.md`.
  - Diga se o item 3 e a bijeção do §1 seguem das peças.
- **R1 (transporte).**
  - Derive T_{±1}(s) a partir de ∂_{u_+}X^{u_-} = −δg_{++}/(2g_{+-}), com X^± = e^{(s−μ)τ}y_±, e confira o
    sinal e o fator de f_∓.
  - Confira a fórmula (*) de `GLOBAL_NULL_GAUGE_TRANSPORT.md`: que resolve a equação, que converge para
    0 < Re s < μ e que preserva a analiticidade num domínio convexo.
  - Confira a unicidade fora da ressonância e a reflexão y_+ = Py_- com a regularidade no eixo, isto é,
    A par e B ímpar.
  - Confira a colagem analítica de f_- com f_+ refletida em ξ = 0, que é a regularidade cartesiana do
    tensor δg.
  - Rode `test_global_null_gauge.py` e acrescente pelo menos um teste próprio, por exemplo um modo com
    fonte analítica não polinomial, resolvido numericamente por (*) e substituído na equação.
- **R2 (ressonância).** Derive passo a passo:
  - T_1(s0)(τe^{im0τ/2}) = e^{im0τ/2};
  - a parte secular de X é κτX0, com a mesma constante nas componentes + e −;
  - δ' = A + τB, com A e B no ansatz, e A Floquet;
  - δω' = e^{s0τ}(h1 + τh0), com h0 = κe^{im0τ/2}g*_perfil;
  - L(s0)h0 = 0 e L(s0)h1 = −h0, usando L(s) = L(0) + s (L-afim);
  - a contradição com a simplicidade algébrica de μ* (S4). Confira em `S4_GAUGE.md` o que foi
    certificado, em qual realização e em qual disco, e se o argumento vale para m0 ≠ 0.

  Procure furos: termos com τ que escapem, R4 aplicado a perturbações seculares, e a realização em que h1
  vive.
- **R3 (extração de δζ, δQ e δφ).** Confira as fórmulas a partir da métrica de RT e a regularidade de δQ
  no eixo. Escreva a prova completa de δζ + δb/(2b) = O(ξ²), função par, ou ache a falha.
- **R4 (equações).**
  - Confira em `GEOMETRIC_FLOQUET_RECONSTRUCTION.md` e no artigo de RT que, para uma perturbação no
    ansatz que resolve as equações linearizadas, todos os δΩ se anulam.
  - Rode `test_geometric_reconstruction.py` e diga o que exatamente ele testa.
  - Diga se R4 vale para perturbações seculares, que R2 usa.
- **R5 (recuperação de raio).**
  - Rode `radius_recovery_geral.py`.
  - Recalcule free_tail(M, N, κ2') de forma independente para pelo menos dois pesos.
  - Confira a cota explícita, isto é, a desigualdade free_tail <= cota, e o limite.
  - Confira D(k) = 2k/(k² − 1) no artigo, onde está o rótulo da divisão regularizada.
  - Confira o argumento de existência forte e unicidade fraca para cadeias.
  - Confira que um perfil analítico Floquet está em algum Y_{κ'}, junto com J_μh.
- **R6 (unicidade módulo gauge).** Confira os três itens, com o uso do Lema G ou G′.

## Relatório

`.codex-runs/2026-10-03-revisao-R/RELATORIO.md`:
- uma tabela R0–R6 com veredicto (OK / OK com ressalva / LACUNA / ERRO), o artefato e o que foi feito;
- a lista de lacunas, com gravidade e conserto proposto;
- o que você não conseguiu verificar;
- a conferência do manifesto no início e no fim.
