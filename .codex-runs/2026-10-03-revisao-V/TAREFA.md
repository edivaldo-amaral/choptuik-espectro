# Revisão independente da volta da reconstrução (03/10/2026)

Você é o revisor independente de uma prova matemática em papel, com testes racionais. O autor é outro
modelo. Seu trabalho é **tentar quebrar** a prova, não confirmá-la. Um "OK" sem artefato conferível não
conta: o artefato é um script com saída, ou uma derivação escrita passo a passo no relatório.
Escreva em português do Brasil.

## O objeto da revisão

`docs/GEOMETRIC_FLOQUET_RECONSTRUCTION.md`, a **volta**. A alegação é que todo h com L(s)h = 0 e
C(s)h = 0, com Re s > 0, vem de uma perturbação física: uma solução e^{sτ}(δQ, δζ, δφ) das equações de
Einstein–campo escalar linearizadas em torno da solução de Reiterer–Trubowitz (RT), na forma do ansatz,
regular no eixo e analítica através do cone.

O documento foi escrito no interlúdio de 14–16/09 por outro modelo e auditado em 16/09, mas nunca
passou por revisão independente. Trate tudo como **alegação**.

**O que está em jogo.** O teorema condicional T2 (`docs/T2_ENUNCIADO.md`) afirma "exatamente 1 modo
instável físico". A ida (`docs/HREC_IDA.md`, Lema R, já revisada) garante que todo modo físico aparece em
ker L ∩ ker C. A volta garante que o modo de ker L ∩ ker C em s ≈ 0,7332 é de fato físico, isto é, que a
contagem 1 é atingida. A parte recíproca usada na ida (R4: perturbação no ansatz que resolve Einstein ⇒
δΩ = 0) já foi revisada. Aqui o foco é a direção oposta.

## Regras

- **Pode ler:**
  - `docs/GEOMETRIC_FLOQUET_RECONSTRUCTION.md`, `docs/GLOBAL_NULL_GAUGE_TRANSPORT.md`, `docs/HREC_IDA.md`,
    `docs/GAUGE_GLOBAL.md`, `docs/CONSTRAINT_ROOT_QUOTIENT.md`, `docs/SHARP_LINEARIZED_IDENTITY.md`,
    `docs/SHARP_DOMAIN.md`, `docs/SPECTRAL_PROBLEM.md`, `docs/SECTOR_B.md`, `docs/RADIUS_RECOVERY.md`,
    `docs/T2_ENUNCIADO.md`, `docs/T2_HIPOTESES.md` e `docs/S4_GAUGE.md`;
  - os scripts em `scripts/*.py` (por exemplo `geometric_reconstruction.py`,
    `test_geometric_reconstruction.py`, `signed_constraints.py`, `signed_operator_L.py`,
    `test_root_constraints.py` e `independent_rt_symbolic.py`, se existir);
  - `build/t2/**` e os dados e o código de RT em `.cache/rt-1203.3766v1/`.
- **O artigo de RT** está em `.cache/rt-1203.3766v1/arXiv-1203.3766v1.tar.gz`, com `choptuik.tex` dentro.
  Extraia para o scratchpad da sua sessão (não para esta pasta: tem 168 MB). As seções relevantes:
  - "From 4D to 2D", com o ansatz, as definições (ekjehee) dos ω, as relações (dhjhfkfhfe1–4), as
    fórmulas de curvatura (dfkhdjhsdshkfd), os Ω (eezuuirzr) e o "Claim";
  - "From 2D to 4D", com a codificação ω_i^+ = −Pω_i e o problema 2Dprob;
  - "The ♯ system", com as identidades diferenciais e a seleção;
  - "Fourier-Chebyshev series", com a base e os espaços.
- **NÃO pode ler:**
  - `docs/C1_REAVALIACAO.md`, `docs/PROOF_OBLIGATIONS.md`, `docs/RETOMADA.md`, `docs/PRIORIDADES.md`;
  - `docs/AUDITORIA_*.md`;
  - nada em `.codex-runs/` além desta pasta;
  - a memória do autor.
- **Não altere nenhum arquivo existente.** Escreva só em `.codex-runs/2026-10-03-revisao-V/`.
- **Ambiente:** `.venv/bin/python`.
  - O `.venv` **não** tem sympy nem pytest. Se precisar de sympy, use `pip install --target <scratchpad>` e
    `PYTHONPATH`.
  - Os testes rodam com `cd scripts && ../.venv/bin/python -m unittest -v <modulo>`.
- **Use só o notebook,** com `OPENBLAS_NUM_THREADS=2` e sem ssh.
- **Snapshot:** `.snapshots/2026-10-03-revisao-V/MANIFEST.sha256`, com o commit em `COMMIT.txt`. Confira-o
  no início e no fim.
- **Orçamento:** até ~2,5 h. Grave o RELATORIO.md aos poucos.

## Alegações e artefatos obrigatórios

- **V0 (enunciado).**
  - Escreva com precisão o que a volta afirma:
    - em qual espaço h vive (D(Y+), ou o toro com sinal, onde as contagens são certificadas, ligados pelo
      lema L-real de `T2_HIPOTESES.md`);
    - a classe da perturbação reconstruída;
    - o sentido de "Floquet", com o fator e^{−4K} do fundo, porque Θ*ḡ = e^{−2K}ḡ.
  - Diga se há hipóteses escondidas.
- **V1 (seleção + constraints = sistema completo).**
  - Mostre que L(s)h = 0 e C(s)h = 0 implicam que os dez resíduos linearizados δΩ_i^σ se anulam.
  - Use as definições de L, que é a seleção SΩ^+ com Ξ⁻¹ regularizado, e de C, em `signed_constraints.py`
    e `CONSTRAINT_ROOT_QUOTIENT.md`, junto com o argumento do sistema ♯ de RT linearizado.
  - Confira que nada se perde na divisão regularizada por ξ nem nas projeções de paridade.
  - O documento diz: "O sistema linearizado COMPLETO significa todos os resíduos delta Omega^sigma nulos,
    não só a seleção." Confira de onde isso sai.
- **V2 (potenciais).**
  - Confira as definições de a_i e b_i e a equivalência D_s^σ z = h_i^σ ⇔ (T_s z = a_i e ∂_ξ z = b_i).
  - Confira a fórmula de T_s⁻¹, com a convergência e a preservação de paridades para Re s > 0, e a
    ausência de obstrução de período.
  - Derive simbolicamente as identidades de compatibilidade −2μξ(T_s b_i − ∂_ξ a_i) = δΩ_j^- − δΩ_j^+,
    para (i, j) = (3, 3), (4, 4) e (2, 2), a partir dos Ω de RT. Diga se são identidades fora das
    soluções.
- **V3 (reconstrução e log W).**
  - Confira as fórmulas q, z e p.
  - Derive a identidade algébrica D^σW − W(ω2^σ − ω3^σ) = (Ω1^σ − Ω2^σ − Ω2♯^σ)/2 e a sua linearização
    no fundo exato.
  - Confira que h1, h2, h3 e h4 são de fato recuperados a partir de (q, z, p).
- **V4 (regularidade).**
  - Confira a paridade: q ímpar no numerador e z, p pares.
  - Confira o centro removível e a analiticidade através do cone, que a volta alega obter sem integral
    radial atravessando ξ = 1.
  - Confira W > 0.
  - Confira que a elipse de Bernstein do peso 5/4 tem semieixo (5/4 + 4/5)/2 = 41/40. Isso significa que
    h ∈ Y+ dá analiticidade em todo o domínio M de RT (0 <= ξ < 41/40)? Confira.
- **V5 (a métrica).** Confira a fórmula de δg⁻¹ pelo frame de RT, com δe_0 = 0 e o perfil de δe_i, e o
  δφ. Confira também o sentido de Floquet com o fator e^{−4K}.
- **V6 (as equações geométricas).** Pelas fórmulas de curvatura (dfkhdjhsdshkfd), δΩ = 0 no fundo exato
  implica δ(Ric − 2dφ⊗dφ) = 0 e δ(□φ) = 0.
  - Confira, inclusive no centro.
  - Rode `test_geometric_reconstruction.py` e diga o que ele testa.
  - Se puder, confira as fórmulas de curvatura de RT com um cálculo do Ricci 4D direto da métrica.
- **V7 (unicidade e cadeias).**
  - Confira que o núcleo da reconstrução é trivial para Re s > 0 e que as constantes invisíveis
    (δζ e δφ constantes) não geram perfis Floquet com Re s > 0.
  - Confira a recorrência para jatos e cadeias.
- **V8 (consequência para T2).** Diga se, com a volta, o vetor de ker L ∩ ker C em s ≈ 0,7332 dá uma
  perturbação física não nula e não de gauge, pelo Lema G de `GAUGE_GLOBAL.md`, na classe H-reg de T2. E
  se a correspondência do §1 de `HREC_IDA.md` fica bijetiva.

## Relatório

`.codex-runs/2026-10-03-revisao-V/RELATORIO.md`:
- uma tabela V0–V8 com veredicto (OK / OK com ressalva / LACUNA / ERRO), o artefato e o que foi feito;
- a lista de lacunas, com gravidade e conserto proposto;
- o que você não conseguiu verificar;
- a conferência do manifesto no início e no fim.
