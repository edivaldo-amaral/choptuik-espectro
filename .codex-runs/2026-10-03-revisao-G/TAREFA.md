# Revisão independente do Lema G: completude global do gauge (03/10/2026)

Você é o revisor independente de um lema matemático com um fato certificado por computador. O autor é
outro modelo. Seu trabalho é **tentar quebrar** o lema, não confirmá-lo. Um "OK" sem artefato conferível
não conta: o artefato é um script com saída, ou uma derivação escrita passo a passo no relatório.
Escreva em português do Brasil.

## O objeto da revisão

`docs/GAUGE_GLOBAL.md`: o enunciado, a prova (Passos 0–8), o fato certificado F1
(`scripts/gauge_global_centro.py`, `build/t2/gauge_global_centro.json`) e as observações.

O lema substitui uma hipótese (H-gauge) no teorema condicional T2 (`docs/T2_ENUNCIADO.md`). O uso em T2
é: o autovetor da raiz física s ≈ 0,7332 não é de gauge, e em s = μ* o espaço de raízes, span(g*), é
todo de gauge. Se o lema estiver errado, ou provar menos do que T2 usa, a contagem física 1 perde um
apoio.

## Regras

- **Pode ler:**
  - `docs/GAUGE_GLOBAL.md`, `docs/GAUGE_RT_ACTION.md`, `docs/GAUGE_DOMAIN_CLASSIFICATION.md`,
    `docs/T2_ENUNCIADO.md`, `docs/T2_HIPOTESES.md`, `docs/S4_GAUGE.md`, `docs/SPECTRAL_PROBLEM.md`,
    `docs/SECTOR_B.md`, `docs/RADIUS_RECOVERY.md`, `docs/GEOMETRIC_FLOQUET_RECONSTRUCTION.md` e
    `docs/CONSTRAINT_ROOT_QUOTIENT.md`;
  - os scripts em `scripts/*.py` (por exemplo `gauge_global_centro.py`, `gauge_refB.py`,
    `gauge_null_profiles.py`, `test_gauge_action.py`);
  - os dados de RT em `.cache/rt-1203.3766v1/` (`sourcecode/RefA.dat`, `RefAplusB.dat`).

  Vários documentos foram escritos por outro modelo em 14–16/09: trate-os como **alegações**.
- **O artigo de RT** (arXiv 1203.3766v1) está em `.cache/rt-1203.3766v1/arXiv-1203.3766v1.tar.gz`, com
  `choptuik.tex` dentro. Extraia para `.codex-runs/2026-10-03-revisao-G/rt/`. As seções relevantes:
  - o "Main result", com o domínio M, Θ e a forma da métrica;
  - §"From 4D to 2D", com o ansatz, as definições (ekjehee) dos ω e a invariância conforme;
  - §"From 2D to 4D", com a codificação ω_i^- = ω_i e ω_i^+ = −Pω_i;
  - §"Fourier-Chebyshev series", com a base, as normas e os espaços;
  - os rótulos [dkjfhk2] e [kdkkkskkskksksks], com a bola.
- **NÃO pode ler:**
  - `docs/C1_REAVALIACAO.md`, `docs/PROOF_OBLIGATIONS.md`, `docs/RETOMADA.md`, `docs/PRIORIDADES.md`;
  - `docs/AUDITORIA_*.md`;
  - nada em `.codex-runs/` além desta pasta;
  - a memória do autor.
- **Não altere nenhum arquivo existente.** Escreva só em `.codex-runs/2026-10-03-revisao-G/`.
- **Ambiente:** `.venv/bin/python`, com sympy e mpmath disponíveis. Use **só o notebook**, com
  `OPENBLAS_NUM_THREADS=2`, sem outras máquinas. Não deve haver cálculo longo.
- **Snapshot:** `.snapshots/2026-10-03-revisao-G/MANIFEST.sha256`, com o commit em `COMMIT.txt`.
  Confira-o no início e no fim.
- **Orçamento:** até ~2 h. Grave o RELATORIO.md aos poucos.

## Alegações e artefatos obrigatórios

- **G1 (definições).**
  - "Gauge" está definido corretamente: campo C¹ em M × S² com o eixo, ação por derivada de Lie em
    (g, φ), e "admissível" como preservação da forma linearizada do ansatz.
  - "Modo de Floquet" (δω∘Θ² = e^{4πs}δω) corresponde ao que T2 usa: modos e^{sτ}h na realização
    espectral, setores A e B.
  - Procure hipóteses escondidas. Por exemplo:
    - μ fixo;
    - o domínio é o M de RT e contém o cone;
    - é legítimo exigir X de classe C¹ através do cone;
    - X complexo;
    - perturbações cujo δω está no espaço de raízes mas que não vêm de X Floquet.
- **G2 (Passo 0).** A média em SO(3), o C¹ e a tangência ao eixo.
- **G3 (Passos 1–2).**
  - Derive (L_X g)_{±±} a partir da métrica de RT em coordenadas nulas.
  - Diga se "forma do ansatz" é exatamente δg_{++} = δg_{--} = 0 ou se exige mais (relação de R com ζ e
    Q, positividade de μ + Qξ², regularidade de Q no eixo). Para o lema basta a condição necessária;
    confira que é necessária.
  - Confira a conexidade das curvas de nível de u_± em M (faça as contas com M do artigo) e a classe C¹
    de f_±.
- **G4 (Passo 3).**
  - O eixo dá f_+ = f_- em (−∞, 0)? Confira que a imagem de u_+ em M é (−∞, 0) e que o eixo cobre
    (−∞, 0).
  - Confira as fórmulas de δτ e δξ, simbolicamente.
- **G5 (Passo 4, a equivariância).**
  - Confira (Θ²)*δω_X = δω_{(Θ²)*X}, com a invariância conforme ζ ↦ ζ + c e o fato de Θ preservar o
    ansatz.
  - Confira f̃(u) = Λ⁻¹f(Λu), com Λ = e^{−4πμ}, e a direção do pullback.
  - Teste simbolicamente com f = u^n, usando ou adaptando `gauge_null_profiles.py`: os expoentes devem
    ser μ(1 − n).
- **G6 (Passo 5, a injetividade).**
  - Derive a fórmula X(φ*) = (e^{μτ}/2μ)[f(u_-)ω4^+ − f(u_+)ω4^-], com os sinais de D^±(u_∓).
  - Confira no eixo: ω4^±(τ, 0) = ∓ω4(τ, 0) pela codificação de RT e φ par em ξ.
  - Confira:
    - a antiperiodicidade de ω4 em τ, que dá um zero;
    - a analiticidade em τ, que dá zeros isolados;
    - o argumento no exterior do cone, isto é, que o aberto {ξ > 1, u_- ∈ I} ∩ M é não vazio e que
      ω4^+ ≡ 0 em M contradiz F1.
  - Confira que δω4^σ = D^σ(X(φ*)), com D^σ fixo.
  - Tente um contraexemplo: um f ≠ 0 com δω_X = 0.
- **G7 (F1, o fato certificado).**
  - Recalcule c_1 = Σ_n v_{1n}i^n de ω4 em ξ = 0 com **código próprio** a partir de `RefA.dat` (não
    reutilize `gauge_global_centro.py`), em racionais exatos.
  - Confira no artigo:
    - a base e^{imτ/2 + inθ} com θ = arccos ξ, de modo que ξ = 0 corresponde a θ = π/2;
    - a simetria v_{m,−n} = v_{mn};
    - que o campo de índice 3 é ω4;
    - a norma da bola e a cota |c_1(ω*) − c_1(RefA)| <= 2⁻²⁵ + 2⁻²⁷⁷.

    Em especial: a bola é para ω* − (RefA + RefB)? ‖RefB‖ <= 2⁻²⁵ está nessa mesma norma? Os pesos são
    >= 1?
- **G8 (Passos 6–7).** Confira a análise de casos para f: ℝ → ℂ, C¹, com f(Λu) = ρf(u), nos três
  regimes de Re s. Os casos estão completos? Confira também:
  - ρ = 1 ⇔ s ≡ μ mod (i/2)ℤ;
  - a cobertura dos setores A e B;
  - o caso Re s = μ com ρ ≠ 1.
- **G9 (Passo 8 e consequências).**
  - Confira que f ≡ 1 dá X0 = e^{μτ}(μ⁻¹∂_τ + ξ∂_ξ) e δω_{X0} = e^{μτ}(Z + 1)ω*, comparando com
    `GAUGE_RT_ACTION.md`, e rode `scripts/test_gauge_action.py`.
  - Confira que o uso em `T2_ENUNCIADO.md` (passo 11b e a "Conta") é exatamente o que o lema prova.
  - Diga se a Observação 1 (o papel do cone) e a Observação 2 (H-cone) estão corretas.

## Relatório

`.codex-runs/2026-10-03-revisao-G/RELATORIO.md`:
- uma tabela G1–G9 com veredicto (OK / OK com ressalva / LACUNA / ERRO), o artefato e o que foi feito;
- a lista de lacunas, com gravidade e conserto proposto;
- o que você não conseguiu verificar;
- a conferência do manifesto no início e no fim.
